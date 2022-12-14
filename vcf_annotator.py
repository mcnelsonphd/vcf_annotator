import requests
import sys, os
import vcf
import json
import csv


vcfjson = {}

def VEP(vlocus):
	"""VEP API handler, built to specifically query the GRCh37 human refernce build and only 
	return desired annotation data. Returning the 
	"""
	anno_data = {}
	endpoint = "http://grch37.rest.ensembl.org/vep/human/hgvs/"
	headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
	var_data = f'{{"hgvs_notations":[{vlocus}]}}'
	r = requests.post(endpoint, 
		headers=headers, 
		data=var_data)

	if not r.ok:
		r.raise_for_status()
		sys.exit()

	decoded = r.json()

	if "colocated_variants" in decoded[0].keys():
		for dataset in decoded[0]["colocated_variants"]:
			try:   # if 'minor_allele' in dataset.keys() we want that data
				anno_data["MinorVar"] = dataset['minor_allele']
				anno_data["MAF"] = dataset['minor_allele_freq']
			except: # handle key error or 'minor_allele' not being present, ie move on.
				pass

	if "transcript_consequences" in decoded[0].keys():
		gene = []
		for cons in decoded[0]["transcript_consequences"]: 
			gene.append(cons["gene_symbol"]) 
		anno_data["Gene"] = max(set(gene), key=gene.count)
	else: 
		anno_data["Gene"] = "N/A"
	
	anno_data["Effect"] = decoded[0]["most_severe_consequence"]			
	return(anno_data)

# Build a dictionary of observed records, splitting multi-allelic, and "correcting" pos value to make VEP g. compliant
for record in vcf.Reader(open(sys.argv[1], 'r')):
	fields = ["Chrom", "Pos", "Ref", "Alt", "VEP_ID", "TCov", "VCov"]
	if len(record.REF) > 1: # For VEP, the locus needs to be provided as the range of bases covered by the variant
		locus_pos = f"{record.POS}_{record.POS+len(record.REF)-1}"
	else:
		locus_pos = record.POS
	if len(record.ALT) == 1: # Standard Biallelic variant
		var_loci=f'"{record.CHROM}:g.{locus_pos}{record.REF}>{record.ALT[0]}"'
		vcfjson[var_loci.strip('"')] = dict(zip(fields, [record.CHROM, record.POS, record.REF, record.ALT[0], var_loci, record.INFO["TC"], record.INFO["TR"][0]]))
	else: # Dealing with a multi-allelic variant, examination of 
		vlocus=f'"{record.CHROM}:g.{locus_pos}{record.REF}>{record.ALT[0]}"'
		vcfjson[vlocus.strip('"')] = dict(zip(fields, [record.CHROM, record.POS, record.REF, record.ALT[0], vlocus, record.INFO["TC"], record.INFO["TR"][0]]))
		vlocus2=f'"{record.CHROM}:g.{locus_pos}{record.REF}>{record.ALT[1]}"'
		vcfjson[vlocus2.strip('"')] = dict(zip(fields, [record.CHROM, record.POS, record.REF, record.ALT[1], vlocus2, record.INFO["TC"], record.INFO["TR"][1]]))
	
# Iterate through the variant dictionary to get VEP annotations and calculate other desired values: variant type, allele freq.
for var in vcfjson.keys():
	vcfjson[var]["AF"] = vcfjson[var]["VCov"] / vcfjson[var]["TCov"]
	if len(vcfjson[var]['Alt']) == len(vcfjson[var]['Ref']):
		vcfjson[var]['VarType'] = "Substitution"
	elif len(vcfjson[var]['Alt']) < len(vcfjson[var]['Ref']):
		vcfjson[var]['VarType'] = "Deletion"
	elif len(vcfjson[var]['Alt']) > len(vcfjson[var]['Ref']):
		vcfjson[var]['VarType'] = "Insertion"
	else:
		vcfjson[var]['VarType'] = "Unknown"
	ret = VEP(vcfjson[var]["VEP_ID"])
	vcfjson[var]["Effect"] = ret["Effect"]
	vcfjson[var]["Gene"] = ret["Gene"]
	vcfjson[var]["MAF"] = ret["MAF"] if "MAF" in ret.keys() else ''
	vcfjson[var]["MinorAllele"] = ret["MinorVar"] if "MinorVar" in ret.keys() else ''

# Write output in tsv format, ordering according to header
with open(sys.argv[2], 'w') as out:
	outheader = ["Chrom", "Pos", "Ref", "Alt", "TCov", "VCov", "AF", "VEP_ID", "Gene", "VarType", "Effect", "MinorAllele", "MAF"]
	writer = csv.DictWriter(out, fieldnames = outheader, delimiter = '\t')
	writer.writeheader() 
	for variant in sorted(vcfjson):
		writer.writerow(vcfjson[variant])
