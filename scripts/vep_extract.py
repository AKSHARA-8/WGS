#! /usr/bin/env python3

import sys
import csv
import re

csv_file=sys.argv[1]	# csv file from the merge_csv process
vep_output = sys.argv[2]# Output file from VEP

map_id = dict()			# Dictionary for the marker occurence
map_cdna_pos = dict()			# Dictionary for cdna_pos
map_cds_pos = dict()			# 
map_protein_pos = dict()		# 
map_amino_acids = dict()		
map_codons = dict()
map_hgvsc = dict()
map_hgvsp = dict()

count = 0
with open (vep_output, 'r') as vep_file:
	vep_handle = csv.reader(vep_file, delimiter = '\t')
	header = next (vep_handle)
	for vep_values in vep_handle:
		chrom = vep_values[1].split(':')[0]
		pos = vep_values[1].split(':')[1].split('-')[0]
		alt = vep_values[2]
		cdna_pos = vep_values[7]
		cds_pos  = vep_values[8]
		protein_pos = vep_values[9]
		amino_acid = vep_values[10]
		codons = vep_values[11]
		hgvsc = vep_values[19]
		hgvsp = vep_values[20]

		if re.search('[a-zA-Z]', str(chrom)):
			chrom = re.sub ("chr","", chrom, flags = re.IGNORECASE)
			chrom = re.sub ("X","23", chrom, flags = re.IGNORECASE)
			chrom = re.sub ("Y","y", chrom, flags = re.IGNORECASE)
		
		vep_id = str(chrom) + ':' +''.join(str(pos)) + ':' +''.join(str(alt))
		map_id[count] = vep_id
		map_cdna_pos[count] = cdna_pos
		map_cds_pos[count] = cds_pos
		map_protein_pos[count] = protein_pos
		map_amino_acids[count] = amino_acid
		map_codons[count] = codons
		map_hgvsc[count] = hgvsc
		map_hgvsp[count] = hgvsp
		count = count + 1
	

with open (csv_file,'r') as vcf:
	csv_handle = csv.reader(vcf)
	header = next (csv_handle)      # Removing the header
	print (*header,"cDNA_position_VEP", "CDS_position_VEP",	"Protein_position_VEP",	"Amino_acids_VEP", "Codons_VEP", "HGVSc_VEP", "HGVSp_VEP", sep=",")
	for str_lines in csv_handle:
		chromosome = str_lines[0]
		position = str_lines[1]
		ref = str_lines[3]
		alt = str_lines[4]
		
		if re.search('[a-zA-Z]', str(chromosome)):
			chromosome = re.sub ("chr","", chromosome, flags = re.IGNORECASE)
			chromosome = re.sub ("X","23", chromosome, flags = re.IGNORECASE)
			chromosome = re.sub ("Y","y", chromosome, flags = re.IGNORECASE)
		
		cdna_pos_list = list()
		cds_pos_list = list()
		protein_pos_list = list()
		amino_acids_list = list()
		codons_list = list()
		hgvsc_list = list()
		hgvsp_list = list()
		csv_id = str(chromosome) + ':' +''.join(str(position)) + ':' +''.join(str(alt))
		for key, values in map_id.items():
			if csv_id in values:
				#print (map_cdna_pos[key], map_cds_pos[key], map_protein_pos[key], map_amino_acids[key], map_codons[key], map_hgvsc[key], map_hgvsp[key], sep=",")
				cdna_pos_list.append(map_cdna_pos[key])
				cds_pos_list.append(map_cds_pos[key])
				protein_pos_list.append(map_protein_pos[key])
				amino_acids_list.append(map_amino_acids[key])
				codons_list.append(map_codons[key])
				hgvsc_list.append(map_hgvsc[key])
				hgvsp_list.append(map_hgvsp[key])


		
		str_cdna_pos_list = str(cdna_pos_list)[1:-1]	# This removes the brackets 
		str_cdna_pos_list = re.sub ("'","", str_cdna_pos_list, flags = re.IGNORECASE)
		str_cds_pos_list = str(cds_pos_list)[1:-1]
		str_cds_pos_list = re.sub ("'","", str_cds_pos_list, flags = re.IGNORECASE)
		str_protein_pos_list = str(protein_pos_list)[1:-1]
		str_protein_pos_list = re.sub ("'","", str_protein_pos_list, flags = re.IGNORECASE)
		str_amino_acids_list = str(amino_acids_list)[1:-1]
		str_amino_acids_list = re.sub ("'","", str_amino_acids_list, flags = re.IGNORECASE)
		str_codons_list = str(codons_list)[1:-1]
		str_codons_list = re.sub ("'","", str_codons_list, flags = re.IGNORECASE)
		str_hgvsc_list = str(hgvsc_list)[1:-1]
		str_hgvsc_list = re.sub ("'","", str_hgvsc_list, flags = re.IGNORECASE)	
		str_hgvsp_list = str(hgvsp_list)[1:-1]
		str_hgvsp_list = re.sub ("'","", str_hgvsp_list, flags = re.IGNORECASE)

		print (*str_lines, str_cdna_pos_list.replace(',',' '), str_cds_pos_list.replace(',',' '), str_protein_pos_list.replace(',',' '), str_amino_acids_list.replace(',',' '), str_codons_list.replace(',',' '), str_hgvsc_list.replace(',',' '), str_hgvsp_list.replace(',',' '), sep=",")
