#Author: Boyang (Betty) Xie
#Date: 11/2015
#email: boyang.xie@vanderbilt.edu
"""This is the main script used for calculations that compare the different possibilities of 5 subpopulations whose protein sequences are different from those
   in the PDB. It gives out the following numbers: for each structure id in the PDB, the different mismatch ratios of different subpopulations; and the overall
   ratio of mismatch for all the structures in different subpopulations (give a overall view of bias)"""
"""Note that this only focuses on SNPs, not including insertions and deletions."""
"""This script count by mismatched snp numbers"""


#!/usr/bin/env python

import sys
import vcf

VCF_FILE = "/dors/capra_lab/xieb2/1kg3_chr2pdb_nosingletons.vcf" #sys.argv[1]   10_indv_1kg3_chr2pdb_nosingletons.vcf
# VCF_FILE = "/dors/capra_lab/xieb2/10_indv_1kg3_chr2pdb_nosingletons.vcf"
PDB_INFO_FILE = "/dors/capra_lab/projects/pdb_ancestry_bias/data/1kg3_chr2pdb_nosingletons_tonyedit.txt" #sys.argv[2]
POPULATION_FILE = "/dors/capra_lab/xieb2/integrated_call_samples_v3.20130502.ALL.panel"

################################################################################
# HELPER FUNCTIONS
################################################################################

complement = {"A": "T", "T": "A", "C": "G", "G": "C"}

subpops = ["EAS", "AFR", "EUR", "SAS", "AMR"]


def diff(a,b):
	return sum(a[i] != b[i] for i in range(len(a)))
	
def same(a,b):
	return sum(a[i] == b[i] for i in range(len(a)))


################################################################################


### Read VCF File

snp_id2individual_id2nt = {}
comvcf_ref = {}
vcf_reader = vcf.Reader(open(VCF_FILE, 'r'))
#print vcf_reader
individual_ids = vcf_reader.samples
for record in vcf_reader:
	chr = record.CHROM
	pos_1 = 1 + int(record.start)
	pos_2 = 1 + pos_1
	snp_id = "chr" + chr + ":" + str(pos_1) + "-" + str(pos_2)
	vcf_ref = record.REF
	vcf_alt = record.ALT
	#print "\n", snp_id

	if len(vcf_ref) == 1 and len(vcf_alt) == 1:
		snp_id2individual_id2nt[snp_id]={}
		comvcf_ref[snp_id] = vcf_ref

		for indiv_id in individual_ids:
			call = record.genotype(indiv_id)
			nucleotides = call.gt_bases
			#print "\t", call, nucleotides
			snp_id2individual_id2nt[snp_id][indiv_id] = [nucleotides]


### Read in PDB Info
pdb_chain2snp_id = {}
pdb_chain2dna = {}

compdb_ref = {}

with open(PDB_INFO_FILE,'r') as f:
	for line in f:
		#print line
		line = line.rstrip()
		columns = line.split("\t")

		if columns[0] == "label":
			continue

		chr = columns[1]
		pos_1 = columns[2]
		pos_2 = columns[3]
		PDB_ID = columns[5]
		CHAIN = columns[6]

		PDB_aa = columns[8]
		ref_aa = columns[9]
		alt_aa = columns[10]

		ref_codon = columns[11]
		alt_codon = columns[12]

		ref_nt = "".join([ntr for ntr in ref_codon if ntr.isupper()])
		alt_nt = "".join([nta for nta in alt_codon if nta.isupper()])
		PDB_ID_CHAIN = PDB_ID + '_' + CHAIN # key
		snp_id = chr + ":" + pos_1 + "-" + pos_2

		if len(ref_nt) == 1 and len(alt_nt) == 1:
			compdb_ref0 = {snp_id:ref_nt}
			compdb_ref.update(compdb_ref0)

			if PDB_aa == ref_aa:
				pdb_chain2snp_id.setdefault(PDB_ID_CHAIN, []).append(snp_id)
				pdb_chain2dna.setdefault(PDB_ID_CHAIN,[]).append(ref_nt)
			elif PDB_aa == alt_aa:
				pdb_chain2snp_id.setdefault(PDB_ID_CHAIN, []).append(snp_id)
				pdb_chain2dna.setdefault(PDB_ID_CHAIN,[]).append(alt_nt)
			else:
				print >> sys.stderr, "Warning: PDB_AA (%s) does not match %s or %s for %s" % (PDB_aa, ref_aa, alt_aa, PDB_ID_CHAIN)
f.close()


###Read in population info
indiv_id2supop = {}

with open(POPULATION_FILE,'r') as pop:
	for line in pop:
		#print line
		line = line.rstrip()
		columns = line.split("\t")

		if columns[0] == "sample":
			continue
		indiv_id = columns[0]
		superpop = columns[2]
		indiv_id2supop_0 = {indiv_id:superpop}
		indiv_id2supop.update(indiv_id2supop_0)

pop.close()
		

#compare pdb_nt to indiv_nt for each structure:
print "\n\n*******"

subpop2allpdb_match = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}
subpop2allpdb_mismatch = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}
subpop2allpdb_total = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}
cal = open("/dors/capra_lab/xieb2/cal_all_snp.txt", "w")

allpdb_nm = 0
allpdb_total = 0

for pdb_id in pdb_chain2snp_id:

	pdb_dna = pdb_chain2dna[pdb_id]
	pdb_snp_ids = pdb_chain2snp_id[pdb_id]

	n_m = 0
	n_nm = 0

	subpop2num_match = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}
	subpop2num_mismatch = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}
	total_pop = {"EAS": 0, "EUR": 0, "AFR": 0, "SAS": 0, "AMR": 0}

	if len(pdb_snp_ids) != len(pdb_dna):
		print "WARNING: Different length:", pdb_id, pdb_dna, pdb_snp_ids
		continue

	pdb_seq = ''.join(pdb_dna)
	print "\n", pdb_id, "\n", pdb_seq
	print >>cal, "\n", pdb_id, "\n", pdb_seq
	
	for indiv_id in individual_ids:
		indiv_dna = []
		allele_1 = []
		allele_2 = []


		for snp_id in pdb_snp_ids:
			
			if snp_id not in snp_id2individual_id2nt.keys():
				#print snp_id
				continue

			alleles = snp_id2individual_id2nt[snp_id][indiv_id]
			for allele in alleles:
				al_1, al_2 = allele.split("|")

				if comvcf_ref[snp_id] != compdb_ref[snp_id]:
					#print al_1, al_2
					al_1 = complement[al_1]
					al_2 = complement[al_2]
					#print al_1, al_2
				allele_1.append(al_1)
				allele_2.append(al_2)

			indiv_dna.append(snp_id2individual_id2nt[snp_id][indiv_id])


		if len(indiv_dna) == len(pdb_dna):
			alseq_1 = ''.join(allele_1)
			alseq_2 = ''.join(allele_2)
			
			print alseq_1,"\t", alseq_2 #indiv_dna
			
			if len(pdb_seq) != len(alseq_2) or len(pdb_seq) != len(alseq_1):
				#print "DUANG!", pdb_id, "\t",pdb_seq, "\t",alseq_1, "\t", alseq_2
				continue
			elif len(pdb_seq) == len(alseq_2) and len(pdb_seq) == len(alseq_1):

				indiv_subpop = indiv_id2supop[indiv_id]
		
				oneindv_nm = diff(pdb_seq,alseq_1) + diff(pdb_seq,alseq_2)
				oneindv_m = same(pdb_seq,alseq_1) + same(pdb_seq,alseq_2)

				subpop2num_match[indiv_subpop] += oneindv_m
				subpop2num_mismatch[indiv_subpop] += oneindv_nm

				n_m += oneindv_m
				n_nm += oneindv_nm
			
			else:
				print "duang"

		else:
			#print "\tWARNING:", indiv_dna
			continue

	total = n_m + n_nm

	for subpop in subpops:
		total_pop[subpop] = subpop2num_match[subpop] + subpop2num_mismatch[subpop]

	print "Overall: ",n_nm,"/",total
	print >>cal, "Overall: ",n_nm,"/",total
	for subpop in subpops:
		print subpop, ": ", subpop2num_mismatch[subpop],"/",total_pop[subpop], "\t", subpop2num_mismatch[subpop]/float(total_pop[subpop])
		print >>cal, subpop, ": ", subpop2num_mismatch[subpop],"/",total_pop[subpop], "\t", subpop2num_mismatch[subpop]/float(total_pop[subpop])

	allpdb_nm += n_nm
	allpdb_total += total

	for subpop in subpop2allpdb_mismatch.keys():
		subpop2allpdb_mismatch[subpop] += subpop2num_mismatch[subpop]
		subpop2allpdb_total[subpop] += total_pop[subpop]

print   "\n","Ratio: ","\n"
print "Overall: ", allpdb_nm/float(allpdb_total), "\n"
print >>cal,  "\n","Ratio: ","\n"
print >>cal, "Overall: ", allpdb_nm/float(allpdb_total), "\n"
for subpop in subpop2allpdb_match.keys():
	print subpop, ": ", subpop2allpdb_mismatch[subpop]/float(subpop2allpdb_total[subpop])
	print >>cal, subpop, ": ", subpop2allpdb_mismatch[subpop]/float(subpop2allpdb_total[subpop])

#print "\n","Ratio: ", "EUR: ", subpop2allpdb_mismatch["EUR"]/float(subpop2allpdb_total["EUR"])
cal.close()