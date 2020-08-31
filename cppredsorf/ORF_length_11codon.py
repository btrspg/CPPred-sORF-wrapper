import sys
try:
	from string import maketrans
except ImportError:
	maketrans=str.maketrans
from . import ORF

def extract_feature_from_seq(seq,stt,stp):
	'''extract features of sequence from fasta entry'''
	
	stt_coden = stt.strip().split(',')
	stp_coden = stp.strip().split(',')
	transtab = maketrans("ACGTNX","TGCANX")
	mRNA_seq = seq.upper()
	mRNA_size = len(seq)
	tmp = ORF.ExtractORF(mRNA_seq)
	(CDS_size1, CDS_integrity, CDS_seq1) = tmp.longest_ORF(start=stt_coden, stop=stp_coden)
	return (mRNA_size, CDS_size1,CDS_integrity,CDS_seq1)

start_codons = 'ATG'
stop_codons = 'TAG,TAA,TGA'
Coverage = 0

def len_cov(seq):

	mRNN_codons = ["TAC","AAC","TAT","ATC","TTC","GAG","AAG","GAT","GAC","AAT","GTG"]
	count_sum = 0
	(mRNA_size, CDS_size,CDS_integrity,CDS_seq) = extract_feature_from_seq(seq = seq, stt = start_codons,stp = stop_codons)
	mRNA_len = mRNA_size
	CDS_len = CDS_size
	Coverage = float(CDS_len)/mRNA_len
	Integrity = CDS_integrity
	CDS_sequence = CDS_seq
	for codon in mRNN_codons:
		count_sum += CDS_sequence.count(codon)
	return(CDS_len,Coverage,Integrity,count_sum)
