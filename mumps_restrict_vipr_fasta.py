from Bio import SeqIO
with open("data/mumps_vipr_full.fasta", 'rU') as fh:
    seqs = [x for x in SeqIO.parse(fh, 'fasta') if 'MuV' in x.name]
for seq in seqs:
    seq.name = seq.name.replace('[','(').replace(']',')')
    seq.description = seq.name
    seq.id = seq.name
with open("data/mumps_vipr.fasta", 'w') as fh:
    SeqIO.write(seqs, fh, "fasta")
