## Upload BCCDC genomes

If you have a FASTA file and CSV metadata, this script will help (with minor modifications as needed)

`python3 scripts/mumps.csv-and-fasta-to-vipr-fasta.py data/input.mumps.raw.fasta data/input.mumps.csv data/input.mumps.vipr.fasta`

Upload to fauna

`python3 vdb/mumps_upload.py -db vdb -v mumps --source bccdc --locus genome --fname mumps.bc.fasta --authors "Gardy et al" --title "Unpublished"`
