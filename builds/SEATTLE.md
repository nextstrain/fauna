# Seattle Flu Pipeline Notes

## Download

    python2 vdb/download.py -db vdb -v seattle --fstem seattle

## Upload

### Fred Hutch sequences

Upload with:

    python2 vdb/seattle_upload.py -db vdb -v seattle --source sfs --authors "Seattle Flu Study Consortium" --url https://github.com/seattleflu --title "Seattle Flu Study" --fname seattle.fasta
