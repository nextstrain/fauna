## the name of the script means "Receiver of the Slain", it's one of the names of Odin: https://en.wikipedia.org/wiki/List_of_names_of_Odin


import re
from Bio import Entrez
from Bio import SeqIO
import datetime as dt

def convertDate(x,start,end):
    """ Converts calendar dates between given formats """
    return dt.datetime.strftime(dt.datetime.strptime(x,start),end)

#define email for entrez login
db           = "nuccore"
Entrez.email = "some_email@somedomain.com"

## define batch size for download
batchSize    = 100
retmax       = 10**9

## define list of accessions to be downloaded
accs=["LN559125","KT592520","KT592529","KT581415","KM015497","KM015504","KM015511","KT592536","JQ922311","KF425658","KF425665","KF425672","KM392474","KM392481","KM392488","KM392495","KM392502","KM392509"]

## apparently Entrez allows fetching via GIs rather than via accessions
query  = " ".join(accs)
handle = Entrez.esearch( db=db,term=query,retmax=retmax )
giList = Entrez.read(handle)['IdList']

#post NCBI query
search_handle     = Entrez.epost(db=db, id=",".join(giList))
search_results    = Entrez.read(search_handle)
webenv,query_key  = search_results["WebEnv"], search_results["QueryKey"] 


#fecth all results in batch of batchSize entries at once
for start in range(0,len(giList),batchSize):
    #fetch entries in batch
    handle = Entrez.efetch(db=db, rettype="gb", retstart=start, retmax=batchSize, webenv=webenv, query_key=query_key)
    
    ## IO the GenBank files
    records=SeqIO.parse(handle,format='gb')
    
    ## iterate over entries
    for record in records:
        ## record.id returns sequence version
        ## record.name returns accession
        
        
        ## fetch authors of associated reference
        authors=[ref.authors for ref in record.annotations["references"]]

        ## get features
        record_features=record.features

        ## the cds part could be replaced just by calling record.seq
        
        print '\n',record.name
        
        CDSs={}
        col_date='XXXX_XX_XX'
        
        ## iterate over sequence features
        for feat in record_features:
            
            qualifiers=feat.qualifiers
            qualifier_keys=qualifiers.keys()

            # coding sequence feature
            if feat.type=='CDS':
                ## name of coding sequence
                cds_name=feat.qualifiers['product'][0]
                ## extract coding sequence
                cds=feat.extract(record.seq)
                ## assign sequence to cds_name
                CDSs[cds_name]=cds
            
            ## source data
            if feat.type=='source':
                information=qualifiers

                ## choose which fields we're interested in
                wanted=['organism','country','strain','collection_date','host']

                ## check which fields are present in GenBank entry
                available=[gimme for gimme in wanted if information.has_key(gimme)]
                
                ## iterate through GenBank information fields
                for field in available:
                    ## for some reason the values of the qualifiers dict are lists with single values
                    entry_value=information[field][0]
                    
                    ## deal with collection date separately
                    if field=='collection_date':
                        ## check date precision, use appropriate conversion format
                        collection_date=entry_value
                        N_fields=len(collection_date.split('-'))

                        if N_fields==1:
                            col_date=convertDate(collection_date,'%Y','%Y_XX_XX')
                        elif N_fields==2:
                            col_date=convertDate(collection_date,'%b-%Y','%Y_%m_XX')
                        elif N_fields==3:
                            col_date=convertDate(collection_date,'%d-%b-%Y','%Y_%m_%d')

                    elif field=='strain':
                        strain_name=entry_value
                    elif field=='organism':
                        cerberus=re.search(' \(([A-Za-z\_\-\/0-9\.]+)\)',entry_value)
                        strain_name=cerberus.group(1)
                    ## non-collection-date field
                    else:
                        pass
#                         print field,information[field]

        ## strain name, collection date, coding sequences, CDS lengths
        print strain_name,col_date,CDSs.keys(),[len(CDSs[product]) for product in CDSs],len(record.seq)