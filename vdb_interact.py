import os, argparse
from vdb_upload import vdb_upload
from vdb_download import vdb_download

parser = argparse.ArgumentParser(description='Interact with the virus data base (VDB) to upload and download viruses')
parser.add_argument('-u', '--upload_file', type = str, default = None, help='fasta file to be uploaded to vdb')
parser.add_argument('-v', '--output_file', type = str, default = None, help='output file to output vdb result')
parser.add_argument('-r', '--raxml_time_limit', type = float, default = 1.0, help='number of hours raxml is run')
parser.add_argument('--interval', nargs = '+', type = float, default = None, help='interval from which to pull sequences')
parser.add_argument('--path', type = str, default = 'data/', help='path of file dumps')
parser.add_argument('--prefix', type = str, default = '', help='prefix of file dumps including auspice')
parser.add_argument('--test', default = False, action="store_true",  help ="don't run the pipeline")
parser.add_argument('--start', default = 'filter', type = str,  help ="start pipeline at specified step")
parser.add_argument('--stop', default = 'export', type=str,  help ="run to end")
parser.add_argument('--skip', nargs='+', type = str,  help ="analysis steps to skip")
parser.add_argument('--ATG', action="store_true", default=False, help ="include full HA sequence starting at ATG")
parser.add_argument('--resolution', type = str,  help ="label for the resolution")
parser.add_argument('--pdb', default = ['1HA0', '2YP7'] type = list,  help ="pdb_structures to run foldx on")


class vdb_interact(object):
    '''
    Interact with vdb user to upload viruses, query vdb or download all sequences
    - option to select uploading, or downloading sequences
    - suboptions within uploading
        - default fasta format
        - specify special fasta format
    - suboptions within downloading
        - specify file to upload, default is ...
        - download all viruses from collection
            - specify what makes best sequence, default is "longest"
        - look for one virus by attribute
        - query for viruses that meet specific attributes
        - count how many viruses meet specific attributes

    '''

