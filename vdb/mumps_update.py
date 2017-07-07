from update import update
from mumps_upload import mumps_upload
from update import parser

class mumps_update(update, mumps_upload):
    def __init__(self, **kwargs):
        update.__init__(self, **kwargs)
        mumps_upload.__init__(self, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = mumps_update(**args.__dict__)
    connVDB.update(**args.__dict__)
