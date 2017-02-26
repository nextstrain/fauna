from update import update
from ebola_upload import ebola_upload
from update import parser

class ebola_update(update, ebola_upload):
    def __init__(self, **kwargs):
        update.__init__(self, **kwargs)
        ebola_upload.__init__(self, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = ebola_update(**args.__dict__)
    connVDB.update(**args.__dict__)
