from update import update
from coronavirus_upload import coronavirus_upload
from update import parser

class coronavirus_update(update, coronavirus_upload):
    def __init__(self, **kwargs):
        update.__init__(self, **kwargs)
        coronavirus_upload.__init__(self, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = coronavirus_update(**args.__dict__)
    connVDB.update(**args.__dict__)
