from update import update
from zika_upload import zika_upload
from update import parser

class zika_update(update, zika_upload):
    def __init__(self, **kwargs):
        update.__init__(self, **kwargs)
        zika_upload.__init__(self, **kwargs)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = zika_update(**args.__dict__)
    connVDB.update(**args.__dict__)
