from update import update
from zika_upload import zika_upload
from upload import parser

class zika_update(update):
    def __init__(self, **kwargs):
        kwargs['virus'] = 'zika'
        update.__init__(self, **kwargs)
        self.zika_upload = zika_upload(**kwargs)

    def fix_name(self, name):
        return self.zika_upload.fix_name(name)

if __name__=="__main__":
    args = parser.parse_args()
    connVDB = zika_update(**args.__dict__)
    connVDB.update(**args.__dict__)
