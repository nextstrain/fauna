import argparse

parser = argparse.ArgumentParser()
parser.add_argument('files', nargs='+', default=[], help="tsvs that will be concatenated")

def concat(files):
    for filename in files:

        if "cdc" in filename.lower():
            source = "cdc"
        elif "crick" in filename.lower():
            source = "crick"
        else:
            source = "unknown"

        if "egg" in filename.lower():
            passage = "egg"
        elif "cell" in filename.lower():
            passage = "cell"
        else:
            passage = "none"

        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                l = "%s\t%s\t%s" % (line, source, passage)
                print l

if __name__=="__main__":
    args = parser.parse_args()
    concat(args.files)
