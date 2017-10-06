import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--files', nargs='*', default=[], help="tsvs that will be concatenated")
parser.add_argument('-o', '--output', type=str, default="data/titers_complete.tsv")

def concat(files,out):
    with open(out, 'w') as o:
        for filename in files:

            print "Concatenating and annotating %s into %s." % (filename, out)

            if "cdc" in filename.lower():
                source = "cdc"
            elif "crick" in filename.lower():
                source = "crick"
            else:
                source = "none"

            if "egg" in filename.lower():
                passage = "egg"
            elif "cell" in filename.lower():
                passage = "egg"
            else:
                passage = "none"

            with open(filename, 'r') as f:
                for line in f.readlines():
                    print line
                    line = line.strip()
                    l = "%s\t%s\t%s\n" % (line, source, passage)
                    o.write(l)

if __name__=="__main__":
    args = parser.parse_args()
    concat(args.files, args.output)
