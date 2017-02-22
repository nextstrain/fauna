import matplotlib.pyplot as plt
import seaborn as sns; sns.set(color_codes=True)
import numpy as np
import math

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--infile', default=None, type=str, help="file to graph")
parser.add_argument('--rtype', default="linear", type=str, help="type of regression: linear or lowess")

def x_y(fname, rtype):
        o = open(fname, 'r')

        x, y = [] , []
        hi, fra = [], []
        hi_v, fra_v = {}, {}

        for line in o.readlines():
            l = line.split('\t')
            tup = (l[0], l[1], l[2])
            if l[5] == 'hi':
                hi.append(tup)
                hi_v[tup] = float(l[4])
            elif l[5] == 'fra':
                fra.append(tup)
                fra_v[tup] = float(l[4])

        paired = []
        for t in hi:
            if t in fra:
                paired.append(t)

        ferret_ref_hi = {}
        ferret_ref_fra = {}
        x = []
        y = []
        for t in paired:
            if t[0] == t[1]:
                ferret_ref_hi[t[2]] = math.log(hi_v[t],2)
                ferret_ref_fra[t[2]] = math.log(fra_v[t],2)

                if rtype == "lowess":
                    x.append(0)
                    y.append(0)
        print("references:", len(x))
        count = 0
        for t in paired:
            if (t[0] != t[1]) and (t[2] in ferret_ref_hi.keys()):
                hi_drop = -( math.log(hi_v[t]) - ferret_ref_hi[t[2]])
                fra_drop = -( math.log(fra_v[t]) - ferret_ref_fra[t[2]])

                x.append(hi_drop)
                if hi_drop < 0:
                    count += 1
                y.append(fra_drop)
        print("total:", len(x))
        print("count:", count)
        return x, y


def plot_it(x,y,rtype):

    x = np.asarray(x)
    y = np.asarray(y)

    if rtype == "lowess":
        ax = sns.regplot(x=x, y=y, color="b", scatter_kws={'alpha':0.3}, x_jitter=.25, y_jitter=.25, lowess=True)
        ax.set(xlabel='HI Titer Drop', ylabel='FRA Titer Drop')
        ax.set_title('HI vs. FRA Titer Drops with Lowess Curve')

        plt.show()
    else:
        ax = sns.regplot(x=x, y=y, color="r", scatter_kws={'alpha':0.3}, x_jitter=.25, y_jitter=.1, fit_reg=True, ci=95)
        ax.set(xlabel='HI Titer Drop', ylabel='FRA Titer Drop')
        ax.set_title('HI vs. FRA Titer Drops with Linear Regression')

        regression = np.polyfit(x, y, 1)
        print(regression)

        plt.show()

if __name__=="__main__":
    args = parser.parse_args()
    f = args.infile
    r = args.rtype
    x, y = x_y(f, r)
    plot_it(x, y, r)
