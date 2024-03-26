import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import math
import pandas as pd
import sys

## This is a script created by Marta Alvarez Presas at the University of Bristol on March 2021 to read the output of the 
## script countingMSAsize.py from E. Ocaña-Pallares, and then generate a normal distribution, selecting the 25 percentile.
## If you use it, please, cite.

print('''Usage: python align_filter_normal_distrib.py Alignmentsizes_file.txt ''')


aln = pd.read_csv(sys.argv[1], sep="\t", header=None) 
aln.columns = ['OG', 'length']
mean = aln.mean()
std = aln.std()

mu = mean
variance = std
sigma = math.sqrt(variance)
x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
plt.plot(x, stats.norm.pdf(x, mu, sigma))
plt.show()

size = aln["length"]
percentile = np.percentile(aln["length"], 25)
selection = aln[size > percentile]

filtered_alignments = selection["OG"].to_list()
selection["OG"].to_csv("filt_ali.csv", index=False)

discard = aln[size < percentile]
discard["OG"].to_csv("discard_ali.csv", index=False)
