---
layout: post
title: "Calculating Mean Conservation Score on Genomic Features"
modified:
excerpt: "Calculating feature conservation from the UCSC conservation track with Pandas"
tags: [code][bioinformatics][genomics][python][pandas]
---


---

This code will calculate the mean conservation score on the features in a BED file from the USCS Genome Browser.

It requires Python and Pandas (`pip install pandas`).

This also requires that your feature file be a BED-like format - chromosome, start, end, and name positions must follow the BED specification, but additional columns can be appended and preserved. Data from the UCSC Conservation Track (downloaded as single base conservation) should also be converted to one BED-like file per chromosome, which can be done trivially by using the `cat -n` command in Terminal, adding the chromosome name (in the same style as the feature file) to each line with `awk`, and ensuring that all whitespace characteres are tabs (`/t`).

A sample `featureFile` and `conservationFile` might look something like this:

{% gist dc7bc62a413cbbc32bbf944125845afd %}

Now, the code:

```python
from __future__ import division
import pandas as pd

FEATUREFILE = 'S2_STARRseq_rep1_vsControl_peaks.bed'
CONSERVATIONFILEDIR = './conservation/'

peakDF = pd.read_csv(str(FEATUREFILE), sep = '\t', header=None, names=['chrom','start','end','name','enrichmentVal'])
#Reject negative peak starts, if they exist (sometimes this can happen w/ MACS)
peakDF.drop(peakDF[peakDF.start <= 0].index, inplace=True)
peakDF.reset_index(inplace=True)
peakDF.drop('index', axis=1, inplace=True)
peakDF['conservation'] = 1.0

chromNames = peakDF.chrom.unique()

for chromosome in chromNames: 
	chromSubset = peakDF[peakDF.chrom == str(chromosome)]
	chromDF = pd.read_csv(str(CONSERVATIONFILEDIR) + str(chromosome)+'.bed', sep='\t', header=None, names=['chrom','start','end','conserveScore'])
	
	for i in xrange(0,len(chromSubset.index)):
		x = chromDF[chromDF.start >= chromSubset['start'][chromSubset.index[i]]]
		featureSubset = x[x.start < chromSubset['end'][chromSubset.index[i]]]
		x=None
		featureConservation = float(sum(featureSubset.conserveScore)/(chromSubset['end'][chromSubset.index[i]]-chromSubset['start'][chromSubset.index[i]]))
		peakDF.set_value(chromSubset.index[i],'conservation',featureConservation)
		featureSubset=None
peakDF.to_csv("featureConservation.td", sep = '\t')
```

This will write a tab-delimited, BED-like file of calculated mean conservation values on each feature in `featureFile`. It might look something like this:

{% gist 1749196f79a23687fe4fa6c6f025c0e8 %}


