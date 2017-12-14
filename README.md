## CSE 549 Project

Professor: [Rob Patro](http://www.robpatro.com)

Course webpage: [Computational Biology](https://rob-p.github.io/CSE549F17/)

### Title: Implementing Long Read Mapping Algorithms

Implement MinHash and Containment Hash for long reads. Long reads are the results of third-generation sequencing of genomes and transcriptomes which, unlike the short read sequencing, (i.e. second or “next”-generation sequencing) methodologies, generate a very long substrings of the reference, but have a higher error rate compared to short reads (~10-15% vs <1% for short reads) and hence should be treated using algorithms, data structures and approaches other than those commonly used in short read analysis. Specifically, approaches that use the idea of Min Hashing or related algorithms have recently drawn a lot of attention in the computational biology community, and have been shown to be useful in assembling and mapping long reads. In this approach, which is borrowed from NLP and estimating document similarity, we use a Min Hash to approximate overlap of two long reads or approximately map a long read to a reference. In addition to this, another approach named containment hashing has been proposed that as an advantage to Min Hashing is addressing potential length differences between the two strings being compared. The basic idea of Min Hashing is well explained in Figure 1 in [this](http://www.nature.com.proxy.library.stonybrook.edu/nbt/journal/v33/n6/pdf/nbt.3238.pdf) paper. 


The goal of this project is to implement the Min Hash and Containment Hash approaches in C++ (one of the implementations is available in Python in [MinHashMetagenomics](https://github.com/dkoslicki/MinHashMetagenomics) repository on github) and using the simulators [SimulatedBiologicalData](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/SimulatedBiologicalData.py) and/or [SimulatedBiologicalDataSmall](https://github.com/dkoslicki/MinHashMetagenomics/blob/master/src/SimulatedBiologicalDataSmall.py), compare the performance of these two approaches on the simulated data.
