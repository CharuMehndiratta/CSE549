import minhash as MH  # import the min hash package
from pybloom import BloomFilter  # import bloom filters for the larger set
import numpy as np

# Define variables
prime = 9999999999971  # taking hashes mod this prime
ksize = 11  # k-mer length
h = 10  # number of hashes to use
p = 0.01  # false positive rate for the bloom filter
len_small_string = 100
len_large_string = 1000

# Create data: small set called A, large set called B

# Random 100 letter string
small_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_small_string))  # small string to form the small set A

# Creating shingles of 11 length
size_A = len(set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)]))  # size of smaller set, used to convert containment index to Jaccard index

 #Reference Genome
large_string = ''.join(np.random.choice(['A', 'C', 'T', 'G'], len_large_string)) + small_string  # large string to form the larger set B

# Populate min hash sketch with smaller set
A_MH = MH.CountEstimator(n=h, max_prime=prime, ksize=ksize)
A_MH.add_sequence(small_string)  # create the min hash of the small string

# Create the bloom filter and populate with the larger set
B_filt = BloomFilter(capacity=1.15*len_large_string, error_rate=p)  # Initialize the bloom filter


size_B_est = 0  # used to count the number of k-mers in B, could do much more intelligently (like with HyperLogLog)

print "large_string", large_string
for i in range(len(large_string) - ksize + 1):
	kmer = large_string[i:i+ksize]
	print "kmer", kmer
	if kmer not in B_filt:
		size_B_est += 1
		B_filt.add(kmer)

# Use the k-mers in the sketch of A and test if they are in the bloom filter of B
int_est = 0  # intersection estimate
for kmer in A_MH._kmers:
	if kmer is not '':  # in case the set "A" was so small the Min Hash was not fully populated
		if kmer in B_filt:
			int_est += 1

int_est -= np.round(p*h)  # adjust for the false positive rate
containment_est = int_est / float(h)  # estimate of the containment index
jaccard_est = size_A * containment_est / (size_A + size_B_est - size_A * containment_est)

# calulate true jaccard for comparison
A = set([small_string[i:i+ksize] for i in range(len(small_string) - ksize + 1)])  # the smaller set
B = set([large_string[i:i+ksize] for i in range(len(large_string) - ksize + 1)])  # the larger set
size_A = len(A)  # number of k-mers in A v
size_B = len(B)  # number of k-mers in B
true_jaccard = len(A.intersection(B)) / float(len(A.union(B)))

print("Containment index estimate: %f" % containment_est)
print("Jaccard index estimate (via the containment approach): %f" % jaccard_est)
print("True Jaccard index: %f" % true_jaccard)
print("Relative error: %f" % (np.abs(jaccard_est-true_jaccard) / true_jaccard))