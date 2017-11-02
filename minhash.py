"""
An implementation of a MinHash bottom sketch, applied to k-mers in DNA.
"""
from __future__ import print_function
import khmer
import h5py
import numpy as np
import os
import re
from blist import *  # note, the import functions import the _mins etc. as lists, and the CE class imports them as blists.
# This shouldn't cause an issue, but will lead to slow performance if a CE is imported, then additional things are added.
# I.e. If you import a CE, don't add new elements, or you might have a bad day (or at least a long one).
import bisect
import warnings
notACTG = re.compile('[^ACTG]')


class CountEstimator(object):
    """
    A simple bottom n-sketch MinHash implementation.
    n is the number of sketches to keep
    Still don't know what max_prime is...
    """

    def __init__(self, n=None, max_prime=9999999999971., ksize=None):
        if n is None:
            raise Exception
        if ksize is None:
            raise Exception

        if ksize % 2 == 0:
            raise Exception("Due to an issue with khmer, only odd ksizes are allowed")

        self.ksize = ksize

        # get a prime to use for hashing
        p = get_prime_lt_x(max_prime)
        self.p = p

        # initialize sketch to size n
        #self._mins = [float("inf")]*n
        self._mins = blist([p]*n)

        # initialize the corresponding counts
        self._counts = blist([0]*n)

        self._kmers = blist(['']*n)
        self.hash_list = None


        # Optional container for the true number of k-mers in the genome used to populate the sketch
        self._true_num_kmers = 0


    def add(self, kmer):
        """
        Add kmer into sketch, keeping sketch sorted, update counts accordingly
        """
        _mins = self._mins
        _counts = self._counts
        _kmers = self._kmers

        h = khmer.hash_murmur3(kmer)
            #h = hash(kmer)

        h = h % self.p
        if self.hash_list:  # If I only want to include hashes that occur in hash_list
            if h not in self.hash_list:  # If the kmer isn't in the hash_list, then break
                return

        if h >= _mins[-1]:
            return

        i = bisect.bisect_left(_mins, h)  # find index to insert h
        if _mins[i] == h:  # if h in mins, increment counts
            _counts[i] += 1
            return
        else:  # otherwise insert h, initialize counts to 1, and insert kmer if necessary
            _mins.insert(i, h)
            _mins.pop()
            _counts.insert(i, 1)
            _counts.pop()
            if _kmers:
                _kmers.insert(i, np.string_(kmer))
                _kmers.pop()
            return

        assert 0, "should never reach this"

    def add_sequence(self, seq, rev_comp=False):
        """
         Sanitize and add a sequence to the sketch.
        """
        # seq = seq.upper().replace('N', 'G')
        seq = notACTG.sub('G', seq.upper())  # more intelligent sanatization?
        print ("seq", seq)
        for kmer in kmers(seq, self.ksize):
            self.add(kmer)


def kmers(seq, ksize):
    """yield all k-mers of len ksize from seq.
    Returns an iterable object
    """
    for i in range(len(seq) - ksize + 1):
        yield seq[i:i+ksize]


# taken from khmer 2.0; original author Jason Pell.
def is_prime(number):
    """Check if a number is prime."""
    if number < 2:
        return False
    if number == 2:
        return True
    if number % 2 == 0:
        return False
    for _ in range(3, int(number ** 0.5) + 1, 2):
        if number % _ == 0:
            return False
    return True

def get_prime_lt_x(target):
    """Backward-find a prime smaller than (or equal to) target.

    Step backwards until a prime number (other than 2) has been
    found.

    Arguments: target -- the number to step backwards from
    """
    if target == 1:
        return 1

    i = int(target)
    if i % 2 == 0:
        i -= 1
    while i > 0:
        if is_prime(i):
            return i
        i -= 2

    if i <= 0:
        raise RuntimeError("unable to find a prime number < %d" % (target))

