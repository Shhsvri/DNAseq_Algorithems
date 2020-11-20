## Fasta File
def readGenome(filename):
    genome = ''
    with open(filename,'r') as fp:
        for line in fp:
            if not line[0] == '>':
                genome += line.rstrip()
    fp.close()
    return genome


def naive(p,t):
    c = 0
    occurences = []
    for i in range(len(t) - len(p) + 1):
        match = True
        for j in range(len(p)):
            c += 1
            if t[i+j] != p[j]:
                match = False
                break
        if match:
            occurences.append(i)
    return c, occurences

def naive_2mm(p,t):
    occurences = []
    for i in range(len(t) - len(p) + 1):
        c = 0
        match = True
        for j in range(len(p)):
            if t[i+j] != p[j]:
                c += 1
                if c > 2:
                    match = False
                    break
        if match:
            occurences.append(i)
    return occurences


# Boyer Moore
import bm_preproc as bm

def boyer_moore(p, p_bm, t):
    i = 0
    occurrences = []
    c = 0
    while i < len(t)-len(p)+1:
        c += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if not p[j] == t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return c, occurrences


# read the genome
genome = readGenome('chr1.GRCh38.excerpt.fasta')


# number of alignments
p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'

print('Naive tries:')
print(len(genome) - len(p) + 1)

print('naive character comparisons:')
print(naive(p, genome))

p_bm = bm.BoyerMoore(p)
print('Boyer Moore:')
print(boyer_moore(p,p_bm, genome))


from kmer_index import Index

index = Index(genome, 8)
print(index)

