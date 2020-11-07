# Algorithms for DNA Seq

import matplotlib.pyplot as plt
import matplotlib

## Fasta File

def readGenome(filename):
	genome = ''
	with open(filename,'r') as fp:
		for line in fp:
			if not line[0] == '>':
				genome += line.rstrip()
	fp.close()
	return genome



genome = readGenome('lambda_virus.fa')

len(genome)

counts = {'A':0,'T':0,'C':0,'G':0}

for base in genome:
	counts[base] += 1

## Fastq Files


def readFastq(filename):
	sequences = []
	qualities = []
	with open(filename,'r') as fp:
		while True:
			fp.readline()
			seq = fp.readline().rstrip()
			fp.readline()
			qual = fp.readline().rstrip()
			if len(seq) == 0:
				break
			sequences.append(seq)
			qualities.append(qual)
	return sequences, qualities

sequences, qualities = readFastq('SRR835775_1.first1000.fastq')


### Convert phred to qual values

def phred33toQ(qual):
	return ord(qual) - 33


phred33toQ('#')


### Create and show a histogram

def createHist(qualities):
	hist = [0] * 50
	for qual in qualities:
		for phred in qual:
			q = phred33toQ(phred)
			hist[q] += 1
	return hist

h = createHist(qualities)

plt.bar(range(len(h)),h)

plt.show()



### Plot GC content by position

def findGCByPos(sequences):
	gc = [0] * 100
	totals = [0] * 100
	for seq in sequences:
		for i in range(len(seq)):
			if seq[i] == 'C' or seq[i] == 'G':
				gc[i] += 1
			totals[i] += 1
	for i in range(len(gc)):
		if totals[i] > 0:
			gc[i] /= float(totals[i])
	return gc


gc = findGCByPos(sequences)
plt.plot(range(len(gc)),gc)
plt.show()



### count gc content
import collections
count = collections.Counter()
for seq in sequences:
	count.update(seq)


print(count)

## Alignment Algorithms

### Practice function to ensure randint is actually generating random integer numbers
def randGen(times):
	dic = [0] * 20
	for _ in range(times):
		num = random.randint(1,20)
		dic[num - 1] += 1
	return dic


### Naive Exact Matching Algorithm
text = 'Welcome to this great Course on Next Generation Sequencing Data Analysis'
word = 'Course'

def naive(p,t):
	occurences = []
	for i in range(len(t) - len(p) + 1):
		match = True
		for j in range(len(p)):
			if t[i+j] != p[j]:
				match = False
				break
		if match:
			occurences.append(i)
	return occurences


genome = readGenome('phix.fa')
string = 'AGCTTAGATAGC'
oc = naive(genome, string)

### Generate Random Reads

import random
def generateReads(genome, numReads, readLen):
	reads = []
	for _ in range(numReads):
		start = random.randint(0, len(genome)-readLen)
		reads.append(genome[start : start+readLen])
	return reads

### naive matches test
reads = generateReads(genome, 100, 100)

numMatched = 0
for r in reads:
	matches = naive(r, genome)
	if len(matches) > 0:
		numMatched += 1

print('%d / %d reads matches exactly!' % (numMatched, len(reads)))


### Implement on actual fastq file

phix_reads, _ = readFastq('ERR266411_1.first1000.fastq')

#### Full match
numMatched = 0
n = 0
for r in phix_reads:
	matches = naive(r, ref)
	n += 1
	if len(matches) > 0:
		numMatched += 1

print('%d / %d matched the genome exactly!' % (numMatched, n))


#### Partial match 30

numMatched = 0
n = 0
for r in phix_reads:
	r = r[:30]
    matches = naive(r, ref)
    n += 1
    if len(matches) > 0:
        numMatched += 1

print('%d / %d matched the genome exactly!' % (numMatched, n))


#### Complement matches

def reverseComplement(seq):
	complement = {'A':'T', 'T':'A','G':'C','C':'G','N':'N'}
	t = ''
	for base in seq:
		t = complement[base] + t
	return t

numMatched = 0
n = 0
for r in phix_reads:
    r = r[:30]
    matches = naive(r, ref)
	matches.extend(naive(reverseComplement(r), ref))	### Note the extend function here
    n += 1
    if len(matches) > 0:
        numMatched += 1

print('%d / %d matched the genome exactly!' % (numMatched, n))


### Allow 2 mismatches
def naive_2mm(p,t):
	occurences = []
	for i in range(len(t) - len(p) + 1):
		match = True
		n = 0
		for j in range(len(p)):
			if t[i+j] != p[j]:
				n += 1
				if n>2:
					match = False
					break
		if match:
			occurences.append(i)
	return occurences


























