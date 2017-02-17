import collections
import sys
import re

if len(sys.argv) < 2:
	print ("python seqTofeat.py seqFilename")
	sys.exit()


	
allseq=dict()

def revcomp(s):

	# find pos of A,C,G,T
	inda = [m.start() for m in re.finditer('A', s )];
	indc = [m.start() for m in re.finditer('C', s )];
	indg = [m.start() for m in re.finditer('G', s )];
	indt = [m.start() for m in re.finditer('T', s )];

	# convert s into char list
	rs = list(s)

	# replace A,C,G,T
	for i in inda:
		rs[i] = 'T'
	for i in indt:
		rs[i] = 'A'
	for i in indc:
		rs[i] = 'G'
	for i in indg:
		rs[i] = 'C'

	#convert rs to string, reverse and return
	rs = ''.join(rs); rs = rs[::-1]

	return rs

def produceAll( S, prefix, k):
	    global allseq
	    if k == 0:
		    allseq[prefix]=0
		    return ['']

	    for i in range(0,len(S)):
		    newprefix=prefix
		    newprefix=newprefix+S[i]
		    produceAll(S, newprefix, k-1)

#generate feature vector indices i.e. all possible k-mers
produceAll(['A','C','G','T'], '', 6)
allseq = collections.OrderedDict(sorted(allseq.items()))

fname=sys.argv[1]
sfile= open(fname, 'r')
ffile= open(fname + '.freq', 'w')

for line in sfile:
	#print line
	line = line.split()
	if len(line) < 2:
		continue
	spname = line[0]
	line = line[1]
	line = line.replace('\n','')
	line = line.replace('(', '').replace(')', '')
	line = line.upper()
	line = line.replace('N', '').replace('U', '').replace('L', '')
	
	if len(line) >= 6:
		#print line
		
		#store word frequency
		for i in range(0, len(line)-5):
			#print line[i:i+6]
			word = line[i:i+6]
			allseq[ word ] = allseq[ word ] + 1
			rword = revcomp( word )
			allseq[ rword ] = allseq[ rword ] + 1

		# store species name
		ffile.write(spname); ffile.write(' ')

		# write word frequencies
		for k,v in allseq.items():
			ffile.write(str(v)); ffile.write(' ')
	
		#get the label
		if spname.find('pregion') == -1:
			label = '-1'
		else:
			label = '1'


		ffile.write(label)
		ffile.write('\n')

		#clear dictionary
		allseq = allseq.fromkeys(allseq, 0)


sfile.close()
ffile.close()


