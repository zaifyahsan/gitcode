import sys, math, random, os, matplotlib.pyplot as plt
import numpy as np

class Cell:
	def __init__(self, seq_1 ):
		self.sequence_1 = seq_1
		#self.structure_1 = struc_1


class Chamber:

	def __init__(self, population, kmerlist):
		self.initpop = population
		self.kmerlist = kmerlist
		return

	def getseq( self, pops ):
		seq = []

		for p in pops:
			seq.append( p.sequence_1 )

		return seq


	def populate( self, target, pop_size):
		population = []
		for i in range(pop_size):
			#get a random sequence to start with
			sequence = "".join([random.choice("ATCG") for _ in range(len(target))]);
			new_cell = Cell(sequence)
			new_cell.id = i
			new_cell.parent = i
			population.append(new_cell)
		return population

	def compute_fitness(self, population, target, beta, nofit):
		"""
			Assigns a fitness and bp_distance value to each cell in the population.
		If either of the k-mer is present, cell is fit otherwise not
		"""

		tot = []

		for cell in population:
			cell.fitness = 0.0
			for kmer in self.kmerlist:
				#kmer = 'GT'
				#print (kmer)
				if kmer in cell.sequence_1 or nofit:
					cell.fitness = cell.fitness + 1
			#print ( cell.sequence_1, cell.fitness)
			#sys.exit()

			tot.append( cell.fitness )

		norm = np.sum( tot )
		#print ( norm )

		#sys.exit()

		for cell in population:
			cell.fitness = cell.fitness / norm

		return None


	def mutate(self, sequence, mutation_rate):
		"""Takes a sequence and mutates bases with probability mutation_rate"""

		#start an empty string to store the mutated sequence
		new_sequence = ""
		#boolean storing whether or not the sequence got mutated
		mutated = False
		#go through every bp in the sequence
		for bp in sequence:
			#generate a random number between 0 and 1
			r = random.random()
			#if r is below mutation rate, introduce a mutation
			if r < mutation_rate:
				#add a randomly sampled nucleotide to the new sequence
				new_sequence = new_sequence + random.choice("ATCG")
				mutated = True
			else:
				#if the mutation condition did not get met, copy the current bp to the new sequence
				#print(bp)
				new_sequence = new_sequence + bp

		return (new_sequence, mutated)


	def selection( self, population, target, mutation_rate, beta, nofit):
		"""
			Returns a new population with offspring of the input population
			"""

		#select the sequences that will be 'parents' and contribute to the next generation
		#look at the documentation for np.random.choice and its optional argument p
		parents = np.random.choice(population, len(population), p=[cell.fitness for cell in population], replace=True)

		#build the next generation using the parents list
		next_generation = []
		for i, p in enumerate(parents):
			new_cell = Cell(p.sequence_1) #, p.structure_1)
			new_cell.id = i
			new_cell.parent = p.id

			next_generation.append(new_cell)

		#introduce mutations in next_generation sequeneces and re-fold when a mutation occurs
		for rna in next_generation:
			mutated_sequence_1, mutated_1 = self.mutate(rna.sequence_1, mutation_rate)
			#mutated_sequence_2, mutated_2 = mutate(rna.sequence_2)

			#if mutation occured assign and fold the new sequence
			if mutated_1:
				rna.sequence_1 = mutated_sequence_1
				#rna.structure_1 = computeStruc(mutated_sequence_1)
			#if mutated_2:
				#rna.sequence_2 = mutated_sequence_2
				#rna.structure_2 = nussinov(mutated_sequence_2)
			else:
				continue

		#update fitness values for the new generation
		self.compute_fitness(next_generation, target, beta, nofit)

		return next_generation


	def evolve( self, target, generations, pop_size, mutation_rate, beta, nofit):
		"""
			Takes target structure and sets up initial population, performs selection and iterates for desired generations.
			"""
		#store list of all populations throughout generations [[cells from generation 1], [cells from gen. 2]...]
		populations = []
		#start a dictionary that will hold some stats on the populations.
		population_stats = {}

		#get a starting population
		#for p in self.initpop:
		#	    print ( p.sequence_1 )

		initial_population = self.initpop #init_populate(target, pop_size)
		#compute fitness of initial population
		self.compute_fitness(initial_population, target, beta, nofit)

		#set current_generation to initial population.
		current_generation = initial_population


		#for p in current_generation:
		#    print ( p.sequence_1 )

		#sys.exit()


		#iterate the selection process over the desired number of generations
		for i in range(generations):

			#let's get some stats on the structures in the populations
			#record_stats(current_generation, population_stats)
			#print (i, '\n', current_generation )
			bpd = []
			#for p in current_generation:
			#    bpd.append( p.bp_distance_1 )

			#population_stats.setdefault('mean_bp_distance_1', []).append( np.mean( bpd ) )


			#add the current generation to our list of populations.
			populations.append(current_generation)

			#select the next generation
			new_gen = self.selection(current_generation, target, mutation_rate, beta, nofit)
			#set current generation to be the generation we just obtained.
			current_generation = new_gen

		return (populations, population_stats)

	def genexamples( self, no_examples, target, pop_size, nofit, tag, wfname ):

		wfid = open( wfname, 'w')

		for n in range(no_examples):

			# starting bubble of positive data
			pop_start = self.populate(target, pop_size)

			# generate 20 generations

			self.initpop = pop_start  # self.populate( target, pop_size)

			pops, pops_stats = self.evolve(target, generations, pop_size, mus[0], beta, nofit)

			# print ( creategen.getseq(pops[ len(pops) - 1]) )

			pop_start_lastgen = pops[len(pops) - 1]

			# print ( creategen.getseq(pop_start_lastgen) )

			# exit()

			# select one dna and replicate it 500 times
			pop_start_lastgen_dna = self.getseq(pop_start_lastgen)[0]

			#print(pop_start_lastgen_dna)

			pop_start_anc = []

			for i in range(pop_size):
				new_cell = Cell(pop_start_lastgen_dna)
				new_cell.id = i
				new_cell.parent = i
				pop_start_anc.append(new_cell)

			#print( tag )
			#print( self.getseq( pop_start_anc ))

			#self.initpop = self.populate( target, pop_size)

			genseq1 = pop_start_lastgen_dna #self.getseq( self.initpop )


			# generate 2 and 3 from 1
			#print ('init gen 1 ' )

			self.initpop = pop_start_anc

			pops, pops_stats = self.evolve( target, generations, pop_size, mus[0], beta, nofit )

			pops2 = pops[ len(pops) - 1 ]

			genseq2 = self.getseq( pops2 )

			#print (' init gen 1 ' )
			pops, pops_stats = self.evolve( target, generations, pop_size, mus[0], beta, nofit )
			pops3 = pops[ len(pops) - 1 ]
			genseq3 = self.getseq( pops3 )



			# generate 4 and 5 from 2
			self.initpop = pops2

			#print ( ' init gen 2 ' )
			pops, pops_stats = self.evolve(target, generations, pop_size, mus[0], beta, nofit)

			pops4 = pops[ len(pops) - 1 ]

			genseq4 = self.getseq( pops4 )

			#print (' init gen 2 ' )
			pops, pops_stats = self.evolve(target, generations, pop_size, mus[0], beta, nofit)
			pops5 = pops[ len(pops) - 1 ]
			genseq5 = self.getseq( pops5 )


			# generate 6 and 7 from 3
			self.initpop = pops3

			# generate 2 and 3 from 1
			#print ('init gen 3 ' )

			pops, pops_stats = self.evolve(target, generations, pop_size, mus[0], beta, nofit)

			pops6 = pops[ len(pops) - 1 ]

			genseq6 = self.getseq( pops6 )

			#print ('init gen 3 ' )
			pops, pops_stats = self.evolve(target, generations, pop_size, mus[0], beta, nofit)
			pops7 = pops[ len(pops) - 1 ]
			genseq7 = self.getseq( pops7 )

			wfid.write ( '1.' + tag + str(n) + ' '+ genseq1 + '\n')
			wfid.write ( '2.' + tag + str(n) + ' '+ genseq2[0] + '\n')
			wfid.write ( '3.' + tag + str(n) + ' '+ genseq3[0] + '\n')
			wfid.write ( '4.' + tag + str(n) + ' '+ genseq4[0] + '\n')
			wfid.write ( '5.' + tag + str(n) + ' '+ genseq5[0] + '\n')
			wfid.write ( '6.' + tag + str(n) + ' '+ genseq6[0] + '\n')
			wfid.write ( '7.' + tag + str(n) + ' '+ genseq7[0] + '\n')



if __name__ == '__main__':

	#if len(sys.argv) < 3:
	#	print ( 'python evolve.py initgen_sequence lastgen_sequence' )
	#	sys.exit()

	no_examples = 2000;

	kmerlist = ['CACACAAC',
				'AGAACCGA',
				 'CCTGGCTA',
				 'TATCCGAA',
				 'TCAGTTTA',
				 'CGAAGGCC',
				 'CTTCAAGA',
				 'AACGGACT',
				 'AGCTGGCT',
				 'ACGGCGAG']

	#print (kmerlist)

	#exit()

	# input settings

	target = np.zeros(200)

	mus = [ 0.005, 0.01, 0.02, 0.05, 0.1 ]

	generations = 20
	pop_size = 500
	beta = -2

	creategen = Chamber(0, kmerlist)

	no_examples = 1000

	# create positive examples

	#creategen.genexamples(no_examples, target, pop_size, False, 'pregion', 'pos')

	# create negative examples

	creategen.genexamples(no_examples, target, pop_size, True, 'nregion', 'neg')

	exit()

	# create positive data

	# starting bubble of positive data
	pop_start = creategen.populate(target, pop_size)

	# print( pop_start)
	# exit()

	# generate 20 generations

	creategen.initpop = pop_start #creategen.populate( target, pop_size)

	pops, pops_stats = creategen.evolve( target, generations, pop_size, mus[0], beta, False) #getseq( creategen.initpop )

	#print ( creategen.getseq(pops[ len(pops) - 1]) )

	pop_start_lastgen = pops[ len(pops) -1 ]

	#print ( creategen.getseq(pop_start_lastgen) )

	#exit()

	# select one dna and replicate it 500 times
	pop_start_lastgen_dna = creategen.getseq(pop_start_lastgen)[0]

	print( pop_start_lastgen_dna)

	pop_start_anc = []

	for i in range(500):
		new_cell = Cell(pop_start_lastgen_dna)
		new_cell.id = i
		new_cell.parent = i
		pop_start_anc.append(new_cell)

	# print ( creategen.getseq( pop_start_anc)[0:2])


	# generate 7 generations 2000 times

	# create negative data
	# starting bubble of negative data
	# generate 20 generations without any pressure
	# select one dna and replicate it 500 times
	# generate 7 generations 2000 times








#
#
##	for line in open('reqTen8mer', 'r'):
##		#print (line )
##		kmerlist.append( line.replace('\n', '') )
#
#	#print (kmerlist)
#	#sys.exit()
#
#
#
#
#
#	creategen = Chamber(0)
#
#
#	for n in range(no_examples):
#
#		creategen.initpop = creategen.populate( targets, pop_size)
#
#		genseq1 = creategen.getseq( creategen.initpop )
#
#
#		# generate 2 and 3 from 1
#		#print ('init gen 1 ' )
#
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#
#		pops2 = pops[ len(pops) - 1 ]
#		
#		genseq2 = creategen.getseq( pops2 )
#
#		#print (' init gen 1 ' )
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#		pops3 = pops[ len(pops) - 1 ]
#		genseq3 = creategen.getseq( pops3 )
#
#
#
#		# generate 4 and 5 from 2
#		creategen.initpop = pops2
#
#		#print ( ' init gen 2 ' )
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#
#		pops4 = pops[ len(pops) - 1 ]
#		
#		genseq4 = creategen.getseq( pops4 )
#
#		#print (' init gen 2 ' )
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#		pops5 = pops[ len(pops) - 1 ]
#		genseq5 = creategen.getseq( pops5 )
#
#
#		# generate 6 and 7 from 3
#		creategen.initpop = pops3
#
#		# generate 2 and 3 from 1
#		#print ('init gen 3 ' )
#
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#
#		pops6 = pops[ len(pops) - 1 ]
#		
#		genseq6 = creategen.getseq( pops6 )
#
#		#print ('init gen 3 ' )
#		pops, pops_stats = creategen.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#		pops7 = pops[ len(pops) - 1 ]
#		genseq7 = creategen.getseq( pops7 )
#
#		print ( '1.pregion' + str(n) + ' '+ genseq1[0] )
#		print ( '2.pregion' + str(n) + ' '+ genseq2[0] )
#		print ( '3.pregion' + str(n) + ' '+ genseq3[0] )
#		print ( '4.pregion' + str(n) + ' '+ genseq4[0] )
#		print ( '5.pregion' + str(n) + ' '+ genseq5[0] )
#		print ( '6.pregion' + str(n) + ' '+ genseq6[0] )
#		print ( '7.pregion' + str(n) + ' '+ genseq7[0] )


		#sys.exit()

#	init_filename = sys.argv[1]
#	lastgen_filename = sys.argv[2]
#
#	
#	chamber = Chamber( init_filename )
#
#
#	
#	pops, pops_stats = chamber.evolve( targets, generations, pop_size, mus[0], beta, kmerlist )
#
#	last_gen = pops[ len(pops) - 1 ]
#
#	#print (last_gen)
#
#	wfile = open( lastgen_filename, 'w' )
#
#	for p in last_gen:
#		print (p.sequence_1)
#		wfile.write( p.sequence_1 + '\n' )
#
#	wfile.close()


#	initial_population = populate( targets, pop_size )
#
#	for p in initial_population:
#		print (p.sequence_1)
#
#	sys.exit()

	#print ( pops )
	#print ( len(pops) )
	#for p in pops:
	#	print (p.sequence_1)
#	a = np.zeros(15)
#	
#	population = populate( a, 10000 )
#
#	#print ( population )
#
#	for p in population:
#		print (p.sequence_1)

	# read required kmer list

#
#
#colors = ['red', 'blue', 'green', 'yellow' ]
#
#generations = 500
#pop_size = 100
#
#beta = -2
#
#for trange in range(0, len(targets)):
#    
#
#    target = targets[ trange ]
#
#    plotname = 'target' + str(trange) + '.png'
#
#    plt.figure()
#    
#    for murange in range( len(mus) ):
#        
#        mu = mus[ murange ]
#        
#        print ( 'target: ', trange, ' mu: ', mu )
#        
#        mutation_rate = mu
#
#        pops, pops_stats = evolve( target , generations, pop_size, mutation_rate, beta)
#        
#        pops_stats = pops_stats['mean_bp_distance_1']
#        #print ( len(pops_stats))
#
#        plt.scatter( range(0,generations), pops_stats, color=colors[murange], label = 'mu = '+str(mu) )
#
#        #break
#
#    plt.legend( loc='upper right')
#
#    plt.xlabel('Generations')
#    plt.ylabel('Average Distance')
#
#    plt.title( 'target = ' + target )
#
#    plt.savefig( plotname )
