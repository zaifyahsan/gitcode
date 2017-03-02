import sys, pickle


orthofilename = sys.argv[1]

spregion = []; xregion = []

xtrain = []; ytrain = []; x_species = []


currdict = dict()
xorthodict = dict()
prev_region = ''
key = 0

for line in open( orthofilename, 'r' ):
	#print line
	line = line.replace('\n', '')

	line = line.split()

	lineid = line[0].split('.')
	spname = lineid[ 0 ].replace('_','')
	region = lineid[0]

	xortho = [ float(i) for i in line[1:4097] ]

	#print spname 

	if not spname == '7':
		spregion.append( spname )
		xregion.append( xortho )
	else:
		#print spname
		xtrain.append( xortho )
		ytrain.append( float( line[4097] ) )
		x_species.append( spname )

	#sys.exit()


	if not region == prev_region:

		#print region, spname, len(xortho), len(spregion), len(xregion)
		#print len(xtrain), len( ytrain )

		prev_region = region

		currdict['species'] = spregion
		currdict['X'] = xregion
		xorthodict[key] = currdict

		spregion = []
		xregion = []
		key = key + 1

	#sys.exit()

print len(xorthodict), len(xtrain), len(x_species), len(ytrain)
pickle.dump( xorthodict, open(orthofilename+'xorthodict.p', 'wb') 	)
pickle.dump( xtrain, open( orthofilename + '.p', 'wb' ) )
pickle.dump( x_species, open( orthofilename + 'species.p', 'wb') )
pickle.dump( ytrain, open( orthofilename + 'y.p', 'wb' ) )

