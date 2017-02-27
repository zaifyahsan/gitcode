'''
Take hyperparameters as input

'''
import logging
import sys, pickle, numpy as np, random, math
from sklearn.metrics import roc_curve, auc

from sklearn.linear_model import LogisticRegression

sys.path.insert(0, '..')
from phyloreg import species, classifiers
from treestructure import getTreeDict
from multiprocessing import Pool


def createXmatrix( x, no_feat):

	X = np.zeros( (len(x), no_feat), dtype = float )

	for i in range(len(x) ):
		item = x[i]
		#print item
		for j in range( len(item) ):
			X[i,j] = item[j]
	return X

def makeLabelOneZero( y ):

	y += np.ones( len(y))
	y = (1/2.0) * y
	return y



def mapper( ranges ):

	# create model
	#model = classifiers.RidgeRegression( ranges[1], ranges[2], normalize_laplacian=True )

	model = classifiers.LogisticRegression( ranges['alpha'], ranges['beta'], opti_learning_rate = 1e-3 ) #, normalize_laplacian=True )

	print 'model created'

	# create distance matrix
	adj_matrix = species.ExponentialAdjacencyMatrixBuilder( ranges['sigma'] )
	species_list, A = adj_matrix( ranges['phylo'] )

	print 'distance matrix created '


	#make labels as zero and one
	ranges['ytrain'] = makeLabelOneZero( ranges['ytrain'])

	ranges['ytest'] = makeLabelOneZero( ranges['ytest'])

	# learn model
	model.fit( ranges['xtrain'], ranges['x_speciestrain'], ranges['ytrain'], ranges['xorthodicttrain'], A, species_list ) # self.reqxtrain, self.reqx_speciestrain, self.reqytrain, self.reqxorthodicttrain, A, species_list )

	print 'model learned'
    #
#	model = LogisticRegression( penalty= 'l2', C = 1e9 )
#	model.fit( ranges['xtrain'], ranges['ytrain'] )
	#test model
	ypred = model.predict( ranges['xtest'] ) #self.reqxvalid )

	print 'model tested'


	print max( ranges['ytrain'] ), ' ', min( ranges['ytrain'] ), ' ', max( ranges['ytest'] ), ' ', min( ranges['ytest'] )

	# return ranges and auc
	#compute auc score
	fpr, tpr, thresholds = roc_curve( ranges['ytest'], ypred )
	roc_auc = auc(fpr, tpr)

	print 'auc: ', roc_auc
	returnranges = [ ranges['sigma'], ranges['alpha'], ranges['beta'], roc_auc ]
	print 'mapper ranges: ', returnranges

	return returnranges



if __name__ == '__main__':

	logging.basicConfig( level=logging.DEBUG,
			format="%(asctime)s.%(msecs)d %(levelname)s %(module)s: %(message)s")
	#sys.stdout = open('threeparvalidaucreport', 'w')


	# read phylo
	# phylo = getTreeDict( 'modtree100' )
    #
	phylo = {'1': (None, None),
			'2': ('1', 0.005),
			'3': ('1', 0.005),
			'4': ('2', 0.005),
			'5': ('2', 0.005),
			'6': ('3', 0.005),
			'7': ('3', 0.005)}

	# print ( phylo )




	# read data
	xfile=sys.argv[4] + 'train.p'; yfile= sys.argv[4] + 'trainy.p'; spfile= sys.argv[4] + 'trainspecies.p';
	orthdict= sys.argv[4] + 'trainxorthodict.p'; xtestfile= sys.argv[4] + 'test.p'; ytestfile= sys.argv[4] + 'testy.p' #):

	# read xtrain
	xtrain = pickle.load( open( xfile, 'rb') )
	X = createXmatrix( xtrain, 4096 )

	# read ytrain
	ytrain = pickle.load( open( yfile, 'rb' ) )

	# make ytrain as 1 or 0
	#ytrain =

	# read xtest
	xtest = pickle.load( open( xtestfile, 'rb' ) )
	Xtest = createXmatrix( xtest, 4096 )

	# read ytest
	ytest = pickle.load( open( ytestfile, 'rb' ) )

	# read x_species
	x_species = pickle.load( open( spfile, 'rb') )

	# read orthologs
	xorthodict = pickle.load( open( orthdict, 'rb') )

	logging.debug('xtrain: %s ytrain: %s x_species: %s xorthodict: %s xtest: %s ytest: %s ',
				  len(xtrain), len(ytrain), len(x_species), len(xorthodict), len(xtest), len(ytest) )

	#print ( xorthodict )

	#exit()

	# create train and valid
	# divide train set into train and validation
	totalrange = X.shape[0]
	indexrange = list( range(totalrange) )

	#random.shuffle( indexrange )

	#print totalrange, 'length of indexrange ', len(indexrange)

	trainindex = indexrange #[0: int( (2.0/3) * totalrange ) ]
	#validindex = indexrange[ int( (2.0/3) * totalrange ) : ]

	#print len(trainindex), len(validindex)

	reqxtrain = X[ trainindex, : ]
	reqytrain = np.asarray( ytrain, dtype = float)[ trainindex ]

	#reqytrain = (1/2.0) * np.add( reqytrain, np.ones( len( reqytrain) ) )


	reqx_speciestrain = x_species[ 0 : len(trainindex) ] #np.asarray(trainindex, dtype=int) ]
	reqxorthodicttrain = dict()

	for k in range( len(trainindex) ):
		reqxorthodicttrain[k] = xorthodict[ trainindex[k] ]
		#reqx_speciestrain = x_species[ trainindex[k] ]
		#print k

	#print ( len( reqx_speciestrain) )
	#print ( reqx_speciestrain )

	reqxvalid = Xtest #X[ validindex, : ]
	reqyvalid = np.asarray(ytest, dtype = float) #np.asarray(ytrain, dtype = float)[ validindex ]

	#reqyvalid = (1/2.0) * np.add( reqyvalid, np.ones( len( reqyvalid) ) )

	sigma_range = float(sys.argv[1]) #np.logspace(-3, 1, 10)
	alpha_range = float(sys.argv[2]) #np.logspace(-6, 6, 10)
	beta_range =  float(sys.argv[3]) #np.logspace(-6, 6, 10)

	ranges = { 'sigma': sigma_range,
			   'alpha': alpha_range,
			   'beta' : beta_range,
			   'xtrain': reqxtrain,
		       'x_speciestrain': reqx_speciestrain,
			   'ytrain': reqytrain,
			   'xorthodicttrain': reqxorthodicttrain,
			   'xtest': reqxvalid,
			   'ytest': reqyvalid,
			   'phylo': phylo
			   }

	print 'xtrain: ',reqxtrain.shape, ' ytrain: ', reqytrain.shape, ' xtest: ', reqxvalid.shape, ' ytest: ', reqyvalid.shape

	print max( reqytrain), ' ', min( reqytrain), ' ', max(reqyvalid), ' ', min(reqyvalid)


	#exit()

	#print sigma_range, alpha_range, beta_range

	#sys.exit()

	#print ranges[0], ranges[1], ranges[2]

	reqranges = mapper( ranges )

	#print reqranges


#	for r in range( len(sigma_range) ):
#		for j in range( len(alpha_range) ):
#			for k in range( len(beta_range) ):
#				ranges.append( ( sigma_range[r], alpha_range[j], beta_range[k], reqxtrain, reqx_speciestrain, reqytrain, reqxorthodicttrain, reqxvalid, reqyvalid, phylo ) )
#
#	print len(ranges)
#	p = Pool( 2 )
#
#	results = p.map( mapper, ranges )

#	reportfilename = 'threeparvalidaucreport'
#
#	report = open ( reportfilename, 'w')
#
#	for r in results:
#		print r
#		r = [ str(ri) for ri in r ]
#		report.write( ' '.join(r) )
#		report.write(' \n')
#
#	report.close()

