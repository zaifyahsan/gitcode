
import sys
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_curve, auc
from sklearn import grid_search #GridSearchCV

class ML:

	def __init__(self):
		self.bestlasso = 0.0
		self.cores = 4
		self.kfold = 10

		return

	def trainLogistic( self, xtrain, ytrain, kfold ):

		parameters = {'C': np.logspace( -6, 6, 20 ) }

		logisticModel = LogisticRegression( penalty = 'l1' )

		if self.kfold > 1:
			clf = grid_search.GridSearchCV( logisticModel, parameters, cv = self.kfold , n_jobs = self.cores )
		else:
			clf = grid_search.GridSearchCV( logisticModel, parameters, n_jobs = self.cores )

		clf.fit( xtrain, ytrain )

		self.bestlasso = clf.best_params_['C']

		#print ( bestparam )

		logisticModel = LogisticRegression( penalty = 'l1', C = self.bestlasso ).fit( xtrain, ytrain)



		return logisticModel

	

if __name__ == '__main__':

	# take train and test data
	n = 10000
	x1 = np.random.normal( 0, 1, (n, 10) )
	x2 = np.random.normal( 3, 2, (n, 10) )

	x3 = np.random.normal( 1, 1, (n, 10) )
	x4 = np.random.normal( 2, 3, (n, 10) )

	#y1 = np.ones(100)
	#y2 = -1*np.ones( 100 )

	xtrain = np.vstack((x1, x2))
	xtest = np.vstack((x3, x4))

	ytrain = np.ones( 2 * n )  #np.vstack((y1, y2))
	#ytest = np.ones( 200 )
	 

	ytrain[ n:] = [ -1 for i in range(n) ]

	ytest = ytrain

	print ( xtrain.shape, ytrain.shape, xtest.shape, ytest.shape )

	#print ( ytrain )

	#sys.exit()

	# learn model on train data with 10 fold cross-validation
	lasso = ML()

	lasso.kfold = 10
	model = lasso.trainLogistic( xtrain, ytrain, 10)
	print ( 'best param', lasso.bestlasso )


	# report accuracy on test data

	ypred = model.predict( xtest )
	
	#fpr, tpr, thresholds = roc_curve(ytest, ypred )
	fpr, tpr, thresholds = roc_curve( ytest, ypred )
	roc_auc = auc(fpr, tpr)
	#print ("Area under the ROC curve : %f" % roc_auc)
	print ('AUC: ', roc_auc)



