import numpy as np, sys, matplotlib.pyplot as plt, pandas as pd, seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages

class Heat:

	def __init__(self, sind, aind, bind, aucind, xvar):
		self.sind = sind
		self.aind = aind
		self.bind = bind
		self.aucind = aucind
		self.xvar = xvar
		return

	def fillval( self, lines ):

		#print (sigmaval)

		aucval = dict() #np.zeros( len(aval) + len(bval) )

		for line in lines:

			if line[:2] == 'ma':
				line = line.replace('\n', ''); line = line.replace('[', ''); line = line.replace(']', ''); line = line.replace(',', '')
				line = line.split()
				k = line[self.sind] + ' ' + line[self.aind] + ' ' + line[self.bind]

				aucval[ k ] = line[ self.aucind ]


		return aucval


	def heatplot( self, aucval ):

		df = pd.DataFrame( columns=['sigma', self.xvar, 'auc'] )

		for k, v in aucval.items():
			k = k.split()
			k = [ float(i) for i in k ]

			if self.xvar == 'alpha':
				reqind = 1
			else:
				reqind = 2

			#if str(k[0]) == str(s) : #(k[0] - s ) < 1e-15:
			df.loc[ df.shape[0] ] = [ round( np.log(k[0])/np.log(10), 2), round( np.log(k[reqind])/np.log(10), 2), v ]
			print('k: %s   reqind: %s', k, reqind)

		#df = df.drop_duplicates( keep=False)

		#print ( df )

		df = df.drop_duplicates( subset = ['sigma', xvar], keep=False)

		print (df)

		df = df.pivot( 'sigma', self.xvar, 'auc' )

		df = df.fillna(0)

		df = df[df.columns].astype(float)  # or int

		fig = plt.figure()

		sns.heatmap( df, annot=True, vmin=0, vmax=1 )

		# plt.title('Beta: 0' )

		return fig




if __name__ == '__main__':

	filename = sys.argv[1]
	xvar = sys.argv[2]

	fid = open( filename, 'r')
	lines = fid.readlines()
	fid.close()

	heatobj = Heat( 2, 3, 4, 5, xvar)

	aucval = heatobj.fillval( lines )

	print (aucval)

	pp = PdfPages( filename + 'heatmap.pdf')
	
	fig = heatobj.heatplot( aucval )

	pp.savefig( fig )

	pp.close()

	sys.exit()

	# #for k, v in aucval.items():
	# #	print ( k, v )
    #
	# #sys.exit()
    #
	# sigma = np.logspace(-6, 2, 10)
	#
	# pp = PdfPages( filename + 'heatmap.pdf')
	#
	# for si, s in enumerate(sigma):
	# 	#print ( s, a, b)
	# 	#aucval = fillval( lines, s, 3, 4 )
	# 	#print ( len(aucval) )
	# 	#s = sigma[1]
	# 	#print ('sigma:', s )
	# 	#for k,v in aucval.items():
	# 	#	print (k, v )
	# 	#sys.exit()
	# 	#if ( si == 2):
	# 	#	sys.exit()
	#
	# 	# plot heatmap
	# 	fig = heatplot( s, aucval )
	# 	pp.savefig( fig )
	# 	#print ( 'Done' )
	# 	#sys.exit()
	#
	# pp.close()
	#sys.exit()
#filename = 'aucreport'
#
#x = np.logspace(-4,4,10)
#y = np.logspace(-4,4,10)
#
#x = np.log(x)/np.log(10)
#y = np.log(y)/np.log(10)
#
#intensity = []
#
#for line in open( filename, 'r'):
#	line = line.replace('\n', '')
#	line = line.split()
#	auc = line[2]
#	#print (auc)
#	intensity.append( float(auc) )
#
#reqintensity = [ intensity[10*i : 10*(i+1)] for i in range(10) ]
#
##setup the 2D grid with Numpy
#x, y = np.meshgrid(x, y)
#
##convert intensity (list of lists) to a numpy array for plotting
#reqintensity = np.array(reqintensity)
#
##print (len(x), len(y), x, y, reqintensity.shape)
##sys.exit()
#
##x = [1, 2, 3, 4, 5]
##y = [0.1, 0.2, 0.3, 0.4, 0.5]
##
##reqintensity = [
##    [5, 10, 15, 20, 25],
##    [30, 35, 40, 45, 50],
##    [55, 60, 65, 70, 75],
##    [80, 85, 90, 95, 100],
##    [105, 110, 115, 120, 125]
##]
#
#plt.figure()
#
##now just plug the data into pcolormesh, it's that easy!
#plt.pcolormesh(x, y, reqintensity)
#plt.colorbar() #need a colorbar to show the intensity scale
##plt.show() #boom
#plt.xlabel('alpha')
#plt.ylabel('beta')
#plt.title('Validation with 2/3 train and 1/3 validation set')
#
#plt.savefig('heatmap'+filename+'.png')



