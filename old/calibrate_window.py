import sys
from numpy import *

fname = sys.argv[1]
data = load(fname)

maxidx = len(data)
width = 150

weights = 1.-((arange(-width, width)/float(width))**2)

# The VCO frequencies are integers so we could dither them
# to avoid quantization error if we wanted to be fancy
# but it seems to make no differece
if False:
    from numpy.random import triangular
    data[:, 1] += triangular(-1.,0.,1.,size=len(data))

# Just fit the whole thing at once, to get a single coefficient
a, b = polyfit(data[:,0], data[:,1], 1)
print "%.1f %u" % (a,b)

# Slide through the data fitting PSL to IMC for data around each sample
coeffs = []
for idx in xrange(maxidx):
    idx1 = max(0, idx-width)
    idx2 = min(idx+width, maxidx)
    coeffs.append(polyfit(data[idx1:idx2,0], data[idx1:idx2,1], 1,
                            w=weights[idx1-idx+width:idx2-idx+width]))
coeffs = array(coeffs)


if False:
    import pylab

    xarr = arange(min(data[:,0]), max(data[:,0]))
    pylab.figure()
    pylab.scatter(data[:,0], data[:,1], c=arange(len(data)), edgecolor='none')
    pylab.plot(xarr, a*xarr+b)
    pylab.title('Overall fit')

    pylab.figure()
    pylab.plot(coeffs[:,0])
    pylab.axhline(a)
    pylab.title('IMC to PSL linear coefficient')

    pylab.figure()
    pylab.plot(coeffs[:,1])
    pylab.axhline(b)
    pylab.title('IMC to PSL offset')

    pylab.figure()
    pylab.plot(data[:,1]-a*data[:,0]-b,c='g')
    pylab.title('Residual')

    pylab.show()

