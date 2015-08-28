import os,sys
import numpy.matlib
import numpy as np
from scipy import interpolate
import pyrap.measures as pm, pyrap.tables as pt, pyrap.quanta as qa
import optparse
import lofar.parmdb as pdb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.linalg as pl
import pylab

me = pm.measures()
parser = optparse.OptionParser()

required = optparse.OptionGroup(parser, "Required Attributes")

required.add_option('--ms', help='Measurement Set', action='store', type='string', dest='inputms')
required.add_option('--in', help='Instrument table', action='store', type='string', dest='inputin')

parser.add_option_group(required)

options, arguments = parser.parse_args()


t = pt.table(options.inputms, readonly=True, ack=False)


#Getting mean position of array

ant = pt.table(options.inputms+'/ANTENNA',readonly=False) #opening antenna table
antpos = ant.getcol('POSITION') #getting positions of all antennas 
meanpos = np.mean(antpos,0) #getting mean position of antennas
frame = ant.getcolkeyword('POSITION','MEASINFO')['Ref']
units = ant.getcolkeyword('POSITION','QuantumUnits')
mpos = me.position(frame,str(meanpos[0])+units[0],str(meanpos[1])+units[1],str(meanpos[2])+units[2])

me.doframe(mpos)

latr=me.measure(mpos,'WGS84')['m1']['value'] #get latitude in radians


#Getting Field direction

field = pt.table(options.inputms+'/FIELD') #opening field table
nfld=field.nrows()
dirs=field.getcol('DELAY_DIR')[:,0,:]
field.close()
print 'Found as many as '+str(nfld)+' fields.' #getting number of fields

#Getting Spectral windows

spec = pt.table(options.inputms+'/SPECTRAL_WINDOW')
nspec = spec.nrows()
bandnames=[x.split('#')[0].split('_')[-1] for x in spec.getcol('NAME')]
spec.close()
print 'Found as many as '+str(nspec)+' spws.' #getting bandnumbers


rah=dirs[0,0]*12.0/np.pi #getting right acesion
decr=dirs[0,1] # getting declination

R=np.zeros((nspec,nfld))
Q=np.zeros((nspec,nfld))


U=np.zeros((nspec,nfld))
mask=np.ones((nspec,nfld),dtype=bool)


QU={}

rah = dirs[0,0]*12.0/np.pi
dech = dirs[0,1]

table = pdb.parmdb(options.inputin)


times = table.getValuesGrid(table.getNames()[0])[table.getNames()[0]]['times'] #only need to get one row of times from the gain table

#time0=86400.0*floor(times[0]/86400.0) **** Why were these lines not included? TC ****
#rtimes=times-time0 

gainsXX = table.getValuesGrid(table.getNames()[0])[table.getNames()[0]]['values'] # **** 00 is XX and 11 is YY TC ****
gainsYY = table.getValuesGrid(table.getNames()[6])[table.getNames()[6]]['values']
freqs = table.getValuesGrid(table.getNames()[0])[table.getNames()[0]]['freqs']

# Getting the range of parangles
parang=np.zeros(len(times))


for itim in range(len(times)):
        tm=me.epoch('UTC',str(times[itim])+'s')
        last=me.measure(tm,'LAST')['m0']['value']
        last-=np.floor(last)  # days
        last*=24.0  # hours
        ha=last-rah  # hours
        har=ha*2.0*np.pi/24.0
        parang[itim]=np.arctan2( (np.cos(latr)*np.sin(har)),(np.sin(latr)*np.cos(decr)-np.cos(latr)*np.sin(decr)*np.cos(har)) )

#parang+=(paoffset*pi/180.) #unknown offset need to check
parangd = parang*(180.0/np.pi)

# Get the amplitude ratio
ampsXX=np.absolute(gainsXX) 
ampsXX[ampsXX==0.0]=1.0   # set any zero gains to 1
ampsYY=np.absolute(gainsYY)
ampsYY[ampsYY==0.0]=1.0
ratio=ampsXX/ampsYY  # get ratio of gains
#ratio.resize(nrows/nants,nants) *** only doing one station at the moment ****


#A=np.ones((nrows/nants,3))
A=np.ones((len(times),3)) # *** Only doing one station at hte moment ***
A[:,1]=np.cos(2*parang[:])
A[:,2]=np.sin(2*parang[:])

fit=pl.lstsq(A,np.square(ratio[:,0])) # (Gx'/Gy') = (Gx/Gy)*(1-2(Q/I)cos2PHI-2(U/I)sin2PHI)^0.5


model=fit[0][0]*A[:,0]+fit[0][1]*A[:,1]+fit[0][2]*A[:,2] #fit[0][0] is the (Gx/Gy)^2 ; fit[0][1] is 2*Q*(Gx/Gy)^2 ; fit[0][2] is 2*U*(Gx/Gy)^2

R[0,0]=fit[0][0]
Q[0,0]=fit[0][1]/R[0,0]/2.0
U[0,0]=fit[0][2]/R[0,0]/2.0
P=np.sqrt(Q[0,0]**2+U[0,0]**2)
X=0.5*np.arctan2(U[0,0],Q[0,0])*180/np.pi

"""ants0=range(nants)
#print ratio.shape
rsum=pl.sum(ratio[:,ants0],1)
#print rsum.shape
rsum/=len(ants0)

fit=pl.lstsq(A,np.square(rsum))

model=fit[0][0]*A[:,0]+fit[0][1]*A[:,1]+fit[0][2]*A[:,2]"""


pylab.subplot(111)
#colors = cm.rainbow(np.linspace(0, 1, len(parangd[:,0])))
rsumsq=np.square(ratio)
parangd_wrpd = parangd[:]
parangd_wrpd[np.where(parangd_wrpd<180.)] = parangd[:]
pylab.scatter(parangd_wrpd,rsumsq)
pylab.plot(parangd_wrpd,model,c='r',label="PKS1057-797")
pylab.plot([-150.,150.],[1.,1.],ls=':')
pylab.xlabel('Parallactic Angle [deg]')
pylab.ylabel(r'$g_{\rm x}/g_{\rm y}$')
pylab.legend()
#pylab.axis([-150.,150.,0.85,1.15])
pylab.show()

