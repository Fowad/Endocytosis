from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy as np
from math import log

fig = plt.figure()
ax = fig.gca(projection='3d')

data = np.genfromtxt('Force_Dist.txt')
#zero coupon maturity dates
x = data[:,0]
#tenor
y = data[:,1]
#rates
z = data[:,2]

##maturity dates chart axis
#uniquemat = np.unique(y)
#end = (np.max(uniquemat))/(uniquemat[0])
##the zc rate maturity axis is arranged in log space
#yi = np.logspace(1, end,len(uniquemat),True,uniquemat[0])

#tenor chart axis
xi = np.unique(x)

yi = np.unique(y)


X, Y = np.meshgrid(xi, yi)
Z = griddata(x, y, z, xi, yi, interp='linear')

# Plot rate surface
print ("Plotting Force Distribution ...")
fig = plt.figure()
fig.suptitle('Force Dist Plot',fontsize=20)
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1,alpha=0.3)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('F')

## Override tenor axis labels
#labels = [item.get_text() for item in ax.get_xticklabels()]
#labels[0] = '1M'
#labels[1] = '3M'
#labels[2] = '6M'
#labels[4] = '1Y'

#ax.set_xticklabels(labels)

#Plot 3D contour
zzlevels = np.linspace(Z.min(),Z.max(),num=4,endpoint=True)
xxlevels = np.linspace(X.min(),X.max(),num=4,endpoint=True)
yylevels = np.linspace(Y.min(),Y.max(),num=4,endpoint=True)

#cset = ax.contour(X, Y, Z, zzlevels, zdir='z',offset=Z.min(), cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, xxlevels, zdir='x',offset=X.min(), cmap=cm.coolwarm)
#cset = ax.contour(X, Y, Z, yylevels, zdir='y',offset=Y.max(), cmap=cm.coolwarm)

#plt.clabel(cset,fontsize=10, inline=1)

ax.set_zlim3d(np.min(Z), np.max(Z))

plt.show()

