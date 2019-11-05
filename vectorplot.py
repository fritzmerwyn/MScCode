# # Plotten von Daten und exakter Funktion der radialen SE
# #==========================================================
# -- Library Imports -- #
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
# from sklearn import preprocessing
# -- File Import -- #
data1=np.genfromtxt('Data/vectordata.dat')
# data1=np.genfromtxt('Data/.dat')
# data2=np.genfromtxt('gluonmass_dressingB.dat', delimiter=' ')
# data3=np.genfromtxt('gluonmass_dressingB78.dat', delimiter=' ')
# data4=np.genfromtxt('gluonmass22m.dat', delimiter=' ')

# data2=np.genfromtxt('vector2_100.dat', delimiter=' ')
# data3=np.genfromtxt('vector3.dat', delimiter=' ')
# -- Processing Data -- #
Y1 = data1 # np.sign changes the sign according to wheather or not the WF is up-side-down or not (due to normalization)
# Y2 = data1[:,2] # np.sign changes the sign according to wheather or not the WF is up-side-down or not (due to normalization)
# Y3 = data1[:,3] # np.sign changes the sign according to wheather or not the WF is up-side-down or not (due to normalization)
# Y4= data4[:,1] # np.sign changes the sign according to wheather or not the WF is up-side-down or not (due to normalization)

# X1= data1[:,0] #np.power((data1[:,0]/1000.0),-2.0) # data[:,1] only gives the second column of the eigvecs.dat file
# X2= data1[:,0] #np.power((data1[:,0]/1000.0),-2.0) # data[:,1] only gives the second column of the eigvecs.dat file
# X3= data1[:,0] #np.power((data2[:,0]/1000.0),-2.0) # data[:,1] only gives the second column of the eigvecs.dat file
# X3= np.power((data3[:,0]/1000.0),-2.0) # data[:,1] only gives the second column of the eigvecs.dat file
# X4= np.power((data4[:,0]/1000.0),-2.0) # data[:,1] only gives the second column of the eigvecs.dat file

# X3= np.sign(data3[1,1])*data3[:,1]
# x1 = np.arange(0., 100. , 100.0/len(X1)) # the division is for displaying the correct x-axis scale of the data
# y2 = np.arange(0., 100. , 100.0/len(X2)) # the division is for displaying the correct x-axis scale of the data
# y3 = np.arange(0., 100. , 100.0/len(X3)) # the division is for displaying the correct x-axis scale of the data
x = np.arange(0.00, 64, 1.0) # this one is just to display the Exact WF as fine as needed.

# plt.xlim(0.0,10.4)
# plt.ylim(-0.1,1.2)
# plt.yscale('logit')
# -- Plotting the exact wavfunctions -- #
# plt.plot(x,2*np.exp(-x),'k', label = '1s - Exact WF')
# a1 = ((6./12.)*((6.0**(0.5))**(-1.)))
# a2 = (-2.0*np.sqrt(0.6666666666666666))
# a3 = 1.0

# plt.plot(x, a1 * np.exp(x*(-0.5))*x,'k', label = 'Exact WF')
# plt.plot(x,a2 * ((-6 + x)*x)/(81.*np.exp(x/3.)),'k',label = '3d - Exact WF')
# plt.plot(x,a3 * (x*(80 - 20*x + x**2))/(256.*np.sqrt(15)*np.exp(x/4.)),'k',label = 'Exact WF')
# plt.plot(x, a1 * np.exp(x*(-0.5))*x,'k', label = '3d - Exact WF')

# plt.plot(x,2*((0.5)**(1.5))*(1-(0.5*x))*np.exp(-(0.5)*x),'k', label = 'Exact WF')
# plt.plot(x,2*((1.0/3.0)**(1.5))*(1-(x*2.0/3.0)+((x**2))*2.0/27.0)*np.exp(-x*(1.0/3.0)),'k',label = 'Exact WF')
# -- Plotting the calculated Finite Element Method Wavefunctions -- #
# plt.plot(X1, Y1,'orange', markersize=4, label=r'Dressing Function A')
# plt.plot(X2, Y2,'r', markersize=4, label=r'Dressing Function B')
plt.figure()
plt.plot(x, Y1,color = 'magenta',marker = '+', markersize=4, label=r'Mass Function M with mc=0.005 GeV')
# plt.plot(X4, Y4,'g', markersize=4, label=r'Effective Mass')


# plt.plot(y2, X2,'bx-.', markersize=4, label=r'4d FDM 1500 ELMTs')
# plt.plot(y3, X3,'gx-.', markersize=4, label=r'FEM 1000 Elmts')
# -- Additional Plot Adjustments -- #
plt.xlim(0.0001,64.4)
# plt.yscale('log')
# plt.xscale('log')
plt.legend()
plt.title('Mass Function M')
plt.xlabel(r'$p^2$')
plt.ylabel(r'$M(p^2)$')
plt.grid(True)

# plt.show()
# plt.figure()
# plt.plot(X1, Y1,'orange', marker = '+',markersize=4, label=r'Dressing Function A')
#
# plt.xscale('log')
# plt.legend()
# plt.title('Dressing Function A')
# plt.xlabel(r'$p^2$')
# plt.ylabel(r'$A(p^2)$')
# plt.grid(True)
#
# # plt.show()
# plt.figure()
# plt.plot(X2, Y2,'r', marker = '+', markersize=4, label=r'Dressing Function B')
#
# plt.xscale('log')
# plt.legend()
# plt.title('Dressing Function B')
# plt.xlabel(r'$p^2$')
# plt.ylabel(r'$B(p^2)$')
# plt.grid(True)
plt.show()
