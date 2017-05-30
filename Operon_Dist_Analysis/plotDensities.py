import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde as kd
import numpy as np 

#Accesses files that have the values for positive and negative controls
posCtrl = open('posCtrlDist.txt', 'r')
negCtrl = open('negCtrlDist.txt', 'r')
posList = []
negList = []

#Adds positive and negative control values to a list
for line in posCtrl:
	posList.append(float(line.strip()))
for line in negCtrl:
	negList.append(float(line.strip()))

#Assembles Plots
fig, ax = plt.subplots(nrows = 1, ncols = 3)
a0, a1, a2 = ax.flatten()

#Plots the first histogram showing positive and negative controls
a0.hist( [posList, negList], 100,  histtype = 'bar', normed = True, label = ['PoitiveControl','NegativeControl'])
a0.set_title( 'Frequency of Operon Distances' )
a0.set_xlim([-75,500])
a0.set_xlabel('Distance(bp.)')
a0.set_ylabel('Frequency')
a0.legend()

#Gets gaussian dist functions from data
pos_kd = kd(posList)
neg_kd = kd(negList)
#Returnes a set of evenly spaced numbers
x_vals = np.linspace(-50, 8000, 12000)

#Analyzing covariance
pos_kd.covariance_factor = lambda : .1
pos_kd._compute_covariance()
neg_kd.covariance_factor = lambda : .1
neg_kd._compute_covariance()

#Adds positive and negative  control plots
a1.plot(x_vals, pos_kd(x_vals),  label = 'Positive Control')
a1.plot(x_vals, neg_kd(x_vals), color = 'green', label = 'Negative Control')
a1.set_xlim([-50, 500])
a1.set_title('Kernal Densities')
a1.set_xlabel('Distance(bp.)')
a1.set_ylabel('Probability')
a1.legend()

print 'Here 1'
probList = []
for i in range(-50,700):
	prb = (pos_kd(i)*.6)/((pos_kd(i)*.6)+(neg_kd(i)*.4))
	probList.append(prb)
a2.plot(range(-50, 700), probList, color = 'black', label = 'Posterior_Probability')
a2.set_title('Posterior Probability')
a2.set_ylabel("Probability")
a2.set_xlabel('Distance(bp.)')
a2.legend()

plt.show()
