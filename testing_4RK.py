import numpy as np
import matplotlib.pyplot as plt 
from amm_functions import statusbar
import time
import pickle

def freq_d(t):
	p = 2
	Bag1TimeConstant = 0.149
	Bag1Frequency = 60
	GammaDynamicGain = 50
	result = ((GammaDynamicGain**p)/(GammaDynamicGain**p + Bag1Frequency**p))*(1 - np.exp(-t/Bag1TimeConstant))
	return(result)
def beta(t):
	Bag1BetaNot = 0.0605
	Bag1Beta1 = 0.2592
	result = Bag1BetaNot + Bag1Beta1*freq_d(t)
	return(result)
def gamma(t):
	Bag1Gamma = 0.0289
	result = Bag1Gamma*freq_d(t)
	return(result)
def C(LengthFirstDeriv):
	result = 0.58*(LengthFirstDeriv >= 0) + 0.42 # when true C = 1 else, C = 0.42  
	return(result)
def f(t,y,T,Length,LengthFirstDeriv,LengthSecondDeriv):
	R = 0.46 #length dependency of the force-velocity relationship
	a = 0.3
	K_SR = 10.4649
	K_PR = 0.15
	M = 0.0002 # intrafusal fiber mass

	LN_SR = 0.0423

	L0_SR = 0.04 #in units of L0
	L0_PR = 0.76 #polar region rest length

	result = float((K_SR/M)*((C(LengthFirstDeriv) * beta(t) * np.sign(LengthFirstDeriv-y/K_SR) \
								* ((abs(LengthFirstDeriv-y/K_SR))**a) * (Length - L0_SR - T/K_SR - R)) \
								+ K_PR * (Length - L0_SR - T/K_SR - L0_PR) \
								+ M*LengthSecondDeriv + gamma(t) - T))
	return(result)
def g(t,y,T):
	result = y
	return(result)
def k_and_l_values(t,y,T,Length,LengthFirstDeriv,LengthSecondDeriv,dt):
	k0 = g(t,y[-1],T[-1])*dt
	l0 = f(t,y[-1],T[-1],Length,LengthFirstDeriv,LengthSecondDeriv)*dt
	k1 = g(t+dt/2,y[-1]+l0/2,T[-1]+k0/2)*dt
	l1 = f(t+dt/2,y[-1]+l0/2,T[-1]+k0/2,Length,LengthFirstDeriv,LengthSecondDeriv)*dt
	k2 = g(t+dt/2,y[-1]+l1/2,T[-1]+k1/2)*dt
	l2 = f(t+dt/2,y[-1]+l1/2,T[-1]+k1/2,Length,LengthFirstDeriv,LengthSecondDeriv)*dt
	k3 = g(t+dt,y[-1]+l2,T[-1]+k2)*dt
	l3 = f(t+dt,y[-1]+l2,T[-1]+k2,Length,LengthFirstDeriv,LengthSecondDeriv)*dt
	return([k0,k1,k2,k3],[l0,l1,l2,l3])
def fourth_order_runge_kutta(t,y,T,Length,LengthFirstDeriv,LengthSecondDeriv,dt):
	"""
	t must be a scalar of value at previous timestep. y and T must be lists with 
	initial values (i.e. most recent values at index [-1])
	"""
	k,l = k_and_l_values(t,y,T,Length,LengthFirstDeriv,LengthSecondDeriv,dt)
	y.append(y[-1] + (l[0] + 2*(l[1] + l[2]) + l[3])/6)
	T.append(T[-1] + (k[0] + 2*(k[1] + k[2]) + k[3])/6)

y = [0]
T = [0]
dt = (1e-6)*5/3
t = np.arange(0,15+dt,dt)
Length = 1.0766863061118987
LengthFirstDeriv = 0
LengthSecondDeriv = 0
StartTime = time.time()
for i in range(len(t)):
	fourth_order_runge_kutta(t[i],y,T,Length,LengthFirstDeriv,LengthSecondDeriv,dt)
	statusbar(i,len(t),StartTime=StartTime,Title='4th Order Runge Kutta')

print('\n')
def bag1_model(CE,Bag1,SamplingPeriod):
	import numpy as np
	Length,LengthFirstDeriv,LengthSecondDeriv = CE['Length'][-1],CE['Velocity'][-1],CE['Acceleration'][-1]
	DynamicSpindleFrequency, GammaDynamicGain = Bag1['DynamicSpindleFrequency'][-1], Bag1['GammaDynamicGain']
	Bag1TensionFirstDeriv, Bag1Tension = Bag1['TensionFirstDeriv'][-1], Bag1['Tension'][-1]
	## Feedback system parameters
	## Spindle Model
	p = 2

	Bag1TimeConstant = 0.149
	Bag1Frequency = 60

	Bag1BetaNot = 0.0605
	Bag1Beta = 0.2592
	Bag1Gamma = 0.0289

	G = 20000  #7000

	C = 0.58*(LengthFirstDeriv >= 0) + 0.42 # when true C = 1 else, C = 0.42       

	DynamicSpindleFrequencyFirstDeriv = (GammaDynamicGain**p/(GammaDynamicGain**p+Bag1Frequency**p) - DynamicSpindleFrequency)\
									/ Bag1TimeConstant
	DynamicSpindleFrequency = (SamplingPeriod)*DynamicSpindleFrequencyFirstDeriv + DynamicSpindleFrequency

	Bag1Beta = Bag1BetaNot + Bag1Beta * DynamicSpindleFrequency
	Bag1Gamma = Bag1Gamma * DynamicSpindleFrequency

	R = 0.46 #length dependency of the force-velocity relationship
	a = 0.3
	K_SR = 10.4649
	K_PR = 0.15
	M = 0.0002 # intrafusal fiber mass

	LN_SR = 0.0423

	L0_SR = 0.04 #in units of L0
	L0_PR = 0.76 #polar region rest length

	Bag1TensionSecondDeriv = float((K_SR/M)*(	(C*Bag1Beta*np.sign(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR) \
										* ((abs(LengthFirstDeriv-Bag1TensionFirstDeriv/K_SR))**a)  \
										* (Length-L0_SR-Bag1Tension/K_SR-R)) \
										+ K_PR*(Length-L0_SR-Bag1Tension/K_SR-L0_PR) \
										+ M*LengthSecondDeriv \
										+ Bag1Gamma \
										- Bag1Tension 	))
	Bag1['TensionFirstDeriv'].append(float(Bag1TensionSecondDeriv*(SamplingPeriod) + Bag1TensionFirstDeriv))
	Bag1['Tension'].append(float(Bag1['TensionFirstDeriv'][-1]*(SamplingPeriod) + Bag1Tension))
	Bag1['DynamicSpindleFrequency'].append(DynamicSpindleFrequency)

CE = {'Length':[Length],'Velocity':[0],'Acceleration':[0]}
Bag1 = {'GammaDynamicGain':50,'DynamicSpindleFrequency':[0],'Tension':[0],'TensionFirstDeriv':[0]}
SamplingPeriod = 1e-6
Time = np.arange(0,15+SamplingPeriod,SamplingPeriod)
StartTime = time.time()
for i in range(len(Time)):
	bag1_model(CE,Bag1,SamplingPeriod)
	statusbar(i,len(Time),StartTime=StartTime,Title='Forward Euler')
Bag1Tension = Bag1['Tension'][:-1]
# Time,Bag1Tension = pickle.load(open('Bag1Tension_FE.pkl','rb'))

plt.plot(Time,Bag1Tension,'r--')
plt.plot(t,T[:-1],'b')
plt.show()
