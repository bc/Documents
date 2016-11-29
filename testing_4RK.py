import numpy as np
import matplotlib.pyplot as plt 
from amm_functions import statusbar,append_dictionary
import time
import pickle


def fourth_order_runge_kutta(f,g,t,x,y,dt):
	def k_and_l_values(f,g,t,x,y,dt):
		l0, k0 = f(t,x,y)*dt, 					g(t,x,y)*dt
		l1, k1 = f(t+dt/2,x+l0/2,y+k0/2)*dt, 	g(t+dt/2,x+l0/2,y+k0/2)*dt
		l2, k2 = f(t+dt/2,x+l1/2,y+k1/2)*dt, 	g(t+dt/2,x+l1/2,y+k1/2)*dt
		l3, k3 = f(t+dt,x+l2,y+k2)*dt, 			g(t+dt,x+l2,y+k2)*dt
		return([l0,l1,l2,l3],[k0,k1,k2,k3])
	l,k = k_and_l_values(f,g,t,x,y,dt)
	next_x = x + (l[0] + 2*(l[1] + l[2]) + l[3])/6
	next_y = y + (k[0] + 2*(k[1] + k[2]) + k[3])/6
	return(next_x,next_y)

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
def define_f(Length,LengthFirstDeriv,LengthSecondDeriv):
	def f(t,y,T):
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
	return(f)

def g(t,y,T):
	result = y
	return(result)
def k_and_l_values(f,g,t,x,y,dt):
	l0, k0 = f(t,x,y)*dt, 					g(t,x,y)*dt
	l1, k1 = f(t+dt/2,x+l0/2,y+k0/2)*dt, 	g(t+dt/2,x+l0/2,y+k0/2)*dt
	l2, k2 = f(t+dt/2,x+l1/2,y+k1/2)*dt, 	g(t+dt/2,x+l1/2,y+k1/2)*dt
	l3, k3 = f(t+dt,x+l2,y+k2)*dt, 			g(t+dt,x+l2,y+k2)*dt
	return([l0,l1,l2,l3],[k0,k1,k2,k3])
def fourth_order_runge_kutta(f,g,t,x,y,dt):
	l,k = k_and_l_values(f,g,t,x,y,dt)
	next_x = x + (l[0] + 2*(l[1] + l[2]) + l[3])/6
	next_y = y + (k[0] + 2*(k[1] + k[2]) + k[3])/6
	return(next_x,next_y)

Bag1_RK = {'TensionFirstDeriv':[0],'Tension':[0]}
dt = 1e-5*(3/5)
t = np.arange(0,15+dt,dt)
Length = 1.0766863061118987
LengthFirstDeriv = 0
LengthSecondDeriv = 0
f = define_f(Length,LengthFirstDeriv,LengthSecondDeriv)
StartTime = time.time()
for i in range(len(t)):
	append_dictionary(Bag1_RK,['TensionFirstDeriv','Tension'],fourth_order_runge_kutta(f,g,t[i],Bag1_RK['TensionFirstDeriv'][-1],Bag1_RK['Tension'][-1],dt))
	statusbar(i,len(t),StartTime=StartTime,Title='4th Order Runge Kutta')
Bag1_RK['Tension'] = Bag1_RK['Tension'][:-1]
t_RK, bag1tension_RK = t[np.arange(0,len(t),10)],np.array(Bag1_RK['Tension'])[np.arange(0,len(t),10)]
del(t,Bag1_RK)
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
Bag1_FE = {'GammaDynamicGain':50,'DynamicSpindleFrequency':[0],'Tension':[0],'TensionFirstDeriv':[0]}
SamplingPeriod = 1e-7
Time = np.arange(0,15+SamplingPeriod,SamplingPeriod)
StartTime = time.time()
for i in range(len(Time)):
	bag1_model(CE,Bag1_FE,SamplingPeriod)
	statusbar(i,len(Time),StartTime=StartTime,Title='Forward Euler')
Bag1_FE['Tension'] = Bag1_FE['Tension'][:-1]
Time_FE, bag1tension_FE = Time[np.arange(0,len(Time),100)],np.array(Bag1_FE['Tension'])[np.arange(0,len(Time),100)]
del(Time,Bag1_FE)
# Time,Bag1Tension = pickle.load(open('Bag1Tension_FE.pkl','rb'))

plt.plot(Time_FE,bag1tension_FE,'r--')
plt.plot(t_RK,bag1tension_RK,'b')
plt.show()
