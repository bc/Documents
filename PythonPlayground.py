import matplotlib.pyplot as plt 

import numpy as np 
"""
mu, sigma = 100, 15
x = mu + sigma*np.random.randn(10000)

n, bins, patches = plt.hist(x, 50, normed = 1, facecolor = 'g', alpha = 0.75)

plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title('Histogram of IQ')
plt.text(60,0.025,r'$\mu=100,\ \sigma=15$')
plt.axis([40,160,0,0.03])
plt.grid(True)
plt.show()
"""
"""
person = raw_input("Enter your name: ")
string = "Hello " + person
print string

"""
#person = raw_input("Enter your name: ")
#string = "Hello " + person
#print string
"""

x = np.array([1.0,5.0,2.0])
y = np.array([7,4,1])
z = x**y
a = []
for i in z:
	a.append(int(i))
print a
=======
x = np.arange(0,3*3.14,0.01)
y = np.sin(x) + np.sin(2*x)
xtick = np.arange(0,3*3.14+0.01,1.57)
xticklabel = [r'0', r'$\frac{\pi}{2}$',r'$\pi$', r'$\frac{3\pi}{2}$',r'$2\pi$', r'$\frac{5\pi}{2}$',r'$3\pi$']
plt.plot(x, y)
plt.xlabel(r'$\theta$',fontsize = 16)
plt.ylabel(r'$\Gamma(\theta)$',fontsize = 16)
plt.xticks(xtick,xticklabel,fontsize = 16)
plt.setp(plt.xticks()[1])
plt.title(r'$\Gamma(\theta) = \mathrm{sin}(\theta) + \mathrm{sin}(2\theta)$', fontsize = 20)
plt.show()
"""
"""
HeightInches = 71
DegreesToRadianFactor = 3.14159/180
ReleaseAngle = 30*DegreesToRadianFactor
Angle1Final = 120*DegreesToRadianFactor
Angle2Final = 84*DegreesToRadianFactor
Angle3Final = -35*DegreesToRadianFactor
Height = 2.54*HeightInches
ShoulderToElbowHeight = 0.186*Height
ForearmLength = 0.146*Height
HandLength = 0.108*Height
FinalPositionInX = ShoulderToElbowHeight*np.sin(Angle1Final) \
		+ ForearmLength*np.sin(Angle1Final + Angle2Final) \
		+ HandLength*np.sin(Angle1Final + Angle2Final - Angle3Final)
print "Final Position in X: " + str(FinalPositionInX)
FinalPositionInY = -ShoulderToElbowHeight*np.cos(Angle1Final) \
		- ForearmLength*np.cos(Angle1Final + Angle2Final) \
		- HandLength*np.cos(Angle1Final + Angle2Final - Angle3Final)
print "Final Position in Y: " + str(FinalPositionInY)
InitialProjectileVelocity = np.sqrt(-490*((434.3+0.152*Height-FinalPositionInX-11.9*np.cos(ReleaseAngle))**2) \
							/((((np.cos(ReleaseAngle))**2)*(304.8-0.87*Height-FinalPositionInY)) \
							-(np.sin(ReleaseAngle)*np.cos(ReleaseAngle)*(434.3+0.152*Height-FinalPositionInX))))
print "Final Velocity is " + str(InitialProjectileVelocity/100) + " m/s"
"""
import pickle

y = pickle.load(open('test.pkl','rb'))