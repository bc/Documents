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