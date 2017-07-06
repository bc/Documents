def planar_3dof_arm(theta1=45,l1=4*2**0.5,l2=4,l3=4,X=[0,8],plot=False):
	import numpy as np
	from numpy import pi 
	from math import cos, sin, acos, asin
	import matplotlib.pyplot as plt 
	import ipdb

	def quick_2D_plot_tool(ax,xlabel,ylabel,title):
		"""
		This will take in the x and y labels as well as a figure title and
		format an already established figure to remove the box, place the
		tick marks outwards, and set aspect ratio to equal.
		"""
		ax = plt.gca()
		ax.spines['left'].set_position('zero')
		ax.spines['right'].set_color('none')
		ax.spines['bottom'].set_position('zero')
		ax.spines['top'].set_color('none')
		ax.xaxis.set_ticks_position('bottom')
		ax.yaxis.set_ticks_position('left')
		ax.set_xlabel(xlabel)
		ax.set_ylabel(ylabel)
		ax.set_title(title)
		ax.set_aspect('equal','datalim')

	def plot_results(ax,X,Y,color):
		ax.plot(X,Y,color,lw = 2)

	def find_theta3(X,theta1,l1,l2,l3):
		theta1 = theta1*pi/180
		x = X[0]
		y = X[1]
		output = acos(((x-l1*cos(theta1))**2 + (y-l1*sin(theta1))**2 - l2**2 - l3**2)/(2*l2*l3))
		return(output*180/pi)

	def find_beta(X,theta1,theta3,l1,l3):
		theta1 = theta1*pi/180
		theta3 = theta3*pi/180
		x = X[0]
		y = X[1]
		output = asin(-l3*sin(theta3)/(((x-l1*cos(theta1))**2 + (y - l1*sin(theta1))**2)**0.5))
		return(output*180/pi)

	def find_theta2(X,theta1,beta,l1):
		theta1 = theta1*pi/180
		beta = beta*pi/180
		x = X[0]
		y = X[1]
		output = acos((x*cos(theta1+beta)+y*sin(theta1+beta)-l1*cos(beta))/(((x-l1*cos(theta1))**2 + (y - l1*sin(theta1))**2)**0.5))+2*beta
		return(output*180/pi)

	theta3 = find_theta3(X,theta1,l1,l2,l3)
	beta = find_beta(X,theta1,theta3,l1,l3)
	theta2 = find_theta2(X,theta1,beta,l1)

	theta3_alt = -theta3
	theta2_alt = theta2-2*beta

	if plot==True:
		limb_x = np.cumsum([0,l1*cos(pi*theta1/180),l2*cos(pi*theta1/180+pi*theta2/180),l3*cos(pi*theta1/180+pi*theta2/180+pi*theta3/180)])
		limb_y = np.cumsum([0,l1*sin(pi*theta1/180),l2*sin(pi*theta1/180+pi*theta2/180),l3*sin(pi*theta1/180+pi*theta2/180+pi*theta3/180)])
		limb_x_alt = np.cumsum([l1*cos(pi*theta1/180),l2*cos(pi*theta1/180+pi*theta2_alt/180),l3*cos(pi*theta1/180+pi*theta2_alt/180+pi*theta3_alt/180)])
		limb_y_alt = np.cumsum([l1*sin(pi*theta1/180),l2*sin(pi*theta1/180+pi*theta2_alt/180),l3*sin(pi*theta1/180+pi*theta2_alt/180+pi*theta3_alt/180)])

		# ipdb.set_trace()
		assert np.sum(np.array([limb_x[-1],limb_y[-1]])-np.array([float(el) for el in X]))<0.0001, "Endpoints do not match."
		assert np.sum(np.array([limb_x_alt[-1],limb_y_alt[-1]])-np.array([float(el) for el in X]))<0.0001, "Endpoints do not match."

		plt.figure()
		plt.plot(limb_x_alt,limb_y_alt,'r',lw=2,marker = 'o',markersize = 10)
		plt.plot(limb_x,limb_y,'k',lw=2,marker = 'o',markersize = 10)
		ax = plt.gca()
		title = "beta: " +"{0:.2f}".format(beta)+", theta 1: " +"{0:.2f}".format(theta1)+ ", theta 2: " +"{0:.2f}".format(theta2)+ ", theta3: " +"{0:.2f}".format(theta3)
		quick_2D_plot_tool(ax,"X","Y",title)
		plt.show()
	return([theta1,theta2,theta3],[theta1,theta2_alt,theta3_alt])
