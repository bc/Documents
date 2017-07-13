def statusbar(i,N,**kwargs):
	"""
	i is the current iteration (must be an int) and N is the length of 
	the range (must be an int). i must also be in [0,N). 

	~~~~~~~~~~~~~~
	**kwargs
	~~~~~~~~~~~~~~

	StartTime should equal time.time() and should be defined before your
	loop to ensure that you get an accurate representation of elapsed time.

	Title should be a str that will be displayed before the statusbar. Title
	should be no longer than 25 characters.

	~~~~~~~~~~~~~~

	NOTE: you should place a print('\n') after the loop to ensure you
	begin printing on the next line.

	"""
	import time
	from scipy import interpolate
	StartTime = kwargs.get("StartTime",False)
	Title = kwargs.get("Title",'')

	assert type(i)==int, "i must be an int"
	assert type(N)==int, "N must be an int"
	assert N>i, "N must be greater than i"
	assert N>0, "N must be a positive integer"
	assert i>=0, "i must not be negative (can be zero)"
	assert type(Title) == str, "Title should be a string"
	assert len(Title) <= 25, "Title should be less than 25 characters"
	if Title != '': Title = ' '*(25-len(Title)) + Title + ': '
	statusbar = Title +'[' + '\u25a0'*int((i+1)/(N/50)) + '\u25a1'*(50-int((i+1)/(N/50))) + '] '
	if StartTime != False:
		global time_array
		global TimeLeft
		if i==0:
			time_array = []
			TimeLeft = '--'
		elif i==int(0.02*N):
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(time_array[-1]*(N/(i+1)))
		elif i%int(0.02*N)==0:
			time_array.append(time.time()-StartTime)
			TimeLeft = '{0:1.1f}'.format(float(interpolate.interp1d(np.arange(len(time_array)),time_array,fill_value='extrapolate')(49)))
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete, ' + '{0:1.1f}'.format(time.time() - StartTime) \
			+ 'sec, (est. ' + TimeLeft + ' sec left)        \r', end='')
	else:
		print(statusbar + '{0:1.1f}'.format((i+1)/N*100) + '% complete           \r',end = '')

def planar_3dof_arm(theta1=45,l1=4*2**0.5,l2=4,l3=4,X=[0,20],plot=False,range_of_motion=[[0,170],[-75,60]],depth=20):
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

	def find_theta3(X,theta1,l1,l2,l3,depth):
		theta1 = theta1*pi/180
		x = X[0]
		y = X[1]
		if ((x-l1*cos(theta1))**2 + (y+depth-l1*sin(theta1))**2)**0.5 - l2 - l3 == 0: 
			output = 0
		else:
			output = acos(((x-l1*cos(theta1))**2 + (y+depth-l1*sin(theta1))**2 - l2**2 - l3**2)/(2*l2*l3))
		output=output*180/pi
		return(output)

	def find_beta(X,theta1,theta3,l1,l3,depth):
		theta1 = theta1*pi/180
		theta3 = theta3*pi/180
		x = X[0]
		y = X[1]
		if theta3 == 0:
			output = 0
		else:
			output = asin(l3*sin(pi-theta3)/(((x-l1*cos(theta1))**2 + (y+depth-l1*sin(theta1))**2)**0.5))
		return(output*180/pi)

	def find_theta2(X,theta1,beta,l1,depth):
		theta1 = theta1*pi/180
		beta = beta*pi/180
		x = X[0]
		y = X[1]
		diagonal = ((x-l1*cos(theta1))**2 + (y-(-depth+l1*sin(theta1)))**2)**0.5
		# output = acos((x*cos(theta1+beta)+y*sin(theta1+beta)-l1*cos(beta))/(((x-l1*cos(theta1))**2 + (y - l1*sin(theta1))**2)**0.5))-2*beta
		# output = asin((y*cos(theta1)-x*sin(theta1)+depth*cos(theta1))/diagonal)-beta
		output = acos((depth*sin(theta1)-l1+y*sin(theta1)+x*cos(theta1))/diagonal)-beta
		output = output*180/pi
		return(output)

	x,y = X[0],X[1]
	diagonal = ((x-l1*cos(pi*theta1/180))**2 + (y+depth-l1*sin(pi*theta1/180))**2)**0.5
	# if diagonal < abs(l2-l3) or diagonal > l2+l3:
	if diagonal > l2+l3:
		return([None,None,None],[None,None,None])
	else:
		try:
			theta3 = find_theta3(X,theta1,l1,l2,l3,depth)
			beta = find_beta(X,theta1,theta3,l1,l3,depth)
			theta2 = find_theta2(X,theta1,beta,l1,depth)

			if abs(l1*cos(pi*theta1/180)+l2*cos(pi*theta1/180+pi*theta2/180)+l3*cos(pi*theta1/180+pi*theta2/180+pi*theta3/180)-x)>0.01:
				theta2=-(theta2+2*beta)
				theta3=-theta3

			if theta3 != 0:
				theta3_alt = -theta3
				theta2_alt = theta2+2*beta

			if plot==True:
				limb_x = np.cumsum([0,l1*cos(pi*theta1/180),l2*cos(pi*theta1/180+pi*theta2/180),l3*cos(pi*theta1/180+pi*theta2/180+pi*theta3/180)])
				limb_y = np.cumsum([-depth,l1*sin(pi*theta1/180),l2*sin(pi*theta1/180+pi*theta2/180),l3*sin(pi*theta1/180+pi*theta2/180+pi*theta3/180)])
				limb_x_alt = np.cumsum([l1*cos(pi*theta1/180),l2*cos(pi*theta1/180+pi*theta2_alt/180),l3*cos(pi*theta1/180+pi*theta2_alt/180+pi*theta3_alt/180)])
				limb_y_alt = np.cumsum([-depth + l1*sin(pi*theta1/180),l2*sin(pi*theta1/180+pi*theta2_alt/180),l3*sin(pi*theta1/180+pi*theta2_alt/180+pi*theta3_alt/180)])

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
			if theta2<range_of_motion[0][0] or theta2>range_of_motion[0][1]:
				Config1 = [None,None,None]
			elif theta3<range_of_motion[1][0] or theta3>range_of_motion[1][1]:
				Config1 = [None,None,None]
			else:
				Config1 = [theta1,theta2,theta3]

			if theta2_alt<range_of_motion[0][0] or theta2_alt>range_of_motion[0][1]:
				Config2 = [None,None,None]
			elif theta3_alt<range_of_motion[1][0] or theta3_alt>range_of_motion[1][1]:
				Config2 = [None,None,None]
			else:
				Config2 = [theta1,theta2_alt,theta3_alt]	

			return(Config1,Config2)
		except ValueError:
			# pass
			return([None,None,None],[None,None,None])

def find_X_values(x_final):
	import numpy as np 
	from numpy import pi
	from math import sin, cos
	t = np.arange(0,1.1,0.1)
	x = -x_final[0]*(15*t**4-6*t**5-10*t**3)
	y = -x_final[1]*(15*t**4-6*t**5-10*t**3)
	X = np.concatenate(([x],[y]),axis=0)
	return(np.array(X))

import numpy as np
from numpy import pi
from math import cos, sin
import random
import matplotlib.pyplot as plt 
import matplotlib.animation as animation
import time
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

X = find_X_values([0,-20])
total2 = []
angle1 = np.arange(0,120,0.1)
movement_length = np.shape(X)[1]
StartTime = time.time()

HeightInInches = 71
Height = HeightInInches*2.54
l1 = 0.186*Height
l2 = 0.146*Height
HandLength = 0.108*Height
l3 = HandLength/2
reaching_depth = 40
for i in range(movement_length):
	# statusbar(i,movement_length,StartTime=StartTime,Title='45 Degree Reach')
	total_in_trial = []
	for j in range(len(angle1)):
		statusbar(len(angle1)*i+j,movement_length*len(angle1),StartTime=StartTime,Title='45 Degree Reach')
		config1,config2 = planar_3dof_arm(theta1=angle1[j],X=X[:,i],l1=l1,l2=l2,l3=l3,depth=reaching_depth)
		if config1 != [None,None,None]:
			config1.append(i/(movement_length-1))
			total_in_trial.append([config1])
		if config2 != [None,None,None]:
			config2.append(i/(movement_length-1))
			total_in_trial.append([config2])
	if total_in_trial != []: total2.append(total_in_trial)

# l1,l2,l3 = 4*2**0.5,4,4
limb_x = lambda theta: np.cumsum([0,l1*cos(pi*theta[0]/180),l2*cos(pi*theta[0]/180+pi*theta[1]/180),l3*cos(pi*theta[0]/180+pi*theta[1]/180+pi*theta[2]/180)])
limb_y = lambda theta: np.cumsum([-reaching_depth,l1*sin(pi*theta[0]/180),l2*sin(pi*theta[0]/180+pi*theta[1]/180),l3*sin(pi*theta[0]/180+pi*theta[1]/180+pi*theta[2]/180)])

# test = np.array(total[3]).T[:,0,:]
# fig, ax = plt.subplots()
# ax = plt.axes(xlim=(0, 5), ylim=(-1, 12))
# line, = ax.plot(limb_x(test[:,0]),limb_y(test[:,0]),'k',lw=2,marker = 'o',markersize = 10)
# # initialization function: plot the background of each frame
# def init():
#     line.set_data([], [])
#     return line,
# def animate(i):
# 	line.set_data(limb_x(test[:,i*400]),limb_y(test[:,i*400]))
# 	return(line,)
# quick_2D_plot_tool(ax,"X","Y",title=None)
# ani = animation.FuncAnimation(fig, animate, init_func=init,frames = int(np.shape(total[3])[0]/400), interval=20,blit=True)
# plt.show()
# import random
# random.seed()

# # testnumber = 0
# testnumber = random.randint(0,len(total))
# plt.figure()
# ax = plt.gca()
# test = np.concatenate(total[testnumber],axis=0)
# [plt.plot(limb_x(test[i,:]),limb_y(test[i,:]),'k',lw=2,marker='o',markersize=10) for i in range(np.shape(test)[0])]
# quick_2D_plot_tool(ax=ax,title="Test for single posture",xlabel='X',ylabel='Y')
# # plt.show()

# concat_test=[np.concatenate(total[i],axis=0) for i in range(np.shape(total)[0])] 
# concat_test=np.concatenate(concat_test)

# from mpl_toolkits.mplot3d import Axes3D

# random.seed()
# sample_index=random.sample(range(np.shape(concat_test)[0]),5000)

# fig2 = plt.figure()
# ax2 = fig2.add_subplot(111,projection= '3d')
# scatterplot = ax2.scatter(concat_test[sample_index,0],concat_test[sample_index,1],concat_test[sample_index,2],c=concat_test[sample_index,3]) 
# ax2.set_ylabel('Elbow')    
# ax2.set_xlabel('Shoulder')   
# ax2.set_zlabel('Wrist') 
# fig2.colorbar(scatterplot)
     
# plt.show()

def find_X_dot_values(x_final):
	import numpy as np 
	from numpy import pi
	from math import sin, cos
	t = np.arange(0,1.1,0.1)
	x = -x_final[0]*(60*t**3-30*t**4-30*t**2)
	y = -x_final[1]*(60*t**3-30*t**4-30*t**2)
	X = np.concatenate(([x],[y]),axis=0)
	return(np.array(X))
def jacobian(theta,l1,l2,l3):
	import numpy as np
	from numpy import pi
	from math import cos, sin
	theta1,theta2,theta3 = [pi*el/180 for el in theta]
	J = np.matrix([[-l1*sin(theta1)-l2*sin(theta1+theta2)-l3*sin(theta1+theta2+theta3),-l2*sin(theta1+theta2)-l3*sin(theta1+theta2+theta3),-l3*sin(theta1+theta2+theta3)],\
					[l1*cos(theta1)+l2*cos(theta1+theta2)+l3*cos(theta1+theta2+theta3),l2*cos(theta1+theta2)+l3*cos(theta1+theta2+theta3),l3*cos(theta1+theta2+theta3)]])
	return(J)
def is_acceptable_joint_velocity(next_X_dot,current_theta,next_theta,l1,l2,l3):
	current_theta,next_theta = np.array(current_theta),np.array(next_theta)
	delta_t = next_theta[0,3]-current_theta[0,3]
	theta_dot = np.array([[next_theta[0,0]-current_theta[0,0]],[next_theta[0,1]-current_theta[0,1]],[next_theta[0,2]-current_theta[0,2]]])*(1/delta_t)
	J = jacobian(next_theta[0,:3],l1,l2,l3)
	X_dot = np.array([[next_X_dot[0]],[next_X_dot[1]]])
	result = ((J*theta_dot-X_dot)<0.001).all()
	return(result)
def return_allowable_jump_indices(current,all_possible_next,X_dot,l1,l2,l3):
	import numpy as np
	from itertools import compress
	is_jump_acceptable = np.array([((np.array(current)-np.array(el))<3.84).all() for el in all_possible_next])
	allowable_jump_indices = list(compress(range(len(is_jump_acceptable)), is_jump_acceptable))
	if allowable_jump_indices != []:
		allowable_jumps = np.array([all_possible_next[el] for el in allowable_jump_indices])
		# ipdb.set_trace()
		with_correct_velocity_index = []
		for i in range(len(allowable_jumps)):
			if is_acceptable_joint_velocity(X_dot,current,allowable_jumps[i],l1,l2,l3):
				with_correct_velocity_index.append(allowable_jump_indices[i])
		return(with_correct_velocity_index)
	else:
		return([])
X_dot = find_X_dot_values([0,-20])
trajectories = []
for i in range(len(total2)-1):
	allowable_jump_indices = []
	statusbar(i,len(total2)-1)
	for j in range(len(total2[i])):
		allowable_jump_indices.append(return_allowable_jump_indices(total2[i][j],total2[i+1],X_dot.T[i+1],l1,l2,l3))
	trajectories.append(allowable_jump_indices)

