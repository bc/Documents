import numpy as np
import matplotlib.pyplot as plt 

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
	ax.set_aspect('equal', 'datalim')
def plot_results(ax,X,Y,color):
	ax.plot(X,Y,color,lw = 2)
def PlanarTransformationMatrix(Angle, Displacement, PreviousTMatrix):
	"""
	This function will take in any angle (radians) of rotation about 
	the k axis and any displacement in the CURRENT FRAME OF REFERENCE

	"""
	NextTMatrix = np.matrix([ 	[np.cos(Angle),	-np.sin(Angle),	0,	Displacement[0]	],\
								[np.sin(Angle),	np.cos(Angle),	0,	Displacement[1]	],\
								[0,				0,				1,				0	],\
								[0,				0,				0,				1	]	])
	T = PreviousTMatrix*NextTMatrix
	Endpoint = np.array(T[:2,3].T)[0]
	return(T,Endpoint)
def rotate_data(X,RotationAngle):
	R = np.matrix([	[np.cos(RotationAngle), -np.sin(RotationAngle)],\
					[np.sin(RotationAngle), np.cos(RotationAngle)]	])
	result = np.array(R*X)
	return(result)
def translate_data(X,Translation):
	result = X + np.array([Translation]*len(X.T)).T
	return(result)
def plot_cap(ax, Radius, Orientation, RotationAngle, Translation, Color):
	"""
	Center is Cartesian point of the center of the cap. Orientation should 
	be either 'up' or 'down' to designate if the cap is for the start or 
	end of the link. RotationAngle rotates the link from the reference posture
	(link oriented downwards). Rotation angle should be given in radians.
	"""
	if Orientation == 'up':
		OrientationFactor = 1
	elif Orientation == 'down':
		OrientationFactor = -1
	Theta = np.arange(0,np.pi,0.001)
	X = OrientationFactor*np.array([Radius*np.cos(Theta), Radius*np.sin(Theta)]) 
	result = translate_data(rotate_data(X,RotationAngle),Translation)
	plot_results(ax,result[0], result[1], Color)
def plot_link_body(ax, Radius, Length, RotationAngle, Center, Color):
	"""
	Length is the length of the link, and RotationAngle (in radians) is the rotation
	from the reference posture (link oriented downwards). Radius is the same as the 
	radius of the cap (will be half of the width of the link body). Center is
	the starting point of this link relative to origin (0,0).
	"""
	RightEdge = np.array([	[Radius, 	Radius],\
							[0,			-Length]])
	LeftEdge = np.array([	[-Radius,	-Radius],\
							[0,			-Length]])
	RightEdge_rot = rotate_data(RightEdge,RotationAngle)
	LeftEdge_rot = rotate_data(LeftEdge,RotationAngle)
	RightResult = RightEdge_rot + np.array([Center,Center]).T
	LeftResult = LeftEdge_rot + np.array([Center,Center]).T
	plot_results(ax,RightResult[0],RightResult[1],Color)
	plot_results(ax,LeftResult[0],LeftResult[1],Color)
def plot_link(ax, Center, Radius, Length, TotalRotation, Translation, Color):
	"""
	Takes in the starting point of the link (Center), the desired radius of the 
	limb (Radius), the length of the link (Length), and the total rotation that
	link has experiences (TotalRotation), as well as the Translation to the new
	endpoint caused by the sum of rotations and the translation it previous
	experienced (Center).
	"""
	plot_cap(ax,Radius,'up',TotalRotation,Center,Color)
	plot_link_body(ax,Radius,Length,TotalRotation,Center,Color)
	plt.plot(Translation[0],Translation[1],'ko')
	plot_cap(ax,Radius,'down',TotalRotation,Translation,Color)
def plot_2_link_planar_model(Angle1,Angle2,Length1,Length2):
	"""
	This takes in two angles and two link lengths and returns a figure of the
	configuration with starting posture (0,0) with both links hanging vertical
	from origin (joint 2 at (0,-Length1) and endpoint at (0,-Length1-Length2))
	"""
	T01,Joint1 = PlanarTransformationMatrix(Angle1,[0,0],np.identity(4))
	T12,Joint2 = PlanarTransformationMatrix(Angle2,[0,-Length1],T01)
	T23,Endpoint = PlanarTransformationMatrix(0,[0,-Length2],T12)
	
	Radius = 0.5
	plt.figure()
	ax = plt.gca()
	plt.plot([0],[0],'ko')
	plt.plot([Joint1[0], Joint2[0], Endpoint[0]],[Joint1[1], Joint2[1], Endpoint[1]],"0.75",lw=3)
	plot_link(ax, Joint1, Radius, Length1, Angle1, Joint2,"0.55")
	plot_link(ax, Joint2, Radius, Length2, Angle1+Angle2, Endpoint, "0.55")
	quick_2D_plot_tool(ax,'x','y','Drawing Rotating Links')
	return(ax)

ax1 = plot_2_link_planar_model(np.pi/6,np.pi/6,5,5)
ax2 = plot_2_link_planar_model(np.pi/2,np.pi/2,4,3)
plt.imshow((ax1,ax2))

from moviepy.editor import VideoClip

def make_frame(t):
    frame_for_time_t = plot_2_link_planar_model(Angle1[t],Angle2[t],5,5)
    return frame_for_time_t # (Height x Width x 3) Numpy array
Angle1 = np.arange(0,np.pi,0.001)
Angle2 = np.arange(0,np.pi/2,0.0005)
animation = VideoClip(make_frame, duration=3) # 3-second clip

# For the export, many options/formats/optimizations are supported
animation.write_videofile("my_animation.mp4", fps=24) # export as video
animation.write_gif("my_animation.gif", fps=24) # export as GIF (slow)