# """
# Matplotlib Animation Example

# author: Jake Vanderplas
# email: vanderplas@astro.washington.edu
# website: http://jakevdp.github.com
# license: BSD
# Please feel free to use and modify this, but keep the above information. Thanks!
# """

# import numpy as np
# from matplotlib import pyplot as plt
# from matplotlib import animation

# # First set up the figure, the axis, and the plot element we want to animate
# fig = plt.figure()
# ax = plt.axes(xlim=(0, 2*np.pi), ylim=(-2, 2))
# line, = ax.plot([], [], lw=2)
# nframes = 200
# # initialization function: plot the background of each frame
# def init():
#     line.set_data([], [])
#     return line,

# # animation function.  This is called sequentially
# def animate(i):
#     x = np.linspace(0, (i/nframes)*2*np.pi, i)
#     y = np.sin(x)
#     line.set_data(x, y)
#     return line,

# # call the animator.  blit=True means only re-draw the parts that have changed.
# anim = animation.FuncAnimation(fig, animate, init_func=init,
#                                frames=200, interval=20, blit=False)

# # save the animation as an mp4.  This requires ffmpeg or mencoder to be
# # installed.  The extra_args ensure that the x264 codec is used, so that
# # the video can be embedded in html5.  You may need to adjust this for
# # your system: for more information, see
# # http://matplotlib.sourceforge.net/api/animation_api.html
# # anim.save('test1.mp4', fps=30, extra_args=['-vcodec', 'libx264'])

# plt.show()

"""
===========================
The double pendulum problem
===========================

This animation illustrates the double pendulum problem.
"""

# Double pendulum formula translated from the C code at
# http://www.physics.usyd.edu.au/~wheat/dpend_html/solve_dpend.c

from numpy import sin, cos
import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import matplotlib.animation as animation

G = 9.8  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum 1 in m
L2 = 1.0  # length of pendulum 2 in m
M1 = 1.0  # mass of pendulum 1 in kg
M2 = 1.0  # mass of pendulum 2 in kg


def derivs(state, t):
	[Angle1, AngularVelocity1, Angle2, AngularVelocity2] = state
	dydx = np.zeros_like(state)
	dydx[0] = AngularVelocity1

	del_ = Angle2 - Angle1
	den1 = (M1 + M2)*L1 - M2*L1*cos(Angle2 - Angle1)*cos(Angle2 - Angle1)
	dydx[1] = (M2*L1*AngularVelocity1**2*sin(Angle2 - Angle1)*cos(Angle2 - Angle1) +
	           M2*G*sin(Angle2)*cos(Angle2 - Angle1) +
	           M2*L2*AngularVelocity2**2*sin(Angle2 - Angle1) -
	           (M1 + M2)*G*sin(Angle1))/den1

	dydx[2] = AngularVelocity2

	den2 = (L2/L1)*den1
	dydx[3] = (-M2*L2*AngularVelocity2**2*sin(Angle2 - Angle1)*cos(Angle2 - Angle1) +
	           (M1 + M2)*G*sin(Angle1)*cos(Angle2 - Angle1) -
	           (M1 + M2)*L1*AngularVelocity1**2*sin(Angle2 - Angle1) -
	           (M1 + M2)*G*sin(Angle2))/den2

	return dydx

# create a time array from 0..100 sampled at 0.05 second steps
dt = 0.05
t = np.arange(0.0, 20, dt)

# Angle1 and Angle2 are the initial angles (degrees)
# AngularVelocity10 and AngularVelocity20 are the initial angular velocities (degrees per second)
Angle1 = 120.0
AngularVelocity1 = 0.0
Angle2 = -10.0
AngularVelocity2 = 0.0

# initial state
state = np.radians([Angle1, AngularVelocity1, Angle2, AngularVelocity2])

# integrate your ODE using scipy.integrate.
y = integrate.odeint(derivs, state, t)

x1 = L1*sin(y[:, 0])
y1 = -L1*cos(y[:, 0])

x2 = L2*sin(y[:, 2]) + x1
y2 = -L2*cos(y[:, 2]) + y1

fig = plt.figure()
ax = fig.add_subplot(111, autoscale_on=False, xlim=(-2, 2), ylim=(-2, 2))
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)


def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text


def animate(i):
    thisx = [0, x1[i], x2[i]]
    thisy = [0, y1[i], y2[i]]

    line.set_data(thisx, thisy)
    time_text.set_text(time_template % (i*dt))
    return line, time_text

ani = animation.FuncAnimation(fig, animate, np.arange(1, len(y)),
                              interval=25, blit=False, init_func=init)

# ani.save('double_pendulum.mp4', fps=15)
plt.show()