import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

N = 10001
Time = np.linspace(0,10,N)
dt = Time[1]-Time[0]

R = lambda t: 2*np.sin(2*np.pi*t)
dR = lambda t: (2*np.pi)*2*np.cos(2*np.pi*t)
d2R = lambda t: -( (2*np.pi)**2)*2*np.sin(2*np.pi*t)
def g1(x1,x2):
    return(x1**2 + x2)
def g2(x1,x2,u):
    return(u)
def update_policy(t,x1,x2,u,dt):
    x2.append(x2[-1] + g2(x1[-1],x2[-1],u)*dt)
    x1.append(x1[-1] + g1(x1[-1],x2[-2])*dt)
def input(t,x1,x2,k1,k2):
    u = (1+k1*k2)*R(t) + (k1 + k2)*dR(t) + d2R(t) \
            - (2*x1[-1] + k1 + k2)*g1(x1[-1],x2[-1]) - k1*k2*x1[-1]
    return(u)

k1,k2 = 5,5
x1 = [0]
x2 = [0]
for t in Time[1:]:
    update_policy(t,x1,x2,input(t,x1,x2,k1,k2),dt)

plt.figure()
plt.title(r'$\dot{x}_{1} = x_{1}^{2} + x_{2}; \hspace{1em} \dot{x}_{2} = u$',\
                fontsize=16,color='gray')
plt.plot(Time,x1,'b',lw=2)
plt.plot(Time,R(Time),'r--')
plt.xlabel("Time (s)")
plt.ylabel("Desired Measure")
plt.legend([r"Output $y = x_{1}$",r"Reference $r(t) = 2\sin(2\pi t)$"])

plt.figure()
plt.title('Error vs. Time')
plt.plot(Time, R(Time)-x1,color='r')
plt.xlabel("Time (s)")
plt.ylabel("Error")

# plt.show()

k1,k2 = 50,50
m1,m2,M = 1,1,1
A,w = 0.10,0.5*np.pi
b1,b2,b3,b4 = 20,20,20,20
CocontractionIndex = 2

def dx1(t,X):
    return(X[1])
def dx2(t,X):
    return(-(k1+k2)/M*X[0] + (k1/M)*X[2] + (k2/M)*X[3])
def dx3(t,X):
    return(X[4])
def dx4(t,X):
    return(X[5])
def dx5(t,X,U):
    return(k1/m1*X[0] -k1/m1*X[2] + U[0]/m1)
def dx6(t,X,U):
    return(k2/m2*X[0] -k2/m2*X[3] - U[1]/m2)

r = lambda t: A*np.sin(w*t)
dr = lambda t: A*w*np.cos(w*t)
d2r = lambda t: -A*w**2*np.sin(w*t)
d3r = lambda t: -A*w**3*np.cos(w*t)
d4r = lambda t: A*w**4*np.sin(w*t)
def z1(t,X):
    return(r(t) - X[0])
def dz1(t,X):
    return(dr(t) - X[1])
def d2z1(t,X):
    return(d2r(t) - dx2(t,X))
def d3z1(t,X):
    return(d3r(t) + (k1+k2)/M*X[1] - k1/M*X[4] - k2/M*X[5])
def a1(t,X):
    return(dr(t) - b1*z1(t,X))
def da1(t,X):
    return(d2r(t) - b1*dz1(t,X))
def z2(t,X):
    return(X[1] - a1(t,X))
def dz2(t,X):
    return(dx2(t,X) - da1(t,X))
def a2(t,X):
    return((k1+k2)/M*X[0] + (1 + b1*b2)*z1(t,X) + (b1+b2)*dz1(t,X) + d2r(t))
def da2(t,X):
    return((k1+k2)/M*X[1] + (1 + b1*b2)*dz1(t,X) + (b1+b2)*d2z1(t,X) + d3r(t))
def d2a2(t,X):
    return((k1+k2)/M*dx2(t,X) + (1 + b1*b2)*d2z1(t,X) + (b1+b2)*d3z1(t,X) + d4r(t))
def z3(t,X):
    return(k1/M*X[2] + k2/M*X[3] - a2(t,X))
def dz3(t,X):
    return(k1/M*X[4] + k2/M*X[5] - da2(t,X))
def a3(t,X):
    return(da2(t,X) - z2(t,X) -b3*z3(t,X))
def da3(t,X):
    return(d2a2(t,X) - dz2(t,X) -b3*dz3(t,X))
def z4(t,X):
    return(k1/M*X[4] + k2/M*X[5] - a3(t,X))
# def dz4(t,X):
#     return(k1/M*dx5(t,X,U) + k2/M*dx6(t,X,U) - da3(t,X))
def c1(t,X):
    return(-(k1**2*m2 + k2**2*m1)/(m1*m2*M)*X[0] + k1**2/(m1*M)*X[2] + k2**2/(m2*M)*X[3] + \
            da3(t,X) - z3(t,X) - b4*z4(t,X))
def return_U(t,X,e,Noise):
    if c1(t,X)<=0:
        u1 = ((m1*m2*M)/(k1*m2-e*k2*m1))*c1(t,X) + Noise[0]
        u2 = e*(u1-Noise[0]) + Noise[1]
    else:
        u2 = ((m1*m2*M)/(e*k1*m2-k2*m1))*c1(t,X) + Noise[1]
        u1 = e*(u2-Noise[1]) + Noise[0]
    return([u1,u2])
def animate_trajectory(response,Time,x1,x3,x4,u1,u2):
    assert type(response)==bool, "Input must be either True or False."

    if response == True:
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
        import matplotlib.animation as animation
        import matplotlib.patches as patches
        from scipy import signal

        fig = plt.figure(figsize=(10,8))
        ax1 = plt.subplot2grid((3,4),(0,0),colspan=4)
        ax2 = plt.subplot2grid((3,4),(1,0),colspan=2)
        ax3 = plt.subplot2grid((3,4),(1,2),colspan=2)
        ax4 = plt.subplot2grid((3,4),(2,0),colspan=3)
        ax5 = plt.subplot2grid((3,4),(2,3))

        plt.suptitle("Underdetermined Mass-Spring System",Fontsize=28,y=0.95)

        # Model Drawing
        IdealBoxScalingFactor = 0.78533496170320571 # Calculated from w = np.pi
        CurrentTrialScalingFactor = max([max(x3)-min(x1),max(x1)-min(x4)])
        StraightLength = 0.05*CurrentTrialScalingFactor/IdealBoxScalingFactor
        RestingLength = max([max(x1)-min(x3),max(x4)-min(x1)])+2*StraightLength\
                        +0.30*CurrentTrialScalingFactor/IdealBoxScalingFactor
        CenterBoxHalfWidth = 0.15*CurrentTrialScalingFactor/IdealBoxScalingFactor
        CenterBoxHalfHeight = 0.2*CurrentTrialScalingFactor/IdealBoxScalingFactor
        SideBoxHalfWidth = 0.1*CurrentTrialScalingFactor/IdealBoxScalingFactor
        SideBoxHalfHeight = 0.075*CurrentTrialScalingFactor/IdealBoxScalingFactor
        ForceScaling = 1*CurrentTrialScalingFactor/IdealBoxScalingFactor

        Spring_array =\
         SideBoxHalfWidth\
            *np.abs(signal.sawtooth(5*2*np.pi*np.linspace(0,1,1001)-np.pi/2))\
                -(1/2)*SideBoxHalfWidth

        Spring1, =\
            ax1.plot(np.linspace(x1[0]+CenterBoxHalfWidth+StraightLength,\
                                    RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,1001),\
                                        Spring_array,'k')
        Spring1_left, = ax1.plot([x1[0]+CenterBoxHalfWidth,x1[0]+CenterBoxHalfWidth+StraightLength],[0,0],'k')
        Spring1_right, = \
            ax1.plot([RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,\
                        RestingLength+x3[0]-SideBoxHalfWidth],\
                            [0,0],'k')

        Spring2, =\
            ax1.plot(np.linspace(-RestingLength+x4[0]+SideBoxHalfWidth+StraightLength,\
                                    x1[0]-CenterBoxHalfWidth-StraightLength,1001),\
                                        Spring_array,'k')
        Spring2_left, = ax1.plot([x1[0]-CenterBoxHalfWidth-StraightLength,x1[0]-CenterBoxHalfWidth],[0,0],'k')
        Spring2_right, = \
            ax1.plot([-RestingLength+x4[0]+SideBoxHalfWidth,\
                        -RestingLength+x4[0]+SideBoxHalfWidth+StraightLength],\
                            [0,0],'k')
        ax1.get_xaxis().set_ticks([])
        ax1.get_yaxis().set_ticks([])
        ax1.set_frame_on(True)
        CenterMass = plt.Rectangle((-CenterBoxHalfWidth,-CenterBoxHalfHeight),\
                                    2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
        ax1.add_patch(CenterMass)
        Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength,-SideBoxHalfHeight),\
                                    2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
        ax1.add_patch(Mass1)
        Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength,-SideBoxHalfHeight),\
                                    2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
        ax1.add_patch(Mass2)

        PositionArrow, = ax1.plot([x1[0],x1[0]],[0,2*CenterBoxHalfHeight],'k')
        PositionArrowHead, = ax1.plot([x1[0]],[2*CenterBoxHalfHeight],'k^')
        PositionArrowTail, = ax1.plot([x1[0]],[0],'ko')

        Scale = ax1.plot([-1.1*A,1.1*A],\
                            [2.75*CenterBoxHalfHeight,2.75*CenterBoxHalfHeight],\
                                '0.60')
        Ticks = np.linspace(-A,A,5)
        TickHeights = [0.3*CenterBoxHalfHeight,\
                        0.15*CenterBoxHalfHeight,\
                        0.3*CenterBoxHalfHeight,\
                        0.15*CenterBoxHalfHeight,\
                        0.3*CenterBoxHalfHeight]
        [ax1.plot([Ticks[i],Ticks[i]],\
                [2.75*CenterBoxHalfHeight-TickHeights[i],2.75*CenterBoxHalfHeight],'0.60') \
                    for i in range(5)]

        Force1Arrow, = ax1.plot([RestingLength+x3[0]+(5/3)*SideBoxHalfWidth,\
                                    RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
                                        +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],\
                                            [0,0],'g')
        Force1ArrowHead, = \
            ax1.plot([RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
                        +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],[0],'g>')
        Force2Arrow, =\
            ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
                        -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:])),\
                            x4[0]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
        Force2ArrowHead, = \
            ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
                        -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:]))],[0],'r<')

        LowerBound = (np.array(x4[5001:])-RestingLength-(5/3)*SideBoxHalfWidth\
                        -ForceScaling*np.array(u2[5000:])/(max(u1[5000:]+u2[5000:]))).min()
        UpperBound = (RestingLength + np.array(x3[5001:])+(5/3)*SideBoxHalfWidth\
                        +ForceScaling*np.array(u1[5000:])/(max(u1[5000:]+u2[5000:]))).max()
        Bound = 1.05*np.array([-LowerBound,UpperBound]).max()
        ax1.set_xlim([-Bound,Bound])
        ax1.set_ylim([-1.5*CenterBoxHalfHeight,3.25*CenterBoxHalfHeight])
        ax1.set_aspect('equal')

        #Force 1

        Force1, = ax3.plot([0],[u1[0]],color = 'g')
        ax3.set_xlim(0,Time[-1])
        ax3.set_xticks(list(np.linspace(0,Time[-1],5)))
        ax3.set_xticklabels([str(0),'','','',str(Time[-1])])
        ax3.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
        if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                        int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
            ax3.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
            ax3.set_yticklabels([""]*(int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))
        else:
            NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
            MaxTick = NumTicks - NumTicks%5
            TickStep = MaxTick/5
            Ticks = list(np.linspace(0,TickStep*5,6))
            ax3.set_yticks(Ticks)
            ax3.set_yticklabels([""]*len(Ticks))
        # ax3.set_yticklabels([str(int(el)) for el in \
        #                         list(np.linspace(0,\
        #                             np.ceil(max(u1[int(len(u1)/2):])*1.1) - \
        #                                 np.ceil(max(u1[int(len(u1)/2):])*1.1)%3,4))],\
        #                                     fontsize=12)
        ax3.spines['right'].set_visible(False)
        ax3.spines['top'].set_visible(False)
        ax3.set_title("Force 1",fontsize=16,fontweight = 4,color = 'g',y = 0.95)
        # ax3.set_xlabel("Time (s)")

        #Force 2

        Force2, = ax2.plot([0],[u2[0]],color = 'r')
        ax2.set_xlim(0,Time[-1])
        ax2.set_xticks(list(np.linspace(0,Time[-1],5)))
        ax2.set_xticklabels([str(0),'','','',str(Time[-1])])
        ax2.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
        ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                        int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
        ax2.set_yticklabels([str(int(el)) for el in \
                                list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                    int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
                                        fontsize=12)
        if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                        int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
            ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
            ax2.set_yticklabels([str(int(el)) for el in \
                                    list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                        int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
                                            fontsize=12)
        else:
            NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
            MaxTick = NumTicks - NumTicks%5
            TickStep = MaxTick/5
            Ticks = list(np.linspace(0,TickStep*5,6))
            ax2.set_yticks(Ticks)
            ax2.set_yticklabels([str(tick) for tick in Ticks])
        ax2.spines['right'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.set_title("Force 2",fontsize=16,fontweight = 4,color = 'r',y = 0.95)
        # ax2.set_xlabel("Time (s)")

        # Trajectory

        Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
        Actual, = ax4.plot([0],[x1[0]],'b')
        ax4.set_xlim(0,Time[-1])
        ax4.set_xticks(list(np.linspace(0,Time[-1],5)))
        ax4.set_xticklabels([str(0),'','','',str(Time[-1])])
        ax4.set_ylim([-1.25*A,1.25*A])
        ax4.set_yticks([-A,0,A])
        ax4.set_xlabel("Time (s)")
        ax4.set_ylabel("Position of Center Mass (m)")
        ax4.spines['right'].set_visible(False)
        ax4.spines['top'].set_visible(False)

        # Error
        ErrorArray = x1-r(Time)
        Error, = ax5.plot([0],[ErrorArray[0]],'k')
        ax5.set_xlim(0,Time[-1])
        ax5.set_xticks(list(np.linspace(0,Time[-1],5)))
        ax5.set_xticklabels([str(0),'','','',str(Time[-1])])
        ax5.set_ylim([ErrorArray.min() - 0.1*(max(ErrorArray)-min(ErrorArray)),\
                        ErrorArray.max() + 0.1*(max(ErrorArray)-min(ErrorArray))])
        ax5.set_xlabel("Time (s)")
        ax5.set_ylabel("Error (m)")
        ax5.yaxis.set_label_position("right")
        ax5.yaxis.tick_right()
        ax5.spines['left'].set_visible(False)
        ax5.spines['top'].set_visible(False)

        def animate(i):
            Spring1.set_xdata(np.linspace(x1[i]+CenterBoxHalfWidth+StraightLength,\
                                    RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,1001))
            Spring1_left.set_xdata([x1[i]+CenterBoxHalfWidth,x1[i]+CenterBoxHalfWidth+StraightLength])
            Spring1_right.set_xdata([RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,\
                                        RestingLength+x3[i]-SideBoxHalfWidth])

            Spring2.set_xdata(np.linspace(-RestingLength+x4[i]+SideBoxHalfWidth+StraightLength,\
                                    x1[i]-CenterBoxHalfWidth-StraightLength,1001))
            Spring2_left.set_xdata([x1[i]-CenterBoxHalfWidth-StraightLength,x1[i]-CenterBoxHalfWidth])
            Spring2_right.set_xdata([-RestingLength+x4[i]+SideBoxHalfWidth,\
                                        -RestingLength+x4[i]+SideBoxHalfWidth+StraightLength])

            CenterMass.xy = (-CenterBoxHalfWidth + x1[i],-CenterBoxHalfHeight)
            Mass1.xy = (-SideBoxHalfWidth+RestingLength + x3[i],-SideBoxHalfHeight)
            Mass2.xy = (-SideBoxHalfWidth-RestingLength + x4[i],-SideBoxHalfHeight)
            PositionArrow.set_xdata([x1[i],x1[i]])
            PositionArrowHead.set_xdata([x1[i]])
            PositionArrowTail.set_xdata([x1[i]])
            Force1Arrow.set_xdata([RestingLength+x3[i]+(5/3)*SideBoxHalfWidth,\
                                    RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
                                        +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))])
            Force1ArrowHead.set_xdata([RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
                                        +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))])
            Force2Arrow.set_xdata([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
                                        -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:])),\
                                            x4[i]-RestingLength-(5/3)*SideBoxHalfWidth])
            Force2ArrowHead.set_xdata([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
                                        -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:]))])

            Force1.set_xdata(Time[:i])
            Force1.set_ydata(u1[:i])

            Force2.set_xdata(Time[:i])
            Force2.set_ydata(u2[:i])

            Actual.set_xdata(Time[:i])
            Actual.set_ydata(x1[:i])

            Error.set_xdata(Time[:i])
            Error.set_ydata(ErrorArray[:i])

            return Spring1,Spring1_left,Spring1_right,Spring2,Spring2_left,Spring2_right,CenterMass,Mass1,Mass2,Force1,Force2,Actual,Error,PositionArrow,PositionArrowHead,PositionArrowTail,Force1Arrow,Force1ArrowHead,Force2Arrow,Force2ArrowHead,

        # Init only required for blitting to give a clean slate.
        def init():
            Spring1, =\
                ax1.plot(np.linspace(x1[0]+CenterBoxHalfWidth+StraightLength,\
                                        RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,1001),\
                                            Spring_array,'k')
            Spring1_left, = \
                ax1.plot([x1[0]+CenterBoxHalfWidth,x1[0]+CenterBoxHalfWidth+StraightLength],\
                    [0,0],'k')
            Spring1_right, = \
                ax1.plot([RestingLength+x3[0]-SideBoxHalfWidth-StraightLength,\
                            RestingLength+x3[0]-SideBoxHalfWidth],[0,0],'k')
            Spring2, =\
                ax1.plot(np.linspace(-RestingLength+x4[0]+SideBoxHalfWidth+StraightLength,\
                                        x1[0]-CenterBoxHalfWidth-StraightLength,1001),\
                                            Spring_array,'k')
            Spring2_left, =\
                ax1.plot([x1[0]-CenterBoxHalfWidth-StraightLength,x1[0]-CenterBoxHalfWidth],\
                            [0,0],'k')
            Spring2_right, = \
                ax1.plot([-RestingLength+x4[0]+SideBoxHalfWidth,\
                            -RestingLength+x4[0]+SideBoxHalfWidth+StraightLength],[0,0],'k')

            CenterMass = plt.Rectangle((-CenterBoxHalfWidth,-CenterBoxHalfHeight),\
                                        2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(CenterMass)
            Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength,-SideBoxHalfHeight),\
                                        2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(Mass1)
            Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength,-SideBoxHalfHeight),\
                                        2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(Mass2)

            PositionArrow, = ax1.plot([x1[0],x1[0]],[0,2*CenterBoxHalfHeight],'k')
            PositionArrowHead, = ax1.plot([x1[0]],[2*CenterBoxHalfHeight],'k^')
            PositionArrowTail, = ax1.plot([x1[0]],[0],'ko')

            Force1Arrow, = ax1.plot([RestingLength+x3[0]+(5/3)*SideBoxHalfWidth,\
                                    RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
                                        +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],\
                                            [0,0],'g')
            Force1ArrowHead, = \
                ax1.plot([RestingLength + x3[0]+(5/3)*SideBoxHalfWidth\
                        +ForceScaling*u1[0]/(max(u1[5000:]+u2[5000:]))],[0],'g<')
            Force2Arrow, = ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
                        -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:])),\
                            x4[0]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
            Force2ArrowHead, = \
                ax1.plot([x4[0]-RestingLength-(5/3)*SideBoxHalfWidth\
                    -ForceScaling*u2[0]/(max(u1[5000:]+u2[5000:]))],[0],'r>')

            Force1, = ax3.plot([0],[u1[0]],color = 'g')
            Force2, = ax2.plot([0],[u2[0]],color = 'r')
            Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
            Actual, = ax4.plot([0],[x1[0]],'b')
            Error, = ax5.plot([0],[ErrorArray[0]],'k')

            Spring1.set_visible(False)
            Spring1_left.set_visible(False)
            Spring1_right.set_visible(False)
            Spring2.set_visible(False)
            Spring2_left.set_visible(False)
            Spring2_right.set_visible(False)
            CenterMass.set_visible(False)
            Mass1.set_visible(False)
            Mass2.set_visible(False)
            PositionArrow.set_visible(False)
            PositionArrowHead.set_visible(False)
            PositionArrowTail.set_visible(False)
            Force1.set_visible(False)
            Force2.set_visible(False)
            Predicted.set_visible(False)
            Actual.set_visible(False)
            Error.set_visible(False)
            Force1Arrow.set_visible(False)
            Force1ArrowHead.set_visible(False)
            Force2Arrow.set_visible(False)
            Force2ArrowHead.set_visible(False)

            return Spring1,Spring1_left,Spring1_right,Spring2,Spring2_left,Spring2_right,CenterMass,Mass1,Mass2,Force1,Force2,Actual,Error,PositionArrow,PositionArrowHead,PositionArrowTail,Force1Arrow,Force1ArrowHead,Force2Arrow,Force2ArrowHead,

        ani = animation.FuncAnimation(fig, animate, np.arange(1, len(Time),10), init_func=init,interval=1, blit=False)
        # if save_as_gif:
        # 	ani.save('test.gif', writer='imagemagick', fps=30)
        plt.show()
def plot_multiple_PDF_frames(response,Time,x1,x3,x4,u1,u2,FileName=None):
    assert type(response)==bool, "Input must be either True or False."
    if FileName != None: assert type(FileName)==str, "FileName must be a string"
    if response == True:
        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib.patches import Ellipse
        import matplotlib.patches as patches
        from scipy import signal
        from matplotlib.backends.backend_pdf import PdfPages
        import os.path

        def return_fig(i):
            fig = plt.figure(figsize=(10,8))
            ax1 = plt.subplot2grid((3,4),(0,0),colspan=4)
            ax2 = plt.subplot2grid((3,4),(1,0),colspan=2)
            ax3 = plt.subplot2grid((3,4),(1,2),colspan=2)
            ax4 = plt.subplot2grid((3,4),(2,0),colspan=3)
            ax5 = plt.subplot2grid((3,4),(2,3))

            plt.suptitle("Underdetermined Mass-Spring System",Fontsize=28,y=0.95)

            # Model Drawing
            IdealBoxScalingFactor = 0.78533496170320571 # Calculated from w = np.pi
            CurrentTrialScalingFactor = max([max(x3)-min(x1),max(x1)-min(x4)])
            StraightLength = 0.05*CurrentTrialScalingFactor/IdealBoxScalingFactor
            RestingLength = max([max(x1)-min(x3),max(x4)-min(x1)])+2*StraightLength\
                            +0.30*CurrentTrialScalingFactor/IdealBoxScalingFactor
            CenterBoxHalfWidth = 0.15*CurrentTrialScalingFactor/IdealBoxScalingFactor
            CenterBoxHalfHeight = 0.2*CurrentTrialScalingFactor/IdealBoxScalingFactor
            SideBoxHalfWidth = 0.1*CurrentTrialScalingFactor/IdealBoxScalingFactor
            SideBoxHalfHeight = 0.075*CurrentTrialScalingFactor/IdealBoxScalingFactor
            ForceScaling = 1*CurrentTrialScalingFactor/IdealBoxScalingFactor

            Spring_array =\
             SideBoxHalfWidth\
                *np.abs(signal.sawtooth(5*2*np.pi*np.linspace(0,1,1001)-np.pi/2))\
                    -(1/2)*SideBoxHalfWidth

            Spring1, =\
                ax1.plot(np.linspace(x1[i]+CenterBoxHalfWidth+StraightLength,\
                                        RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,1001),\
                                            Spring_array,'k')
            Spring1_left, = \
                ax1.plot([x1[i]+CenterBoxHalfWidth,x1[i]+CenterBoxHalfWidth+StraightLength],\
                            [0,0],'k')
            Spring1_right, = \
                ax1.plot([RestingLength+x3[i]-SideBoxHalfWidth-StraightLength,\
                            RestingLength+x3[i]-SideBoxHalfWidth],\
                                [0,0],'k')

            Spring2, =\
                ax1.plot(np.linspace(-RestingLength+x4[i]+SideBoxHalfWidth+StraightLength,\
                                        x1[i]-CenterBoxHalfWidth-StraightLength,1001),\
                                            Spring_array,'k')
            Spring2_left, = \
                ax1.plot([x1[i]-CenterBoxHalfWidth-StraightLength,x1[i]-CenterBoxHalfWidth],\
                            [0,0],'k')
            Spring2_right, = \
                ax1.plot([-RestingLength+x4[i]+SideBoxHalfWidth,\
                            -RestingLength+x4[i]+SideBoxHalfWidth+StraightLength],\
                                [0,0],'k')
            ax1.get_xaxis().set_ticks([])
            ax1.get_yaxis().set_ticks([])
            ax1.set_frame_on(True)
            CenterMass = plt.Rectangle((-CenterBoxHalfWidth + x1[i],-CenterBoxHalfHeight),\
                                        2*CenterBoxHalfWidth,2*CenterBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(CenterMass)
            Mass1 = plt.Rectangle((-SideBoxHalfWidth+RestingLength + x3[i],-SideBoxHalfHeight),\
                                        2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(Mass1)
            Mass2 = plt.Rectangle((-SideBoxHalfWidth-RestingLength + x4[i],-SideBoxHalfHeight),\
                                        2*SideBoxHalfWidth,2*SideBoxHalfHeight,Color='#4682b4')
            ax1.add_patch(Mass2)

            PositionArrow, = ax1.plot([x1[i],x1[i]],[0,2*CenterBoxHalfHeight],'k')
            PositionArrowHead, = ax1.plot([x1[i]],[2*CenterBoxHalfHeight],'k^')
            PositionArrowTail, = ax1.plot([x1[i]],[0],'ko')

            Scale = ax1.plot([-1.1*A,1.1*A],\
                                [2.75*CenterBoxHalfHeight,2.75*CenterBoxHalfHeight],\
                                    '0.60')
            Ticks = np.linspace(-A,A,5)
            TickHeights = [0.3*CenterBoxHalfHeight,\
                            0.15*CenterBoxHalfHeight,\
                            0.3*CenterBoxHalfHeight,\
                            0.15*CenterBoxHalfHeight,\
                            0.3*CenterBoxHalfHeight]
            [ax1.plot([Ticks[i],Ticks[i]],\
                    [2.75*CenterBoxHalfHeight-TickHeights[i],2.75*CenterBoxHalfHeight],'0.60') \
                        for i in range(5)]

            Force1Arrow, = ax1.plot([RestingLength+x3[i]+(5/3)*SideBoxHalfWidth,\
                                        RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
                                            +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))],\
                                                    [0,0],'g')
            Force1ArrowHead, = \
                ax1.plot([RestingLength + x3[i]+(5/3)*SideBoxHalfWidth\
                            +ForceScaling*u1[i]/(max(u1[5000:]+u2[5000:]))],[0],'g>')
            Force2Arrow, =\
                ax1.plot([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
                            -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:])),\
                                x4[i]-RestingLength-(5/3)*SideBoxHalfWidth],[0,0],'r')
            Force2ArrowHead, = \
                ax1.plot([x4[i]-RestingLength-(5/3)*SideBoxHalfWidth\
                            -ForceScaling*u2[i]/(max(u1[5000:]+u2[5000:]))],[0],'r<')

            LowerBound = (np.array(x4[5001:])-RestingLength-(5/3)*SideBoxHalfWidth\
                            -ForceScaling*np.array(u2[5000:])/(max(u1[5000:]+u2[5000:]))).min()
            UpperBound = (RestingLength + np.array(x3[5001:])+(5/3)*SideBoxHalfWidth\
                            +ForceScaling*np.array(u1[5000:])/(max(u1[5000:]+u2[5000:]))).max()
            Bound = 1.05*np.array([-LowerBound,UpperBound]).max()
            ax1.set_xlim([-Bound,Bound])
            ax1.set_ylim([-1.5*CenterBoxHalfHeight,3.25*CenterBoxHalfHeight])
            ax1.set_aspect('equal')

            #Force 1

            Force1, = ax3.plot(Time[:i],u1[:i],color = 'g')
            ax3.set_xlim(0,Time[-1])
            ax3.set_xticks(list(np.linspace(0,Time[-1],5)))
            ax3.set_xticklabels([str(0),'','','',str(Time[-1])])
            ax3.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
            if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
                ax3.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
                ax3.set_yticklabels([""]*(int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))
            else:
                NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
                MaxTick = NumTicks - NumTicks%5
                TickStep = MaxTick/5
                Ticks = list(np.linspace(0,TickStep*5,6))
                ax3.set_yticks(Ticks)
                ax3.set_yticklabels([""]*len(Ticks))
            # ax3.set_yticklabels([str(int(el)) for el in \
            #                         list(np.linspace(0,\
            #                             np.ceil(max(u1[int(len(u1)/2):])*1.1) - \
            #                                 np.ceil(max(u1[int(len(u1)/2):])*1.1)%3,4))],\
            #                                     fontsize=12)
            ax3.spines['right'].set_visible(False)
            ax3.spines['top'].set_visible(False)
            ax3.set_title("Force 1",fontsize=16,fontweight = 4,color = 'g',y = 0.95)
            # ax3.set_xlabel("Time (s)")

            #Force 2

            Force2, = ax2.plot(Time[:i],u2[:i],color = 'r')
            ax2.set_xlim(0,Time[-1])
            ax2.set_xticks(list(np.linspace(0,Time[-1],5)))
            ax2.set_xticklabels([str(0),'','','',str(Time[-1])])
            ax2.set_ylim(0,1.15*max(u1[5000:]+u2[5000:]))
            ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
            ax2.set_yticklabels([str(int(el)) for el in \
                                    list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                        int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
                                            fontsize=12)
            if np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1).shape[0] < 5:
                ax2.set_yticks(list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1)))
                ax2.set_yticklabels([str(int(el)) for el in \
                                        list(np.linspace(0,np.floor(1.15*max(u1[5000:]+u2[5000:])),\
                                            int(np.floor(1.15*max(u1[5000:]+u2[5000:])))+1))],\
                                                fontsize=12)
            else:
                NumTicks = np.floor(1.15*max(u1[5000:]+u2[5000:]))
                MaxTick = NumTicks - NumTicks%5
                TickStep = MaxTick/5
                Ticks = list(np.linspace(0,TickStep*5,6))
                ax2.set_yticks(Ticks)
                ax2.set_yticklabels([str(tick) for tick in Ticks])
            ax2.spines['right'].set_visible(False)
            ax2.spines['top'].set_visible(False)
            ax2.set_title("Force 2",fontsize=16,fontweight = 4,color = 'r',y = 0.95)
            # ax2.set_xlabel("Time (s)")

            # Trajectory

            Predicted, = ax4.plot(Time,r(Time),'0.60',linestyle='--')
            Actual, = ax4.plot(Time[:i],x1[:i],'b')
            ax4.set_xlim(0,Time[-1])
            ax4.set_xticks(list(np.linspace(0,Time[-1],5)))
            ax4.set_xticklabels([str(0),'','','',str(Time[-1])])
            ax4.set_ylim([-1.25*A,1.25*A])
            ax4.set_yticks([-A,0,A])
            ax4.set_xlabel("Time (s)")
            ax4.set_ylabel("Position of Center Mass (m)")
            ax4.spines['right'].set_visible(False)
            ax4.spines['top'].set_visible(False)

            # Error
            ErrorArray = x1-r(Time)
            Error, = ax5.plot(Time[:i],ErrorArray[:i],'k')
            ax5.set_xlim(0,Time[-1])
            ax5.set_xticks(list(np.linspace(0,Time[-1],5)))
            ax5.set_xticklabels([str(0),'','','',str(Time[-1])])
            ax5.set_ylim([ErrorArray.min() - 0.1*(max(ErrorArray)-min(ErrorArray)),\
                            ErrorArray.max() + 0.1*(max(ErrorArray)-min(ErrorArray))])
            ax5.set_xlabel("Time (s)")
            ax5.set_ylabel("Error (m)")
            ax5.yaxis.set_label_position("right")
            ax5.yaxis.tick_right()
            ax5.spines['left'].set_visible(False)
            ax5.spines['top'].set_visible(False)

            return(fig)
        i = 1
        if FileName == None:
            FileName = "ReferenceTrackingTest.pdf"
        else:
            FileName = FileName + ".pdf"
        if os.path.exists(FileName) == True:
            while os.path.exists(FileName) == True:
                i += 1
                FileName = FileName[:-4]
                FileName = FileName	+ "_" + "{:0>2d}".format(i) +".pdf"
        PDFFile = PdfPages(FileName)

        for i in np.linspace(0,len(Time)-1,201)[:-1]:
            t_i = int(i)
            fig = return_fig(t_i)
            PDFFile.savefig(fig)
            plt.close("all")
        PDFFile.close()

N = 10001
Time = np.linspace(0,10,N)
x1,x2,x3,x4,x5,x6 = [A],[-1],[0],[0],[0],[0]
u1,u2 = [],[]

AddNoise = True
if AddNoise == True:
    np.random.seed(seed=None)
    NoiseArray = np.random.normal(loc=0.0,scale=0.2,size=(2,len(Time)))
else:
    NoiseArray = np.zeros((2,len(Time)))

def update_policy(t,x1,x2,x3,x4,x5,x6,dt,NoiseArray,e=2):
    import numpy as np
    X = [x1[-1],x2[-1],x3[-1],x4[-1],x5[-1],x6[-1]]
    U = return_U(t,X,e,NoiseArray[:,int(t/dt)])
    u1.append(U[0])
    u2.append(U[1])
    x6.append(x6[-1] + dx6(t,X,U)*dt)
    x5.append(x5[-1] + dx5(t,X,U)*dt)
    x4.append(x4[-1] + dx4(t,X)*dt)
    x3.append(x3[-1] + dx3(t,X)*dt)
    x2.append(x2[-1] + dx2(t,X)*dt)
    x1.append(x1[-1] + dx1(t,X)*dt)

for t in Time[1:]:
    update_policy(t,x1,x2,x3,x4,x5,x6,dt,NoiseArray,e=CocontractionIndex)

plt.figure()
# plt.title(r'$\dot{x}_{1} = x_{1}^{2} + x_{2}; \hspace{1em} \dot{x}_{2} = u$',\
#                 fontsize=16,color='gray')
plt.title("Underdetermined Spring Example",\
                fontsize=16,color='gray')
plt.plot(Time,x1,'b',lw=2)
plt.plot(Time,r(Time),'r--')
plt.xlabel("Time (s)")
plt.ylabel("Desired Measure")
plt.legend([r"Output $y = x_{1}$",r"Reference $r(t) = " + str(A) + "\sin(" + str(w) + "t)$"])

plt.figure()
plt.title('Error vs. Time')
plt.plot(Time, r(Time)-x1,color='r')
plt.xlabel("Time (s)")
plt.ylabel("Error")
#
# plt.show()

plt.close('all')
plot_multiple_PDF_frames(False,Time,x1,x3,x4,u1,u2)
animate_trajectory(True,Time,x1,x3,x4,u1,u2)
