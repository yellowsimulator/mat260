import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from math import sqrt, sin, cos
"""
Animation of the zombi trajectory from problem 3 in homwork 1
"""

    
def eulers_method(h,sz):
    """
    solve the system of ode:
    
    xz' = sz*(R*cos(w*t) - xz)*1/alpha
    yz' = sz*(R*sin(w*t) - yz)*1/alpha
    alpha = sqrt((R*cos(w*t)-xz)**2 + (R*sin(w*t)-yz)**2)
    xz(0) = a
    yz(0) = -a
    0< a < sqrt(0.5)
    by Euler's method
    Y[n+1] = Y[n] + h*F[n]
    """
    T = 8*np.pi # simulation time t in [0, T=1]
    N = int(T/h) 
    t = np.linspace(0, T, N+1)
    xz = np.zeros(N+1)
    yz = np.zeros(N+1)
    
    #initial conditions
    a = sqrt(0.5)/2
    w = 1; R = 1
    xz[0] = a; yz[0] = -a
    
    for n in range(0, N):
        alpha = sqrt( (R*cos(w*t[n])-xz[n])**2 + (R*sin(w*t[n])-yz[n])**2   )
        xz[n+1] = xz[n] + ((sz*h)/alpha)*(R*cos(w*t[n])-xz[n])
        yz[n+1] = yz[n] + ((sz*h)/alpha)*(R*sin(w*t[n])-yz[n])
    return xz, yz, t
     
def particle_path(x,y,z):
    path = np.vstack((x, y, z))
    return path

def update_lines(num, dataLines, lines) :
    for line, data in zip(lines, dataLines) :
        line.set_data(data[0:2, num-1:num])
        line.set_3d_properties(data[2,num-1:num])
    return lines


if __name__ == "__main__":
    
    #set 3d plot and particle path properties 
    fig = plt.figure()
    ax = p3.Axes3D(fig)
    sz = 1.5
    h = 0.001
    T = 8*np.pi # simulation time t in [0, T=1]
    N = int(T/h)
    zz = np.zeros(N+1)
    xz, yz, t = eulers_method(h,sz)
    
    data = [particle_path(xz,yz,zz)]
    lines = [ax.plot(data[0][0,0:1], data[0][1,0:1], data[0][2,0:1], 'o')[0]]


    ax.set_xlim3d([-1.5, 1.5])
    ax.set_xlabel('x')
    ax.set_ylim3d([-1.5, 1.5])
    ax.set_ylabel('y')
    ax.set_zlim3d([-0.06, 0.06])
    ax.set_zlabel('z')
    

    # Create the animation object
    ani = animation.FuncAnimation(fig, update_lines, N, fargs=(data, lines),
                              interval=1, blit=False)
    plt.show()