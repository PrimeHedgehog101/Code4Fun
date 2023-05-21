import numpy as np 
import matplotlib.pyplot as plt 
import math 

##########################USER-CHANGE AS YOU WISH#############################

# Set stopping time and step length: 

Nstep = 5000
h     = 0.05

# Parameters for the mesh grid used to plot multiple particle paths in Fig1 
# Flow travels away from (0,0) so recomended initial conditions
# are (0,0) +- 0.05

x0_min = -0.05
x0_max = 0.05

y0_min = -0.05
y0_max = 0.05

# Choose how many blue flow lines (squared) are calculated:

Num_of_blue_lines_squared=5

# Initial conditions for single particle path in red: 

x0=-0.2
y0=0

# Different h values for comparason in Figure 3

h2=[0.005,0.01,0.1,0.2,0.3,0.5,1]  

# Set figure axis for figures 1,4,5:        (Figures 2,3 are auto)

Figure_1_4_5_Axis = [-5,5,-5,10]

#Toggle if you want to see chaos like behavoir from (1,0):

chaos='off'
#chaos='on'


#######################SETUP##################################################

# Creates two arrays of equally spaced points, used for plotting flow lines:

listx=np.linspace(x0_min,x0_max,Num_of_blue_lines_squared) 
listy=np.linspace(y0_min,y0_max,Num_of_blue_lines_squared)  


# Define function to compute velocity field:

def velfun(x,y,t): 
    
# Defines values to make things shorter:

    L1=x**2+y**2 
    L2=x**2+y**2-2*x+1 
    L3=x**2+y**2-x-y*np.sqrt(3)+1 
    
# Each P variable is a fraction, P1,P2,P3 are real, P4,P5,P6 are imagenary

    P1=(x*math.sqrt(3)-y)/(4*math.pi*L1) 
    P2=(y-x*math.sqrt(3)+math.sqrt(3))/(4*math.pi*L2) 
    P3=(y-0.5*math.sqrt(3))/(2*math.pi*L3) 

    P4=(x+y*math.sqrt(3))/(4*math.pi*L1) 
    P5=(-x+1-y*math.sqrt(3))/(4*math.pi*L2) 
    P6=(-x+0.5)/(2*math.pi*L3) 

    xvelocity = P1+P2+P3
    yvelocity = P4+P5+P6

    return xvelocity, yvelocity 

# Modified Runge-cutter method:

def loop(x0,y0,h): 
    
# Creates arrays of zeros:
    
    t     = np.linspace(0.,Nstep*h,Nstep) 
    x     = np.zeros(Nstep) 
    y     = np.zeros(Nstep) 
    
# Uses initial conditions set by user:
    
    x[0]=x0 
    y[0]=y0 
    
# Runge-cutter method:
    
    for k in range(1,Nstep): 
            
        # Provisional Euler step 
        tt    = t[k-1] 
        xx    = x[k-1] 
        yy    = y[k-1] 
        ux,uy = velfun(xx,yy,tt) 
        xp    = x[k-1] + h*ux 
        yp    = y[k-1] + h*uy 
        
        # Compute velocity at provisional point 
        uxp,uyp = velfun(xp,yp,tt+h) 
        
        # Compute average velocity 
        uxa = 0.5*(ux + uxp) 
        uya = 0.5*(uy + uyp)  
        
        # Make final Euler step 
        # using average velocity 
        x[k] = x[k-1] + h*uxa 
        y[k] = y[k-1] + h*uya 
        
# Stops chaos like behavour by restricting computation near (1,0):
        if chaos=='off':
            if 0.9<x[k]<1.1 and -0.1<y[k]<0.1:
                x=np.trim_zeros(x,'b')
                y=np.trim_zeros(y,'b')
                #x=np.delete(x,k)
                #y=np.delete(y,k)
                break
    
    return x,y

##############################GRAPHING########################################

# Fig1() is a plot of the blue flow lines using different initial conditions
# Stops division by zero at (0,0)
# Plots using a 2x2 grid of initial conditions around (0,0) 
  
def fig1(): 
    fig1=plt.figure()
    for i in listx: 
        for j in listy: 
            if not i==j==0:
                x,y=loop(i,j,h) 
                plt.plot(x,y,'b') 
    plt.xlabel('x') 
    plt.ylabel('y') 
    plt.grid() 
    plt.axis(Figure_1_4_5_Axis) 
    plt.title('Streamlines of flow')
    plt.savefig("Flow.pdf", bbox_inches='tight')
    return fig1 

# Fig2() is a plot of a single particle path in red using initial conditions
# x0 and y0 set by user

def fig2(): 
    fig2=plt.figure() 
    x,y=loop(x0,y0,h)   
    plt.plot(x,y,'r') 
    plt.xlabel('x') 
    plt.ylabel('y') 
    plt.grid() 
    plt.axis([-3,3,-2,3])
    plt.title('Figure 2 : Flow of idividual particle starting at ({},{})'.format(x0,y0))
    return fig2 

# Fig3() is a plot comparing different particle paths depending on different 
# h values set by user, using initial conditions x0,y0.

def fig3(): 
    fig3=plt.figure() 
    L=np.linspace(0,len(h2),len(h2)+1,dtype=int) 
    colour= ['r','g','c','k','m','y','b'] 
    for i,j in zip(h2,L): 
        x,y=loop(x0,y0,i) 
        plt.plot(x,y,colour[j],label='h={}'.format(i)) 
    #plt.xlabel('x') 
    #plt.ylabel('y') 
    plt.grid() 
    plt.legend() 
    plt.axis([1.8,2.05,0.4,1.4])
    #plt.title('Plot of different h values for x0=-0.1, y0=0.2')
    fig3.savefig("Hvaules3new2.pdf", bbox_inches='tight')
    return fig3 

# Figcomb1() is a combination of fig1 and fig2

def figcomb1(): 
    figcomb1=plt.figure() 
    fig1() 
    x,y=loop(x0,y0,h)    
    plt.plot(x,y,'r',label='Trajectory of particle starting at ({},{})'.format(x0,y0)) 
    plt.legend() 
    #plt.title('Flow of source independent of others')
    return figcomb1 

# Figcomb2() is a combination of fig1 and fig3

def figcomb2(): 
    figcomb2=plt.figure() 
    fig1() 
    L=np.linspace(0,len(h2),len(h2)+1,dtype=int) 
    colour= ['r','g','c','k','m','y','b'] 
    for i,j in zip(h2,L): 
       x,y=loop(x0,y0,i) 
       plt.plot(x,y,colour[j],label='h={}'.format(i)) 
    plt.legend() 
    plt.title('Figure 5 : Flow of idividual particle starting at ({},{})'.format(x0,y0))
    return figcomb2 




####################TOGGLE GRAPHS AS REQUIRED#################################

#fig1() 
#fig2() 
#fig3() 
figcomb1() 
#figcomb2() 


##################EXTRA STUFF###########################################


def figexact1():
    figexact1=plt.figure()
    a_values = np.linspace(-5, 5, 10)

    # Create a meshgrid for x and y values
    x = np.linspace(-10, 10, 1000)
    y = np.linspace(-10, 10, 1000)
    X, Y = np.meshgrid(x, y)

    # Iterate through each value of a and plot the equation
    for a in a_values:
        Z = np.log(np.abs(X)) - np.sqrt(3) * np.arctan(Y / X) + 0.5 * np.log((Y**2/X**2)+1)-a
        
        plt.contour(X, Y, Z, levels=[0])
    
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis([-3,3,-2,3])
    plt.title('Exact solution for (0,0)')
    return figexact1
    
    
def figexact2():
    figexact2=plt.figure()
    a_values = np.linspace(-20, 20, 12)

    x = np.linspace(-10, 10, 1000)
    y = np.linspace(-10, 10, 1000)
    X, Y = np.meshgrid(x, y)

    for a in a_values:
        a1 = np.sqrt(3)*Y+X-1
        a2 = np.sqrt(3)*(X-1)-Y
        Z = 6*np.sqrt(3) * np.arctan(a1/a2) - 3*np.log(np.abs(0.25*(Y**2/(X-1)**2+1))) - 6*np.log(np.abs(X-1))-a
        plt.contour(X, Y, Z, levels=[0])
    
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis([-3,3,-2,3])
    plt.title('Exact solution for (1,0)')
    return figexact2
  
def figexact3():
    figexact3=plt.figure()
    a_values = np.linspace(-1, 0, 6)

    x = np.linspace(-10, 10, 1000)
    y = np.linspace(-10, 10, 1000)
    X, Y = np.meshgrid(x, y)

    for a in a_values:
        Z = X**2+Y**2-np.sqrt(3)*Y-X-a
        plt.contour(X, Y, Z, levels=[0])
    
    plt.grid()
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis([-3,3,-2,3])
    plt.title('Exact solution for (0.5,0.866)')
    return figexact3
    
    
def figcomb3():
    figcomb3=plt.figure()
    fig1()
    a_values = np.linspace(-5, 5, 10)

    x = np.linspace(-10, 10, 1000)
    y = np.linspace(-10, 10, 1000)
    X, Y = np.meshgrid(x, y)

    for a in a_values:
        Z = np.log(np.abs(X)) - np.sqrt(3) * np.arctan(Y / X) + 0.5 * np.log((Y**2/X**2)+1)-a
        
        plt.contour(X, Y, Z, levels=[0])
    
    plt.xlabel('x')
    plt.ylabel('y')
    return figcomb3






#figexact1()
#figexact2() 
#figexact3()  

#figcomb3()
