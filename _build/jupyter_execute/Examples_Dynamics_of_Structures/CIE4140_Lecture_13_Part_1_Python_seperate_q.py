#!/usr/bin/env python
# coding: utf-8

# # CIE4140_Lecture_13_Part_1_Python_seperate_q #

# ## Consider the steady-state vibrations of a string subjected to distributed viscous damping under a harmonic force ##

# In[1]:


import sympy as sp


# We represent the harmonic time-dependence by $e^{i \Omega t}$. The real part of the response to this load will give the response to the corresponding cosinusoidal load , whereas the imaginary part will give the response to the corresponding sinusoidal load

# Let us write the equation of motion

# In[2]:


w = sp.symbols('w',cls=sp.Function)
q1 = sp.symbols('q1',cls=sp.Function)
x,t = sp.symbols('x t')
nd, c, Omega= sp.symbols('nd c Omega',real=True)


# In[3]:


q1_second_normal_mode = sp.sin(2 * sp.pi * x / 1)
q1_arbitrary = 14.39 * x **5 * (x - 1)
q1_constant = 1
q1 = q1_arbitrary


# In[4]:


EQM = sp.diff(w(x,t),t,2) + 2 * nd * sp.diff(w(x,t),t) - c**2 * sp.diff(w(x,t),x,2) - q1*sp.exp(sp.I*Omega*t)
display(EQM)


# In[5]:


sp.plot(q1, (x, 0, 1));


# We search for the steady-state solution in the form:

# In[6]:


W = sp.symbols('W',cls=sp.Function)
w_form = W(x) * sp.exp(sp.I * Omega * t)


# Substitution of this equation in the equation of motion gives:

# In[7]:


EQM_fr = sp.simplify(EQM.subs(w(x,t),w_form))
display(EQM_fr)


# $e^{i \Omega t}$ is not dropped, so can dropped manually, not necessary:

# In[8]:


display(EQM_fr.subs(sp.exp(sp.I * Omega * t),1))


# Let us find the soluton to this equation that satisfies the boundary conditions of a fixed-fixed string

# In[9]:


L = sp.symbols('L',real=True)
W_sol = sp.dsolve(EQM_fr.subs(sp.exp(sp.I * Omega * t),1),W(x),ics={W(0):0,W(L):0}).rhs


# We introduce numerical values for the parameters

# In[10]:


W_sol = W_sol.subs([(L,1),(c,2),(nd,1/2)])


# Let us animate the reponses to these two loads assuming that the time-depandence is sinusoidal and introducing a numerical value of the load frequency: $\Omega=5$

# In[11]:


W_sol_sin = sp.im(W_sol*sp.exp(sp.I*Omega*t))


# In[12]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[13]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# In[14]:


W_sol_func = sp.lambdify((x,t,Omega),W_sol_sin)


# In[15]:


fig, ax = plt.subplots()
xdata = np.linspace(0,1,100)
line1, = ax.plot([], [])
ax.set_xlim(0, 1)
ax.set_ylim(-0.085,0.085)
frame = 0.1

def update(frame):
    ydata1 = W_sol_func(x=xdata,t=frame,Omega=5)
    ax.set_title("Displacement for t = "+str(np.round(frame,2)))
    line1.set_data(xdata, ydata1)

ani = FuncAnimation(fig, update, frames=np.linspace(0,2*np.pi/5,100),interval = 100)
plt.show()


# Now we plot the amplitude-frequency response functions at two locations along the string for the two load shapes

# In[16]:


AFRF_sol = sp.Abs(W_sol.subs(x,3/4))
sp.plot(AFRF_sol,(Omega,0.1,50),adaptive=False);


# Let us animate the reponses to these two loads assuming that the time-depandence is sinusoidal and introducing a numerical value of the load frequency: $\Omega=13$

# In[17]:


fig, ax = plt.subplots()
xdata = np.linspace(0,1,100)
line1, = ax.plot([], [])
ax.set_xlim(0, 1)
ax.set_ylim(-0.085,0.085)
frame = 0.1

def update(frame):
    ydata1 = W_sol_func(x=xdata,t=frame,Omega=13)
    ax.set_title("Displacement for t = "+str(np.round(frame,2)))
    line1.set_data(xdata, ydata1)

ani = FuncAnimation(fig, update, frames=np.linspace(0,2*np.pi/5,100),interval = 100)
plt.show()


# In[ ]:




