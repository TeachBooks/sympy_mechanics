#!/usr/bin/env python
# coding: utf-8

# # CIE_4140_Lecture_2_1_Python #

# In[1]:


import sympy as sp


# In[2]:


u = sp.symbols('u',cls=sp.Function)
t, omega_0 = sp.symbols('t omega_0')
C, s = sp.symbols('C, s')


# Equation governing free motion of an SDOF

# In[3]:


Equation_of_Motion= sp.diff(u(t),t,2)+omega_0 **2 * u(t)
display(Equation_of_Motion)


# Characteristic equation

# In[4]:


CharacteristicEquation = sp.simplify(Equation_of_Motion.subs(u(t), C*sp.exp(s*t)) / (C * sp.exp(s*t)))
display(CharacteristicEquation)


# The Eigenvalues

# In[5]:


EigenValues = sp.solve(CharacteristicEquation,s)
display(EigenValues)


# The Eigenfrequencies

# In[6]:


EigenFrequency = []
for k in [0,1]:
    EigenFrequency.append(EigenValues[k]/ sp.I)
    display(EigenFrequency[k])


# The General Solution

# In[7]:


C_1, C_2 = sp.symbols('C_1, C_2')
General_Solution_Complex = C_1*sp.exp(EigenValues[0]*t) + C_2*sp.exp(EigenValues[0]*t)
display(General_Solution_Complex)


# In[8]:


A_1, A_2 = sp.symbols('A_1, A_2')
General_Solution_Real= A_1*sp.sin(EigenFrequency[0]*t)+A_2*sp.cos(EigenFrequency[1]*t)
display(General_Solution_Real)


# In[9]:


A, phi = sp.symbols('A phi')
General_Solution_Real_Compact =A*sp.cos(EigenFrequency[1]*t+phi)
display(General_Solution_Real_Compact)


# Visualisation

# In[10]:


sp.plot(General_Solution_Real_Compact.subs([(A,1),(phi,0),(omega_0, 1.3)]),(t,0,10));


# In[11]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[12]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# In[13]:


General_Solution_Real_Compact_func = sp.lambdify((t,A,phi,omega_0),General_Solution_Real_Compact)


# In[14]:


fig, ax = plt.subplots()
tdata = np.linspace(0,10,100)
line, = ax.plot([], [])
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)

def update(frame):
    ydata = General_Solution_Real_Compact_func(t=tdata,A=1,phi=frame,omega_0=1.3)
    ax.set_title("$\phi$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0, 2*np.pi, 100),interval = 30)
plt.show()


# In[15]:


fig, ax = plt.subplots()
tdata = np.linspace(0,10,100)
line, = ax.plot([], [])
ax.set_xlim(0, 10)
ax.set_ylim(-2, 2)

def update(frame):
    ydata = General_Solution_Real_Compact_func(t=tdata,A=frame,phi=0,omega_0=1.3)
    ax.set_title("$A$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(1, 2, 100),interval = 30)
plt.show()


# In[16]:


fig, ax = plt.subplots()
tdata = np.linspace(0,10,100)
line, = ax.plot([], [])
ax.set_xlim(0, 10)
ax.set_ylim(-1, 1)

def update(frame):
    ydata = General_Solution_Real_Compact_func(t=tdata,A=1,phi=0,omega_0=frame)
    ax.set_title("$\omega_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0.3, np.pi, 100),interval = 30)
plt.show()


# In[ ]:




