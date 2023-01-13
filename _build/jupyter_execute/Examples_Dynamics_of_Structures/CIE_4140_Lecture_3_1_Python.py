#!/usr/bin/env python
# coding: utf-8

# # CIE_4140_Lecture_3_1_Python #

# ## Free vibration of an SDOF with viscous damping ##

# In[1]:


import sympy as sp


# In[2]:


u = sp.symbols('u',cls=sp.Function)
t, omega_n, m, zeta = sp.symbols('t, omega_n, m, zeta',real=True,positive=True)
u_0, v_0 = sp.symbols('u_0, v_0')


# The equation of motion:

# In[3]:


Equation_of_Motion= sp.Eq(sp.diff(u(t),t,2)+2*zeta*omega_n*sp.diff(u(t),t)+omega_n **2 * u(t),0)
display(Equation_of_Motion)


# In[4]:


u_sol = sp.dsolve(Equation_of_Motion, u(t), ics={u(0): u_0 , u(t).diff(t , 1).subs(t , 0) : v_0}).rhs
display(u_sol)


# ## Super-critically damed system $1 < \xi $ ##
# 

# Example plot

# In[5]:


sp.plot(u_sol.subs([(u_0,1),(v_0,1),(omega_n,1),(zeta,1.2)]),(t,0,12));


# Effect of the initial displacement:

# In[6]:


get_ipython().run_line_magic('matplotlib', 'notebook')


# In[7]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# In[8]:


u_func = sp.lambdify((u_0,v_0,omega_n,zeta,t),u_sol)


# In[9]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_func(u_0=frame,v_0=1,omega_n=1,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $u_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of the initial velocity:

# In[10]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_func(u_0=1,v_0=frame,omega_n=1,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $v_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of the undamped natural frequency:

# In[11]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_func(u_0=1,v_0=0.1,omega_n=frame,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $\omega_n$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0.5,2,100),interval = 100)
plt.show()


# Effect of the damping ratio:

# In[12]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_func(u_0=1,v_0=0.1,omega_n=1,zeta=frame,t=tdata)
    ax.set_title("Displacement versus time for $\zeta$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(1.1,2,100),interval = 100)
plt.show()


# ## Critically damed system $\xi = 1 $ ##
# 

# In[13]:


u_crit = sp.limit(u_sol,zeta,1)
display(u_crit)


# Example plot

# In[14]:


sp.plot(u_crit.subs([(u_0,1),(v_0,1),(omega_n,1),(zeta,1.2)]),(t,0,12));


# Effect of the initial displacement

# In[15]:


u_crit_func = sp.lambdify((u_0,v_0,omega_n,zeta,t),u_crit)


# In[16]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_crit_func(u_0=frame,v_0=1,omega_n=1,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $u_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of the initial velocity

# In[17]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_crit_func(u_0=1,v_0=frame,omega_n=1,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $v_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of undamped natural frequency

# In[18]:


fig, ax = plt.subplots()
tdata = np.linspace(0,12,100)
line, = ax.plot([], [])
ax.set_xlim(0, 12)
ax.set_ylim(0, 1.2)

def update(frame):
    ydata = u_crit_func(u_0=1,v_0=0.1,omega_n=frame,zeta=1.2,t=tdata)
    ax.set_title("Displacement versus time for $\omega_n$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0.5,2,100),interval = 100)
plt.show()


# ## Sub-critically damped system $\xi < 1 $ ##

# In[19]:


u_sol_sub = sp.exp(-zeta*omega_n*t)*(u_0*sp.cos(omega_n**sp.sqrt(1-zeta**2)*t)+(v_0+zeta*omega_n*u_0)/(omega_n*sp.sqrt(1-zeta**2))*sp.sin(omega_n*sp.sqrt(1-zeta**2)*t))
display(u_sol_sub)


# Example plot

# In[20]:


sp.plot(u_sol_sub.subs([(u_0,1),(v_0,1),(omega_n,1),(zeta,0.05)]),(t,0,42));


# Effect of the initial displacement:

# In[21]:


u_sub_func = sp.lambdify((u_0,v_0,omega_n,zeta,t),u_sol_sub)


# In[22]:


fig, ax = plt.subplots()
tdata = np.linspace(0,42,500)
line, = ax.plot([], [])
ax.set_xlim(0, 42)
ax.set_ylim(-1.5, 1.5)

def update(frame):
    ydata = u_sub_func(u_0=frame,v_0=1,omega_n=1,zeta=0.05,t=tdata)
    ax.set_title("Displacement versus time for $u_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of the initial velocity:

# In[23]:


fig, ax = plt.subplots()
tdata = np.linspace(0,42,500)
line, = ax.plot([], [])
ax.set_xlim(0, 42)
ax.set_ylim(-1.5, 1.5)

def update(frame):
    ydata = u_sub_func(u_0=1,v_0=frame,omega_n=1,zeta=0.05,t=tdata)
    ax.set_title("Displacement versus time for $v_0$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0,1,100),interval = 100)
plt.show()


# Effect of the undamped natural frequency:

# In[24]:


fig, ax = plt.subplots()
tdata = np.linspace(0,42,500)
line, = ax.plot([], [])
ax.set_xlim(0, 42)
ax.set_ylim(-2.2,2.2)

def update(frame):
    ydata = u_sub_func(u_0=1,v_0=0.1,omega_n=frame,zeta=0.05,t=tdata)
    ax.set_title("Displacement versus time for $\omega_n$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0.5,2,100),interval = 100)
plt.show()


# Effect of the damping ratio

# In[25]:


fig, ax = plt.subplots()
tdata = np.linspace(0,42,500)
line, = ax.plot([], [])
ax.set_xlim(0, 42)
ax.set_ylim(-1.5,1.5)

def update(frame):
    ydata = u_sub_func(u_0=1,v_0=1,omega_n=1,zeta=frame,t=tdata)
    ax.set_title("Displacement versus time for $\zeta$ = "+str(np.round(frame,2)))
    line.set_data(tdata, ydata)

ani = FuncAnimation(fig, update, frames=np.linspace(0.01,0.9,100),interval = 100)
plt.show()


# In[ ]:




