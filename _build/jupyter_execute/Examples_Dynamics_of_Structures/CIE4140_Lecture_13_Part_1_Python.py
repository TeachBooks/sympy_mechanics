#!/usr/bin/env python
# coding: utf-8

# # CIE4140_Lecture_13_Part_1_Python #

# ## Consider the steady-state vibrations of a string subjected to distributed viscous damping under a harmonic force ##

# In[1]:


import sympy as sp


# We represent the harmonic time-dependence by $e^{i \Omega t}$. The real part of the response to this load will give the response to the corresponding cosinusoidal load , whereas the imaginary part will give the response to the corresponding sinusoidal load

# Let us write the equation of motion

# In[2]:


w = sp.symbols('w',cls=sp.Function)
q1 = sp.symbols('q1',cls=sp.Function,real=True)
x,t = sp.symbols('x t',real=True)
nd, c, Omega= sp.symbols('nd c Omega',real=True)


# In[3]:


EQM = sp.diff(w(x,t),t,2) + 2 * nd * sp.diff(w(x,t),t) - c**2 * sp.diff(w(x,t),x,2) - q1(x)*sp.exp(sp.I*Omega*t)
display(EQM)


# We search for the steady-state solution in the form:

# In[4]:


W = sp.symbols('W',cls=sp.Function)
w_form = W(x) * sp.exp(sp.I * Omega * t)


# Substitution of this equation in the equation of motion gives:

# In[5]:


EQM_fr = sp.simplify(EQM.subs(w(x,t),w_form))
display(EQM_fr)


# $e^{i \Omega t}$ is not dropped, so can dropped manually, not necessary:

# In[6]:


display(EQM_fr.subs(sp.exp(sp.I * Omega * t),1))


# Let us find the soluton to this equation that satisfies the boundary conditions of a fixed-fixed string

# In[7]:


L = sp.symbols('L')
W_sol = sp.dsolve(EQM_fr,W(x),ics={W(0):0,W(L):0}).rhs
display(W_sol)


# We introduce numerical values for the parameters

# In[8]:


W_sol = W_sol.subs([(L,1),(c,2),(nd,1/2)])
display(W_sol)


# Let us consider two different shapes of the exernal force:

# In[9]:


q1_second_normal_mode = sp.sin(2 * sp.pi * x / 1)
q1_arbitrary = 14.39 * x **5 * (x - 1)
q1_constant = 1


# In[10]:


p0 = sp.plotting.plot(q1_second_normal_mode, (x, 0, 1),label='$q1_{second\ normal\ mode}$'    ,legend=True,show=False)
p1 = sp.plotting.plot(q1_arbitrary         , (x, 0, 1),label='$q1_{arbitrary}$' ,show=False)
p2 = sp.plotting.plot(q1_constant          , (x, 0, 1),label='$q1_{constant}$'  ,show=False)
p0.append(p1[0])
p0.append(p2[0])
p0.show()


# In[11]:


W_to_second_normal_mode = W_sol.subs(q1(x),q1_second_normal_mode).subs(q1(L),q1_second_normal_mode.subs(x,L)).subs(q1(0),q1_second_normal_mode.subs(x,0))
display(W_to_second_normal_mode)
W_to_arbitrary = W_sol.subs(q1(x),q1_arbitrary).subs(q1(L),q1_second_normal_mode.subs(x,L)).subs(q1(0),q1_second_normal_mode.subs(x,0))
display(W_to_arbitrary)
W_to_constant = W_sol.subs(q1(x),q1_constant).subs(q1(L),q1_second_normal_mode.subs(x,L)).subs(q1(0),q1_second_normal_mode.subs(x,0))
display(W_to_constant)


# Let us animate the reponses to these two loads assuming that the time-depandence is sinusoidal and introducing a numerical value of the load frequency: $\Omega=5$

# In[12]:


W_to_second_normal_mode_sin = sp.im(W_to_second_normal_mode*sp.exp(sp.I*Omega*t)).subs(Omega,5)
W_to_arbitrary_sin = sp.im(W_to_arbitrary*sp.exp(sp.I*Omega*t)).subs(Omega,5)
W_to_constant_sin = sp.im(W_to_constant*sp.exp(sp.I*Omega*t)).subs(Omega,5)


# In[13]:


display(W_to_second_normal_mode_sin.subs([(x,0.25),(t,3)]).evalf())


# In[14]:


W_to_second_normal_mode_sin_func = sp.lambdify((x,t),W_to_second_normal_mode_sin)
W_to_arbitrary_sin_func = sp.lambdify((x,t),W_to_arbitrary_sin)
W_to_constant_sin_func = sp.lambdify((x,t),W_to_constant_sin)


# Integral cannot be evaluated numerically, see `CIE4140_Lecture_13_Part_1_Python_seperate_q`

# Now we plot the amplitude-frequency response functions at two locations along the string for the two load shapes

# In[121]:


AFRF_second_normal_mode = sp.Abs(W_to_second_normal_mode.subs(x,3/4))
display(AFRF_second_normal_mode.evalf())


# AFRF can neither be evaluated

# In[ ]:




