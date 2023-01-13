#!/usr/bin/env python
# coding: utf-8

# # CIE_4140_Lecture_18_Part_1_Python #

# ## Excitation of a rod by a point pulse-load $P_0 Dirac(x) Dirac(t)$ ##

# In[1]:


import sympy as sp


# In[2]:


u = sp.symbols('u',cls=sp.Function,real=True)
x, t, k= sp.symbols('x, t, k',real=True)
c, P_0 = sp.symbols('c P_0',real=True)


# In[3]:


EM = sp.diff(u(x,t),t,2)-c**2*sp.diff(u(x,t),x,2)-P_0 * sp.DiracDelta(x)*sp.DiracDelta(t)
display(EM)


# In[4]:


EM_fourier = sp.fourier_transform(EM,x,k)
display(EM_fourier)


# In[5]:


U = sp.symbols('U',cls=sp.Function)
EM_fourier = EM_fourier.subs(sp.FourierTransform(sp.Derivative(u(x, t), (x, 2)), x, k),-k**2*U(t))
EM_fourier = EM_fourier.subs(sp.FourierTransform(sp.Derivative(u(x, t), (t, 2)), x, k),sp.Derivative(U(t), (t, 2)))
display(EM_fourier)


# In[6]:


C1, C2 = sp.symbols('C1, C2')
Solution_fourier = sp.dsolve(EM_fourier,U(t))
display(Solution_fourier)


# Why diracdelta not taken into account?

# In[7]:


Solution_fourier = Solution_fourier.subs([(C1,0),(C2,0)]).rhs
display(sp.simplify(Solution_fourier))


# In[8]:


Solution = sp.integrate(Solution_fourier.subs([(P_0,1),(c,1)]),(k,0,sp.oo))/sp.pi
display(Solution)


# In[ ]:




