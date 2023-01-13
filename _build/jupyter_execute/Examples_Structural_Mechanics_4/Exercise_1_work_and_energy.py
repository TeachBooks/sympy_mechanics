#!/usr/bin/env python
# coding: utf-8

# # Exercise 1 work and energy #

# In[1]:


import sympy as sp


# In[2]:


x, q, L = sp.symbols('x q L')


# In[3]:


M1 = 0
M2 = -q*x**2/2
M3 = -q*(L-x)**2/2
m1 = -x/2
m2 = -L/2 - x/2
m3 = -L + x


# In[4]:


w = sp.integrate(M1*m1,(x,0,L))+ sp.integrate(M2*m2,(x,0,L)) + sp.integrate(M3*m3,(x,0,L))
display(w)


# In[ ]:




