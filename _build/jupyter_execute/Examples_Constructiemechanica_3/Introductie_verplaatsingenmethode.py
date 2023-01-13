#!/usr/bin/env python
# coding: utf-8

# # Introductie verplaatsingenmethode

# ![image.png](attachment:image.png)

# In[1]:


import sympy as sp


# In[2]:


a, b, c, k1, k2, k3 = sp.symbols('a b c k1 k2 k3')


# In[3]:


Bg = sp.Matrix([[0,-1,-a],[0,-1,b],[-1,0,c]])
BgT = sp.transpose(Bg)
Dg = sp.Matrix([[k1,0,0],[0,k2,0],[0,0,k3]])
K = BgT*Dg*Bg
display(K)
F = K.inv()
display(F)
load = sp.Matrix([50,150,-5])
disp = F*load
disp_sol = disp.subs([(a,3),(b,2),(c,1),(k1,1000),(k2,2000),(k3,3000)])
display(disp_sol)
display(disp_sol.evalf())

