#!/usr/bin/env python
# coding: utf-8

# # Example Extension week 1 strain energy#

# In[1]:


import sympy as sp


# In[2]:


q, L, x, F, N, EA = sp.symbols('q L x F N EA')


# In[3]:


N = q * (L-x) + F
eps = N  / EA
Ps = sp.integrate(EA*eps**2/2,(x,0,L))
display(sp.simplify(Ps))
Ps = sp.integrate(1/EA*N**2/2,(x,0,L))
display(sp.simplify(Ps))
sp.plot(N.subs([(L,10),(F,25e3),(q,5e3),(EA,2500e3)]),(x,0,10),title="Normal force $N(x)$ in N")
display(Ps.subs([(L,10),(F,25e3),(q,5e3),(EA,2500e3)]).evalf())


# In[ ]:




