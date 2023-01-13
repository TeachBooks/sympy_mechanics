#!/usr/bin/env python
# coding: utf-8

# # DV-SB-ligger #

# In[1]:


import sympy as sp


# In[2]:


w = sp.symbols('w', cls=sp.Function)
q0, x = sp.symbols('q0 x')
L, EI = sp.symbols('L EI')
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')


# In[3]:


DV = sp.Eq(EI*sp.diff(w(x),x,4),q0) 
display(DV)


# In[4]:


w = sp.dsolve(DV, w(x)) 
w = w.rhs 
display(w)


# In[5]:


phi = -sp.diff(w, x)
kappa = sp.diff(phi, x)
M = EI * kappa
V = sp.diff(M, x)


# In[6]:


Eq1 = sp.Eq(w.subs(x, 0), 0) 
Eq2 = sp.Eq(w.subs(x, L), 0)
Eq3 = sp.Eq(M.subs(x, 0), 0)
Eq4 = sp.Eq(M.subs(x, L), 0)


# In[7]:


sol = sp.solve((Eq1,Eq2,Eq3,Eq4),(C1,C2,C3,C4))
display(sol)


# In[8]:


w_sol = w.subs(sol)
phi_sol = phi.subs(sol)
M_sol = M.subs(sol)
V_sol = V.subs(sol)
display(phi_sol.subs(x,0))
display(phi_sol.subs(x,L))
display(w_sol.subs(x,L/2))


# In[9]:


w_subs = w_sol.subs([(EI,1000),(q0,5),(L,8)])
phi_subs = phi_sol.subs([(EI,1000),(q0,5),(L,8)])
M_subs = M_sol.subs([(EI,1000),(q0,5),(L,8)])
V_subs = V_sol.subs([(EI,1000),(q0,5),(L,8)])


# In[10]:


sp.plot(-w_subs,(x,0,8),title='$w$');
sp.plot(-phi_subs,(x,0,8),title='$\phi$');
sp.plot(-M_subs,(x,0,8),title='$M$');
sp.plot(-V_subs,(x,0,8),title='$V$');

