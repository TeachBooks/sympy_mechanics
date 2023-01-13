#!/usr/bin/env python
# coding: utf-8

# # Beam with varying EI and q #

# In[1]:


import sympy as sp
w = sp.symbols('w', cls=sp.Function)
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
x = sp.symbols('x')
L = 10
EI = sp.nsimplify(1000)
q = sp.nsimplify(-500*w(x)+10+400*sp.diff(w(x),x,2))

sp.plotting.plot(EI,(x,0,L))
#sp.plotting.plot(q,(x,0,L))

diffeq = sp.Eq(EI*sp.diff(w(x),x,4),q)
display(diffeq)

w = sp.dsolve(diffeq)
w = w.rhs
display(w)

phi = -sp.diff(w, x)
kappa = sp.diff(phi, x)
M = EI * kappa
V = sp.diff(M, x)

eq1  = sp.Eq(w.subs(x , 0) , 0)
eq2  = sp.Eq(M.subs(x , 0) , 0)
eq3  = sp.Eq(w.subs(x , L) , 0)
eq4  = sp.Eq(phi.subs(x , L) , 0)

sol = sp.solve((eq1,eq2,eq3,eq4) ,
               (C1 ,C2 ,C3 ,C4))
w_sol = w.subs(sol)
M_sol = M.subs(sol)

sp.plot(w_sol,(x,0,L))
#sp.plot(M_sol,(x,0,L));


# In[ ]:




