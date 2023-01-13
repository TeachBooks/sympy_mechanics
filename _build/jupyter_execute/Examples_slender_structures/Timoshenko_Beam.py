#!/usr/bin/env python
# coding: utf-8

# # Timoshenko beam #

# In[1]:


import sympy as sp
w, phi = sp.symbols('w phi', cls=sp.Function)
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
x, EI, GA, L, q = sp.symbols('x EI GA L q')

diffeq1 = sp.Eq(EI*sp.diff(phi(x),x,2)-GA*(sp.diff(w(x),x,1)+phi(x)),0)
diffeq2 = sp.Eq(GA*(sp.diff(w(x),x,2)+sp.diff(phi(x),x,1)),-q)
display(diffeq1)
display(diffeq2)

w,phi = sp.dsolve([diffeq1, diffeq2],[w(x),phi(x)])
w = w.rhs
phi = phi.rhs


# In[2]:


gamma_shear = phi + sp.diff(w, x)
kappa = sp.diff(phi, x)
M = EI * kappa
V = GA*gamma_shear

eq1  = sp.Eq(w.subs(x , 0) , 0)
eq2  = sp.Eq(M.subs(x , 0) , 0)
eq3  = sp.Eq(w.subs(x , L) , 0)
eq4  = sp.Eq(M.subs(x , L) , 0)

sol = sp.solve((eq1,eq2,eq3,eq4) ,
               (C1 ,C2 ,C3 ,C4))
w_sol = w.subs(sol)
M_sol = M.subs(sol)

display(sp.simplify(w_sol.subs(x,L/2)))

