#!/usr/bin/env python
# coding: utf-8

# # Exercise 6 Work and Energy #

# In[1]:


import sympy as sp


# In[2]:


rho, EI, L, q, x, Fv, Q = sp.symbols('rho EI L q x Fv Q')
numset = {EI : 10000, L : 8, q : 5, Q : 0}
k = rho * EI / L**3
M = q * (L-x)**2 / 2 - (Q-Fv)*(L-x)
Ev1 = sp.integrate(M**2/(2*EI),(x,0,L))
Ev2 = Fv**2 / (2 * k)
w1 = sp.diff(Ev1, Q)
w2 = sp.diff(Ev2, Fv)
eq1 = sp.Eq(w1.subs(Q,0),w2.subs(Q,0))
Fv_sol = sp.solve(eq1,Fv)[0]
M_sol = M.subs(Fv,Fv_sol)
V = sp.diff(M_sol,x)
xmax = sp.solve(V,x)[0]
display(sp.simplify(xmax.subs(numset)))
Mmax = M_sol.subs(x,xmax)
display(Mmax.subs(numset))
Mi = M_sol.subs(x,0)
Eq2 = sp.Eq(Mi,-3*Mmax)
rho_sol = sp.solve(Eq2,rho)[1]
display(rho_sol.subs(numset))
Mi_sol = Mi.subs(rho,rho_sol)
Mmax_sol = Mmax.subs(rho,rho_sol)
display(sp.simplify(Mi_sol))
display(sp.simplify(3*Mmax_sol))
Mi_sub = Mi_sol.subs(numset)
Mmax_sub = Mmax_sol.subs(numset)
display(Mi_sub.evalf())
display(Mmax_sub.evalf())
M_sub = M_sol.subs(rho,rho_sol).subs(numset)
display(M_sub)
sp.plot(-M_sub,(x,0,8));


# In[ ]:





# In[ ]:




