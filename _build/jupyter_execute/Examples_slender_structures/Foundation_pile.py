#!/usr/bin/env python
# coding: utf-8

# # Foundation pile #

# In[1]:


import sympy as sp
u = sp.symbols('u', cls=sp.Function)
C1, C2 = sp.symbols('C1 C2')
x = sp.symbols('x')

b = 1/4
h = 1/4
Acircum = 2*b + 2*h
A = b*h
L = 20
EA = sp.nsimplify(5e9)
tau = 3e4
Fo = 2500e3
c = 16e8
q0 = sp.nsimplify(-tau*Acircum)
k = A * c
diffeq = sp.Eq(EA*sp.diff(u(x),x,2),-q0)

u = sp.dsolve(diffeq,u(x))
u = u.rhs
display(u)

N = EA * sp.diff(u,x)

eq1  = sp.Eq(N.subs(x , 0) , -Fo)
eq2  = sp.Eq(N.subs(x , 0) , -k*u.subs(x,0))

sol = sp.solve((eq1,eq2),
               (C1 ,C2 ))
u_sol = u.subs(sol)
N_sol = N.subs(sol)

display(u_sol)
display(N_sol)
sp.plot(u_sol,(x,0,L))
sp.plot(N_sol,(x,0,L));

