#!/usr/bin/env python
# coding: utf-8

# # Wrinkler element definition

# In[1]:


import sympy as sp
w, phi = sp.symbols('w phi', cls=sp.Function)
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
x, EI, L, q, beta = sp.symbols('x EI L q beta',positive=True)
w1, w2, phi1, phi2 = sp.symbols('w1 w2 phi1 phi2')

k = sp.S(4)*EI*beta**4

ODE = sp.Eq(EI*sp.diff(w(x),x,4)+k*w(x),q)
display(ODE)

w = sp.dsolve(ODE,w(x))
w = w.rhs

display(w)


# In[2]:


phi = -sp.diff(w, x)
kappa = sp.diff(phi, x)
M = EI * kappa
V = sp.diff(M, x)

eq1  = sp.Eq(w.subs(x , 0) , w1)
eq2  = sp.Eq(phi.subs(x , 0) , phi1)
eq3  = sp.Eq(w.subs(x , L) , w2)
eq4  = sp.Eq(phi.subs(x , L) , phi2)

display(eq1)
display(eq2)
display(eq3)
display(eq4)

sol = sp.solve((eq1,eq2,eq3,eq4) ,
               (C1 ,C2 ,C3 ,C4))
w_sol = w.subs(sol)
V_sol = V.subs(sol)
M_sol = M.subs(sol)

Fz1 = sp.expand(sp.simplify(-V_sol.subs(x,0)))
Fz2 = sp.expand(sp.simplify(V_sol.subs(x,L)))
Ty1 = sp.expand(sp.simplify(-M_sol.subs(x,0)))
Ty2 = sp.expand(sp.simplify(M_sol.subs(x,L)))


# In[ ]:





# In[50]:


k11 = sp.simplify(Fz1.coeff(w1))
k12 = sp.simplify(Fz1.coeff(phi1))
k13 = sp.simplify(Fz1.coeff(w2))
k14 = sp.simplify(Fz1.coeff(phi2))
k21 = sp.simplify(Ty1.coeff(w1))
k22 = sp.simplify(Ty1.coeff(phi1))
k23 = sp.simplify(Ty1.coeff(w2))
k24 = sp.simplify(Ty1.coeff(phi2))
k31 = sp.simplify(Fz2.coeff(w1))
k32 = sp.simplify(Fz2.coeff(phi1))
k33 = sp.simplify(Fz2.coeff(w2))
k34 = sp.simplify(Fz2.coeff(phi2))
k41 = sp.simplify(Ty2.coeff(w1))
k42 = sp.simplify(Ty2.coeff(phi1))
k43 = sp.simplify(Ty2.coeff(w2))
k44 = sp.simplify(Ty2.coeff(phi2))

Ksys = sp.Matrix([[k11,k12,k13,k14],
        [k21,k22,k23,k24],
        [k31,k32,k33,k34],
        [k41,k42,k43,k44]])
display(Ksys)


# In[48]:


f1 = sp.simplify(-Fz1.coeff(q))*q
f2 = sp.simplify(-Ty1.coeff(q))*q
f3 = sp.simplify(-Fz2.coeff(q))*q
f4 = sp.simplify(-Ty2.coeff(q))*q
Fsys = sp.Matrix([f1,f2,f3,f4])
display(Fsys)


# In[49]:


k11#.rewrite(sp.cos).simplify()


# In[60]:


K11=4*EI*beta**3*(sp.sin(2*beta*L)+sp.sinh(2*beta*L))/((sp.cos(2*beta*L)+sp.cosh(2*beta*L)-2))
display(K11)
check = K11-k11
check.subs([(L,2),(EI,500),(beta,sp.pi**2)]).evalf()


# In[ ]:





# In[ ]:




