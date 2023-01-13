#!/usr/bin/env python
# coding: utf-8

# # Variation of load cable #

# In[1]:


import sympy as sp


# In[2]:


q, p, H, dH, q1, q2, L = sp.symbols('q p H dH q1 q2 L',positive=True,real=True)
x = sp.symbols('x')
C1, C2, C3, C4 = sp.symbols('C1 C2 C3 C4')
q1 = q + p
q2 = q - p


# In[3]:


z1 = sp.nsimplify((-1/(H+dH))*((1/2)*q1*x**2+C1+C2*x))
z2 = sp.nsimplify((-1/(H+dH))*((1/2)*q2*x**2+C3+C4*x))


# In[4]:


V1 = (H+dH)*sp.diff(z1,x)
V2 = (H+dH)*sp.diff(z2,x)
eq1= sp.Eq(z1.subs(x,0),0)
eq2= sp.Eq(z1.subs(x,L/2),z2.subs(x,L/2))
eq3= sp.Eq(V1.subs(x,L/2),V2.subs(x,L/2))
eq4= sp.Eq(z2.subs(x,L),0)
sol = sp.solve((eq1,eq2,eq3,eq4),(C1,C2,C3,C4))
display(sol)


# In[5]:


z1 = z1.subs(sol)
z2 = z2.subs(sol)
display(z1)
display(z2)


# In[6]:


eq5 = sp.Eq(sp.integrate(1+sp.diff(z1,x)**2/2,(x,0,L/2)) + sp.integrate(1+sp.diff(z2,x)**2/2,(x,L/2,L)),L+q**2*L**3/(24*H**2))
dH_sol = sp.solve(eq5,dH)[0]
dH_sol


# In[7]:


beta=dH_sol/H
sp.simplify(beta)


# In[8]:


beta.subs(p,1/4*q)


# In[9]:


zk1 = z1
pp1 = zk1.subs([(x,L/4),(dH,dH_sol),(p,q/4)])/((3/32)*q*L**2/H)
sp.simplify(pp1.evalf())


# In[ ]:




