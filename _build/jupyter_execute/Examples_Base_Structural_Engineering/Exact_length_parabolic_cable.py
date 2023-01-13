#!/usr/bin/env python
# coding: utf-8

# # Parabolic cable #

# In[1]:


import sympy as sp


# In[2]:


q, L, H = sp.symbols('q L H',positive=True,real=True,nonzero=True)
C1, C2 = sp.symbols('C1 C2')
x = sp.symbols('x',real=True,positive=True,nonzero=True)
z = sp.symbols('z',cls=sp.Function)


# In[3]:


ODE = sp.Eq(H * sp.diff(z(x),x,2), -q)


# In[4]:


z = sp.dsolve(ODE,z(x)).rhs


# In[5]:


eq1 = sp.Eq(z.subs(x,0),0)
eq2 = sp.Eq(z.subs(x,L),0)
sol = sp.solve((eq1,eq2),(C1,C2))
z = z.subs(sol)
display(z)


# Length based on Taylor approximation

# In[6]:


Ltaylor = sp.integrate(sp.nsimplify(1+1/2*sp.diff(z,x)**2),(x,0,L))
display(sp.simplify(Ltaylor))


# Length based on archlength

# In[7]:


display(z)
Lexact = sp.integrate(sp.sqrt(1+sp.diff(z,x)**2),(x,0,L))
display(Lexact) #cannot evaluate integral, sympy needs some manual help here? dz/dx = sinh z'
Lexact = sp.integrate(sp.sqrt(1+sp.diff(z.subs([(q,5),(L,10),(H,60)]),x)**2),(x,0,10))
display(Lexact)


# In[8]:


display(Ltaylor.subs([(q,5),(L,10),(H,60)]).evalf())
display(Lexact.subs([(q,5),(L,10),(H,60)]).evalf())


# In[ ]:




