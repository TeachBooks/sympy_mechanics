#!/usr/bin/env python
# coding: utf-8

# # VGM-q-last #

# In[1]:


import sympy as sp


# In[2]:


w1, w2, w3 = sp.symbols('w1 w2 w3', cls=sp.Function)
q1, q2, q3, x = sp.symbols('q1 q2 q3 x')
dq = sp.symbols('dq')
a, b, c, EI = sp.symbols('a b c EI')
C1, C2, C3, C4, C5, C6, C7, C8, C9, C10, C11, C12 = sp.symbols('C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12')


# In[3]:


b = a
c = a
q2 = q1 + dq * (x-a)/b
L = a + b + c


# In[4]:


DV1 = sp.Eq(EI*sp.diff(w1(x),x,4),q1)
DV2 = sp.Eq(EI*sp.diff(w2(x),x,4),q2)
DV3 = sp.Eq(EI*sp.diff(w3(x),x,4),q3)


# In[5]:


w1, w2, w3 = sp.dsolve([DV1, DV2, DV3],[w1(x),w2(x),w3(x)])
w1 = w1.rhs
w2 = w2.rhs
w3 = w3.rhs
display(w1)
display(w2)
display(w3)


# In[6]:


phi1 = -sp.diff(w1, x)
kappa1 = sp.diff(phi1, x)
M1 = EI * kappa1
V1 = sp.diff(M1, x)
phi2 = -sp.diff(w2, x)
kappa2 = sp.diff(phi2, x)
M2 = EI * kappa2
V2 = sp.diff(M2, x)
phi3 = -sp.diff(w3, x)
kappa3 = sp.diff(phi3, x)
M3 = EI * kappa3
V3 = sp.diff(M3, x)


# In[7]:


Eq1 =  sp.Eq(w1.subs(x, 0), 0) 
Eq2 =  sp.Eq(M1.subs(x, 0), 0)
Eq3 =  sp.Eq(  w1.subs(x, a),   w2.subs(x, a))
Eq4 =  sp.Eq(  V1.subs(x, a),   V2.subs(x, a))
Eq5 =  sp.Eq(phi1.subs(x, a), phi2.subs(x, a)) 
Eq6 =  sp.Eq(  M1.subs(x, a),   M2.subs(x, a)) 
Eq7 =  sp.Eq(  w2.subs(x, a + b),   w3.subs(x, a + b))
Eq8 =  sp.Eq(  M2.subs(x, a + b),   M3.subs(x, a + b))
Eq9 =  sp.Eq(phi2.subs(x, a + b), phi3.subs(x, a + b))
Eq10 = sp.Eq(  V2.subs(x, a + b),   V3.subs(x, a + b))
Eq11 = sp.Eq(  w3.subs(x, a + b + c), 0)
Eq12 = sp.Eq(  M3.subs(x, a + b + c), 0)


# In[8]:


sol = sp.solve((Eq1,Eq2,Eq3,Eq4,Eq5,Eq6,Eq7,Eq8,Eq9,Eq10,Eq11,Eq12),(C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12))
display(sol)


# In[9]:


w1_sol = w1.subs(sol)
phi1_sol = phi1.subs(sol)
M1_sol = M1.subs(sol)
V1_sol = V1.subs(sol)
w2_sol = w2.subs(sol)
phi2_sol = phi2.subs(sol)
M2_sol = M2.subs(sol)
V2_sol = V2.subs(sol)
w3_sol = w3.subs(sol)
phi3_sol = phi3.subs(sol)
M3_sol = M3.subs(sol)
V3_sol = V3.subs(sol)
display(phi1_sol.subs(x,0))
display(phi3_sol.subs(x,L))


# example with some numbers: activate a line by removing the # and de-activate a line by inserting a # at the begin

# In[10]:


numset = {dq : 2 , a : 2, q1 : 2, q3 : 4,  EI: 5000} 
w1_subs = w1_sol.subs(numset)
phi1_subs = phi1_sol.subs(numset)
M1_subs = M1_sol.subs(numset)
V1_subs = V1_sol.subs(numset)
w2_subs = w2_sol.subs(numset)
phi2_subs = phi2_sol.subs(numset)
M2_subs = M2_sol.subs(numset)
V2_subs = V2_sol.subs(numset)
w3_subs = w3_sol.subs(numset)
phi3_subs = phi3_sol.subs(numset)
M3_subs = M3_sol.subs(numset)
V3_subs = V3_sol.subs(numset)


# In[11]:


sp.plot((-2,(x,0,2)),(-2 - 2 * (x-2)/2,(x,2,4)),(-4,(x,4,6)),title='$q$');


# In[12]:


sp.plot((-w1_subs,(x,0,2)),(-w2_subs,(x,2,4)),(-w3_subs,(x,4,6)),title='$w$')
sp.plot((-phi1_subs,(x,0,2)),(-phi2_subs,(x,2,4)),(-phi3_subs,(x,4,6)),title='$\phi$')
sp.plot((-M1_subs,(x,0,2)),(-M2_subs,(x,2,4)),(-M3_subs,(x,4,6)),title='$M$')
sp.plot((-V1_subs,(x,0,2)),(-V2_subs,(x,2,4)),(-V3_subs,(x,4,6)),title='$V$');


# In[13]:


display(phi1_subs.subs(x,0).evalf())
display(phi3_subs.subs(x,6).evalf())

