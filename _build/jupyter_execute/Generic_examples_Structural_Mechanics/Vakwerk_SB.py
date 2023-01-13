#!/usr/bin/env python
# coding: utf-8

# # Vakwerk_SB #

# Stelsel vergelijkingen (knoopevenwicht van een vakwerk)

# In[1]:


import sympy as sp
Ah, Av, Bv, N1, N2, N3, N4, N5, F = sp.symbols('Ah Av Bv N1 N2 N3 N4 N5 F')


# In[2]:


eq1 = sp.Eq(sp.nsimplify(Ah+N1+(1/sp.sqrt(5))*N2),0)
eq2 = sp.Eq(sp.nsimplify(Av+(2/sp.sqrt(5))*N2),0)
eq3 = sp.Eq(sp.nsimplify(-N1-(1/sp.sqrt(5))*N3+(1/sp.sqrt(5))*N4),0)
eq4 = sp.Eq(sp.nsimplify(Bv+(2/sp.sqrt(5))*N3+(2/sp.sqrt(5))*N4),0)
eq5 = sp.Eq(sp.nsimplify((-1/sp.sqrt(5))*N2+(1/sp.sqrt(5))*N3+N5),0)
eq6 = sp.Eq(sp.nsimplify((-2/sp.sqrt(5))*N2-(2/sp.sqrt(5))*N3),0)
eq7 = sp.Eq(sp.nsimplify((-1/sp.sqrt(5))*N4-N5+(3/5)*F),0)
eq8 = sp.Eq(sp.nsimplify((-2/sp.sqrt(5))*N4+(4/5)*F),0)
sol = sp.solve((eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8),(N1,N2,N3,N4,N5,Ah,Av,Bv))
display(sol)


# In[3]:


display(sol[N1].subs(F,50))
display(sol[N2].subs(F,50))
display(sol[N3].subs(F,50))
display(sol[N4].subs(F,50))
display(sol[N5].subs(F,50))
display(sol[Ah].subs(F,50))
display(sol[Av].subs(F,50))
display(sol[Bv].subs(F,50))

