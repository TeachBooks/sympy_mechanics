#!/usr/bin/env python
# coding: utf-8

# # Sympy voorbeeld 1

# ![image.png](attachment:image.png)

# In[1]:


import sympy as sp


# In[2]:


Ah, Av, Bv, R, q, L = sp.symbols('Ah Av Bv R q L')


# De resultante van de q-last:

# In[3]:


R = q * L


# Som van de horizontale krachten is nul:

# In[4]:


eq1 = sp.Eq(Ah,0)


# Som van de verticale krachten is nul:

# In[5]:


eq2 = sp.Eq(Av+Bv-R,0)


# Som van de momenten om punt A is nul (linksom positief aangenomen):

# In[6]:


eq3 = sp.Eq(-R/2*L+Bv*L,0)


# Laat sympy de onbekende oplegreacties berekenen (3 onbedenken en 3 vergelijkingen, dus dit moet lukken):

# In[7]:


sol = sp.solve((eq1,eq2,eq3),(Ah,Av,Bv))
Ah_sol = sol[Ah]
Av_sol = sol[Av]
Bv_sol = sol[Bv]
display(Ah_sol)
display(Av_sol)
display(Bv_sol)


# Substitueer de waardes van Q en L:

# In[8]:


Ah_sub = Ah_sol.subs([(q,5),(L,10)])
Av_sub = Av_sol.subs([(q,5),(L,10)])
Bv_sub = Bv_sol.subs([(q,5),(L,10)])


# Print de uitkomsten van de oplegreacties met de juiste waardes van q en L:

# In[9]:


display(Ah_sub)
display(Av_sub)
display(Bv_sub)

