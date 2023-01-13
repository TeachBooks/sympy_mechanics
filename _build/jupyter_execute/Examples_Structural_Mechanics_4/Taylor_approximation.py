#!/usr/bin/env python
# coding: utf-8

# # Taylor approximation #

# In[1]:


import sympy as sp


# In[2]:


x = sp.symbols('x')


# In[3]:


z = 1 - sp.sqrt(1-x**2)


# In[4]:


display(z.series(x, n=3))

