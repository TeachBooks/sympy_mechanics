#!/usr/bin/env python
# coding: utf-8

# # Taylor series #

# In[1]:


import sympy as sp


# In[2]:


x = sp.symbols('x')


# In[3]:


f = sp.sqrt(1+x**2)


# In[4]:


display(f.series(x, n=5))


# In[ ]:




