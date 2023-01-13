#!/usr/bin/env python
# coding: utf-8

# # CIE_4140_Lecture_3_4_Python #

# ## Integral Fourier Transform and frequency-domain analysis ##

# In[1]:


import sympy as sp


# In[2]:


tau= sp.symbols('tau',real=True,positive = True)
t = sp.symbols('t',real=True)


# Example Function of Time

# In[3]:


f = sp.Heaviside(t-tau)-sp.Heaviside(t-2*tau)


# In[4]:


sp.plot(f.subs(tau,2),(t,0,10));


# In[5]:


omega = sp.symbols('omega',real=True,positive=True)
phi = sp.symbols('phi',real=True,positive=True)
f_omega = sp.integrate(f*sp.exp(-sp.I*omega*t), (t,-sp.oo,sp.oo))
f_phi = sp.fourier_transform(f,t,phi) 
f_omega = sp.simplify(f_omega)
f_phi = sp.simplify(f_phi)
display(f_omega)
display(f_phi)


# Another definition for the fourier transformation is used by sympy

# In[6]:


p0 = sp.plotting.plot(sp.re(f_omega.subs([(tau,2)])),(omega,-10,10),label='Real part' ,legend=True,show=False,adaptive=False,nb_of_points=3000)
p1 = sp.plotting.plot(sp.im(f_omega.subs([(tau,2)])),(omega,-10,10),label='Imaginary part',legend=True,show=False,adaptive=False,nb_of_points=3000)
p0.append(p1[0])
p0.show()
p0 = sp.plotting.plot(sp.re(f_phi.subs([(tau,2)])),(phi,-10/2/sp.pi,10/2/sp.pi),label='Real part' ,legend=True,show=False,adaptive=False,nb_of_points=3000)
p1 = sp.plotting.plot(sp.im(f_phi.subs([(tau,2)])),(phi,-10/2/sp.pi,10/2/sp.pi),label='Imaginary part',legend=True,show=False,adaptive=False,nb_of_points=3000)
p0.append(p1[0])
p0.show()


# In[7]:


x = sp.symbols('x',cls = sp.Function)
m, c, k = sp.symbols('m, c, k',real=True,positive=True)
Equation_of_Motion = m*sp.diff(x(t),t,2)+c*sp.diff(x(t),t)+k*x(t)-f
display(Equation_of_Motion)


# In[8]:


sp.fourier_transform(Equation_of_Motion,t,phi)


# So Homogeneous equation of motion cannot be Fourier transformed

# In[131]:


Equation_of_Motion_in_frequency_domain_homogeneous = sp.fourier_transform(m*sp.diff(x(t),t,2)+c*sp.diff(x(t),t)+k*x(t),t,phi)
Equation_of_Motion_in_frequency_domain = Equation_of_Motion_in_frequency_domain_homogeneous - f_phi
display(Equation_of_Motion_in_frequency_domain)


# Derivatives of ${d \over {dt}}x\left( t \right)$ not evaluated

# In[132]:


Equation_of_Motion_in_frequency_domain = Equation_of_Motion_in_frequency_domain.subs(sp.fourier_transform(x(t).diff(t),t,phi),-sp.I*c*2*sp.pi*phi*sp.fourier_transform(x(t),t,phi))
Equation_of_Motion_in_frequency_domain = Equation_of_Motion_in_frequency_domain.subs(sp.fourier_transform(x(t).diff(t,2),t,phi),(-sp.I*c*2*sp.pi*phi)**2*sp.fourier_transform(x(t),t,phi))
display(Equation_of_Motion_in_frequency_domain)


# In[133]:


solution_in_frequency_domain = sp.solve(sp.Eq(Equation_of_Motion_in_frequency_domain,0),sp.FourierTransform(x(t), t, phi))[0]
display(solution_in_frequency_domain)


# In[134]:


p0 = sp.plotting.plot(sp.re(solution_in_frequency_domain.subs([(tau,2),(k,4),(m,1),(c,0.75)])),(phi,0.001,15/2/sp.pi),label='Real part' ,legend=True,show=False,adaptive=False,nb_of_points=3000)
p1 = sp.plotting.plot(sp.im(solution_in_frequency_domain.subs([(tau,2),(k,4),(m,1),(c,0.75)])),(phi,0.001,15/2/sp.pi),label='Imaginary part',legend=True,show=False,adaptive=False,nb_of_points=3000)
p0.append(p1[0])
p0.show()


# In[135]:


solution = sp.inverse_fourier_transform(solution_in_frequency_domain, phi,t)
display(solution)


# Takes quite a while

# In[115]:


solution_in_frequency_domain = solution
solution_in_frequency_domain_numeric = solution_in_frequency_domain.subs([(tau,2),(k,4),(m,1),(c,3/4)])
display(solution_in_frequency_domain_numeric)


# In[117]:


sp.plot(solution_in_frequency_domain_numeric,(t,0,15))


# In[119]:


solution2 = sp.inverse_fourier_transform(sp.re(solution_in_frequency_domain), phi,t)
display(solution2)


# In[ ]:




