#!/usr/bin/env python
# coding: utf-8

# # CIE4140_Lecture_13_Part_2_Python #

# ## Let us calculate the natural and resonanse frequencies of a string with arbitrary (though linear) boundary conditions ##

# The universal recipe is to substitute $w(x,t)=W(x) e^{i \omega t})

# Strictly speaking, the natural frequencies are defined for undamped systems only. Therefore, we assume that there is no damping along the string nor in the boundaries.

# In[1]:


import sympy as sp


# The equation of motion in this case reads

# In[2]:


w = sp.symbols('w',cls=sp.Function,real=True)
x,t = sp.symbols('x,t',real=True)
c = sp.symbols('c',real=True)
EQM = sp.diff(w(x,t),t,2)-c**2*sp.diff(w(x,t),x,2)
display(EQM)


# Substituting $w(x,t)=W(x) e^{i \omega t}$, we obtain the following ordinary differential equation, in which $\omega$ is one of the natural frequencies
# 

# In[3]:


W = sp.symbols('W',cls=sp.Function,real=True)
omega = sp.symbols('omega',real=True)
EQM_Fr = sp.simplify(EQM.subs(w(x,t),W(x)*sp.exp(sp.I*omega*t)))
display(EQM_Fr)


# The general solution to this equation reads

# In[4]:


C1, C2 = sp.symbols('C1, C2',real=True)
W_sol = sp.dsolve(EQM_Fr,W(x)).rhs
display(W_sol)
W_sol = C1*sp.sin(omega*x/c)+C2*sp.cos(omega*x/c)
display(W_sol)


# Now we will consider different boundary conditions

# 1: Fixed-fixed string: $W(0)=0$, $W(L)=L$

# In[5]:


L = sp.symbols('L',real=True)
BC_left =  W_sol.subs(x,0)
BC_right = W_sol.subs(x,L)
display(BC_left)
display(BC_right)


# When $\omega$ equals one of the natural frequencies, the determinant of the coefficient matrix composed of the two boundary conditions (i.e., the frequency equation) should be zero. 

# Let us compose the coefficient matrix

# In[6]:


Coeff_Matrix= sp.Matrix([[BC_left.diff(C1),BC_left.diff(C2)],
                         [BC_right.diff(C1),BC_right.diff(C2)]])
display(Coeff_Matrix)


# In[7]:


Frequency_Equation = sp.det(Coeff_Matrix)
display(Frequency_Equation)


# Now we search for the roots of the frequency equation, which, by definition, are the natural frequencies (in the undamped case these frequencies coincide with the resonance frequencies) 

# Often, especially in the case of more complex boundary conditions, we introduce $\Omega_{dl}= {\Omega \over L} C$. For this dimensionless parameter, the roots of the frequency equation can be always found numerically 

# In[8]:


Omega_dl = sp.symbols('Omega_dl',real=True)
Frequency_Equation_DL = Frequency_Equation.subs(omega,Omega_dl/L*c)
display(Frequency_Equation_DL)


# Obviously, independently of the wave speed and the string length $\Omega_{dl}=\pi n$, where $n$ is an integer
# 

# In[9]:


sol = sp.solveset(sp.Eq(Frequency_Equation_DL,0),Omega_dl)
display(sol)


# Leta us find values of $\beta$ also numerically

# In[10]:


Dimensionless_Natural_Frequency = []
step = 0.01
for i in range(5000):
    xi = i * step
    left_value = Frequency_Equation_DL.subs(Omega_dl,xi)
    right_value = Frequency_Equation_DL.subs(Omega_dl,xi+step)
    if left_value*right_value < 0:
        Dimensionless_Natural_Frequency.append(xi+step/2)
        
display(Dimensionless_Natural_Frequency)


# 2: Fixed-free string: $W(0)=0, {\left. {{{dW\left( x \right)} \over {dx}}} \right|_{x = L}} = 0$

# In[11]:


BC_left =  W_sol.subs(x,0)
BC_right = W_sol.diff(x).subs(x,L)
display(BC_left)
display(BC_right)


# In[12]:


Coeff_Matrix= sp.Matrix([[BC_left.diff(C1),BC_left.diff(C2)],
                         [BC_right.diff(C1),BC_right.diff(C2)]])
display(Coeff_Matrix)


# In[13]:


Frequency_Equation = sp.det(Coeff_Matrix)
display(Frequency_Equation)


# In[14]:


Frequency_Equation_DL = Frequency_Equation.subs(omega,Omega_dl/L*c)
display(Frequency_Equation_DL)


# In[15]:


sol = sp.solveset(sp.Eq(Frequency_Equation_DL,0),Omega_dl)
display(sol)


# 3: Fixed-{Elastically Supported Inertial End}: 
# - ${w\left( {0,t} \right) = 0}$
# - ${\left. {T{{\delta w\left( {x,t} \right)} \over {\delta x}}} \right|_{x = L}} =  - m{\left. {{{{\delta ^2}w\left( {x,t} \right)} \over {\delta {x^2}}}} \right|_{x = L}} + kw\left( {x,t} \right)$

# In[16]:


T, m, k = sp.symbols('T, m, k',real=True)
BC_left =  W_sol.subs(x,0)
BC_right = T*sp.diff(w(x,t),x)+m*sp.diff(w(x,t),t,2)+k*w(x,t)
display(BC_left)
display(BC_right)


# In[17]:


BC_right = sp.simplify(BC_right.subs(w(x,t),W_sol*sp.exp(sp.I*omega*t))).subs(x,L)
display(BC_right)


# $e^{i \Omega t}$ is not dropped, so can dropped manually, not necessary:

# In[18]:


BC_right = BC_right.subs(sp.exp(sp.I*omega*t),1)
display(BC_right)


# In[19]:


Coeff_Matrix= sp.Matrix([[BC_left.diff(C1),BC_left.diff(C2)],
                         [BC_right.diff(C1),BC_right.diff(C2)]])
display(Coeff_Matrix)


# In[20]:


Frequency_Equation = sp.det(Coeff_Matrix)
display(Frequency_Equation)


# This equation should be solved numerically. The natural frequencies depend both on the string parameters and the mass and stiffness of the boundary

# In[21]:


Frequency_Equation_subs = Frequency_Equation.subs([(L,1),(c,2),(T,10),(k,20),(m,7)])


# In[22]:


sol = sp.solveset(sp.Eq(Frequency_Equation_subs,0),omega)
display(sol)


# In[23]:


Dimensionless_Natural_Frequency = []
step = 0.01
for i in range(5000):
    xi = i * step
    left_value = Frequency_Equation_subs.subs(omega,xi)
    right_value = Frequency_Equation_subs.subs(omega,xi+step)
    if left_value*right_value < 0:
        Dimensionless_Natural_Frequency.append(xi+step/2)
        
display(Dimensionless_Natural_Frequency)


# In[ ]:




