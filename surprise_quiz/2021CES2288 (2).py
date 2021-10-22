#!/usr/bin/env python
# coding: utf-8

# In[7]:


import math
import numpy as np
E=210000000
v=0.3
t=0.025
w=3000


# In[4]:


xi,yi = (0,0)
xj,yj = (0.5,0.25)
xm,ym = (0,0.25)

xp,yp = (0,0)
xq,yq = (0.5,0)
xr,yr = (0.5,0.25)

A1 = (xi*(yj-ym)+xj*(ym-yi)+xm*(yi-yj))/2
A2 = (xp*(yq-yr)+xq*(yr-yp)+xr*(yp-yq))/2


# In[8]:


B1 = (np.array([yj-ym,0,ym-yi,0,yi-yj,0,0,xm-xj,0,xi-xm,0,xj-xi,xm-xj,yj-ym,xi-xm,ym-yi,xj-xi,yi-yj]).reshape((3,6)))/(2*A1)
B2 = (np.array([yq-yr,0,yr-yp,0,yp-yq,0,0,xr-xq,0,xp-xr,0,xq-xp,xr-xq,yq-yr,xp-xr,yr-yp,xq-xp,yp-yq]).reshape((3,6)))/(2*A2)


# In[9]:


D = (E/(1-v**2)) * (np.array([1,v,0,v,1,0,0,0,(1-v)/2]).reshape((3,3)))


# In[10]:


B1_T = B1.transpose()
k1 = t*A1*B1_T
k11=np.matmul(k1,D)
k111=np.matmul(k11,B1)

B2_T = B2.transpose()
k2 = t*A2*B2_T
k22=np.matmul(k2,D)
k222=np.matmul(k22,B2)


# In[11]:


tsm = np.zeros((8,8))
tsm[0:6,0:6] += k111
tsm[]


# In[ ]:




