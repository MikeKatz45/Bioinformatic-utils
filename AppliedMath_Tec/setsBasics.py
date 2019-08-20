#!/usr/bin/env python
# coding: utf-8

# In[15]:


emptyset = set() #Initialize empty set (just braces defaults to dictionary)
numberset = {1, 3, 2} #Initialize a set of numbers
stringset = {'hola', 'hi', 'konichiwa', 'hi'} #Initialize set of strings

#elements may print in a different order compared to imput because order
#in sets is not relevant. Duplicates are removed automatically
print(emptyset)
print(numberset)
print(stringset)


# In[12]:


A = {1, 2, 4}
A.add(3) #add elements to a set
print(A)


# In[14]:


B = {'bat', 'raquet', 'ball'}
B.discard('bat') #remove elements from a set
print(B)

B.clear() #remove all elements from a set
print(B)


# In[18]:


A = {2, 1, 3}
B = B = {'bat', 'raquet', 'ball'}

#Order/sorte ascendingly
print(sorted(A))
print(sorted(B))

#Order/sort descendingly
print(sorted(A, reverse=True))
print(sorted(B, reverse=True))


# In[21]:


A = {1, 2, 3}
B = {4, 5, 6}

A.union(B) #perform union of two sets
A | B # alternative


# In[23]:


A = {1, 2, 3}
B = {3, 4, 5}

A.intersection(B) #perform intersection of two sets
A & B # alternative


# In[25]:


A = {1, 2}
B = {2}

A.difference(B) #perform difference of two sets
A - B # alternative


# In[26]:


A = {1, 2, 3}
B = {3, 4, 5}

A.symmetric_difference(B) #perform symmetric difference of two sets
A ^ B # alternative


# In[27]:


A = {'swim', 'run', 'play'}
b = 'run'

b in A #membership test of element in a set


# In[28]:


A = {'swim', 'run', 'play'}
B = {'swim', 'play'}

B.issubset(A) #test of subset relationship between 2 sets

