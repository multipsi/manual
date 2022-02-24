#!/usr/bin/env python
# coding: utf-8

# # Configuration interaction
# 
# Multipsi provides a very general configuration code, allowing a wide variety of CI expansions, including FullCI and truncated CI (for example CI with single and double excitations CISD). Note however that the code is mostly designed around Full CI (for use in particular in complete active space CAS methods, see following chapters) and is thus not as efficient as a dedicated CIS or CISD code.

# ## Full CI
# 
# Let's start with a simple FullCI of water.
# 
# MultiPsi relies on Veloxchem to define the molecule and basis set. For more informations about how to run Veloxchem, please go to the Veloxchem webpage.
# 
# > <https://veloxchem.org>
# 
# Here we will start from a simple Hartree-Fock calculation in Veloxchem:

# In[1]:


import veloxchem as vlx

water_str='''
O   0.0   0.0   0.0
H   0.0   1.4   1.1
H   0.0  -1.4   1.1
'''
water=vlx.Molecule.read_str(water_str)
basis = vlx.MolecularBasis.read(water,"6-31G")
scfdrv = vlx.ScfRestrictedDriver()
scfdrv.compute(water, basis)


# We then define an orbital space object to store the starting orbitals and define the type of CI expansion:

# In[2]:


import multipsi as mtp

space=mtp.OrbSpace(water,scfdrv.mol_orbs)
print(space)


# By default, OrbSpace defines a restricted Hartree-Fock wavefunction (or restricted open-shell if the molecule is not a singlet). To define the full CI, we simply use the command FCI().

# In[3]:


space.FCI()
print(space)


# As we can see, internally the code assumes a complete active space, but since all orbitals are included in the active space, it is indeed a full CI.
# 
# We can also choose to freeze the core orbitals, as the correlation for core electrons tend to be constant and thus cancel out in relative energies.

# In[4]:


space.FCI(nfrozen=1)
print(space)


# This reduces significantly the number of determinants and thus the calculation time.
# 
# We can now run the configuration interaction:

# In[5]:


CIdrv=mtp.CIDriver(water,basis,space)
CIdrv.compute(1)


# If we are interested in several states and not just the ground state, we simply provide the number of states to the compute function:

# In[6]:


CIdrv.compute(3) #We want 3 states now


# At the end of the calculation, the CI Driver prints a few informations, but the key results can also be accessed directly in python.

# In[7]:


print("Ground state energy:",CIdrv.getEnergy(0)) #Get the energy of a specific state, here state number 0
print("Array of all state energies:",CIdrv.getEnergies()) #Get a python list of energies
print()
print("Ground state natural orbitals occupation numbers:",CIdrv.getNON(0))


# ## Truncated CI

# To define a truncated CI, one simply needs to define the truncated expansion in the OrbSpace object.
# 
# CIS and CISD expansions have a specific command, and higher order truncations can be done by using CI(n) with n the excitation order.

# In[8]:


space=mtp.OrbSpace(water,scfdrv.mol_orbs)
space.CISD(nfrozen=1) #Equivalent to space.CI(2,nfrozen=1)

CIdrv.compute(1)


# In[ ]:




