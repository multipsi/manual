#!/usr/bin/env python
# coding: utf-8

# # MCSCF
# The core functionality of MultiPsi is multiconfigurational methods, with multiconfigurational self-consistent field (MCSCF) as its cornerstone.
# 
# ## CASSCF
# The most commonly used form of MCSCF nowadays is the complete active space self-consistent field (CASSCF), which is a form of MCSCF where the configurations are generated from a full CI expansion within an active space consisting of a few orbitals. Because of the factorial growth of full CI with the number of active orbitals, the active space should usually be limited to less than 20 electrons in 20 orbitals, with the current record being a CAS(22,22) on the triplet state of hydrogenase. However, this large calculation was only possible on a large cluster using the entire distributed memory available.
# 
# To define a CASSCF, we simply need to define our OrbSpace as a full CI expansion (see previous chapters) and provide it to our MCSCF module. For example, we can start with a simple Hartree-Fock, here for furan:

# In[1]:


import veloxchem as vlx
import multipsi as mtp
furan_xyz="""9

 C     -0.86213    -0.90784     0.00007
 H     -1.63433    -1.64264    -0.00003
 C      0.50727    -0.90524     0.00007
 C      0.92057     0.47886    -0.00003
 C     -0.22323     1.23186    -0.00003
 O     -1.35123     0.40376    -0.00013
 H      1.17117    -1.74724     0.00017
 H      1.93767     0.81866     0.00007
 H     -0.46573     2.26986    -0.00013
"""

molecule = vlx.Molecule.from_xyz_string(furan_xyz)
basis = vlx.MolecularBasis.read(molecule,"def2-sv(p)")

scfdrv = vlx.ScfRestrictedDriver()
scfdrv.compute(molecule, basis)


# Then we choose all $\pi$ orbitals (so 5 orbitals total with 6 electrons):

# In[2]:


space=mtp.OrbSpace(molecule,scfdrv.mol_orbs)
space.CAS(6,5)

mcscfdrv=mtp.McscfDriver(molecule,basis,space)
mcscfdrv.compute(1)


# You see the convergence is very poor. The reason for it is that we started with the wrong active orbitals and the MCSCF had to make up for this. Let's visualize the resulting active orbitals:

# In[3]:


space.act2cube(basis)

import py3Dmol

viewer = py3Dmol.view(viewergrid=(2,1),width=800,height=400)

#for iorb in range(0):
#    filename="cube_"+str(iorb+1)+".cube"
#    with open(filename, "r") as f: cube = f.read()
#    # Plot strick structures
#    viewer.addModel(cube, "cube",viewer=(0,iorb)); viewer.setStyle({'stick':{}},viewer=(0,iorb))
#    # Negative and positive nodes
#    #viewer.addVolumetricData(cube, "cube", {"isoval": -0.02, "color": "blue", "opacity": 0.75},viewer=(0,iorb))
    #viewer.addVolumetricData(cube, "cube", {"isoval":  0.02, "color":  "red", "opacity": 0.75},viewer=(0,iorb))

filename="cube_1.cube"
with open(filename, "r") as f: cube = f.read()
# Plot strick structures
viewer.addModel(cube, "cube",viewer=(1,0)); #viewer.setStyle({'stick':{}},viewer=(0,0))

cube2=cube
#filename="cube_1.cube"
#with open(filename, "r") as f: cube = f.read()
# Plot strick structures
viewer.addModel(cube2, "cube",viewer=(0,0)); viewer.setStyle({'stick':{}},viewer=(0,0))

'''
iorb=0
filename="cube_"+str(iorb+1)+".cube"
with open(filename, "r") as f: cube = f.read()
# Plot strick structures
viewer.addModel(cube, "cube",viewer=(0,iorb)); viewer.setStyle({'stick':{}},viewer=(0,iorb))
# Negative and positive nodes
#viewer.addVolumetricData(cube, "cube", {"isoval": -0.02, "color": "blue", "opacity": 0.75},viewer=(0,iorb))
#viewer.addVolumetricData(cube, "cube", {"isoval":  0.02, "color":  "red", "opacity": 0.75},viewer=(0,iorb))
'''
viewer.show()


# We can see that we got 3 $\pi$ orbitals correctly, but 2 were not. By default, when only specifying the number of electrons and orbitals like we did here, the program chooses those around the HOMO-LUMO gap. While it is sometimes fine, most often the orbitals we want are not those, and we thus need to resort to visual inspection.
# 
# In this case, it is simple to just look at the MO coefficients. Indeed, the $\pi$ orbitals have a different symmetry than the $\sigma$ and thus show in the coefficients. Here, the p0 orbitals are those forming the $\pi$ system.

# In[4]:


space.molorb.print_orbitals(molecule,basis,True)


# We can thus recognize that the orbitals we want are 12, 17, 18, 19 and 21. We can specify them using the CASOrb keyword, which both defines a CAS expansion and allows to select the active orbitals:

# In[5]:


space.CASOrb([11,16,17,18,20])
mcscfdrv=mtp.McscfDriver(molecule,basis,space)
mcscfdrv.compute(1)


# The calculation converged much faster, and you can also notice that the occupation numbers are a bit further away from 2 and 0, which is common of $\pi$ bonds. Let us visualize them to confirm:

# In[7]:


space.act2cube(basis)

viewer = py3Dmol.view(viewergrid=(1,1),width=800,height=400)

filename="cube_1.cube"
with open(filename, "r") as f: cube = f.read()
# Plot strick structures
viewer.addModel(cube, "cube",viewer=(0,0)); viewer.setStyle({'stick':{}},viewer=(0,0))
# Negative and positive nodes
viewer.addVolumetricData(cube, "cube", {"isoval": -0.02, "color": "blue", "opacity": 0.75},viewer=(0,0))
viewer.addVolumetricData(cube, "cube", {"isoval":  0.02, "color":  "red", "opacity": 0.75},viewer=(0,0))

viewer.show()


# It is always recommended to inspect the orbitals after your calculations, as the optimisation may sometimes swap our active orbitals, especially if their occupations are close to 2 or 0.

# ## State-specific and state-averaged CASSCF
# 
# When one is interested in more than just the ground state, for example in spectroscopy or photochemistry, we can specify this in the mcscf driver:

# In[8]:


nstates=3
mcscfdrv.compute(nstates)


# By default, the program then uses a state-averaged optimisation. This means the orbitals are found to minimize the average of all the selected states (here 3). By doing so, we obtain a compromise and all states are treated at a similar accuracy, favouring energy cancellation in the energy differences (i.e. excitation energies). On the other hand, the orbitals are not optimal for any specific state, and here we can see for example that the ground state energy is slightly higher than the one found in the previous section.
# 
# Alternatively, one can decide to optimize the orbitals for a specific state, or biasing the average for a specific state. For this, we use the state-average optimizer but with custom weights instead of the default equal weights:

# In[9]:


mcscfdrv.compute(nstates, [0,1,0])


# Here we selected 0 weights on all states except the second one, and as expected, that state is now lower in energy than in the state-averaged case.
# 
# Computing each state this way (called the state-specific optimization) provides the best description of each state but suffers from 2 major issues. First, it is not always easy to optimize excited states this way, as the order of the states may change depending on the orbitals and it is thus difficult to truly follow a specific state. Additionally, the orthonormality of the states wavefunction is not enforced, and thus the excited states may appear lower in energy than they should be. This can be fixed a posteriori by computing the overlap and Hamiltonian between the different states and 
