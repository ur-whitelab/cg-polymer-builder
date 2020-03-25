# User variables

bond_length = 0.05
scaling_method = 'grid'
degree_of_polymerization = 10
N_peg = 5
box_size = 5

# Calculations

MW_polymer = 133.11*degree_of_polymerization+44.05*N_peg

def density_to_n(density, molecular_weight, box_dimensions):
  """
  Description: output the number of molecules of molecular_weight and density
  that are needed to fill a box with size box_dimensions
  density: density of the molecules in g/ml
  molecular_weight: molecular weight of the full in amu or g/mol
  box_dimensions: dimensions of box to fill in nm
  """
  box_volume = np.prod(box_dimensions)
  N = density / molecular_weight * (1e-21) * (6.02e23) * box_volume
  return int(N)

box_dimensions = [box_size, box_size, box_size]
# N_water = density_to_n(1, 18.02, box_dimensions)
# N_polymer = density_to_n(25e-3, MW_polymer, box_dimensions)
N_water = 200
N_polymer = 10

# Start of script

import numpy as np
import mbuild as mb
import random

class Monomer(mb.Compound):
  """
  Create a monomer that is read from a file.
  filepath: path of file on system that describes the monomer.
  bond_indices: indices of atoms that other monomers will connect to.
  """
  def __init__(self, filepath, bond_indices):
    super(Monomer, self).__init__()
    mb.load(filepath, compound=self)
    # mb.Compound.translate(self, -self[0].pos)  # center molecule at origin
    # self.add_bond((self[1], self[3]))
    self.add(mb.Port(anchor=self[bond_indices[0]]), label='up')
    self.add(mb.Port(anchor=self[bond_indices[1]]), label='down')
    mb.Compound.translate(self['up'], [bond_length, 0, 0])
    mb.Compound.translate(self['down'], [-bond_length, 0, 0])

class Polymer(mb.Compound):
  """
  Create a custom polymer block.
  monomer_chain: a list containing the monomers the block should comprise of in the specified order
  """
  def __init__(self, monomer_chain):
    super(Polymer, self).__init__()
    #
    for i in range(len(monomer_chain)):
      next_monomer = mb.clone(monomer_chain[i])
      if i == 0:
        # Initialize first monomer if no monomers exist yet
        current_monomer = next_monomer
      else:
        # Move monomer over from up port to down port
        mb.force_overlap(move_this=next_monomer,
                         from_positions=next_monomer['up'],
                         to_positions=current_monomer['down'])
      self.add(next_monomer)
      current_monomer = next_monomer
      mb.Compound.translate(self[0], [0, 0, 0]) # center molecule at origin

def PolymerBlock(monomers, monomer_amounts, number_repeated = 1):
  """
  Description: creates and returns a custom-built polymer block.
  monomers: a list containing the monomers the block should comprise of in the correct order.
  monomer_amounts: a list containing the number of times the respective monomers in monomers will be added.
  number_repeated: the number of times the above arrangement should be repeated. Default 1
  Example: To create a polymer A-A-A-A-B-B-B-A-A-C-C-C-C-D, we use monomers=[A,B,A,C,D] and monomer_list=[4,3,2,4,1].
  """
  if type(monomers) != list or type(monomer_amounts) != list:
    raise TypeError('Your inputs must be lists. Please see the docstring for more information.')
  if len(monomers) != len(monomer_amounts):
    raise ValueError('Please ensure your input lists are the same size.')
  monomer_names = []
  for i in range(number_repeated): # repeat block number_repeated times
    for j in range(len(monomers)): # loop through all monomers
      for k in range(monomer_amounts[j]): # keep adding monomer_amounts of that monomer to the chain
        monomer_names.append(monomers[j])
  return Polymer(monomer_names)

def RandomizedPolymerBlock(monomers, monomer_ratios, chain_length):
  """
  Description: creates and returns polymer block with a random arrangement of monomers.
  monomers: a list containing the monomers which the polymer should be made out of
  monomer_ratios: a list comprised of the respective abundances of the monomers
  chain_length: the number of monomers the chain should be made of
  Example: To create a polymer with A and B in a 4:1 ratio, we use monomer = [A,B]
  and monomer_ratios = [4,1] or [0.4, 0.1] or [0.80, 0.20].
  """
  probability_ranges = []
  monomer_names = []
  # Make probability array to randomly choose a monomer using rand
  for i in range(len(monomer_ratios)):
    probability_ranges.append(monomer_ratios[i] * (1 / sum(monomer_ratios))) # calculate probability that ith monomer will show up
    if i != 0:
      probability_ranges[i] += probability_ranges[i-1] # keep adding probabilities to get list of cumulative probabilities
  # Build up randomized arrays of monomers according to ratio
  for i in range(chain_length):
    rand = random.random() # random number between 0 and 1
    num = min(i for i in probability_ranges if i > rand) # check which monomer this corresponds to by comparing to 'probability_ranges'
    pos = probability_ranges.index(num) # get index of respective monomer
    monomer_names.append(monomers[pos])
  return Polymer(monomer_names)

# Creation of polymers

# Monomers and molecules
alpha = Monomer('molecules/alpha_acetic_acid.pdb', [2, 0])
beta = Monomer('molecules/beta_acetic_acid.pdb', [2, 0])
peg = Monomer('molecules/PEG.pdb', [1, 0])
water = mb.load('molecules/water.pdb')

# Create polymers
polymer_block_A = RandomizedPolymerBlock([alpha, beta], [1, 3], degree_of_polymerization)
polymer_block_PEG = PolymerBlock([peg], [1], N_peg)
polymer = PolymerBlock([polymer_block_A, polymer_block_PEG], [1, 1])
# polymer = peg

# system = mb.packing.fill_box([polymer, water], box=box_dimensions, n_compounds=[N_polymer, N_water])
system = polymer_block_PEG
system.save('polymer_test.gsd', overwrite=True)
















# if scaling_method == 'grid':
#   polymer_positions = mb.Grid3DPattern(1, 8, 3) # [x, y, z]
#   polymer_positions.scale([3, 2, 1.5]) # [x, y, z]; initial (ratio) scaling
#   polymer_positions.scale(3) # absolute scaling
# else:
#   polymer_positions = mb.Random3DPattern(20)
#   polymer_positions.scale(3)

# polymers = mb.Compound()
# for position in polymer_positions:
#   current_polymer = mb.clone(polymer1)
#   mb.Compound.translate(current_polymer, position)
#   polymers.add(current_polymer)

# Solvate system

# solvate_box_extension = np.array([0, 0, 0])
# polymer_box = polymers.boundingbox.maxs
# solvate_box = polymer_box + solvate_box_extension
# solvate_box = np.around(solvate_box, decimals=1).tolist() # needs to be a python list

# solvate_box = polymers.boundingbox

# print(solvate_box)

# system = mb.packing.solvate(
#   polymers, water, 200, box = solvate_box,
#   # overlap=0.2, seed=12345, edge=0.2, fix_orientation=False,
#   # temp_file=None, update_port_locations=False
# )

