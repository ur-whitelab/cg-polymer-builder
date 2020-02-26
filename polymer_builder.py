import numpy as np
import mbuild as mb
import random

# VARIABLES
bond_length = 0.05

class Monomer(mb.Compound):
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
  monomers: a list containing the monomers the block should comprise of in the correct order
  monomer_list: a list containing the number of times the respective monomers in monomers will be added
  Example: To create a polymer A-A-A-A-B-B-B-A-A-C-C-C-C, we use monomers=[A,B,A,C] and monomer_list=[4,3,2,4].
  No number_repeated value is necessary.
  """
  def __init__(self, monomer_chain):
    super(Polymer, self).__init__()
    #
    for i in range(len(monomer_chain)):
      next_monomer = mb.clone(monomer_chain[i])
      if i == 0:
        # Initialize the first monomer properly
        current_monomer = next_monomer
      else:
        # Move monomer over and add it to the chain
        mb.force_overlap(move_this=next_monomer,
                         from_positions=next_monomer['up'],
                         to_positions=current_monomer['down'])
      self.add(next_monomer)
      current_monomer = next_monomer
      mb.Compound.translate(self[0], [0, 0, 0]) # center molecule at origin

def PolymerBlock(monomers, monomer_amounts, number_repeated = 1):
  if type(monomers) != list or type(monomer_amounts) != list: raise TypeError('Your inputs must be lists. Please see the docstring for more information.')
  if len(monomers) != len(monomer_amounts): raise ValueError('Please ensure your input lists are the same size.')
  monomer_names = []
  for i in range(number_repeated): # repeat block number_repeated times
    for j in range(len(monomers)): # loop through all monomers
      for k in range(monomer_amounts[j]): # keep adding monomer_amounts of that monomer to the chain
        monomer_names.append(monomers[j])
  return Polymer(monomer_names)

def RandomizedPolymerBlock(monomers, monomer_ratios, chain_length):
  probability_ranges = []
  monomer_names = []
  # Make probability array to randomly choose a monomer using rand
  for i in range(len(monomer_ratios)):
    probability_ranges.append(monomer_ratios[i] * (1 / sum(monomer_ratios))) # calculate probability that ith monomer will show up
    if i != 0: probability_ranges[i] += probability_ranges[i-1] # keep adding probabilities to get list of cumulative probabilities
  # Build up randomized arrays of monomers according to ratio
  for i in range(chain_length):
    rand = random.random() # random number between 0 and 1
    num = min(i for i in probability_ranges if i > rand) # check which monomer this corresponds to by comparing to 'probability_ranges
    pos = probability_ranges.index(num) # get index of respective monomer
    monomer_names.append(monomers[pos])
  return Polymer(monomer_names)

# Creation of system

alpha = Monomer('alpha_acetic_acid.pdb', [2, 0])
beta = Monomer('beta_acetic_acid.pdb', [2, 0])

polymer = RandomizedPolymerBlock([alpha, beta], [1, 3], 20)
polymer2 = PolymerBlock([alpha, beta], [1, 3], 5)

# SYSTEM

# system = mb.Compound()

# pattern_disk = mb.DiskPattern(50)
# pattern_disk.scale(5)

# # now clone the polymer and move it to the points in the pattern
# for pos in pattern_disk:
#     current_polymer = mb.clone(polymer)
#     mb.translate(current_polymer, pos)
#     system.add(current_polymer)

polymer2.save('mbuild_test.gsd', overwrite=True)