from typing import List

def pseudoangle(v: vector):
  """
  Compute the pseudoangle of a vector in R^2
  """
  if v[0] == 0:
    if v[1] > 0:
      return float('inf')
    else:
      return float('-inf')
  return v[1]/v[0]

class CSPolygon:
  """
  A centrally symmetric polygon in R^2 represented by a list of segments (vectors)
  that are summed together using the Minkowski sum.
  """

  def __init__(self, vectors: List[vector]):
    self.vectors = []
    for v in vectors:
      assert v[0] != 0 or v[1] != 0, "Cannot have zero vector"
      if v[0] < 0:
        self.vectors.append(-v)
      else:
        self.vectors.append(v)
    self.vectors.sort(key=pseudoangle)

  def __repr__(self):
    return f"CSPolygon({self.vectors})"

  def __add__(self, rhs):
    return CSPolygon(self.vectors + rhs.vectors)
  
  def volume(self):
    vol = 0
    sweap = vector([0,0])
    for v in self.vectors:
      vol += abs(sweap[0]*v[1] - sweap[1]*v[0])
      sweap += v
    return vol

  def mixed_volume(self, rhs):
    return (self + rhs).volume() - self.volume() - rhs.volume()