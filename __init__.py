#Copyright (c) 2008, Riccardo De Maria
#All rights reserved.

from polmap import *
from funset import *
from dop import *
#from lintrack import *
#from lie import *


if __name__=='__main__':
  import doctest
  mods='polmap pol'.split()
  for mod in mods:
      m=__import__(mod)
      doctest.testmod(m)
