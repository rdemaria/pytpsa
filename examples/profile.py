import cProfile

cProfile.run('from compose import *;slow()',sort='time')
cProfile.run('from compose import *;fast()',sort='time')
