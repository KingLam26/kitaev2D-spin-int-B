from core import *
from utils import *
from datetime import datetime
import multiprocessing as mp
import sys, time


bp = [1,1,1,1,2,2,2,2]
bp1 = [1,1,1,1,3,3,2,2]
bp2 = [1,1,1,1,2,3,3,2]
bp3 = [1,1,1,1,2,3,2,3]

bp = [1,1,1,1,2,2,3,3]
bp1 = [0,0,1,1,2,2,3,3]
bp2 = [1,0,0,1,2,2,3,3]
bp3 = [0,1,0,1,2,2,3,3]


"""
bp = [1,1,1,2,1,2,2,2]
bp1 = [1,1,1,3,1,3,2,2]
bp2 = [1,1,1,2,1,3,3,2]
bp3 = [1,1,1,2,1,3,2,3]

bp = [1,1,2,1,1,2,2,2]
bp1 = [1,1,2,1,1,2,3,3]
bp2 = [1,1,2,1,1,3,3,2]
bp3 = [1,1,2,1,1,3,2,3]

bp = [1,2,1,2,1,2,1,2]
bp1 = [1,2,1,2,1,3,1,3]
bp2 = [1,2,1,3,1,3,1,2]
bp3 = [1,2,1,3,1,2,1,3]
"""


_, a = process_permute(bp)
print(abs(a))
_, a1 = process_permute(bp1)
print(a1)
_, a2 = process_permute(bp2)
print(a2)
_, a3 = process_permute(bp3)
print(a3)

b = a1+a2+a3
print(b*2 / abs(a))