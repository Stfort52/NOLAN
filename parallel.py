import multiprocessing as mp
import numpy as np
from functools import reduce
from operator import mul
import math

class SharedArray:
    def __init__(self, dtype, shape):
        self.shape = shape
        self.nptype = np.dtype(dtype)
        ctype = np.ctypeslib.as_ctypes_type(self.nptype)
        self.shm = mp.RawArray(ctype, reduce(mul, shape))

    def to_numpy(self):
        return np.frombuffer(self.shm, dtype=self.nptype).reshape(self.shape)

    def to_shm(self):
        return self.shm
        
    def dims(self):
        return self.shape
    
    def type(self):
        return self.nptype
    
    
    ## someday add subscription support lol

def roundDown2sqrt(x: int) -> int:
    """
    Round down to the nearest square root.
    """
    return int(math.floor(math.sqrt(x)))

def splitRanges(x: int, n: int) -> list[tuple]:
    """
    Split `range(x)` into n tuple containing new ranges with balanced length.
    """
    if x <= n:
        for i in range(x):
            yield i, i + 1
    else:
        k, m = divmod(x, n)
        for i in range(n):
            yield i * k + min(i, m), (i + 1) * k + min(i + 1, m)