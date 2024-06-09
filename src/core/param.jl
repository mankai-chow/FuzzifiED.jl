
"""
    SilentStd :: bool = false 

a flag to determine whether logs of the FuzzifiED functions should be turned off. False by default. If you want to evaluate without log, put `FuzzifiED.SilentStd = true`. This parameter can be defined for each process separately. 
"""
SilentStd = false

"""
    NumThreads :: Int = Threads.nthreads()

an integer to define how many threads OpenMP uses. By default, it is the same as the number of threads in Julia. If you use Jupyter notebooks, which by default uses one core only, you may need to define this by hand, _e.g._, `FuzzifiED.NumThreads = 8`. This parameter can be defined for each process separately. 
"""
NumThreads = Threads.nthreads()

"""
    Libpath :: Int = FuzzifiED_jll.LibpathFuzzifiED

define where the Fortran library are compiled. You do not need to modify that by yourself. However, if you compile the Fortran codes by yourself, you need to point this to your compiled library. 
"""
Libpath = FuzzifiED_jll.LibpathFuzzifiED