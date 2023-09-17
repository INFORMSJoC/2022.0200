This folder consists of three Python scripts:

QCQP.py is the main script that runs the robust quadratically constrained quadratic program. The file include:
- functions necessary for efficient setting up and manipulating the data used in such problems
- the proper script which takes as command line arguments:
- - m+1=4: (meaning m=3 uncertain constraint and an objective, for running unconstrained problems choose 1)
- - n=10 (dimension of variable x)
- l=10 (number of rows in linear transformation of x)
- k=10 (the dimension of the z_i in each constraint i)
- seed - calculated automatically
- Time limit=600  
