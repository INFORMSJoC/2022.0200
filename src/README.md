This folder consists of three Python scripts:

**QCQP.py** is the main script that runs the robust quadratically constrained quadratic program. The file includes:
- functions necessary for efficient setting up and manipulating the data used in such problems
- the proper script which takes as command line arguments:
  - m (meaning m=3 uncertain constraint and an objective, for running unconstrained problems choose 1)
  - n=10 (dimension of variable x)
  - l=10 (number of rows in linear transformation of x)
  - k=10 (the dimension of the z_i in each constraint i)
  - instance number: used to initiate the seed that generates data instance
  - max_time: maximum number of seconds to run the algorithms


  In the latter part of the script, all the data and relevant functions (data transformers, gradient computers etc.) are initialized and wrapped into a dict Data which serves as the main vehicle of carrying all components of the problem.
  
  Later, we have the loop that runs the first-order methods we compare ourselves to, where the steering parameters are the booleans i_FOP and i_FOD, which indicated whether one uses first-order algorithms on the primal or dual side of the problem. 
  All the possible combinations of these parameters are ran. 

  At the end of the script, we have the call to our SGSP algorithm to solve the same problem.

**PrimalDualOCO.py** which stores a single-function implementation of the methods we compare ourselves to.

**SGSP.py** which stores the implementation of our SGSP algorithm.
