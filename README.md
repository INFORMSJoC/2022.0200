[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# CacheTest

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data
that were used in the research reported on in the paper 
First-order algorithms for robust optimization problems via
convex-concave saddle-point Lagrangian reformulation (https://doi.org/10.1287/ijoc.2022.0200) by K. Postek and S. Shtern. 
The snapshot is based on 
(https://github.com/tkralphs/JoCTemplate/commit/f7f30c63adbcb0811e5a133e1def696b74f3ba15) 
in the development repository. 

**Important: This code is being developed on an on-going basis at 
https://github.com/tkralphs/JoCTemplate. Please go there if you would like to
get a more recent version or would like support**

## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2022.0200

https://doi.org/10.1287/ijoc.2022.0200.cd

Below is the BibTex for citing this snapshot of the respoitory.

```
@article{CacheTest,
  author =        {K. Postek and S. Shtern},
  publisher =     {INFORMS Journal on Computing},
  title =         {First-order algorithms for robust optimization problems via convex-concave saddle-point Lagrangian reformulation},
  year =          {2022},
  doi =           {10.1287/ijoc.2022.0200.cd},
  url =           {https://github.com/INFORMSJoC/2022.0200},
}  
```

## Description

The goal of this software is to generate QCQP problems with ball uncertainty and compare the performance of the approach suggested in the paper using the Lagrangian formulation with the existing approach based on online convex optimization.

## Running

In Linux, the experiment can be conducted using the script Run_QCQP_script.sh from the script folder

```
chmod +x Run_QCQP_script.sh
./Run_QCQP_script.sh
```

Note that script will run a total of 50 realizations in batches of 10.
The parameter used in the script (inputed in this order) are
m+1=4 (meaning 3 uncertain constraint and an objective, for running unconstrained problems choose 1)
n=10 (dimension of variable x)
l=10 (number of rows in linear transformation of x)
k=10 (the dimension of the z_i in each constraint i)
seed - calculated automatically
Time limit=600  

The scripts runs the Python script QCQP.py located in the src folder. This script prepares the data and elementary functions and calls both SGSP and OCO with several parameter combination for solving the robuts problem.


## Results

Figure 1 in the paper shows the results of the multiplication test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/mult-test.png)

Figure 2 in the paper shows the results of the sum test with different
values of K using `gcc` 7.5 on an Ubuntu Linux box.

![Figure 1](results/sum-test.png)

## Replicating

To replicate the results in the paper the Run_QCQP_script.sh script should be changed as follows:

The small examples in the paper are ran with 
```
cpulimit -l 100 -i python3 ../src/QCQP.py 1 10 10 10 $seed 600&>output$seed.txt&
cpulimit -l 100 -i python3 ../src/QCQP.py 4 10 10 10 $seed 600&>output$seed.txt&

```
The medium examples are ran with
```
cpulimit -l 100 -i python3 ../src/QCQP.py 1 600 15 25 $seed 1200&>output$seed.txt&
cpulimit -l 100 -i python3 ../src/QCQP.py 4 600 15 25 $seed 1200&>output$seed.txt&
```

The medium examples are ran with
```
cpulimit -l 100 -i python3 ../src/QCQP.py 1 3600 16 30 $seed 3600&>output$seed.txt&
cpulimit -l 100 -i python3 ../src/QCQP.py 4 3600 16 30 $seed 3600&>output$seed.txt&
```

To replicate the results in [Figure 1](results/mult-test), do either

```
make mult-test
```
or
```
python test.py mult
```
To replicate the results in [Figure 2](results/sum-test), do either

```
make sum-test
```
or
```
python test.py sum
```
