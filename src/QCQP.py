import numpy as np
import scipy
import scipy.linalg
import scipy.stats
import gurobipy
from gurobipy import *
import casadi as ca
import pandas as pd
import time
import sys
from PrimalDualOCO import PrimalDualOCO
from SGSP import SGSP
#from Ronny import Ronny

# function for the worst-case z given x
def Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, i_constraint):

    # this function, given an x,
    # transforms the problems data into a form that
    # gives a quadratic function in terms of the uncertain parameter

    Px = np.matrix(np.zeros((l, K)))
    for k in range(0, K):
        Px[:, k] = np.matmul(Pjk[i_constraint][k], x)

    sth = np.matmul(Pjk[i_constraint][-1], x)
    rx = 2 * np.matmul(Px.T, sth)
    sx = cjk[i_constraint][-1] + np.matmul(bjk[i_constraint][-1].T, x) + np.linalg.norm(sth) ** 2
    Qx = np.matmul(Px.T, Px)

    for k in range(0, K):
        rx[k, 0] = rx[k, 0] + np.matmul(bjk[i_constraint][k].T, x) + cjk[i_constraint][k]

    return Qx, rx, sx

def Transform_to_Qz_rz_sz(z, Pjk, bjk, cjk, n, m, l, K, i_constraint):

    P_z = Pjk[i_constraint][-1]
    r_z = bjk[i_constraint][-1]
    s_z = cjk[i_constraint][-1]

    for k in range(K):
        P_z = P_z + np.multiply(Pjk[i_constraint][k], np.asscalar(z[k]))
        r_z = r_z + np.multiply(bjk[i_constraint][k], np.asscalar(z[k]))
        s_z = s_z + np.multiply(cjk[i_constraint][k], np.asscalar(z[k]))

    Q_z = np.matmul(P_z.T, P_z)#np.matmul(P_z.T, P_z)

    return Q_z, r_z, s_z


def Transform_to_Qzx_rz_sz(z, x,Pjk, bjk, cjk, n, m, l, K, i_constraint):

    P_z = Pjk[i_constraint][-1]
    r_z = bjk[i_constraint][-1]
    s_z = cjk[i_constraint][-1]

    for k in range(K):
        P_z = P_z + np.multiply(Pjk[i_constraint][k], np.asscalar(z[k]))
        r_z = r_z + np.multiply(bjk[i_constraint][k], np.asscalar(z[k]))
        s_z = s_z + np.multiply(cjk[i_constraint][k], np.asscalar(z[k]))

    Q_zx = np.matmul(P_z.T, np.matmul(P_z,x))#P_z.T @ (P_z @ x)

    return Q_zx, r_z, s_z
##### make it separable w.r.t. x and z

def ComputeGradients(x, z, Pjk, bjk, cjk, n, m, l, K):
    # this function computes gradients of all functions returning them as a list, the input
    # z should be a list of z's corresponding to different constraints

    gx = []
    gz = []

    for j in range(0, m):
        [Qzx, rz, _] = Transform_to_Qzx_rz_sz(z[j], x,Pjk, bjk, cjk, n, m, l, K, j)
        [Qx, rx, _] = Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, j)
        [lambda_max, v] = scipy.linalg.eigh(Qx, eigvals = tuple([K-1, K-1]))
        gradientz = 2 * np.matmul(Qx  - np.asscalar(lambda_max) * np.eye(K),z[j]) + rx
        gz.append(gradientz)

        # computing the term related to lambda_max
        P_x = np.matrix(np.zeros((l, n)))

        for k in range(K):
            P_x = P_x + np.asscalar(v[k]) * Pjk[j][k]

        gradientx = 2 * Qzx + rz + 2 * (1 - np.linalg.norm(z[j]) ** 2) * np.matmul(P_x.T, np.matmul(P_x, x))


        gx.append(gradientx)

    return gx, gz

def ComputeGradientX(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint):
    # this function computes gradients of all functions returning them as a list, the input
    # z is a single realization
    [Qzx, rz, _] = Transform_to_Qzx_rz_sz(z, x, Pjk, bjk, cjk, n, m, l, K, i_constraint)
    [Qx, _, _] = Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, i_constraint)
    [_, v] = scipy.linalg.eigh(Qx, eigvals = tuple([K-1, K-1]))

    # computing the term related to lambda_max
    P_x = np.matrix(np.zeros((l, n)))
    for k in range(K):
        P_x = P_x + np.asscalar(v[k]) * Pjk[i_constraint][k]

    gradientx = 2 * Qzx + rz + 2 * (1 - np.linalg.norm(z) ** 2) * np.matmul(P_x.T, np.matmul(P_x, x))

    return gradientx

def ComputeGradientZ(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint):
    # this function computes gradients of all functions returning them as a list, the input
    # z should be a list of z's corresponding to different constraints

    [Qx, rx, _] = Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, i_constraint)
    #print(Qx, rx)
    lambda_max = scipy.linalg.eigvalsh(Qx, eigvals = tuple([K-1, K-1]))
    gradientz = 2 * np.matmul(Qx  - np.asscalar(lambda_max) * np.eye(K), z) + rx

    return gradientz

def NominalModel(z, length_z, Pjk, bjk, cjk, n, m, l, K, model = None):
    # we assume that the first constraint (index 0) is the objective function, the input
    # z should be a list of z's corresponding to different constraints

    if not (model == None):
        i_history = length_z - 1

        varnames = [('xvar[' + str(i) + ']') for i in range(n)]
        #print(varnames)
        xvar = pd.Series([model.getVarByName(v) for v in varnames])
        objvar = model.getVarByName('objvar')
        #print(model.getVars())

        # objective
        [Qz, rz, sz] = Transform_to_Qz_rz_sz(z[i_history][0], Pjk, bjk, cjk, n, m, l, K, 0)
        model.addConstr(quicksum(float(Qz[i, j]) * xvar[i] * xvar[j] for i  in range(n) for j  in range(n))\
                + quicksum(float(rz[i]) * xvar[i] for i  in range(n))\
                    + float(sz) <= objvar)
        #
        for i_constraint in range(1, m):
                    [Qz, rz, sz] = Transform_to_Qz_rz_sz(z[i_history][i_constraint], Pjk, bjk, cjk, n, m, l, K, i_constraint)
                    model.addConstr(quicksum(float(Qz[i, j]) * xvar[i] * xvar[j] for i  in range(n) for j  in range(n))\
                        + quicksum(float(rz[i]) * xvar[i] for i  in range(n))\
                        + float(sz) <= 0.0)
    else:
        model = Model('Nominal problem')
        model.setParam( 'OutputFlag', False)
        model.setParam( 'Threads', 1)
        xvar = model.addVars(n, lb=-GRB.INFINITY, name='xvar')
        objvar = model.addVar(lb=-GRB.INFINITY, name='objvar')
        model.setObjective(objvar, GRB.MINIMIZE)
        model.addConstr(quicksum(xvar[i] * xvar[i] for i in range(n)) <= 1.0)

        for i_history in range(length_z):
            [Qz, rz, sz] = Transform_to_Qz_rz_sz(z[i_history][0], Pjk, bjk, cjk, n, m, l, K, 0)
            # possibly all of that can be accelerated by referring directly to gurobi API from the matrix level
            model.addConstr(quicksum(float(Qz[i, j]) * xvar[i] * xvar[j] for i  in range(n) for j  in range(n))\
                + quicksum(float(rz[i]) * xvar[i] for i  in range(n))\
                    + float(sz) <= objvar)

            if(m > 1):
                for i_constraint in range(1, m):
                    [Qz, rz, sz] = Transform_to_Qz_rz_sz(z[i_history][i_constraint], Pjk, bjk, cjk, n, m, l, K, i_constraint)
                    model.addConstr(quicksum(float(Qz[i, j]) * xvar[i] * xvar[j] for i  in range(n) for j  in range(n))\
                        + quicksum(float(rz[i]) * xvar[i] for i  in range(n))\
                        + float(sz) <= 0.0)

    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        infeasible = False
        xvar_ = model.getAttr('X', xvar)
        xvar_output = np.matrix(np.zeros((n, 1)))
        for j in range(n):
            xvar_output[j] = xvar_[j]

        return_val = model.objVal
        xvar_ = xvar_output
    else:
        infeasible = True
        xvar_ = 0
        return_val = 0

    return xvar_, return_val, infeasible, model

### Make it callable with index - input is the entire list + index
### for projection on U make lambda_max
### Make a vector of projection operators because for the objective it will be without lambda_max

def RQ_ProjectionZ_All(z_new):
    #project onto unit balls with a given radius
    for i in range(0, len(z_new)):
        z_new[i] = RQ_ProjectionZ_Single(z_new[i])

    return z_new[:]

def RQ_ProjectionZ_Single(z_input):
    #project onto unit balls with a given radius
    if np.linalg.norm(z_input) > 1:
        z_output = z_input / np.linalg.norm(z_input)
    else:
        z_output = z_input

    return z_output

def RQ_ProjectionU_Single(u_input, lambda_bar):
    #project onto unit balls with a given radius
    z_input = u_input[:-1]
    lambda_input = np.asscalar(u_input[-1])
    lambda_output = lambda_input
    z_output = z_input
    norm_z=np.linalg.norm(z_input)
    if norm_z>lambda_input:
        lambda_output = np.maximum(np.minimum((norm_z + lambda_input)/2, lambda_bar), 0)
    elif lambda_input>lambda_bar:
            lambda_output=lambda_bar

    if norm_z>lambda_output:
        z_output = lambda_output * z_input / norm_z


    return z_output, lambda_output

def RQ_ProjectionX(x):
    #project onto unit balls with a given radius
    if np.linalg.norm(x) > 1:
        x_new = x / np.linalg.norm(x)
    else:
        x_new = x

    return x_new

def RQ_WorstCaseZ(x, Pjk, bjk, cjk, n, m, l, K):
    #solves the generalized eigenvalue problems
    output_z = [] # for the list storing the worst case z's
    slack_z = [] # for the w-c values of the constraints

    for i_constraint in range(m):
        [Qx, rx, sx] = Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, i_constraint)
        [lambda_max, _] = scipy.linalg.eigh(Qx, eigvals = tuple([K-1, K-1]))
        #[lambda_max, _] = scipy.sparse.linalg.eigsh(Qx, k = 1)
        Qx -= np.eye(K) * (lambda_max)# + 10 ** -16)
        sx += lambda_max

        # scipy implementation
        z = ca.MX.sym("z", K)
        f = - ca.dot(ca.mtimes(z.T, Qx).T, z)  - ca.dot(rx, z) - sx
        g = ca.dot(z, z) - 1
        lbX = ca.repmat(-np.inf, K)
        ubX = ca.repmat(np.inf, K)
        lbg = np.full(1, -np.inf)
        ubg = np.full(1, 0.0)
        X = ca.veccat(z)
        x0 = np.zeros((K, 1))

        nlp = {"f": f, "g": g, "x": X}
        solver = ca.nlpsol(
            "nlpsol",
            "ipopt",
            nlp,
            {
                "ipopt": {
                    "print_level": 0
                },
                "print_time": 0
            }
        )

        solution = solver(lbx=lbX, ubx=ubX, lbg=lbg, ubg=ubg, x0=x0)
        z_res = solution["x"].full()

        #output_z.append(output_zvar)
        #slack_z.append(model.objVal)

        output_z.append(z_res)
        slack_z.append(-np.asscalar(solution["f"].full()))

    return output_z, slack_z

def TRSP(Qx, rx, sx):
    K = rx.shape[0]
    z = ca.MX.sym("z", K)
    f = - ca.dot(ca.mtimes(z.T, Qx).T, z)  - ca.dot(rx, z) - sx
    g = ca.dot(z, z) - 1
    lbX = ca.repmat(-np.inf, K)
    ubX = ca.repmat(np.inf, K)
    lbg = np.full(1, -np.inf)
    ubg = np.full(1, 0.0)
    X = ca.veccat(z)
    x0 = np.zeros((K, 1))

    nlp = {"f": f, "g": g, "x": X}
    solver = ca.nlpsol(
        "nlpsol",
        "ipopt",
        nlp,
        {
            "ipopt": {
                "print_level": 0
            },
            "print_time": 0
        }
    )

    solution = solver(lbx=lbX, ubx=ubX, lbg=lbg, ubg=ubg, x0=x0)
    z_res = solution["x"].full()

    return z_res, -np.asscalar(solution["f"].full())

def RQ_FunctionValue(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint):

    [Qx, rx, sx] = Transform_to_Qx_rx_sx(x, Pjk, bjk, cjk, n, m, l, K, i_constraint)
    [lambda_max, _] = scipy.linalg.eigh(Qx, eigvals = tuple([K-1, K-1]))
    Qx -= np.identity(K) * lambda_max # + 10 ** -16)
    sx += lambda_max

    return (z.T * Qx * z + rx.T * z + sx)

###### This is the preparation of the instances
#m - number of constraints
#n - dimension of variable x
#l - number of rows in linear transformation of x
#K - the dimension of the z_{j}
#i_instance - the generated instance
#max_time - maximum time for running methods

m = int(sys.argv[1])
n = int(sys.argv[2])
l = int(sys.argv[3])
K = int(sys.argv[4])
i_instance = int(sys.argv[5])
max_time = int(sys.argv[6])

# Prepare the data to run the function
#the matrices construction ensure that 0 is in the interior of the feasible set
np.random.seed(i_instance)
matrices_scaling = 1

# preallocations for lists of matrices
Pjk = [] # list of matrices P^k_j j = 1...m
for j in range(m):
    Pjk.append([])
    for k in range(K + 1):
        Pjk[j].append(matrices_scaling * np.matrix(2 * np.random.rand(l, n) - 1))

# function definitions as in Ho-Nguyen, only that the sign in front of b_i and c_i is +, not -
# comment: last entry is the nominal value, the first K are multiplied by subsequent perturbations

bjk = []
for j in range(m):
    bjk.append([])
    for k in range(K):
        bjk[j].append(np.matrix(np.zeros((n, 1))))

    bjk[j].append(matrices_scaling * np.matrix(- 2 * np.random.rand(n, 1) + 1))

# comment: last entry is the nominal value, the first K are multiplied by subsequent perturbations
cjk = []
for j in range(m):
    cjk.append([])
    for k in range(K):
        cjk[j].append(0.0)

    cjk[j].append(-0.05)


# normalizing the numbers in the constraints
for j in range(m):
    max_bjk_norm = scipy.linalg.norm(bjk[j][-1], 2)
    total_Pjk = np.zeros((1, Pjk[j][-1].shape[1]))

    for k in range(K):
        total_Pjk = np.vstack((total_Pjk, Pjk[j][k]))

        if(scipy.linalg.norm(bjk[j][k], 2) > max_bjk_norm):
                max_bjk_norm = scipy.linalg.norm(bjk[j][k], 2)

    normalizer_0 = scipy.linalg.norm(Pjk[j][-1], 2) # normalize the P0 matrix by its norm
    normalizer_1 = scipy.linalg.norm(total_Pjk, 2) # normalize all the pertrurbation matrices by the total norm of the matrix obtained by stacking them on top of each other
    normalizer_2 = max_bjk_norm #normalize the b0 vector by its norm

    # normalizing them
    for k in range(K):
        Pjk[j][k] = Pjk[j][k] / normalizer_1

    Pjk[j][-1] = Pjk[j][-1] / normalizer_0
    bjk[j][-1] = bjk[j][-1] / normalizer_2

# comment: last entry is the nominal value, the first K are multiplied by subsequent perturbations

# Prepare for OCO
def inline_grds(x, z):
    return ComputeGradients(x, z, Pjk, bjk, cjk, n, m, l, K)

def inline_Nominal(z, length_z, model):
    return NominalModel(z, length_z, Pjk, bjk, cjk, n, m, l, K, model)

def inline_WCz(x):
    return RQ_WorstCaseZ(x, Pjk, bjk, cjk, n, m, l, K)

def inline_Pz(z):
    return RQ_ProjectionZ_All(z)

def inline_Pu(z,bar_lambda):
    return RQ_ProjectionU_Single(z,bar_lambda)

def inline_Px(x):
    return RQ_ProjectionX(x)

# computing the parameters needed in the problem
#see equation (26) in [Ho-Nguyen and Kılın¸c-Karzan, 2017]
sigma = 0
rho = 0
chi = 0
beta = 0

for j in range(m):
    candidate_sigma = 0

    for k in range(K):
        candidate_sigma = candidate_sigma + np.linalg.norm(Pjk[j][k], 'fro') ** 2

        P_max_norm = np.linalg.norm(Pjk[j][k], 2)
        if(P_max_norm > chi):
            chi = P_max_norm

    if(candidate_sigma > sigma ** 2):
        sigma = np.sqrt(candidate_sigma)

    A_norm = scipy.linalg.norm(Pjk[j][-1], 2)

    if(A_norm > rho):
        rho = A_norm

    b_norm = np.linalg.norm(bjk[j][-1], 2)
    if(b_norm > beta):
        beta = b_norm

constraint_magnitudes = np.zeros((m, 1))

for j in range(m):
    max_bjk_norm = scipy.linalg.norm(bjk[j][-1], 2)
    total_Pjk = Pjk[j][-1]

    for k in range(K):
        total_Pjk = np.vstack((total_Pjk, Pjk[j][k]))

        if(scipy.linalg.norm(bjk[j][k], 2) > max_bjk_norm):
                max_bjk_norm = scipy.linalg.norm(bjk[j][k], 2)

    constraint_magnitudes[j] = (2 * scipy.linalg.norm(total_Pjk, 2) ** 2) \
        + max_bjk_norm + np.abs(cjk[j][-1])

# computing Gz and Gx
Omegaz = 0.5
Omegax =  np.power(2, 0.5)
Gz = 2 * (sigma ** 2 + rho * sigma)
Gx = 4 * (rho + np.sqrt(K) * sigma) ** 2 + beta

#the functions themselves
def f_function_builder(i_constraint):
    def function(x, z):
        return RQ_FunctionValue(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint)

    return function

f = []
for i_constraint in range(m):
    f.append(f_function_builder(i_constraint))

#get x gradient for ith constraint
def gx_function_builder(i_constraint):
    def function(x, z):
        return ComputeGradientX(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint)

    return function

gx = []
for i_constraint in range(m):
    gx.append(gx_function_builder(i_constraint))

#get z_i gradient for the ith constraint
def gz_function_builder(i_constraint):
    def function(x, z):
        return ComputeGradientZ(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint)

    return function

gz = []
for i_constraint in range(m):
    gz.append(gz_function_builder(i_constraint))
    #myfun = lambda x, z: ComputeGradientZ(x, z, Pjk, bjk, cjk, n, m, l, K, i_constraint)
    #gz.append(myfun)

#accuracy
epsilon=0.001

#put everything in the data structure
Data = {"gx": gx, #subgradients in x
"gz": gz, #the subgradient in z
"Gx": Gx, #upper bound on the norm of the subgradients in x
"Gz": Gz * np.ones((m, 1)), #upper bound on the norm of the subgradients in z
"Omegax": Omegax, #bound on the diameter of the set X
"Omegaz": Omegaz * np.ones((m, 1)), #bound on the diameter of the set Z
"ub_z": np.power(2 * Gz ** 2 * Omegaz, 0.5), # and upper bound on the convergence coefficient of z
"nb_of_constraints": m, #number of constraints
"nb_of_variables": n, #dimension of x
"nb_of_z": K * np.ones((m, 1)), #dimension of z for each constraint
"Nominal": inline_Nominal, #nominal problem
"WCz": inline_WCz, #problem of computing the worst case z for constraint i
"f": f, #constraints functions
"Pz": inline_Pz, #projection of z
"Px": inline_Px} #projection of x

logical_range = [False, True]

z0 = []
for i in range(0, m):
    z0.append(np.zeros((K, 1)))

# starting point correction
x0,_,_,model = inline_Nominal([z0], 1, None)
#x0=np.random.rand(n,1)
#x0=x0/np.linalg.norm(x0)

file_tag = "m=" + repr(m) + "_n=" + repr(n) + "_l=" + repr(l) + "_K=" + repr(K) + "_i=" + repr(i_instance)+"_"
output_folder = "."


for i_FOP in logical_range: #[True]
    for i_FOD in logical_range:
        if i_FOP==False:
            print_every_iter=1
        else:
            print_every_iter=100

        if((i_FOP)):# or (not i_FOD)):
        #if((not i_FOP) and (not i_FOD)):
		#if (i_FOP):
        #if (False):
            [x_res, x_curr, max_iter, total_iter, tot_time] = \
                    PrimalDualOCO(Data,x0, model, z0, epsilon, i_FOP, i_FOD, constraint_magnitudes[0],  -constraint_magnitudes[0], max_time, print_every_iter,output_folder,file_tag)

#Updating Data structure for SGSP
Data = {
"grds": inline_grds,
"gx": gx,
"gz": gz,
"Gx": Gx,
"Gz": Gz * np.ones((m, 1)),
"Omegax": Omegax,
"Omegaz": Omegaz * np.ones((m, 1)),
"ub_z": np.power(2 * Gz * Omegaz, 0.5), 
"nb_of_constraints": m,
"nb_of_variables": n,
"nb_of_z": K * np.ones((m, 1)),
"Nominal": inline_Nominal,
"WCz": inline_WCz,
"f": f,
"f_max": constraint_magnitudes,
"Pu": np.repeat(inline_Pu, m),
"Px": inline_Px,
"lb":-1}

print_every_iter=100

[x_Slater,x_res, x_curr, max_iter, Slater_time,total_iter, tot_time] = SGSP(Data, x0, epsilon, max_time, print_every_iter,output_folder, file_tag)
