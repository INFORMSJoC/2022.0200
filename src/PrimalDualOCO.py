import numpy as np
import scipy
import scipy.linalg
import scipy.stats
import time

def PrimalDualOCO(Data, x0, model, z0, epsilon, First_Order_Primal, First_Order_Dual, tau_upper, tau_lower, max_time, print_every_iter=100, output_folder=None, file_tag=None):
    #Implementation of OCO code from ``Online First-Order Framework for Robust Convex Optimization'' by Nam Ho-Nguyen and Fatma Kılınc-Karzan
    #Data - parameters structure
    #x0 - initial point
    #z0 - initial point for dual problem
    #epsilon - feasibility and optimality tolerance
    #First_Order_Primal - bolean indicating if first order method is used on primal problem
    #First_Order_Dual - bolean indicating if first order method is used on dual variables
    #tau_upper - an upper bound for the value of the robust problem
    #tau_lower - a lower bound for the value of the robust problem
    #max_time - maximum running time allowed
    #print_every_iter - printing to file one in every print_every_iter iterations
    #output_folder - output folder for results
    #file_tag - file name prefix given by user

    use_gx = True
    gz_exists = True

    if not file_tag == None:
        if(First_Order_Primal):
            str1 = "T"
        else:
            str1 = "F"

        if(First_Order_Dual):
            str2 = "T"
        else:
            str2 = "F"

        print('Start OCO '+str1+str2)

        filestring1 = output_folder + "/" + file_tag + "PDOCO"+str1 + str2+"_Iteration.txt"
        file_object1 = open(filestring1, "w")
        file_object1.close()

        filestring2 = output_folder + "/" + file_tag + "PDOCO"+ str1 + str2 +"_Final.txt"
        file_object2 = open(filestring2, "w")

    ###################################################
    if 'gx' in Data: # this needs to be remade
        gx = Data['gx']
        gz = Data['gz']
        use_gx = True
        gz_exists = True

    #copy data
    f = Data['f'] #function values
    Px = Data['Px'] #projection of primal
    Pz = Data['Pz'] #projection of dual
    WCz = Data['WCz'] #compute worst case dual
    Nominal = Data['Nominal'] #Nominal problem for input z
    Omegaz = Data['Omegaz'] #upper bound on diameter of dual sets
    Gz = Data['Gz'] #upper bounds on subgradients in z
    ub_z = Data['ub_z'] #upper bound
    Omegax = Data['Omegax'] #upper bound on diameter of primal set
    Gx = Data['Gx'] #upper bounds on subgradients in x
    nb_of_variables = Data['nb_of_variables'] #number of primal variables
    nb_of_z = Data['nb_of_z'].astype(int) #list of number of variables z
    nb_of_constraints = Data['nb_of_constraints'] #number of constraints

    infeasible = True
    tau = (tau_upper + tau_lower) / 2

    #initialize
    x = x0
    x_new = x
    x_avg = x
    sumstepsize=1
    print(np.linalg.norm(x_avg))

    z=z0
    #z = Pz(z[:])
    z_new = z[:]

    #set max number of iterations
    const = np.zeros((nb_of_constraints, 1))

    #calculate number of iterations needed
    if First_Order_Primal and First_Order_Dual:
        max_iter = np.uint64(np.ceil(((np.power(2 * Omegax * (Gx ** 2), 0.5) + ub_z) / epsilon) ** 2))
    else:
        if First_Order_Dual:
            max_iter = np.uint64(np.ceil((ub_z / epsilon) ** 2))
        else:
            if First_Order_Primal:
                max_iter = np.uint64(np.ceil((np.power(2 * Omegax * Gx ** 2, 0.5) / epsilon) ** 2))
            else:
                max_iter = 10 ** 7

    #run bi-section on objective value
    outer_iter = 0
    x_res = x
    total_iter = 0
    tot_time = 0
    stepsize_x = np.power(2 * Omegax / (Gx ** 2), 0.5)
    stepsize_z = np.zeros((nb_of_constraints, 1))
    not_time_out = True

    #calculates step sizes
    for i in range(0, nb_of_constraints):
        stepsize_z[i] = np.power(2 * Omegaz[i] / (Gz[i] ** 2), 0.5)

    n_hist = 1
    z_hist = []
    z_hist.append(z[:])
    while (np.absolute(tau_upper - tau_lower) > epsilon or infeasible) and not_time_out:
        outer_iter = outer_iter + 1
        print_output_OCO(filestring1,x_avg,x_new,sumstepsize,outer_iter,-1,tot_time,tau_lower,tau_upper,tau,nb_of_constraints,WCz)

        #only relevant for first_order_primal=TRUE
        max_const=-np.Inf
        const = np.zeros((nb_of_constraints, 1))
        if First_Order_Primal:
            const[i] = f[i](x, z[i]) - tau_upper * (i == 0) # careful if i==0 then computes value
            #finds the maximum index
            if max_const < const[i]:
                max_const = const[i]
                max_i = i
        avg_const = const
        sum_f = np.amax(avg_const)


        if First_Order_Dual or First_Order_Primal:
            n_hist = 0
            z_hist = [] #cell(max_iter,1); #initialize z_hist if needed # careful here

        periodic_check = 1
        early_stop=False
        for iter in range(max_iter):
            #computes step for x
            tic = time.time()
            total_iter = total_iter + 1
            if (First_Order_Primal or First_Order_Dual) and not use_gx:
                [gradx, gradz] = Data['grds'](x, z[:]) # check if lists

            #we start by checking if we are doing first order primal and doing the step
            if First_Order_Primal:
                if use_gx: #computed gradient for x with respect to old constraints
                    gradx = gx[max_i](x, z[max_i])

                if np.linalg.norm(gradx, 2) > Gx + 10 ** -5:
                    print('Error: Oops')
                    1/0


                # updates x value (avg is updates later)
                if  np.linalg.norm(gradx)>0:
                    stepsize_x = np.power(2 * Omegax , 0.5)* np.power(1 / (iter + 1), 0.5)/ np.linalg.norm(gradx)
                else:
                    stepsize_x=0

                x_new = x - stepsize_x  * gradx
                x_new = Px(x_new)
                K_x = np.power(2 * Omegax * Gx ** 2 / (iter + 1), 0.5) / epsilon #The value of kappa_t for x

            #checks if first order dual
            if First_Order_Dual:
                max_norm_gradz = 0

                for i in range(0, nb_of_constraints):
                    if gz_exists: #not sure what this is
                        gradz_i = gz[i](x, z[i]) #updates gradient for z
                    else:
                        gradz_i = gradz[i]

                    if np.linalg.norm(gradz_i) > Gz[i]:
                        print('Error: Oops gradient z too large')
                        1/0

                    if  np.linalg.norm(gradz_i)>0:
                        stepsize_z[i] = np.power(2 * Omegaz[i] , 0.5)* np.power(1 / (iter + 1), 0.5)/ np.linalg.norm(gradz_i)
                    else:
                        stepsize_z[i]=0

                    max_norm_gradz = np.maximum(np.linalg.norm(gradz_i), max_norm_gradz)
                    z_new[i] = z[i] + np.asscalar(stepsize_z[i]) * gradz_i # updates according to stepsize

                #projecting z (avg updates later)
                z_new = Pz(z_new[:])
                K_z = ub_z / np.power(iter + 1, 0.5) / epsilon #upper bound on number of iterations in z

                if not First_Order_Primal: #solve nominal with respect to z
                    [x_new, tau_lower, infeasible, model] = Nominal([z[:]], 1, model)
                    if infeasible:#checks feasibility with regard to new point (why?)
                        tau_lower = np.Inf
                        print('Infeasible iteration')
                        break

                    K_x = 0

                if max_norm_gradz < (10 ** -20) and iter > 0:
                    print("OMG")
                    break

            else: #compute worst case dual
                if not First_Order_Primal: #solve nominal with respect to z historical realization
                    # this is the case where z is updated by pasimization oracle
                    # and x by an extended nominal solution using all the
                    # historical data

                    if (n_hist == 1):
                        x_new = x0
                        infeasible = False
                        tau_lower = model.objVal
                    else:
                        [x_new, tau_lower, infeasible, model] = Nominal(z_hist, n_hist, model)

                    if infeasible:
                        tau_lower = np.Inf
                        early_stop=True
                        toc = time.time()
                        tot_time = tot_time + (toc - tic)
                        print('Infeasible iteration')
                        break

                    #print("Lower bound update tau=%f,tau_lower=%f" % (tau,tau_lower))
                    K_x = 0

                [z_new, slacks_z] = WCz(x_new)

                if not First_Order_Primal: # When doing dual cuts and extended nominal
                    # Adds the last worst cases to history
                    n_hist = n_hist + 1
                    z_hist.append(z_new)

                if nb_of_constraints == 1 or np.amax(slacks_z[1:]) < epsilon:
                    infeasible = False
                    tau_upper = min(slacks_z[0],tau_upper)
                    #print("Upper bound update %f" % (tau_upper))
                    if slacks_z[0]<=tau:
                        early_stop=True
                        toc = time.time()
                        tot_time = tot_time + (toc - tic)
                        break

            if First_Order_Primal: #important update for next step
                max_const=-np.Inf
                for i in range(0, nb_of_constraints): #checks constraints for current point (not new one)
                    const[i] = f[i](x_new, z_new[i]) - tau * (i == 0) # careful if i==0 then computes value
                    #finds the maximum index
                    if max_const < const[i]:
                        max_const = const[i]
                        max_i = i
                avg_const = avg_const + const
                sum_f = np.amax(avg_const)
            toc = time.time()
            tot_time = tot_time + (toc - tic)

            if((iter + 1) == periodic_check):
                if First_Order_Dual:
                    if First_Order_Primal:
                        x_avg_tentative = (x_avg + x_new*stepsize_x) / (sumstepsize+stepsize_x)
                    else:
                        x_avg_tentative = (x_avg + x_new) / (iter+2)

                    [_, z_slack_tentative] = WCz(x_avg_tentative)
                    if(((nb_of_constraints > 1) and (np.amax(z_slack_tentative[1:]) <= 0)) or nb_of_constraints == 1) and (z_slack_tentative[0] <= tau):
                        tau_upper = z_slack_tentative[0]
                        early_stop=True
                        break
                    else:
                        periodic_check *= 2

            x_curr = x
            x = x_new
            z = z_new
            if First_Order_Primal:#update within for-loop
                x_avg = (x_avg + x_new*stepsize_x)
                sumstepsize=sumstepsize+stepsize_x #x_avg over the step sizes
            else:
                x_avg = (x_avg + x_new)
                sumstepsize= sumstepsize + 1
            if np.linalg.norm(x_avg/sumstepsize)>1+1e-5:
                print(iter,np.linalg.norm(x_avg/sumstepsize))
                1/0


            #if np.linalg.norm(z_new[0], 2) < 1 - 1e-4 and iter > 0:
                #print('Worst-case z is bad')
                #print(np.linalg.norm(z_new[0], 2))

            if (not file_tag == None) and ((iter + 1) % print_every_iter) == 0:
                print_output_OCO(filestring1,x_avg,x_new,sumstepsize,outer_iter,iter,tot_time,tau_lower,tau_upper,tau,nb_of_constraints,WCz)


            if (tau_upper - tau_lower) < epsilon: #and not First_Order_Primal:
                break

            if tot_time > max_time:
                not_time_out = False
                break

        # after the for loop is finished:
        if early_stop: #if early stop need to update x_avg
            if First_Order_Primal:
                x_avg=x_avg+stepsize_x*x_new
                sumstepsize=sumstepsize+stepsize_x
            else:
                x_avg=x_avg+x_new
                sumstepsize=sumstepsize+1
            if np.linalg.norm(x_avg/sumstepsize)>1:
                print(iter,np.linalg.norm(x_avg/sumstepsize))
                1/0


        # printing the final output of the iteration
        if early_stop or ((iter + 1) % print_every_iter) > 0:
            print_output_OCO(filestring1,x_avg,x_new,sumstepsize,outer_iter,iter,tot_time,tau_lower,tau_upper,tau,nb_of_constraints,WCz)

        x_avg = x_avg / sumstepsize

        #print("Before update - tau_lower= %f, tau_upper= %f, tau=%f" % (tau_lower,tau_upper,tau))
        if First_Order_Primal:
            sum_f = sum_f / sumstepsize

            #if First_Order_Dual:
            if not early_stop and sum_f > K_x * epsilon:
                infeasible = True
                tau_lower = tau
            else:
                x_res = x_avg
                infeasible = False
                tau_upper = min(tau,tau_upper)

        tau = (tau_upper + tau_lower) / 2
        #print("After update - tau_lower= %f, tau_upper= %f, tau=%f" % (tau_lower,tau_upper,tau))


    if not file_tag == None:
        [_, z_avg_slack_tentative] = WCz(x_res)
        [_, z_new_slack_tentative] = WCz(x_new)

        file_object2.write(' {:10d}\n'.format(max_iter))
        file_object2.write(' {:10d}\n'.format(total_iter))
        file_object2.write(' {:10.5f}\n'.format(tot_time))
        file_object2.write(' {:10.5f}\n'.format(z_avg_slack_tentative[0]))

        if(nb_of_constraints > 1):
            file_object2.write(' {:10.5f}\n'.format(np.amax(z_avg_slack_tentative[1:])))
        else:
            file_object2.write(' {:10.5f}\n'.format(0))

        file_object2.write(' {:10.5f}\n'.format(z_new_slack_tentative[0]))

        if(nb_of_constraints > 1):
            file_object2.write(' {:10.5f}\n'.format(np.amax(z_new_slack_tentative[1:])))
        else:
            file_object2.write(' {:10.5f}\n'.format(0))

        file_object2.write('\n')

        for i in range(x_res.shape[0]):
            file_object2.write(' {:10.5f}\n'.format(np.asscalar(x_res[i])))

        file_object2.write('\n')

        for i in range(x_res.shape[0]):
            file_object2.write(' {:10.5f}\n'.format(np.asscalar(x_new[i])))

        file_object2.close()

    outputfunction = [x_res, x_new, max_iter, total_iter, tot_time]

    return outputfunction


def print_output_OCO(filestring1,x_avg,x_new,sumstepsize,outer_iter,iter,tot_time,tau_lower,tau_upper,tau,nb_of_constraints,WCz):
    file_object1 = open(filestring1, "a")
    x_avg_tentative = x_avg / max(sumstepsize,1)
    [_, z_avg_slack_tentative] = WCz(x_avg_tentative)
    [_, z_new_slack_tentative] = WCz(x_new)

    # computing the slack
    if(nb_of_constraints > 1):
        violation_avg = np.amax(z_avg_slack_tentative[1:])
        violation_new = np.amax(z_new_slack_tentative[1:])
    else:
        violation_avg = 0
        violation_new = 0

    file_object1.write(' {:10d}'.format(outer_iter))
    file_object1.write(' {:10d}'.format(iter + 1))
    file_object1.write(' {:10.5f}'.format(tot_time))
    file_object1.write(' {:10.5f}'.format(z_avg_slack_tentative[0]))
    file_object1.write(' {:10.5f}'.format(violation_avg))
    file_object1.write(' {:10.5f}'.format(z_new_slack_tentative[0]))
    file_object1.write(' {:10.5f}'.format(violation_new))
    file_object1.write(' {:10.5f}'.format(np.asscalar(np.asmatrix(tau_lower))))
    file_object1.write(' {:10.5f}'.format(np.asscalar(np.asmatrix(tau_upper))))
    file_object1.write(' {:10.5f}\n'.format(np.asscalar(np.asmatrix(tau))))
    file_object1.close()
