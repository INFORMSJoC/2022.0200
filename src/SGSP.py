# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 13:21:30 2020

@author: Shim
"""
import numpy as np
import time

def SGSP(Data, x0, accuracy, max_time=None, print_every_iter=100, output_folder=None, file_tag=None):
    #input : 
    # Data - problem data including functions, variables dimension, gradients,
    #        projections, upper bounds on function values and norms of variables and gradients
    # x0 - Initial point for x
    # accuracy - Accuracy we are trying to obtain for feasibility and optimality
    # max_time - a limit on the time the program is allowed to run.
    # print_every_iter - print status for optimality and feasibility bounds
    # output_folder - the folder output files will be written to
    # file_tag - the prefix for the file name given by the user, if empty results not printed
    if not file_tag == "":
        #file name for saving SGSP iterations
        filestring1 = output_folder + "/" + file_tag + "SGSP_Iteration.txt"
        file_object1 = open(filestring1, "w")
        file_object1.close()
        
        #file name for final solution
        filestring2 = output_folder + "/" + file_tag + "SGSP_Final.txt"
        file_object2 = open(filestring2, "w")

    lb = Data['lb'] #lower bound on the function values for the objective
    f = Data['f'] #functions for each constraint
    gx = Data['gx'] #gradient function with respect to x
    gz = Data['gz'] #gradient function with respect to z
    grads = Data['grds'] #all gradients
    Px = Data['Px'] #projection with respect to x
    Pu = Data['Pu'] #projections with respect to u
    WCz = Data['WCz'] #calculates worst case z and slack for each constraint
    Omegaz = Data['Omegaz'] #upper bound on the norm of z
    Gz = Data['Gz'] #upper bound on norm of the gradients gz
    ub_z = Data['ub_z'] # upper bound on the norm of each z
    Omegax = Data['Omegax'] #upper bound on the norm of x under the appropriate norm
    f_max = Data['f_max'] #upper bound on the magnitude of all f on their domain
    Gx = Data['Gx'] #upper bound on norm of the gradients gx
    nb_of_variables = Data['nb_of_variables'] #dimension of x
    nb_of_z = Data['nb_of_z'].astype(int) #number of z
    nb_of_z = np.squeeze(np.asarray(nb_of_z),1)
    nb_of_constraints = Data['nb_of_constraints'] #number of constraints in problem
    Nominal = Data['Nominal'] # 

    n=x0.shape[0] #dimension of variable x

    x_res = x0
    x_curr = x0

    z_curr = []
    u_curr = []
    for i in range(nb_of_constraints):
        z_curr.append(np.zeros((np.asscalar(nb_of_z[i]), 1)))
        u_curr.append(np.zeros((np.asscalar(nb_of_z[i])+1, 1))) #lifted z variable

    mu_curr=np.zeros(nb_of_constraints)
    mu_curr[0]=1

    total_iter = 0
    stepsize_x = np.power(2 * Omegax / (Gx ** 2), 0.5) #theoretical coeffient for stepsize of x
    stepsize_z = np.zeros((nb_of_constraints, 1))
    for i in range(nb_of_constraints):
        stepsize_z[i] = np.power(2 * Omegaz[i] / (Gz[i] ** 2), 0.5) #theoretical coeffient for stepsize of u_i

    init_time=time.time()
    total_time = 0 #calculates time
    sumstepsize = 1 #for avareging x
    #Calculate Slater point if there are constraints and xSlater not provided
    #Since we are using an epigraph constraint for the objective this corresponds with nb_of_constraints>=2.
    if not 'xSlater' in Data and not nb_of_constraints==1: #if there is only one constraint it is the objective
        outer_stage=0 #calculating slater
        K_prev=0
        K=2 #how many iteration to try
        #calculates worst case z and slack for each constraint
        [z_curr, slack]=WCz(x_curr)
        x_avg=x_curr
        t_curr=np.max(slack[1:nb_of_constraints])
        t_avg=t_curr
        delta=t_curr #approximate Slater point for current problem
        t_bar=t_curr+delta #upper and lower bound for t
        t_ubar=-t_bar
        lambda_bar=(t_bar-t_ubar)/delta #upper bound on dual solution
        epsSlater=0
        if not file_tag == None:
            #print to file
            print_output_SGSP(filestring1,x_avg,x_curr,1,outer_stage,-1,total_time,nb_of_constraints,WCz)
        print('Calculate Slater Point')
        if t_bar<0:
            xSlater=x_curr
            epsSlater=-t_curr
            fSlater=slack[0]
            print('Found Slater Point')
        #Run SGSP exits when we find Slater point
        while total_time<max_time and t_bar>=0:
            start_time=time.time()
            #how many iteration to try for the current guess for epsilon-Slater
            for k in range(1+K_prev,K+K_prev+1):#
                ggx,ggz = grads(x_curr,z_curr)
                grad_L_x=np.zeros((nb_of_variables,1))
                grad_L_t=1-sum(mu_curr[1:nb_of_constraints])#gradient of Lagrangian with respect to t
                for i in range(1,nb_of_constraints): #this is only with respect to the actual constraints not objective
                    grad_L_x+=(mu_curr[i])*(np.asarray(ggx[i]))#gradient of Lagrangian with respect to x for constraint i
                x_curr_vector=(np.asarray(x_curr))
                if np.linalg.norm(grad_L_x)>0:
                    stepsize_x=np.power(2*Omegax,0.5)/(np.sqrt(k+1)*np.linalg.norm(grad_L_x)) #try normalized step size
                else:
                    stepsize_x=0
                x_curr_vector=Px(x_curr_vector-stepsize_x*grad_L_x) #projected subgradient for x
                x_curr=x_curr_vector.reshape(n,1)
                x_avg=(x_avg+stepsize_x*x_curr)

                t_curr=min(max(t_curr-stepsize_x*grad_L_t,t_ubar),t_bar) #projected subgradient for t
                
                sumstepsize+=stepsize_x
                glambda=0
                for i in range(1,nb_of_constraints):
                    u_curr[i][0:np.asscalar(nb_of_z[i])]=np.asarray(z_curr[i])*mu_curr[i]
                    u_curr[i][np.asscalar(nb_of_z[i])]=mu_curr[i]
                    if i==0: #this is only with respect to the actual constraints not objective
                        print("should not happen here")
                    else:
                        #computing gradient with respoct to lambda_i
                        glambda= (np.asarray(t_curr-f[i](x_curr,z_curr[i])+np.transpose(ggz[i]).dot(z_curr[i]))) 
                    gu = np.zeros((np.asscalar(nb_of_z[i])+1,1))
                    gu[0:np.asscalar(nb_of_z[i])]=-ggz[i]
                    gu[np.asscalar(nb_of_z[i])]=glambda
                    if np.linalg.norm(gu)==0:
                        stepsize_z[i]=0
                    else:
                        stepsize_z[i]=np.power(2*Omegaz[i],0.5)/(np.sqrt(k+1)*np.linalg.norm(gu)) #try normalized step size
                    u_step_matrix=(u_curr[i]-stepsize_z[i]*gu).reshape(nb_of_z[i]+1,1)

                    tilde_z_curr,mu_curr[i]=Pu[i](u_step_matrix,lambda_bar) #projetion
                    #converting z_tilde to z
                    if mu_curr[i]>0:
                        z_curr[i]=tilde_z_curr.reshape((nb_of_z[i],1))/mu_curr[i]
                    else:
                        z_curr[i]=z_curr[i]*0

            [z_temp, slack]=WCz(x_curr)
            [z_temp_avg, slack_avg]=WCz(x_avg/sumstepsize)
            t_temp=np.max(slack[1:nb_of_constraints])
            t_temp_avg=np.max(slack_avg[1:nb_of_constraints])
            finish_time=time.time()
            iter_time=finish_time-start_time
            total_time+=iter_time
            print_output_SGSP(filestring1,x_avg,x_curr,sumstepsize,outer_stage,k,total_time,nb_of_constraints,WCz)
            print(k,t_temp,t_temp_avg)
            if t_temp<0: #if current point is Slater point
                xSlater=x_curr
                epsSlater=-t_curr
                fSlater=slack[0]
                print('Found Slater Point')
                break
            elif t_temp_avg<0: #if average point is Slater point
                xSlater=x_avg/sumstepsize
                epsSlater=-t_temp_avg
                fSlater=slack_avg[0]
                print('Found Slater Point')
                break
            else: #if we have not obtained an epsilon-Slater point with epsSlater=t_bar
                if 2*t_temp<t_bar: #update epsSlater=t_temp
                    delta=t_temp
                    t_bar=min(t_temp+delta,t_bar)
                    t_ubar=-t_bar
                    lambda_bar=(t_bar-t_ubar)/delta
                    t_curr=t_temp
                K_prev=K #increase number of iteration by factor of 2
                K=K_prev*2
    else:
        if 'epsSlater' in Data:
            #slater point given
            epsSlater=Data['epsSlater']
            xSlater=Data['xSlater']
            fSlater=WCz(xSlater)[1][0]
        else:
            #slater point is not given but nb_of_constraints==1 so no real constraints
            epsSlater=0

    if epsSlater>0:
         ub_lambda=(fSlater-lb)/epsSlater
         x_avg=xSlater
         x_curr=xSlater
    else:
        ub_lambda=0
        fSlater=WCz(x0)[1][0]
        xSlater=x0
        x_avg=x0
        x_curr=x0


    #update bounds
    GGx=Gx*(1+ub_lambda)
    Gi=np.sqrt(Gz**2+(f_max+Omegaz*Gz)**2)
    #compute upper bound on the number of iterations needed
    max_iter = np.uint64(np.ceil((( Omegax* GGx+1*Gi[0]+sum(ub_lambda*Gi[1:nb_of_constraints])) / accuracy) ** 2))[0]
    now_time=time.time()
    Slater_time=now_time-init_time
    total_time=Slater_time #update time used to find Slater point
    print('Start SGSP')
    outer_stage=1 #starting actual iterations
    mu_curr[0]=1
    max_violation=[]
    function_value=[fSlater]
    sumstepsize=1
    iter=-1
    print_output_SGSP(filestring1,x_avg,x_curr,sumstepsize,outer_stage,-1,total_time,nb_of_constraints,WCz)
    for iter in range(max_iter):
        if total_time > max_time:
            iter=iter-1
            break
        total_iter += 1
        start_time=time.time()

        # CHUNK 1
        #compute gradient of x
        ggx,ggz=grads(x_curr,z_curr)
        grad_L_x=(ggx[0])#gx[0](x_curr,z_curr[0])))
        for i in range(1,nb_of_constraints):
            grad_L_x+=mu_curr[i]*(ggx[i])#gx[i](x_curr,z_curr[i])))
        x_curr_vector=x_curr #HERE

        if np.linalg.norm(grad_L_x)>0:
            stepsize_x=np.power(2*Omegax,0.5)/(np.linalg.norm(grad_L_x)*np.sqrt(iter+1)) #try normalized step size
        else:
            stepsize_x=0
        x_curr_vector=Px(x_curr_vector-stepsize_x*grad_L_x)
        x_curr=x_curr_vector.reshape(n,1)
        x_avg=(x_avg+stepsize_x*x_curr)
        sumstepsize+=stepsize_x

        # CHUNK 2
        glambda=0
        for i in range(nb_of_constraints):
            u_curr[i][0:(nb_of_z[i])] = np.asarray(z_curr[i])*mu_curr[i]
            u_curr[i][(nb_of_z[i])] = mu_curr[i]
            #ggz = ((gz[i](x_curr,z_curr[i])))
            if i>0:
                glambda= (np.asarray(-f[i](x_curr,z_curr[i])+np.transpose(ggz[i]).dot(z_curr[i])))
            gu = np.zeros(((nb_of_z[i])+1,1))
            gu[0:(nb_of_z[i])] = -ggz[i]
            gu[nb_of_z[i]] = glambda

            if np.linalg.norm(gu)==0:
                stepsize_z[i]=0
            else:
                stepsize_z[i]=np.power(2*Omegaz[i],0.5)/(np.linalg.norm(gu)*np.sqrt(iter+1)) #try normalized step size
            #print(u_curr[i][20],gu[20])
            #print(stepsize_z[i])
            u_step_matrix=(u_curr[i]-stepsize_z[i]*gu).reshape(nb_of_z[i]+1,1)

            if i>0:
                tilde_z_curr, mu_curr[i] = Pu[i](u_step_matrix, ub_lambda)
            else: # when i=0 mu should always be 1
                tilde_z_curr, mu_curr[i] = Pu[i](u_step_matrix,1)
                if mu_curr[i]<1: #this is another safe guard
                    print('oops should not be smaller than 1')
                    print(i, mu_curr[i])
                    1/0

            #update z for z_tilde
            if mu_curr[i]>0:
                z_curr[i]=tilde_z_curr.reshape((nb_of_z[i],1))/mu_curr[i]
            else:
                z_curr[i]=z_curr[i]*0

        finish_time=time.time()
        iter_time=finish_time-start_time
        total_time+=iter_time


        #printing results to file
        if not file_tag == None and (iter + 1) % print_every_iter == 0:
            print_output_SGSP(filestring1,x_avg,x_curr,sumstepsize,outer_stage,iter,total_time,nb_of_constraints,WCz)





    #prints last iteration
    if not file_tag == None and (iter + 1) % print_every_iter != 0:
        print_output_SGSP(filestring1,x_avg,x_curr,sumstepsize,outer_stage,iter,total_time,nb_of_constraints,WCz)


    if not file_tag == None:
        x_avg=x_avg/sumstepsize
        [_, z_avg_slack_tentative] = WCz(x_avg)
        [_, z_new_slack_tentative] = WCz(x_curr)

        file_object2.write(' {:10d}\n'.format(max_iter))
        file_object2.write(' {:10d}\n'.format(total_iter))
        file_object2.write(' {:10.5f}\n'.format(total_time))
        file_object2.write(' {:10.5f}\n'.format(Slater_time))
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
            file_object2.write(' {:10.5f}\n'.format(np.asscalar(x_avg[i])))

        file_object2.write('\n')

        for i in range(x_res.shape[0]):
            file_object2.write(' {:10.5f}\n'.format(np.asscalar(x_curr[i])))

        file_object2.write('\n')

        for i in range(x_res.shape[0]):
            file_object2.write(' {:10.5f}\n'.format(np.asscalar(xSlater[i])))

        file_object2.close()

    outputfunction = [xSlater,x_curr, x_avg/sumstepsize, max_iter, iter + 1, Slater_time,total_time]

    return outputfunction


#printing output to file
def print_output_SGSP(filestring1,x_avg,x_curr,sumstepsize,outer_stage,iter,total_time,nb_of_constraints,WCz):
    file_object1 = open(filestring1, "a")
    x_avg=x_avg/sumstepsize
    z_avg_slack_tentative = WCz(x_avg[:])[1]
    z_new_slack_tentative = WCz(x_curr[:])[1]
    if(nb_of_constraints > 1):
        violation_avg = np.amax(z_avg_slack_tentative[1:])
        violation_new = np.amax(z_new_slack_tentative[1:])
    else:
        violation_avg = 0
        violation_new = 0

    file_object1.write(' {:10d}'.format(outer_stage))
    file_object1.write(' {:10d}'.format(iter + 1))
    file_object1.write(' {:10.5f}'.format(total_time))
    file_object1.write(' {:10.5f}'.format(z_avg_slack_tentative[0]))
    file_object1.write(' {:10.5f}'.format(violation_avg))
    file_object1.write(' {:10.5f}'.format(z_new_slack_tentative[0]))
    file_object1.write(' {:10.5f}\n'.format(violation_new))
    #file_object1.write(' {:10.5f}\n'.format(np.asscalar(mu_curr)))
    file_object1.close()
