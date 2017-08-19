import numpy as np
from scipy import special
import time
#==============================================================================
# get matrix [j,d]=n^{(d)}_j from pi.01mat.list (
# checked; but without colnames
#==============================================================================
#J = 2
#D = 3
#pi_01mat_list = [np.matrix('1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0'), 
#                 np.matrix('0 1 0; 0 1 0; 1 0 0; 0 0 1')]

def getNj_JD(J, D, pi_01mat_list):
    nj_JD = np.zeros(shape=(J, D))
    for d in range(D):
        nj_JD[:,d] = np.sum(pi_01mat_list[d], axis=0)
    return (nj_JD)
    
#getNj_JD(J, D, pi_01mat_list)


#==============================================================================
# extract 0-1 matrix (n_d by J) for pi from vector pi.1toJvec for replicate d
# checked
#==============================================================================
#nd = 5
#J = 3

#pi_1toJvec = np.array([0, 1, 0, 2, 1])
def getPi_01mat_fromPi_1toJvec(nd, J, pi_1toJvec):
    pi_01mat = np.zeros(shape=(nd, J))
    pi_01mat[range(nd), pi_1toJvec] = 1
    return (pi_01mat)

#getPi_01mat_fromPi_1toJvec(nd, J, pi_1toJvec)

#==============================================================================
# extract vector pi_1toJvec from 0-1 matrix (n_d by J) for pi
# checked
#==============================================================================
#J = 3
#pi_01mat = np.matrix('1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0')

#def getPi_1toJvec_fromPi_01mat(J, pi_01mat):
#    pi_1toJvec = np.squeeze(np.array(np.where(pi_01mat)[1]+1))
#    return (pi_1toJvec)

#getPi_1toJvec_fromPi_01mat(J, pi_01mat)






#==============================================================================
# update Pi_01mat_list, Pi_1toJvec_list, Nj_JD
# may need double check (because random multinomial) 
#==============================================================================
#J = 3
#D = 2
#
#currPi_01mat_list = [np.matrix('1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0'), 
#                     np.matrix('0 1 0; 0 1 0; 1 0 0; 0 0 1')]
#n_D = np.array([5, 4])
#E_D = np.array([1, 1])
#lambda_J = np.array([0.3, 0.5, 0.2])
#logPrReads_jMinus1_list = [np.matrix('0 0 0; 0 -2 -4; 0 3 -1; 0 1 -1; 0 0 0', dtype='f'),
#                           np.matrix('0 0 0; 0 -1 -4; 0 0 0; 0 1 -2', dtype='f')]
                           
def updatePi_list(J, D, n_D, currPi_01mat_list, E_D, lambda_J, logPrReads_jMinus1_list):
    currNj_JD = getNj_JD(J, D, currPi_01mat_list)
    pi_Ed1_temp1 = lambda_J + np.sum(currNj_JD[:,E_D==1],axis=1)
    
    newPi_01mat_list = [None]*D
#    newPi_1toJvec_list = [None]*D
    newNj_JD = np.zeros(shape=(J, D))
    
    for d in range(D):
        logPrReads_jMinus1_mat = logPrReads_jMinus1_list[d]
        pi_d_temp = currPi_01mat_list[d]
        if E_D[d] == 0:
            pi_Ed0_temp1 = (lambda_J+currNj_JD[:,d]) - pi_d_temp
            pi_Ed0_temp2 = pi_Ed0_temp1[:,0].reshape(-1,1)
            logPiProb_ndByJ = ( special.gammaln(pi_Ed0_temp1 + 1) - special.gammaln(pi_Ed0_temp1) - 
                            special.gammaln(pi_Ed0_temp2+1) + special.gammaln(pi_Ed0_temp2) + 
                            logPrReads_jMinus1_mat )
        else:
            pi_Ed1_temp2 = pi_Ed1_temp1 - pi_d_temp
            pi_Ed1_temp3 = pi_Ed1_temp2[:,0].reshape(-1,1)
            logPiProb_ndByJ = ( special.gammaln(pi_Ed1_temp2 + 1) - special.gammaln(pi_Ed1_temp2) - 
                            special.gammaln(pi_Ed1_temp3 + 1) + special.gammaln(pi_Ed1_temp3) + 
                            logPrReads_jMinus1_mat )
        
        piProb_ndByJ = np.exp( logPiProb_ndByJ )
        piProb_ndByJ = piProb_ndByJ/piProb_ndByJ.sum(axis=1)
#        newPi_01mat_list[d] = np.zeros(shape=(n_D[d], J))
#        for x in range(n_D[d]):
#            newPi_01mat_list[d][x, :] = np.random.multinomial(1, np.squeeze(np.array(piProb_ndByJ[x, :])), size=1)
#            newPi_1toJvec_list[d] = getPi_1toJvec_fromPi_01mat( J=J, pi_01mat=newPi_01mat_list[d] )
        prob_1tond =np.random.uniform(0, 1, size=n_D[d]).reshape(n_D[d],1)
        prob_cdf = np.cumsum(piProb_ndByJ,axis=1)
        bin = np.argmin(prob_1tond > prob_cdf, axis=1)
        pi_1toJvec = np.squeeze(np.asarray(bin))
        newPi_01mat_list[d] = getPi_01mat_fromPi_1toJvec(n_D[d], J, pi_1toJvec)
        newNj_JD[:, d] = newPi_01mat_list[d].sum(axis=0)
#    return([newPi_01mat_list, newPi_1toJvec_list, newNj_JD])
    return([newPi_01mat_list, newNj_JD])
        
#updatePi_list(J, D, n_D, currPi_01mat_list, E_D, lambda_J, logPrReads_jMinus1_list)



#==============================================================================
# update E_d (given d)
# checked
#==============================================================================
#d = 0
#D = 2
#J = 3
#n_D = np.array([5, 4])
#currE_D = np.array([1, 1])
#lambda_J = np.array([0.3, 0.5, 0.2])
#currPi_01mat_list = [np.matrix('1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0'), 
#                     np.matrix('0 1 0; 0 1 0; 1 0 0; 0 0 1')]
#currNj_JD = getNj_JD(J, D, currPi_01mat_list)
#gamma = 0.8

def updateEd( d, J, D, n_D, currNj_JD, currE_D, gamma, lambda_J):
    lambdaSum= np.sum(lambda_J)
    gammaOdds = gamma/(1-gamma)
    nj_d_J = currNj_JD[:,d]
    
    E_minus_d = np.delete(currE_D, d)
    nj_JD_minus_d = np.delete(currNj_JD, d, axis=1)
    temp1 = lambda_J + nj_JD_minus_d[:, E_minus_d==1].sum(axis=1)
    temp2 = np.sum(temp1)
    
    num1_index = np.sum(special.gammaln(temp1+nj_d_J))-special.gammaln(temp2 + n_D[d])
    denom1_index = np.sum(special.gammaln(temp1)) - special.gammaln(temp2)
    
    num2_index = np.sum(special.gammaln(lambda_J)) - special.gammaln(lambdaSum)
    denom2_index = np.sum(special.gammaln(lambda_J+nj_d_J)) - special.gammaln(lambdaSum+n_D[d])
    
    exponent = min(100, num1_index + num2_index - denom1_index - denom2_index)
    odds = np.exp(exponent) * gammaOdds
    Ed = np.random.binomial(1, p = odds/(1+odds), size=1)
    return(Ed)

#updateEd( d, J, D, n_D, currNj_JD, currE_D, gamma, lambda_J)



#==============================================================================
# update \bm E_d 
# checked
#==============================================================================
#D = 2
#J = 3
#n_D = np.array([5, 4])
#currE_D = np.array([1, 1])
#lambda_J = np.array([0.3, 0.5, 0.2])
#currPi_01mat_list = [np.matrix('1 0 0; 0 1 0; 1 0 0; 0 0 1; 0 1 0'), 
#                     np.matrix('0 1 0; 0 1 0; 1 0 0; 0 0 1')]
#currNj_JD = getNj_JD(J, D, currPi_01mat_list)
#gamma = 0.8

def updateE_D( J, D, n_D, currNj_JD, currE_D, gamma, lambda_J):
    gammaOdds = gamma/(1-gamma)
    lambdaSum = sum(lambda_J) 
    
    for d in range(D):
        currEd = updateEd( d=d, J=J, D=D, n_D=n_D, 
                          currNj_JD=currNj_JD, 
                          currE_D=currE_D, 
                          gamma=gamma, 
                          lambda_J=lambda_J)
        currE_D[d] = currEd
    return(currE_D)

#updateE_D( J, D, n_D, currNj_JD, currE_D, gamma, lambda_J)




#==============================================================================
# update parameter \gamma
# checked
#==============================================================================
#D = 4
#D1 = 3
#a = 8
#b = 2

def updateGamma(D, D1, a, b):
    newGamma = np.random.beta(a= D1+a, b=D-D1+b, size=1)
    return(newGamma)
    
    
    
#==============================================================================
# main function
#==============================================================================
def Gibbs_Model(maxIter, burnIn, 
               J, D, n_D,
               A, B, lambda_J,
               initE_D, 
               logPrReads_jMinus1_list, 
               seed):
    if seed > 0:
        np.random.seed(seed=seed)
        
    Es_Iter = np.zeros(shape=(maxIter, D))
    Pi_lists_Iter = [None]*maxIter
    Gammas_Iter = np.array([None]*maxIter) 
    
    currE_D = initE_D
    currGamma = np.random.beta(a=A, b=B, size=1)
    
    initPi_01mat_list = [None]*D
    for d in range(D):
        initPi_01mat_list[d] = np.zeros(shape=(n_D[d], J))
        for x in range(n_D[d]):
            initPi_01mat_list[d][x, :] = np.random.multinomial(1, [1/J]*J, size=1)
    currPi_01mat_list = initPi_01mat_list
    
#    currPi_1toJvec_list = [None]*D
#    for d in range(D):
#        currPi_1toJvec_list[d] = getPi_1toJvec_fromPi_01mat(J, currPi_01mat_list[d])
    
    currNj_JD = getNj_JD(J=J, D=D, pi_01mat_list=currPi_01mat_list)
#    time1=[0]*maxIter
#    time2=[0]*maxIter
    for i in range(maxIter):
#        if i%1000 == 0:
#            print("Iter", i)
#        start_time1 = time.clock()
        newE_D = updateE_D( J=J, D=D, n_D=n_D, currNj_JD=currNj_JD, currE_D=currE_D, 
                         gamma=currGamma, lambda_J=lambda_J)
#        time1[i] = time.clock()-start_time1
        Es_Iter[i, :] = newE_D
        currE_D = newE_D
        currD1 = sum(currE_D)
        
#        start_time2 = time.clock()
        currPi_list = updatePi_list( J, D, n_D,  currPi_01mat_list=currPi_01mat_list, 
                                  E_D=currE_D, lambda_J=lambda_J, 
                                  logPrReads_jMinus1_list = logPrReads_jMinus1_list )
#        time2[i] = time.clock()-start_time2
        Pi_lists_Iter[i] = currPi_list
        currPi_01mat_list =  currPi_list[0]
#        currPi_1toJvec_list = currPi_list[1]
#       currNj_JD =  currPi_list[2]
        currNj_JD =  currPi_list[1]
    
        currGamma = updateGamma( D, D1=currD1, a=A, b=B )
        Gammas_Iter[i] = currGamma
        
    Pi_lists_Iter = Pi_lists_Iter[-(maxIter - burnIn):]
    result = { "Es_Iter": Es_Iter[-(maxIter - burnIn):],
            "Pi_lists_Iter": Pi_lists_Iter,
            "Gammas_Iter": Gammas_Iter[-(maxIter - burnIn):] }
    
    Gammas_Iter = result["Gammas_Iter"]
    estGamma = np.mean(Gammas_Iter)
    Es_Iter = result["Es_Iter"]
    estE_D = np.mean(Es_Iter,axis=0)
    estE_D[estE_D >= 0.5] = 1
    estE_D[estE_D < 0.5] = 0

    Pi_lists_Iter = result["Pi_lists_Iter"]
    ind = np.where(np.sum(Es_Iter,axis=1)!=0)[0]
    iter_ef = len(ind)
    Es_1s = Es_Iter[ind,:]
#    Pi_count_Iter = np.zeros(shape=(iter_ef, J)) 
    Pi_count = np.array([0]*J)
    for i in range(iter_ef):
        Nj_JD = Pi_lists_Iter[ind[i]][1]
 #       Pi_1toJvec_list = Pi_lists_Iter[ind[i]][1]
        indexEd = np.where(Es_1s[i, :]==1)[0]
    
#        a = np.array([], dtype="int64")
#        for j in indexEd:
#            a = np.hstack((a,Pi_1toJvec_list[j]))
        
        a = Nj_JD[:,indexEd]
        Pi_count = Pi_count + np.sum(a, axis=1)
#        tauTemp = np.zeros(shape=(1, J))
#        count = np.bincount(a)
#        tauTemp[0, 0 : (np.size(count)-1)] = count[1:]
#        Pi_count_Iter[i, :] = tauTemp

#    tau_Iter_E1 = Pi_count_Iter/Pi_count_Iter.sum(axis=1, keepdims=True)

#    tau_JSLIDE = tau_Iter_E1.mean(axis=0)
    tau_JSLIDE = Pi_count/(sum(Pi_count))
    return{"tau_JSLIDE":tau_JSLIDE,  "estE_D":estE_D, "estGamma":estGamma}

