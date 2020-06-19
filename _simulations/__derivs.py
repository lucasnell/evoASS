

""" 

dVi / dVi 

reticulate::use_python('/usr/local/bin/python3')
reticulate::source_python('_simulations/__derivs.py')


For these, each row gives the dynamics of a particular component of the system, 
while columns give the responses to other components of the system 
(see http://www.scholarpedia.org/article/Equilibrium).


"""

import sympy
import theano
theano.config.cxx = ""
import theano.tensor as T
import numpy as np
import pandas as pd
from tqdm import tqdm
import math
pd.options.display.max_columns = 10



def no_negative(X):
    """Equivalent to ifelse(X < 0, 0, X)"""
    out = np.where(X < 0, 0, X)
    return out


def dVi_dVi(i, V, O, C, f, a0, s2):
    """Automatic differentiation of dVi/dVi using theano pkg"""
    i = int(i)
    Vi = T.dvector('Vi')
    Vhat = np.absolute(Vi + 2 * s2 * (
        ( a0 * O * T.exp(-1 * T.dot(Vi, Vi.T)) * Vi) - 
        ( f * T.dot(Vi, C) )
    ))
    J, updates = theano.scan(lambda i, Vhat, Vi : T.grad(Vhat[i], Vi), 
                         sequences=T.arange(Vhat.shape[0]), non_sequences=[Vhat, Vi])
    num_fun = theano.function([Vi], J, updates=updates)
    out_array = num_fun(V[i])
    return out_array


def dVi_dVk(i, k, N, V, D, f, a0, C, s2):
    """Automatic differentiation of dVi/dVk using theano pkg"""
    i = int(i)
    k = int(k)
    Vi = V[i]
    Ni = N[i]
    Nk = N[k]
    P = [np.exp(-1 * np.dot(np.dot(V[j], D), V[j].T)) * N[j] 
         for j in range(0, len(N)) if j != i and j != k]
    P = np.sum(P) + Ni
    Vk_ = T.dvector('Vk_')
    Vhat = np.absolute(Vi + 2 * s2 * ( T.dot(Nk * T.exp(-1 * T.dot(T.dot(Vk_, D), Vk_.T)) + P, 
                              a0 * T.dot(T.exp(-1 * T.dot(Vi, Vi.T)), Vi)) -
                       f * T.dot(Vi, C) ))
    J, updates = theano.scan(lambda i, Vhat, Vk_ : T.grad(Vhat[i], Vk_), 
                         sequences=T.arange(Vhat.shape[0]), non_sequences=[Vhat, Vk_])
    num_fun = theano.function([Vk_], J, updates=updates)
    out_array = num_fun(V[k])
    return out_array


