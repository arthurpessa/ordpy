"""
A Python package for data analysis with permutation entropy 
and ordinal networks methods.

Copyright (C) 2020 Haroldo V. Ribeiro (hvr@dfi.uem.br)
                   Arthur A. B. Pessa (arthur_pessa@hotmail.com)

If you have used ordpy in a scientific publication, we would appreciate 
citations to the following reference:

A. A. B. Pessa, H. V. Ribeiro, "ordpy: A Python package for data 
analysis with permutation entropy and ordinal networks methods", 
Journal ??, ???? (2021).

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or 
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import numpy as np
import itertools


def np_setdiff(a, b):
    """
    Finds diferent elements (subarrays) in a that are not 
    contained in b. 
    
    From: 
    https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays
    
    Parameters
    ----------    
    a, b: tuple, list or array [[x_11, x_12, x_13, ..., x_1m], 
          object in the format  [x_21, x_22, x_23, ..., x_2m], ..., 
                                [x_n1, x_n2, x_n3, ..., x_nm]] (n x m)
    ----------
    Returns a np.ndarray containing the elements in a that are not 
    contained in b.
    
    Examples
    --------
    >>> a = np.asarray(((0,1,2), (0,1,2), (1,0,2), (2,0,1)))
    >>> b = np.asarray([[0,2,1], [0,1,2], [0,1,2]])    
    >>> np_setdiff(a, b)
    array([[1, 0, 2],
           [2, 0, 1]])
    """
    a = np.asarray(a)
    b = np.asarray(b)

    _, ncols = a.shape

    dtype={'names':['f{}'.format(i) for i in range(ncols)],
           'formats':ncols * [a.dtype]}

    C = np.setdiff1d(a.view(dtype), b.view(dtype))
    C = C.view(a.dtype).reshape(-1, ncols)

    return(C)


def symbolic_distribution(data, dx=3, dy=1, tau_x=1, tau_y=1, missing=False):
    """
    Applies the Bandt and Pompe symbolization process to extract the probability 
    distribution of ordinal patterns (permutations) from data.
    
    Parameters
    ----------
    data    : array object in the format [x_1, x_2, x_3, ..., x_n]
              or [[x_11, x_12, x_13, ..., x_1m], 
                 [x_21, x_22, x_23, ..., x_2m], ..., 
                 [x_n1, x_n2, x_n3, ..., x_nm]] (n x m);
    dx     : embedding dimension (horizontal axis) (default: 3).
    dy     : embedding dimension (vertical axis); must be 1 for time series 
             (default: 1).
    tau_x  : embedding delay (horizontal axis) (default: 1).
    tau_y  : embedding delay (vertical axis) (default: 1).
    missing: if True, permutations that do not appear in the symbolic sequence 
             obtained from data are shown; if False, they are ommited. 
             (default: False)
    ----------
    Returns a list of lists containing the occurring permutations and 
    their corresponding probabilities.

    Examples
    --------
    >>> symbolic_distribution([4,7,9,10,6,11,3], dx=2)
    [[[0, 1], [1, 0]], [0.6666666666666666, 0.3333333333333333]]
    >>>
    >>> symbolic_distribution([4,7,9,10,6,11,3], dx=3, missing=True)
    [[[0, 1, 2], [1, 0, 2], [2, 0, 1], [0, 2, 1], [1, 2, 0], [2, 1, 0]],
     [0.4, 0.2, 0.4, 0.0, 0.0, 0.0]]
    >>>
    >>> symbolic_distribution([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    [[[0, 1, 3, 2], [1, 0, 2, 3], [1, 2, 3, 0]], [0.5, 0.25, 0.25]]
    >>>
    >>> symbolic_distribution([[1,2,1,4],[8,3,4,5],[6,7,5,6]], 
                              dx=2, dy=2, tau_x=2)
    [[[0, 1, 3, 2], [0, 2, 1, 3], [1, 3, 2, 0]], [0.5, 0.25, 0.25]]
    
    References
    ----------
    For notation or context, please see: Citar o paper aqui.
    """
    try:
        ny, nx = np.shape(data)
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])
        
    partitions = np.concatenate(
        [
            [np.concatenate(data[j:j+dy*tau_y:tau_y,i:i+dx*tau_x:tau_x]) for i in range(nx-(dx-1)*tau_x)] 
            for j in range(ny-(dy-1)*tau_y)
        ]
    )

    symbols = np.apply_along_axis(np.argsort, 1, partitions)
    symbols, symbols_count = np.unique(symbols, return_counts=True, axis=0)

    probabilities          = symbols_count/len(partitions)

    if missing==False:
        return [symbols, probabilities]
    
    else:
        all_symbols   = list(map(list,list(itertools.permutations(np.arange(dx*dy)))))
        miss_symbols  = np_setdiff(all_symbols, symbols)
        symbols       = np.concatenate((symbols, miss_symbols))
        probabilities = np.concatenate((probabilities, np.zeros(miss_symbols.__len__())))
        
        return [symbols, probabilities]


def permutation_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1, base='e', normalized=True, probs=False):
    """
    Calculates Shannon's entropy using the symbolic ditribution extracted from
    data.
    
    Parameters
    ----------
    data      : data array object in the format [x_1, x_2, x_3, ..., x_n] 
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m)
                or the ordinal probabilities obtained from symbolic_distribution
                (if probs=True).
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    base      : logarithm base in Shannon's entropy. Either 'e' or 2 (default: 'e').
    normalized: if True, permutation entropy lies in the interval [0,1]; if False, 
                it does not (default: True).
    probs     : if True, assumes data imput to be an ordinal probability distribution.
                (default: False)

    ----------
    Returns the value of permutation entropy.
    
    Examples
    --------
    >>> permutation_entropy([4,7,9,10,6,11,3], dx=2, base=2, normalized=True)
    0.9182958340544896
    >>>
    >>> permutation_entropy([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, base=2, normalized=True)
    0.32715643797829735
    >>>
    >>> permutation_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, tau_x=2, normalized=False)
    1.0397207708399179
    
    References
    ----------
    For notation or context, please see: Citar o paper aqui.
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, missing=False)
    else:
        probabilities = data

    if normalized==True and base in [2, '2']:        
        smax = np.log2(np.math.factorial(dx*dy))
        s    = -np.sum(probabilities*np.log2(probabilities))

        return s/smax
    elif normalized==True and base=='e':        
        smax = np.log(np.math.factorial(dx*dy))
        s    = -np.sum(probabilities*np.log(probabilities))
        return s/smax
    elif normalized==False and base in [2, '2']:
        return -np.sum(probabilities*np.log2(probabilities))
    else:
        return -np.sum(probabilities*np.log(probabilities))


def complexity_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1):
    """
    Calculates permutation entropy and statistical complexity using 
    the symbolic ditribution extracted from data.
    
    Parameters
    ----------
    data      : tuple, list or array object in the format [x_1, x_2, x_3, ..., x_n] 
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m);
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    ----------
    Returns the values of normalized permutation entropy and statistical complexity.
    
    Examples
    --------
    >>> complexity_entropy([4,7,9,10,6,11,3], dx=2)
    [0.9182958340544894, 0.06112816548804511]
    >>>
    >>> complexity_entropy([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    [0.3271564379782973, 0.2701200547320647]
    >>>
    >>> complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]],dx=3, dy=2)
    [0.21070701155008006, 0.20704765093242872]
    
    References
    ----------
    For notation or context, please see: Citar o paper aqui.    
    """
    _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, missing=True)   
    h                = permutation_entropy(probabilities[probabilities>0], dx, dy, tau_x, tau_y, probs=True)

    n            = np.math.factorial(dx*dy)
    uniform_dist = np.full(n, 1/n)

    p_plus_u_over_2      = (uniform_dist + probabilities)/2  
    s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

    probabilities = probabilities[probabilities!=0]
    s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
    s_of_u_over_2 = np.log(n)/2.

    js_div_max = -0.5*(((n+1)/n)*np.log(n+1) + np.log(n) - 2*np.log(2*n))    
    js_div     = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

    return [h, h*js_div/js_div_max]


def logq(x, q=1):
    """
    Calculates the q-logarithm of x.

    Parameters
    ----------
    x: real number or array containing real numbers.
    q: Tsallis q parameter (default: 1).
    ----------
    Returns the value/array of the 
    q-logarithm of x.
    
    Examples
    --------
    >>> logq(np.math.e)
    1.0
    >>>    
    >>> logq([np.math.e for i in range(5)])
    array([1., 1., 1., 1., 1.])
    
    References
    ----------
    For notation or context, please see: Citar o paper aqui.     
    """
    x = np.asarray(x, dtype=float)

    if q==1:
        return np.log(x) 
    else:       
        return (x**(1-q) - 1)/(1-q)


def tsallis_entropy(data, q=1, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates the normalized Tsallis's entropy using the symbolic ditribution 
    extracted from data.
    
    Parameters
    ----------
    data      : tuple, list or array object in the format [x_1, x_2, x_3, ..., x_n]
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m)
                or the ordinal probabilities obtained from symbolic_distribution
                (if probs=True).
    q         : Tsallis q parameter, it can be an array of values (default: 1).
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    probs     : if True, assumes data imput to be an ordinal probability distribution.
                (default: False)
    ----------
    Returns the normalized values of Tsallis's entropy for each parameter q.
    
    Examples
    --------
    >>> tsallis_entropy([4,7,9,10,6,11,3], dx=2) 
    0.9182958340544894
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=[1,2], dx=2) 
    array([0.91829583, 0.88888889])
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=2, dx=2)
    0.888888888888889
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=2, dx=3)
    0.768

    References
    ----------
    For notation or context, please see: Physical Review E 95, 062106 (2017).    
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y)
    else:
        probabilities = data

    if isinstance(q,(list,np.ndarray)):
        s = []
        for q_ in q:
            smax          = logq(np.math.factorial(dx*dy),q_)
            lnq_1_over_p  = logq(1./probabilities,q_)
            s             += [np.sum(probabilities*lnq_1_over_p)/smax]
        s = np.asarray(s)
    else:
        smax          = logq(np.math.factorial(dx*dy),q)
        lnq_1_over_p  = logq(1./probabilities,q)
        s             = np.sum(probabilities*lnq_1_over_p)/smax

    return s


def jensen_tsallis_divergence_max(n_states, q):
    """
    Returns the maximum value of the Jensen Tsallis divergence.

    Parameters
    ----------
    n_states: number of ordinal states.
    q: Tsallis q parameter.
    ----------
    Returns the maximum divergence value.
    """
    if q==1:
        return(-0.5*(((n_states+1)/n_states)*np.log(n_states+1) + np.log(n_states) - 2*np.log(2*n_states)))
    else:
        # Equation 11 of Physical Review E 95, 062106 (2017).
        return(
               ((2**(2-q))*n_states - (1+n_states)**(1-q) - n_states*(1+1/n_states)**(1-q) - n_states + 1)/ 
               ((1-q)*(2**(2-q))*n_states)
              )


def tsallis_complexity_entropy(data, q=1, dx=3, dy=1, tau_x=1, tau_y=1):
    """
    Calculates Tsallis generalized complexity-entropy plane using the symbolic 
    ditribution extracted from data.
    
    Parameters
    ----------
    data      : tuple, list or array object in the format [x_1, x_2, x_3, ..., x_n] 
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m);
    q         : Tsallis q parameter, it can be an array of values (default: 1).
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    ----------
    Returns the values of Tsallis's generalized normalized permutation entropy 
    and Tsallis's generalized statistical complexity for each parameter q.

    Examples
    --------
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], dx=2, q=[1,2])
    array([[0.91829583, 0.06112817],
       [0.88888889, 0.07619048]])
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], q=2, dx=2)
    array([0.88888889, 0.07619048])
    >>> tsallis_complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], q=3, dx=3, dy=2)
    array([0.93750181, 0.92972165])
    
    References
    ----------
    For notation or context, please see: Physical Review E 95, 062106 (2017).
    """

    _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, missing=True)
    h_q              = tsallis_entropy(probabilities[probabilities>0], q, dx, dy, 
                                       tau_x, tau_y, probs=True)

    p             = probabilities[probabilities!=0]
    n             = np.math.factorial(dx*dy)
    uniform_dist  = np.full(n, 1/n)

    if isinstance(q, (list,np.ndarray)):
        jt_div     = []
        jt_div_max = []
        for q_ in q:
            #first and second terms of Jensen-Tsallis divergence: 
            #equation 10 of Physical Review E 95, 062106 (2017).
            first_term  = (uniform_dist[:len(p)] + p)/(2*p)
            first_term  = -0.5*np.sum(p*logq(first_term, q_)) 

            second_term = n*(uniform_dist + probabilities)/2
            second_term = -(0.5/n)*np.sum(logq(second_term, q_))

            jt_div      += [first_term + second_term]
            jt_div_max  += [jensen_tsallis_divergence_max(n,q_)]

        jt_div     = np.asarray(jt_div)
        jt_div_max = np.asarray(jt_div_max)
    else:
        first_term  = (uniform_dist[:len(p)] + p)/(2*p)
        first_term  = -0.5*np.sum(p*logq(first_term, q)) 

        second_term = n*(uniform_dist + probabilities)/2
        second_term = -(0.5/n)*np.sum(logq(second_term, q))

        jt_div      = first_term + second_term
        jt_div_max  = jensen_tsallis_divergence_max(n,q)


    return np.asarray([h_q, h_q*jt_div/jt_div_max]).T


def renyi_entropy(data, alpha=1, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates the normalized Rényi's entropy using the symbolic ditribution 
    extracted from data.
    
    Parameters
    ----------
    data      : tuple, list or array object in the format [x_1, x_2, x_3, ..., x_n]
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m)
                or the ordinal probabilities obtained from symbolic_distribution
                (if probs=True).
    alpha     : Rényi alpha parameter, it can be an array of values (default: 1).
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    probs     : if True, assumes data imput to be an ordinal probability distribution.
                (default: False)
    ----------
    Returns the normalized value of Rényi's entropy for each parameter alpha.
    
    Examples
    --------
    >>> renyi_entropy([4,7,9,10,6,11,3], dx=2)
    0.9182958340544894
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=2, dx=2)
    0.84799690655495
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=[1,2], dx=2)
    array([0.91829583, 0.84799691])
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=2, dx=3)
    0.5701944178769374

    References
    ----------
    For notation or context, please see: Citar o paper aqui.    
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y)
    else:
        probabilities = data

    smax = np.log(np.math.factorial(dx*dy))

    if isinstance(alpha, (list,np.ndarray)):
        s = []
        for alpha_ in alpha:
            if alpha_ !=1:
                s += [(1/(1-alpha_))*np.log(np.sum(probabilities**alpha_))/smax]
            else:
                s += [-np.sum(probabilities*np.log(probabilities))/smax]
        s = np.asarray(s)
    else:
        if alpha !=1:
            s = (1/(1-alpha))*np.log(np.sum(probabilities**alpha))/smax
        else:
            s = -np.sum(probabilities*np.log(probabilities))/smax

    return s


def jensen_renyi_divergence_max(n_states, q):
    """
    Returns the maximum value of the Jensen Renyi divergence.
    
    Parameters
    ----------
    n_states: number of ordinal states.
    q: Renyi q parameter.
    ----------
    Returns the maximum divergence value.
    """
    if q==1:
        return(-0.5*(((n_states+1)/n_states)*np.log(n_states+1) + np.log(n_states) - 2*np.log(2*n_states)))
    else:
        # Equation 5 in Physica A 498 (2018) 74-85.
        return(
                (np.log(((n_states + 1.)**(1. - q) + n_states - 1.)/(2.**(1. - q)*n_states)) + (1. - q)*
                 np.log((n_states + 1.)/(2.*n_states)))/(2.*(q - 1))
              )


def renyi_complexity_entropy(data, alpha=1, dx=3, dy=1, tau_x=1, tau_y=1):
    """
    Calculates Rényi generalized complexity-entropy plane using a symbolic 
    distribution obtained from data.
    
    Parameters
    ----------
    data      : tuple, list or array object in the format [x_1, x_2, x_3, ..., x_n] 
                or [[x_11, x_12, x_13, ..., x_1m], 
                    [x_21, x_22, x_23, ..., x_2m], ..., 
                    [x_n1, x_n2, x_n3, ..., x_nm]] (n x m);
    alpha     : Rényi alpha parameter (default: 1).
    dx        : embedding dimension (horizontal axis) (default: 3).
    dy        : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x     : embedding delay (horizontal axis) (default: 1).
    tau_y     : embedding delay (vertical axis) (default: 1).
    ----------

    Returns the values of Rényi's generalized normalized permutation entropy 
    and Rényi's generalized statistical complexity for each parameter alpha.

    Examples
    --------
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], alpha=2, dx=2)
    array([0.84799691, 0.08303895])
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], alpha=[1, 2], dx=2)
    array([[0.91829583, 0.06112817],
       [0.84799691, 0.08303895]])
    >>> renyi_complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], alpha=3, dx=3, dy=2)
    array([0.21070701, 0.20975673])

    References
    ----------
    For notation or context, please see: Citar o paper aqui.    
    """

    _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, missing=True)
    h_a              = renyi_entropy(probabilities, alpha, dx, dy, tau_x, tau_y, probs=True)

    n             = np.math.factorial(dx*dy)
    uniform_dist  = np.full(n, 1/n)

    if isinstance(alpha, (list,np.ndarray)):
        jr_div = []
        jr_div_max = []
        for alpha_ in alpha:
            if alpha_==1:
                p_plus_u_over_2      = (uniform_dist + probabilities)/2  
                s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

                probabilities = probabilities[probabilities!=0]
                s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
                s_of_u_over_2 = np.log(n)/2.

                jr_div       += [s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2]
            else:
                # Equation 4 in Physica A 498 (2018) 74-85.
                first_term  = ((probabilities + uniform_dist)/2)**(1-alpha_)
                second_term = (1/(n**alpha_))*first_term

                first_term  = np.log(np.sum(first_term*probabilities**alpha_))
                second_term = np.log(np.sum(second_term))

                jr_div     += [(1/(2*(alpha_-1)))*(first_term + second_term)]

            jr_div_max += [jensen_renyi_divergence_max(n, alpha_)]

        jr_div = np.asarray(jr_div)
        jr_div_max = np.asarray(jr_div_max)
    else:
        if alpha==1:
            p_plus_u_over_2      = (uniform_dist + probabilities)/2  
            s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

            probabilities = probabilities[probabilities!=0]
            s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
            s_of_u_over_2 = np.log(n)/2.

            jr_div        = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2
        else:
            first_term  = ((probabilities + uniform_dist)/2)**(1-alpha)
            second_term = (1/(n**alpha))*first_term

            first_term  = np.log(np.sum(first_term*probabilities**alpha))
            second_term = np.log(np.sum(second_term))

            jr_div      = (1/(2*(alpha-1)))*(first_term + second_term)
            
        jr_div_max = jensen_renyi_divergence_max(n, alpha)

    return np.asarray([h_a,h_a*jr_div/jr_div_max ]).T


def ordinal_network(data, dx=3, dy=1, tau_x=1, tau_y=1, normalized=True, overlapping=True, connections="all"):
    """
    Generates the elements (nodes, edges and edge weights) necessary to 
    obtain an ordinal network from data.
    
    Parameters
    ----------
    data        : array object in the format [x_1, x_2, x_3, ..., x_n] 
                  or [[x_11, x_12, x_13, ..., x_1m], 
                      [x_21, x_22, x_23, ..., x_2m], ..., 
                      [x_n1, x_n2, x_n3, ..., x_nm]] (nxm);
    dx          : embedding dimension (horizontal axis) (default: 3).
    dy          : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x       : embedding delay (horizontal axis) (default: 1).
    tau_y       : embedding delay (vertical axis) (default: 1).
    overlapping : if True, partitions data into overlapping sliding windows; if False, it does not. 
                  (Also notice that definitions of tau_x and tau_y become somewhat strange if 
                  partitions do not overlap. So tau_x = tau_y = 1 if False). Parameter only
                  valid for time series data. (default: True). 
    normalized  : if True, edge weights represent transition probabilities between permutations.
                  If False, edge weights are transition counts.
    connections : ordinal network is constructed using "all" permutation successions in symbolic 
                  sequence or only "horizontal" or "vertical" successions. Parameter only valid 
                  for image data (default: "all"). 
    ----------
    Returns a list of arrays containing the nodes, edges and edge weights
    corresponding to an ordinal network mapped from data.
    
    Examples
    --------
    >>> ordinal_network([4,7,9,10,6,11,8,3,7], dx=2, normalized=False)
    (array(['0|1', '1|0'], dtype='<U3'),
     array([['0|1', '0|1'],
            ['0|1', '1|0'],
            ['1|0', '0|1'],
            ['1|0', '1|0']], dtype='<U3'),
     array([2, 2, 2, 1]))
    >>>
    >>> ordinal_network([4,7,9,10,6,11,8,3,7], dx=2, overlapping=False, normalized=False)
    (array(['0|1', '1|0'], dtype='<U3'),
     array([['0|1', '0|1'],
            ['0|1', '1|0']], dtype='<U3'),
     array([2, 1]))
    >>>
    >>> ordinal_network([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, normalized=False)
    (array(['0|1|3|2', '1|0|2|3', '1|2|3|0'], dtype='<U7'),
     array([['0|1|3|2', '1|0|2|3'],
            ['0|1|3|2', '1|2|3|0'],
            ['1|0|2|3', '0|1|3|2'],
            ['1|2|3|0', '0|1|3|2']], dtype='<U7'),
     array([1, 1, 1, 1]))
    >>> 
    >>> ordinal_network([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, normalized=False, connections='horizontal')
    (array(['0|1|3|2', '1|0|2|3', '1|2|3|0'], dtype='<U7'),
     array([['0|1|3|2', '1|0|2|3'],
            ['1|2|3|0', '0|1|3|2']], dtype='<U7'),
     array([1, 1]))
    
    References
    ----------
    For notation or context, please see: Citar o paper aqui.    
    """     
    try:
        ny, nx = np.shape(data)
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])

    #TIME SERIES DATA
    if ny==1:
        if overlapping == True:
            partitions = np.concatenate(
                [
                    [np.concatenate(data[j:j+dy*tau_y:tau_y,i:i+dx*tau_x:tau_x]) for i in range(nx-(dx-1)*tau_x)] 
                    for j in range(ny-(dy-1)*tau_y)
                ]
            )

            states = np.apply_along_axis(np.argsort, 1, partitions)
            links  = np.concatenate([states[i:(i-1) or None, None] for i in range(2)], axis=1)

            unique_links, occurrences = np.unique(links, return_counts=True, axis=0)
            if normalized == True: 
                occurrences = occurrences/len(links)

            unique_links = np.apply_along_axis(np.array2string, 2, unique_links, separator="|")
            for char in ['[', ']']:
                unique_links = np.char.replace(unique_links, char, '')   
            
            return (np.unique(unique_links), unique_links, occurrences)

        else:
            partitions = np.concatenate(
                [
                    [np.concatenate(data[j:j+dy*1:1, i:i+dx:1]) for i in range(0, nx-(dx-1), dx)] 
                    for j in range(ny-(dy-1))
                ]
            )
            
            states = np.apply_along_axis(np.argsort, 1, partitions)
            links  = np.concatenate([states[i:(i-1) or None, None] for i in range(2)], axis=1)

            unique_links, occurrences = np.unique(links, return_counts=True, axis=0)
            if normalized == True: 
                occurrences = occurrences/len(links)

            unique_links = np.apply_along_axis(np.array2string, 2, unique_links, separator="|")
            for char in ['[', ']']:
                unique_links = np.char.replace(unique_links, char, '')   

            return (np.unique(unique_links), unique_links, occurrences)

    #IMAGE DATA
    else:
        partitions = np.concatenate(
            [
                [[np.concatenate(data[j:j+dy*tau_y:tau_y,i:i+dx*tau_x:tau_x]) for i in range(nx-(dx-1)*tau_x)]]
                for j in range(ny-(dy-1)*tau_y)
            ]
        )
        #horizontal and vertical successions        
        if not connections in ["horizontal", "vertical"]:
            states = np.apply_along_axis(np.argsort, 2, partitions)
        #horizontal successions of permutation states
            hlinks = np.concatenate([[row[i:i+2]] for row in states for i in range(len(row)-1)])
        #vertical successions of permutation states
            states = states.transpose(1,0,2)
            vlinks = np.concatenate([[column[i:i+2]] for column in states for i in range(len(column)-1)])
        #all occurring transitions
            links  = np.concatenate((hlinks, vlinks))

            unique_links, occurrences = np.unique(links, return_counts=True, axis=0)
            if normalized == True: 
                occurrences = occurrences/len(links)

        #only horizontal successions of permutation states
        elif connections == 'horizontal':
            states = np.apply_along_axis(np.argsort, 2, partitions)
            hlinks = np.concatenate([[line[i:i+2]] for line in states for i in range(len(line)-1)], axis=0)
            
            unique_links, occurrences = np.unique(hlinks, return_counts=True, axis=0)
            if normalized == True: 
                occurrences = occurrences/len(hlinks)

        #only vertical successions of permutation states        
        else: #connections == 'vertical':
            states = np.apply_along_axis(np.argsort, 2, partitions)
            states = states.transpose(1,0,2)
            vlinks = np.concatenate([[line[i:i+2]] for line in states for i in range(len(line)-1)])
            
            unique_links, occurrences = np.unique(vlinks, return_counts=True, axis=0)
            if normalized == True: 
                occurrences = occurrences/len(vlinks)
                
        unique_links = np.apply_along_axis(np.array2string, 2, unique_links, separator="|")
        for char in ['[', ']']:
            unique_links = np.char.replace(unique_links, char, '')   
            
        return [np.unique(unique_links), unique_links, occurrences]


def global_network_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1, normalized=True, overlapping=True, connections="all"):
    """
    Calculates global node entropy from data.

    Parameters
    ----------
    data        : array object in the format [x_1, x_2, x_3, ..., x_n] 
                  or [[x_11, x_12, x_13, ..., x_1m], 
                      [x_21, x_22, x_23, ..., x_2m], ..., 
                      [x_n1, x_n2, x_n3, ..., x_nm]] (n x m).
                  or the ordinal network returned by ordinal_network().
    dy          : embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x       : embedding delay (horizontal axis) (default: 1).
    tau_y       : embedding delay (vertical axis) (default: 1).
    overlapping : if True, partitions data into overlapping sliding windows; if False, it does not. 
                  (Also notice that definitions of tau_x and tau_y become somewhat strange if 
                  partitions do not overlap. So tau_x = tau_y = 1 if False). Parameter only
                  valid for time series data. (default: True). 
    normalized  : if True, edge weights represent transition probabilities between permutations.
                  If False, edge weights are transition counts.
    connections : ordinal network is constructed using "all" permutation successions in symbolic 
                  sequence or only "horizontal" or "vertical" successions. Parameter only valid 
                  for image data (default: "all"). 
    ----------
    Returns the value of global node entropy.
    
    Examples
    --------
    
    """
    if len(data)==3 and type(data[0][0])==numpy.str_:
        nodes, links, weights = data
    else:
        nodes, links, weights = ordinal_network(data, dx, dy, tau_x, tau_y, normalized, overlapping, connections)

    links_source = links.transpose()[0]
    links_target = links.transpose()[1]
    h_gn  = 0
    for node in nodes:
        args           = np.argwhere(links_source==node).flatten()
        renorm_weights = weights[args]/np.sum(weights[args])

        args_in        = np.argwhere(links_target==node).flatten()
        p_in           = np.sum(weights[args_in])
        
        h_i            = -np.sum(renorm_weights*np.log2(renorm_weights))
        h_gn          += p_in*h_i
        
    return h_gn


def random_ordinal_network(dx=3, dy=1):
    """
    Creates a random ordinal network representing the mapping of a random time 
    series or a random bidimensional field. The result assume overlapping window
    partitions and embbeding delays unitaries.
    
    Parameters
    ----------
    dx     : horizontal embedding dimension (default: 3)
    dy     : vertical embedding dimension (default: 1 [time series data])
    ----------
    Returns a list of arrays containing the nodes, edges and edge weights
    corresponding to a random ordinal network.
    
    """

#THEORETICAL RESULTS FOR IMAGE DATA
    if dx>1 and dy>1:
        allowed_links     = []
        theoretical_probs = []

    #this loop has the job of finding all allowed transitions starting from an ordinal pattern and calculating 
    #the expected probability of such transitions for a random bidimensional scalar field.
        for pattern in it.permutations(np.arange(dx*dy).astype('int')):  
            pattern            = list(pattern)
            hlinks, vlinks     = pattern.copy(), pattern.copy()
            hnumbers, vnumbers = [], []            

            #horizontal moves
            for i in range(len(pattern)):
                if hlinks[i]%dx==0:
                    hnumbers.append(hlinks[i]+dx-1)
                    hlinks[i] = dx*dy+2
                else:
                    hlinks[i] = hlinks[i]-1 

            #vertical moves           
            for i in range(len(pattern)):
                if vlinks[i]<dx:
                    vnumbers.append(vlinks[i]+dx*(dy-1))
                    vlinks[i] = dx*dy+2
                else:
                    vlinks[i] = vlinks[i]-dx            

            hlinks, vlinks = [hlinks], [vlinks]     

            #allowed horizontal transitions
            for number, k in zip(hnumbers, range(len(hnumbers))):
                d = []
                for permutation in hlinks:
                    for index in range(len(permutation),-1,-1):
                        b = permutation.copy()
                        b.insert(index, number)
                        if k==len(hnumbers)-1:
                            for _ in range(len(hnumbers)):
                                b.remove(dx*dy+2)
                        d.append(b)
                    hlinks = d.copy()        

            #allowed vertical transitions
            for number, k in zip(vnumbers, range(len(vnumbers))):
                d = []
                for permutation in vlinks:
                    for index in range(len(permutation),-1,-1):
                        b = permutation.copy()
                        b.insert(index, number)
                        if k==len(vnumbers)-1:
                            for _ in range(len(vnumbers)):
                                b.remove(dx*dy+2)
                        d.append(b)   
                    vlinks = d.copy()                       

            unique_links, frequencies = np.unique(hlinks+vlinks, axis=0, return_counts=True)
            links                     = np.concatenate((np.tile(pattern,len(unique_links)).reshape(-1,dx*dy),
                                                        unique_links), axis=1).reshape(len(unique_links),2,-1)
            frequencies               = frequencies/frequencies.sum()/np.math.factorial(dx*dy)

    #transforming lists of permutations to strings to save the transitions        
            links                     = np.apply_along_axis(np.array2string, 2, links, separator='|')
            for char in ['[', ']']:
                links = np.char.replace(links, char, '')    

            allowed_links.append(links.tolist())
            theoretical_probs.append(frequencies.tolist())

        edge_array, weight_array = np.concatenate(allowed_links), np.concatenate(theoretical_probs)
        vertices_names           = np.unique(edge_array)

        return [vertices_names, edge_array, weight_array]
    
#THEORETICAL RESULTS FOR TIME SERIES DATA
    else:
        vertices_names = []
        edge_array = []   
        weight_array    = []

        for j in it.permutations(np.arange(dx)):
            sub_perm = np.asarray(j)-1
            sub_perm = np.extract(sub_perm > -1, sub_perm)
            source   = '|'.join(map(str, j))
            index    = source.find('0')

            for k in range(len(sub_perm), -1, -1):
                target = '|'.join(map(str, np.insert(sub_perm, k, dx-1)))

                if not target.find(str(dx-1)) == index:
                    weight_array.append(1/((dx+1)*np.math.factorial(dx)))
                else:
                    weight_array.append(2/((dx+1)*np.math.factorial(dx)))

                edge_array.append([source, target])
        
        vertices_names = np.unique(edge_array)    

        return [vertices_names, np.asarray(edge_array), np.asarray(weight_array)]