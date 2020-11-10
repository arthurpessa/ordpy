"""
ordpy: A Python Package for Data Analysis with Permutation Entropy and Ordinal Networks Methods
===============================================================================================

``ordpy`` is pure Python module\\ [#pessa2021]_ that implements data analysis methods based
on Band and Pompe\\ [#bandt_pompe]_ symbolic encoding scheme.

.. note::

   If you have used ``ordpy`` in a scientific publication, we would appreciate 
   citations to the following reference\\ [#pessa2021]_:

   - A. A. B. Pessa, H. V. Ribeiro, `ordpy: A Python package for data 
     analysis with permutation entropy and ordinal networks methods 
     <https://ourpaper_url>`_, ????, ??-?? (2021).

    .. code-block:: bibtex

       @article{pessa2021ordpy,
         title    = {ordpy: A Python Package for Data Analysis with Permutation 
                     Entropy and Ordinal Networks Methods},
         author   = {Pessa, Arthur A. B. and Ribeiro, Haroldo V.},
         journal  = {?},
         volume   = {?},
         number   = {?},
         pages    = {?},
         year     = {2021},
         doi      = {?}
       }

``ordpy`` implements the following data analysis methods:

- permutation entropy for time series\\ [#bandt_pompe]_ and images\\ [#ribeiro_2012]_;
- complexity-entropy plane for time series\\ [#lopezruiz]_\\ :sup:`,`\\ [#rosso]_ and 
  images\\ [#ribeiro_2012]_;
- multiscale complexity-entropy plane for time series\\ [#zunino2012]_ and 
  images\\ [#zunino2016]_;
- Tsallis\\ [#ribeiro2017]_ and Rényi\\ [#jauregui]_ generalized complexity-entropy
  curves for time series and images;
- ordinal networks for time series\\ [#small]_\\ :sup:`,`\\ [#pessa2019]_ and 
  images\\ [#pessa2020]_;
- global node entropy of ordinal networks for 
  time series\\ [#McCullough]_\\ :sup:`,`\\ [#pessa2019]_ and images\\ [#pessa2020]_.

References
----------

.. [#pessa2021] Pessa, A. A., & Ribeiro, H. V. (2020). ordpy: A Python package
   for data analysis with permutation entropy and ordinal networks methods. 
   arXiv preprint arXiv:2007.03090.

.. [#bandt_pompe] Bandt, C., & Pompe, B. (2002). Permutation entropy: A Natural 
   Complexity Measure for Time Series. Physical Review Letters, 88, 174102.

.. [#ribeiro_2012] Ribeiro, H. V., Zunino, L., Lenzi, E. K., Santoro, P. A., &
   Mendes, R. S. (2012). Complexity-Entropy Causality Plane as a Complexity
   Measure for Two-Dimensional Patterns. PLOS ONE, 7, e40689.

.. [#lopezruiz] Lopez-Ruiz, R., Mancini, H. L., & Calbet, X. (1995). A Statistical
   Measure of Complexity. Physics Letters A, 209, 321-326.

.. [#rosso] Rosso, O. A., Larrondo, H. A., Martin, M. T., Plastino, A., &
   Fuentes, M. A. (2007). Distinguishing Noise from Chaos. Physical Review 
   Letters, 99, 154102.

.. [#zunino2012] Zunino, L., Soriano, M. C., & Rosso, O. A. (2012). 
   Distinguishing Chaotic and Stochastic Dynamics from Time Series by Using 
   a Multiscale Symbolic Approach. Physical Review E, 86, 046210.

.. [#zunino2016] Zunino, L., & Ribeiro, H. V. (2016). Discriminating Image 
   Textures with the Multiscale Two-Dimensional Complexity-Entropy Causality 
   Plane. Chaos, Solitons & Fractals, 91, 679-688.

.. [#ribeiro2017] Ribeiro, H. V., Jauregui, M., Zunino, L., & Lenzi, E. K. 
   (2017). Characterizing Time Series Via Complexity-Entropy Curves. 
   Physical Review E, 95, 062106.

.. [#jauregui] Jauregui, M., Zunino, L., Lenzi, E. K., Mendes, R. S., &
   Ribeiro, H. V. (2018). Characterization of Time Series via Rényi 
   Complexity-Entropy Curves. Physica A, 498, 74-85.

.. [#small] Small, M. (2013). Complex Networks From Time Series: Capturing 
   Dynamics. In 2013 IEEE International Symposium on Circuits and Systems
   (ISCAS2013) (pp. 2509-2512). IEEE.

.. [#pessa2019] Pessa, A. A., & Ribeiro, H. V. (2019). Characterizing Stochastic 
   Time Series With Ordinal Networks. Physical Review E, 100, 042304.

.. [#pessa2020] Pessa, A. A., & Ribeiro, H. V. (2020). Mapping Images Into
   Ordinal Networks. arXiv preprint arXiv:2007.03090.

.. [#McCullough] McCullough, M., Small, M., Iu, H. H. C., & Stemler, T. (2017).
   Multiscale Ordinal Network Analysis of Human Cardiac Dynamics.
   Philosophical Transactions of the Royal Society A, 375, 20160292.

.. [#amigó] Amigó, J. M., Zambrano, S., & Sanjuán, M. A. F. (2007).
   True and False Forbidden Patterns in Deterministic and Random Dynamics.
   Europhysics Letters, 79, 50001.

.. [#rosso_curvas] Martin, M. T., Plastino, A., & Rosso, O. A. (2006). 
   Generalized Statistical Complexity Measures: Geometrical and 
   Analytical Properties, Physica A, 369, 439–462.


Installing
==========

To install the ordpy use

.. code-block:: console

   pip install ordpy

or

.. code-block:: console

   git clone https://gitlab.com/hvribeiro/ordpy.git
   cd ordpy
   pip install -e .


Basic usage
===========

We provide a `notebook <https://github.com/hvribeiro/ordpy/blob/master/examples/sample_notebook.ipynb>`_
illustrating how to use ``ordpy``. This notebook reproduces all figures of our
article\\ [#pessa2021]_. The codes below show simple usages of ``ordpy``.

**Complexity-entropy plane for logistic map and Gaussian noise**

.. code-block:: python
   
    import numpy as np
    import ordpy
    from matplotlib import pylab as plt

    def logistic(a=4, n=100000, x0=0.4):
        x = np.zeros(n)
        x[0] = x0
        for i in range(n-1):
            x[i+1] = a*x[i]*(1-x[i])
        return(x)

    time_series = [logistic(a) for a in [3.05, 3.55, 4]]
    time_series += [np.random.normal(size=100000)]

    HC = [ordpy.complexity_entropy(series, dx=4) for series in time_series]


    f, ax = plt.subplots(figsize=(9.1,7))

    for HC_, label_ in zip(HC, ['Simple periodic (a=3.05)', 
                                '4-period (a=3.55)', 
                                'Chaotic (a=4)', 
                                'Gaussian noise']):
        ax.scatter(*HC_, label=label_, s=100)
        
    ax.set_xlabel('Permutation entropy, $H$')
    ax.set_ylabel('Statistical complexity, $C$')

    plt.legend()

.. figure:: ../examples/figs/sample_fig.png
   :height: 489px
   :width: 633px
   :scale: 80 %
   :align: center

**Ordinal networks for logistic map and Gaussian noise**

.. code-block:: python

    import numpy as np
    import igraph
    import ordpy
    from matplotlib import pylab as plt

    vertex_list, edge_list, edge_weight_list = list(), list(), list()

    for series in time_series:
        v_, e_, w_ = ordpy.ordinal_network(series, dx=4)
        vertex_list += [v_]
        edge_list += [e_]
        edge_weight_list += [w_]

    def create_ig_graph(vertex_list, edge_list, edge_weight):
        
        G = igraph.Graph(directed=True)
        
        for v_ in vertex_list:
            G.add_vertex(v_)
        
        for [in_, out_], weight_ in zip(edge_list, edge_weight):
            G.add_edge(in_, out_, weight=weight_)
            
        return G

    graphs = []

    for v_, e_, w_ in zip(vertex_list, edge_list, edge_weight_list):
        graphs += [create_ig_graph(v_, e_, w_)]

    def igplot(g):
        f = igraph.plot(g,
                        layout=g.layout_circle(),
                        bbox=(500,500),
                        margin=(40, 40, 40, 40),
                        vertex_label = [s.replace('|','') for s in g.vs['name']],
                        vertex_label_color='#202020',
                        vertex_color='#969696',
                        vertex_size=20,
                        vertex_font_size=6,
                        edge_width=(1 + 8*np.asarray(g.es['weight'])).tolist(),
                       )
        return f

    from IPython.core.display import display, SVG

    for graph_, label_ in zip(graphs, ['Simple periodic (a=3.05)', 
                                       '4-period (a=3.55)', 
                                       'Chaotic (a=4)', 
                                       'Gaussian noise']):
        print(label_)
        display(SVG(igplot(graph_)._repr_svg_()))

.. figure:: ../examples/figs/sample_net.png
   :height: 1648px
   :width: 795px
   :scale: 50 %
   :align: center    




List of functions
=================
"""

import numpy as np
import itertools


#(Arthur): changed the name of the function from np_setdiff to setdiff. (18/08/2020)
def setdiff(a, b):
    """
    Searches for elements (subarrays) in `a` that are not contained in `b` [*]_. 

    Parameters
    ----------    
    a : tuples, lists or arrays
        Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
        \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]` (:math:`n \\times m`).
    b : tuples, lists or arrays
        Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
        \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]` (:math:`n \\times m`).
    
    Returns
    -------
      : array
        An array containing the elements in `a` that are not contained in `b`.

    Notes
    -----
    .. [*] This function was adapted from https://stackoverflow.com/questions/8317022/get-intersecting-rows-across-two-2d-numpy-arrays

    Examples
    --------
    >>> a = ((0,1,2), (0,1,2), (1,0,2), (2,0,1))
    >>> b = [[0,2,1], [0,1,2], [0,1,2]] 
    >>> setdiff(a, b)
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


def symbolic_distribution(data, dx=3, dy=1, tau_x=1, tau_y=1, return_missing=False):
    """
    Applies the Bandt and Pompe\\ [#bandt_pompe]_ symbolization process to extract a
    probability distribution of ordinal patterns (permutations) from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    return_missing: boolean
                    If `True`, permutations that do not appear in the symbolic sequence 
                    obtained from data are shown. If `False`, they are omitted. 
                    (default: `False`)

    Returns
    -------
     : tuple
       A tuple of arrays containing the permutations occurring in data 
       and their corresponding probabilities.

    Examples
    --------
    >>> symbolic_distribution([4,7,9,10,6,11,3], dx=2)
    (array([[0, 1],
            [1, 0]]),
     array([0.66666667, 0.33333333]))
    >>>
    >>> symbolic_distribution([4,7,9,10,6,11,3], dx=3, return_missing=True)
    (array([[0, 1, 2],
            [1, 0, 2],
            [2, 0, 1],
            [0, 2, 1],
            [1, 2, 0],
            [2, 1, 0]]),
     array([0.4, 0.2, 0.4, 0. , 0. , 0. ]))
    >>>
    >>> symbolic_distribution([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    (array([[0, 1, 3, 2],
            [1, 0, 2, 3],
            [1, 2, 3, 0]]),
     array([0.5 , 0.25, 0.25]))
    >>>
    >>> symbolic_distribution([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, tau_x=2)
    (array([[0, 1, 3, 2],
            [0, 2, 1, 3],
            [1, 3, 2, 0]]),
     array([0.5 , 0.25, 0.25]))
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

    probabilities = symbols_count/len(partitions)

    if return_missing==False:
        return symbols, probabilities
    
    else:
        all_symbols   = list(map(list,list(itertools.permutations(np.arange(dx*dy)))))
        miss_symbols  = setdiff(all_symbols, symbols)
        symbols       = np.concatenate((symbols, miss_symbols))
        probabilities = np.concatenate((probabilities, np.zeros(miss_symbols.__len__())))
        
        return symbols, probabilities


def permutation_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1, base='e', normalized=True, probs=False):
    """
    Calculates Shannon's entropy using a symbolic ditribution extracted from
    data\\ [#bandt_pompe]_\\ :sup:`,`\\ [#ribeiro_2012]_.
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.symbolic_distribution`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or the ordinal probabilities returned by :func:`ordpy.symbolic_distribution`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
             Embedding delay (vertical axis) (default: 1).
    base : str, int
           Logarithm base in Shannon's entropy. Either 'e' or 2 (default: 'e').
    normalized: boolean
                If `True`, permutation entropy is normalized by its maximum value.
                If `False`, it does not (default: `True`).
    probs : boolean
            If `True`, assumes data input to be an ordinal probability
            distribution (default: `False`). 

    Returns
    -------
     : float
       The value of permutation entropy.
    
    Examples
    --------
    >>> permutation_entropy([4,7,9,10,6,11,3], dx=2, base=2, normalized=True)
    0.9182958340544896
    >>>
    >>> permutation_entropy([.5,.5], dx=2, base=2, normalized=False, probs=True)
    1.0
    >>>
    >>> permutation_entropy([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, base=2, normalized=True)
    0.32715643797829735
    >>>
    >>> permutation_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, tau_x=2, normalized=False)
    1.0397207708399179
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=False)
#(Arthur): added a convertion to array and a filter to remove 0's in the probability distribution
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]

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


#(Arthur): added the 'probs' parameter. (24/08/2020)
def complexity_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates permutation entropy\\ [#bandt_pompe]_ and statistical
    complexity\\ [#lopezruiz]_ (measures which define the complexity-entropy 
    causality plane [#rosso]_\\ :sup:`,`\\ [#ribeiro_2012]_) using a 
    symbolic ditribution obtained from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    probs: boolean
           If `True`, assumes the data input to be an ordinal probability distribution
           (default: `False`)
    
    Returns
    -------
     : tuple
       Values of the normalized permutation entropy and statistical complexity.
    
    Examples
    --------
    >>> complexity_entropy([4,7,9,10,6,11,3], dx=2)
    (0.9182958340544894, 0.06112816548804511)
    >>>
    >>> p = symbolic_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
    >>> complexity_entropy(p, dx=2, probs=True)
    (0.9182958340544894, 0.06112816548804511)
    >>> 
    >>> complexity_entropy([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    (0.3271564379782973, 0.2701200547320647)
    >>>
    >>> complexity_entropy([1/3, 1/15, 4/15, 2/15, 1/5, 0], dx=3, probs=True)
    (0.8314454838586238, 0.16576716623440763)
    >>>
    >>> complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]],dx=3, dy=2)
    (0.21070701155008006, 0.20704765093242872)
    """
#checking if 'data' is a probability distribution or not
    if probs==False:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=True)   
        h                = permutation_entropy(probabilities[probabilities>0], dx, dy, tau_x, tau_y, probs=True)
    else:
        if len(data)==np.math.factorial(dx*dy):
            probabilities = np.asarray(data)
            h             = permutation_entropy(probabilities[probabilities>0], dx, dy, tau_x, tau_y, probs=True)
        else:
            raise Exception("The data provided as input is not adequate. The length of" \
                            " the ordinal distribution array does not match the number of" \
                            " possible permutations for dx = {} and dy = {}." \
                            " Please run symbolic_function(..., missing=True) instead of" \
                            " symbolic_function(..., missing=False) to use part of its" \
                            " return as input in complexity_entropy() or change the" \
                            " values of the embedding dimensions.".format(dx, dy))

    n            = np.math.factorial(dx*dy)
    uniform_dist = np.full(n, 1/n)

    p_plus_u_over_2      = (uniform_dist + probabilities)/2  
    s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

    probabilities = probabilities[probabilities!=0]
    s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
    s_of_u_over_2 = np.log(n)/2.

    js_div_max = -0.5*(((n+1)/n)*np.log(n+1) + np.log(n) - 2*np.log(2*n))    
    js_div     = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

    return h, h*js_div/js_div_max


def logq(x, q=1):
    """
    Calculates the `q`-logarithm of `x`.

    Parameters
    ----------
    x : float or array
        Real number or array containing real numbers.
    q : float
        Tsallis `q` parameter (default: 1).

    Returns
    -------
     : float or array
       Value or array of values with the `q`-logarithm of `x`.

    Notes
    -----
    The `q`-logarithm of `x` is difined as [*]_

    .. math::

       \\log_q (x) = \\frac{x^{1-q} - 1}{1-q}~~\\text{for}~~q\\neq 1

    and :math:`\\log_q (x) = \\log (x)` for :math:`q=1`.

    .. [*] Tsallis, C. (2009). Introduction to Nonextensive Statistical 
       Mechanics: Approaching a Complex World. Springer.

    Examples
    --------
    >>> logq(np.math.e)
    1.0
    >>>
    >>> logq([np.math.e for i in range(5)])
    array([1., 1., 1., 1., 1.])
    """
    x = np.asarray(x, dtype=float)

    if q==1:
        return np.log(x) 
    else:       
        return (x**(1-q) - 1)/(1-q)


def tsallis_entropy(data, q=1, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates the normalized Tsallis's entropy\\ [#ribeiro2017]) using a
    symbolic ditribution obtained from data.
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.symbolic_distribution`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`) or the ordinal probabilities returned
           by :func:`ordpy.symbolic_distribution`.
    q : float or array
        Tsallis `q` parameter; it can be an array of values (default: 1).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    probs: boolean
           If `True`, assumes the data input to be an ordinal probability distribution
           (default: `False`)

    Returns
    -------
     : float or array
       The normalized values of Tsallis's entropy for each parameter `q`.
    
    Examples
    --------
    >>> tsallis_entropy([4,7,9,10,6,11,3], dx=2) 
    0.9182958340544894
    >>>
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=[1,2], dx=2) 
    array([0.91829583, 0.88888889])
    >>>
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=2, dx=2)
    0.888888888888889
    >>>
    >>> tsallis_entropy([1/3, 1/15, 4/15, 2/15, 1/5, 0], dx=3, probs=True)
    0.8314454838586238
    >>>
    >>> tsallis_entropy([4,7,9,10,6,11,3], q=2, dx=3)
    0.768 
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y)
    else:
#(Arthur): added a convertion to array and a filter to remove 0's in the probability distribution
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]

#(Arthur): changed to add tuple. (19/08/2020)
    if isinstance(q, (tuple, list, np.ndarray)):
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


#(Arthur): added the 'probs' parameter. (24/08/2020)
def tsallis_complexity_entropy(data, q=1, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates permutation entropy\\ [#bandt_pompe]_ and statistical
    complexity\\ [#lopezruiz]_ using Tsallis's entropy (these measures define 
    a generalized complexity-entropy causality plane\\ [#ribeiro2017]_ ) 
    using a symbolic ditribution obtained from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`).
    q : float
        Tsallis `q` parameter; it can be an array of values (default: 1).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    probs: boolean
           If `True`, assumes the data input to be an ordinal probability distribution
           (default: `False`)

    Returns
    -------
     : array
       Values of Tsallis's generalized normalized permutation entropy and 
       Tsallis's generalized statistical complexity for each parameter `q`.

    Examples
    --------
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>>
    >>> p = symbolic_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
    >>> tsallis_complexity_entropy(p, dx=2, probs=True)
    array([0.91829583, 0.06112817])
    >>>
    >>> tsallis_complexity_entropy([1/3, 1/15, 4/15, 2/15, 1/5, 0], dx=3, probs=True)
    array([0.83144548, 0.16576717])
    >>>
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], dx=2, q=[1,2])
    array([[0.91829583, 0.06112817],
           [0.88888889, 0.07619048]])
    >>>
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], q=2, dx=2)
    array([0.88888889, 0.07619048])
    >>>
    >>> tsallis_complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], q=3, dx=3, dy=2)
    array([0.93750181, 0.92972165])
    """
#(Arthur): trouxe essa função pra dentro da função principal tsallis_complexity_entropy (25/08/2020)
    def jensen_tsallis_divergence_max(n_states, q):
        """
        Estimates the maximum value of the Jensen Tsallis divergence\\ [#ribeiro2017]_.
    
        Parameters
        ----------
    
        n_states : int
                  Number of ordinal states.
        q : float
            Renyi `alpha` parameter.
        
        Returns
        -------
         : float
           The maximum divergence value.
        """
        if q==1:
            return(-0.5*(((n_states+1)/n_states)*np.log(n_states+1) + np.log(n_states) - 2*np.log(2*n_states)))
        else:
            # Equation 11 of Physical Review E 95, 062106 (2017).
            return(
                   ((2**(2-q))*n_states - (1+n_states)**(1-q) - n_states*(1+1/n_states)**(1-q) - n_states + 1)/ 
                   ((1-q)*(2**(2-q))*n_states)
                  )

##################################################################################

    if probs==False:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=True)
        h_q              = tsallis_entropy(probabilities[probabilities>0], q, dx, dy, 
                                           tau_x, tau_y, probs=True)
#(Arthur): added this conditional to check whether data is a probability distribution or not (24/08/2020)
    else:
        if len(data)==np.math.factorial(dx*dy):
            probabilities = np.asarray(data)
            h_q           = tsallis_entropy(probabilities[probabilities>0], q, dx, dy, 
                                            tau_x, tau_y, probs=True)
        else:
            raise Exception("The data provided as input is not adequate. The length of" \
                            " the ordinal distribution array does not match the number of" \
                            " possible permutations for dx = {} and dy = {}." \
                            " Please run symbolic_function(..., missing=True) instead of" \
                            " symbolic_function(..., missing=False) to use part of its" \
                            " return as input in tsallis_complexity_entropy() or change the" \
                            " values of the embedding dimensions.".format(dx, dy))

    p             = probabilities[probabilities!=0]
    n             = np.math.factorial(dx*dy)
    uniform_dist  = np.full(n, 1/n)

#(Arthur): changed to add tuple. (19/08/2020)
    if isinstance(q, (tuple, list, np.ndarray)):
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
    Calculates the normalized Rényi's entropy\\ [#jauregui]_ using a
    symbolic ditribution obtained from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`).
    alpha : float
            Rényi `alpha` parameter; it can be an array of values (default: 1).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    probs: boolean
           If `True`, assumes the data input to be an ordinal probability distribution
           (default: `False`)

    Returns
    -------
     : float, array
       The normalized value of Rényi's entropy for each parameter alpha.
    
    Examples
    --------
    >>> renyi_entropy([4,7,9,10,6,11,3], dx=2)
    0.9182958340544894
    >>>
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=2, dx=2)
    0.84799690655495
    >>>
    >>> renyi_entropy([1/3, 1/15, 4/15, 2/15, 1/5, 0], dx=3, probs=True)
    0.8314454838586238
    >>>
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=[1,2], dx=2)
    array([0.91829583, 0.84799691])
    >>>
    >>> renyi_entropy([4,7,9,10,6,11,3], alpha=2, dx=3)
    0.5701944178769374 
    """
    if not probs:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y)
#(Arthur): added a conversion to array and a filter to eliminate 0's from the distribution
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]

    smax = np.log(np.math.factorial(dx*dy))
#(Arthur): changed to add tuple. (19/08/2020)
    if isinstance(alpha, (tuple, list, np.ndarray)):
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


#(Arthur): added the 'probs' parameter. (24/08/2020)
def renyi_complexity_entropy(data, alpha=1, dx=3, dy=1, tau_x=1, tau_y=1, probs=False):
    """
    Calculates permutation entropy\\ [#bandt_pompe]_ and statistical
    complexity\\ [#lopezruiz]_ using Rényi's entropy (these measures define 
    a generalized complexity-entropy causality plane\\ [#jauregui]_) using 
    a symbolic ditribution obtained from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`).
    alpha : float
            Rényi `alpha` parameter; it can be an array of values (default: 1).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    probs: boolean
           If `True`, assumes the data input to be an ordinal probability distribution
           (default: `False`)

    Returns
    -------
     : array
       Values of Rényi's generalized normalized permutation entropy and 
       Rényi's generalized statistical complexity for each parameter alpha.

    Examples
    --------
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>>
    >>> p = symbolic_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
    >>> renyi_complexity_entropy(p, dx=2, probs=True)
    array([0.91829583, 0.06112817])
    >>>
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], alpha=2, dx=2)
    array([0.84799691, 0.08303895])
    >>>
    >>> renyi_complexity_entropy([1/3, 1/15, 4/15, 2/15, 1/5, 0], dx=3, probs=True)
    array([0.83144548, 0.16576717])
    >>>
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], alpha=[1, 2], dx=2)
    array([[0.91829583, 0.06112817],
           [0.84799691, 0.08303895]])
    >>>
    >>> renyi_complexity_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], alpha=3, dx=3, dy=2)
    array([0.21070701, 0.20975673])   
    """
#(Arthur): trouxe essa função pra dentro da função principal renyi_complexity_entropy (25/08/2020)
    def jensen_renyi_divergence_max(n_states, q):
        """
        Estimates the maximum value of the Jensen Renyi divergence\\ [#jauregui]_.
        
        Parameters
        ----------
    
        n_states : int
                  Number of ordinal states.
        q : float
            Renyi `alpha` parameter.
        
        Returns
        -------
         : float
           The maximum divergence value.
        """
        if q==1:
            return(-0.5*(((n_states+1)/n_states)*np.log(n_states+1) + np.log(n_states) - 2*np.log(2*n_states)))
        else:
            # Equation 5 in Physica A 498 (2018) 74-85.
            return(
                    (np.log(((n_states + 1.)**(1. - q) + n_states - 1.)/(2.**(1. - q)*n_states)) + (1. - q)*
                     np.log((n_states + 1.)/(2.*n_states)))/(2.*(q - 1))
                  )

#############################################################################################        

    if probs==False:
        _, probabilities = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=True)
        h_a              = renyi_entropy(probabilities, alpha, dx, dy, tau_x, tau_y, probs=True)
    else:
        if len(data)==np.math.factorial(dx*dy):
            probabilities = np.asarray(data)
            h_a           = renyi_entropy(probabilities, alpha, dx, dy, tau_x, tau_y, probs=True)
#(Arthur): added this conditional to check whether data is a probability distribution or not (24/08/2020)
        else:
            raise Exception("The data provided as input is not adequate. The length of" \
                            " the ordinal distribution array does not match the number of" \
                            " possible permutations for dx = {} and dy = {}." \
                            " Please run symbolic_function(..., missing=True) instead of" \
                            " symbolic_function(..., missing=False) to use part of its" \
                            " return as input in renyi_complexity_entropy() or change the" \
                            " values of the embedding dimensions.".format(dx, dy))

    n             = np.math.factorial(dx*dy)
    uniform_dist  = np.full(n, 1/n)

#(Arthur): changed to add tuple. (19/08/2020)
    if isinstance(alpha, (tuple, list, np.ndarray)):
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


def ordinal_network(data, dx=3, dy=1, tau_x=1, tau_y=1, normalized=True, overlapping=True, directed=True, connections="all"):
    """
    Generates the elements (nodes, edges and edge weights) necessary to obtain
    an ordinal network from data\\ [#small]_\\ :sup:`,`\\ [#pessa2019]_\\ :sup:`,`\\ [#pessa2020]_.
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.ordinal_network`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`) or the nodes, edges and edge weights 
           returned by :func:`ordpy.ordinal_network`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
             Embedding delay (vertical axis) (default: 1).
    normalized : boolean
                 If `True`, edge weights represent transition probabilities 
                 between permutations; if `False`, edge weights are transition
                 counts.
    overlapping : boolean
                  If `True`, data is partitioned into overlapping sliding windows 
                  (default: `True`); if `False`, it does not. 
    directed : boolean
               If `True`, ordinal network links are directed (default: `True`); if `False`, edges
               correspond to an undirected network.
    connections : str
                  The ordinal network is constructed using `'all'` permutation
                  successions in symbolic sequence or only `'horizontal'` or 
                  `'vertical'` successions. Parameter only valid for image data
                  (default: `'all'`). 

    Returns
    -------
     : tuple
       A tuple of arrays containing the nodes, edges and edge weights corresponding
       to an ordinal network mapped from data.
    
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
    """     
#(Arthur): wrote this "inside" function to allow for undirected ordinal networks to be built for time series and images.
    def undirected_ordinal_network(unique_links, occurrences):
        """
        Removes edges that are duplicated in case edge directionality 
        is not relevant to an ordinal network.
        
        Parameters
        ----------
        unique_links : array
                       Edges of a directed ordinal network.
        occurrences : array
                      Edge weights a directed ordinal network.
        
        Returns
        -------
         : tuple
           A tuple of arrays containing the nodes, edges and edge weights
           corresponding to an UNDIRECTED ordinal network mapped from data.
        """
        rev_links    = unique_links[:,[1,0]]
        strlinks     = np.apply_along_axis(np.array2string, 1, unique_links, separator='')
        strrev_links = np.apply_along_axis(np.array2string, 1, rev_links, separator='')

        checked, und_occurrences = [], np.copy(occurrences)                

        for str_rev_ed, cnt in zip(strrev_links, range(strlinks.__len__())):
            #scaping potential double checking
            if not str_rev_ed in checked:
                checked.append(strlinks[cnt])        
                args = np.argwhere(strlinks==str_rev_ed).flatten()
                if args.__len__()==0:
                    pass
                else:
                    for arg_ in args:
                        #checking for autoloops
                        if rev_links[arg_][0]==rev_links[arg_][1]:
                            pass
                        else:
                            und_occurrences[cnt]  += occurrences[arg_]
                            und_occurrences[arg_] = 0
            else:
                pass
            
        #the order of these operations is important
        und_unique_links = np.delete(unique_links, np.argwhere(und_occurrences==0), axis=0)
        und_occurrences  = np.delete(und_occurrences, np.argwhere(und_occurrences==0))

        return (np.unique(und_unique_links), und_unique_links, und_occurrences)
    
#############################################################################################        

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
        
        else: #non overlapping
            partitions = np.concatenate(
                [
                    [np.concatenate(data[j:j+dy*tau_y:tau_y, i:i+dx*tau_x:tau_x]) for i in range(0, nx-(dx-1)*tau_x, dx+(dx-1)*(tau_x-1))] 
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
#
        if directed == True:
            return (np.unique(unique_links), unique_links, occurrences)

        else:                
            return undirected_ordinal_network(unique_links, occurrences)

    #IMAGE DATA
    else:
        if overlapping == True:
            partitions = np.concatenate(
                [
                    [[np.concatenate(data[j:j+dy*tau_y:tau_y,i:i+dx*tau_x:tau_x]) for i in range(nx-(dx-1)*tau_x)]]
                    for j in range(ny-(dy-1)*tau_y)
                ]
            )
            
        else: #non overlapping
            partitions = np.concatenate(
                [
                    [[np.concatenate(data[j:j+dy*tau_y:tau_y, i:i+dx*tau_x:tau_x]) for i in range(0, nx-(dx-1)*tau_x, dx+(dx-1)*(tau_x-1))]] 
                    for j in range(0, ny-(dy-1)*tau_y, dy+(dy-1)*(tau_y-1))
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
#
        if directed == True:
            return (np.unique(unique_links), unique_links, occurrences)

        else:                
            return undirected_ordinal_network(unique_links, occurrences)


def global_node_entropy(data, dx=3, dy=1, tau_x=1, tau_y=1, overlapping=True, connections="all"):
    """
    Calculates global node entropy\\ [#McCullough]_\\ :sup:`,`\\ [#pessa2019]_ for an ordinal
    network obtained from data. (Assumes directed )

    Parameters
    ----------
    data : array or return of :func:`ordpy.ordinal_network`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`) or ordinal network returned by 
           :func:`ordpy.ordinal_network`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
             Embedding delay (vertical axis) (default: 1).
    overlapping : boolean
                  If `True`, data is partitioned into overlapping sliding windows
                  (default: True). If `False`, it does not. 
    connections : str
                  The ordinal network is constructed using `'all'` permutation
                  successions in symbolic sequence or only `'horizontal'` or 
                  `'vertical'` successions. Parameter only valid for image data
                  (default: `'all'`). 
    Returns
    -------
     : float
       The value of global node entropy.

    Examples
    --------
    >>> global_node_entropy([1,2,3,4,5,6,7,8,9], dx=2)
    0.0
    >>>
    >>> global_node_entropy(ordinal_network([1,2,3,4,5,6,7,8,9], dx=2))
    0.0
    >>>
    >>> global_node_entropy(np.random.uniform(size=100000), dx=3)
    1.4988332319747597
    >>>
    >>> global_node_entropy(random_ordinal_network(dx=3))
    1.5
    """
    if len(data)==3 and type(data[0][0])==np.str_:
        nodes, links, weights = data
    else:
#(Arthur): Assumes 'normalized==True' and 'directed==True'.
        nodes, links, weights = ordinal_network(data, dx, dy, tau_x, tau_y, True, 
                                                overlapping, True, connections)

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


def random_ordinal_network(dx=3, dy=1, overlapping=True):
    """
    Generates the elements (nodes, edges and edge weights) necessary to obtain
    a random ordinal network, that is, a  network representing the mapping of a 
    random time series\\ [#pessa2019]_ or a random bidimensional field\\ 
    [#pessa2020]_. Assumes directed links and unitary embbeding 
    delays.
    
    Parameters
    ----------
    dx : int
         Horizontal embedding dimension (default: 3).
    dy : int
         Vertical embedding dimension (default: 1 [time series data])
    overlapping : boolean
                  If `True`, data is partitioned into overlapping sliding windows
                  (default: True). If `False`, it does not. 
    
    Returns
    -------
     : tuple
       A tuple of arrays containing the nodes, edges and edge weights of a
       random ordinal network.

    Examples
    --------
    >>> random_ordinal_network(dx=2)
    (array(['0|1', '1|0'], dtype='<U3'),
     array([['0|1', '0|1'],
            ['0|1', '1|0'],
            ['1|0', '0|1'],
            ['1|0', '1|0']], dtype='<U3'),
     array([0.16666667, 0.33333333, 0.33333333, 0.16666667]))
    >>>
    >>> random_ordinal_network(dx=2, overlapping=False)
    (array(['0|1', '1|0'], dtype='<U3'),
     array([['0|1', '0|1'],
            ['0|1', '1|0'],
            ['1|0', '0|1'],
            ['1|0', '1|0']], dtype='<U3'),
     array([0.25, 0.25, 0.25, 0.25]))
    """

#THEORETICAL RESULTS FOR IMAGE DATA
    if overlapping==True:
        if dx>1 and dy>1:
            allowed_links     = []
            theoretical_probs = []

        #this loop has the job of finding all allowed transitions starting from an ordinal pattern and calculating 
        #the expected probability of such transitions for a random bidimensional scalar field.
            for pattern in itertools.permutations(np.arange(dx*dy).astype('int')):  
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
                links = np.apply_along_axis(np.array2string, 2, links, separator='|')
                for char in ['[', ']']:
                    links = np.char.replace(links, char, '')    

                allowed_links.append(links.tolist())
                theoretical_probs.append(frequencies.tolist())

            edge_array, weight_array = np.concatenate(allowed_links), np.concatenate(theoretical_probs)
            vertices_names           = np.unique(edge_array)

            return vertices_names, edge_array, weight_array
        
    #THEORETICAL RESULTS FOR TIME SERIES DATA
        else:
            vertices_names = []
            edge_array     = []   
            weight_array   = []

            for j in itertools.permutations(np.arange(dx)):
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

            return vertices_names, np.asarray(edge_array), np.asarray(weight_array)
#(Arthur): added more recently the lines below to return the random ordinal network for the case of nonoverlapping partitions.
    else: #nonoverlapping partitions
        permutations = list(itertools.permutations(np.arange(dx*dy).astype('int')))
        edge_array   = [] 

        for j in permutations:
            for i in permutations:
                edge_array.append([j, i])

        edge_array = np.apply_along_axis(np.array2string, 2, edge_array, separator='|')
        for char in ['[', ']']:
            edge_array = np.char.replace(edge_array, char, '')    
            
        weight_array = np.full(shape=len(edge_array), fill_value=1/len(edge_array))

        return np.unique(edge_array), edge_array, weight_array


def missing_states(data, dx=3, dy=1, tau_x=1, tau_y=1, return_fraction=True, return_states=True):
    """
    Searches for ordinal patterns (permutations) which do not occur 
    in data\\ [#amigó]_.
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.symbolic_distribution`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`) or the symbols and probabilities 
           returned by :func:`ordpy.symbolic_distribution`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series (default: 1).
    tau_x : int
            Embedding delay (horizontal axis) (default: 1).
    tau_y : int
            Embedding delay (vertical axis) (default: 1).
    return_fraction : boolean
                      if `True`, returns the fraction of missing ordinal patterns relative 
                      to the total number of ordinal patterns given choices of dx and dy
                      (default: `True`); if `False`, returns the raw number of missing 
                      states.
    return_states : boolean 
                    if `True`, returns the missing ordinal patterns not found in data
                    (default: `True`); if `False`, it only returns the fraction/number 
                    of these missing states.
    Returns
    -------
     : tuple
       A tuple containing an array indicating missing ordinal patterns
       and the relative faction (or raw number) of these missing patterns 
       in data.
    
    Examples
    --------
    >>> missing_states([4,7,9,10,6,11,3], dx=2, return_fraction=False)
    (array([], shape=(0, 2), dtype=int64), 0)
    >>>
    >>> missing_states(symbolic_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True), dx=2)
    (array([], shape=(0, 2), dtype=int64), 0.0)
    >>>
    >>> missing_states([4,7,9,10,6,11,3,5,6,2,3,1], dx=3, return_fraction=False)
    (array([[0, 2, 1],
            [2, 1, 0]]),
     2)
    >>>
    >>> missing_states([4,7,9,10,6,11,3,5,6,2,3,1], dx=3, return_fraction=True, return_states=False)    
    0.3333333333333333
    """
    if not return_states==False:
        if len(data)==2 and type(data[0])==np.ndarray:
            states, probs = data
# inferring dx and dy from data. (It is not possible in the case of missing_links() for the sets of possible links are different 
# for time series and images.
            dx = len(states[0])
            dy = 1

            if len(states)==np.math.factorial(dx*dy):
                pass        
            else:
                raise Exception("The data provided as input is not adequate. The length of" \
                                " the symbolic array does not match the number of" \
                                " possible ordinal patterns for dx = {} and dy = {}." \
                                " Please run symbolic_function(..., missing=True) instead of" \
                                " symbolic_function(..., missing=False) to use its" \
                                " return as input in missing_states() or change the" \
                                " values of the embedding dimensions.".format(dx, dy))

        else:
            states, probs = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=True)
            
        missing_args   = np.argwhere(probs==0).flatten()
        missing_states = states[missing_args]
        
        if return_fraction==True:
            return missing_states, len(missing_args)/np.math.factorial(dx*dy)
        else:
            return missing_states, len(missing_args) 

#not return missing ordinal patterns.
    else:
        if len(data)==2 and type(data[0])==np.ndarray:
            states, _ = data
        else:
            states, _ = symbolic_distribution(data, dx, dy, tau_x, tau_y, return_missing=False)

        if return_fraction==True:
            n = np.math.factorial(dx*dy)
            return (n-len(states))/n
        else:
            return np.math.factorial(dx*dy)-len(states)


def missing_links(data, dx=3, dy=1, return_fraction=True, return_links=True):
    """
    Searches for successions between ordinal patterns (permutations) 
    which do not occur in data. (These successions correspond 
    to directed links in ordinal networks\\ [#pessa2019]_.) Assumes 
    overlapping windows, tau_x = tau_y = 1. In case dx>1 and dy>1, 
    'horizontal' and 'vertical' connections are considered.
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.ordinal_network`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           (:math:`n \\times m`) or the nodes, edges and edge weights
           returned by :func:`ordpy.ordinal_network`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series (default: 1).
    return_fraction : boolean
                      if `True`, returns the fraction of missing links between ordinal 
                      patterns relative to the total number of possible links 
                      (successions) given choices of dx and dy (default: True). If 
                      `False`, returns the raw number of missing links.
    return_links : boolean 
                   if True, returns the missing links between ordinal patterns not found 
                   in data; if False, it only returns the fraction/number of these
                    missing links.

    Returns
    -------
     : tuple
       A tuple containing an array indicating missing links (successions) 
       between ordinal patterns and the relative faction (or raw number) 
       of these missing links in data.
    
    Examples
    --------
    >>> missing_links([4,7,9,10,6,11,3], dx=2, return_fraction=False)
    (array([['1|0', '1|0']], dtype='<U3'), 1)
    >>>
    >>> missing_links(ordinal_network([4,7,9,10,6,11,3], dx=2), dx=2, return_fraction=True)
    (array([['1|0', '1|0']], dtype='<U3'), 0.25)
    >>>
    >>> missing_links([4,7,9,10,6,11,3,5,6,2,3,1], dx=3, return_fraction=False)
    (array([['0|1|2', '0|2|1'],
        ['0|2|1', '1|0|2'],
        ['0|2|1', '1|2|0'],
        ['0|2|1', '2|1|0'],
        ['1|0|2', '0|1|2'],
        ['1|0|2', '0|2|1'],
        ['1|2|0', '0|2|1'],
        ['2|0|1', '2|1|0'],
        ['2|1|0', '1|0|2'],
        ['2|1|0', '1|2|0'],
        ['2|1|0', '2|1|0']], dtype='<U5'),
     11)

    References
    ----------
    For notation or context, please see: Citar o paper aqui.
    """
    if not return_links==False:
        if len(data)==3 and type(data[0][0])==np.str_:
            _, data_links, _ = data
        else:
            _, data_links, _ = ordinal_network(data, dx, dy)
            
        _, all_links, _ = random_ordinal_network(dx=dx, dy=dy)
        missing_links   = setdiff(all_links, data_links)

        if return_fraction==True:
            return missing_links, len(missing_links)/len(all_links)
        else:
            return missing_links, len(missing_links)

#not return missing links
    else:
        if len(data)==3 and type(data[0][0])==np.str_:
            _, data_links, _ = data
        else:
            _, data_links, _ = ordinal_network(data, dx, dy)

#assumes dx==1 means time series data, although the algorithm works for images too
#in this case and te results are different.
        if dy==1:
            all_links = np.math.factorial(dx)*dx

            if return_fraction==True:
                return (all_links-len(data_links))/all_links
            else:
                return all_links-len(data_links)

#not return missing links for image data.
        else:
            _, all_links, _ = random_ordinal_network(dx=dx, dy=dy)
            missing_links   = setdiff(all_links, data_links)

            if return_fraction==True:
                return len(missing_links)/len(all_links)
            else:
                return len(missing_links)


def minimum_complexity_entropy(dx=3, dy=1, size=100):
    """
    Generates data corresponding to the lower bound limit curve of
    the complexity-entropy causality plane\\ [#rosso_curvas]_ [*]_.
    
    Parameters
    ----------
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series (default: 1).
    size : int
           number of returned pairs of values of permutation entropy and 
           statistical complexity.    
           
    Returns
    -------
     : array
       Pairs of values of permutation entropy and statistical complexity 
       which delimit the lower boundary in the complexity-entropy 
       causality plane.

    Notes
    -----
    .. [*] This function was adapted from Sebastian Sippel et al., 
    statcomp: Statistical Complexity and Information Measures for 
    Time Series Analysis, version 0.1.0. (Computer software in R).
    
    Examples
    --------
    >>> minimum_complexity_entropy(dx=3, size=10)
    array([[-0.        , -0.        ],
           [ 0.25534537,  0.17487933],
           [ 0.43376898,  0.21798915],
           [ 0.57926767,  0.21217508],
           [ 0.70056313,  0.18118876],
           [ 0.80117413,  0.1378572 ],
           [ 0.88242522,  0.09085847],
           [ 0.94417789,  0.04723328],
           [ 0.98476216,  0.01396083],
           [ 1.        ,  0.        ]])
    >>>
    >>> minimum_complexity_entropy(dx=4, size=5)
    array([[-0.00000000e+00, -0.00000000e+00],
           [ 4.09625322e-01,  2.10690777e-01],
           [ 6.90580872e-01,  1.84011590e-01],
           [ 8.96072453e-01,  8.64190054e-02],
           [ 1.00000000e+00, -3.66606083e-16]])
    """
    size        += 1
    N            = np.math.factorial(dx*dy)
    prob_params  = np.linspace(1/N, 1, num=size-1)
    uniform_dist = np.full(N, 1/N)
    
    hc_ = []
    for i in range(size-1):
        probabilities    = np.full(shape=N, fill_value=(1-prob_params[i])/(N-1))
        probabilities[0] = prob_params[i]
#                
        h = permutation_entropy(probabilities, dx, dy, probs=True)
#        
        p_plus_u_over_2      = (uniform_dist + probabilities)/2  
        s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

        probabilities = probabilities[probabilities!=0]
        s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
        s_of_u_over_2 = np.log(N)/2.

        js_div_max = -0.5*(((N+1)/N)*np.log(N+1) + np.log(N) - 2*np.log(2*N))    
        js_div     = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

        hc_.append([h, h*js_div/js_div_max])
        
    return np.flip(hc_, axis=0)


def maximum_complexity_entropy(dx=3, dy=1, m=1):
    """
    Generates data corresponding to the upper bound limit curve of
    the complexity-entropy causality plane\\ [#rosso_curvas]_ [*]_.
    
    Parameters
    ----------
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series (default: 1).
    m : int
        The number of pairs of values of permutation entropy and statistical 
        complexity is equal to :math: `[(d_x \\times d_y)!-1] \\times m`.
           
    Returns
    -------
     : array
       Pairs of values of permutation entropy and statistical complexity 
       which delimit the upper boundary in the complexity-entropy 
       causality plane.

    Notes
    -----
    .. [*] This function was adapted from Sebastian Sippel et al., 
    statcomp: Statistical Complexity and Information Measures for 
    Time Series Analysis, version 0.1.0. (Computer software in R).
    
    Examples
    --------
    >>> maximum_complexity_entropy(dx=3, dy=1, m=1)
    array([[-0.        , -0.        ],
           [ 0.38685281,  0.27123863],
           [ 0.61314719,  0.29145164],
           [ 0.77370561,  0.22551573],
           [ 0.8982444 ,  0.12181148]])
    >>>
    >>> maximum_complexity_entropy(dx=3, dy=1, m=2)
    array([[-0.        , -0.        ],
           [ 0.251463  ,  0.19519391],
           [ 0.38685281,  0.27123863],
           [ 0.57384034,  0.28903864],
           [ 0.61314719,  0.29145164],
           [ 0.76241899,  0.2294088 ],
           [ 0.77370561,  0.22551573],
           [ 0.89621768,  0.12320718],
           [ 0.8982444 ,  0.12181148],
           [ 1.        ,  0.        ]])
    """
    
    N              = np.math.factorial(dx*dy)
    hlist_, clist_ = np.zeros(shape=(N-1,m)), np.zeros(shape=(N-1,m))

    for i in range(N-1):
    #     print('i: ', i)
        p             = np.zeros(shape=N)
        uniform_dist  = np.full(N, 1/N)
        prob_params   = np.linspace(0, 1/N, num=m)
    #
        for k in range(len(prob_params)):
    #         print('\m: ', k)
            p[0] = prob_params[k]
            for j in range(1,N-i):
    #             print('j: ', j)
                p[j] = (1-prob_params[k])/(N-i-1)

    #         print(p)
            h = permutation_entropy(p, dx, dy, probs=True)

            p_plus_u_over_2      = (uniform_dist + p)/2  
            s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2))

            p_non_zero    = p[p!=0]
            s_of_p_over_2 = -np.sum(p_non_zero*np.log(p_non_zero))/2
            s_of_u_over_2 = np.log(N)/2.

            js_div_max = -0.5*(((N+1)/N)*np.log(N+1) + np.log(N) - 2*np.log(2*N))    
            js_div     = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

            hlist_[i, k] = h
            clist_[i, k] = h*js_div/js_div_max
            
#flatenning the arrays and ordering the pairs of values.
    hlist_ = hlist_.flatten()
    clist_ = clist_.flatten()
    args   = np.argsort(hlist_)
    
    return np.asarray((hlist_[args], clist_[args])).T

#     ranges_      = np.linspace(0,1,m+1)
#     final_hlist_ = []
#     final_clist_ = []

#     for i in range(m+1)[:-1]:
#         try:
#             x, y = np.argwhere((hlist_>ranges_[i]) & (hlist_<ranges_[i+1]))[0]
#     #         print('location: ', x, y)
#             final_hlist_.append(hlist_[x, y])
#     #         print(hlist_[x, y])
#             final_clist_.append(clist_[x, y])
#         except IndexError:
#             pass
        
#     return np.asarray((final_hlist_, final_clist_)).T