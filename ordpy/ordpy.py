"""
ordpy: A Python Package for Data Analysis with Permutation Entropy and Ordinal Network Methods
===============================================================================================

``ordpy`` is a pure Python module [#pessa2021]_ that implements data analysis methods based
on Bandt and Pompe's [#bandt_pompe]_ symbolic encoding scheme.

.. note::

   If you have used ``ordpy`` in a scientific publication, we would appreciate 
   citations to the following reference [#pessa2021]_:

   - A. A. B. Pessa, H. V. Ribeiro, `ordpy: A Python package for data 
     analysis with permutation entropy and ordinal network methods 
     <https://doi.org/10.1063/5.0049901>`_, Chaos 31, 063110 (2021). `arXiv:2102.06786 <https://arxiv.org/abs/2102.06786>`_

    .. code-block:: bibtex
    
       @misc{pessa2021ordpy,
         title         = {ordpy: A Python package for data analysis with permutation entropy and ordinal network methods}, 
         author        = {Arthur A. B. Pessa and Haroldo V. Ribeiro},
         journal       = {Chaos: An Interdisciplinary Journal of Nonlinear Science},
         volume        = {31},
         number        = {6},
         pages         = {063110},
         year          = {2021},
         doi           = {10.1063/5.0049901},
        }

``ordpy`` implements the following data analysis methods:

Released on version 1.0.0 (February 2021):

- Permutation entropy for time series [#bandt_pompe]_ and images [#ribeiro_2012]_;
- Complexity-entropy plane for time series [#lopezruiz]_, [#rosso]_ and 
  images [#ribeiro_2012]_;
- Multiscale complexity-entropy plane for time series [#zunino2012]_ and 
  images [#zunino2016]_;
- Tsallis [#ribeiro2017]_ and Rényi [#jauregui]_ generalized complexity-entropy
  curves for time series and images;
- Ordinal networks for time series [#small]_, [#pessa2019]_ and 
  images [#pessa2020]_;
- Global node entropy of ordinal networks for 
  time series [#McCullough]_, [#pessa2019]_ and images [#pessa2020]_.
- Missing ordinal patterns [#amigo]_ and missing transitions between ordinal
  patterns [#pessa2019]_ for time series and images.

Released on version 1.1.0 (January 2023):

- Weighted permutation entropy for time series [#fadlallah]_ and images;
- Fisher-Shannon plane for time series [#olivares]_ and images;
- Permutation Jensen-Shannon distance for time series [#zunino2022]_ and images;
- Four pattern permutation contrasts (up-down balance, persistence, 
  rotational-asymmetry, and up-down scaling.) for time series [#bandt]_; 
- Smoothness-structure plane for images [#bandt_wittfeld]_.

References
----------

.. [#pessa2021] Pessa, A. A., & Ribeiro, H. V. (2021). ordpy: A Python package
   for data analysis with permutation entropy and ordinal networks methods. 
   Chaos, 31, 063110.

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
   Ordinal Networks. Physical Review E, 102, 052312.

.. [#McCullough] McCullough, M., Small, M., Iu, H. H. C., & Stemler, T. (2017).
   Multiscale Ordinal Network Analysis of Human Cardiac Dynamics.
   Philosophical Transactions of the Royal Society A, 375, 20160292.

.. [#amigo] Amigó, J. M., Zambrano, S., & Sanjuán, M. A. F. (2007).
   True and False Forbidden Patterns in Deterministic and Random Dynamics.
   Europhysics Letters, 79, 50001.

.. [#fadlallah] Fadlallah B., Chen, B., Keil A. & Príncipe, J. (2013). 
   Weighted-permutation entropy: a complexity measure for time series 
   incorporating amplitude information. Physical Review E, 97, 022911.

.. [#olivares] Olivares, F., Plastino, A., & Rosso, O. A. (2012). 
   Contrasting chaos with noise via local versus global 
   information quantifiers. Physics Letters A, 376, 1577–1583.

.. [#zunino2022] Zunino L., Olivares, F., Ribeiro H. V. & Rosso, O. A. (2022). 
   Permutation Jensen-Shannon distance: A versatile and fast symbolic tool
   for complex time-series analysis. Physical Review E, 105, 045310.

.. [#bandt] Bandt, C. (2022). Statistics and modelling of order patterns in
   univariate series. arXiv:2212.14386.

.. [#bandt_wittfeld] Bandt, C., & Wittfeld, K. (2022). Two new parameters for 
   the ordinal analysis of images. arXiv:2212.14643.

.. [#rosso_curvas] Martin, M. T., Plastino, A., & Rosso, O. A. (2006). 
   Generalized Statistical Complexity Measures: Geometrical and 
   Analytical Properties, Physica A, 369, 439–462.

Installing
==========

Ordpy can be installed via the command line using

.. code-block:: console

   pip install ordpy

or you can directly clone its git repository:

.. code-block:: console

   git clone https://github.com/arthurpessa/ordpy.git
   cd ordpy
   pip install -e .


Basic usage
===========

We provide a `notebook <https://github.com/arthurpessa/ordpy/blob/master/examples/ordpy.ipynb>`_
illustrating how to use ``ordpy``. This notebook reproduces all figures of our
article [#pessa2021]_. The code below shows simple applications of ``ordpy``.

.. code-block:: python

    #Complexity-entropy plane for logistic map and Gaussian noise.

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


    f, ax = plt.subplots(figsize=(8.19, 6.3))

    for HC_, label_ in zip(HC, ['Period-2 (a=3.05)', 
                                'Period-8 (a=3.55)', 
                                'Chaotic (a=4)', 
                                'Gaussian noise']):
        ax.scatter(*HC_, label=label_, s=100)
        
    ax.set_xlabel('Permutation entropy, $H$')
    ax.set_ylabel('Statistical complexity, $C$')

    ax.legend()

.. figure:: https://raw.githubusercontent.com/arthurpessa/ordpy/master/examples/figs/sample_fig.png
   :height: 489px
   :width: 633px
   :scale: 80 %
   :align: center

.. code-block:: python

    #Ordinal networks for logistic map and Gaussian noise.

    import numpy as np
    import igraph
    import ordpy
    from matplotlib import pylab as plt
    from IPython.core.display import display, SVG

    def logistic(a=4, n=100000, x0=0.4):
        x = np.zeros(n)
        x[0] = x0
        for i in range(n-1):
            x[i+1] = a*x[i]*(1-x[i])
        return(x)

    time_series = [logistic(a=4), np.random.normal(size=100000)]

    vertex_list, edge_list, edge_weight_list = list(), list(), list()
    for series in time_series:
        v_, e_, w_   = ordpy.ordinal_network(series, dx=4)
        vertex_list += [v_]
        edge_list   += [e_]
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

    for graph_, label_ in zip(graphs, ['Chaotic (a=4)', 
                                       'Gaussian noise']):
        print(label_)
        display(SVG(igplot(graph_)._repr_svg_()))

.. figure:: https://raw.githubusercontent.com/arthurpessa/ordpy/master/examples/figs/sample_net.png
   :height: 1648px
   :width: 795px
   :scale: 50 %
   :align: center


List of functions
=================
"""
import itertools
import math
import warnings

import numpy as np


def ordinal_sequence(data, dx=3, dy=1, taux=1, tauy=1, overlapping=True, tie_precision=None):
    """
    Applies the Bandt and Pompe\\ [#bandt_pompe]_ symbolization approach to obtain 
    a sequence of ordinal patterns (permutations) from data.
    
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
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    overlapping : boolean
                  If `True`, **data** is partitioned into overlapping sliding 
                  windows (default: `True`). If `False`, adjacent partitions are
                  non-overlapping.
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    number of decimals (default: `None`).

    Returns
    -------
     : array
       Array containing the sequence of ordinal patterns.

    Examples
    --------
    >>> ordinal_sequence([4,7,9,10,6,11,3], dx=2)
    array([[0, 1],
           [0, 1],
           [0, 1],
           [1, 0],
           [0, 1],
           [1, 0]])
    >>> 
    >>> ordinal_sequence([4,7,9,10,6,11,3], dx=2, taux=2)
    array([[0, 1],
           [0, 1],
           [1, 0],
           [0, 1],
           [1, 0]])
    >>>
    >>> ordinal_sequence([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2)
    array([[[0, 1, 3, 2],
            [1, 0, 2, 3],
            [0, 1, 2, 3]], 
           [[1, 2, 3, 0],
            [0, 1, 3, 2],
            [0, 1, 2, 3]]])
    >>>
    >>> ordinal_sequence([1.3, 1.2], dx=2), ordinal_sequence([1.3, 1.2], dx=2, tie_precision=0)
    (array([[1, 0]]), array([[0, 1]]))
    """
    try:
        ny, nx = np.shape(data)
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])

    if tie_precision is not None: 
        data = np.round(data, tie_precision)
    
    if ny==1: #time series
        if overlapping == True:
            partitions = np.apply_along_axis(func1d       = np.lib.stride_tricks.sliding_window_view, 
                                             axis         = 1, 
                                             arr          = data, 
                                             window_shape = (dx+(dx-1)*(taux-1),)
                                            )
            partitions = partitions[::, ::, ::taux].reshape(-1, dx) 
                                                                    
        else: #non overlapping
            partitions = np.concatenate(
                [
                    [np.concatenate(data[j:j+dy*tauy:tauy, i:i+dx*taux:taux]) for i in range(0, nx-(dx-1)*taux, dx+(dx-1)*(taux-1))] 
                    for j in range(ny-(dy-1)*tauy)
                ]
            )
        symbols = np.apply_along_axis(np.argsort, 1, partitions)
        
    else: #image
        if overlapping == True:
            partitions        = np.lib.stride_tricks.sliding_window_view(x            = data, 
                                                                         window_shape = (dy+(dy-1)*(tauy-1), dx+(dx-1)*(taux-1))
                                                                        )
            rows, columns, *_ = np.shape(partitions)
            partitions        = partitions[::,::,::tauy,::taux].reshape(rows, columns, dx*dy)

        else: #non overlapping
            partitions = np.concatenate(
                [
                    [[np.concatenate(data[j:j+dy*tauy:tauy, i:i+dx*taux:taux]) for i in range(0, nx-(dx-1)*taux, dx+(dx-1)*(taux-1))]] 
                    for j in range(0, ny-(dy-1)*tauy, dy+(dy-1)*(tauy-1))
                ]
            )            
        symbols = np.apply_along_axis(np.argsort, 2, partitions)
    
    return symbols


def ordinal_distribution(data, dx=3, dy=1, taux=1, tauy=1, return_missing=False, tie_precision=None):
    """
    Applies the Bandt and Pompe\\ [#bandt_pompe]_ symbolization approach to obtain 
    a probability distribution of ordinal patterns (permutations) from data.
    
    Parameters
    ----------
    data : array 
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    return_missing: boolean
                    If `True`, it returns ordinal patterns not appearing in the 
                    symbolic sequence obtained from **data**. If `False`, these
                    missing patterns (permutations) are omitted (default: `False`).
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).
    Returns
    -------
     : tuple
       Tuple containing two arrays, one with the ordinal patterns occurring in data 
       and another with their corresponding probabilities.

    Examples
    --------
    >>> ordinal_distribution([4,7,9,10,6,11,3], dx=2)
    (array([[0, 1],
            [1, 0]]),
     array([0.66666667, 0.33333333]))
    >>>
    >>> ordinal_distribution([4,7,9,10,6,11,3], dx=3, return_missing=True)
    (array([[0, 1, 2],
            [1, 0, 2],
            [2, 0, 1],
            [0, 2, 1],
            [1, 2, 0],
            [2, 1, 0]]),
     array([0.4, 0.2, 0.4, 0. , 0. , 0. ]))
    >>>
    >>> ordinal_distribution([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    (array([[0, 1, 3, 2],
            [1, 0, 2, 3],
            [1, 2, 3, 0]]),
     array([0.5 , 0.25, 0.25]))
    >>>
    >>> ordinal_distribution([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, taux=2)
    (array([[0, 1, 3, 2],
            [0, 2, 1, 3],
            [1, 3, 2, 0]]),
     array([0.5 , 0.25, 0.25]))
    """
    def setdiff(a, b):
        """
        Searches for elements (subarrays) in `a` that are not contained in `b` [*]_. 

        Parameters
        ----------    
        a : tuples, lists or arrays
            Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
            \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]`.
        b : tuples, lists or arrays
            Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
            \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]`.
        
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
        a = np.asarray(a).astype('int64')
        b = np.asarray(b).astype('int64')

        _, ncols = a.shape

        dtype={'names':['f{}'.format(i) for i in range(ncols)], 'formats':ncols * [a.dtype]}

        C = np.setdiff1d(a.view(dtype), b.view(dtype))
        C = C.view(a.dtype).reshape(-1, ncols)

        return(C)


    symbols = ordinal_sequence(data, dx=dx, dy=dy, taux=taux, tauy=tauy, overlapping=True, tie_precision=tie_precision)

    if np.shape(symbols).__len__() == 3:
        symbols = symbols.reshape(-1, dx*dy) #we reshape the symbols array to facilitate the application of np.unique below.
    else:
        pass
    
    symbols, symbols_count = np.unique(symbols, return_counts=True, axis=0)
    probabilities          = symbols_count/symbols_count.sum()

    if return_missing==False:
        return symbols, probabilities
    
    else:
        all_symbols   = np.asarray(list(itertools.permutations(range(dx*dy))), dtype='int')
        miss_symbols  = setdiff(all_symbols, symbols)
        symbols       = np.concatenate((symbols, miss_symbols)).astype('int')
        probabilities = np.concatenate((probabilities, np.zeros(miss_symbols.__len__())))
        
        return symbols, probabilities


def permutation_entropy(data, dx=3, dy=1, taux=1, tauy=1, base=2, normalized=True, probs=False, tie_precision=None):
    """
    Calculates the Shannon entropy using an ordinal distribution obtained from
    data\\ [#bandt_pompe]_\\ :sup:`,`\\ [#ribeiro_2012]_.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
           :func:`ordpy.ordinal_distribution`).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    base : str, int
           Logarithm base in Shannon's entropy. Either 'e' or 2 (default: 2).
    normalized: boolean
                If `True`, permutation entropy is normalized by its maximum value 
                (default: `True`). If `False`, it is not.
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : float
       Value of permutation entropy.
    
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
    >>> permutation_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, taux=2, base='e', normalized=False)
    1.0397207708399179
    """
    if not probs:
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=False, tie_precision=tie_precision)
    else:
        probabilities    = np.asarray(data)
        probabilities    = probabilities[probabilities>0]

    logfunc = np.log if base=='e' else np.log2
    s       = -np.sum(probabilities*logfunc(probabilities))

    if normalized==True:        
        smax = logfunc(float(math.factorial(dx*dy)))
        return s/smax
    else:
        return s 


def complexity_entropy(data, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates permutation entropy\\ [#bandt_pompe]_ and statistical
    complexity\\ [#lopezruiz]_\\ :sup:`,`\\ [#rosso]_ using an ordinal 
    distribution obtained from data.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
           :func:`ordpy.ordinal_distribution`).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : tuple
       Values of normalized permutation entropy and statistical complexity.
    
    Examples
    --------
    >>> complexity_entropy([4,7,9,10,6,11,3], dx=2)
    (0.9182958340544894, 0.06112816548804511)
    >>>
    >>> p = ordinal_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
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

    n = float(math.factorial(dx*dy)) 

    if probs==False:
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=False, tie_precision=tie_precision)   
        h                = permutation_entropy(probabilities, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]
        h             = permutation_entropy(probabilities, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)

    n_states_not_occuring = n-len(probabilities)
    uniform_dist          = 1/n

    p_plus_u_over_2      = (uniform_dist + probabilities)/2  
    s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2)) - (0.5*uniform_dist)*np.log(0.5*uniform_dist)*n_states_not_occuring

    s_of_p_over_2 = -np.sum(probabilities*np.log(probabilities))/2
    s_of_u_over_2 = np.log(n)/2.

    js_div_max = -0.5*(((n+1)/n)*np.log(n+1) + np.log(n) - 2*np.log(2*n))    
    js_div     = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2

    return h, h*js_div/js_div_max


def logq(x, q=1):
    """
    Calculates the `q`-logarithm of x.

    Parameters
    ----------
    x : float, array
        Real number or array containing real numbers.
    q : float
        Tsallis's `q` parameter (default: 1).

    Returns
    -------
     : float or array
       Value or array of values containing the `q`-logarithm of x.

    Notes
    -----
    The `q`-logarithm of `x` is defined as\\ [*]_

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


def tsallis_entropy(data, q=1, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates the normalized Tsallis permutation entropy\\ [#ribeiro2017]_ using 
    an ordinal distribution obtained from data.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
       :func:`ordpy.ordinal_distribution`).
    q : float or array
        Tsallis's `q` parameter (default: 1); an array of values is also 
        accepted for this parameter.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
            Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : float, array
       Value(s) of the normalized Tsallis permutation entropy.
    
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
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, tie_precision=tie_precision)
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]

    if isinstance(q, (tuple, list, np.ndarray)):
        s = []
        for q_ in q:
            smax          = logq(float(math.factorial(dx*dy)), q_)
            lnq_1_over_p  = logq(1./probabilities, q_)
            s            += [np.sum(probabilities*lnq_1_over_p)/smax]
        s = np.asarray(s)
    else:
        smax          = logq(float(math.factorial(dx*dy)), q)
        lnq_1_over_p  = logq(1./probabilities, q)
        s             = np.sum(probabilities*lnq_1_over_p)/smax

    return s


def tsallis_complexity_entropy(data, q=1, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates the Tsallis normalized permutation entropy and statistical
    complexity\\ [#ribeiro2017]_ using an ordinal distribution obtained from 
    data.

    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
       :func:`ordpy.ordinal_distribution`).
    q : float, array
        Tsallis's `q` parameter (default: 1); an array of values is also 
        accepted for this parameter.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : array
       Value(s) of normalized permutation entropy and statistical complexity in 
       Tsallis's formalism.

    Examples
    --------
    >>> tsallis_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>>
    >>> p = ordinal_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
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


    n = float(math.factorial(dx*dy))

    if probs==False:
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=False, tie_precision=tie_precision)
        h_q              = tsallis_entropy(probabilities, q, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)
    else:
        probabilities = np.asarray(data)      
        probabilities = probabilities[probabilities!=0]
        h_q           = tsallis_entropy(probabilities, q, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)

    n_states_not_occur = n-len(probabilities)
    uniform_dist       = 1/n
    p                  = probabilities

    if isinstance(q, (tuple, list, np.ndarray)):
        jt_div     = []
        jt_div_max = []
        for q_ in q:
            #first and second terms of Jensen-Tsallis divergence: equation 10 of Physical Review E 95, 062106 (2017).
            first_term  = (uniform_dist + p)/(2*p)
            first_term  = -0.5*np.sum(p*logq(first_term, q_)) 

            second_term = n*(uniform_dist + p)/2
            second_term = -(0.5/n) * ( np.sum(logq(second_term, q_)) + logq(0.5, q_)*n_states_not_occur )

            jt_div      += [first_term + second_term]
            jt_div_max  += [jensen_tsallis_divergence_max(n,q_)]

        jt_div     = np.asarray(jt_div)
        jt_div_max = np.asarray(jt_div_max)

    else:
        first_term  = (uniform_dist + p)/(2*p)
        first_term  = -0.5*np.sum(p*logq(first_term, q)) 

        second_term = n*(uniform_dist + p)/2
        second_term = -(0.5/n) * ( np.sum(logq(second_term, q)) + logq(0.5, q)*n_states_not_occur )

        jt_div      = first_term + second_term
        jt_div_max  = jensen_tsallis_divergence_max(n,q)

    return np.asarray([h_q, h_q*jt_div/jt_div_max]).T
    

def renyi_entropy(data, alpha=1, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates the normalized Rényi permutation entropy\\ [#jauregui]_ using an 
    ordinal distribution obtained from data.

    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
          :func:`ordpy.ordinal_distribution`).
    alpha : float, array
            Rényi's `alpha` parameter (default: 1); an array of values is also 
            accepted for this parameter.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : float, array
       Value(s) of the normalized Rényi permutation entropy.

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
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, tie_precision=tie_precision)
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities>0]

    smax = np.log(float(math.factorial(dx*dy)))

    if isinstance(alpha, (tuple, list, np.ndarray)):
        s = []
        for alpha_ in alpha:
            if alpha_!=1:
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


def renyi_complexity_entropy(data, alpha=1, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates the Rényi normalized permutation entropy and statistical
    complexity\\ [#jauregui]_ using an ordinal distribution obtained from 
    data.

    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
           :func:`ordpy.ordinal_distribution`).
    alpha : float, array
            Rényi's `alpha` parameter (default: 1); an array of values is also 
            accepted for this parameter.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : array
       Value(s) of normalized permutation entropy and statistical complexity in 
       Rényi's formalism.

    Examples
    --------
    >>> renyi_complexity_entropy([4,7,9,10,6,11,3], dx=2)
    array([0.91829583, 0.06112817])
    >>>
    >>> p = ordinal_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
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



    n = float(math.factorial(dx*dy))

    if probs==False:
        _, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=False, tie_precision=tie_precision)
        h_a              = renyi_entropy(probabilities, alpha, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)
    else:
        probabilities = np.asarray(data)
        probabilities = probabilities[probabilities!=0]
        h_a           = renyi_entropy(probabilities, alpha, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)

    n_states_not_occur = n - len(probabilities)
    uniform_dist       = 1/n
    p                  = probabilities
    
    if isinstance(alpha, (tuple, list, np.ndarray)):
        jr_div     = []
        jr_div_max = []
        for alpha_ in alpha:
            if alpha_==1: #limit of Shannon entropy
                p_plus_u_over_2      = (uniform_dist + p)/2
                s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2)) - (0.5*uniform_dist)*np.log(0.5*uniform_dist)*n_states_not_occur

                s_of_p_over_2 = -np.sum(p*np.log(p))/2
                s_of_u_over_2 = np.log(n)/2.

                jr_div       += [s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2]
            else:
                # Equation 4 in Physica A 498 (2018) 74-85.
                first_term  = ((p + uniform_dist)/2)**(1-alpha_)
                first_term  = np.log(np.sum(first_term * p**alpha_))
                second_term = np.log(np.sum((1/(n**alpha_)) * ((p + uniform_dist)/2)**(1-alpha_)) + ( n_states_not_occur * (1/n**alpha_) * (1/(2*n))**(1-alpha_) ))
                
                jr_div     += [(1/(2*(alpha_-1)))*(first_term + second_term)]

            jr_div_max += [jensen_renyi_divergence_max(n, alpha_)]

        jr_div     = np.asarray(jr_div)
        jr_div_max = np.asarray(jr_div_max)

    else:
        if alpha==1:
            p_plus_u_over_2      = (uniform_dist + p)/2  
            s_of_p_plus_u_over_2 = -np.sum(p_plus_u_over_2*np.log(p_plus_u_over_2)) - (0.5*uniform_dist)*np.log(0.5*uniform_dist)*n_states_not_occur

            s_of_p_over_2 = -np.sum(p*np.log(p))/2
            s_of_u_over_2 = np.log(n)/2.

            jr_div        = s_of_p_plus_u_over_2 - s_of_p_over_2 - s_of_u_over_2
        else:
            first_term  = ((p + uniform_dist)/2)**(1-alpha)
            first_term  = np.log(np.sum(first_term * p**alpha))    
            second_term = np.log(np.sum((1/(n**alpha)) * ((p + uniform_dist)/2)**(1-alpha)) + ( n_states_not_occur * (1/n**alpha) * (1/(2*n))**(1-alpha) ))
        
            jr_div = (1/(2*(alpha-1)))*(first_term + second_term)
            
        jr_div_max = jensen_renyi_divergence_max(n, alpha)

    return np.asarray([h_a,h_a*jr_div/jr_div_max ]).T


def ordinal_network(data, dx=3, dy=1, taux=1, tauy=1, normalized=True, overlapping=True, directed=True, connections="all", tie_precision=None):
    """
    Maps a data set into the elements (nodes, edges and edge weights) of its
    corresponding ordinal network representation\\ [#small]_\\ :sup:`,`\\ 
    [#pessa2019]_\\ :sup:`,`\\ [#pessa2020]_.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    normalized : boolean
                 If `True`, edge weights represent transition probabilities between 
                 permutations (default: `True`). If `False`, edge weights are 
                 transition counts.
    overlapping : boolean
                  If `True`, **data** is partitioned into overlapping sliding 
                  windows (default: `True`). If `False`, adjacent partitions are
                  non-overlapping.
    directed : boolean
               If `True`, ordinal network edges are directed (default: `True`). If 
               `False`, edges are undirected.
    connections : str
                  The ordinal network is constructed using `'all'` permutation
                  successions in a symbolic sequence or only `'horizontal'` or 
                  `'vertical'` successions. Parameter only valid for image data
                  (default: `'all'`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : tuple
       Tuple containing three arrays corresponding to nodes, edges and edge weights 
       of an ordinal network.
    
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
    >>> ordinal_network([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, normalized=True, connections='horizontal')
    (array(['0|1|3|2', '1|0|2|3', '1|2|3|0'], dtype='<U7'),
     array([['0|1|3|2', '1|0|2|3'],
            ['1|2|3|0', '0|1|3|2']], dtype='<U7'),
     array([0.5, 0.5]))
    """     

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
           Tuple of arrays containing the nodes, edges and edge weights
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



    try:
        ny, nx = np.shape(data)
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])

    if tie_precision is not None:
        data = np.round(data, tie_precision)

    if ny==1: #time series
        if ny==1:
            if overlapping == True:
                partitions = np.apply_along_axis(func1d       = np.lib.stride_tricks.sliding_window_view, 
                                                axis         = 1, 
                                                arr          = data, 
                                                window_shape = (dx+(dx-1)*(taux-1),)
                                                )
                partitions = partitions[::, ::, ::taux].reshape(-1, dx) 
                                                                        
            else: #non overlapping
                partitions = np.concatenate(
                    [
                        [np.concatenate(data[j:j+dy*tauy:tauy, i:i+dx*taux:taux]) for i in range(0, nx-(dx-1)*taux, dx+(dx-1)*(taux-1))] 
                        for j in range(ny-(dy-1)*tauy)
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

        if directed == True:
            return (np.unique(unique_links), unique_links, occurrences)
        else:                
            return undirected_ordinal_network(unique_links, occurrences)

    else: #image
        if overlapping == True:
            partitions        = np.lib.stride_tricks.sliding_window_view(x            = data, 
                                                                         window_shape = (dy+(dy-1)*(tauy-1), dx+(dx-1)*(taux-1))
                                                                        )
            rows, columns, *_ = np.shape(partitions)
            partitions        = partitions[::,::,::tauy,::taux].reshape(rows, columns, dx*dy)

        else: #non overlapping
            partitions = np.concatenate(
                [
                    [[np.concatenate(data[j:j+dy*tauy:tauy, i:i+dx*taux:taux]) for i in range(0, nx-(dx-1)*taux, dx+(dx-1)*(taux-1))]] 
                    for j in range(0, ny-(dy-1)*tauy, dy+(dy-1)*(tauy-1))
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

        if directed == True:
            return (np.unique(unique_links), unique_links, occurrences)
        else:                
            return undirected_ordinal_network(unique_links, occurrences)


def global_node_entropy(data, dx=3, dy=1, taux=1, tauy=1, overlapping=True, connections="all", tie_precision=None):
    """
    Calculates global node entropy\\ [#pessa2019]_\\ :sup:`,`\\ [#McCullough]_ for an
    ordinal network obtained from data. (Assumes directed and weighted edges).

    Parameters
    ----------
    data : array, return of :func:`ordpy.ordinal_network`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]` 
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           or an ordinal network returned by :func:`ordpy.ordinal_network`\\ [*]_.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    overlapping : boolean
                  If `True`, **data** is partitioned into overlapping sliding 
                  windows (default: `True`). If `False`, adjacent partitions are
                  non-overlapping.
    connections : str
                  The ordinal network is constructed using `'all'` permutation
                  successions in a symbolic sequence or only `'horizontal'` or 
                  `'vertical'` successions. Parameter only valid for image data
                  (default: `'all'`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : float
       Value of global node entropy.

    Notes
    -----
    .. [*] In case **data** is an ordinal network returned by 
           :func:`ordpy.ordinal_network`, the parameters of 
           :func:`ordpy.global_node_entropy` are inferred from the network.

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
    >>>
    >>> global_node_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, connections='horizontal')
    0.25
    >>>
    >>> global_node_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, connections='vertical')
    0.0
    """
    if len(data)==3 and type(data[0][0])==np.str_:
        nodes, links, weights = data
    else:
        #assumes 'normalized==True' and 'directed==True'.
        nodes, links, weights = ordinal_network(data, dx, dy, taux, tauy, True, overlapping, True,
                                                connections, tie_precision=tie_precision)

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
    Generates the nodes, edges and edge weights of a random ordinal network: the 
    theoretically expected network representation of a random time series\\ 
    [#pessa2019]_ or a random two-dimensional array\\ [#pessa2020]_. 
    (Assumes directed edges and unitary embedding delays.)
    
    Parameters
    ----------
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    overlapping : boolean
                  If `True`, **data** is partitioned into overlapping sliding 
                  windows (default: `True`). If `False`, adjacent partitions are
                  non-overlapping.
    
    Returns
    -------
     : tuple
       Tuple containing three arrays corresponding to nodes, edges and edge weights 
       of the random ordinal network.

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

    if overlapping==True: #theoretical results for image
        if dx>1 and dy>1:
            allowed_links     = []
            theoretical_probs = []

            #this loop finds all allowed transitions starting from an ordinal pattern and calculates the 
            #expected probability of such transitions for a random bidimensional scalar field.
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
                frequencies               = frequencies/frequencies.sum()/math.factorial(dx*dy)

                #transforming lists of permutations to strings to save the transitions        
                links = np.apply_along_axis(np.array2string, 2, links, separator='|')
                for char in ['[', ']']:
                    links = np.char.replace(links, char, '')    

                allowed_links.append(links.tolist())
                theoretical_probs.append(frequencies.tolist())

            edge_array, weight_array = np.concatenate(allowed_links), np.concatenate(theoretical_probs)
            vertices_names           = np.unique(edge_array)

            return vertices_names, edge_array, weight_array
        
        else: #theoretical results for time series
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
                        weight_array.append(1/((dx+1)*math.factorial(dx)))
                    else:
                        weight_array.append(2/((dx+1)*math.factorial(dx)))

                    edge_array.append([source, target])
            
            vertices_names = np.unique(edge_array)    

            return vertices_names, np.asarray(edge_array), np.asarray(weight_array)

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


def missing_patterns(data, dx=3, dy=1, taux=1, tauy=1, return_fraction=True, return_missing=True, probs=False, tie_precision=None):
    """
    Identifies ordinal patterns (permutations) not occurring in data\\ [#amigo]_.

    Parameters
    ----------
    data : array, return of :func:`ordpy.ordinal_distribution`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           or the ordinal patterns and probabilities 
           returned by :func:`ordpy.ordinal_distribution` with `return_missing=True`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    return_fraction : boolean
                      if `True`, returns the fraction of missing ordinal patterns 
                      relative to the total number of ordinal patterns for given 
                      values of **dx** and **dy** (default: `True`); if `False`, 
                      returns the raw number of missing patterns.
    return_missing : boolean 
                     if `True`, returns the missing ordinal patterns in **data** 
                     (default: `True`); if `False`, it only returns the 
                     fraction/number of these missing patterns.
    probs : boolean
            If `True`, assumes **data** to be the return of 
            :func:`ordpy.ordinal_distribution` with `return_missing=True`.
            If `False`, **data** is expected to be a one- or two-dimensional 
            array (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : tuple
       Tuple containing an array and a float indicating missing 
       ordinal patterns and their relative fraction (or raw number).
       
    Examples
    --------
    >>> missing_patterns([4,7,9,10,6,11,3], dx=2, return_fraction=False)
    (array([], shape=(0, 2), dtype=int64), 0)
    >>> 
    >>> missing_patterns([4,7,9,10,6,11,3,5,6,2,3,1], dx=3, return_fraction=True, return_missing=False)
    0.3333333333333333
    >>>
    >>> missing_patterns(ordinal_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True), dx=2, probs=True)
    (array([], shape=(0, 2), dtype=int64), 0.0)
    >>>
    >>> missing_patterns(ordinal_distribution([4,7,9,10,6,11,3], dx=3, return_missing=True), dx=3, probs=True)
    (array([[0, 2, 1],
            [1, 2, 0],
            [2, 1, 0]]), 0.5)
    """

    if probs == True:
        states, probs = data
    else:
        states, probs = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=True, tie_precision=tie_precision)

    missing_args   = np.argwhere(probs==0).flatten()
    missing_states = states[missing_args]

    if return_fraction==True:
        if return_missing:
            return missing_states, len(missing_args)/math.factorial(dx*dy)
        else:
            return len(missing_args)/math.factorial(dx*dy)
    else:
        if return_missing:
            return missing_states, len(missing_args)
        else:
            return len(missing_args)


def missing_links(data, dx=3, dy=1, return_fraction=True, return_missing=True, tie_precision=None):
    """
    Identifies transitions between ordinal patterns (permutations) not occurring 
    in data. (These transtitions correspond to directed edges in ordinal 
    networks\\ [#pessa2019]_.) Assumes overlapping windows and unitary embedding 
    delay. In case dx>1 and dy>1, both 'horizontal' and 'vertical' connections are 
    considered (see :func:`ordpy.ordinal_network`).
    
    Parameters
    ----------
    data : array, return of :func:`ordpy.ordinal_network`
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]` 
           or the nodes, edges and edge weights
           returned by :func:`ordpy.ordinal_network`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3)
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series 
         (default: 1).
    return_fraction : boolean
                      if `True`, returns the fraction of missing links among 
                      ordinal patterns relative to the total number of possible 
                      links (transitions) for given values of **dx** and **dy** 
                      (default: `True`). If `False`, returns the raw number of 
                      missing links.
    return_missing : boolean 
                     if True, returns the missing links in **data**; if False, it 
                     only returns the fraction/number of these missing links.
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : tuple
       Tuple containing an array and a float indicating missing links and their 
       relative fraction (or raw number).
    
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
            ['2|1|0', '2|1|0']], dtype='<U5'), 11)
    """
    def setdiff(a, b):
        """
        Searches for elements (subarrays) in `a` that are not contained in `b` [*]_. 

        Parameters
        ----------    
        a : tuples, lists or arrays
            Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
            \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]`.
        b : tuples, lists or arrays
            Array in the format :math:`[[x_{21}, x_{22}, x_{23}, \\ldots, x_{2m}], 
            \\ldots, [x_{n1}, x_{n2}, x_{n3}, ..., x_{nm}]]`.
        
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



    if not return_missing==False:
        if len(data)==3 and type(data[0][0])==np.str_:
            _, data_links, _ = data
        else:
            _, data_links, _ = ordinal_network(data, dx, dy, tie_precision=tie_precision)
            
        _, all_links, _ = random_ordinal_network(dx=dx, dy=dy)
        missing_links   = setdiff(all_links, data_links)

        if return_fraction==True:
            return missing_links, len(missing_links)/len(all_links)
        else:
            return missing_links, len(missing_links)

    else: #not return missing links
        if len(data)==3 and type(data[0][0])==np.str_:
            _, data_links, _ = data
        else:
            _, data_links, _ = ordinal_network(data, dx, dy, tie_precision=tie_precision)

        #assumes dx==1 means time series data, although the algorithm works
        #for images too in this case and te results are different.
        if dy==1:
            all_links = math.factorial(dx)*dx

            if return_fraction==True:
                return (all_links-len(data_links))/all_links
            else:
                return all_links-len(data_links)

        else: #not return missing links for image
            _, all_links, _ = random_ordinal_network(dx=dx, dy=dy)
            missing_links   = setdiff(all_links, data_links)

            if return_fraction==True:
                return len(missing_links)/len(all_links)
            else:
                return len(missing_links)


def minimum_complexity_entropy(dx=3, dy=1, size=100):
    """
    Generates data corresponding to values of normalized permutation entropy and
    statistical complexity which delimit the lower boundary in the 
    complexity-entropy causality plane\\ [#rosso_curvas]_:sup:`,`\\ [*]_.
    
    Parameters
    ----------
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series 
         (default: 1).
    size : int
           The length of the array returned containing pairs of values of 
           permutation entropy and statistical complexity.    
           
    Returns
    -------
     : array
       Values of normalized permutation entropy and statistical complexity 
       belonging to the lower limiting curve of the complexity-entropy 
       causality plane.

    Notes
    -----
    .. [*] This function was adapted from Sebastian Sippel et al., 
           statcomp: Statistical Complexity and Information Measures for 
           Time Series Analysis, version 0.1.0. (Computer Software in R).
    
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
    N            = math.factorial(dx*dy)
    prob_params  = np.linspace(1/N, 1, num=size-1)
    uniform_dist = np.full(N, 1/N)
    
    hc_ = []
    for i in range(size-1):
        probabilities    = np.full(shape=N, fill_value=(1-prob_params[i])/(N-1))
        probabilities[0] = prob_params[i]

        h = permutation_entropy(probabilities, dx, dy, probs=True)

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
    Generates data corresponding to values of normalized permutation entropy and 
    statistical complexity which delimit the upper boundary in the 
    complexity-entropy  causality plane\\ [#rosso_curvas]_:sup:`,`\\ [*]_.
    
    Parameters
    ----------
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series 
         (default: 1).
    m : int
        The length of the returned array containing values of permutation entropy
        and statistical complexity is given by 
        :math:`[(d_x \\times d_y)!-1] \\times m`. 
           
    Returns
    -------
     : array
       Values of normalized permutation entropy and statistical complexity 
       belonging to the upper limiting curve of the complexity-entropy causality 
       plane.

    Notes
    -----
    .. [*] This function was adapted from Sebastian Sippel et al., 
           statcomp: Statistical Complexity and Information Measures for 
           Time Series Analysis, version 0.1.0. (Computer Software in R).
    
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
    
    N              = math.factorial(dx*dy)
    hlist_, clist_ = np.zeros(shape=(N-1,m)), np.zeros(shape=(N-1,m))

    for i in range(N-1):

        p             = np.zeros(shape=N)
        uniform_dist  = np.full(N, 1/N)
        prob_params   = np.linspace(0, 1/N, num=m)

        for k in range(len(prob_params)):
            p[0] = prob_params[k]
            for j in range(1,N-i):
                p[j] = (1-prob_params[k])/(N-i-1)

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


def weighted_permutation_entropy(data, dx=3, dy=1, taux=1, tauy=1, base=2, normalized=True, tie_precision=None):
    """
    Calculates Shannon entropy using a weighted ordinal distribution
    obtained from data\\ [#fadlallah]_.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    base : str, int
           Logarithm base in Shannon's entropy. Either 'e' or 2 (default: 'e').
    normalized: boolean
                If `True`, weighted permutation entropy is normalized by its 
                maximum value. If `False`, it does not (default: `True`).
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).
    Returns
    -------
     : float
       The value of weighted permutation entropy.
    
    Examples
    --------
    >>> weighted_permutation_entropy([4,7,9,10,6,11,3], dx=2, base=2, normalized=False)
    0.9125914261094841
    >>>
    >>> weighted_permutation_entropy([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2, base=2)
    0.2613186822347165
    >>>
    >>> weighted_permutation_entropy([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=2, dy=2, taux=2, base='e')
    0.22725901135766047
    """
    try:
        ny, nx = np.shape(data)
        data   = np.array(data)
    except:
        nx     = np.shape(data)[0]
        ny     = 1
        data   = np.array([data])

    if tie_precision is not None: 
        data = np.round(data, tie_precision)

    if ny==1: #time series
        partitions = np.apply_along_axis(func1d       = np.lib.stride_tricks.sliding_window_view, 
                                         axis         = 1, 
                                         arr          = data, 
                                         window_shape = (dx+(dx-1)*(taux-1),)
                                        )
        partitions = partitions[::, ::, ::taux].reshape(-1, dx*dy) 
    else: #image
        partitions = np.lib.stride_tricks.sliding_window_view(x            = data, 
                                                              window_shape = (dy+(dy-1)*(tauy-1), dx+(dx-1)*(taux-1))
                                                             )
        partitions = partitions[::, ::, ::tauy, ::taux].reshape(-1, dx*dy)
    
    if dx*dy>6 and nx*ny>1000: #a special approach is employed when dealing with large arrays and large embedding dimensions.
        weights                   = partitions.var(axis=1)
        flat_symbolic_sequence    = np.apply_along_axis(np.argsort, 1, partitions)

        str_symbols               = np.apply_along_axis(np.array2string, 1, flat_symbolic_sequence, separator='')
        idx_sort                  = np.argsort(str_symbols)
        _, idx_start              = np.unique(str_symbols[idx_sort], return_index=True)
        sets_same_symbol          = np.split(idx_sort, idx_start[1:])

        total_variance_per_symbol = np.fromiter((weights[symbol_].sum() for symbol_ in sets_same_symbol), dtype='float')        
    
    else: #simpler approach employed to small series (less than 1000 observations) and small embedding dimensions (dx*dy < 7, < ~5000 permutations).
        weights                   = partitions.var(axis=1)
        flat_symbolic_sequence    = np.apply_along_axis(np.argsort, 1, partitions)
        symbols                   = np.unique(flat_symbolic_sequence, return_counts=False, axis=0)
        total_variance_per_symbol = np.fromiter((weights[(flat_symbolic_sequence==symbol_).all(axis=1)].sum() for symbol_ in symbols), dtype='float')
        
    probs = total_variance_per_symbol/total_variance_per_symbol.sum()

    logfunc = np.log if base=='e' else np.log2
    s       = -np.sum(probs[probs>0]*logfunc(probs[probs>0]))

    if normalized==True:        
        smax = logfunc(float(math.factorial(dx*dy)))
        return s/smax
    else:
        return s


def fisher_shannon(data, dx=3, dy=1, taux=1, tauy=1, probs=False, tie_precision=None):
    """
    Calculates permutation entropy\\ [#bandt_pompe]_ and Fisher 
    information\\ [#olivares]_ using an ordinal distribution obtained 
    from data.
    
    Parameters
    ----------
    data : array
           Array object in the format :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or  :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`
           or an ordinal probability distribution (such as the ones returned by 
           :func:`ordpy.ordinal_distribution`).
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).
    probs : boolean
            If `True`, it assumes **data** is an ordinal probability distribution.
            If `False`, **data** is expected to be a one- or two-dimensional array
            (default: `False`). 
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).

    Returns
    -------
     : tuple
       Values of normalized permutation entropy and Fisher information.

    Examples
    --------
    >>> fisher_shannon([4,7,9,10,6,11,3], dx=2)
    (0.9182958340544896, 0.028595479208968325)
    >>>
    >>> p = ordinal_distribution([4,7,9,10,6,11,3], dx=2, return_missing=True)[1]
    >>> fisher_shannon(p, dx=2, probs=True)
    (0.9182958340544896, 0.028595479208968325)
    >>>
    >>> fisher_shannon([[1,2,1],[8,3,4],[6,7,5]], dx=2, dy=2)
    (0.32715643797829735, 1.0)
    >>>
    >>> fisher_shannon([1/3, 1/15, 2/15, 3/15, 4/15, 0], dx=3, probs=True)
    (0.8314454838586239, 0.19574180681374548)
    >>>
    >>> fisher_shannon([[1,2,1,4],[8,3,4,5],[6,7,5,6]], dx=3, dy=2)
    (0.2107070115500801, 1.0)    
    """  
    if not probs:
        if dx*dy>10:
            warnings.warn("Large values of embedding dimension may lead to memory allocation problems.")                
        else: pass
        
        symbols, probabilities = ordinal_distribution(data, dx, dy, taux, tauy, return_missing=False, tie_precision=tie_precision)
        h                      = permutation_entropy(probabilities, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)        
        
        #In case not all possible permutations occur in **data**, we assign zero probability the non-occurring states. 
        #This is necessary to calculate the Fisher information measure
        if (probabilities==1.).sum()==1.:
            return h, 1.
        elif len(probabilities)==dx*dy:
            all_probs = probabilities
        else:
            all_symbols           = np.asarray(list(itertools.permutations(range(dx*dy))), dtype='int')
            all_probs             = np.full(math.factorial(dx*dy), 0.)

            _, ncols              = all_symbols.shape
            dtype                 = {'names':['f{}'.format(i) for i in range(ncols)], 'formats': ncols*[all_symbols.dtype]}
            logifilter            = np.isin(all_symbols.view(dtype), symbols.view(dtype)).flatten()
            all_probs[logifilter] = probabilities            

    else:
        warnings.warn("Be mindful the correct calculation of Fisher information depends on all possible permutations " + 
                      "given the embedding dimensions dx and dy, not only on the permutations occurring in data.")

        all_probs = np.asarray(data)
        h         = permutation_entropy(all_probs, dx, dy, taux, tauy, probs=True, tie_precision=tie_precision)
        
        if (all_probs==1.).sum()==1.:
            return h, 1.

    return h, np.sum(np.diff(np.sqrt(all_probs[::-1]))**2.)/2  #F = F_0 \sum_{i = 1}^{N - 1} (\sqrt{p_{i+1}} - \sqrt{p_{i}})^2, with F_0=1/2

    
def permutation_js_distance(data, dx=3, dy=1, taux=1, tauy=1, base='e', normalized=True, tie_precision=None):
    """
    Calculates the permutation Jensen-Shannon distance\\ [#zunino2022]_ between 
    multiple time series (or channels of a multivariate time series) or 
    multiple matrices (images) using ordinal distributions obtained from data.
    
    Parameters
    ---------
    data : array 
           Array object containing  arrays in the format 
           :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`
           or :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`.
    dx : int
         Embedding dimension (horizontal axis) (default: 3).
    dy : int
         Embedding dimension (vertical axis); it must be 1 for time series 
         (default: 1).
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).    
    base : str, int
           Logarithm base in Shannon's entropy. Either 'e' or 2 (default: 2).
    normalized: boolean
                If `True`, the permutation Jensen-Shannon distance is normalized by  
                the logarithm of 2. (default: `True`). If `False`, it is not.
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).
                    
    Returns
    -------
    : float
      Value of permutation Jensen-Shannon distance.
    
    Examples
    --------
    >>> permutation_js_distance([[1,2,6,5,4], [1,2,1]], dx=3, base=2, normalized=False)
    0.677604543245723
    >>>
    >>> permutation_js_distance([[1,2,6,5,4], [1,2,1]], dx=3, base='e', normalized=False)
    0.5641427870206323
    >>>
    >>> permutation_js_distance([[[1,2,6],[4,2,6]], [[1,2,6],[1,2,4]]], dx=2, dy=2)
    0.7071067811865476
    >>>
    >>> permutation_js_distance([[[1,2,6,5,4],[4,2,6,3,4]], [[1,2,6],[1,2,4]]], dx=2, dy=2)
    0.809715420521041
    """
    if dx*dy>10:
        warnings.warn("Large values of embedding dimension may lead to memory allocation problems.")                
    else: pass
    
    logfunc                         = np.log if base=="e" else np.log2
    num_of_channels                 = len(data)
    sum_of_entropies, probabilities = 0., []

    for channel_ in data:
        symbols, probs     = ordinal_distribution(channel_, dx=dx, dy=dy, taux=taux, tauy=tauy, return_missing=True, tie_precision=tie_precision)

        symbols_as_strings = np.apply_along_axis(np.array2string, 1, symbols, separator='')
        mask_probs         = np.apply_along_axis(np.argsort,      0, symbols_as_strings)

        probabilities.append(probs[mask_probs])
        sum_of_entropies  += -np.sum(probs[probs>0]*logfunc(probs[probs>0]))/num_of_channels

    sum_of_distros = np.sum(probabilities, axis=0)/num_of_channels
    entropy_of_sum = -np.sum(sum_of_distros[sum_of_distros>0]*logfunc(sum_of_distros[sum_of_distros>0]))
    distance       = entropy_of_sum - sum_of_entropies
    is_it_nan      = np.isnan(np.sqrt(distance))        

    if is_it_nan==True or num_of_channels==1: #added `is_it_nan` after having problems calculating the square root of very small negative numbers caused by floating-point representation problems in Python.
        return 0                              #the second boolean check basically covers the case only one series is provided.
    else:
        if normalized==True:
            return np.sqrt(distance/logfunc(num_of_channels))
        else:
            return np.sqrt(distance)


def permutation_contrasts(data, taux=1, tie_precision=None):
    """
    Calculates the four pattern (permutation) contrasts [#bandt]_ 
    of a time series (with dx=3) using an ordinal distribution 
    obtained from data.
    
    Parameters
    ---------
    data : array 
           Array object containing a pair of arrays in the format 
           :math:`[x_{1}, x_{2}, x_{3}, \\ldots ,x_{n}]`.
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).
                    
    Returns
    -------
     : tuple
       Tuple containing values corresponding to the pattern (permutation) 
       contrasts: up-down balance, persistence, rotational-asymmetry,
       and up-down scaling.
       
    Examples
    --------
    >>> permutation_contrasts([4,7,9,10,6,11,3])
    (0.4, 0.4, 0.6, -0.2)
    >>>
    >>> permutation_contrasts([4,7,9,10,6,11,3], taux=2)
    (0.0, 0.666, -0.333, 0.333)
    >>>
    >>> permutation_contrasts([5,2,3,4,2,7,4])
    (0.2, 0.2, 0.0, 0.0)
    >>>
    >>> permutation_contrasts([5, 2, 3, 4, 2, 7, 4], taux=2)
    (0.0, 0.666, 0.333, 0.333)
    """
    dict_symbols = dict([('[012]', 0),('[021]', 1),('[102]', 2),
                         ('[120]', 3),('[201]', 4),('[210]', 5)])
    
    symbols                = ordinal_sequence(data, dx=3, dy=1, taux=taux, tauy=1, overlapping=True, tie_precision=tie_precision)
    symbols, symbols_count = np.unique(symbols, return_counts=True, axis=0)
    probabilities          = symbols_count/symbols_count.sum()
    symbols                = np.apply_along_axis(np.argsort, 1, symbols) #sort the permutations once again so they correspond to the ones used by Bandt and Wittfeld.


    symbols         = np.apply_along_axis(np.array2string, 1, symbols, separator='')
    indices         = np.fromiter(map(dict_symbols.get, symbols), dtype='int')
    all_probs       = np.fromiter((probabilities[np.argwhere(indices==type_).flatten()].sum() for type_ in range(6)), dtype='float')            

    up_down_balance = all_probs[0] - all_probs[5] #p012 - p210
    persistence     = all_probs[0] + all_probs[5] #p012 + p210
    rot_asymetry    = all_probs[2] + all_probs[3] - all_probs[1] - all_probs[4] #p102 + p120 - p021 - p201
    updown_scaling  = all_probs[1] + all_probs[2] - all_probs[3] - all_probs[4] #p021 + p102 - p120 - p201
    
    return up_down_balance, persistence, rot_asymetry, updown_scaling


def smoothness_structure(data, taux=1, tauy=1, tie_precision=None):
    """
    Calculates the smoothness and curve structure\\ [#bandt_wittfeld]_ 
    of an image (with embedding parameters dx=dy=2) using an ordinal
    distribution obtained from data.
    
    Parameters
    ---------
    data : array
           Array object in the format 
           :math:`[[x_{11}, x_{12}, x_{13}, \\ldots, x_{1m}],
           \\ldots, [x_{n1}, x_{n2}, x_{n3}, \\ldots, x_{nm}]]`.
    taux : int
           Embedding delay (horizontal axis) (default: 1).
    tauy : int
           Embedding delay (vertical axis) (default: 1).    
    tie_precision : None, int
                    If not `None`, **data** is rounded with `tie_precision`
                    decimal numbers (default: `None`).
                    
    Returns
    -------
     : tuple
       Values of smoothness and curve structure.
       
    Examples
    --------
    >>> smoothness_structure([[1,2,1,4],[8,3,4,5],[6,7,5,6]])
    (0.0, 0.333)
    >>>
    >>> smoothness_structure([[1,2,1,4],[8,3,4,5],[6,7,5,6]], taux=2, tauy=1)
    (-0.0833, 0.75)
    >>>
    >>> smoothness_structure([[1,2,1,4],[8,3,4,5],[6,7,5,6]], taux=2, tauy=2)
    (-0.333, 1.0)
    """
    #These permutations correspond to a different sorting process laid out in the paper by Bandt and Wittfeld
    dict_types  = dict([('[0123]', 0),('[0132]', 1),('[0213]', 0),('[0231]', 2),
                        ('[0312]', 1),('[0321]', 2),('[1023]', 1),('[1032]', 0),
                        ('[1203]', 1),('[1230]', 2),('[1302]', 0),('[1320]', 2),
                        ('[2013]', 2),('[2031]', 0),('[2103]', 2),('[2130]', 1),
                        ('[2301]', 0),('[2310]', 1),('[3012]', 2),('[3021]', 1),
                        ('[3102]', 2),('[3120]', 0),('[3201]', 1),('[3210]', 0)])

    symbols                = ordinal_sequence(data, dx=2, dy=2, taux=taux, tauy=tauy, overlapping=True, tie_precision=tie_precision)
    symbols                = symbols.reshape(-1, 4) #reshape the symbols array to facilitate the application of np.unique below
    symbols, symbols_count = np.unique(symbols, return_counts=True, axis=0)
    probabilities          = symbols_count/symbols_count.sum()
    symbols                = np.apply_along_axis(np.argsort, 1, symbols) #resort the permutations so they correspond to the ones used by Bandt and Wittfeld.

    probs_types = np.full(shape=3, fill_value=0.)
    symbols     = np.apply_along_axis(np.array2string, 1, symbols, separator='')
    indices     = np.fromiter(map(dict_types.get, symbols), dtype='int')
    probs_types = np.fromiter((probabilities[np.argwhere(indices==type_).flatten()].sum() for type_ in range(3)), dtype='float')

    smoothness  = probs_types[0] - 1/3
    structure   = probs_types[1] - probs_types[2]

    return smoothness, structure

