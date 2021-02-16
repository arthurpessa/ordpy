.. |logo1| image:: https://img.shields.io/pypi/v/ordpy?style=plastic   :alt: PyPI 
   :target: https://pypi.org/project/ordpy/
   :scale: 100%
.. |logo2| image:: https://img.shields.io/github/license/arthurpessa/ordpy?style=plastic   :alt: GitHub 
   :target: https://github.com/arthurpessa/ordpy/blob/master/LICENSE
   :scale: 100%
.. |logo3| image:: https://img.shields.io/pypi/dm/ordpy?style=plastic   :alt: PyPI - Downloads
   :target: https://pypi.org/project/ordpy/
   :scale: 100%

|logo1| |logo2| |logo3|

ordpy: A Python Package for Data Analysis with Permutation Entropy and Ordinal Network Methods
===============================================================================================

``ordpy`` is a pure Python module [#pessa2021]_ that implements data analysis methods based
on Bandt and Pompe's [#bandt_pompe]_ symbolic encoding scheme.

If you have used ``ordpy`` in a scientific publication, we would appreciate citations to the following reference [#pessa2021]_:

- A. A. B. Pessa, H. V. Ribeiro, `ordpy: A Python package for data analysis with permutation entropy and ordinal 
  network methods <https://arxiv.org/abs/2102.06786>`_, arXiv: 2102.06786 (2021).

.. code-block:: bibtex
    
   @misc{pessa2021ordpy,
    title         = {ordpy: A Python module implementing permutation entropy and ordinal network methods for data analysis}, 
    author        = {Arthur A. B. Pessa and Haroldo V. Ribeiro},
    year          = {2021},
    eprint        = {2102.06786},
    archivePrefix = {arXiv},
   }

``ordpy`` implements the following data analysis methods:

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

For more detailed information about the methods implemented in ``ordpy``, please 
consult its `documentation <https://arthurpessa.github.io/ordpy/_build/html/index.html>`_.

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



References
==========

.. [#pessa2021] Pessa, A. A., & Ribeiro, H. V. (2021). ordpy: A Python package
   for data analysis with permutation entropy and ordinal networks methods. 
   arXiv: 2102.06786.

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

.. [#rosso_curvas] Martin, M. T., Plastino, A., & Rosso, O. A. (2006). 
   Generalized Statistical Complexity Measures: Geometrical and 
   Analytical Properties, Physica A, 369, 439–462.
