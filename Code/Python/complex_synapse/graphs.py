"""Tools for plotting models as graphs
"""
from __future__ import annotations
from numbers import Number
from typing import Dict, List, Optional, Tuple

# from typing import List, Optional, Sequence, Tuple

import matplotlib as mpl
# import matplotlib.animation as mpla
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
# import sl_py_tools.containers as cn
import sl_py_tools.iter_tricks as it
import sl_py_tools.numpy_tricks.markov as ma
# import sl_py_tools.matplotlib_tricks as mpt

from . import synapse_mem as sm
from . import optimise as opt
from . import identify as idfy
# =============================================================================


def  model_to_graph(model: sm.SynapseModel, thresh: float) -> nx.DiGraph:
    """Create a directed graph from a synapse model

    Parameters
    ----------
    model : SynapseModel
        Synapse model

    Returns
    -------
    graph : DiGraph
        Graph describing model
    """
    peq = model.peq()
    graph = nx.DiGraph()
    for i, weight, prob in it.zenumerate(model.weight, peq):
        graph.add_node(i, kind=weight, value=prob)
    inds = np.unravel_index(ma.indices.offdiag_inds(model.nstate),
                            model.plast.shape[-2:])
    for k, sig in enumerate(model.signal):
        for i, j in zip(*inds):
            val = model.plast[k, i, j]
            if val > thresh:
                graph.add_edge(i, j, kind=sig, value=val)
    return graph


def id_model_to_graph(model: idfy.SynapseIdModel) -> nx.DiGraph:
    """Create a directed graph from a synapse identification model

    Parameters
    ----------
    model : SynapseIdModel
        Synapse model

    Returns
    -------
    graph : DiGraph
        Graph describing model
    """
    peq = model.peq()
    graph = nx.DiGraph()
    for i, weight, prob in it.zenumerate(model.readout, peq):
        graph.add_node(i, kind=weight, value=prob)
    inds = np.unravel_index(ma.indices.offdiag_inds(model.nstate),
                            model.plast.shape[-2:])
    for k in range(model.nplast):
        for i, j in zip(*inds):
            graph.add_edge(i, j, kind=k, value=model.plast[k, i, j])
    return graph


def param_model_to_graph(model: opt.SynapseParamModel) -> nx.DiGraph:
    """Create a directed graph from a parameterised synapse model

    Parameters
    ----------
    model : SynapseParamModel
        Synapse model

    Returns
    -------
    graph : DiGraph
        Graph describing model
    """
    peq = model.peq()
    param = model.get_params().reshape((model.nplast, -1))
    return param_to_graph(param, peq, model.weight, model.topology)


def param_to_graph(param: la.lnarray, peq: la.lnarray, weight: la.lnarray,
                   topology: opt.TopologyOptions) -> nx.DiGraph:
    """Create a directed graph from a parameterised synapse model

    Parameters
    ----------
    param : la.lnarray (P,k), k in [M(M-1), M, M-1, 1]
        Independent parameters of model
    peq : la.lnarray (M,)
        Steady state distribution
    weight : la.lnarray (M,)
        Synaptic weight in each state
    topology : opt.TopologyOptions
        Encapsulation of model class

    Returns
    -------
    graph : DiGraph
        Graph describing model
    """
    nstate = len(peq)
    nplast = param.shape[0]
    graph = nx.DiGraph()
    for i, wgt, prob in it.zenumerate(weight, peq):
        graph.add_node(i, kind=wgt, value=prob)
    # (P,)(2,)(Q,)
    inds = [ma.indices.param_subs(nstate, **topology.directed(k))
            for k in range(nplast)]
    for plast in it.zenumerate(inds, param):
        for i, j, rate in zip(*plast[1], plast[2]):
            graph.add_edge(i, j, kind=plast[0], value=rate)
    return graph


def collect_node_values(graph: nx.DiGraph, key: str) -> np.ndarray:
    """Collect values of node attributes"""
    return np.array([node[key] for node in graph.nodes.values()])


def collect_edge_values(graph: nx.DiGraph, key: str) -> np.ndarray:
    """Collect values of edge attributes"""
    return np.array([edge[key] for edge in graph.edges.values()])


def collect_node_colours(graph: nx.DiGraph, key: str) -> Dict[str, np.ndarray]:
    """Collect values of node attributes for the colour"""
    vals = collect_node_values(graph, key)
    vmin, vmax = vals.min(), vals.max()
    return {'node_color': vals, 'vmin': vmin, 'vmax': vmax}


def collect_edge_colours(graph: nx.DiGraph, key: str) -> Dict[str, np.ndarray]:
    """Collect values of edge attributes for the colour"""
    vals = collect_edge_values(graph, key)
    vmin, vmax = vals.min(), vals.max()
    return {'edge_color': vals, 'edge_vmin': vmin, 'edge_vmax': vmax}


def draw_serial(graph: nx.DiGraph, axs: Optional[mpl.axes.Axes] = None
                ) -> Tuple[mpl.collections.PathCollection,
                           List[mpl.patches.FancyArrowPatch]]:
    """draw a serial model"""
    axs = ag.default_eval(axs, plt.gca)
    cmap = plt.get_cmap('coolwarm')
    pos = {n: (n, 0) for n in graph.nodes}

    node_col = collect_node_colours(graph, 'kind')
    node_siz = np.sqrt(collect_node_values(graph, 'value')) * 300
    edge_col = collect_edge_colours(graph, 'kind')
    edge_wid = collect_edge_values(graph, 'value') * 3

    node_col.update(ax=axs, cmap=cmap, edgecolors='k')
    edge_col.update(ax=axs, edge_cmap=cmap)
    nodes = nx.draw_networkx_nodes(graph, pos, node_size=node_siz, **node_col)
    edges = nx.draw_networkx_edges(graph, pos, width=edge_wid, **edge_col)
    for edge in edges:
        edge.set_connectionstyle('Arc3', rad=-.5)
    return nodes, edges
