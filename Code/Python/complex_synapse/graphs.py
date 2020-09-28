"""Tools for plotting models as graphs
"""
from __future__ import annotations

import operator as opr
from numbers import Number
from typing import (Any, Callable, Dict, Hashable, Iterable, List, Optional,
                    Tuple, TypeVar, Union)

import matplotlib as mpl
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import numpy_linalg as la
import sl_py_tools.arg_tricks as ag
import sl_py_tools.containers as cn
# import sl_py_tools.iter_tricks as it
import sl_py_tools.numpy_tricks.markov as ma
import sl_py_tools.graph_tricks as gt
# import sl_py_tools.graph_plots as gp
import sl_py_tools.options_classes as op

import complex_synapse.identify as idfy
import complex_synapse.optimise as opt
# import complex_synapse.options as op
import complex_synapse.synapse_base as sb
import complex_synapse.synapse_mem as sm

# import sl_py_tools.matplotlib_tricks as mpt


# =============================================================================
# Graph builders
# =============================================================================


def _base_to_graph(model: sb.SynapseBase,
                   node_id: Getter = opr.attrgetter('weight'),
                   edge_id: Getter = opr.attrgetter('signal'),
                   node_val: Getter = opr.methodcaller('peq')
                   ) -> gt.MultiDiGraph:
    """Create a directed graph from a synapse model.

    Parameters
    ----------
    model : SynapseModel
        Synapse model.

    Returns
    -------
    graph : nx.DiGraph
        Graph describing model.
    """
    params = ma.params.gen_mat_to_params(model.plast, (0,) * model.nplast)
    return gt.param_to_graph(params, node_val(model),
                             node_id(model), edge_id(model))


def model_to_graph(model: sm.SynapseModel) -> gt.MultiDiGraph:
    """Create a directed graph from a synapse model.

    Parameters
    ----------
    model : SynapseModel
        Synapse model.

    Returns
    -------
    graph : nx.DiGraph
        Graph describing model.
    """
    return _base_to_graph(model, node_id=opr.attrgetter('weight'),
                          edge_id=opr.attrgetter('signal'),
                          node_val=opr.methodcaller('peq'))


def idmodel_to_graph(model: idfy.SynapseIdModel) -> gt.MultiDiGraph:
    """Create a directed graph from a synapse model.

    Parameters
    ----------
    model : SynapseModel
        Synapse model.

    Returns
    -------
    graph : nx.DiGraph
        Graph describing model.
    """
    return _base_to_graph(model, node_id=opr.attrgetter('readout'),
                          edge_id=lambda x: np.arange(x.nplast),
                          node_val=opr.attrgetter('initial'))


def param_model_to_graph(model: opt.SynapseParamModel) -> gt.MultiDiGraph:
    """Create a directed graph from a parameterised synapse model.

    Parameters
    ----------
    model : SynapseParamModel
        Synapse model.

    Returns
    -------
    graph : nx.DiGraph
        Graph describing model.
    """
    peq = model.peq()
    param = model.get_params()
    return gt.param_to_graph(param, peq, model.weight, None, model.topology)


def param_to_graph(param: la.lnarray, values: la.lnarray, kinds: la.lnarray,
                   topology: opt.TopologyOptions) -> gt.MultiDiGraph:
    """Create a directed graph from a parameters of a synapse model.

    Parameters
    ----------
    param : la.lnarray (PQ,), Q in [M(M-1), M, M-1, 1]
        Independent parameters of model.
    values : la.lnarray (M,)
        Value assigned to each state.
    kinds : la.lnarray (M,)
        Type label for each state.
    topology : opt.TopologyOptions
        Encapsulation of model class.

    Returns
    -------
    graph : DiGraph
        Graph describing model.
    """
    return gt.param_to_graph(param, values, kinds, None, topology)


# =============================================================================
# Graph attributes
# =============================================================================


def has_node_attr(graph: nx.DiGraph, key: str) -> bool:
    """Test for existence of node attributes.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph with nodes whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    has : bool
        `True` if every node has the attribute, `False` otherwise.
    """
    return all(key in node for node in graph.nodes.values())


def has_edge_attr(graph: nx.DiGraph, key: str) -> bool:
    """Test for existence of edge attributes.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph with edges whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    has : bool
        `True` if every edge has the attribute, `False` otherwise.
    """
    return all(key in edge for edge in graph.edges.values())


def node_attrs(graph: nx.DiGraph, key: str) -> Dict[Node, Any]:
    """Dictionary of node attribute values.

    Parameters
    ----------
    graph : nx.DiGraph (N,E)
        Graph with nodes whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    attrs : Dict[Node, Any]
        Dictionary of node attribute values, keyed by node.
    """
    return {node: values[key] for node, values in graph.nodes.items()}


def edge_attrs(graph: nx.DiGraph, key: str) -> Dict[Edge, Any]:
    """Dictionary of edge attribute values.

    Parameters
    ----------
    graph : nx.DiGraph (N,E)
        Graph with edges whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    attrs : Dict[Tuple[Node, Node], Any]
        Dictionary of edge attribute values, keyed by edge (node tuple).
    """
    return {edge: values[key] for edge, values in graph.edges.items()}


def node_attr_vec(graph: nx.DiGraph, key: str) -> np.ndarray:
    """Collect values of node attributes.

    Parameters
    ----------
    graph : nx.DiGraph (N,E)
        Graph with nodes whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    vec : np.ndarray (N,)
        Vector of node attribute values.
    """
    return np.array([node[key] for node in graph.nodes.values()])


def edge_attr_vec(graph: nx.DiGraph, key: str) -> np.ndarray:
    """Collect values of edge attributes.

    Parameters
    ----------
    graph : nx.DiGraph (N,E)
        Graph with edges whose attributes we want.
    key : str
        Name of attribute.

    Returns
    -------
    vec : np.ndarray (E,)
        Vector of edge attribute values.
    """
    return np.array([edge[key] for edge in graph.edges.values()])


def edge_attr_matrix(graph: nx.DiGraph, key: str, fill: Number = 0.
                        ) -> np.ndarray:
    """Collect values of edge attributes in a matrix.

    Parameters
    ----------
    graph : nx.DiGraph (N,E)
        Graph with edges whose attributes we want.
    key : str
        Name of attribute to use for matrix elements.
    fill : Number
        Value given to missing edges.

    Returns
    -------
    mat : np.ndarray (N,N)
        Matrix of edge attribute values.
    """
    nodes = list(graph.nodes)
    mat = np.full((len(nodes),) * 2, fill)
    for edge, vals in graph.edges.items():
        ind = tuple(nodes.index(u) for u in edge)
        mat[ind] = vals[key]
    return mat


def collect_node_colours(graph: nx.DiGraph, key: str) -> Dict[str, np.ndarray]:
    """Collect values of node attributes for the colour

    Parameters
    ----------
    graph : nx.DiGraph
        Graph with nodes whose attributes we want.
    key : str
        Name of attribute to map to colour.

    Returns
    -------
    kwargs : Dict[str, np.ndarray]
        Dictionary of keyword arguments for `nx.draw_networkx_nodes` related to
        colour values: `{'node_color', 'vmin', 'vmax'}`.
    """
    vals = node_attr_vec(graph, key)
    vmin, vmax = vals.min(), vals.max()
    return {'node_color': vals, 'vmin': vmin, 'vmax': vmax}


def collect_edge_colours(graph: nx.DiGraph, key: str) -> Dict[str, np.ndarray]:
    """Collect values of edge attributes for the colour

    Parameters
    ----------
    graph : nx.DiGraph
        Graph with edges whose attributes we want.
    key : str
        Name of attribute to map to colour.

    Returns
    -------
    kwargs : Dict[str, np.ndarray]
        Dictionary of keyword arguments for `nx.draw_networkx_edges` related to
        colour values: `{'edge_color', 'edge_vmin', 'edge_vmax'}`.
    """
    vals = edge_attr_vec(graph, key)
    vmin, vmax = vals.min(), vals.max()
    return {'edge_color': vals, 'edge_vmin': vmin, 'edge_vmax': vmax}


# =============================================================================
# Options
# =============================================================================


# pylint: disable=too-many-ancestors
class GraphOptions(op.Options):
    """Options for drawing graphs.

    Parameters
    ----------
    topology : TopologyOptions
        Topology specifying options, for creating graphs.
    layout : Callable[DiGraph -> Dict[Node, ArrayLike]]
        Function to compute node positions.
    node_cmap : str|mpl.colors.Colormap
        Map `node['kind']` to node colour. `str` passed to `mpl.cm.get_cmap`.
    edge_cmap : str|mpl.colors.Colormap
        Map `edge['kind']` to edge colour. `str` passed to `mpl.cm.get_cmap`.
    size : float
        Scale factor between `node['value']` and node area.
    width : float
        Scale factor between `edge['value']` and edge thickness.
    rad : float
        Curvature of edges: ratio betweeen max perpendicular distance from
        straight line to curve and length of straight line.
    """
    map_attributes: op.Attrs = ('topology',)
    prop_attributes: op.Attrs = ('node_cmap', 'edge_cmap')
    # topology specifying options
    topology: opt.TopologyOptions
    layout: Layout
    _node_cmap: mpl.colors.Colormap
    _edge_cmap: mpl.colors.Colormap
    size: float
    width: float
    rad: float

    def __init__(self, *args, **kwds) -> None:
        self.topology = opt.TopologyOptions(serial=True)
        self.layout = linear_layout
        self.size = 600
        self.width = 5
        self.rad = 0.7
        self._node_cmap = mpl.cm.get_cmap('coolwarm')
        self._edge_cmap = mpl.cm.get_cmap('seismic')
        super().__init__(*args, **kwds)

    def set_node_cmap(self, value: Union[str, mpl.colors.Colormap]) -> None:
        """Set the colour map for nodes.

        Does nothing if `value` is `None`. Converts to `Colormap` if `str`.
        """
        if value is None:
            pass
        elif isinstance(value, str):
            self._node_cmap = mpl.cm.get_cmap(value)
        else:
            self._node_cmap = value

    def set_edge_cmap(self, value: Union[str, mpl.colors.Colormap]) -> None:
        """Set the colour map for edges.

        Does nothing if `value` is `None`. Converts to `Colormap` if `str`.
        """
        if value is None:
            pass
        elif isinstance(value, str):
            self._edge_cmap = mpl.cm.get_cmap(value)
        else:
            self._edge_cmap = value

    @property
    def node_cmap(self) -> mpl.colors.Colormap:
        """Get the colour map for nodes.
        """
        return self._node_cmap

    @property
    def edge_cmap(self) -> mpl.colors.Colormap:
        """Get the colour map for nodes.
        """
        return self._edge_cmap
# pylint: enable=too-many-ancestors


# =============================================================================
# Plot graph
# =============================================================================


def linear_layout(graph: nx.DiGraph, sep: sb.ArrayLike = (1., 0.)) -> NodePos:
    """Layout graph nodes in a line.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph whose nodes need laying out.
    sep : ArrayLike, optional
        Separation of nodes along line, by default `(1.0, 0.0)`.

    Returns
    -------
    pos : Dict[Node, np.ndarray]
        Dictionary of node ids -> position vectors.
    """
    sep = np.array(sep)
    return {node: pos * sep for pos, node in enumerate(graph.nodes)}


def draw_graph(graph: nx.DiGraph,
               pos: Union[NodePos, Layout, None] = None,
               axs: Optional[mpl.axes.Axes] = None,
               opts: Optional[GraphOptions] = None,
               **kwds) -> Tuple[NodePlots, DirectedEdgeCollection]:
    """Draw a synapse model's graph.

    Parameters
    ----------
    graph : nx.DiGraph
        The graph we are drawing.
    pos : NodePos|Layout|None, optional
        Dictionary of node positions, or a function to compute them,
        by default `None -> linear_layout`.
    axs : Axes|None, optional
        Axes on which we draw, by default `None -> plt.gca()`.
    opts : GraphOptions|None, optional
        Options for graph plot, by default `None -> GraphOptions()`.
    Other keywords passed to `opts`, if accepted, or to both
    `nx.draw_networkx_nodes` and `nx.draw_networkx_edges`.

    Returns
    -------
    nodes : PathCollection
        Collection of for the plots of the nodes.
    edges : DirectedEdgeCollection
        Collection of objects for the plots of the edges.
    """
    opts = ag.default_eval(opts, GraphOptions)
    opts.pop_my_args(kwds)
    axs = ag.default_eval(axs, plt.gca)
    pos = ag.default(pos, opts.layout)
    if callable(pos):
        pos = pos(graph)

    node_col = collect_node_colours(graph, 'kind')
    # note: PathCollection.size is area
    node_siz = node_attr_vec(graph, 'value') * opts.size
    edge_col = collect_edge_colours(graph, 'kind')
    edge_wid = edge_attr_vec(graph, 'value') * opts.width

    node_col.update(kwds, ax=axs, cmap=opts.node_cmap, edgecolors='k')
    edge_col.update(kwds, ax=axs, edge_cmap=opts.edge_cmap, node_size=node_siz,
                    connectionstyle=f'arc3,rad={-opts.rad}')

    nodes = nx.draw_networkx_nodes(graph, pos, node_size=node_siz, **node_col)
    edges = nx.draw_networkx_edges(graph, pos, width=edge_wid, **edge_col)

    return nodes, DirectedEdgeCollection(edges, graph)


def update_graph(nodes: NodePlots, edges: DirectedEdgeCollection,
                 node_siz: np.ndarray, edge_wid: np.ndarray,
                 opts: Optional[GraphOptions] = None, **kwds) -> None:
    """Update a synapse model's graph plot.

    Parameters
    ----------
    nodes : PathCollection
        The objects for the plots of the nodes.
    edges : List[FancyArrowPatch]
        The objects for the plots of the edges.
    node_siz : np.ndarray
        Sizes of the nodes (proportional to area)
    edge_wid : np.ndarray
        Widths of edges.
    opts : GraphOptions|None, optional
        Options for graph plot, by default `None -> GraphOptions()`.
    Other keywords passed to `opts`.
    """
    opts = ag.default_eval(opts, GraphOptions)
    opts.pop_my_args(kwds)
    nodes.set_sizes(node_siz * opts.size)
    edges.set_widths(edge_wid * opts.width)
    edges.set_node_sizes(node_siz * opts.size)


# =============================================================================
# Edge collection
# =============================================================================


class DirectedEdgeCollection:
    """A collection of edge plots"""
    _edges: Dict[Edge, EdgePlot]
    _node_ids: List[Node]

    def __init__(self, edges: Iterable[EdgePlot], graph: nx.DiGraph) -> None:
        self._edges = dict(zip(graph.edges, edges))
        self._node_ids = list(graph.nodes)

    def __len__(self) -> int:
        return len(self._edges)

    def __getitem__(self, key: Edge) -> EdgePlot:
        return self._edges[key]

    def __iter__(self) -> Iterable[Edge]:
        return iter(self._edges)

    def keys(self) -> Iterable[Edge]:
        """A view of edge dictionary keys"""
        return self._edges.keys()

    def values(self) -> Iterable[EdgePlot]:
        """A view of edge dictionary values"""
        return self._edges.values()

    def items(self) -> Iterable[Tuple[Edge, EdgePlot]]:
        """A view of edge dictionary items"""
        return self._edges.items()

    def set_widths(self, edge_wid: sb.ArrayLike) -> None:
        """Set line widths of edges"""
        edge_wid = cn.tuplify(edge_wid, len(self))
        for edge, wid in zip(self._edges.values(), edge_wid):
            edge.set_linewidth(wid)

    def set_node_sizes(self, node_siz: sb.ArrayLike) -> None:
        """Set line widths of edges"""
        node_siz = cn.tuplify(node_siz, len(self._node_ids))
        for edge_id, edge_plot in self._edges.items():
            src, dst = [self._node_ids.index(node) for node in edge_id]
            edge_plot.shrinkA = _to_marker_edge(node_siz[src], 'o')
            edge_plot.shrinkB = _to_marker_edge(node_siz[dst], 'o')


def _to_marker_edge(marker_size, marker):
    if marker in "s^>v<d":  # `large` markers need extra space
        return np.sqrt(2 * marker_size) / 2
    return np.sqrt(marker_size) / 2


class GraphPlots:
    """Class for plotting model as a graph.

    Parameters
    ----------
    graph : nx.DiGraph
        Graph object describing model. Nodes have attributes `kind` and
        `value`. Edges have attributes `kind`, `value` and `pind` (if the
        model was a `SynapseParamModel`).
    opts : GraphOptions|None, optional
        Options for plotting the graph, by default `None -> GraphOptions()`.
    Other keywords passed to `opt` or `complex_synapse.graph.draw_graph`.
    """
    nodes: NodePlots
    edges: DirectedEdgeCollection
    pinds: np.ndarray
    opts: GraphOptions

    def __init__(self, graph: nx.DiGraph,
                 opts: Optional[GraphOptions] = None, **kwds) -> None:
        """Class for plotting model as a graph.

        Parameters
        ----------
        graph : nx.DiGraph
            Graph object describing model. Nodes have attributes `kind` and
            `value`.  Edges have attributes `kind`, `value` and `pind` (if the
            model was a `SynapseParamModel`).
        opt : GraphOptions|None, optional
            Options for drawing the graph, by default `None -> GraphOptions()`.
        Other keywords passed to `opt` or `complex_synapse.graph.draw_graph`.
        """
        self.opts = ag.default_eval(opts, GraphOptions)
        self.opts.pop_my_args(kwds)
        self.nodes, self.edges = draw_graph(graph, opts=self.opts, **kwds)
        if has_edge_attr(graph, 'pind'):
            self.pinds = edge_attr_vec(graph, 'pind')
        else:
            self.pinds = np.arange(len(self.edges))

    def update(self, edge_vals: la.lnarray, node_vals: Optional[la.lnarray]
               ) -> None:
        """Update plots.

        Parameters
        ----------
        edge_vals : la.lnarray (E,)
            Transition probabilities, for edge line widdths.
        node_vals : None|la.lnarray (N,), optional
            Equilibrium distribution,for nodes sizes (area), by default `None`
            -> calculate from `params`.
        """
        edge_vals = la.asarray(edge_vals).ravel()
        # peq = ag.default_eval(peq, lambda: serial_eqp(params))
        self.nodes.set_sizes(node_vals * self.opts.size)
        self.edges.set_widths(edge_vals[self.pinds] * self.opts.width)
        self.edges.set_node_sizes(node_vals * self.opts.size)

    def update_from(self, graph: nx.DiGraph) -> None:
        """Update plots using a graph object.

        Parameters
        ----------
        graph : nx.DiGraph
            Graph object describing model. Nodes have attributes `kind` and
            `value`.  Edges have attributes `kind`, `value`.
        """
        params = edge_attr_vec(graph, 'value')
        peq = node_attr_vec(graph, 'value')
        self.update(params, peq)


# =============================================================================
# Aliases
# =============================================================================
Node = TypeVar('Node', int, str, Hashable)
Edge = Tuple[Node, Node]
NodePlots = mpl.collections.PathCollection
EdgePlot = mpl.patches.FancyArrowPatch
NodePos = Dict[Node, sb.ArrayLike]
Layout = Callable[[nx.DiGraph], NodePos]
Getter = Callable[[sb.SynapseBase], np.ndarray]
