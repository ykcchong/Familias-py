"""Pedigree plotting via matplotlib + networkx (replacing R's kinship2).

Usage::

    from familias import FamiliasPedigree, plot_familias_pedigree
    p = FamiliasPedigree(id=..., dadid=..., momid=..., sex=...)
    plot_familias_pedigree(p)
"""
from __future__ import annotations
from typing import Optional


def plot_familias_pedigree(ped, ax=None, show: bool = True, title: Optional[str] = None):
    try:
        import matplotlib.pyplot as plt
        import networkx as nx
    except ImportError as e:  # pragma: no cover
        raise ImportError(
            "Plotting requires matplotlib and networkx. "
            "Install with `pip install familias[plot]`."
        ) from e

    G = nx.DiGraph()
    for pid, sex in zip(ped.id, ped.sex):
        G.add_node(pid, sex=sex)
    for i, pid in enumerate(ped.id):
        if ped.findex[i] > 0:
            G.add_edge(ped.id[ped.findex[i] - 1], pid)
        if ped.mindex[i] > 0:
            G.add_edge(ped.id[ped.mindex[i] - 1], pid)

    # Compute generations top-down (founders -> generation 0)
    gen = {n: 0 for n in G.nodes if G.in_degree(n) == 0}
    changed = True
    while changed:
        changed = False
        for u, v in G.edges:
            new = gen.get(u, 0) + 1
            if gen.get(v, -1) < new:
                gen[v] = new
                changed = True
    # Group by generation, lay out horizontally
    layers = {}
    for n, g in gen.items():
        layers.setdefault(g, []).append(n)
    pos = {}
    max_per_layer = max(len(layers[g]) for g in layers)
    for g, nodes in layers.items():
        for k, n in enumerate(sorted(nodes)):
            pos[n] = (k - len(nodes) / 2.0, -g)

    if ax is None:
        _, ax = plt.subplots(figsize=(max(6, max_per_layer * 1.3),
                                       max(4, len(layers) * 1.3)))

    # Draw nodes by sex (squares for male, circles for female)
    males = [n for n, d in G.nodes(data=True) if d["sex"] == "male"]
    females = [n for n, d in G.nodes(data=True) if d["sex"] == "female"]
    nx.draw_networkx_nodes(G, pos, nodelist=males, node_shape="s",
                            node_color="lightblue", node_size=1200, ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=females, node_shape="o",
                            node_color="lightpink", node_size=1200, ax=ax)
    nx.draw_networkx_edges(G, pos, arrows=False, ax=ax)
    nx.draw_networkx_labels(G, pos, font_size=8, ax=ax)

    if title:
        ax.set_title(title)
    ax.set_axis_off()
    if show:
        import matplotlib.pyplot as plt
        plt.show()
    return ax
