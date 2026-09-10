"""
greedy_color.py -- greedy vertex colouring of the cluster adjacency graph, so bordering clusters on a map get different colours.

Called from _inc/conus_map.do and _inc/state_hangover_map.do as:
    python greedy_color.py <nodes_csv> <edges_csv> <out_csv>
Reads:  <nodes_csv> (column id); <edges_csv> (columns a, b; may be empty)
Writes: <out_csv> with id, color (1-based)
Notes:  uses igraph's vertex_coloring_greedy; the number of colours is whatever the greedy order needs, and the Stata caller maps it onto the tableau palette.
"""

import sys

import pandas as pd
import igraph as ig

def main(nodes_csv: str, edges_csv: str, out_csv: str) -> None:
    nodes = pd.read_csv(nodes_csv)
    ids = nodes["id"].tolist()
    pos = {v: i for i, v in enumerate(ids)}

    try:
        edges = pd.read_csv(edges_csv)
        elist = [(pos[a], pos[b]) for a, b in zip(edges["a"], edges["b"])
                 if a in pos and b in pos and a != b]
    except (pd.errors.EmptyDataError, KeyError):
        elist = []

    g = ig.Graph(n=len(ids), edges=elist, directed=False)
    g.simplify(multiple=True, loops=True)

    colors = list(g.vertex_coloring_greedy())
    ncol = (max(colors) + 1) if colors else 1

    pd.DataFrame({"id": ids, "color": [c + 1 for c in colors]}).to_csv(
        out_csv, index=False)
    print(f"greedy_color: {len(ids)} clusters colored with {ncol} colors "
          f"({len(elist)} adjacency edges)")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit("usage: greedy_color.py <nodes_csv> <edges_csv> <out_csv>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
