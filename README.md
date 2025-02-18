# Customizable Contraction Hierarchies: Metric-Independent Contraction Hierarchy

 (though the metric-independent contraction hierarchy is actually not customizable)


## General design approach
Rather than storing each "core graph" separately and then taking the union to obtain the Contraction Hierarchy, we use a single representation which wraps the underlying Node/Edge types with corresponding enums that annotate contraction and shortcut status respectively, and iteration of contraction/shortcut, in such a way that any previous core graph can be computed in a constant memory footprint.


## Relevant sections of the CCH paper


### src/lib.rs (contraction hierarchy)
> Given a fixed weight w, one can exploit that in many applications it is sufficient to only preserve all shortest path distances [21]. Weighted vertex contraction of a vertex v in the graph G is defined as the operation of removing v and inserting (a minimum number of shortcuts) among the neighbors of v to obtain a graph G′ such that distG(x, y) = distG′ (x, y) for all vertices x != v and y != v. To compute G′, one iterates over all pairs of neighbors x, y of v increasing by distG(x, y). For each pair one checks whether a xy-path of length distG(x, y) exists in G\{v}, i. e., one checks whether removing v destroys the xy-shortest path. This check is called witness search [21] and the xy-path is called witness, if it exists. If a witness is found, the considered vertex pair is skipped and no shortcut added. Otherwise, if an edge {x, y} already exists, its weight is decreased to distG(x, y), or a new shortcut edge with that weight is added to G. This new shortcut edge is considered in witness searches for subsequent neighbor pairs as part of G. If shortest paths are not unique, it is important to iterate over the pairs increasing by distG(x, y), because otherwise more edges than strictly necessary can be inserted: Shorter shortcuts can make longer shortcuts superfluous.

src/nested_dissection.rs
> To support metric-independence, we therefore use nested dissection orders as suggested in [6] or ND-orders for short. An order π for G is computed recursively by determining a balanced separator S of minimum cardinality that splits G into two parts induced by the vertex sets A and B. The vertices of S are assigned to π(n − |S|) … π(n) in an arbitrary order. Orders πₐ and π_B are computed recursively and assigned to π(1) … π(|A|) and π(|A| + 1) … π(|A| + |B|), respectively. The base case of the recursion is reached when the subgraphs are empty.
