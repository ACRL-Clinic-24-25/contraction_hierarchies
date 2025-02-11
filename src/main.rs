#![allow(dead_code)]

use contraction_hierarchies::CH;
use petgraph::graph::NodeIndex;
use petgraph::Graph;

fn main() {
    let mut graph = Graph::<(), ()>::new();

    // Create 3 nodes in a path: 0 -> 1 -> 2
    let n0 = graph.add_node(());
    let n1 = graph.add_node(());
    let n2 = graph.add_node(());

    graph.add_edge(n0, n1, ());
    graph.add_edge(n1, n2, ());

    let contraction_order = [1].map(NodeIndex::new).into_iter();

    println!("Starting graph:\n{:#?}", graph);
    let mut ch = CH::new(graph, contraction_order);

    let core_graph = ch.core_graph(1);
    let contraction_hierarchy = ch.contraction_hierarchy();

    println!("core graph:\n{:#?}", core_graph);
    println!("contraction hierarchy:\n{:#?}", contraction_hierarchy);

    // let mut graph = Graph::<(), ()>::new();
    // let mut nodes = Vec::new();
    // let n = 10;
    // for _ in 0..n {
    //     nodes.push(graph.add_node(()));
    // }
    // for (a, b) in nodes
    //     .iter()
    //     .enumerate()
    //     .flat_map(|(i, a)| nodes.iter().skip(i).map(|b| (*a, *b)).take(2))
    // {
    //     graph.add_edge(a, b, ());
    // }
    // let contraction_order = nested_dissection_contraction_order(graph.clone()).into_iter();
    //
    // let core_graph = ch.core_graph(contraction_order.len()-1);
    // let contraction_hierarchy = ch.contraction_hierarchy();
    // println!("core graph:\n{:#?}", core_graph);
    // println!("contraction hierarchy:\n{:#?}", contraction_hierarchy);
}
