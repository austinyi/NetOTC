# GraphOTC
 
This is the official repository for the paper titled "Graph Optimal Transport with Transition Couplings of Random Walks" by O'Connor et al. It contains a Matlab implementation of GraphOTC as well as the other Graph OT methods used as baselines. Moreover, code for reproducing experimental results and examples from the paper is also included.

## Running GraphOTC
In order to run the GraphOTC algorithm on weighted adjacency matrices `A1` and `A2` with cost matrix `c`, first transform to transition matrices as follows:
```
P1 = adj_to_trans(A1);
P2 = adj_to_trans(A2);
```
Then, we can run GraphOTC using the ExactOTC algorithm by calling
```
[gotc_cost, gotc, node_alignment] = exact_otc(P1, P2, c);
```
For moderately sized graphs, one may want to use the more efficient EntropicOTC algorithm by calling
```
[gotc_cost, gotc, node_alignment] = entropic_otc(P1, P2, c, L, T, xi, sink_iter, get_sd);
```
The parameters `L` and `T` are used in the ApproxTCE step that approximately evaluates the current transition coupling; the parameters `xi` and `sink_iter` are used in the EntropicTCI step that uses Sinkhorn's algorithm to improve the transition coupling. Larger values of each of these parameters will give a better approximation of ExactOTC at the expense of increased runtime. Finally, the parameter `get_sd` indicates whether the node alignment should be computed and returned. Setting this to 0 when the node alignment is not needed can save some time for larger graphs.

In order to the mass assigned to a pair of edges `(x1, x2)` and `(y1, y2)`, compute `node_alignment(x1, y1)*gotc((x1, y1), (x2, y2))`. For an example of this, please refer to `pc_align_exp.m`.

## Running Experiments and Examples
Code for reproducing the examples and experimental results described in the paper may be found in the folders `examples` and `experiments`. Note that you will have to modify any directories in the scripts before you run them. Moreover, the classification experiment code should be run with `longleaf=0`.

## Citing this Repository
If you wish to cite our work, please use the following BibTeX code:
```
@article{o2021graph,
  title={Graph Optimal Transport with Transition Couplings of Random Walks},
  author={O'Connor, Kevin and Yi, Bongsoo and McGoff, Kevin and Nobel, Andrew B},
  journal={arXiv preprint arXiv:2106.07106},
  year={2021}
}
```
