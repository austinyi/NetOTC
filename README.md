# NetOTC
 
This is the official repository for the paper titled "Alignment and Comparison of Directed Networks via
Transition Couplings of Random Walks" by Bongsoo Yi, Kevin O'Connor, Kevin McGoff, Andrew B. Nobel. It contains a Matlab implementation of NetOTC as well as the other Network OT methods used as baselines. Moreover, code for reproducing experimental results and examples from the paper is also included.

## Running NetOTC
In order to run the NetOTC algorithm on weighted adjacency matrices `A1` and `A2` with cost matrix `c`, first transform to transition matrices as follows:
```
P1 = adj_to_trans(A1);
P2 = adj_to_trans(A2);
```
Then, we can run NetOTC using the ExactOTC algorithm by calling
```
[gotc_cost, gotc, node_alignment] = exact_otc(P1, P2, c);
```
For moderately sized networks, one may want to use the more efficient EntropicOTC algorithm by calling
```
[gotc_cost, gotc, node_alignment] = entropic_otc(P1, P2, c, L, T, xi, sink_iter, get_sd);
```
The parameters `L` and `T` are used in the ApproxTCE step that approximately evaluates the current transition coupling; the parameters `xi` and `sink_iter` are used in the EntropicTCI step that uses Sinkhorn's algorithm to improve the transition coupling. Larger values of each of these parameters will give a better approximation of ExactOTC at the expense of increased runtime. Finally, the parameter `get_sd` indicates whether the node alignment should be computed and returned. Setting this to 0 when the node alignment is not needed can save some time for larger networks.

In order to compute the mass assigned to a pair of edges `(x1, x2)` and `(y1, y2)`, compute `node_alignment(x1, y1)*gotc((x1, y1), (x2, y2))`. For an example of this, please refer to `pc_align_exp.m`.

## Running Experiments and Examples
Code for reproducing the examples and experimental results described in the paper may be found in the folder `experiments`. Note that you will have to modify any directories in the scripts before you run them. Moreover, the isomorphism and classification experiment code should be run with `longleaf=0`.

## Citing this Repository
If you wish to cite our work, please use the following BibTeX code:
```
@article{yi2025alignment,
  title={Alignment and comparison of directed networks via transition couplings of random walks},
  author={Yi, Bongsoo and O'Connor, Kevin and McGoff, Kevin and Nobel, Andrew B},
  journal={Journal of the Royal Statistical Society Series B: Statistical Methodology},
  pages={qkae085},
  year={2024},
  doi = {10.1093/jrsssb/qkae085}
}
```
