# Notes on cgranges

These notes should help to grok Heng Li's [cgranges](https://github.com/lh3/cgranges) interval tree data structure, which stores the intervals in a flat array (sorted by begin position) and interprets this array as an [implicit](https://en.wikipedia.org/wiki/Implicit_data_structure) binary search tree. That interpretation has some tricky details, but once we work them out, we can easily treat the binary tree as a [textbook interval tree](https://en.wikipedia.org/wiki/Interval_tree#Augmented_tree) in which each node is augmented with the maximum interval end position amongst itself and its descendants.

We begin by reproducing the [explanatory code comment](https://github.com/lh3/cgranges/blob/master/cpp/IITree.h):

```
/* Suppose there are N=2^(K+1)-1 sorted numbers in an array a[]. They
 * implicitly form a complete binary tree of height K+1. We consider leaves to
 * be at level 0. The binary tree has the following properties:
 *
 * 1. The lowest k-1 bits of nodes at level k are all 1. The k-th bit is 0.
 *    The first node at level k is indexed by 2^k-1. The root of the tree is
 *    indexed by 2^K-1.
 *
 * 2. For a node x at level k, its left child is x-2^(k-1) and the right child
 *    is x+2^(k-1).
 *
 * 3. For a node x at level k, it is a left child if its (k+1)-th bit is 0. Its
 *    parent node is x+2^k. Similarly, if the (k+1)-th bit is 1, x is a right
 *    child and its parent is x-2^k.
 *
 * 4. For a node x at level k, there are 2^(k+1)-1 nodes in the subtree
 *    descending from x, including x. The left-most leaf is x&~(2^k-1) (masking
 *    the lowest k bits to 0).
 *
 * When numbers can't fill a complete binary tree, the parent of a node may not
 * be present in the array. The implementation here still mimics a complete
 * tree, though getting the special casing right is a little complex. There may
 * be alternative solutions.
 *
 * As a sorted array can be considered as a binary search tree, we can
 * implement an interval tree on top of the idea. We only need to record, for
 * each node, the maximum value in the subtree descending from the node.
 */
 ```

Let's look at a couple of these binary search trees of integers starting from zero. Consider the number of each node as its rank in the sorted array, so the item with the smallest interval begin position has rank 0.

N=7, K=2:
```
       3
     1   5
    0 2 4 6
```

N=15, K=3:
```
            7
       3        11
     1   5    9    13
    0 2 4 6  8 10 12 14
```

Notice that each parent node in the sorted array is flanked by its left and right subtree contiguously; the node ranked one lower is the rightmost descendant of its left subtree, and the node ranked one higher is the leftmost descendant of its right subtree. Its immediate left and right children reside in the middle of their respective stretches of the array, and so on fractally. The code comment provides closed formulae for calculating the ranks of these related nodes, so that there's no need to store pointers between them.

In cgranges, the implicit tree layout is determined solely by the *number* of sorted items indexed, otherwise independent of their actual interval positions. This causes at least two complications. First, given a node with interval begin position B, other nodes with that same begin B can exist in *either or both* of its subtrees, even if its *immediate* children have different begin positions. Its parent or any other relatives could potentially share B as well: just work backwards from the extreme where the whole "sorted" array has the same begin. This can be a thorn for various tree algorithms, although it ends up not really affecting cgranges itself.

* Items with colliding begin positions are contiguous in the flat sorted array, so it's much easier to find them there than by tree traversal; just scan items adjacent to the node of interest. Rank calculations can then distinguish whether they're inside or outside the node's subtree based on the rank difference. The iitii code has examples.

The other tricky complication...

### Imaginary nodes

Consider again the K=3 tree, but suppose we have only 10 items in our sorted array instead of the full 15.

```
            7
       3        ?
     1   5    9    ?
    0 2 4 6  8  ? ?  ?
```

The predetermined [full binary tree](https://web.cecs.pdx.edu/~sheard/course/Cs163/Doc/FullvsComplete.html) structure doesn't rearrange to accommodate the "odd" number of nodes, so we have to deal with ranks 10-14 as *imaginary* nodes. Almost half of the tree could be imaginary (when the root is the highest-ranked real node). Keep "rank" and "array index" separate: the commented rank traversal formulae work fine on imaginary node ranks, so long as we're careful to avoid attempting to access the array beyond its length.

Or how about N=21:

```
                       15
            7                    ?
       3        11           19     ?
     1   5    9    13     17    ?      ?
    0 2 4 6  8 10 12 14 16 18 20  ?  ?   ?
```

Notice that, if we permit ourselves to traverse only between real nodes, this tree has *three* distinct, disconnected components, with node 20 orphaned unto itself.

The order in which the binary tree fills up (left, parent, right) implies useful facts about imaginary nodes.

* A real node can have imaginary descendants in its right subtree only.
* An imaginary node can have real descendants in its left subtree only.
* The border of the real and imaginary nodes can be traced by starting from the rightmost real leaf and following parents up to the root. This path may include imaginary nodes, and fully connects the tree.

cgranges uses the last fact to propagate the interval end bounds from the orphaned parts of the tree up to the root.

### "There may be alternative solutions."

[Brodal, Fagerberg & Jacob (2001) §3.3](https://tidsskrift.dk/brics/article/download/21696/19132) give a different construction for implicit search trees which avoids the complication of imaginary nodes. Briefly, write *N* as a sum of powers of two, *N* = 2<sup>*p*<sub>0</sub></sup> + 2<sup>*p*<sub>1</sub></sup> + 2<sup>*p*<sub>2</sub></sup> + ... (where the *p*'s are the places of *N*'s binary 1 digits). Then, consider item 0 in the array to be an index node for the implicit, full tree of items [1, 2<sup>*p*<sub>0</sub></sup>); item 2<sup>*p*<sub>0</sub></sup> indexes the full tree of items [2<sup>*p*<sub>0</sub></sup>+1, 2<sup>*p*<sub>0</sub></sup>+2<sup>*p*<sub>1</sub></sup>), and so on. Search proceeds by first finding the rightmost index node less than the query, then searching the adjacent full tree (and all trees whose index nodes are equal to the query). This approach could be used for interval queries by augmenting each index node with the maximum interval end position of itself and its adjacent tree.
