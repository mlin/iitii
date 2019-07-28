#include <vector>
#include <limits>
#include <algorithm>
#include <assert.h>

// Base template for the internal representation of a node within an implicit interval tree
// User should not care about this; subclass instantiations may add more members for more-
// exotically augmented classes of implicit interval trees
template<
    typename Pos,              // numeric position type (e.g. uint32_t, double)
    typename Item,             // arbitrary type of items to be indexed
    Pos get_beg(const Item&),  // function to get interval begin position from item
    Pos get_end(const Item&)   // function to get interval end position from item
>   // get_beg & get_end should just access members of Item (i.e. cache-local fetch)
struct iit_node_base {
    static const Pos npos = std::numeric_limits<Pos>::max();    // reserved constant for invalid Pos

    Item item;
    Pos inside_max_end;   // max end of this & subtree (as in textbook augmented interval tree)

    iit_node_base(const Item& item_)
        : item(item_)
        , inside_max_end(get_end(item))
        {}
    inline Pos beg() const {
        return get_beg(item);
    }
    inline Pos end() const {
        return get_end(item);
    }
    bool operator<(const iit_node_base<Pos, Item, get_beg, get_end>& rhs) const {
        auto lbeg = beg(), rbeg = rhs.beg();
        if (lbeg == rbeg) {
            return end() < rhs.end();
        }
        return lbeg < rbeg;
    }
};

// Base template for an implicit interval tree, with internal repr Node<Pos, Item, ...>
// Users should not use this directly, but instantiate iit (below) or other subclasses
template<typename Pos, typename Item, typename Node>
class iit_base {
protected:
    // aliases to help keep the Pos, Rank, and Level concepts straight
    typedef std::size_t Rank;   // rank of a node, its index in the sorted array (or beyond, if dark)
    typedef std::size_t Level;  // level in tree

    static const Rank nrank = std::numeric_limits<Rank>::max();

    std::vector<Node> nodes;  // array of Nodes sorted by beginning position
    size_t full_size;         // size of the full binary tree containing the nodes; liable to be
                              // as large as 2*nodes.size()-1 including implicit "dark" nodes

    Rank root;
    Level root_level;         //  = K in cgranges

    // compute a node's level, the # of 1 bits below the lowest 0 bit
    inline Level level(Rank node) const {
        assert(node < full_size);
        Level ans = __builtin_ctzll(~node);  // bitwise-negate & count trailing zeroes
        #ifndef NDEBUG
        Level chk;
        for (chk=0; node&1; chk++, node>>=1);
        assert(ans == chk);
        #endif
        return ans;
    }

    // get node's parent, or nrank if called on the root
    inline Rank parent(Rank node) const {
        assert(node < full_size);
        if (node == root) {
            return nrank;
        }
        Level lv = level(node);
        Rank ofs = Rank(1) << lv;
        assert(node >= ofs-1);
        if (((node>>(lv+1)) & 1)) {
            assert(node >= ofs);
            return node-ofs;
        }
        return node+ofs;
    }

    // get node's left child, or nrank if called on a leaf
    inline Rank left(Rank node) const {
        Level lv = level(node);
        return lv > 0 ? node - (Rank(1) << (lv-1)) : nrank;
    }

    // get node's right child, or nrank if called on a leaf
    inline Rank right(Rank node) const {
        Level lv = level(node);
        return lv > 0 ? node + (Rank(1) << (lv-1)) : nrank;
    }

    // top-down overlap scan for [qbeg,qend). return # of nodes visited
    size_t scan(Rank subtree, Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        if (subtree == nrank) {
            return 0;
        }
        assert(subtree < full_size);
        if (subtree >= nodes.size()) {
            // When we arrive at a dark node, its right subtree must also be completely dark, owing
            // to the order in which the binary tree fills up recursively (left, parent, right). So
            // we just need to descend into the left branch, which must eventually lead to real
            // node(s) perhaps via a chain of dark ones.
            // TODO: there should be a closed formula for the length of this left chain of dark
            //       nodes given subtree, n.size(), and full_size. or at very least we could
            //       precompute a lookup table.
            return 1 + scan(left(subtree), qbeg, qend, ans);
        }

        size_t cost = 1;
        const Node& n = nodes[subtree];
        if (n.inside_max_end > qbeg) {  // something in current subtree extends into/over query
            cost += scan(left(subtree), qbeg, qend, ans);
            Pos nbeg = n.beg();
            if (nbeg < qend) {          // this node isn't already past query
                if (n.end() > qbeg) {   // this node overlaps query
                    ans.push_back(n.item);
                }
                cost += scan(right(subtree), qbeg, qend, ans);
            }
        }
        return cost;
    }

public:
    template<typename InputIterator>
    iit_base(InputIterator first, InputIterator last)
        : root_level(0)
        , root(std::numeric_limits<Rank>::max())
    {
        // store the items
        std::for_each(first, last, [&](const Item& it) { nodes.push_back(Node(it)); });

        // compute the implied tree geometry
        for (root_level = 0, full_size = 0; full_size < nodes.size();
             ++root_level, full_size = (size_t(1)<<(root_level+1)) - 1);
        root = (Rank(1) << root_level) - 1;

        if (nodes.size()) {
            // sort the nodes by interval
            std::sort(nodes.begin(), nodes.end());

            // Memoize the path from the rightmost leaf up to the root. This will trace the border
            // between the nodes present and dark nodes (if any), which we'll refer to in indexing
            // below. Some of these border nodes may be dark.
            std::vector<Rank> right_border_nodes({
                nodes.size() - (2 - nodes.size() % 2)  // rightmost leaf in nodes (not dark)
            });
            while (right_border_nodes.back() != root) {
                right_border_nodes.push_back(parent(right_border_nodes.back()));
            }

            // bottom-up indexing
            Pos right_border_ime = nodes[right_border_nodes[0]].inside_max_end;
            for (Level lv=1; lv <= root_level; ++lv) {
                // for each in nodes on this level
                size_t x = size_t(1)<<(lv-1), step = x<<2;
                for (Rank n = (x<<1)-1; n < nodes.size(); n += step) {
                    // figure inside_max_end
                    Pos ime = nodes[n].end();
                    ime = std::max(ime, nodes[left(n)].inside_max_end);
                    if (right(n) < nodes.size()) {
                        ime = std::max(ime, nodes[right(n)].inside_max_end);
                    } else {
                        // right child is dark; take the last border observation
                        ime = std::max(ime, right_border_ime);
                    }
                    assert(ime != Node::npos);
                    nodes[n].inside_max_end = ime;

                    if (n == right_border_nodes[lv]) {
                        // track inside_max_end of the non-dark nodes on the border
                        right_border_ime = ime;
                    }
                }
            }
        }
    }

    // overlap query; fill ans and return query cost (number of tree nodes visited)
    virtual size_t overlap(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        ans.clear();
        return scan(root, qbeg, qend, ans);
    }

    // overlap query, return vector of results
    std::vector<Item> overlap(Pos qbeg, Pos qend) const {
        std::vector<Item> ans;
        overlap(qbeg, qend, ans);
        return ans;
    }
};

// Basic implicit interval tree
template<
    typename Pos,              // numeric position type (e.g. uint32_t, double)
    typename Item,             // arbitrary type of items to be indexed
    Pos get_beg(const Item&),  // function to get interval begin position from item
    Pos get_end(const Item&)   // function to get interval end position from item
>   // get_beg & get_end should just access members of Item (i.e. cache-local fetch)
class iit : public iit_base<Pos, Item, iit_node_base<Pos, Item, get_beg, get_end>> {
public:
    template<typename InputIterator>
    iit(InputIterator first, InputIterator last)
        : iit_base<Pos, Item, iit_node_base<Pos, Item, get_beg, get_end>>(first, last)
        {}
};

template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
struct iitii_node : public iit_node_base<Pos, Item, get_beg, get_end> {
    // Additional augment values for iitii nodes, which help us know when we can stop climbing in
    // the bottom-up search.
    //
    // outside_max_end of node n is the maximum m.end() of all nodes m outside of n & its subtree
    //     with m.beg() < n.beg(); -infinity if there are no such nodes.
    // outside_min_beg of node n is the minimum m.beg() of all nodes m outside of n & its subtree
    //     with m.beg() >= n.beg(); infinity if there are no such nodes.
    //
    // Suppose during a query for [qbeg, qend) we climb up to a node n with,
    //   (i) n.outside_max_end <= qbeg; AND
    //  (ii) qend <= n.outside_min_beg
    // By (i), any node outside n & subtree with beg < n's cannot overlap the query. By (ii), any
    // node outside n & subtree with beg >= n's cannot overlap the query. This exhausts all nodes
    // outside of n & subtree, so we don't need to climb past n.
    Pos outside_max_end;
    Pos outside_min_beg;

    iitii_node(const Item& item_)
        : iit_node_base<Pos, Item, get_beg, get_end>(item_)
        , outside_max_end(std::numeric_limits<Pos>::min())
        , outside_min_beg(std::numeric_limits<Pos>::max())
        {}
};

template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
class iitii : public iit_base<Pos, Item, iitii_node<Pos, Item, get_beg, get_end>> {
    using Node = iitii_node<Pos, Item, get_beg, get_end>;
    using super = iit_base<Pos, Item, Node>;
    using typename super::Rank, typename super::Level;
    using super::left, super::right, super::parent;
    using super::nodes;
    using super::full_size, super::nrank;
    using super::level, super::root, super::root_level;

    // leaf prediction model: the leaves are the even-ranked nodes; so predict,
    //   2*floor(w[1]*qbeg + w[0])
    float w[2] = { 0.0f, 0.0f };
    Rank predict_leaf(Pos qbeg) const {
        float halfrank = w[1]*float(qbeg) + w[0];
        assert(halfrank == halfrank);
        Rank r = 2*Rank(std::max(0.0f,halfrank));
        const auto nsz = nodes.size();
        return r < nsz ? r : (nsz - (2 - nsz%2));
    }

    void train() {
        // Simple linear regression of rank/2 ~ node.beg() over all the leaf nodes
        // Sure, let's write it ourselves, that sounds like a great idea =)

        double sum_beg = 0, sum_halfrank = 0;
        size_t leaf_count = 0;
        for (Rank rank = 0; rank < nodes.size(); rank += 2) {
            sum_beg += nodes[rank].beg();
            sum_halfrank += rank/2;
            ++leaf_count;
        }

        if (leaf_count == 0) {
            w[0] = w[1] = 0.0f;
            return;
        }

        const double mean_beg = sum_beg/leaf_count,
                     mean_halfrank = sum_halfrank/leaf_count;

        double cov = 0, var = 0;
        for (Rank rank = 0; rank < nodes.size(); rank += 2) {
            const double beg_err = nodes[rank].beg() - mean_beg;
            cov += beg_err*(rank/2 - mean_halfrank);
            var += beg_err*beg_err;
        }

        const double m = cov / var;
        w[1] = float(m);
        w[0] = float(mean_halfrank - m*mean_beg);
        assert(w[0] == w[0] && w[1] == w[1]);
    }

    inline Rank leftmost_child(Rank subtree) const {
        Level k = level(subtree);
        auto ofs = (Rank(1)<<k) - 1;
        assert(subtree >= ofs);
        return subtree - ofs;
    }

    inline Rank rightmost_child(Rank subtree) const {
        Level k = level(subtree);
        auto ofs = (Rank(1)<<k) - 1;
        assert(subtree + ofs < full_size);
        return subtree + ofs;
    }

public:
    template<typename InputIterator>
    iitii(InputIterator first, InputIterator last)
        : super(first, last)
    {
        if (nodes.size()) {
            // compute running max_end along the sorted array, which we'll look up while computing
            // outside_max_end below
            std::vector<Pos> running_max_end { nodes[0].end() };
            for (Rank n = 1; n < nodes.size(); ++n) {
                running_max_end.push_back(std::max(running_max_end[n-1], nodes[n].end()));
            }

            // Fill outside_min_beg and outside_max_end. Neatly, with the cgranges beg-sorted array
            // and running_max_end, we can find each value independently in constant time (modulo
            // beg position collisions)
            for (Rank n = 0; n < nodes.size(); ++n) {
                Node& node = nodes[n];
                Rank l = leftmost_child(n), r = rightmost_child(n);

                // outside_min_beg is beg() of the node ranked one higher than n's rightmost child
                node.outside_min_beg = r < nodes.size()-1 ? nodes[r+1].beg()
                                                          : std::numeric_limits<Pos>::max();
                if (l>0) {
                    if (nodes[l-1].beg() == node.beg()) {
                        // corner case: nodes to the left of n's subtree can have the same beg as n
                        // and outside_min_beg is defined on nodes with beg >= n's.
                        node.outside_min_beg = node.beg();
                    }

                    // outside_max_end is the running_max_end of the highest-ranked node ranked
                    // below n's leftmost child & has beg strictly below n's
                    Rank leq = l-1;
                    while (nodes[leq].beg() == node.beg()) {
                        if (leq == 0) {
                            break;
                        }
                        --leq;
                    }
                    assert(nodes[leq].beg() <= node.beg());
                    node.outside_max_end = nodes[leq].beg() < node.beg()
                                                ? running_max_end[leq]
                                                : std::numeric_limits<Pos>::min();
                }

                assert(node.outside_min_beg >= node.beg());
            }

            // train the model for leaf rank ~ beg
            train();
        }
    }

    size_t overlap(Pos qbeg, Pos qend, std::vector<Item>& ans) const override {
        // ask model which leaf we should begin our bottom-up climb at
        Rank subtree = predict_leaf(qbeg);
        assert(level(subtree) == 0);
        size_t climb_cost = 0;
        // climb until we hit the root or our necessary & sufficient criteria are met
        while (subtree != root) {
            if (subtree < nodes.size() &&
                nodes[subtree].outside_max_end <= qbeg &&
                qend <= nodes[subtree].outside_min_beg) {
                break;
            }
            subtree = parent(subtree);
            // TODO: scheme to skip through chains of dark nodes along the border
            ++climb_cost;
        }

        // scan the subtree for query results 
        ans.clear();
        return super::scan(subtree, qbeg, qend, ans) + climb_cost;
    }

    using super::overlap;
};
