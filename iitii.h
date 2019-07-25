#include <vector>
#include <limits>
#include <algorithm>
#include <assert.h>

// base template for the internal representation of a node within an implicit interval tree
//   Pos      : integral position type (e.g. uint32_t, uint64_t)
//   Item     : arbitrary type of items to be indexed
//   get_beg  : function to get interval begin position from Item
//   get_end  : function to get interval end position from Item
// For good performance, beg & end should just access members of Item (i.e. O(1) cache-local fetch)
// Subclass instantiations may add more augment values.
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
struct iit_node_base {
    Item item;
    Pos inside_max_end;   // max end of this & subtree (as in textbook augmented interval tree)

    iit_node_base(const Item& item_)
        : item(item_)
        , inside_max_end(std::numeric_limits<Pos>::max())
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

// base template for an implicit interval tree, with internal repr Node<Pos, Item, ...>
template<typename Pos, typename Item, typename Node>
class iit_base {
    // aliases to help keep the Pos, Rank, and Level concepts straight
    typedef std::size_t Rank;   // rank of a node, its index in the sorted array (or beyond, if dark)
    typedef std::size_t Level;  // level in tree

    static const Pos npos = std::numeric_limits<Pos>::max();    // reserved constant for invalid Pos
    static const Rank nrank = std::numeric_limits<Rank>::max();

    std::vector<Node> nodes;  // array of Nodes sorted by beginning position
    size_t full_size;         // size of the full binary tree containing the nodes; liable to be
                              // as large as 2*nodes.size()-1 including implicit "dark" nodes

    Rank root;
    Level root_level;          //  = K in cgranges


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
            // When we visit a dark node, the right child must also be dark, due to the order in
            // which the binary tree fills up (left, parent, right). So we just need to descend
            // into the left branch, which must eventually lead to real node(s) perhaps via a chain
            // of dark ones. If there were nothing there then the current dark node wouldn't exist.
            return 1 + scan(left(subtree), qbeg, qend, ans);
        }

        size_t cost = 1;
        const Node& n = nodes[subtree];
        if (n.inside_max_end > qbeg) {     // something in current subtree extends into/over query
            cost += scan(left(subtree), qbeg, qend, ans);
            Pos nbeg = n.beg();
            if (nbeg < qend) {             // this node isn't already past query
                if (n.end() > qbeg) {  // this node overlaps query
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

            // memoize the rightmost leaf and its ancestors, which we'll use during indexing to
            // skip dark areas of the tree
            std::vector<Rank> rightmost_anc({
                nodes.size() - (2 - nodes.size() % 2)  // rightmost leaf in nodes (not dark)
            });
            while (rightmost_anc.back() != root) {
                rightmost_anc.push_back(parent(rightmost_anc.back()));
            }

            // bottom-up indexing
            Pos rightmost_ime;
            for (Level lv=1; lv <= root_level; ++lv) {
                if (rightmost_anc[lv-1] < nodes.size()) {
                    // inside_max_end of the rightmost leaf's ancestor on lv-1.
                    // if the ancestor is dark, then it adopts the previous level value.
                    rightmost_ime = nodes[rightmost_anc[lv-1]].inside_max_end;
                }

                // for each in nodes on this level
                size_t x = size_t(1)<<(lv-1), step = x<<2;
                for (Rank n = (x<<1)-1; n < nodes.size(); n += step) {
                    // figure inside_max_end
                    Pos ime = nodes[n].end();
                    ime = std::max(ime, nodes[left(n)].inside_max_end);
                    if (right(n) < nodes.size()) {
                        ime = std::max(ime, nodes[right(n)].inside_max_end);
                    } else {
                        ime = std::max(ime, rightmost_ime);
                    }
                    nodes[n].inside_max_end = ime;
                }
            }
        }
    }

    // overlap query; fill ans and return query cost (number of tree nodes visited)
    size_t overlap(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
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

// basic implicit interval tree
//   Pos      : integral position type (e.g. uint32_t, uint64_t)
//   Item     : arbitrary type of items to be indexed
//   get_beg  : function to get interval begin position from Item
//   get_end  : function to get interval end position from Item
// For good performance, beg & end should just access members of Item (i.e. O(1) cache-local fetch)
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
class iit : public iit_base<Pos, Item, iit_node_base<Pos, Item, get_beg, get_end>> {
public:
    template<typename InputIterator>
    iit(InputIterator first, InputIterator last)
        : iit_base<Pos, Item, iit_node_base<Pos, Item, get_beg, get_end>>(first, last)
        {}
};
