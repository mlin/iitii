#include <vector>
#include <limits>
#include <algorithm>
#include <assert.h>

// class iitii template parameters:
//   Pos  : integral position type (e.g. uint32_t, uint64_t)
//   Item : arbitrary type of items to be indexed
//   beg  : function to get interval begin position from Item
//   end  : function to get interval end position from Item
// For good performance, beg & end should just access members of Item (i.e. O(1) cache-local fetch)
template<typename Pos, typename Item, Pos beg(const Item&), Pos end(const Item&)>
class iitii {
    static const Pos npos = std::numeric_limits<Pos>::max();  // reserved constant for invalid Pos

    struct Node {  // augmented item
        Item item;
        Pos inside_max_end;   // max end of this & subtree (as in textbook augmented interval tree)
        Pos outside_max_end;  // max end of all nodes n with A[n].beg <= this->beg, excluding this & subtree
        Pos outside_min_beg;  // min beg of all nodes n with A[n].beg >= this->beg, excluding this & subtree

        Node(const Item& item_)
            : item(item_)
            , inside_max_end(end(item_))
            , outside_max_end(std::numeric_limits<Pos>::min())
            , outside_min_beg(npos)
            {}
    };
    std::vector<Node> A;       // sorted array of Nodes

    // beg & end extreme values in A
    Pos min_beg, max_beg;

    // type aliases Rank and Level, just to help keep these concepts straight (& Pos)
    typedef std::size_t Rank;  // rank of a node, or its index in A
    Rank root;
    static const Rank nrank = std::numeric_limits<Rank>::max();  // reserved constant for invalid Rank
    typedef std::size_t Level;
    Level root_level;          //  = K in cgranges

    // compute a node's level, the # of 1 bits below the lowest 0 bit
    inline Level level(Rank node) const {
        assert(node < (Rank(1)<<(root_level+1))-1);
        Level ans = __builtin_ctzll(~node);
        #ifndef NDEBUG
        Level chk;
        for (chk=0; node&1; chk++, node>>=1);
        assert(ans == chk);
        #endif
        return ans;
    }

    inline Rank parent(Rank node) const {
        assert(node < ((Rank(1)<<(root_level+1))-1));
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

    inline Rank left(Rank node) const {
        Level lv = level(node);
        return lv > 0 ? node - (Rank(1) << (lv-1)) : nrank;
    }

    inline Rank right(Rank node) const {
        Level lv = level(node);
        return lv > 0 ? node + (Rank(1) << (lv-1)) : nrank;
    }

    size_t scan(Rank subtree, Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        if (subtree == nrank) {
            return 0;
        }
        if (subtree >= A.size()) {
            return 1 + scan(left(subtree), qbeg, qend, ans);
        }
        size_t cost = 1;
        const Node& n = A[subtree];
        if (n.inside_max_end > qbeg) {
            cost += scan(left(subtree), qbeg, qend, ans);
            Pos nbeg = beg(n.item);
            if (nbeg < qend) {
                if (qbeg < end(n.item)) {
                    ans.push_back(n.item);
                }
                cost += scan(right(subtree), qbeg, qend, ans);
            }
        }
        return cost;
    }

    Rank interpolate(Pos q) const {
        if (q <= min_beg) {
            return 0;
        }
        if (q >= max_beg) {
            return A.size()-1;
        }
        return 2*Rank(((q-min_beg)/double(max_beg-min_beg))*(A.size()/2));
    }

public:
    template<typename InputIterator>
    iitii(InputIterator first, InputIterator last)
        : root_level(0)
        , root(std::numeric_limits<Rank>::max())
        , min_beg(npos)
        , max_beg(std::numeric_limits<Pos>::min())
    {
        std::for_each(first, last, [&](const Item& it) {
            Pos beg_it = beg(it);
            min_beg = std::min(min_beg, beg_it);
            max_beg = std::max(max_beg, beg_it);
            A.push_back(Node(it));
        });
        if (A.size()) {
            std::sort(A.begin(), A.end(),
                  [](const Node& lhs, const Node& rhs) { return beg(lhs.item) < beg(rhs.item); });

            for (root_level = 0; (Rank(1)<<(root_level+1))-1 < A.size(); root_level++);
            root = (Rank(1) << root_level) - 1;
            
            Rank rightmost_leaf = A.size() - (2 - A.size() % 2);
            Rank rightmost_anc = rightmost_leaf;
            Pos rightmost_ime = A[rightmost_anc].inside_max_end;
            for (Level lv=1; lv <= root_level; ++lv) {
                size_t x = size_t(1)<<(lv-1), step = x<<2;
                for (Rank n = (x<<1)-1; n < A.size(); n += step) {
                    Pos ime = A[n].inside_max_end;
                    ime = std::max(ime, A[left(n)].inside_max_end);
                    if (right(n) < A.size()) {
                        ime = std::max(ime, A[right(n)].inside_max_end);
                    } else {
                        ime = std::max(ime, rightmost_ime);
                    }
                    A[n].inside_max_end = ime;
                }
                rightmost_anc = parent(rightmost_anc);
                if (rightmost_anc < A.size()) {
                    rightmost_ime = A[rightmost_anc].inside_max_end;
                }
            }
        }
        // FIXME populate outside_max_end & outside_min_beg
    }

    size_t overlap(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        ans.clear();
        return scan(root, qbeg, qend, ans);
    }

    std::vector<Item> overlap(Pos qbeg, Pos qend) const {
        std::vector<Item> ans;
        scan(root, qbeg, qend, ans);
        return ans;
    }
    /*
    size_t overlap_interpolate(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        Rank subtree = interpolate(qbeg);
        Level level = 0;
        while (level < root_level && !(effective_outside_max_end(subtree) <= qbeg &&
                                       qend <= effective_outside_min_beg(subtree))) {
            subtree = parent(subtree, level++);
        }

        ans.clear();
        return scan(qbeg, qend, ans);
    } */
};
