#include <vector>
#include <limits>
#include <algorithm>
#include <assert.h>

template<typename Pos, typename Item, Pos beg(const Item&), Pos end(const Item&)>
class iitii {
    static const Pos npos = std::numeric_limits<Pos>::max();
    struct Node {
        Item item;
        Pos interior_max_end;  // max end of this & subtree (as in standard augmented interval tree)
        Pos exterior_max_end;  // max end of all nodes n with A[n].beg <= this->beg, excluding this & subtree
        Pos exterior_min_beg;  // min beg of all nodes n with A[n].beg >= this->beg, excluding this & subtree

        Node(const Item& item_)
            : item(item_)
            , interior_max_end(std::numeric_limits<Pos>::min())
            , exterior_max_end(std::numeric_limits<Pos>::min())
            , exterior_min_beg(npos)
            {}
    };
    std::vector<Node> A;
    typedef std::size_t Rank;
    typedef std::size_t Level;
    size_t nsz = std::numeric_limits<size_t>::max();
    Rank root;
    Level root_level;          // K
    Pos max_end;
    Pos min_beg, max_beg;      // for interpolation

    inline Level level(Rank node) const {
        Level ans;
        for (ans=0; node&1; ans++, node>>=1);
        return ans;
    }

    inline Rank parent(Rank node, Level level) const {
        assert(node != root);
        Rank ofs = Rank(1) << level;
        assert(node >= ofs-1);
        if (((node>>(level+1)) & 1)) {
            assert(node >= ofs);
            return node-ofs;
        }
        return node+ofs;
    }

    inline Rank left(Rank node, Level lv) const {
        assert(lv && lv == level(node));
        return node - (Rank(1) << (lv-1));
    }

    inline Rank left(Rank node) const {
        return left(node, level(node));
    }

    inline Rank right(Rank node, Level lv) const {
        assert(lv && lv == level(node));
        return node + (Rank(1) << (lv-1));
    }

    inline Rank right(Rank node) const {
        return right(node, level(node));
    }

    inline Pos effective_interior_max_end(Rank node) const {
        if (node >= A.size()) {
            return max_end;
        }
        return A[node].interior_max_end;
    }

    inline Pos effective_exterior_max_end(Rank node) const {
        if (node >= A.size()) {
            // FIXME what to put here?
        }
        return A[node].exterior_max_end;
    }

    inline Pos effective_exterior_min_beg(Rank node) const {
        if (node >= A.size()) {
            // FIXME what to put here?
        }
        return A[node].exterior_min_beg;
    }

    void scan(Rank subtree, Level lv, Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        assert(lv == level(subtree));
        if (effective_interior_max_end(subtree) > qbeg) {
            if (lv) {
                scan(left(subtree, lv), lv-1, qbeg, qend, ans);
            }
            if (subtree < A.size()) {
                const Node& n = A[subtree];
                Pos begn = beg(n.item);
                if (begn < qend && qbeg < end(n.item)) {
                    ans.push_back(n.item);
                }
                if (lv && begn < qend) {
                    scan(right(subtree, lv), lv-1, qbeg, qend, ans);
                }
            } else if (lv) {
                // TODO: ghost node stopping condition
                scan(right(subtree, lv), lv-1, qbeg, qend, ans);
            }
        }
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
        , max_end(std::numeric_limits<Pos>::min())
    {
        std::for_each(first, last, [&](const Item& it) {
            Pos beg_it = beg(it), end_it = end(it);
            min_beg = std::min(min_beg, beg_it);
            max_beg = std::max(max_beg, beg_it);
            max_end = std::max(max_end, end_it);
            A.push_back(Node(it));
            (A.end()-1)->interior_max_end = end_it;
        });
        std::sort(A.begin(), A.end(), [](const Node& lhs, const Node& rhs) { return beg(lhs.item) < beg(rhs.item); });
        if (A.size()) {
            for (root_level = 0; (Rank(1)<<(root_level+1))-1 < A.size(); root_level++);
            root = (Rank(1) << root_level) - 1;
        }
        for (Level lv=1; lv <= root_level; ++lv) {
            size_t x = size_t(1)<<(lv-1), step = x<<2;
            for (Rank n = (x<<1)-1; n < A.size(); n += step) {
                A[n].interior_max_end = std::max(A[n].interior_max_end, effective_interior_max_end(left(n, lv)));
                A[n].interior_max_end = std::max(A[n].interior_max_end, effective_interior_max_end(right(n, lv)));
            }
        }
        // FIXME populate exterior_max_end & exterior_min_beg
    }

    void overlap(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        ans.clear();
        scan(root, root_level, qbeg, qend, ans);
    }

    void overlap_interpolate(Pos qbeg, Pos qend, std::vector<Item>& ans) const {
        Rank subtree = interpolate(qbeg);
        Level level = 0;
        while (level < root_level && !(effective_exterior_max_end(subtree) <= qbeg &&
                                       qend <= effective_exterior_min_beg(subtree))) {
            subtree = parent(subtree, level++);
        }

        ans.clear();
        scan(qbeg, qend, ans);
    }
};
