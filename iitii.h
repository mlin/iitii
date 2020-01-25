/*
Implicit Interval Tree with Interpolation Index (iitii)

This file provides two template classes,
    iit: implicit interval tree, a reimplementation of cgranges by Heng Li
    iitii: iit + interpolation index, experimental extension to speed up queries on very large
           datasets

Both classes take the folllowing four template parameters:
    Pos     :  numeric position type (e.g. uint32_t, double)
    Item    :  arbitrary type of the items to be indexed
    get_beg :  a function (const Item& -> Pos) which accesses the interval begin position of item
    get_end :         "              "                   "                 end          "
For good performance, get_beg and get_end should just access members of Item (cache-local fetch).

Example:

    using intpair = std::pair<int,int>;
    int p_get_beg(const intpair& p) { return p.first; }
    int p_get_end(const intpair& p) { return p.second; }
    using p_iit = iit<int, intpair, p_get_beg, p_get_end>; // first arg is position type

    p_iit::builder br;
    br.add(intpair(12,34));
    br.add(intpair(0,23));
    br.add(intpair(34,56));
    p_iit db = br.build();
    // alternative: p_iit db = p_iit::builder(container.begin(), container.end()).build();

    std::vector<const intpair*> results = db.overlap(22, 25);
    // alternative: db.overlap(22, 25, results);

Building iitii works the same way, except build() takes a size_t argument giving the number of
model domains.

This header file has other helper template classes that allow most code to be shared between iit
and iitii (without burdening the former with baggage from the latter -- important for fair
comparative benchmarking). The template structure has gotten a little out of hand, which always
seems to happen.
*/

#include <vector>
#include <limits>
#include <algorithm>
#include <cmath>
#include <assert.h>

// Base template for the internal representation of a node within an implicit interval tree
// User should not care about this; subclass instantiations may add more members for more-
// exotically augmented classes of implicit interval trees
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
struct iit_node_base {
    static const Pos npos = std::numeric_limits<Pos>::max();  // reserved constant for invalid Pos

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

// Base template for an implicit interval tree, with internal repr
//     Node<Pos, Item, ...> : iit_node_base<Pos, Item, ...>
// User should not deal with this directly, but instantiate sub-templates iit or iiitii (below)
template<typename Pos, typename Item, typename Node>
class iit_base {
protected:
    // aliases to help keep the Pos, Rank, and Level concepts straight
    typedef std::size_t Rank;   // rank of a node, its index in the sorted array (or beyond, if
                                //                                                imaginary)
    typedef std::size_t Level;  // level in tree

    static const Rank nrank = std::numeric_limits<Rank>::max();  // invalid Rank

    std::vector<Node> nodes;  // array of Nodes sorted by beginning position
    size_t full_size;         // size of the full binary tree containing the nodes; liable to be
                              // as large as 2*nodes.size()-1, including imaginary nodes.

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
        if (((node>>(lv+1)) & 1)) {  // node is right child
            assert(node >= ofs);
            return node-ofs;
        }
        // node is left child
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
    // recursion depth limited to tree height
    // TODO: nip&tuck optimizations (eliminate recursion, unroll traversal of low levels, etc.)
    size_t scan(Rank subtree, Pos qbeg, Pos qend, std::vector<const Item*>& ans) const {
        if (subtree == nrank) {
            return 0;
        }
        assert(subtree < full_size);
        if (subtree >= nodes.size()) {
            // When we arrive at an imaginary node, its right subtree must be all imaginary, so we
            // only need to descend left.
            // TODO: should be able to work out a closed formula for the next real node we'll find,
            //       given subtree, n.size(), and full_size.
            return 1 + scan(left(subtree), qbeg, qend, ans);
        }

        size_t cost = 1;
        const Node& n = nodes[subtree];
        if (n.inside_max_end > qbeg) {  // something in current subtree extends into/over query
            cost += scan(left(subtree), qbeg, qend, ans);
            Pos nbeg = n.beg();
            if (nbeg < qend) {          // this node isn't already past query
                if (n.end() > qbeg) {   // this node overlaps query
                    ans.push_back(&(n.item));
                }
                cost += scan(right(subtree), qbeg, qend, ans);
            }
        }
        return cost;
    }

public:
    iit_base(std::vector<Node>& nodes_)
        : nodes(std::move(nodes_))
        , root_level(0)
        , root(std::numeric_limits<Rank>::max())
    {
        // compute the implied tree geometry
        for (root_level = 0, full_size = 0; full_size < nodes.size();
             ++root_level, full_size = (size_t(1)<<(root_level+1)) - 1);
        root = (Rank(1) << root_level) - 1;

        if (nodes.size()) {
            // sort the nodes by interval beg (then end)
            std::sort(nodes.begin(), nodes.end());

            // Memoize the path from the rightmost leaf up to the root. This will trace the border
            // between the real and imaginary nodes (if any), which we'll refer to in indexing
            // below. Some of these border nodes may be imaginary.
            std::vector<Rank> right_border_nodes({
                nodes.size() - (2 - nodes.size() % 2)  // rightmost real leaf
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
                        // right child is imaginary; take the last border observation
                        ime = std::max(ime, right_border_ime);
                    }
                    assert(ime != Node::npos);
                    nodes[n].inside_max_end = ime;

                    if (n == right_border_nodes[lv]) {
                        // track inside_max_end of the real nodes on the border
                        right_border_ime = ime;
                    }
                }
            }
        }
    }

    // overlap query; fill ans and return query cost (number of tree nodes visited)
    virtual size_t overlap(Pos qbeg, Pos qend, std::vector<const Item*>& ans) const {
        ans.clear();
        return scan(root, qbeg, qend, ans);
    }

    // overlap query, return vector of results
    std::vector<const Item*> overlap(Pos qbeg, Pos qend) const {
        std::vector<const Item*> ans;
        overlap(qbeg, qend, ans);
        return ans;
    }
};

// template for the builder class exposed by each user-facing class, which takes in items either
// all at once from InputIterator, or streaming one-by-one
template<class iitT, typename Item, typename Node>
class iit_builder_base {
    std::vector<Node> nodes_;

public:
    iit_builder_base() = default;
    template<typename InputIterator>
    iit_builder_base(InputIterator begin, InputIterator end) {
        add(begin, end);
    }

    void add(const Item& it) {
        nodes_.push_back(Node(it));
    }

    template<typename InputIterator>
    void add(InputIterator begin, InputIterator end) {
        std::for_each(begin, end, [&](const Item& it) { add(it); });
    }

    template<typename... Args>
    iitT build(Args&&... args) {
        return iitT(nodes_, std::forward<Args>(args)...);
    }
};

// Basic implicit interval tree (a reimplementation of cgranges)
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
class iit : public iit_base<Pos, Item, iit_node_base<Pos, Item, get_beg, get_end>> {
    using Node = iit_node_base<Pos, Item, get_beg, get_end>;

    iit(std::vector<Node>& nodes_)
        : iit_base<Pos, Item, Node>(nodes_)
        {}

public:
    using builder = iit_builder_base<iit<Pos, Item, get_beg, get_end>, Item, Node>;
    friend builder;
};


// iitii-specialized node type
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
struct iitii_node : public iit_node_base<Pos, Item, get_beg, get_end> {
    // Additional augment value for iitii nodes, which helps us prove when we can stop climbing in
    // the bottom-up search for a subtree root which must contain all query results beneath it.
    Pos outside_max_end;
    // outside_max_end of node n is the maximum m.end() of all nodes m outside of n & its subtree
    //     with m.beg() < n.beg(); -infinity if there are no such nodes.
    //
    // Furthermore,
    // outside_min_beg of node n is the minimum m.beg() of all nodes m outside of n & its subtree
    //     with m.beg() >= n.beg(); infinity if there are no such nodes.
    // But we don't need to store outside_min_beg, because we can compute it in constant time using
    // rank offsets in the beg-sorted node array.
    //
    // Suppose during a query for [qbeg, qend) we climb up to a node n with,
    //   (i) n.outside_max_end <= qbeg; AND
    //  (ii) qend <= n.outside_min_beg
    // By (i), any node outside n & subtree with beg < n's cannot overlap the query. By (ii), any
    // node outside n & subtree with beg >= n's cannot overlap the query. This exhausts all nodes
    // outside of n & subtree, so we need not climb past n.

    iitii_node(const Item& item_)
        : iit_node_base<Pos, Item, get_beg, get_end>(item_)
        , outside_max_end(std::numeric_limits<Pos>::min())
        {}
};

// simple linear regression of y ~ x given points [(x,y)], returning (intercept, slope)
template<typename XT, typename YT>
std::pair<double,double> regress(const std::vector<std::pair<XT,YT>>& points) {
    if (points.size() < 1) {
        return std::make_pair(0.0, 0.0);
    }
    double sum_x, sum_y, cov, var;
    sum_x = sum_y = cov = var = 0.0;
    for (const auto& pt : points) {
        sum_x += double(pt.first);
        sum_y += double(pt.second);
    }
    const double mean_x = sum_x/points.size(), mean_y = sum_y/points.size();
    for (const auto& pt : points) {
        const double x_err = pt.first - mean_x;
        cov += x_err*(pt.second - mean_y);
        var += x_err*x_err;
    }
    if (var == 0.0) {
        return std::make_pair(0.0, 0.0);
    }
    const double m = cov / var;
    return std::make_pair(mean_y - m*mean_x, m);
}

// here it is
template<typename Pos, typename Item, Pos get_beg(const Item&), Pos get_end(const Item&)>
class iitii : public iit_base<Pos, Item, iitii_node<Pos, Item, get_beg, get_end>> {
    using Node = iitii_node<Pos, Item, get_beg, get_end>;
    using super = iit_base<Pos, Item, Node>;
    using typename super::Rank;
    using typename super::Level;
    using super::left;
    using super::right;
    using super::parent;
    using super::nodes;
    using super::full_size;
    using super::nrank;
    using super::level;
    using super::root;
    using super::root_level;

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

    inline Pos outside_min_beg(Rank subtree) const {
        // constant-time computation of outside_min_beg: beg() of the node ranked one higher than
        // subtree's rightmost child
        const Rank r = rightmost_child(subtree);
        __builtin_prefetch(&(nodes[r+1]));
        const Pos beg = nodes[subtree].beg();
        const Rank l = leftmost_child(subtree);
        if (l && nodes[l-1].beg() == beg) {
            // corner case: nodes to the left of the subtree can have the same beg as subroot
            // and outside_min_beg is defined on nodes with beg >= subroot's.
            return beg;
        }
        return r < nodes.size()-1 ? nodes[r+1].beg() : std::numeric_limits<Pos>::max();
    }

    // Additional tree navigation concept, LevelRank: the rank of a node **within its level**
    // e.g. a level-k node with LevelRank=1 is the second-lowest (second-leftmost) node on level k
    typedef std::size_t LevelRank;

    inline Rank rank_of_levelrank(Level lv, LevelRank ofs) const {
        return (size_t(1)<<lv)*(2*ofs+1)-1;
    }

    inline LevelRank levelrank_of_rank(Rank r) const {
        return ((r+1)/(size_t(1)<<level(r))-1)/2;
    }

    // Rank prediction model: the [min_beg, max_beg] range is partitioned into a number C of
    // domains, each covering an equal-sized portion of that range. The domain pertaining to a
    // position beg is d(beg) = floor((beg-min_beg)*C/(max_beg-min_beg)), bounded to [0,C).
    //
    // For each domain d, we store three parameters: a Level l[d] ∈ [0,root_level] into which we
    // will jump, and linear weights w[d,0] and w[d,1] for the regression of LevelRank on Pos,
    //   lr(beg) ~ w[d(beg),0] + w[d(beg),1]*beg
    //
    // To start a query for qbeg, jump to the node: rank_of_levelrank(l[d(qbeg)], lr(qbeg))

    typedef std::size_t Domain;
    Domain domains;               // C
    Pos min_beg = std::numeric_limits<Pos>::max(),
        domain_size = Node::npos;;
    std::vector<float> parameters;  // C rows of three parameters (row-major storage): w[0,d],
                                    // w[1,d] and l[d]. NB: the third is a Level stored as a float.

    inline Domain which_domain(Pos beg) const {
        if (beg < min_beg) {
            return 0;
        }
        return std::min(domains-1, Domain((beg-min_beg)/domain_size));
    }

    inline Rank interpolate(Level lv, float w0, float w1, Pos qbeg) const {
        // given model parameters within a domain, return the node to start searching for qbeg
        const float ofs_f = w0 + w1*float(qbeg);
        assert(std::isfinite(ofs_f));
        const Rank r = rank_of_levelrank(lv, LevelRank(std::max(0.0f, roundf(ofs_f))));
        assert(r >= nodes.size() || level(r) == lv);

        // detail: if rank is imaginary (qbeg is off-scale high), start from rightmost real leaf
        const auto nsz = nodes.size();
        return r < nsz ? r : (nsz - (2 - nsz%2));
    }

    void train() {
        // scan the nodes to extract <Pos,Rank> points partitioned by domain
        std::vector<std::vector<std::pair<Pos,Rank>>> points(domains);
        for (Rank r = 0; r < nodes.size(); ++r) {
            points.at(which_domain(nodes[r].beg())).push_back(std::make_pair(nodes[r].beg(), r));
        }
        // train each domain-specific model
        for (Domain domain = 0; domain < domains; ++domain) {
            // partition the domain points by tree level, converting Ranks to LevelRanks
            std::vector<std::vector<std::pair<Pos,LevelRank>>> points_by_level(root_level+1);
            for (const auto& p : points[domain]) {
                Level lv = level(p.second);
                points_by_level.at(lv)
                               .push_back(std::make_pair(p.first, levelrank_of_rank(p.second)));
            }
            // for each level,
            double lowest_cost = std::numeric_limits<double>::max();
            for (Level lv = 0; lv <= 3*root_level/4; ++lv) {
                // regress points on this level
                auto w = regress<Pos,LevelRank>(points_by_level[lv]);
                if (w.second) {
                    // calculate estimate of search cost (average over all domain points)
                    double cost = 0.0;
                    for (const auto& p : points[domain]) {
                        Rank fb = interpolate(lv, float(w.first), float(w.second), p.first);
                        double error = fabs(double(fb) - double(p.second))/double(size_t(1)<<lv);
                        double error_penalty = error > 1.0 ? 2.0*log2(1+error) : 0.0;
                        double overlap_penalty = nodes[fb].outside_max_end > p.first ? 0.5*(root_level-lv) : 0;
                        cost += lv + std::max(error_penalty, overlap_penalty);
                    }
                    cost /= points[domain].size();
                    assert(std::isfinite(cost) && cost >= 0.0);
                    // store parameters if cost estimate is lower than top-down search and lower
                    // than previous levels
                    if (cost < root_level && cost < lowest_cost) {
                        lowest_cost = cost;
                        float *pp = &(parameters[3*domain]);
                        pp[0] = float(w.first);
                        pp[1] = float(w.second);
                        pp[2] = float(lv);
                    }
                }
            }
            points[domain].clear(); // free a little memory
            /*
            std::cout << "domain = " << domain << " level = " << Level(parameters[3*domain+2])
                      << " E[cost] = " << lowest_cost << std::endl;
            */
        }
    }

    // Given qbeg, select domain and predict search start node
    Rank predict(Pos qbeg) const {
        auto which = which_domain(qbeg);
        assert(which < domains);
        const float *pp = &(parameters[3*which]);

        const float lv_f = pp[2];
        if (lv_f < 0) {
            return nrank;
        }
        assert(lv_f >= 0 && lv_f <= root_level);
        const Level lv = (Level) lv_f;

        return interpolate(lv, pp[0], pp[1], qbeg);
    }

    iitii(std::vector<Node>& nodes_, Domain domains_)
        : super(nodes_)
        , domains(std::max(Domain(1),domains_))
        , domain_size(std::numeric_limits<Pos>::max())
    {
        parameters.resize(domains*3, -1.0f);

        if (nodes.size()) {
            // equal size (in Pos units) of each domain
            min_beg = nodes[0].beg();
            domain_size = 1 + (nodes.back().beg()-min_beg)/domains;

            // compute running max_end along the sorted array, which we'll look up while computing
            // outside_max_end below
            std::vector<Pos> running_max_end { nodes[0].end() };
            for (Rank n = 1; n < nodes.size(); ++n) {
                running_max_end.push_back(std::max(running_max_end[n-1], nodes[n].end()));
            }

            // fill outside_max_end
            for (Rank n = 0; n < nodes.size(); ++n) {
                Node& node = nodes[n];
                Rank l = leftmost_child(n);

                if (l>0) {
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
            }

            // train the rank prediction models
            train();
        }
    }

public:
    // iitii::builder::build() takes a size_t argument giving the number of domains to model
    using builder = iit_builder_base<iitii<Pos, Item, get_beg, get_end>, Item, Node>;
    friend builder;

    size_t overlap(Pos qbeg, Pos qend, std::vector<const Item*>& ans) const override {
        // ask model which leaf we should begin our bottom-up climb at
        Rank prediction = predict(qbeg);
        if (prediction == nrank) {
            // the model did not make a prediction for some reason, so just go to the root
            return super::overlap(qbeg, qend, ans);
        }
        assert(level(prediction) <= root_level);
        __builtin_prefetch(&(nodes[parent(prediction)]));

        // climb until our necessary & sufficient criteria are met, or the root
        size_t climb_cost = 0;
        Rank subtree = prediction;
        while (subtree != root &&                           // stop at root
                (subtree >= nodes.size() ||                 // continue climb through imaginary
                 qbeg < nodes[subtree].outside_max_end ||   // possible outside overlap from left
                 outside_min_beg(subtree) < qend)) {        // possible outside overlap from right
            subtree = parent(subtree);
            __builtin_prefetch(&(nodes[parent(subtree)]));
            // TODO: scheme to skip through chains of imaginary nodes along the border
            ++climb_cost;
        }

        auto self = const_cast<iitii<Pos, Item, get_beg, get_end>*>(this);  // getting around const
        self->queries++;
        self->total_climb_cost += climb_cost;

        // scan the subtree for query results.
        // pessimistically, we triple the climbing cost when adding it to the top-down search cost,
        // because the outside_min_beg() lookup may incur two additional cache misses.
        ans.clear();
        return super::scan(subtree, qbeg, qend, ans) + 3*climb_cost;
    }

    size_t queries = 0;
    size_t total_climb_cost = 0;

    using super::overlap;
};
