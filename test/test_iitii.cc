#include <iostream>
#include <utility>
#include <random>
#include <math.h>
#include "iitii.h"
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

using namespace std;

typedef uint32_t pos;
typedef pair<pos, pos> pospair;

pos get_beg(const pospair& p) {
    return p.first;
}

pos get_end(const pospair& p) {
    return p.second;
}

auto build_iit(const vector<pospair>& examples) {
    return iit<pos, pospair, &get_beg, &get_end>(examples.begin(), examples.end());
}

auto build_iitii(const vector<pospair>& examples) {
    return iitii<pos, pospair, &get_beg, &get_end>(examples.begin(), examples.end());
}

TEST_CASE("cgranges example") {
    auto tree = build_iit({ { 12, 34 }, { 0, 23 }, { 34, 56 } });

    auto results = tree.overlap(22, 25);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 12);
}

TEST_CASE("cgranges example with iitii") {
    auto tree = build_iitii({ { 12, 34 }, { 0, 23 }, { 34, 56 } });

    auto results = tree.overlap(22, 25);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 12);
}

TEST_CASE("dark nodes (N=5)") {
    // the three rank levels are:
    // 2.       3
    // 1.     1   5
    // 0.    0 2 4 6
    // with N=5 (ranks 0-4), nodes 5 and 6 are "dark"; dark node 5 is the parent of 4
    //
    // lets set up an example where node 4 shall be part of the result set

    vector<pospair> examples = { { 0, 7 }, {1, 2}, {2, 4}, {3, 6}, {4, 9} };
    auto tree = build_iit(examples);

    auto results = tree.overlap(6, 10);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 4);

    // to answer the following query, interval tree algo should visit nodes 1, 3, 4, and 5
    REQUIRE(tree.overlap(7, 10, results) == 4);
}

TEST_CASE("dark nodes (N=5) with iitii") {
    vector<pospair> examples = { { 0, 7 }, {1, 2}, {2, 4}, {3, 6}, {4, 9} };
    auto tree = build_iitii(examples);

    auto results = tree.overlap(6, 10);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 4);

    // with the interpolation index, we can answer this query in one step
    REQUIRE(tree.overlap(7, 10, results) == 1);
}

TEST_CASE("fuzz") {
    default_random_engine R(42);

    uniform_int_distribution<uint32_t> begD(1, 420000);
    geometric_distribution<uint16_t> lenD(0.01);
    const vector<uint32_t> spike_positions = {100, 1000, 10000, 100000, 420000};

    for (int N = 3; N < 2000000; N *= 3) {  // base != 2 provides varying tree fullness patterns
        // generate random intervals with beg ~ begD and length ~ lenD
        vector<pospair> examples;
        for (int i = 0; i < N; ++i) {
            auto beg = begD(R);
            examples.push_back({ beg, beg+lenD(R) });
        }
        // also spike in a bunch of intervals starting at a few selected positions, since colliding
        // beg positions trigger certain corner cases
        for (int i = 0; i < N/10; ++i) {
            auto beg = spike_positions[i % spike_positions.size()];
            examples.push_back({ beg, beg+lenD(R) });
        }

        // build trees
        auto tree = build_iit(examples);
        auto treeii = build_iitii(examples);

        std::sort(examples.begin(), examples.end(), [](const pospair& lhs, const pospair& rhs) {
            auto begl = get_beg(lhs), begr = get_beg(rhs);
            if (begl == begr) {
                return get_end(lhs) < get_end(rhs);
            }
            return begl < begr;
        });

        // run random queries and check that the result sets are correct
        const size_t Q = 1000;
        size_t results = 0, cost = 0, costii = 0;
        for (size_t i = 0; i < Q; ++i) {
            auto qbeg = begD(R);
            auto qend = qbeg + 42;
            vector<pospair> ans;
            cost += tree.overlap(qbeg, qend, ans);

            vector<pospair> naive;
            for (const auto& p : examples) {
                if (qbeg < p.second && p.first < qend) {
                    naive.push_back(p);
                }
            }

            REQUIRE(ans.size() == naive.size());
            bool alleq = true;
            for (auto p1 = ans.begin(), p2 = naive.begin(); p1 != ans.end(); ++p1, ++p2) {
                alleq = alleq && (*p1 == *p2);
            }
            REQUIRE(alleq);
            results += ans.size();

            costii += treeii.overlap(qbeg, qend, ans);
            REQUIRE(ans.size() == naive.size());
            alleq = true;
            for (auto p1 = ans.begin(), p2 = naive.begin(); p1 != ans.end(); ++p1, ++p2) {
                alleq = alleq && (*p1 == *p2);
            }
            REQUIRE(alleq);
        }

        cout << "fuzz N = " << N << ": results = " << results << ", cost = " << cost << ", costii = " << costii << endl;

        // guess bound on cost per query: 2*(lg(N) + results)
        size_t cost_bound = size_t(2*Q*(log2(N) + results/Q + 1));
        REQUIRE(cost < cost_bound);

        // iitii should show a cost advantage unless query cost is dominated by result set size
        if (log2(N) > results/Q + 2) {
            REQUIRE(costii < cost);
        }
    }
}
