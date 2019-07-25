#include <iostream>
#include <utility>
#include <random>
#include "iitii.h"
#define CATCH_CONFIG_MAIN
#include "catch2/catch.hpp"

using namespace std;

typedef uint32_t pos;
typedef pair<pos, pos> pospair;

pos iitii_beg(const pospair& p) {
    return p.first;
}

pos iitii_end(const pospair& p) {
    return p.second;
}

auto build_tree(const vector<pospair>& examples) {
    return iit<pos, pospair, &iitii_beg, &iitii_end>(examples.begin(), examples.end());
}

TEST_CASE("cgranges example") {
    auto tree = build_tree({ { 12, 34 }, { 0, 23 }, { 34, 56 } });

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
    auto tree = build_tree(examples);

    auto results = tree.overlap(6, 10);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 4);
}

TEST_CASE("fuzz") {
    default_random_engine R(42);
    uniform_int_distribution<uint32_t> begD(1, 42000);
    geometric_distribution<uint16_t> lenD(0.01);
    size_t count = 0, cost = 0;

    for (int N = 3; N < 1000000; N *= 3) {  // base != 2 provides varying tree fullness patterns
        vector<pospair> examples;
        for (int i = 0; i < N; ++i) {
            auto beg = begD(R);
            examples.push_back({ beg, beg+lenD(R) });
        }

        auto tree = build_tree(examples);

        std::sort(examples.begin(), examples.end(), [](const pospair& lhs, const pospair& rhs) {
            auto begl = iitii_beg(lhs), begr = iitii_beg(rhs);
            if (begl == begr) {
                return iitii_end(lhs) < iitii_end(rhs);
            }
            return begl < begr;
        });

        for (int i = 0; i < 1000; ++i) {
            auto qbeg = begD(R);
            auto qend = qbeg + 100;
            vector<pospair> results;
            cost += tree.overlap(qbeg, qend, results);

            vector<pospair> results2;
            for (const auto& p : examples) {
                if (qbeg < p.second && p.first < qend) {
                    results2.push_back(p);
                }
            }

            REQUIRE(results.size() == results2.size());
            bool alleq = true;
            for (auto p1 = results.begin(), p2 = results2.begin(); p1 != results.end(); ++p1, ++p2) {
                alleq = alleq && (*p1 == *p2);
            }
            REQUIRE(alleq);
            count += results.size();
        }
    }
    cout << "fuzz count = " << count << ", cost = " << cost << endl;
}
