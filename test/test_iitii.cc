#include <iostream>
#include <utility>
#include <random>
#include "iitii.h"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
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
    return iitii<pos, pospair, &iitii_beg, &iitii_end>(examples.begin(), examples.end());
}

TEST_CASE("cgranges example") {
    auto tree = build_tree({ { 12, 34 }, { 0, 23 }, { 34, 56 } });

    auto results = tree.overlap(22, 25);
    REQUIRE(results.size() == 2);
    REQUIRE(results[0].first == 0);
    REQUIRE(results[1].first == 12);
}

TEST_CASE("ghost nodes (N=5)") {
    // the three rank levels are:
    // 2.       3
    // 1.     1   5
    // 0.    0 2 4 6
    // with N=5 (ranks 0-4), nodes 5 and 6 are "ghosts"; nb ghost node 5 is the parent of 4
    //
    // lets set up an example where node 4 shall be part of the result set

    vector<pospair> examples = { { 0, 7 }, {1, 2}, {2, 4}, {3, 6}, {4, 9} };
    iitii<pos, pospair, &iitii_beg, &iitii_end> tree(examples.begin(), examples.end());

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

    for (int N = 3; N < 500000; N *= 3) {
        vector<pospair> examples;
        for (int i = 0; i < N; ++i) {
            auto beg = begD(R);
            examples.push_back({ beg, beg+lenD(R) });
        }
        iitii<pos, pospair, &iitii_beg, &iitii_end> tree(examples.begin(), examples.end());
        for (int i = 0; i < 10000; ++i) {
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
            count += results.size();
        }
    }
    cout << "fuzz count = " << count << ", cost = " << cost << endl;
}
