#include <iostream>
#include <utility>
#include "iitii.h"
#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch2/catch.hpp"

using namespace std;

unsigned int Factorial( unsigned int number ) {
    return number <= 1 ? number : Factorial(number-1)*number;
}

uint32_t iitii_beg(const pair<uint32_t, uint32_t>& p) {
    return p.first;
}

uint32_t iitii_end(const pair<uint32_t, uint32_t>& p) {
    return p.second;
}

TEST_CASE("cgranges example") {
    vector<pair<uint32_t, uint32_t>> examples = { { 12, 34 }, { 0, 23 }, { 34, 56 } };
    iitii<uint32_t, pair<uint32_t, uint32_t>, &iitii_beg, &iitii_end> tree(examples.begin(), examples.end());

    examples.clear();
    tree.overlap(22, 25, examples);
    REQUIRE(examples.size() == 2);
}
