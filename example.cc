// g++ -std=c++11 -o example example.cc
#include <iostream>
#include <vector>
#include "iitii.h"

using intpair = std::pair<int,int>;
int p_get_beg(const intpair& p) { return p.first; }
int p_get_end(const intpair& p) { return p.second; }
using p_iitii = iitii<int, intpair, p_get_beg, p_get_end>;  // first arg is position type

int main() {
    p_iitii::builder br;
    br.add(intpair(12,34));
    br.add(intpair(0,23));
    br.add(intpair(34,56));
    p_iitii db = br.build(1);  // 1 = # domains
    // alternative: p_iitii db = p_iitii::builder(container.begin(), container.end()).build(1);

    std::vector<intpair> results = db.overlap(22, 25);
    // alternative: db.overlap(22, 25, results);

    for (const auto& p : results)
        std::cout << p.first << "\t" << p.second << std::endl;
    return 0;
}
