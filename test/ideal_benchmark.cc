// benchmark iitii with synthetic ideal data (beg positions can be interpolated perfectly)

#include "util.h"
#include <fstream>
#include <random>

struct ideal_item {
    uint32_t beg;
    uint32_t end;
};

uint32_t ideal_beg(const ideal_item& it) { return it.beg; }
uint32_t ideal_end(const ideal_item& it) { return it.end; }

using ideal_iit = iit<uint32_t, ideal_item, ideal_beg, ideal_end>;
using ideal_iitii = iitii<uint32_t, ideal_item, ideal_beg, ideal_end>;

vector<ideal_item> generate(size_t N) {
    // each item ranked i has begin position 10*i, with geometrically distributed length mean 20.
    default_random_engine R(42);
    geometric_distribution<uint32_t> lenD(0.05);
    vector<ideal_item> ans;

    for (uint32_t i = 0; i < N; i++) {
        ideal_item it;
        it.beg = i*10;
        it.end = it.beg + lenD(R);
        ans.push_back(it);
    }

    return ans;
}

template <class tree>
size_t run_queries(const tree& t, uint32_t max_end, size_t queries, size_t& cost) {
    default_random_engine R(42);
    uniform_int_distribution<uint32_t> begD(0, max_end);
    geometric_distribution<uint32_t> lenD(0.025);
    size_t ans = 0;
    cost = 0;
    for (size_t i = 0; i < queries; i++) {
        auto qbeg = begD(R);
        auto qend = qbeg+lenD(R);
        vector<const ideal_item*> results;
        cost += t.overlap(qbeg, qend, results);
        ans += results.size();
    }
    return ans;
}

template <class tree, typename... Args>
size_t run_experiment(vector<ideal_item> items, uint64_t& build_ms, uint64_t& queries_ms, size_t& cost, Args&&... args) {
    unique_ptr<tree> ptree;
    build_ms = milliseconds_to([&](){
        auto t = typename tree::builder(items.begin(), items.end()).build(forward<Args>(args)...);
        ptree.reset(new tree(move(t)));
    });

    cost = 0;
    size_t result_count = 0;
    queries_ms = milliseconds_to([&](){
        result_count = run_queries<tree>(*ptree, items.size()*10, 10000000, cost);
    });

    return result_count;
}

int main(int argc, char** argv) {
    cout << "#tree_type\tN\tbuild_ms\tqueries_ms\tqueries_cost\tresult_count" << endl;
    for (size_t s = 20; s <= 28; s+=2) {
        size_t N = 1 << s;
        auto items = generate(N);
        uint64_t build_ms, queries_ms;
        size_t cost;
        size_t result_count = run_experiment<ideal_iit>(items, build_ms, queries_ms, cost);
        cout << "iit\t" << N << "\t" << build_ms << "\t" << queries_ms << "\t" << cost << "\t" << result_count << endl;
        if (result_count != run_experiment<ideal_iitii>(items, build_ms, queries_ms, cost, 1)) {
            throw runtime_error("RED ALERT: inconsistent results");
        }
        cout << "iitii\t" << N << "\t" << build_ms << "\t" << queries_ms << "\t" << cost << "\t" << result_count << endl;
    }

    return 0;
}
