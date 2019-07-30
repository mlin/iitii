#include "util.h"
#include <fstream>
#include <random>
#include <chrono>
#include <functional>

uint32_t milliseconds_to(function<void()> f) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    f();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    return chrono::duration_cast<chrono::milliseconds>(end - begin).count();
}

template <class tree>
size_t run_queries(const vector<variant>& variants, const tree& t, int max_end, int queries, size_t& cost) {
    default_random_engine R(42);
    uniform_int_distribution<uint32_t> begD(0, max_end);
    uniform_int_distribution<size_t> vtD(0, variants.size()-1);
    size_t ans = 0;
    cost = 0;
    for (int i = 0; i < queries; i++) {
        // 50% queries for the interval of a random existing variant (results will include itself)
        // and 50% for 10bp intervals with a uniform random begin position
        auto qbeg = begD(R);
        auto qend = qbeg+10;
        if (i % 2 == 1) {
            const auto& vt = variants.at(vtD(R));
            qbeg = vt.beg;
            qend = vt.end;
        }
        vector<variant> results;
        cost += t.overlap(qbeg, qend, results);
        ans += results.size();
    }
    return ans;
}

template <class tree, typename... Args>
size_t run_experiment(const vector<variant>& variants, const size_t N,
                      uint32_t& build_ms, uint32_t& queries_ms, size_t& cost, Args&&... args) {
    vector<variant> variantsN(variants.begin(), variants.begin()+N);
    int max_end = -1;
    for (const auto& vt : variantsN) {
        max_end = std::max(max_end, vt.end);
    }
    unique_ptr<tree> ptree;
    build_ms = milliseconds_to([&](){
        auto t = typename tree::builder(variantsN.begin(), variantsN.end()).build(forward<Args>(args)...);
        ptree.reset(new tree(move(t)));
    });

    cost = 0;
    size_t result_count = 0;
    queries_ms = milliseconds_to([&](){
        result_count = run_queries<tree>(variantsN, *ptree, max_end, 10000000, cost);
    });
    // cout << "mean climbing per iitii query: " << double(ptree->total_climb_cost)/ptree->queries << endl;

    return result_count;
}

int main(int argc, char** argv) {
    // As of this writing (2019-07-29) newer gnomAD versions have far larger files but not many
    // additional variants (a lot more metadata)
    const string filename = "/tmp/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz";
    const string url = "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz";
    #ifdef NDEBUG
    const int megabases = 244;
    #else
    const int megabases = 24;
    #endif

    ifstream vcf(filename), tbi(filename + ".tbi");
    if (!(vcf.good() && tbi.good())) {
        cerr << "This program requires " << filename
             << " and .tbi to be present. Download them to that location from " << url << endl;
        return 1;
    }

    auto variants = load_variants_parallel(filename, 0, megabases);
    int max_len = -1, max_end = -1;
    for (const auto& vt : variants) {
        max_len = max(max_len, vt.end - vt.beg);
        max_end = max(max_end, vt.end);
    }
    cerr << variants.size() << " variants, max END = " << max_end
         << ", max rlen = " << max_len << endl;

    std::sort(variants.begin(), variants.end(), [](const variant& lhs, const variant& rhs) {
        auto begl = lhs.beg, begr = rhs.beg;
        if (begl == begr) {
            return lhs.end < rhs.end;
        }
        return begl < begr;
    });

    cout << "#tree_type\tnum_variants\tbuild_ms\tqueries_ms\tqueries_cost\tresult_count" << endl;
    for(size_t N = variants.size(); N >= 10000; N /= 4) {
        uint32_t build_ms, queries_ms;
        size_t cost;
        size_t result_count = run_experiment<variant_iit>(variants, N, build_ms, queries_ms, cost);
        cout << "iit\t" << N << "\t" << build_ms << "\t" << queries_ms << "\t" << cost << "\t" << result_count << endl;
        for (unsigned domains = 1; domains <= 4096; domains *= 16) {
            if (result_count != run_experiment<variant_iitii>(variants, N, build_ms, queries_ms, cost, domains)) {
                throw runtime_error("RED ALERT: inconsistent results");
            }
            cout << "iitii(" << domains << ")\t" << N << "\t" << build_ms << "\t" << queries_ms << "\t" << cost << "\t" << result_count << endl;
        }
    }

    return 0;
}
