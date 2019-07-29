#include "util.h"
#include <fstream>
#include <random>
#include <chrono>
#include <functional>

size_t milliseconds_to(function<void()> f) {
    chrono::steady_clock::time_point begin = chrono::steady_clock::now();
    f();
    chrono::steady_clock::time_point end = chrono::steady_clock::now();
    return chrono::duration_cast<chrono::milliseconds>(end - begin).count();
}

template <class tree>
size_t run_queries(const vector<variant>& variants, const tree& t, int max_end, int queries) {
    default_random_engine R(42);
    uniform_int_distribution<uint32_t> begD(0, max_end);
    uniform_int_distribution<size_t> vtD(0, variants.size()-1);
    size_t ans = 0;
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
        auto results = t.overlap(qbeg, qend);
        ans += results.size();
    }
    return ans;
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
    cerr << "Loaded " << variants.size() << " variants, max END = " << max_end
         << ", max rlen = " << max_len << endl;

    unique_ptr<variant_iit> tree;
    cout << "milliseconds to build iit: " << milliseconds_to([&](){
        auto t = variant_iit::builder(variants.begin(), variants.end()).build();
        tree.reset(new variant_iit(move(t)));
    }) << endl;

    unique_ptr<variant_iitii> treeii;
    cout << "milliseconds to build iitii: "  << milliseconds_to([&](){
        // configure iitii with one interpolation domain per 100kbp
        auto t = variant_iitii::builder(variants.begin(), variants.end()).build(megabases*10);
        treeii.reset(new variant_iitii(move(t)));
    }) << endl;

    const size_t trials = 10000000;
    size_t result_count = 0;
    cout << "milliseconds for iit queries: " << milliseconds_to([&](){
        result_count = run_queries<variant_iit>(variants, *tree, max_end, trials);
    }) << endl;

    cout << "milliseconds for iitii queries: " << milliseconds_to([&](){
        if (run_queries<variant_iitii>(variants, *treeii, max_end, trials) != result_count) {
            throw runtime_error("PANIC: result sets differed!");
        }
    }) << endl;

    cout << "mean climbing per iitii query: " << double(treeii->total_climb_cost)/treeii->queries << endl;

    return 0;
}
