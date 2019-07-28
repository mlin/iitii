#include <iostream>
#include <utility>
#include <random>
#include <math.h>
#include <string>
#include <thread>
#include <assert.h>
#include "iitii.h"
#include "ctpl_stl.h"
#include "kstring.h"
#include "tbx.h"
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

// split s on delim & return strlen(s). s is damaged by side-effect
template<typename Out>
size_t split(char *s, char delim, Out result, uint64_t maxsplit = ULLONG_MAX) {
    string delims("\n");
    delims[0] = delim;
    char *cursor = s;
    char *token = strsep(&cursor, delims.c_str());
    char *last_token = token;
    uint64_t i = 0;
    while (token) {
        *(result++) = last_token = token;
        if (++i < maxsplit) {
            token = strsep(&cursor, delims.c_str());
        } else {
            if (cursor) {
                *result = last_token = cursor;
            }
            break;
        }
    }
    return last_token ? (last_token + strlen(last_token) - s) : 0;
}

template<typename Out>
size_t split(string& s, char delim, Out result, uint64_t maxsplit = ULLONG_MAX) {
    return split(&s[0], delim, result, maxsplit);
}

class TabixIterator {
    htsFile *fp_;
    tbx_t *tbx_;

    hts_itr_t *it_;
    kstring_t str_ = {0,0,0};
    bool valid_ = false;

    TabixIterator(htsFile *fp, tbx_t *tbx, hts_itr_t *it) {
        fp_ = fp;
        tbx_ = tbx;
        it_ = it;
        assert(tbx_ && fp_ && it_);
    }

public:
    ~TabixIterator() {
        tbx_itr_destroy(it_);
        if (str_.s) {
            free(str_.s);
        }
    }

    static unique_ptr<TabixIterator> Open(htsFile *fp, tbx_t *tbx, int tid, int beg, int end) {
        if (!tbx) {
            return nullptr;
        }
        hts_itr_t *it = tbx_itr_queryi(tbx, tid, beg, end);
        if (!it) {
            return nullptr;
        }
        auto ans = unique_ptr<TabixIterator>(new TabixIterator(fp, tbx, it));
        ans->Next();
        return ans;
    }

    bool Valid() const {
        return valid_ && str_.s;
    }

    const char* Line() const {
        if (Valid()) {
            assert(ks_len((kstring_t*)&str_) == strlen(str_.s));
            return str_.s;
        }
        return nullptr;
    }

    void Next() {
        valid_ = (tbx_itr_next(fp_, tbx_, it_, &str_) >= 0);
        assert(!valid_ || str_.s);
    }
};

struct variant {
    int rid;
    int beg;
    int end;
    string id;
    string ref;
    vector<string> alt;

    string str() const {
        string ans = to_string(beg) + "/" + ref + "/";
        bool first = true;
        for (const auto& one_alt : alt) {
            if (!first) {
                ans += ",";
            }
            ans += one_alt;
        }
        if (id.size()) {
            ans += "/" + id;
        }
        return ans;
    }

    bool operator==(const variant& rhs) const {
        return rid == rhs.rid && beg == rhs.beg && ref == rhs.ref;
    }
};

int variant_beg(const variant& v) { return v.beg; }
int variant_end(const variant& v) { return v.end; }

vector<variant> load_variants(const string& filename, int rid, int beg, int end) {
    auto fp = shared_ptr<htsFile>(hts_open(filename.c_str(), "r"),
                                  [] (htsFile* f) { if (f && hts_close(f)) throw runtime_error("hts_close"); });
    if (!fp) {
        throw runtime_error("Failed to open " + filename);
    }
    auto tbx = shared_ptr<tbx_t>(tbx_index_load(filename.c_str()), [] (tbx_t *t) { if (t) tbx_destroy(t); });
    if (!tbx) {
        throw runtime_error("Failed to open .tbi/.csi index of " + filename);
    }
    auto itr = TabixIterator::Open(fp.get(), tbx.get(), rid, beg, end);
    if (!itr) {
        throw runtime_error("Failed to seek position " + to_string(beg) + " using index of " + filename);
    }

    vector<char*> tokens;
    vector<variant> variants;
    for (; itr->Valid(); itr->Next()) {
        string linecpy(itr->Line());
        tokens.clear();
        split(linecpy, '\t', back_inserter(tokens), 6);
        if (tokens.size() < 6) {
            throw runtime_error("invalid VCF line: " + string(itr->Line()));
        }

        variant vt;
        vt.beg = atoi(tokens[1]);
        if (vt.beg >= beg && vt.beg < end) {  // exclude danglers
            vt.rid = rid;
            vt.id = strcmp(tokens[2], ".") ? string(tokens[2]) : "";
            vt.ref = string(tokens[3]);
            vt.end = vt.beg + vt.ref.size();
            string altcpy(tokens[4]);
            tokens.clear();
            split(altcpy, ',', back_inserter(tokens));
            for (auto alt : tokens) {
                vt.alt.push_back(string(alt));
            }
            variants.push_back(move(vt));
        }
    }
    //cout << beg << " " << end << " " << variants.size() << endl;
    return variants;
}

vector<variant> load_variants_parallel(const string& filename, int rid, int megabases) {
    ctpl::thread_pool pool(thread::hardware_concurrency());
    mutex mu;
    vector<future<void>> futures;
    vector<variant> ans;

    for (int mb = 0; mb < megabases; ++mb) {
        auto fut = pool.push([&, mb](int tid) {
            auto variants_i = load_variants(filename, rid, mb*1000000, (mb+1)*1000000);
            lock_guard<mutex> lock(mu);
            ans.insert(ans.end(), variants_i.begin(), variants_i.end());
        });
        futures.push_back(move(fut));
    }

    for (auto& fut : futures) {
        fut.get();
    }

    return ans;
}

TEST_CASE("gnomAD chr2") {
    const int rid = 0;
    #ifdef NDEBUG
    const int megabases = 240;
    #else
    const int megabases = 24;
    #endif
    const string filename = "/tmp/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz";
    const string url = "https://storage.googleapis.com/gnomad-public/release/2.0.2/vcf/genomes/gnomad.genomes.r2.0.2.sites.chr2.vcf.bgz";

    ifstream vcf(filename), tbi(filename + ".tbi");
    if (vcf.bad() || tbi.bad()) {
        WARN("Skipping test because " + filename + " and .tbi aren't present. Download them to that location from " + url);
    } else {
        auto variants = load_variants_parallel(filename, rid, megabases);
        int max_len = -1, max_end = -1;
        for (const auto& vt : variants) {
            max_len = std::max(max_len, vt.end - vt.beg);
            max_end = std::max(max_end, vt.end);
        }
        cout << "Loaded " << variants.size() << " variants from first " << megabases << "Mbp, max len =  " << max_len << endl;

        default_random_engine R(42);
        uniform_int_distribution<uint32_t> begD(0, max_end);
        const size_t trials = 10000;

        iit<int, variant, variant_beg, variant_end> tree(variants.begin(), variants.end());
        iitii<int, variant, variant_beg, variant_end> treeii(variants.begin(), variants.end());
        size_t cost = 0, costii=0, count=0;
        vector<variant> results, resultsii;

        for (size_t i = 0; i < trials; i++) {
            auto qbeg = begD(R);
            auto qend = qbeg+10;
            cost += tree.overlap(qbeg, qend, results);
            costii += treeii.overlap(qbeg, qend, resultsii);
            REQUIRE(results.size() == resultsii.size());
            count += results.size();
            bool alleq = true;
            for (auto p1 = results.begin(), p2 = resultsii.begin(); p1 != results.end(); ++p1, ++p2) {
                alleq = alleq && (*p1 == *p2);
            }
            REQUIRE(alleq);
        }

        cout << count << " " << cost << " " << costii << endl;
        cout << "mean absolute error of rank prediction = " << double(treeii.rank_error)/treeii.queries << endl;

        std::sort(variants.begin(), variants.end(), [](const variant& lhs, const variant& rhs) {
            auto begl = lhs.beg, begr = rhs.beg;
            if (begl == begr) {
                return lhs.end < rhs.end;
            }
            return begl < begr;
        });
        ofstream variants_txt("variants.txt");
        for (const auto& vt : variants) {
            variants_txt << vt.str() << endl;
        }
    }
}
