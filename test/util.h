#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <utility>
#include <assert.h>
#include "kstring.h"
#include "tbx.h"
#include "ctpl_stl.h"
#include "iitii.h"

using namespace std;

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

// sugar for seeking into a htsFile with tabix & iterating through the relevant lines
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

// iit Item for testing with VCF variants
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

// use tabix to open VCF file & load variants overlapping given region
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
    return variants;
}

// read the first 'megabases' Mbp worth of variants from the file quickly, by parallelizing each
// megabase on a thread pool. results are unordered
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

using variant_iit = iit<int, variant, variant_beg, variant_end>;
using variant_iitii = iitii<int, variant, variant_beg, variant_end>;
