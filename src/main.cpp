#include <iostream>
#include <vector>
#include <random>
#include "args.hxx"
#include "iitii.h"

using intpair = std::pair<uint64_t,uint64_t>;
uint64_t p_get_beg(const intpair& p) { return p.first; }
uint64_t p_get_end(const intpair& p) { return p.second; }
using p_iitii = iitii::iitii<uint64_t, intpair, p_get_beg, p_get_end>;  // first arg is position type

int main(int argc, char** argv) {

    args::ArgumentParser parser("memmapped interpolated implicit interval tree");
    args::HelpFlag help(parser, "help", "display this help summary", {'h', "help"});
    args::ValueFlag<std::string> test_file(parser, "FILE", "test mmmultimap with random data in this file", {'T', "test-file"});
    args::ValueFlag<uint64_t> test_size(parser, "N", "test this many pairs", {'s', "test-size"});
    args::ValueFlag<uint64_t> max_val(parser, "N", "generate test data in the range [1,max_value]", {'M', "max-value"});
    args::ValueFlag<uint64_t> range_mean(parser, "N", "the mean length for intervals (under gaussian distribution)", {'m', "range-mean"});
    args::ValueFlag<uint64_t> range_stdev(parser, "N", "the standard deviation for intervals (under gaussian distribution)", {'D', "range-stdev"});
    args::ValueFlag<uint64_t> threads(parser, "N", "number of threads to use", {'t', "threads"});
    args::ValueFlag<uint64_t> domains(parser, "N", "number of domains for interpolation", {'d', "domains"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }
    
    assert(!args::get(test_file).empty());
    assert(args::get(test_size));
    assert(args::get(max_val));
    
    std::remove(args::get(test_file).c_str());
    p_iitii::builder bb = p_iitii::builder(args::get(test_file));

    //bb.add(intpair(12,34));
    //bb.add(intpair(0,23));
    //bb.add(intpair(34,56));

    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    uint64_t max_value = args::get(max_val);
    std::uniform_int_distribution<uint64_t> dis(0, max_value);
    std::normal_distribution<> dlen(args::get(range_mean),args::get(range_stdev));
    uint64_t x_len = args::get(test_size);
#pragma omp parallel for
    for (int n=0; n<x_len; ++n) {
        uint64_t q = dis(gen);
        uint64_t r = std::min(q + std::max((uint64_t)0, (uint64_t)std::round(dlen(gen))), max_value);
        //uint64_t a = std::min(q, r);
        //uint64_t b = std::max(q, r);
        //std::cerr << q << ", " << r << std::endl;
        bb.add(intpair(q, r));
    }
    //tree.index();
    uint64_t n_domains = std::max((uint64_t)1, (uint64_t)args::get(domains));
    p_iitii db = bb.build(n_domains);
//#pragma omp parallel for
    for (int n=0; n<max_value; ++n) {
        std::vector<intpair> ovlp = db.overlap(n, n+1);
        if (n % 1000 == 0) std::cerr << n << "\r";
        //std::cerr << n << " has " << ovlp.size() << " overlaps" << std::endl;
        for (auto& s : ovlp) {
            if (s.first > n || s.second < n) {
                std::cerr << "tree broken at " << n << std::endl;
            }
        }
    }

    //std::vector<intpair> results = db.overlap(22, 25);
    // alternative: db.overlap(22, 25, results);

    //for (const auto& p : results)
    //std::cout << p.first << "\t" << p.second << std::endl;
    return 0;
}
