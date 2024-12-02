#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include <iomanip>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <unordered_map>
#include <openfhe.h>
#include <random>
#include <algorithm>
#include <thread>

#pragma once

using namespace std;
using namespace lbcrypto;
using namespace chrono;

inline string machine_name = "c5.24xlarge";
inline size_t TOTAL_MACHINE_THREAD = 48;
inline thread* MACHINE_THREADS = new thread[TOTAL_MACHINE_THREAD];
inline vector<size_t> start_end_threads;
inline vector<size_t> start_end_expected_threads;
inline uint32_t PRINT_SIDES = 5;
inline uint32_t PRINT_TOTAL = 9;
inline double fx_coeff_x3 = -0.5;
inline double fx_coeff_x = 1.5;
inline double gx_coeff_x3 = -1359/1024.0;
inline double gx_coeff_x = 2126/1024.0;
inline uint32_t multDepth = 50;
inline uint32_t ringDimBits = 16;
inline uint32_t total_rows = 32768;

// only for the purpose of collecting data deleting this variable after
inline uint32_t p_comps_var = 0;

// set this value to a minimum of 3 bits for multiple column evaluation
// any number of bits that is lower offers no multiplicative depth advantage
// increasing the amount of compositions and the multiplicative depth also improves performance
// For the most part, the fastest performance will come from single column compositions
inline uint32_t numDimBits = 14;
inline uint32_t numSlotBits = ringDimBits - 1;
inline uint32_t ringDim = 1 << ringDimBits;
inline uint32_t numDim = 1 << numDimBits;
inline uint32_t numSlots = 1 << numSlotBits;

// can choose between FIXEDMANUAL or any AUTO technique
// FIXEDMANUAL is less accurate than FLEXIBLE AUTO
inline ScalingTechnique rescaleTech = FLEXIBLEAUTO; // FIXEDMANUAL OR ANY AUTO TECHNIQUE
inline uint32_t dcrtBits = 40;
inline uint32_t firstMod = 41;

// estimate of 2^17 ring size 128 bit security num bits to allocate
inline double ring_2pow17_128bit_sec_est = 3544.0;
inline double ring_2pow16_128bit_sec = 1772.0;
inline double ring_2pow15_128bit_sec = 881.0;
inline int tpchq6_pred_mults_reduction = 2;
inline uint32_t depth;
inline SecretKeyDist secretKeyDist = UNIFORM_TERNARY;
inline mt19937_64 mt(time(nullptr));
inline string prev_ID = "0";
inline int cur_q_index = 0;
inline double min_prec = 1.0 / (1 << (numDimBits + 1));
inline string tableFilePath = "/home/ubuntu/";
inline string queryFilePath = "/home/ubuntu/OPENFHE_CKKS/";
inline int add_count_flag = 0;

inline CCParams<CryptoContextCKKSRNS>* parameters;
inline CryptoContext<DCRTPoly>* cc;
inline KeyPair<lbcrypto::DCRTPoly>* keys;

inline unordered_map<string, int> priority = 
{{"AND", 2}, {"OR", 1}, {"+", 1}, {"-", 1}, {"*", 2}};

inline vector<string> ops = {">=", "<=", ">", "<", "==", "!="};

// for added precision add one or two compositions
// especially when using multiple columns
// to g or f but these are the minimum requirements
// for sub 1 percent error for single column compositions
// fp and gp are for added precision in compositions
inline uint32_t fp = 0;
inline uint32_t gp = 0;
inline uint32_t p_comps = 0;
inline unordered_map<int, pair<uint32_t, uint32_t>> comp_ref;
inline bool depth_optimization = true;

inline uint64_t count_time, pred_comb_time, query_comp_time, query_xz_time, query_decomp_time, query_expand_time;
inline vector<uint64_t> agg_ops_time;
inline vector<uint64_t> selection_ops_time;
inline vector<uint64_t> client_agg_time;
inline unordered_map<string, uint64_t> pred_ops_time;
inline uint64_t end_to_end_time, total_pred_time, total_agg_time;
inline string query_results = "";
inline uint64_t resp_comp_time_count, resp_decomp_time_count;
inline vector<uint64_t> resp_comp_time_aggs;
inline vector<uint64_t> resp_decomp_time_aggs;

// vector formatting function taken from Microsoft SEAL examples
template <typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if (vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}

inline void define_compositions() {
    comp_ref = {{0, make_pair(0, 0)},
                {1, make_pair(2 + fp, 1 + gp)},
                {2, make_pair(2 + fp, 2 + gp)},
                {3, make_pair(2 + fp, 3 + gp)},
                {4, make_pair(3 + fp, 3 + gp)},
                {5, make_pair(3 + fp, 4 + gp)},
                {6, make_pair(3 + fp, 5 + gp)},
                {7, make_pair(3 + fp, 6 + gp)},
                {8, make_pair(3 + fp, 7 + gp)},
                {9, make_pair(3 + fp, 8 + gp)},
                {10, make_pair(3 + fp, 9 + gp)},
                {11, make_pair(3 + fp, 10 + gp)},
                {12, make_pair(3 + fp, 11 + gp)},
                {13, make_pair(3 + fp, 12 + gp)},
                {14, make_pair(3 + fp, 13 + gp)},
                {15, make_pair(3 + fp, 14 + gp)},
                {16, make_pair(3 + fp, 15 + gp)},
                {17, make_pair(3 + fp, 16 + gp)},
                {18, make_pair(3 + fp, 17 + gp)},
                {19, make_pair(3 + fp, 18 + gp)},
                {20, make_pair(3 + fp, 19 + gp)},
                {21, make_pair(4 + fp, 19 + gp)},
                {22, make_pair(4 + fp, 20 + gp)},
                {23, make_pair(4 + fp, 21 + gp)},
                {24, make_pair(4 + fp, 22 + gp)},
                {25, make_pair(4 + fp, 23 + gp)},
                {26, make_pair(4 + fp, 24 + gp)},
                {27, make_pair(4 + fp, 25 + gp)},
                {28, make_pair(4 + fp, 26 + gp)},
                {29, make_pair(4 + fp, 27 + gp)},
                {30, make_pair(4 + fp, 28 + gp)},
                {31, make_pair(4 + fp, 29 + gp)},
                {32, make_pair(4 + fp, 30 + gp)},
                {33, make_pair(4 + fp, 31 + gp)}};
}

inline void define_params(bool security128 = true, uint32_t desired_mult_depth = 0) {
    parameters = new CCParams<CryptoContextCKKSRNS>();
    (*parameters).SetSecretKeyDist(secretKeyDist);

    if (!security128) {
        (*parameters).SetSecurityLevel(HEStd_NotSet);
        (*parameters).SetRingDim(ringDim);

    } else {
        (*parameters).SetSecurityLevel(HEStd_128_classic);
        (*parameters).SetRingDim(ringDim);
    }

    multDepth = desired_mult_depth;
    auto approx = 0;
    if (ringDimBits == 17) {
        approx = min(ceil(ring_2pow17_128bit_sec_est / multDepth), 60.0);
        if (approx + (approx - 1) * (multDepth - 1) > ring_2pow17_128bit_sec_est) approx--;
    } else if (ringDimBits == 16) {
        approx = min(ceil(ring_2pow16_128bit_sec / multDepth), 60.0);
        if (approx + (approx - 1) * (multDepth - 1) > ring_2pow16_128bit_sec) approx--;
    } else if (ringDimBits == 15) {
        approx = min(ceil(ring_2pow15_128bit_sec / multDepth), 60.0);
        if (approx + (approx - 1) * (multDepth - 1) > ring_2pow15_128bit_sec) approx--;
    }
    firstMod = approx;
    dcrtBits = approx - 1;

    (*parameters).SetScalingModSize(dcrtBits);
    (*parameters).SetFirstModSize(firstMod);

    (*parameters).SetScalingTechnique(rescaleTech);

    cout << "\nComputed or desired multDepth: " << multDepth << endl;
    cout << "Newly computed firstMod: " << firstMod << endl;
    cout << "Newly computed dcrtbits: " << dcrtBits << endl << endl;

    (*parameters).SetMultiplicativeDepth(multDepth);
    (*parameters).SetBatchSize(numSlots);

    cc = new CryptoContext<DCRTPoly>(GenCryptoContext((*parameters)));

    (*cc)->Enable(PKE);
    (*cc)->Enable(KEYSWITCH);
    (*cc)->Enable(LEVELEDSHE);                   

    cout << "CKKS scheme is using ring dimension " << ringDim << endl << endl;

    keys = new KeyPair<lbcrypto::DCRTPoly>((*cc)->KeyGen());
    (*cc)->EvalMultKeyGen(keys->secretKey);

    vector<int32_t> indexList;
    for (size_t i = 0; i < numSlotBits; i++) {
        indexList.push_back((1 << i));
        indexList.push_back(-(1 << i));
    }

    (*cc)->EvalRotateKeyGen(keys->secretKey, indexList);   
}

inline Ciphertext<DCRTPoly> encrypt_vector(vector<double> &vec) {
    Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(vec);
    x_plain->SetLength(numSlots);
    return (*cc)->Encrypt(keys->publicKey, x_plain);
}

inline void print_ciphertext_and_expected(Ciphertext<DCRTPoly> &enc, vector<double> &input) {
    Plaintext result;
    cout.precision(8);

    (*cc)->Decrypt(keys->secretKey, enc, &result);
    result->SetLength(numSlots);
    cout << "expected = " << endl;
    print_vector(input, PRINT_SIDES, PRINT_TOTAL);
    vector<double> realNum;
    for (auto el: result->GetCKKSPackedValue()) {
        realNum.push_back(static_cast<double>(el.real()));
    }
    cout << "x_enc = " << endl;
    print_vector(realNum, PRINT_SIDES, PRINT_TOTAL);
    cout << "Estimated precision in bits: " << result->GetLogPrecision() << endl;
    auto levels_left = enc->GetLevel();
    cout << "Number of levels remaining: " << multDepth - levels_left << endl;
}

inline void print_ciphertext(Ciphertext<DCRTPoly>& enc) {
    Plaintext result;
    cout.precision(8);

    (*cc)->Decrypt(keys->secretKey, enc, &result);
    result->SetLength(numSlots);
    vector<double> realNum;
    for (auto el: result->GetCKKSPackedValue()) {
        realNum.push_back(static_cast<double>(el.real()));
    }
    cout << "x_enc = " << endl;
    print_vector(realNum, PRINT_SIDES, PRINT_TOTAL);
    cout << "Estimated precision in bits: " << result->GetLogPrecision() << endl;
    auto levels_left = enc->GetLevel();
    cout << "Number of levels remaining: " << multDepth - levels_left << endl;
}

inline void comp(uint32_t f_x, uint32_t g_x, Ciphertext<DCRTPoly>& enc,
 vector<double>& pub, bool transform, bool verbose = false) {
    // setting coefficients up for the multiplication
    double x3_coeff, x_coeff;
    if (f_x == 0 && g_x == 0) {
        if (transform) {
            enc = (*cc)->EvalAdd(enc, 0.5);
            for (uint64_t i = 0; i < pub.size(); i++) pub[i] += 0.5;
        }
        return;
    } else if (g_x != 0) x3_coeff = gx_coeff_x3, x_coeff = gx_coeff_x;
    else if (f_x == 1 && transform) x3_coeff = fx_coeff_x3 / 2.0, x_coeff = fx_coeff_x / 2.0;
    else x3_coeff = fx_coeff_x3, x_coeff = fx_coeff_x;

    auto coeff3_x = (*cc)->EvalMult(enc, x3_coeff);
    if (rescaleTech == FIXEDMANUAL) coeff3_x = (*cc)->Rescale(coeff3_x);
    
    auto x2 = (*cc)->EvalMult(enc, enc);
    if (rescaleTech == FIXEDMANUAL) x2 = (*cc)->Rescale(x2);

    auto coeff3_x3 = (*cc)->EvalMult(x2, coeff3_x);
    if (rescaleTech == FIXEDMANUAL) coeff3_x3 = (*cc)->Rescale(coeff3_x3);

    auto coeff_x = (*cc)->EvalMult(enc, x_coeff);
    if (rescaleTech == FIXEDMANUAL) coeff_x = (*cc)->Rescale(coeff_x);

    // set_to_same_level(coeff3_x3, coeff_x);
    enc = (*cc)->EvalAdd(coeff3_x3, coeff_x);

    for (uint64_t i = 0; i < pub.size(); i++) {
        pub[i] = (x3_coeff * pub[i] *  pub[i] + x_coeff) * pub[i];
    }

    if (verbose) print_ciphertext_and_expected(enc, pub);

    if (g_x != 0) comp(f_x, g_x - 1, enc, pub, transform, verbose);
    else comp(f_x - 1, 0, enc, pub, transform, verbose);
}

inline void pred_comp_ENC(uint32_t iterations, Ciphertext<DCRTPoly>& pred, bool verbose = false) {

    if (iterations == 0) return;

    auto neg2_x = (*cc)->EvalMult(pred, -2);
    if (rescaleTech == FIXEDMANUAL) neg2_x = (*cc)->Rescale(neg2_x);

    auto neg2_x_plus3 = (*cc)->EvalAdd(neg2_x, 3);

    auto x2 = (*cc)->EvalMult(pred, pred);
    if (rescaleTech == FIXEDMANUAL) x2 = (*cc)->Rescale(x2);

    pred = (*cc)->EvalMult(neg2_x_plus3, x2);
    if (rescaleTech == FIXEDMANUAL) pred = (*cc)->Rescale(pred);

    if (verbose) print_ciphertext(pred);

    pred_comp_ENC(iterations - 1, pred, verbose);
}

inline void pred_comp_PUB(uint32_t iterations, vector<double>& pred, bool verbose = false) {

    if (iterations == 0) return;

    for (uint64_t i = 0; i < pred.size(); i++) {
        pred[i] = (-2 * pred[i] + 3) * pred[i] * pred[i];
    }

    if (verbose) print_vector(pred);

    pred_comp_PUB(iterations - 1, pred, verbose);
}

inline void comp_ENC(uint32_t f_x, uint32_t g_x, Ciphertext<DCRTPoly>& enc,
 bool transform, bool verbose = false) {
    // setting coefficients up for the multiplication
    double x3_coeff, x_coeff;
    if (f_x == 0 && g_x == 0) {
        if (transform) {
            enc = (*cc)->EvalAdd(enc, 0.5);
        }
        return;
    } else if (g_x != 0) x3_coeff = gx_coeff_x3, x_coeff = gx_coeff_x;
    else if (f_x == 1 && transform) x3_coeff = fx_coeff_x3 / 2.0, x_coeff = fx_coeff_x / 2.0;
    else x3_coeff = fx_coeff_x3, x_coeff = fx_coeff_x;

    auto coeff3_x = (*cc)->EvalMult(enc, x3_coeff);
    if (rescaleTech == FIXEDMANUAL) coeff3_x = (*cc)->Rescale(coeff3_x);

    auto x2 = (*cc)->EvalMult(enc, enc);
    if (rescaleTech == FIXEDMANUAL) x2 = (*cc)->Rescale(x2);

    auto coeff3_x3 = (*cc)->EvalMult(x2, coeff3_x);
    if (rescaleTech == FIXEDMANUAL) coeff3_x3 = (*cc)->Rescale(coeff3_x3);

    auto coeff_x = (*cc)->EvalMult(enc, x_coeff);
    if (rescaleTech == FIXEDMANUAL) coeff_x = (*cc)->Rescale(coeff_x);

    // set_to_same_level(coeff3_x3, coeff_x);
    enc = (*cc)->EvalAdd(coeff3_x3, coeff_x);

    if (verbose) print_ciphertext(enc);

    if (g_x != 0) comp_ENC(f_x, g_x - 1, enc, transform, verbose);
    else comp_ENC(f_x - 1, 0, enc, transform, verbose);
}

inline void comp_PUB(uint32_t f_x, uint32_t g_x,
 vector<double>& pub, bool transform, bool verbose = false) {
    // setting coefficients up for the multiplication
    double x3_coeff, x_coeff;
    if (f_x == 0 && g_x == 0) {
        if (transform) {
            for (uint64_t i = 0; i < pub.size(); i++) pub[i] += 0.5;
        }
        return;
    } else if (g_x != 0) x3_coeff = gx_coeff_x3, x_coeff = gx_coeff_x;
    else if (f_x == 1 && transform) x3_coeff = fx_coeff_x3 / 2.0, x_coeff = fx_coeff_x / 2.0;
    else x3_coeff = fx_coeff_x3, x_coeff = fx_coeff_x;

    for (uint64_t i = 0; i < pub.size(); i++) {
        pub[i] = (x3_coeff * pub[i] *  pub[i] + x_coeff) * pub[i];
    }

    if (verbose) print_vector(pub);

    if (g_x != 0) comp_PUB(f_x, g_x - 1, pub, transform, verbose);
    else comp_PUB(f_x - 1, 0, pub, transform, verbose);
}

inline void min_precision() {
    cout << "lowest fractional measurement for comparisons: " << 1 << "/" << pow(2, numDimBits + 1) << endl;
    cout << "lowest fractional measurement for comparisons containing both positive and negative numbers: " << 1 << "/" << pow(2, numDimBits + 2) << endl;
}

inline void test_comp(int f, int g, int log_range) {
    vector<double> pub;
    // pushes back the minimum fractional differences between numbers
    // when they are mapped to [0, 1] to ensure that there are a 
    // sufficient amount of compositions to evaluate a and b
    // within a certain range of numbers
    auto l_min_prec = 1.0 / (1 << log_range);
    cout << "min precision is: " << l_min_prec << endl;
    for (uint32_t i = 0; i < numSlots; i++) {
        // if a - b == 1
        if (i % 5 == 0) pub.push_back(l_min_prec + 2 * l_min_prec);
        // if a == b
        else if (i % 5 == 1) pub.push_back(l_min_prec + 0 * l_min_prec);
        // max pos val in ciphertext
        else if (i % 5 == 2) pub.push_back(1 - l_min_prec);
        // max neg val in ciphertext
        else if (i % 5 == 3) pub.push_back(l_min_prec - 1);
        // if a - b == -1
        else pub.push_back(l_min_prec - 2 * l_min_prec);
    }
    Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(pub);
    x_plain->SetLength(numSlots);

    auto x_enc = (*cc)->Encrypt(keys->publicKey, x_plain);
    std::cout << "cur ciphertext level: " << x_enc->GetLevel() << std::endl;
    print_ciphertext_and_expected(x_enc, pub);

    comp(f, g, x_enc, pub, true);

    std::cout << "number of levels remaining: " << multDepth - x_enc->GetLevel() << std::endl;

    print_ciphertext_and_expected(x_enc, pub);
}

inline void test_comp() {
    vector<double> pub;
    // pushes back the minimum fractional differences between numbers
    // when they are mapped to [0, 1] to ensure that there are a 
    // sufficient amount of compositions to evaluate a and b
    // within a certain range of numbers
    for (uint32_t i = 0; i < numSlots; i++) {
        // if a - b == 1
        if (i % 3 == 0) pub.push_back(min_prec + 2 * min_prec);
        // if a == b
        else if (i % 3 == 1) pub.push_back(min_prec + 0 * min_prec);
        // if a - b == -1
        else pub.push_back(min_prec - 2 * min_prec);
    }
    Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(pub);
    x_plain->SetLength(numSlots);

    auto x_enc = (*cc)->Encrypt(keys->publicKey, x_plain);
    std::cout << "cur ciphertext level: " << x_enc->GetLevel() << std::endl;
    auto comps = comp_ref[numDimBits + 1];
    comp(comps.first, comps.second, x_enc, pub, true);

    std::cout << "number of levels remaining: " << multDepth - x_enc->GetLevel() << std::endl;

    print_ciphertext_and_expected(x_enc, pub);
}

inline string next_ID() {
    prev_ID = to_string(stoi(prev_ID) + 1);
    return prev_ID;
}

// https://www.geeksforgeeks.org/cpp-string-to-vector-using-delimiter/
inline vector<string> split(string str, string delimiter)
{
    vector<string> v;
    if (!str.empty()) {
        int start = 0;
        do {
            auto idx = str.find(delimiter, start);
            if (idx == string::npos) {
                break;
            }

            int length = idx - start;
            v.push_back(str.substr(start, length));
            start += (length + delimiter.size());
        } while (true);
        v.push_back(str.substr(start));
    }
 
    return v;
}

// from https://leetcode.com/problems/number-of-days-between-two-dates/solutions/621582/c-short-simple-modular-solution/
inline bool is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
}

inline int days_in_month(int m, int year) { 
    if(m==1 || m==3 || m==5 || m==7 || m==8 || m==10 || m==12 ) return 31;
    if(m==2) return is_leap_year(year) ? 29 : 28;
    return 30;
}

inline int date_to_int(int cur_date) {
    int Y = cur_date / 10000;
    int M = (cur_date / 100) % 100;
    int D = cur_date % 100;
    
    int date = 0;
    for(int y = 1970; y < Y; ++y) date += is_leap_year(y) ? 366 : 365;
    for(int m = 1; m < M; ++m) date += days_in_month(m, Y);
    return date + D;
}

inline int daysSinceJan1st1970(int date) {
    return date_to_int(date) - date_to_int(19700101);
}

inline void reset_params() {
    MACHINE_THREADS = new thread[TOTAL_MACHINE_THREAD];
    start_end_threads.clear();
    start_end_expected_threads.clear();
    // delete parameters;
    delete cc;
    delete keys;
    count_time = pred_comb_time = query_comp_time = query_xz_time = query_decomp_time = query_expand_time = 0;
    agg_ops_time.clear();
    selection_ops_time.clear();
    client_agg_time.clear();
    pred_ops_time.clear();
    end_to_end_time = total_pred_time = total_agg_time = 0;
    query_results = "";
    resp_comp_time_count = resp_decomp_time_count = 0;
    resp_comp_time_aggs.clear();
    resp_decomp_time_aggs.clear();
    prev_ID = "0";
    cur_q_index = 0;
}