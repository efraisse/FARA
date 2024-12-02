#include "globals.h"
#include "db.h"
#include "q.h"
#include "header.h"

void test_TPCH_6(int num_rows, bool pub_res = true, bool exp_res = true, bool depth_opt = true) {
    total_rows = num_rows;
    depth_optimization = depth_opt;
    cout << "Testing TPC-H 6 query with " << num_rows << " rows" << endl;
    cout << "Testing on Machine " << machine_name << endl;
    cout << "Using " << TOTAL_MACHINE_THREAD << " threads" << endl;
    numDimBits = 14;
    numDim = 1 << numDimBits;
    ringDimBits = 16;
    ringDim = 1 << ringDimBits;
    numSlotBits = ringDimBits - 1;
    numSlots = 1 << numSlotBits;
    fp = 0;
    gp = 1;
    p_comps = p_comps_var;
    min_precision();

    header quanlity_data("quanlity_data", 0, 63);
    header extended_price("extended_price", 0, 10500000, 2);
    header discount_data("discount_data", 0, 15, 2);
    header ship_date("ship_date", daysSinceJan1st1970(19920101), daysSinceJan1st1970(19981231));

    vector<header> headers = {quanlity_data, extended_price, discount_data, ship_date};

    subQuery sub1 = subQuery(ship_date, "<", daysSinceJan1st1970(19950101));
    subQuery sub2 = subQuery(ship_date, ">=", daysSinceJan1st1970(19940101));
    subQuery sub3 = subQuery(discount_data, "<=", 12);
    subQuery sub4 = subQuery(discount_data, ">=", 6);
    subQuery sub5 = subQuery(quanlity_data, "<", 15);

    vector<subQuery> subs = {sub1,
                            subQuery("AND"),
                            sub2,
                            subQuery("AND"),
                            sub3,
                            subQuery("AND"),
                            sub4,
                            subQuery("AND"),
                            sub5};

    vector<vector<subQuery>> t_aggs;
    vector<string> final_ops;

    vector<subQuery> aggs1 = {subQuery(extended_price),
                            subQuery("*"),
                            subQuery(discount_data)};

    final_ops = {"SUM"};
    t_aggs = {aggs1};

    query q = query(final_ops, subs, t_aggs);

    define_compositions();

    // experimental min depth calculator
    // auto estimated_depth = q.calculate_min_depth(total_rows);
    // estimated_depth -= tpchq6_pred_mults_reduction;
    // define_params(false, estimated_depth);

    // define_params(false, 40);
    define_params(false, 36 + p_comps_var * 2);
    p_comps_var++;

    db TPCHdb(num_rows, headers, pub_res, exp_res);
    TPCHdb.fill_db_TPCH({4,5,6,10});

    // TPCHdb.print_db();
    cout << "formulating the query to send to the server" << endl;

    q.formulate_query();

    auto q_comp_time = chrono::high_resolution_clock::now();
    auto return_code = Serial::SerializeToFile("comp_query.txt", q.getCompressedQuery(), SerType::BINARY);
    query_comp_time = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - q_comp_time)).count();
    cout << "The size of the compressed query is: " << std::filesystem::file_size(queryFilePath + "comp_query.txt") << " bytes.\n";
    assert(return_code == 1);

    // xz compression provides marginal query size improvement
    // for some configurations can provide 20-30 percent additional size compression
    // not very worthwile when scaling the database to millions of rows
    // 
    // auto q_xz_time = chrono::high_resolution_clock::now();
    // return_code = system("xz comp_query.txt");
    // assert(return_code == 0);
    // query_xz_time =  (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - q_xz_time)).count();

    // cout << "The size of the compressed xz query is: " << std::filesystem::file_size(queryFilePath + "comp_query.txt.xz") << " bytes.\n";

    auto q_decomp_time = chrono::high_resolution_clock::now();
    Ciphertext<DCRTPoly> decomp_query;

    // return_code = system("xz -d comp_query.txt.xz");
    // assert(return_code == 0);

    return_code = Serial::DeserializeFromFile("comp_query.txt", decomp_query, SerType::BINARY);
    assert(return_code == 1);
    q.setCompressedQuery(decomp_query);
    query_decomp_time =  (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - q_decomp_time)).count();

    auto q_exp_time = chrono::high_resolution_clock::now();
    q.query_expansion();
    query_expand_time = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - q_exp_time)).count();
    
    cout << "starting query evaluation" << endl;
    TPCHdb.evaluate_query(q);
}

void test_4bit_comp(int num_rows, bool pub_res = true, bool exp_res = true) {
    total_rows = num_rows;
    cout << "Testing 4 bits with " << num_rows << " rows" << endl;
    cout << "Testing on Machine " << machine_name << endl;
    cout << "Using " << TOTAL_MACHINE_THREAD << " threads" << endl;
    numDimBits = 4;
    numDim = 1 << numDimBits;
    fp = 0;
    gp = 0;
    // p_comps = 1;
    min_precision();

    header test_data("test_data", 0, (1 << 4) - 1);

    vector<header> headers = {test_data};

    subQuery sub1 = subQuery(test_data, "==", (1 << 3) - 1);
    vector<subQuery> subs = {sub1};

    vector<vector<subQuery>> t_aggs;
    vector<string> final_ops;

    vector<subQuery> aggs1 = {subQuery(test_data)};

    final_ops = {"SUM"};
    t_aggs = {aggs1};

    query q = query(final_ops, subs, t_aggs);
    
    define_compositions();
    define_params(false, 19);

    db TPCHdb(num_rows, headers, pub_res, exp_res);
    TPCHdb.fill_db_random_values();

    // TPCHdb.print_db();
    cout << "formulating the query to send to the server" << endl;

    q.formulate_query();
    q.query_expansion();
    
    cout << "starting query evaluation" << endl;
    TPCHdb.evaluate_query(q);
}

void test_16bit_comp(int num_rows, bool pub_res = true, bool exp_res = true) {
    total_rows = num_rows;
    cout << "Testing 16 bits with " << num_rows << " rows" << endl;
    cout << "Testing on Machine " << machine_name << endl;
    cout << "Using " << TOTAL_MACHINE_THREAD << " threads" << endl;
    numDimBits = 16;
    numDim = 1 << numDimBits;
    fp = 0;
    gp = 0;
    min_precision();

    header test_data("test_data", -(1 << 15), (1 << 15) - 1);

    vector<header> headers = {test_data};

    subQuery sub1 = subQuery(test_data, "==", 0);

    vector<subQuery> subs = {sub1};

    vector<vector<subQuery>> t_aggs;
    vector<string> final_ops;

    vector<subQuery> aggs1 = {subQuery(test_data)};

    final_ops = {"SUM"};
    t_aggs = {aggs1};

    query q = query(final_ops, subs, t_aggs);
    
    define_compositions();
    define_params(false, 43);

    db TPCHdb(num_rows, headers, pub_res, exp_res);
    TPCHdb.fill_db_random_values();

    TPCHdb.print_db();
    cout << "formulating the query to send to the server" << endl;

    q.formulate_query();
    q.query_expansion();
    
    cout << "starting query evaluation" << endl;
    TPCHdb.evaluate_query(q);
}

void test_32bit_comp(int num_rows, bool pub_res = true, bool exp_res = true) {
    total_rows = num_rows;
    cout << "Testing 32 bits with " << num_rows << " rows" << endl;
    cout << "Testing on Machine " << machine_name << endl;
    cout << "Using " << TOTAL_MACHINE_THREAD << " threads" << endl;
    numDimBits = 11;
    numDim = 1 << numDimBits;
    ringDimBits = 16;
    ringDim = 1 << ringDimBits;
    numSlotBits = ringDimBits - 1;
    numSlots = 1 << numSlotBits;
    min_precision();

    header test_data("test_data", -(long long)(pow(2, 31)), (long long)(pow(2, 31) - 1));

    vector<header> headers = {test_data};

    subQuery sub1 = subQuery(test_data, "==", 0);

    vector<subQuery> subs = {sub1};

    vector<string> final_ops;

    final_ops = {};

    query q = query(final_ops, subs);

    // multiplying by plaintext for aggregation does not work
    // because the plaintext is large so DECRYPTION FAILS
    // due to accuracy standards from OpenFHE library
    // therefore, we probably need a ring size of 2^17 to actually be able to perform this.

    // vector<subQuery> subs = {sub1};

    // vector<vector<subQuery>> t_aggs;
    // vector<string> final_ops;

    // vector<subQuery> aggs1 = {subQuery(test_data)};

    // final_ops = {"SUM"};
    // t_aggs = {aggs1};

    // query q = query(final_ops, subs, t_aggs);
    
    define_compositions();
    define_params(false, 45);

    db TPCHdb(num_rows, headers, pub_res, exp_res);
    TPCHdb.fill_db_random_values();

    // TPCHdb.print_db();
    cout << "formulating the query to send to the server" << endl;

    q.formulate_query();
    q.query_expansion();
    
    cout << "starting query evaluation" << endl;
    TPCHdb.evaluate_query(q);
}

void test_64bit_comp(int num_rows, bool pub_res = true, bool exp_res = true) {
    total_rows = num_rows;
    cout << "Testing 32 bits with " << num_rows << " rows" << endl;
    cout << "Testing on Machine " << machine_name << endl;
    cout << "Using " << TOTAL_MACHINE_THREAD << " threads" << endl;
    numDimBits = 13;
    numDim = 1 << numDimBits;
    ringDimBits = 17;
    ringDim = 1 << ringDimBits;
    numSlotBits = ringDimBits - 1;
    numSlots = 1 << numSlotBits;
    min_precision();

    header test_data("test_data", -(long long)(pow(2, 63)), (long long)(pow(2, 63) - 1));

    vector<header> headers = {test_data};

    subQuery sub1 = subQuery(test_data, "==", 0);

    vector<subQuery> subs = {sub1};

    vector<string> final_ops;

    final_ops = {};

    query q = query(final_ops, subs);

    // multiplying by plaintext for aggregation does not work
    // because the plaintext is large so DECRYPTION FAILS
    // due to accuracy standards from OpenFHE library
    // therefore, we probably need a ring size of 2^17 to actually be able to perform this.

    // vector<subQuery> subs = {sub1};

    // vector<vector<subQuery>> t_aggs;
    // vector<string> final_ops;

    // vector<subQuery> aggs1 = {subQuery(test_data)};

    // final_ops = {"SUM"};
    // t_aggs = {aggs1};

    // query q = query(final_ops, subs, t_aggs);
    
    define_compositions();
    define_params(false, 53);

    db TPCHdb(num_rows, headers, pub_res, exp_res);
    TPCHdb.fill_db_random_values();

    // TPCHdb.print_db();
    cout << "formulating the query to send to the server" << endl;

    q.formulate_query();
    q.query_expansion();
    
    cout << "starting query evaluation" << endl;
    TPCHdb.evaluate_query(q);
}

int main(int, char **)
{
    // int num_tests = 4;
    // for (int i = 0; i < num_tests; i++) {
    //     cout << "Running Test #" << i + 1 << ": \n\n";
    //     test_TPCH_6((1 << 15) * 1, true, true, true);
    //     reset_params();
    //     cout << "\n\n";
    // }
    test_4bit_comp((1 << 15) * 1, true, true);
}