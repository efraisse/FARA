#include "globals.h"
#include "q.h"

#pragma once

class db {
    public:
        db(uint32_t num_rows, vector<header> headers, bool pub_res, bool exp_res);

        void print_initial_db();

        void print_db();

        void print_expected_db();

        void fill_db_random_values();

        void evaluate_query(query& q);

        vector<Ciphertext<DCRTPoly>> evaluate_sub_query_ENC(subQuery& sq, query& q);

        vector<vector<double>> evaluate_sub_query_PUB(subQuery sq);

        vector<Ciphertext<DCRTPoly>> evaluate_agg_query_ENC(subQuery sq, vector<Ciphertext<DCRTPoly>>& countRes);

        vector<vector<double>> evaluate_agg_query_PUB(subQuery sq, vector<vector<double>>& countRes); 

        pair<vector<vector<vector<double>>>, vector<double>> generate_db_columns(header &h);

        void evaluate_expected_count(query& q);

        void evaluate_expected_sum(query &q, int index);

        vector<long long> count(query &q);

        vector<double> sum(query &q, int index);

        void fill_db_TPCH(vector<int> indices);

        vector<vector<vector<double>>> generate_db_TPCH(header& h, vector<long long>& vals);

        void printReport(query& q);

        void combineSubQueriesENC(query& q);

        void combineSubQueriesPUB(query &q);

        void combineExpectedSubQueries(query &q);

        void combineExpectedAggQueries(int index, query& q);

    private:
        uint32_t num_rows;
        uint32_t ct_rows;
        vector<header> headers;
        unordered_map<string, vector<vector<Plaintext>>> plain_data;
        unordered_map<string, vector<double>> expected_data;
        unordered_map<string, vector<vector<vector<double>>>> public_data;
        bool pub_res;
        bool exp_res;

        vector<Ciphertext<DCRTPoly>> enc_pred_res;
        vector<vector<double>> pub_pred_res;
        vector<long long> exp_pred_res;
        vector<Ciphertext<DCRTPoly>> enc_agg_res;
        vector<vector<double>> pub_agg_res;
        vector<double> exp_agg_res;
};