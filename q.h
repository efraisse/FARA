#include "globals.h"
#include "header.h"

#pragma once

class subQuery {
    public:
        subQuery(const header &h, string op, long long cmp_val);
        subQuery(const header &h);
        subQuery(const subQuery &sq);
        subQuery(double scalar);
        subQuery(string op);
        string ID;
        string op;
        header h;
        long long cmp_val;
        double scalar = 0.0;
};


class query {
    public:
        query(vector<string> final_operations, vector<subQuery> encSubQueries);

        query(vector<string> final_operations, vector<subQuery> encSubQueries, vector<vector<subQuery>> encAggQueries);

        void formulate_query();

        void query_expansion();

        void mapSubQueryResultENC(string subQueryID, const vector<Ciphertext<DCRTPoly>>& res);

        void mapSubQueryResultPUB(string subQueryID, const vector<vector<double>>& res);

        void mapAggQueryResultENC(string subQueryID, const vector<Ciphertext<DCRTPoly>>& res);

        void mapAggQueryResultPUB(string subQueryID, const vector<vector<double>>& res);

        void PredCompsENC(vector<Ciphertext<DCRTPoly>>& pred);

        void PredCompsPUB(vector<vector<double>>& pred);

        void combineSubQueriesENC();

        void combineSubQueriesPUB();

        void combineAggSubQueriesENC(int index);

        void combineAggSubQueriesPUB(int index);

        void combineExpectedSubQueries();

        void combineExpectedAggQueries(int index);

        void mapSubQueryExpectedResult(string subQueryID, const vector<long long>& res);

        void mapAggQueryExpectedResult(string subQueryID, const vector<double> res);

        vector<subQuery> getEncSubQueries();

        vector<vector<subQuery>> getEncAggQueries();

        vector<string> getFinalOperations();

        vector<Ciphertext<DCRTPoly>> getQueryResultENC();

        vector<vector<double>> getQueryResultPUB();

        vector<Ciphertext<DCRTPoly>> getAggQueryResultENC();

        vector<vector<double>> getAggQueryResultPUB();

        vector<long long> getExpectedQueryResult();

        vector<double> getExpectedAggQueryResult();

        vector<Ciphertext<DCRTPoly>> getExpandedQuery();

        unordered_map<string, pair<size_t, size_t>> getSubQueryQIndices();

        Ciphertext<DCRTPoly> getCompressedQuery();

        unordered_map<string, vector<Ciphertext<DCRTPoly>>> getEncRows();

        unordered_map<string, vector<Ciphertext<DCRTPoly>>> getAggEncRows();

        unordered_map<string, vector<vector<double>>> getPubRows();

        unordered_map<string, vector<vector<double>>> getAggPubRows();

        unordered_map<string, vector<long long>> getExpectedRows();

        unordered_map<string, vector<double>> getAggExpectedRows();

        void setCompressedQuery(Ciphertext<DCRTPoly> compressed_query);

        void combineResultsENC(string prev_res_ID, string sq_ID, string op);

        void combineResultsPUB(string prev_res_ID, string sq_ID, string op);

        void combineExpResults(string prev_res_ID, string sq_ID, string op);

        void combineAggResultsENC(vector<Ciphertext<DCRTPoly>>& pred, vector<double>& agg, vector<Ciphertext<DCRTPoly>>& res);

        void combineAggResultsPUB(vector<vector<double>>& pred, vector<double>& agg, vector<vector<double>>& res);

        void combineAggExpResults(string prev_res_ID, string sq_ID, string op);

        uint32_t calculate_min_depth(int num_rows);

        vector<subQuery> infix_to_postfix(vector<subQuery> infix);

    private:
        long long numRows;
        vector<string> final_operations;
        vector<subQuery> encSubQueries;
        vector<vector<subQuery>> encAggQueries;
        vector<Ciphertext<DCRTPoly>> query_result_ENC;
        vector<vector<double>> query_result_PUB;
        unordered_map<string, vector<Ciphertext<DCRTPoly>>> enc_rows;
        unordered_map<string, vector<vector<double>>> pub_rows;
        unordered_map<string, vector<long long>> expected_rows;
        vector<long long> expected_query_result;
        vector<Ciphertext<DCRTPoly>> agg_query_result_ENC;
        vector<vector<double>> agg_query_result_PUB;
        unordered_map<string, vector<Ciphertext<DCRTPoly>>> agg_enc_rows;
        unordered_map<string, vector<vector<double>>> agg_pub_rows;
        unordered_map<string, vector<double>> agg_expected_rows;
        vector<double> expected_agg_query_result;

        // variables related to query expansion
        Ciphertext<DCRTPoly> compressed_query;
        vector<Ciphertext<DCRTPoly>> expanded_query;
        unordered_map<string, pair<size_t, size_t>> subQuery_q_indices;
        size_t q_slots;
        size_t num_vals;
};