#include "q.h"
#include "header.h"

subQuery::subQuery(const header &h, string op, long long cmp_val) {
    if (op == "<=" || op == ">=") {
        assert(cmp_val > h.low && cmp_val < h.high);
    }
    if (op == "<" || op == ">" || op == "==" || op == "!=") {
        assert(cmp_val >= h.low && cmp_val <= h.high);
    }
    this->op = op;
    this->h = h;
    this->ID = next_ID();
    this->cmp_val = cmp_val;
}

subQuery::subQuery(const header &h) {
    this->op = "";
    this->h = h;
    this->ID = next_ID();
}

subQuery::subQuery(double scalar) {
    this->op = "";
    this->scalar = scalar;
    this->ID = next_ID();
}

subQuery::subQuery(string op) {
    this->op = op;
}

subQuery::subQuery(const subQuery &sq) {
    this->op = sq.op;
    this->h = sq.h;
    this->ID = sq.ID;
    this->cmp_val = sq.cmp_val;
    this->scalar = sq.scalar;
}

query::query(vector<string> final_operations, vector<subQuery> encSubQueries) {
    assert(encSubQueries.size() > 0);
    this->final_operations = final_operations;
    this->encSubQueries = encSubQueries;
}

query::query(vector<string> final_operations, vector<subQuery> encSubQueries, vector<vector<subQuery>> encAggQueries) {
    assert(encSubQueries.size() > 0);
    this->final_operations = final_operations;
    this->encSubQueries = encSubQueries;
    this->encAggQueries = encAggQueries;
}

void query::mapSubQueryResultENC(string subQueryID, const vector<Ciphertext<DCRTPoly>>& res) {
    this->enc_rows.insert({subQueryID, res});
}

void query::mapSubQueryResultPUB(string subQueryID, const vector<vector<double>>& res) {
    this->pub_rows.insert({subQueryID, res});
}

void query::mapAggQueryResultENC(string subQueryID, const vector<Ciphertext<DCRTPoly>>& res) {
    this->agg_enc_rows.insert({subQueryID, res});
}

void query::mapAggQueryResultPUB(string subQueryID, const vector<vector<double>>& res) {
    this->agg_pub_rows.insert({subQueryID, res});
}

void query::mapSubQueryExpectedResult(string subQueryID, const vector<long long>& res) {
    this->expected_rows.insert({subQueryID, res});
}

void query::mapAggQueryExpectedResult(string subQueryID, const vector<double> res) {
    this->agg_expected_rows.insert({subQueryID, res});
}

vector<subQuery> query::getEncSubQueries() {
    return this->encSubQueries;
}

vector<vector<subQuery>> query::getEncAggQueries() {
    return this->encAggQueries;
}

vector<string> query::getFinalOperations() {
    return this->final_operations;
}

vector<Ciphertext<DCRTPoly>> query::getQueryResultENC() {
    return this->query_result_ENC;
}

vector<vector<double>> query::getQueryResultPUB() {
    return this->query_result_PUB;
}

vector<Ciphertext<DCRTPoly>> query::getAggQueryResultENC() {
    return this->agg_query_result_ENC;
}

vector<vector<double>> query::getAggQueryResultPUB() {
    return this->agg_query_result_PUB;
}

vector<long long> query::getExpectedQueryResult() {
    return this->expected_query_result;
}

vector<double> query::getExpectedAggQueryResult() {
    return this->expected_agg_query_result;
}

vector<Ciphertext<DCRTPoly>> query::getExpandedQuery() {
    return this->expanded_query;
}

unordered_map<string, pair<size_t, size_t>> query::getSubQueryQIndices() {
    return this->subQuery_q_indices;
}

Ciphertext<DCRTPoly> query::getCompressedQuery() {
    return this->compressed_query;
}

unordered_map<string, vector<Ciphertext<DCRTPoly>>> query::getEncRows() {
    return this->enc_rows;
}

unordered_map<string, vector<Ciphertext<DCRTPoly>>> query::getAggEncRows() {
    return this->agg_enc_rows;
}

unordered_map<string, vector<vector<double>>> query::getPubRows() {
    return this->pub_rows;
}

unordered_map<string, vector<vector<double>>> query::getAggPubRows() {
    return this->agg_pub_rows;
}

unordered_map<string, vector<long long>> query::getExpectedRows() {
    return this->expected_rows;
}

unordered_map<string, vector<double>> query::getAggExpectedRows() {
    return this->agg_expected_rows;
}

void query::setCompressedQuery(Ciphertext<DCRTPoly> compressed_query) {
    this->compressed_query = compressed_query;
}

void query::formulate_query() {
    vector<long long> q_values;
    for (auto &sq : this->encSubQueries) {
        if (std::find(ops.begin(), ops.end(), sq.op) != ops.end()) {
            auto max_val = sq.h.high;
            auto sq_val = sq.cmp_val;
            auto start = cur_q_index;
            while (max_val) {
                long long cur_sq_col_val = abs(sq_val) % numDim;
                if (sq_val < 0) cur_sq_col_val *= -1;
                sq_val /= numDim;
                max_val /= numDim;
                q_values.push_back(cur_sq_col_val);
                cur_q_index += 1;
            }
            auto end = cur_q_index;
            subQuery_q_indices.insert({sq.ID, make_pair(start, end)});
        }
    }

    num_vals = q_values.size();
    q_slots = 1 << (numSlotBits - int(ceil(log2(q_values.size()))));
    q_values.resize(numSlots / q_slots);
    vector<double> q_plain(numSlots);
    for (size_t i = 0; i < q_plain.size(); i++) {
        q_plain[i] = q_values[i / q_slots];
    }
    compressed_query = encrypt_vector(q_plain);

    return;
}

void expand_slots(vector<Ciphertext<DCRTPoly>>& expanded_query, size_t i, size_t q_slots,
 Ciphertext<DCRTPoly>& compressed_query) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    for (size_t i = start; i < end; i++) {

        vector<double> extract_slots(numSlots);
        for (size_t j = i * q_slots; j < (i + 1) * q_slots; j++) {
            extract_slots[j] = 1;
        }
        auto enc_slots = encrypt_vector(extract_slots);
        expanded_query[i] = (*cc)->EvalMult(enc_slots, compressed_query);
        if (rescaleTech == FIXEDMANUAL) expanded_query[i] = (*cc)->Rescale(expanded_query[i]);

        for (size_t l = log2(numSlots / q_slots); l > 0; l--) {
            int32_t step = 1 << (numSlotBits - l);
            Ciphertext<DCRTPoly> rot = (*cc)->EvalRotate(expanded_query[i], step);
            expanded_query[i] = (*cc)->EvalAdd(expanded_query[i], rot);
        }
    }
    return;
}

void query::query_expansion() {
    start_end_threads = vector<size_t>(min(TOTAL_MACHINE_THREAD, this->num_vals));
    for (size_t i = 0; i < this->num_vals; i++) {
        start_end_threads[i % TOTAL_MACHINE_THREAD] += 1;
    }
    start_end_threads.insert(start_end_threads.begin(), 0);
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        start_end_threads[i] += start_end_threads[i - 1];
    }

    expanded_query = vector<Ciphertext<DCRTPoly>>(num_vals);

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(expand_slots, ref(expanded_query), i, q_slots, ref(compressed_query));
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }

    return;
}

// cite from geeks for geeks algorithm example
vector<subQuery> query::infix_to_postfix(vector<subQuery> infix) {
    stack<subQuery> s;
    vector<subQuery> postfix;
    for (auto &sq : infix) {
        if (std::find(ops.begin(), ops.end(), sq.op) != ops.end() ||
         sq.op == "(" || sq.op == "") postfix.push_back(sq);

        else if (sq.op == ")") {
            while (s.top().op != "(") {
                postfix.push_back(s.top());
                s.pop();
            }
            s.pop();

        } else {
            while (!s.empty() && priority[sq.op] < priority[s.top().op]) {
                postfix.push_back(s.top());
                s.pop();
            }
            s.push(sq.op);
        }
    }
    while (!s.empty()) {
        postfix.push_back(s.top());
        s.pop();
    }
    return postfix;
}

void combineResultsENCHelper(vector<Ciphertext<DCRTPoly>>& cur, vector<Ciphertext<DCRTPoly>>& prev, string op, size_t i) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    if (op == "AND") {
        for (size_t i = start; i < end; i++) {
            cur[i] = (*cc)->EvalMult(prev[i], cur[i]);
            if (rescaleTech == FIXEDMANUAL) cur[i] = (*cc)->Rescale(cur[i]);
        }
    } else {
        for (size_t i = start; i < end; i++) {
            cur[i] = (*cc)->EvalSub(1, cur[i]);
            prev[i] = (*cc)->EvalSub(1, prev[i]);
            cur[i] = (*cc)->EvalMult(prev[i], cur[i]);
            if (rescaleTech == FIXEDMANUAL) cur[i] = (*cc)->Rescale(cur[i]);
            cur[i] = (*cc)->EvalSub(1, cur[i]);
        }
    }
}

void query::combineResultsENC(string prev_res_ID, string sq_ID, string op) {
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineResultsENCHelper, ref(enc_rows[sq_ID]),
         ref(enc_rows[prev_res_ID]), op, i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void PredCompsENCHelper(vector<Ciphertext<DCRTPoly>>& pred, size_t i) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    for (size_t i = start; i < end; i++) {
        pred_comp_ENC(p_comps, pred[i]);
    }
}

void query::PredCompsENC(vector<Ciphertext<DCRTPoly>>& pred) {
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(PredCompsENCHelper, ref(pred), i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void PredCompsPUBHelper(vector<vector<double>>& pred, size_t i) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    for (size_t i = start; i < end; i++) {
        pred_comp_PUB(p_comps, pred[i]);
    }
}

void query::PredCompsPUB(vector<vector<double>>& pred) {
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(PredCompsPUBHelper, ref(pred), i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void combineResultsPUBHelper(vector<vector<double>>& cur, vector<vector<double>>& prev, string op, size_t i) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    if (op == "AND") {
        for (size_t i = start; i < end; i++) {
            for (size_t j = 0; j < prev[i].size(); j++) {
                cur[i][j] *= prev[i][j];
            }
        }
    } else {
        for (size_t i = start; i < end; i++) {
            for (size_t j = 0; j < prev[i].size(); j++) {
                cur[i][j] = 1 - (1 - prev[i][j]) * (1 - cur[i][j]);
            }
        }
    }
}

void query::combineResultsPUB(string prev_res_ID, string sq_ID, string op) {

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineResultsPUBHelper, ref(pub_rows[sq_ID]),
         ref(pub_rows[prev_res_ID]), op, i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void combineExpResultsHelper(vector<long long>& cur, vector<long long>& prev, string op, size_t i) {

    if (i >= start_end_expected_threads.size()) return;
    size_t start = start_end_expected_threads[i - 1];
    size_t end = start_end_expected_threads[i];

    if (op == "AND") {
        for (size_t i = start; i < end; i++) {
            cur[i] = prev[i] && cur[i];
        }
        
    } else {
        for (size_t i = start; i < end; i++) {
            cur[i] = prev[i] || cur[i];
        }
    }
}

void query::combineExpResults(string prev_res_ID, string sq_ID, string op) {

    for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineExpResultsHelper, ref(expected_rows[sq_ID]),
         ref(expected_rows[prev_res_ID]), op, i);
    }

    for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void combineAggResultsENCHelper(vector<Ciphertext<DCRTPoly>>& pred, vector<double>& agg,
 vector<Ciphertext<DCRTPoly>>& res, size_t i) {

    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    for (size_t i = start; i < end; i++) {
        vector<double> plain_vec(numSlots, 0);
        for (size_t j = 0; j < numSlots && i * numSlots + j < agg.size(); j++) {
            plain_vec[j] = agg[i * numSlots + j];
        }
        Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(plain_vec);
        res[i] = (*cc)->EvalMult(pred[i], x_plain);
        if (rescaleTech == FIXEDMANUAL) res[i] = (*cc)->Rescale(res[i]);
    }
}

void query::combineAggResultsENC(vector<Ciphertext<DCRTPoly>>& pred, vector<double>& agg,
 vector<Ciphertext<DCRTPoly>>& res) {
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineAggResultsENCHelper, ref(pred),
        ref(agg), ref(res), i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void combineAggResultsPUBHelper(vector<vector<double>>& pred,
  vector<double>& agg, vector<vector<double>>& res, size_t i) {
    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];

    for (size_t i = start; i < end; i++) {
        for (size_t j = 0; j < pred[i].size() && i * numSlots + j < agg.size(); j++) {
            res[i][j] = pred[i][j] * agg[i * numSlots + j];
        }
    }

}

void query::combineAggResultsPUB(vector<vector<double>>& pred, vector<double>& agg, vector<vector<double>>& res) {
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineAggResultsPUBHelper, ref(pred), ref(agg), ref(res), i);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

void combineAggExpResultsHelper(vector<double>& cur, vector<double>& prev, string op, size_t i) {
    
    if (i >= start_end_expected_threads.size()) return;
    size_t start = start_end_expected_threads[i - 1];
    size_t end = start_end_expected_threads[i];

    if (op == "+") {
        for (size_t i = start; i < end; i++) {
            cur[i] += prev[i];
        }
    } else if (op == "-") {
        for (size_t i = start; i < end; i++) {
            cur[i] = prev[i] - cur[i];
        }
    } else if (op == "*") {
        for (size_t i = start; i < end; i++) {
            cur[i] *= prev[i];
        }
    } else {
        for (size_t i = start; i < end; i++) {
            cur[i] = prev[i] / cur[i];
        }
    }
}

void query::combineAggExpResults(string prev_res_ID, string sq_ID, string op) {

    for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(combineAggExpResultsHelper, ref(agg_expected_rows[sq_ID]),
         ref(agg_expected_rows[prev_res_ID]), op, i);
    }

    for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }
}

// calculator that estimates min mult depth required to perform operations
// It will not always be very accurate and give the best multiplicative depth because it 
// assumes one level is consumed per predicate subquery evaluation (which is not always the case)
// for the worst case scenario.
// calculator may also not work for 64 bit values due to long long restrictions
uint32_t query::calculate_min_depth(int num_rows) {

    uint32_t new_mult_depth = 0;
    uint32_t numBits = 0;
    auto n_cols = 1;
    auto pred_mults = 0;
    auto neg_flag = 0;
    auto equality_flag = 0;
    auto agg_mults = 0;
    auto multiple_rows_flag = 0;
    auto query_expansion_cost = 1;

    // if rows are not divisible by numSlots, need an extra multiplication to get rid of dud rows in last ciphertext
    auto uneven_rows_flag = 0;
    if (num_rows % numSlots != 0) {
        uneven_rows_flag = 1;
    }

    for (auto &s : encSubQueries) {
        if (std::find(ops.begin(), ops.end(), s.op) == ops.end()) continue;
        if (s.op == "==" || s.op == "!=") equality_flag = 2;
        if (s.scalar == 0.0)
            n_cols = max(n_cols, int(ceil(log2(max(abs(s.h.low), abs(s.h.high))) / log2(numDim))));
        if (abs(s.h.high - s.h.low) > abs(s.h.high)) neg_flag = 1;

        numBits = max(numBits, uint32_t(ceil(log2(s.h.high - s.h.low) + 1)));
        cout << s.h.high - s.h.low << endl;
        pred_mults++;
    }
    pred_mults -= 1;

    if (encAggQueries.size() > 0) agg_mults = 1;
    if (n_cols > 1) multiple_rows_flag = 1;

    cout << "NUMBITS: " << numBits << endl;

    numBits = min(numBits, numDimBits) + neg_flag;
    cout << "\nMIN DEPTH CALCULATOR STATS:" << endl;
    cout << "number of bits to accomodate largest range: " << numBits << endl;
    cout << "Max Header columns in query: " << n_cols << endl;
    cout << "Max multiplications seen in agg queries: " << agg_mults << endl;
    cout << "Number of predicate multiplications: " << pred_mults << endl;
    cout << "Number of predicate compositions: " << p_comps << endl;
    cout << "Equality flag status: " << bool(equality_flag) << endl;
    cout << "Uneven rows flag status: " << bool(uneven_rows_flag) << endl;
    cout << "Multiple rows flag status: " << bool(multiple_rows_flag) << endl << endl;
    auto pred_comps = comp_ref[numBits];
    auto col_comps = comp_ref[n_cols];
    if (n_cols == 1) col_comps = comp_ref[0];

    cout << pred_comps.first << " " << pred_comps.second << endl;
    cout << col_comps.first << " " << col_comps.second << endl;

    new_mult_depth = query_expansion_cost + (2 * (pred_comps.first + pred_comps.second) + 1) + 
                (2 * (col_comps.first + col_comps.second)) + (2 * p_comps) +
                pred_mults + agg_mults + equality_flag + uneven_rows_flag + multiple_rows_flag;

    return new_mult_depth + 1;
}