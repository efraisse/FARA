#include "db.h"

db::db(uint32_t num_rows, vector<header> headers, bool pub_res, bool exp_res) {
    this->num_rows = num_rows;
    this->ct_rows = (uint32_t)(ceil(this->num_rows / (double)(numSlots)));
    this->headers = headers;
    this->pub_res = pub_res;
    this->exp_res = exp_res;
}

pair<vector<vector<vector<double>>>, vector<double>> db::generate_db_columns(header& h) {

    vector<double> expected_values;
    vector<vector<vector<double>>> vector_db;
    for (size_t i = 0; i < this->ct_rows; i++) {
        vector<vector<double>> row;
        vector<long long> nums;
        auto min_val_placeholder = h.low;
        if (h.col_name.find("date") != std::string::npos) {
            for (size_t n = 0; n < numSlots; n++) {
                auto date_low = daysSinceJan1st1970(h.low);
                auto date_high = daysSinceJan1st1970(h.high);
                long long random_val = (mt() % (date_high - date_low + 1)) + date_low;
                nums.push_back(random_val);
                expected_values.push_back((double)(random_val) * (1.0 / pow(10, h.precision)));
            }
        } else {
            for (size_t n = 0; n < numSlots; n++) {
                long long random_val = (mt() % (h.high - h.low + 1)) + h.low;
                nums.push_back(random_val);
                expected_values.push_back((double)(random_val) * (1.0 / pow(10, h.precision)));
            }
        }
        auto temp = max(abs(h.high), abs(h.low));
        auto f = 1.0;
        while (temp >>= 1) f += 1;
        h.num_ct = (int)(ceil(f / numDimBits));
        for (size_t j = 0; j < h.num_ct; j++) {
            vector<double> ct;
            for (size_t k = 0; k < numSlots; k++) {
                if (i * numSlots + k < this->num_rows) {
                    int sign = nums[k] > 0 ? 1 : -1;
                    ct.push_back(sign * (abs(nums[k]) % numDim));
                    nums[k] /= numDim;
                } else {
                    int sign = min_val_placeholder > 0 ? 1 : -1;
                    ct.push_back(sign * (abs(min_val_placeholder) % numDim));
                }
            }
            min_val_placeholder /= numDim;
            row.insert(row.begin(), ct);
        }
        vector_db.push_back(row);
    }
    cout << "generated db columns for header " << h.col_name << endl;
    return pair(vector_db, expected_values);
}

vector<vector<vector<double>>> db::generate_db_TPCH(header& h, vector<long long>& vals) {
    vector<vector<vector<double>>> vector_db;
    for (size_t i = 0; i < this->ct_rows; i++) {
        vector<vector<double>> row;
        auto temp = max(abs(h.high), abs(h.low));
        auto min_val_placeholder = h.low;
        auto f = 1.0;
        while (temp >>= 1) f += 1;
        h.num_ct = (int)(ceil(f / numDimBits));
        for (size_t j = 0; j < h.num_ct; j++) {
            vector<double> ct;
            for (size_t k = 0; k < numSlots; k++) {
                if (i * numSlots + k < this->num_rows) {
                    int sign = vals[i * numSlots + k] > 0 ? 1 : -1;
                    ct.push_back(sign * (abs(vals[i * numSlots + k]) % numDim));
                    vals[i * numSlots + k] /= numDim;
                } else {
                    int sign = min_val_placeholder > 0 ? 1 : -1;
                    ct.push_back(sign * (abs(min_val_placeholder) % numDim));
                }
            }
            min_val_placeholder /= numDim;
            row.insert(row.begin(), ct);
        }
        vector_db.push_back(row);
    }
    cout << "generated db columns for header " << h.col_name << endl;
    return vector_db;
}

pair<vector<vector<long long>>, vector<vector<double>>> grab_TPCH_data(int num_rows, vector<int> indices, vector<header> &headers)
{ 
	ifstream file(tableFilePath + "lineitem_10GB.tbl");
	string line;
    vector<vector<long long>> TPCH_data(indices.size(), vector<long long>(num_rows, 0));
    vector<vector<double>> expected_data(indices.size(), vector<double>(num_rows, 0.0));
    int index = 0;
	if (file.is_open()) { 
		while (getline(file, line) && index < num_rows) { 
            vector<string> res = split(line, "|");
            for (size_t i = 0; i < indices.size(); i++) {
                if (headers[i].col_name.find("date") != std::string::npos) {
                    vector<string> date = split(res[indices[i]], "-");
                    int format_date = stoi(date[2]) + stoi(date[1]) * 100 + stoi(date[0]) * 10000;
                    TPCH_data[i][index] = daysSinceJan1st1970(format_date);
                    expected_data[i][index] = daysSinceJan1st1970(format_date);
                } else {
                    TPCH_data[i][index] = (long long)(stod(res[indices[i]]) * (pow(10, headers[i].precision)));
                    expected_data[i][index] = stod(res[indices[i]]);
                }
            }
            index += 1;
		}

		file.close();
	} 
	else {
		cerr << "Error while processing tpc-h data" << endl; 
	}

    return make_pair(TPCH_data, expected_data);
}

void db::fill_db_TPCH(vector<int> indices) {
    auto data = grab_TPCH_data(this->num_rows, indices, headers);
    int index = 0;
    for (auto &h : headers) {
        auto col = this->generate_db_TPCH(h, data.first[index]);
        vector<vector<Plaintext>> plain_col;
        for (auto &v1: col){
            vector<Plaintext> plain_row;
            for(auto &v2: v1) {
                Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(v2);
                plain_row.push_back(x_plain);
            }
            plain_col.push_back(plain_row);
        }
        this->plain_data.insert({h.col_name, plain_col});
        if (pub_res) this->public_data.insert({h.col_name, col});
        this->expected_data.insert({h.col_name, data.second[index]});
        index += 1;
    }
    cout << "Successfully filled the db with TPC-H values " << endl;
}

void db::fill_db_random_values() {
    for (auto &h : headers) {
        auto cols = this->generate_db_columns(h);
        vector<vector<vector<double>>> col = cols.first;
        vector<double> expected_cols = cols.second;
        vector<vector<Plaintext>> plain_col;
        for (auto &v1: col){
            vector<Plaintext> plain_row;
            for(auto &v2: v1) {
                Plaintext x_plain = (*cc)->MakeCKKSPackedPlaintext(v2);
                plain_row.push_back(x_plain);
            }
            plain_col.push_back(plain_row);
        }
        this->plain_data.insert({h.col_name, plain_col});
        if (pub_res) this->public_data.insert({h.col_name, col});
        this->expected_data.insert({h.col_name, expected_cols});
    }
    cout << "Successfully filled the db with random values " << endl;
}

void db::print_initial_db() {
    if (!pub_res) {
        cout << "public rows were not stored" << endl;
        return;
    }
    cout.precision(8);
    cout << "printing the contents of the initial db: " << endl; 
    for (auto &h : this->headers) {
        cout << "header " << h.col_name << "'s rows and columns:" << endl;
        for (size_t i = 0; i < this->public_data[h.col_name].size(); i++) {
            for (size_t j = 0; j < this->public_data[h.col_name][i].size(); j++) {
                cout << "Row: " << to_string(i) << ": " << "Col " << to_string(j) << ": " << endl;
                print_vector(this->public_data[h.col_name][i][j], PRINT_SIDES, PRINT_TOTAL);
            }
        }
    }
}

void db::print_expected_db() {
    cout.precision(8);
    cout << "printing the contents of the expected db: " << endl; 
    for (auto &h : this->headers) {
        cout << "header " << h.col_name << "'s rows and columns:" << endl;
        for (size_t i = 0; i < this->expected_data[h.col_name].size(); i++) {
            print_vector(this->expected_data[h.col_name], 4, 9);
        }
    }
}

void db::print_db() {
    cout.precision(8);
    cout << "printing the contents of the encrypted db: " << endl; 
    for (auto &h : this->headers) {
        cout << "header " << h.col_name << "'s rows and columns:" << endl;
        for (size_t i = 0; i < this->plain_data[h.col_name].size(); i++) {
            for (size_t j = 0; j < this->plain_data[h.col_name][i].size(); j++) {
                cout << "Row: " << to_string(i) << ": " << "Col " << to_string(j) << ": " << endl;
                vector<double> realNum;
                for (auto el: this->plain_data[h.col_name][i][j]->GetCKKSPackedValue()) {
                    realNum.push_back(static_cast<double>(el.real()));
                }
                print_vector(realNum, PRINT_SIDES, PRINT_TOTAL);
            }
        }
    }
}

// NOTE: depth calculator below is only tested on tpc-h q6 to make it run faster
// when using the FIXEDMANUAL rescaling technique.
// calculator is not the most reliable, it is always better to
// observe and then fill in the minimum multiplicative depth
// The reason why is due to the predicate multiplications
// the calculator assumes that one multiplicative depth will be consumed per
// predicate subquery operation but that is not always the case
uint32_t min_depth_requirements_TPCHQ6(query& q, subQuery& sq, bool verbose = false) {

    uint32_t new_mult_depth = 0;
    uint32_t numBits = 0;
    auto n_cols = 1;
    auto pred_mults = 0;
    auto neg_flag = 0;
    auto equality_flag = 0;
    auto agg_mults = 0;
    auto multiple_rows_flag = 0;

    // if rows are not divisible by numSlots, need an extra multiplication to get rid of dud rows in last ciphertext
    auto uneven_rows_flag = 0;
    if (total_rows % numSlots != 0) {
        uneven_rows_flag = 1;
    }

    for (auto &s : q.getEncSubQueries()) {
        if (std::find(ops.begin(), ops.end(), s.op) == ops.end()) continue;

        pred_mults++;
    }
    numBits = max(numBits, uint32_t(ceil(log2(abs(sq.h.high - sq.h.low) + 1))));
    n_cols = max(n_cols, int(ceil(log2(max(abs(sq.h.low), abs(sq.h.high))) / log2(numDim))));
    if (sq.op == "==" || sq.op == "!=") equality_flag = 2;
    if (abs(sq.h.high - sq.h.low) > abs(sq.h.high)) neg_flag = 1;
    pred_mults -= 1;

    if (q.getEncAggQueries().size() > 0) agg_mults = 1;
    if (n_cols > 1) multiple_rows_flag = 1;

    numBits = min(numBits, numDimBits) + 2 + neg_flag;
    if (verbose) {
        cout << "number of bits to accomodate largest range: " << numBits << endl;
        cout << "Max Header columns in query: " << n_cols << endl;
        cout << "Max multiplications seen in agg queries: " << agg_mults << endl;
        cout << "Number of predicate multiplications: " << pred_mults << endl;
        cout << "Number of predicate compositions: " << p_comps << endl;
        cout << "Equality flag status: " << bool(equality_flag) << endl;
        cout << "Uneven rows flag status: " << bool(uneven_rows_flag) << endl;
        cout << "Multiple rows flag status: " << bool(multiple_rows_flag) << endl;
    }
    auto pred_comps = comp_ref[numBits];
    auto col_comps = comp_ref[n_cols];
    if (n_cols == 1) col_comps = comp_ref[0];

    new_mult_depth = 2 * (pred_comps.first + pred_comps.second) +
                2 * (col_comps.first + col_comps.second) + (2 * p_comps) + multiple_rows_flag +
                (pred_mults - tpchq6_pred_mults_reduction) + agg_mults + equality_flag + uneven_rows_flag;
    
    new_mult_depth = min(new_mult_depth, multDepth - 1);
    if (verbose) cout << "NEW DEPTH FOR THIS COLUMN: " << new_mult_depth << endl;

    return new_mult_depth;
}

void greater_than_equal_ENC(vector<Ciphertext<DCRTPoly>>& h_enc_data, query& q,
 subQuery& sq, long long total_h_low, long long total_h_high, bool transform) {
    
    auto mult_cols = 1;
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    auto ciph_index = q.getSubQueryQIndices()[sq.ID].first;

    for (size_t j = h_enc_data.size(); j --> 0;) {

        if (j != h_enc_data.size() - 1) {
            mult_cols = 0;
        }

        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                if (h_enc_data.size() == 1) {
                    vector<double> vec_copy(numSlots, 1.0);
                    h_enc_data[j] = encrypt_vector(vec_copy);
                } else {
                    vector<double> vec_copy(numSlots, 0.5);
                    h_enc_data[j] = encrypt_vector(vec_copy);
                }
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }
        
        Ciphertext<DCRTPoly> ct_low = (*cc)->EvalSub(q.getExpandedQuery()[ciph_index], 0.5 * mult_cols);
        if (mult_cols && sq.op == ">") ct_low = (*cc)->EvalAdd(ct_low, 1);
        h_enc_data[j] = (*cc)->EvalSub(h_enc_data[j], ct_low);
        h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 1 / (double)(chosen_dim));
        if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);

        if (depth_optimization && rescaleTech == FIXEDMANUAL) {
            auto new_depth = min_depth_requirements_TPCHQ6(q, sq);
            h_enc_data[j] = (*cc)->Compress(h_enc_data[j], new_depth);
        }

        if (chosen_dim + (mult_cols && sq.op == "<") > 1) {
            comp_ENC(cur_comps.first, cur_comps.second, h_enc_data[j], transform);
        } else {
            if (transform) {
                h_enc_data[j] = (*cc)->EvalAdd(h_enc_data[j], 1);
                h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 0.5);
                if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
        ciph_index++;
    }
}

void greater_than_equal_PUB(vector<vector<double>>& h_pub_data, long long total_q_range_low,
        long long total_h_low, long long total_h_high, bool transform) {
    
    auto mult_cols = 1;
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    for (size_t j = h_pub_data.size(); j --> 0;) {

        int cur_q_range_low = abs(total_q_range_low) % numDim;
        if (total_q_range_low < 0) cur_q_range_low *= -1;
        total_q_range_low /= numDim;

        if (j != h_pub_data.size() - 1) {
            mult_cols = 0;
        }

        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                if (h_pub_data.size() == 1) {
                    for (size_t k = 0; k < h_pub_data[j].size(); k++)
                        h_pub_data[j][k] = 1;
                    vector<double> vec_copy(numSlots, 1.0);
                } else {
                    for (size_t k = 0; k < h_pub_data[j].size(); k++)
                        h_pub_data[j][k] = 0.5;
                    vector<double> vec_copy(numSlots, 0.5);
                }
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }

        for (size_t k = 0; k < h_pub_data[j].size(); k++)
            h_pub_data[j][k] = ((2 * h_pub_data[j][k]) - (2 * cur_q_range_low - mult_cols)) / (double)(chosen_dim * 2);
        

        if (chosen_dim > 1)
            comp_PUB(cur_comps.first, cur_comps.second, h_pub_data[j], transform);
        else {
            if (transform) {
                for (uint64_t i = 0; i < h_pub_data[j].size(); i++) h_pub_data[j][i] = (h_pub_data[j][i] + 1) / 2;
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
    }
}

void less_than_equal_ENC(vector<Ciphertext<DCRTPoly>>& h_enc_data, query& q,
 subQuery& sq, long long total_h_low, long long total_h_high, bool transform) {
    
    auto mult_cols = 1;
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    auto ciph_index = q.getSubQueryQIndices()[sq.ID].first;

    for (size_t j = h_enc_data.size(); j --> 0;) {

        if (j != h_enc_data.size() - 1) {
            mult_cols = 0;
        }

        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                if (h_enc_data.size() == 1) {
                    vector<double> vec_copy(numSlots, 1.0);
                    h_enc_data[j] = encrypt_vector(vec_copy);
                } else {
                    vector<double> vec_copy(numSlots, 0.5);
                    h_enc_data[j] = encrypt_vector(vec_copy);
                }
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }

        Ciphertext<DCRTPoly> ct_high = (*cc)->EvalAdd(q.getExpandedQuery()[ciph_index], 0.5 * mult_cols);
        if (mult_cols && sq.op == "<") ct_high = (*cc)->EvalSub(ct_high, 1);
        h_enc_data[j] = (*cc)->EvalSub(ct_high, h_enc_data[j]);
        h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 1 / (double)(chosen_dim));
        if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);

        if (depth_optimization && rescaleTech == FIXEDMANUAL) {
            auto new_depth = min_depth_requirements_TPCHQ6(q, sq);
            h_enc_data[j] = (*cc)->Compress(h_enc_data[j], new_depth);
        }

        if (chosen_dim + (mult_cols && sq.op == "<") > 1)
            comp_ENC(cur_comps.first, cur_comps.second, h_enc_data[j], transform);
        else {
            if (transform) {
                h_enc_data[j] = (*cc)->EvalAdd(h_enc_data[j], 1);
                h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 0.5);
                if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
        ciph_index++;
    }
}

void less_than_equal_PUB(vector<vector<double>>& h_pub_data, long long total_q_range_high,
        long long total_h_low, long long total_h_high, bool transform) {
    
    auto mult_cols = 1;
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    for (size_t j = h_pub_data.size(); j --> 0;) {

        int cur_q_range_high = abs(total_q_range_high) % numDim;
        if (total_q_range_high < 0) cur_q_range_high *= -1;
        total_q_range_high /= numDim;

        if (j != h_pub_data.size() - 1) {
            mult_cols = 0;
        }

        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                if (h_pub_data.size() == 1) {
                    for (size_t k = 0; k < h_pub_data[j].size(); k++)
                        h_pub_data[j][k] = 1;
                    vector<double> vec_copy(numSlots, 1.0);
                } else {
                    for (size_t k = 0; k < h_pub_data[j].size(); k++)
                        h_pub_data[j][k] = 0.5;
                    vector<double> vec_copy(numSlots, 0.5);
                }
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }

        for (size_t k = 0; k < h_pub_data[j].size(); k++)
            h_pub_data[j][k] = ((2 * cur_q_range_high + mult_cols) - (2 * h_pub_data[j][k])) / (double)(chosen_dim * 2);
        
        if (chosen_dim > 1)
            comp_PUB(cur_comps.first, cur_comps.second, h_pub_data[j], transform);
        else {
            if (transform) {
                for (uint64_t i = 0; i < h_pub_data[j].size(); i++) h_pub_data[j][i] = (h_pub_data[j][i] + 1) / 2;
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
    }
}

// same as greater than equal and less than equal but with no equality check.
// POTENTIALLY JUST DO A FLOOR ON CHOSEN_DIM TO REDUCE IT BY ONE
void equal_ENC(vector<Ciphertext<DCRTPoly>>& h_enc_data, query& q,
 subQuery& sq, long long total_h_low, long long total_h_high, bool transform) {
            
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    auto ciph_index = q.getSubQueryQIndices()[sq.ID].first;

    for (size_t j = h_enc_data.size(); j --> 0;) {

        // do not have to worry about equality check operations
        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                vector<double> vec_copy(numSlots, 0.5);
                h_enc_data[j] = encrypt_vector(vec_copy);
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }

        h_enc_data[j] = (*cc)->EvalSub(h_enc_data[j], q.getExpandedQuery()[ciph_index]);
        h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 1 / (double)(chosen_dim));
        if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);

        if (depth_optimization && rescaleTech == FIXEDMANUAL) {
            auto new_depth = min_depth_requirements_TPCHQ6(q, sq);
            h_enc_data[j] = (*cc)->Compress(h_enc_data[j], new_depth);
        }

        if (chosen_dim > 1)
            comp_ENC(cur_comps.first, cur_comps.second, h_enc_data[j], transform);
        else {
            if (transform) {
                h_enc_data[j] = (*cc)->EvalAdd(h_enc_data[j], 1);
                h_enc_data[j] = (*cc)->EvalMult(h_enc_data[j], 0.5);
                if (rescaleTech == FIXEDMANUAL) h_enc_data[j] = (*cc)->Rescale(h_enc_data[j]);
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
        ciph_index++;
    }
}

void equal_PUB(vector<vector<double>>& h_pub_data, long long total_q_range_low, 
 long long total_h_low, long long total_h_high, bool transform) {
            
    auto chosen_dim = numDim;
    if (total_h_low < 0 && total_h_high > 0) chosen_dim *= 2;
    auto cur_comps = comp_ref[numDimBits];
    if (total_h_low < 0 && total_h_high > 0) cur_comps = comp_ref[numDimBits + 1];
    for (size_t j = h_pub_data.size(); j --> 0;) {

        int cur_q_range_low = abs(total_q_range_low) % numDim;
        if (total_q_range_low < 0) cur_q_range_low *= -1;
        total_q_range_low /= numDim;

        if (j == 0) {
            chosen_dim = abs(total_h_high - total_h_low) + 1;
            if (chosen_dim == 1) {
                for (size_t k = 0; k < h_pub_data[j].size(); k++)
                    h_pub_data[j][k] = 0.5;
                vector<double> vec_copy(numSlots, 0.5);
                break;
            }
            cur_comps = comp_ref[int(ceil(log2(chosen_dim)))];
        }

        for (size_t k = 0; k < h_pub_data[j].size(); k++)
            h_pub_data[j][k] = (h_pub_data[j][k] - cur_q_range_low) / (double)(chosen_dim);
        
        if (chosen_dim > 1)
            comp_PUB(cur_comps.first, cur_comps.second, h_pub_data[j], transform);
        else {
            if (transform) {
                for (uint64_t i = 0; i < h_pub_data[j].size(); i++) h_pub_data[j][i] = (h_pub_data[j][i] + 1) / 2;
            }
        }

        total_h_low /= numDim;
        total_h_high /= numDim;
    }
}

void combine_multiple_columns_ENC(vector<Ciphertext<DCRTPoly>>& enc_col) {
    
    // aggregate all results in the last column
    size_t agg_i = enc_col.size() - 1;
    for (size_t i = enc_col.size(); i --> 0;) {
        auto res_col = (*cc)->EvalMult(enc_col[i], pow(2, -((int)(i) + 1)));
        if (rescaleTech == FIXEDMANUAL) res_col = (*cc)->Rescale(res_col);

        if (i == agg_i) enc_col[agg_i] = res_col;
        else { 
            enc_col[agg_i] = (*cc)->EvalAdd(enc_col[agg_i], res_col);
        }
    }

    auto factors = comp_ref[enc_col.size()];
    comp_ENC(factors.first, factors.second, enc_col[agg_i], true);
}

void combine_multiple_columns_PUB(vector<vector<double>>& pub_col) {
    
    // aggregate all results in the last column
    size_t agg_i = pub_col.size() - 1;

    for (size_t i = pub_col.size(); i --> 0;) {
        for (size_t j = 0; j < pub_col[i].size(); j++) {
            auto res_col = pub_col[i][j] * (double)(pow(2, -((int)(i) + 1)));
            if (i == agg_i) pub_col[agg_i][j] = res_col;
            else pub_col[agg_i][j] += res_col;
        }
    }

    auto factors = comp_ref[pub_col.size()];
    comp_PUB(factors.first, factors.second, pub_col[agg_i], true);
}

void evaluate_sub_query_helper_ENC(vector<vector<Ciphertext<DCRTPoly>>>& header_enc_data, query& q,
 subQuery& sq, size_t i, uint32_t num_rows) {
    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];
    for (size_t i = start; i < end; i++) {
        auto transform = true;
        long long total_h_low = sq.h.low;
        long long total_h_high = sq.h.high;

        if (header_enc_data[i].size() > 1) {
            transform = false;
        }

        if (sq.op == ">=") 
            greater_than_equal_ENC(header_enc_data[i], q, sq, total_h_low, total_h_high, transform);

        else if (sq.op == "<=")
            less_than_equal_ENC(header_enc_data[i], q, sq, total_h_low, total_h_high, transform);

        else if (sq.op == ">")
            greater_than_equal_ENC(header_enc_data[i], q, sq, total_h_low, total_h_high, transform);

        else if (sq.op == "<")
            less_than_equal_ENC(header_enc_data[i], q, sq, total_h_low, total_h_high, transform);

        else
            equal_ENC(header_enc_data[i], q, sq, total_h_low, total_h_high, transform);

        if (header_enc_data[i].size() > 1) {
            combine_multiple_columns_ENC(header_enc_data[i]);
        }

        auto last_index = header_enc_data[i].size() - 1;

        if (sq.op == "==") {
            auto enc_comp = (*cc)->EvalSub(1, header_enc_data[i][last_index]);
            enc_comp = (*cc)->EvalMult(enc_comp, 4);
            if (rescaleTech == FIXEDMANUAL) enc_comp = (*cc)->Rescale(enc_comp);
            header_enc_data[i][last_index] = (*cc)->EvalMult(header_enc_data[i][last_index], enc_comp);
            if (rescaleTech == FIXEDMANUAL) header_enc_data[i][last_index] = (*cc)->Rescale(header_enc_data[i][last_index]);
            header_enc_data[i][last_index] = (*cc)->EvalMult(header_enc_data[i][last_index], header_enc_data[i][last_index]);
            if (rescaleTech == FIXEDMANUAL) header_enc_data[i][last_index] = (*cc)->Rescale(header_enc_data[i][last_index]);

        } else if (sq.op == "!=") {
            auto enc_comp = (*cc)->EvalSub(1, header_enc_data[i][last_index]);
            enc_comp = (*cc)->EvalMult(enc_comp, 4);
            if (rescaleTech == FIXEDMANUAL) enc_comp = (*cc)->Rescale(enc_comp);
            header_enc_data[i][last_index] = (*cc)->EvalMult(header_enc_data[i][last_index], enc_comp);
            if (rescaleTech == FIXEDMANUAL) header_enc_data[i][last_index] = (*cc)->Rescale(header_enc_data[i][last_index]);
            header_enc_data[i][last_index] = (*cc)->EvalMult(header_enc_data[i][last_index], header_enc_data[i][last_index]);
            if (rescaleTech == FIXEDMANUAL) header_enc_data[i][last_index] = (*cc)->Rescale(header_enc_data[i][last_index]);
            header_enc_data[i][last_index] = (*cc)->EvalSub(1, header_enc_data[i][last_index]);
        }

        // if slot count does not divide the number of rows evenly
        // then we must disregard the extra elements in the last ciphertext
        // Default value to minimum (i.e. header.low) so as to avoid
        // polynomial composition errors if the value exists outside of bounds
        if (i == header_enc_data.size() - 1 && num_rows % numSlots != 0) {
            vector<double> omit_last_els(numSlots, 0);
            for (size_t j = 0; j < (num_rows % numSlots); j++) {
                omit_last_els[j] = 1;
            }
            auto enc_omit = encrypt_vector(omit_last_els);
            header_enc_data[i][last_index] = (*cc)->EvalMult(header_enc_data[i][last_index], enc_omit);
            if (rescaleTech == FIXEDMANUAL) header_enc_data[i][last_index] = (*cc)->Rescale(header_enc_data[i][last_index]);
        }
    }
}

vector<Ciphertext<DCRTPoly>> db::evaluate_sub_query_ENC(subQuery& sq, query& q) {
    vector<vector<Ciphertext<DCRTPoly>>> header_enc_data;

    for (size_t i = 0; i < this->plain_data[sq.h.col_name].size(); i++) {
        vector<Ciphertext<DCRTPoly>> header_enc_row;
        for (size_t j = 0; j < this->plain_data[sq.h.col_name][i].size(); j++) {
            Ciphertext<DCRTPoly> enc = (*cc)->Encrypt(keys->publicKey, this->plain_data[sq.h.col_name][i][j]);
            header_enc_row.push_back(enc);
        }
        header_enc_data.push_back(header_enc_row);
    }

    vector<Ciphertext<DCRTPoly>> enc_response;

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(evaluate_sub_query_helper_ENC, ref(header_enc_data), ref(q), ref(sq), i, num_rows);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }

    for (size_t i = 0; i < header_enc_data.size(); i++) {
        enc_response.push_back(header_enc_data[i][header_enc_data[i].size() - 1]);
    }

    return enc_response;
}

void evaluate_sub_query_helper_PUB(vector<vector<vector<double>>>& header_pub_data, subQuery& sq, size_t i, uint32_t num_rows) {
    if (i >= start_end_threads.size()) return;
    size_t start = start_end_threads[i - 1];
    size_t end = start_end_threads[i];
    for (size_t i = start; i < end; i++) {
        auto transform = true;
        long long cmp_val = sq.cmp_val;
        long long total_h_low = sq.h.low;
        long long total_h_high = sq.h.high;

        if (header_pub_data[i].size() > 1) {
            transform = false;
        }

        if (sq.op == ">=") 
            greater_than_equal_PUB(header_pub_data[i], cmp_val, total_h_low, total_h_high, transform);

        else if (sq.op == "<=")
            less_than_equal_PUB(header_pub_data[i], cmp_val, total_h_low, total_h_high, transform);

        else if (sq.op == ">")
            greater_than_equal_PUB(header_pub_data[i], cmp_val + 1, total_h_low, total_h_high, transform);

        else if (sq.op == "<")
            less_than_equal_PUB(header_pub_data[i], cmp_val - 1, total_h_low, total_h_high, transform);

        else
            equal_PUB(header_pub_data[i], cmp_val, total_h_low, total_h_high, transform);

        if (header_pub_data[i].size() > 1) {
            combine_multiple_columns_PUB(header_pub_data[i]);
        }

        auto last_index = header_pub_data[i].size() - 1;

        if (sq.op == "==") {

            for (size_t j = 0; j < header_pub_data[i][last_index].size(); j++) {
                auto pub_comp = 1 - header_pub_data[i][last_index][j];
                header_pub_data[i][last_index][j] = header_pub_data[i][last_index][j] * pub_comp;
                header_pub_data[i][last_index][j] = header_pub_data[i][last_index][j] * 4;
                header_pub_data[i][last_index][j] *= header_pub_data[i][last_index][j];
            }

        } else if (sq.op == "!=") {

            for (size_t j = 0; j < header_pub_data[i][last_index].size(); j++) {
                auto pub_comp = 1 - header_pub_data[i][last_index][j];
                header_pub_data[i][last_index][j] = header_pub_data[i][last_index][j] * pub_comp;
                header_pub_data[i][last_index][j] = header_pub_data[i][last_index][j] * 4;
                header_pub_data[i][last_index][j] *= header_pub_data[i][last_index][j];
                header_pub_data[i][last_index][j] = 1 - header_pub_data[i][last_index][j];
            }
        }

        // if slot count does not divide the number of rows evenly
        // then we must disregard the extra elements in the last ciphertext
        // Default value to minimum (i.e. header.low) so as to avoid
        // polynomial composition errors if value exists outside of bounds
        if (i == header_pub_data.size() - 1 && num_rows % numSlots != 0) {
            vector<double> omit_last_els(numSlots, 0);
            for (size_t j = 0; j < (num_rows % numSlots); j++) {
                omit_last_els[j] = 1;
            }
            for (size_t j = 0; j < omit_last_els.size(); j++) {
                header_pub_data[i][last_index][j] *= omit_last_els[j];
            }
        }
    }
}

vector<vector<double>> db::evaluate_sub_query_PUB(subQuery sq) {
    vector<vector<vector<double>>> header_pub_data = this->public_data[sq.h.col_name];
    vector<vector<double>> pub_response;

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1] = thread(evaluate_sub_query_helper_PUB, ref(header_pub_data), ref(sq), i, num_rows);
    }

    for (size_t i = 1; i < start_end_threads.size(); i++) {
        MACHINE_THREADS[i - 1].join();
    }

    for (size_t i = 0; i < header_pub_data.size(); i++) {
        pub_response.push_back(header_pub_data[i][header_pub_data[i].size() - 1]);
    }

    return pub_response;
}

vector<long long> db::count(query& q) {
    auto count_eval_time = chrono::high_resolution_clock::now();
    long long encrypted_result = 0;
    double real_encrypted_result = 0;

    auto res_enc = enc_pred_res[0];

    for (size_t i = 1; i < enc_pred_res.size(); i++) {
        res_enc = (*cc)->EvalAdd(res_enc, enc_pred_res[i]);
    }

    Plaintext result;
    (*cc)->Decrypt(keys->secretKey, res_enc, &result);
    result->SetLength(numSlots);
    vector<double> realNum;
    for (auto el: result->GetCKKSPackedValue()) {
        encrypted_result += (long long)(round(el.real()));
        real_encrypted_result += el.real();
    }

    count_time = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - count_eval_time)).count();

    query_results += "Remaining levels after COUNT operation: " + to_string(multDepth - res_enc->GetLevel() - 1) + "\n\n";
    res_enc = (*cc)->Compress(res_enc, 1);

    auto resp_comp = chrono::high_resolution_clock::now();
    auto return_code = Serial::SerializeToFile("comp_resp.txt", res_enc, SerType::BINARY);
    query_results += "The size of the compressed response is: " + to_string(std::filesystem::file_size(queryFilePath + "comp_resp.txt")) + " bytes." + "\n\n";
    assert(return_code == 1);
    resp_comp_time_count = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - resp_comp)).count();

    auto resp_decomp = chrono::high_resolution_clock::now();
    return_code = Serial::DeserializeFromFile("comp_resp.txt", res_enc, SerType::BINARY);
    assert(return_code == 1);
    resp_decomp_time_count = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - resp_decomp)).count();

    query_results += "Estimated response compression time for COUNT operation: " + to_string(resp_comp_time_count / 1000.0) + " seconds\n";
    query_results += "Estimated response decompression time for COUNT operation: " + to_string(resp_decomp_time_count / 1000.0) + " seconds\n\n";

    double real_public_result = 0;
    long long public_result = 0;
    
    if (pub_res) {
        auto res_pub = pub_pred_res[0];

        for (size_t i = 1; i < pub_pred_res.size(); i++) {
            for (size_t j = 0; j < pub_pred_res[i].size(); j++) {
                res_pub[j] += pub_pred_res[i][j];
            }
        }

        for (auto el: res_pub) {
            public_result += (long long)(round(el));
            real_public_result += el;
        }
    }

    long long expected_result = 0;

    if (exp_res) {
        combineExpectedSubQueries(q);

        for (size_t i = 0; i < exp_pred_res.size(); i++) {
            expected_result += exp_pred_res[i];
        }
    }

    if (real_public_result != 0) {
        stringstream ss;
        ss << std::fixed << std::setprecision(10) << ((real_encrypted_result - real_public_result) / real_public_result) * 100;
        query_results += "Percent error between encrypted count and expected count from ciphertext operations: " 
            + ss.str() + " %\n";
    } else query_results += "Percent error could not be calculated for the encrypted count result since the expected count result is 0\n";
    if (expected_result != 0) {
        stringstream ss;
        ss << std::fixed << std::setprecision(10) << ((real_encrypted_result - expected_result) / expected_result) * 100;
        query_results += "Percent error between encrypted count and actual expected count: " 
            + ss.str() + " %\n";
    } else query_results += "Percent error could not be calculated for the encrypted count result since the actual count result is 0\n";
    if (real_public_result != 0 && expected_result != 0) {
        stringstream ss;
        ss << std::fixed << std::setprecision(10) << ((real_public_result - expected_result) / expected_result) * 100;
        query_results += "Percent error between expected count from ciphertext operations and actual expected count: " 
            + ss.str() + " %\n";
    } else query_results += "Percent error could not be calculated for the expected count from ciphertext operations result since the actual count result is 0\n";
    
    query_results += "\n";

    return {encrypted_result, public_result, expected_result};
}

vector<double> db::sum(query &q, int index) {
    auto sel_time = chrono::high_resolution_clock::now();
    auto EncAggQueryRow = q.getEncAggQueries()[index];
    combineExpectedAggQueries(index, q);
    enc_agg_res = vector<Ciphertext<DCRTPoly>>(this->ct_rows);
    q.combineAggResultsENC(enc_pred_res, exp_agg_res, enc_agg_res);
    selection_ops_time.push_back((chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - sel_time)).count());

    auto res_enc = enc_agg_res[0];
    double agg_enc_result = 0;

    auto agg_time = chrono::high_resolution_clock::now();
    // compress everything down to a single row of ciphertexts
    for (size_t i = 1; i < enc_agg_res.size(); i++) {
        res_enc = (*cc)->EvalAdd(res_enc, enc_agg_res[i]);
    }
    agg_ops_time.push_back((chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - agg_time)).count());

    Plaintext result;
    auto client_agg = chrono::high_resolution_clock::now();
    (*cc)->Decrypt(keys->secretKey, res_enc, &result);
    result->SetLength(numSlots);
    vector<double> realNum;
    for (auto el: result->GetCKKSPackedValue()) {
        agg_enc_result += el.real();
    }
    client_agg_time.push_back((chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - client_agg)).count());

    auto resp_comp = chrono::high_resolution_clock::now();
    query_results += "Remaining levels after SUM operation: " + to_string(multDepth - res_enc->GetLevel() - 1) + "\n\n";
    res_enc = (*cc)->Compress(res_enc, 1);

    auto return_code = Serial::SerializeToFile("comp_resp.txt", res_enc, SerType::BINARY);
    query_results += "The size of the compressed response is: " + to_string(std::filesystem::file_size(queryFilePath + "comp_resp.txt")) + " bytes." + "\n\n";
    assert(return_code == 1);
    resp_comp_time_aggs.push_back((chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - resp_comp)).count());

    auto resp_decomp = chrono::high_resolution_clock::now();
    return_code = Serial::DeserializeFromFile("comp_resp.txt", res_enc, SerType::BINARY);
    assert(return_code == 1);
    resp_decomp_time_aggs.push_back((chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - resp_decomp)).count());

    query_results += "Estimated response compression time for the aggregation operation: " + to_string(resp_comp_time_aggs[resp_comp_time_aggs.size() - 1] / 1000.0) + " seconds\n";
    query_results += "Estimated response decompression time for the aggregation operation: " + to_string(resp_decomp_time_aggs[resp_decomp_time_aggs.size() - 1] / 1000.0) + " seconds\n\n";

    double agg_pub_result = 0;
    if (pub_res) {

        pub_agg_res = vector<vector<double>>(this->ct_rows, vector<double>(numSlots, 0));
        q.combineAggResultsPUB(pub_pred_res, exp_agg_res, pub_agg_res);
        auto res_pub = pub_agg_res[0];
        
        // compress everything down to a single row of ciphertexts
        for (size_t i = 1; i < pub_agg_res.size(); i++) {
            for (size_t j = 0; j < pub_agg_res[0].size(); j++) {
                res_pub[j] += pub_agg_res[i][j];
            }
        }

        for (size_t i = 0; i < res_pub.size(); i++) {
            agg_pub_result += res_pub[i];
        }
    }

    double expected_result = 0;

    if (exp_res) {
        for (size_t i = 0; i < exp_agg_res.size(); i++) {
            if (exp_pred_res[i]) expected_result += exp_agg_res[i];
        }
    }

    return {agg_enc_result, agg_pub_result, expected_result};
}

void db::printReport(query& q) {
    if (!pub_res) cout << "PUBLIC CIPHERTEXT OPERATIONS WERE NOT CALCULATED" << endl;
    if (!exp_res) cout << "EXPECTED CIPHERTEXT OPERATIONS WERE NOT CALCULATED" << endl;
    cout << "\nTime taken for query compression: " << query_comp_time / 1000.0 << " seconds" << endl;
    cout << "Time taken for xz compression: " << query_xz_time / 1000.0 << " seconds" << endl;
    cout << "Time taken for query decompression: " << query_decomp_time / 1000.0 << " seconds" << endl;
    cout << "Time taken for query expansion: " << query_expand_time / 1000.0 << " seconds" << endl << endl;
    auto cur_index = 0;
    for (auto sq : q.getEncSubQueries()) {
        if (sq.op != "(" && sq.op != ")" && sq.op != "AND" && sq.op != "OR") {
            cout << "Time taken for subquery " << sq.h.col_name << " " << sq.op << 
                " " << sq.cmp_val * (1.0 / pow(10, sq.h.precision)) << ": " << 
                pred_ops_time[sq.ID] / 1000.0 << " seconds" << endl;
            total_pred_time += pred_ops_time[sq.ID];
            cur_index++;
        }
    }
    cout << "\nTime taken to evaluate predicate expression: " << pred_comb_time / 1000.0 << " seconds" << endl;
    total_pred_time += pred_comb_time;
    cout << "\nTotal time taken for predicate evaluation: " << total_pred_time / 1000.0 << " seconds" << endl;
    cout << "\nTime taken for the count operation: " << count_time / 1000.0 << " seconds" << endl;

    for (size_t i = 0; i < q.getFinalOperations().size(); i++) {
        cout << "Time taken for client-side aggregation operation " << i + 1 << ": " 
            << client_agg_time[i] / 1000.0 << " seconds" << endl;
        cout << "Time taken for selection operation " << i + 1 << ": " 
            << selection_ops_time[i] / 1000.0 << " seconds" << endl;
        cout << "Time taken for aggregation operation " << i + 1 << ": " 
            << agg_ops_time[i] / 1000.0 << " seconds" << endl;
        cout << "Response compression time for aggregation operation: " << i + 1 << ": "
            << resp_comp_time_aggs[i] / 1000.0 << " seconds" << endl;
        cout << "Response decompression time for aggregation operation: " << i + 1 << ": "
            << resp_decomp_time_aggs[i] / 1000.0 << " seconds" << endl;
        total_agg_time += client_agg_time[i] + agg_ops_time[i] + resp_comp_time_aggs[i] + resp_decomp_time_aggs[i] + selection_ops_time[i];
    }
    cout << "\nTotal time taken to evaluate all aggregation expressions: " << total_agg_time / 1000.0 << " seconds" << endl;
    end_to_end_time = total_pred_time + total_agg_time;
    end_to_end_time += query_comp_time + query_xz_time + query_decomp_time + query_expand_time;
    if (add_count_flag) end_to_end_time += count_time;
    cout << "\nTotal end-to-end query time: " << end_to_end_time / 1000.0 << " seconds" << endl << endl;
    cout << "Query results: " << endl << endl;
    cout << query_results << endl;
}

void db::combineSubQueriesENC(query& q) {
    auto comb_time = chrono::high_resolution_clock::now();
    stack<subQuery> sub;
    for (auto &sq: q.infix_to_postfix(q.getEncSubQueries())) {
        if (std::find(ops.begin(), ops.end(), sq.op) != ops.end()) {
            sub.push(sq);

            if (pred_ops_time.find(sq.ID) == pred_ops_time.end()) {
                auto cur_pred_time = chrono::high_resolution_clock::now();
                q.mapSubQueryResultENC(sq.ID, evaluate_sub_query_ENC(sq, q));
                pred_ops_time.insert({sq.ID, (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - cur_pred_time)).count()});
            }

        } else {

            auto cur_sq = sub.top();
            sub.pop();
            auto prev_sq = sub.top();
            sub.pop();
            q.combineResultsENC(prev_sq.ID, cur_sq.ID, sq.op);
            q.getEncRows().erase(prev_sq.ID);
            sub.push(cur_sq);
        }
    }
    enc_pred_res = q.getEncRows()[sub.top().ID];
    q.getEncRows().erase(sub.top().ID);
    // add a function here that performs predicate comps
    q.PredCompsENC(enc_pred_res);
    pred_comb_time = (chrono::duration_cast<milliseconds>(chrono::high_resolution_clock::now() - comb_time)).count();
    for (auto &t : pred_ops_time) pred_comb_time -= t.second;
}

void db::combineSubQueriesPUB(query &q) {
    stack<subQuery> sub;
    unordered_map<string, bool> isEvaluated;
    for (auto &sq: q.infix_to_postfix(q.getEncSubQueries())) {
        if (std::find(ops.begin(), ops.end(), sq.op) != ops.end()) {

            sub.push(sq);

            if (isEvaluated.find(sq.ID) == isEvaluated.end()) {
                q.mapSubQueryResultPUB(sq.ID, evaluate_sub_query_PUB(sq));
                isEvaluated.insert({sq.ID, true});
            }

        } else {
            auto cur_sq = sub.top();
            sub.pop();
            auto prev_sq = sub.top();
            sub.pop();
            q.combineResultsPUB(prev_sq.ID, cur_sq.ID, sq.op);
            q.getPubRows().erase(prev_sq.ID);
            sub.push(cur_sq);
        }
    }

    pub_pred_res = q.getPubRows()[sub.top().ID];
    q.getPubRows().erase(sub.top().ID);
    q.PredCompsPUB(pub_pred_res);
}

void evaluate_expected_count_helper(vector<long long>& result, vector<double>& expected_data, subQuery& sq, size_t i) {
    if (i >= start_end_expected_threads.size()) return;
    size_t start = start_end_expected_threads[i - 1];
    size_t end = start_end_expected_threads[i];

    for (size_t i = start; i < end; i++) {
        long long cor_val = expected_data[i] * (long long)(pow(10, sq.h.precision));
        if (sq.op == ">=" && cor_val >= sq.cmp_val) result[i] = 1;
        else if (sq.op == "<=" && cor_val <= sq.cmp_val) result[i] = 1;
        else if (sq.op == ">" && cor_val > sq.cmp_val) result[i] = 1;
        else if (sq.op == "<" && cor_val < sq.cmp_val) result[i] = 1;
        else if (sq.op == "==" && cor_val == sq.cmp_val) result[i] = 1;
        else if (sq.op == "!=" && cor_val != sq.cmp_val) result[i] = 1;
        else result[i] = 0;
    }
}

void db::combineExpectedSubQueries(query &q) {
    stack<subQuery> sub;
    unordered_map<string, bool> isEvaluated;
    for (auto &sq: q.infix_to_postfix(q.getEncSubQueries())) {
        if (std::find(ops.begin(), ops.end(), sq.op) != ops.end()) {
            sub.push(sq);

            if (isEvaluated.find(sq.ID) == isEvaluated.end()) {
                vector<long long> result(num_rows);
                for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
                    MACHINE_THREADS[i - 1] = thread(evaluate_expected_count_helper, ref(result), ref(expected_data[sq.h.col_name]), ref(sq), i);
                }

                for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
                    MACHINE_THREADS[i - 1].join();
                }
                q.mapSubQueryExpectedResult(sq.ID, result); 
                isEvaluated.insert({sq.ID, true});
            }

        } else {
            auto cur_sq = sub.top();
            sub.pop();
            auto prev_sq = sub.top();
            sub.pop();
            q.combineExpResults(prev_sq.ID, cur_sq.ID, sq.op);
            q.getExpectedRows().erase(prev_sq.ID);
            sub.push(cur_sq);
        }
    }
    exp_pred_res = q.getExpectedRows()[sub.top().ID];
    q.getExpectedRows().erase(sub.top().ID);
}

void db::combineExpectedAggQueries(int index, query& q) {
    stack<subQuery> sub;
    unordered_map<string, bool> isEvaluated;
    for (auto &sq: q.infix_to_postfix(q.getEncAggQueries()[index])) {
        if (sq.op == "") {
            sub.push(sq);

            if (isEvaluated.find(sq.ID) == isEvaluated.end()) {
                q.mapAggQueryExpectedResult(sq.ID, expected_data[sq.h.col_name]);
                isEvaluated.insert({sq.ID, true});
            }

        } else {
            auto cur_sq = sub.top();
            sub.pop();
            auto prev_sq = sub.top();
            sub.pop();
            q.combineAggExpResults(prev_sq.ID, cur_sq.ID, sq.op);
            q.getAggExpectedRows().erase(prev_sq.ID);
            sub.push(cur_sq);
        }
    }
    exp_agg_res = q.getAggExpectedRows()[sub.top().ID];
    q.getAggExpectedRows().erase(sub.top().ID);
}

void db::evaluate_query(query& q) {
    start_end_threads = vector<size_t>(min((uint32_t)(TOTAL_MACHINE_THREAD), ct_rows));
    for (size_t i = 0; i < this->ct_rows; i++) {
        start_end_threads[i % TOTAL_MACHINE_THREAD] += 1;
    }
    start_end_threads.insert(start_end_threads.begin(), 0);
    for (size_t i = 1; i < start_end_threads.size(); i++) {
        start_end_threads[i] += start_end_threads[i - 1];
    }

    start_end_expected_threads = vector<size_t>(min((uint32_t)(TOTAL_MACHINE_THREAD), num_rows));
    for (size_t i = 0; i < this->num_rows; i++) {
        start_end_expected_threads[i % TOTAL_MACHINE_THREAD] += 1;
    }
    start_end_expected_threads.insert(start_end_expected_threads.begin(), 0);
    for (size_t i = 1; i < start_end_expected_threads.size(); i++) {
        start_end_expected_threads[i] += start_end_expected_threads[i - 1];
    }

    combineSubQueriesENC(q);
    if (pub_res) combineSubQueriesPUB(q);

    auto res_count = count(q);

    long long enc_count = res_count[0];
    long long pub_count = res_count[1];
    long long exp_count = res_count[2];

    query_results += "ENC_COUNT: " + to_string(enc_count) + "\n";
    query_results += "PUB_COUNT: " + to_string(pub_count) + "\n";
    query_results += "EXP_COUNT: " + to_string(exp_count) + "\n\n";

    for (size_t i = 0; i < q.getFinalOperations().size(); i++) {
        query_results += "Results for aggregation query at index " + to_string(i + 1) + ":\n";
        if (q.getFinalOperations()[i] == "COUNT") {
            query_results += "Estimated encrypted count: " + to_string(enc_count) + "\n";
            query_results += "Expected count from ciphertext operations: " + to_string(enc_count) + "\n";
            query_results += "Actual expected count: " + to_string(enc_count) + "\n";
            add_count_flag = 1;
        }

        else if (q.getFinalOperations()[i]  == "SUM") {
            auto agg_res = sum(q, i);

            query_results += "Estimated encrypted sum: " + to_string(agg_res[0]) + "\n";
            query_results += "Expected sum from ciphertext operations: " + to_string(agg_res[1]) + "\n";
            query_results += "Actual expected sum: " + to_string(agg_res[2]) + "\n";
            if (agg_res[1] != 0) {
                stringstream ss;
                ss << std::fixed << std::setprecision(10) << ((agg_res[0] - agg_res[1]) / agg_res[1]) * 100;
                query_results += "Percent error between encrypted sum and expected sum from ciphertext operations: " 
                    + ss.str() + " %\n";
            } else query_results += "Percent error could not be calculated for the encrypted sum result since the expected sum result is 0\n";

            if (agg_res[2] != 0) {
                stringstream ss;
                ss << std::fixed << std::setprecision(10) << ((agg_res[0] - agg_res[2]) / agg_res[2]) * 100;
                query_results += "Percent error between encrypted sum and actual expected sum: " 
                    + ss.str() + " %\n";
            } else query_results += "Percent error could not be calculated for the encrypted sum result since the actual sum result is 0\n";

            if (agg_res[1] != 0 && agg_res[2] != 0) {
                stringstream ss;
                ss << std::fixed << std::setprecision(10) << ((agg_res[1] - agg_res[2]) / agg_res[2]) * 100;
                query_results += "Percent error between expected sum from ciphertext operations and actual expected sum: " 
                    + ss.str() + " %\n";
            } else query_results += "Percent error could not be calculated for the expected sum from ciphertext operations result since the actual sum result is 0\n";

        }

        else if (q.getFinalOperations()[i] == "AVG") {
            auto agg_res = sum(q, i);
            add_count_flag = 1;

            if (enc_count == 0) {
                query_results += "No elements were selected to compute the average from the encrypted operations\n";
            } else {
                auto enc_avg = agg_res[0] / (double)(enc_count);
                auto est_avg = agg_res[1] / (double)(pub_count);            
                query_results += "Estimated encrypted average" + to_string(enc_avg) + "\n";
                query_results += "Expected average from ciphertext operations" + to_string(enc_avg) + "\n";
                if (est_avg != 0) {
                    stringstream ss;
                    ss << std::fixed << std::setprecision(10) << ((enc_avg - est_avg) / est_avg) * 100;
                    query_results += "Percent error between encrypted average and expected average from ciphertext operations: "
                    + ss.str() + " %\n";
                } else query_results += "Percent error could not be calculated for the encrypted average because the expected average from ciphertext operations is 0\n";

            }

            if (exp_count != 0) {
                auto enc_avg = agg_res[0] / (double)(enc_count);
                auto act_avg = agg_res[2] / (double)(exp_count);
                query_results += "Actual expected average " + to_string(act_avg) + "\n";
                if (act_avg != 0) {
                    stringstream ss;
                    ss << std::fixed << std::setprecision(10) << ((enc_avg - act_avg) / act_avg) * 100;
                    query_results += "Percent error between encrypted average and actual expected average: "
                    + ss.str() + " %\n";
                } else query_results += "Percent error could not be calculated for the encrypted average because the actual expected average is 0\n";

            } else {
                query_results += "No elements were selected to compute the average from the expected operations\n";
            }

            if (exp_count != 0) {
                auto est_avg = agg_res[1] / (double)(enc_count);
                auto act_avg = agg_res[2] / (double)(exp_count);
                if (act_avg != 0) {
                    stringstream ss;
                    ss << std::fixed << std::setprecision(10) << ((est_avg - act_avg) / act_avg) * 100;
                    query_results += "Percent error between expected average from ciphertext operations and actual expected average: "
                    + ss.str() + " %\n";
                } else query_results += "Percent error could not be calculated for the expected average from ciphertext operations because the actual expected average is 0\n";

            } else {
                query_results += "No elements were selected to compute the average from the expected operations\n";
            }

            cout << endl;
        } 

        // MIN and MAX return the minimum and maximum element found in a header,
        // Range query support for these operations are not supported in this work.
        // These functionalities are purely to identify the range of values that are
        // actually present in a header as opposed to its own lower and upper bounds
        else {
            combineExpectedSubQueries(q);
            combineExpectedAggQueries(i, q);

            auto exp_count = 0;
            for (size_t j = 0; j < exp_pred_res.size(); j++) exp_count += exp_pred_res[j];

            if (exp_count == 0) {
                query_results += "No elements have been selected for evaluating MIN or MAX\n";
                return;
            }

            else if (q.getFinalOperations()[i] == "MIN") {
                size_t j = 0;
                while (exp_pred_res[j] != 1) j += 1;

                double Min = exp_agg_res[j];
                for (size_t i = j + 1; i < num_rows; i++) {
                    if (exp_pred_res[i] == 1) Min = min(Min, exp_agg_res[i]);
                }

                query_results += "Min element: " + to_string(Min) + "\n";
            }

            else {
                size_t j = 0;
                while (exp_pred_res[j] != 1) j += 1;

                double Max = exp_agg_res[j];
                for (size_t i = j + 1; i < num_rows; i++) {
                    if (exp_pred_res[i] == 1) Max = max(Max, exp_agg_res[i]);
                }

                query_results += "Max element: " + to_string(Max) + "\n";
            }
        }
        query_results += "\n";
    }
    printReport(q);
}