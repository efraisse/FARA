#include "globals.h"

#pragma once

class header {
    public:
        header(string col_name, long long low = 0, long long high = 0, uint32_t precision = 0) {
            this->col_name = col_name;
            this->num_ct = 1;
            this->low = low;
            this->high = high;
            this->precision = precision;
        }

        header() {}

        header(const header &h) {
            this->col_name = h.col_name;
            this->num_ct = h.num_ct;
            this->low = h.low;
            this->high = h.high;
            this->precision = h.precision;
        }
        
        string col_name;
        uint32_t num_ct;
        long long low;
        long long high;
        uint32_t precision;
};