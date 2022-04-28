#pragma once
#include "cmf.h"
#include <vector>
template <typename data_t> struct avg_qty_t
{
    int numEntries;
    std::vector<data_t> accumulated_data;
    std::vector<data_t> inst_data;
    std::vector<avg_qty_t<double>*>* registry;
    avg_qty_t(int dsize, std::vector<avg_qty_t<double>*>& reg_in)
    {
        accumulated_data.resize(dsize, 0.0);
        inst_data.resize(dsize);
        registry = &reg_in;
        numEntries = 0;
        reg_in.push_back(this);
    }
    std::size_t size(void) const { return accumulated_data.size(); }
    const data_t& operator [] (const std::size_t i) const {return accumulated_data[i];}
    data_t& operator [] (const std::size_t i) {return accumulated_data[i];}
    void accum(void)
    {
        double beta  = 1.0/(numEntries+1);
        double alpha = 1.0 - beta;
        for (std::size_t i = 0; i < this->size(); i++)
        {
            accumulated_data[i] = alpha*accumulated_data[i] + beta*inst_data[i];
        }
        numEntries++;
    }
};