#pragma once
#include <iostream>
template <const std::size_t ar_dim, class data_t> struct val_grad
{
    data_t data[ar_dim];
    data_t& operator [] (const std::size_t& idx) {return data[idx];}
    const data_t& operator [] (const std::size_t& idx) const {return data[idx];}
};