#pragma once

template <typename data_t> struct m9
{
    data_t v[9];
    data_t& operator()(const std::size_t& row, const std::size_t& col)
    {
        return v[col + 3*row];
    }
    const data_t& operator()(const std::size_t& row, const std::size_t& col) const
    {
        return v[col + 3*row];
    }
};