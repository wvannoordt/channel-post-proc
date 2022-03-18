#pragma once

#include <string>
#include <iostream>
#include <fstream>


template <typename ltype, typename rtype> auto min(const ltype& a, const rtype& b)
{
    return a<b?a:b;
}

namespace detail
{
    template <class indexable_t> std::size_t get_min_size_r(const indexable_t& i)
    {
        return i.size();
    }
    template <class indexable_t, class... indexables_t> std::size_t get_min_size_r(const indexable_t& i, indexables_t... is)
    {
        return min(i.size(), get_min_size_r(is...));
    }

    template <class indexable_t> void write_r(std::ofstream& myfile, const std::size_t idx, const indexable_t& i)
    {
        myfile << i[idx] << "\n";
    }
    template <class indexable_t, class... indexables_t> void write_r(std::ofstream& myfile, const std::size_t idx, const indexable_t& i, indexables_t... is)
    {
        myfile << i[idx] << ", ";
        write_r(myfile, idx, is...);
    }
}
template <class... indexable_t> static void save_csv(
    const std::string& filename,
    // const std::vector<std::string> names,
    indexable_t... vecs)
{
    std::ofstream myfile(filename);
    // for (int j = 0; j < names.size(); j++)
    // {
    //     myfile << names[j];
    //     if (j<names.size()-1) myfile << ", ";
    // }
    // myfile << "\n";
    std::size_t minsize = detail::get_min_size_r(vecs...);
    for (std::size_t i = 0; i < minsize; i++)
    {
        detail::write_r(myfile, i, vecs...);
    }
}