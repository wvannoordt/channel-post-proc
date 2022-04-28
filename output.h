#pragma once

#include <string>
#include <iostream>
#include <fstream>


template <typename ltype, typename rtype> auto min(const ltype& a, const rtype& b)
{
    return a<b?a:b;
}
template <typename ltype, typename rtype> auto max(const ltype& a, const rtype& b)
{
    return a>b?a:b;
}

std::string str_pad(const std::string& str, const char& pad_char, const std::size_t& len)
{
    return str;
    std::string output = str;
    while (output.length() < len) output += pad_char;
    return output;
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

    template <class indexable_t> void write_r(const std::size_t& pad_size, std::ofstream& myfile, const std::size_t idx, const indexable_t& i)
    {
        myfile << str_pad(std::to_string(i[idx]), ' ', pad_size) << "\n";
    }
    template <class indexable_t, class... indexables_t> void write_r(const std::size_t& pad_size, std::ofstream& myfile, const std::size_t idx, const indexable_t& i, indexables_t... is)
    {
        myfile << str_pad(std::to_string(i[idx])+',', ' ', pad_size);
        write_r(pad_size, myfile, idx, is...);
    }
}
template <class... indexable_t> static void save_csv(
    const std::string& filename,
    const std::vector<std::string> names,
    indexable_t... vecs)
{
    std::ofstream myfile(filename);
    std::size_t max_size = 0;
    for (int j = 0; j < names.size(); j++)
    {
        max_size = max(max_size, names[j].length());
    }
    std::size_t field_size = 2+max_size;
    for (int j = 0; j < names.size(); j++)
    {
        myfile << str_pad(names[j] + ((j<names.size()-1)?',':' '), ' ', field_size);
    }
    myfile << "\n";
    std::size_t minsize = detail::get_min_size_r(vecs...);
    for (std::size_t i = 0; i < minsize; i++)
    {
        detail::write_r(field_size, myfile, i, vecs...);
    }
}