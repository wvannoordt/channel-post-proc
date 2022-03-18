#pragma once

#include <concepts>
#include <array>

#include "cmf.h"
#include "fluid_state.h"
#include "range.h"

template <class T, typename dtype> concept void_callable = requires(T t)
{
    { t() } -> std::same_as<prim_t<dtype>>;
};

template <class T, typename dtype> concept prim_callable = requires(T t, const prim_t<dtype>& p)
{
    { t(p) } -> std::same_as<prim_t<dtype>>;
};

namespace detail
{
    template <typename dtype, void_callable<dtype> func_t> prim_t<dtype> forward_prim_callable_args(
        const func_t& func, cmf::Vec3<double>& xyz, cmf::Vec3<int>& ijk, const prim_t<dtype>& prev_data)
    {
        return func();
    }
    template <typename dtype, prim_callable<dtype> func_t> prim_t<dtype> forward_prim_callable_args(
        const func_t& func, cmf::Vec3<double>& xyz, cmf::Vec3<int>& ijk, const prim_t<dtype>& prev_data)
    {
        return func(prev_data);
    }
}

template <class callable, typename dtype=double> void transform_inplace(cmf::CartesianMeshArray& ar, const callable& func)
{
    for (auto& lb: ar)
    {
        cmf::BlockArray<dtype,1> data = ar[lb];
        cmf::cell_t i0 = data.imin-data.exchangeI;
        cmf::cell_t i1 = data.imax+data.exchangeI;
        cmf::cell_t j0 = data.jmin-data.exchangeJ;
        cmf::cell_t j1 = data.jmax+data.exchangeJ;
        cmf::cell_t k0 = data.kmin-data.exchangeK;
        cmf::cell_t k1 = data.kmax+data.exchangeK;
        cmf::BlockInfo info = ar.GetBlockInfo(lb);
        for (auto i: range(i0,i1,j0,j1,k0,k1))
        {
            cmf::Vec3<dtype> xyz(info.blockBounds[0]+i[0]*info.dx[0],info.blockBounds[2]+i[1]*info.dx[1],info.blockBounds[4]+i[2]*info.dx[2]);
            cmf::Vec3<int> ijk(i[0],i[1],i[2]);
            prim_t<dtype> dataPrimLoc;
            for (auto n: range(0,5)) dataPrimLoc[n[0]] = data(n[0], i[0], i[1], i[2]);
            prim_t<dtype> dataNew = detail::forward_prim_callable_args(func, xyz, ijk, dataPrimLoc);
            for (auto n: range(0,5)) data(n[0], i[0], i[1], i[2]) = dataNew[n[0]];
        }
    }
}