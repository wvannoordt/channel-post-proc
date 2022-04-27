#pragma once

#include <concepts>

#include "cmf.h"

#include "range.h"
#include "fluid_state.h"
#include "val_grad.h"

template <class T, typename dtype> concept void_reduce_callable = requires(T t)
{
    { t()  } -> std::same_as<dtype>;
};

template <class T, typename dtype> concept prim_reduce_callable = requires(T t, const prim_t<dtype>& p)
{
    { t(p) } -> std::same_as<dtype>;
};

template <class T, typename dtype> concept prim_qgrad_reduce_callable = requires(T t, const prim_t<dtype>& p, const val_grad<3, prim_t<dtype>>& grad)
{
    { t(p, grad) } -> std::same_as<dtype>;
};

template <class T, typename dtype> concept xyz_reduce_callable = requires(T t, const cmf::Vec3<double>& p)
{
    { t(p) } -> std::same_as<dtype>;
};

template <class T, typename dtype> concept primi_reduce_callable = requires(T t, const prim_t<dtype>& p, const std::size_t& j)
{
    { t(p,j) } -> std::same_as<dtype>;
};

namespace detail
{
    template <typename dtype, void_reduce_callable<dtype> func_t> dtype forward_prim_callable_reduce_args(
        const func_t& func,
        const cmf::Vec3<double>& xyz,
        cmf::Vec3<int>& ijk,
        const prim_t<dtype>& prev_data,
        const val_grad<3, prim_t<dtype>>& grad,
        const std::size_t& jblock)
    {
        return func();
    }
    template <typename dtype, prim_reduce_callable<dtype> func_t> dtype forward_prim_callable_reduce_args(
        const func_t& func,
        const cmf::Vec3<double>& xyz,
        cmf::Vec3<int>& ijk,
        const prim_t<dtype>& prev_data,
        const val_grad<3, prim_t<dtype>>& grad,
        const std::size_t& jblock)
    {
        return func(prev_data);
    }
    
    template <typename dtype, prim_qgrad_reduce_callable<dtype> func_t> dtype forward_prim_callable_reduce_args(
        const func_t& func,
        const cmf::Vec3<double>& xyz,
        cmf::Vec3<int>& ijk,
        const prim_t<dtype>& prev_data,
        const val_grad<3, prim_t<dtype>>& grad,
        const std::size_t& jblock)
    {
        return func(prev_data, grad);
    }
    
    template <typename dtype, xyz_reduce_callable<dtype> func_t> dtype forward_prim_callable_reduce_args(
        const func_t& func,
        const cmf::Vec3<double>& xyz,
        cmf::Vec3<int>& ijk,
        const prim_t<dtype>& prev_data,
        const val_grad<3, prim_t<dtype>>& grad,
        const std::size_t& jblock)
    {
        return func(xyz);
    }
    template <typename dtype, primi_reduce_callable<dtype> func_t> dtype forward_prim_callable_reduce_args(
        const func_t& func,
        const cmf::Vec3<double>& xyz,
        cmf::Vec3<int>& ijk,
        const prim_t<dtype>& prev_data,
        const val_grad<3, prim_t<dtype>>& grad,
        const std::size_t& jblock)
    {
        return func(prev_data, jblock);
    }
}

template <class callable, typename dtype=double> void compute_average(cmf::CartesianMeshArray& ar, std::vector<dtype>& output, const callable& func)
{
    for (auto& d:output) d = 0.0;
    int numBlocksY = ar.Mesh()->GetBlockDim()[1];
    auto vb = ar.Mesh()->Blocks()->GetBlockBounds();
    double dxDom = (vb[3]-vb[2])/numBlocksY;
    std::vector<std::size_t> counts;
    counts.resize(output.size(), 0);
    for (auto& lb: ar)
    {
        cmf::BlockArray<dtype,1> data = ar[lb];
        cmf::cell_t i0 = data.imin;
        cmf::cell_t i1 = data.imax;
        cmf::cell_t j0 = data.jmin;
        cmf::cell_t j1 = data.jmax;
        cmf::cell_t k0 = data.kmin;
        cmf::cell_t k1 = data.kmax;
        cmf::BlockInfo info = ar.GetBlockInfo(lb);
        for (auto i: range(i0,i1,j0,j1,k0,k1))
        {
            cmf::Vec3<dtype> xyz(
                info.blockBounds[0]+((dtype)i[0]+0.5)*info.dx[0],
                info.blockBounds[2]+((dtype)i[1]+0.5)*info.dx[1],
                info.blockBounds[4]+((dtype)i[2]+0.5)*info.dx[2]);
            cmf::Vec3<int> ijk(i[0],i[1],i[2]);
            prim_t<dtype> q;
            for (auto n: range(0,5)) q[n[0]] = data(n[0], i[0], i[1], i[2]);
            val_grad<3, prim_t<dtype>> qgrad;
            for (auto n: range(0,5))
            {
                qgrad[0] = (data(n[0], i[0]+1, i[1],   i[2])   - data(n[0], i[0]-1, i[1],   i[2]))  /(2.0*info.dx[0]);
                qgrad[1] = (data(n[0], i[0],   i[1]+1, i[2])   - data(n[0], i[0],   i[1]-1, i[2]))  /(2.0*info.dx[1]);
                qgrad[2] = (data(n[0], i[0],   i[1],   i[2]+1) - data(n[0], i[0],   i[1],   i[2])-1)/(2.0*info.dx[2]);
            }
            auto c = lb->GetBlockCenter();
            std::size_t jblock = (int)((c[1]-vb[2]) / (dxDom));
            int idx = i[1] + (j1-j0)*jblock;
            dtype dataNew = detail::forward_prim_callable_reduce_args(func, xyz, ijk, q, qgrad, idx);
            output[idx] += dataNew;
            counts[idx] ++;
        }
    }
    for (std::size_t i = 0; i < counts.size(); i++)
    {
        output[i] /= counts[i];
    }
}