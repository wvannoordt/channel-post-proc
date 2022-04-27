#pragma once
#include "cmf.h"
#include "compute_average.h"
using cmf::print;
template <class callable, typename dtype=double>
void extract_profile(
    cmf::CartesianMeshArray& ar,
    const std::size_t& j,
    const std::size_t& k,
    std::vector<dtype>& output,
    const callable& func)
{
    auto pp  = ar.Mesh()->GetBlockDim();
    auto nxb = ar.Mesh()->GetMeshDataDim();
    cmf::Vec3<int> numBlocks(pp[0], pp[1], pp[2]);
    auto vb = ar.Mesh()->Blocks()->GetBlockBounds();
    cmf::Vec3<double> dxDom((vb[1]-vb[0])/numBlocks[0], (vb[3]-vb[2])/numBlocks[1], (vb[5]-vb[4])/numBlocks[2]);
    for (auto lb: ar)
    {
        cmf::BlockArray<dtype,1> data = ar[lb];
        cmf::BlockInfo info = ar.GetBlockInfo(lb);
        auto c = lb->GetBlockCenter();
        cmf::Vec3<int> idx_block;
        for (auto i: range(0, cmf::Dim())) idx_block[i[0]] = (int)((c[i[0]]-vb[2*i[0]]) / (dxDom[i[0]]));
        bool block_contains_profile = (j/nxb[1]==idx_block[1]) && (k/nxb[2]==idx_block[2]);
        if (block_contains_profile)
        {
            int jloc = j-idx_block[1]*nxb[1];
            int kloc = k-idx_block[2]*nxb[2];
            for (auto iloc: range(0, nxb[0]))
            {
                cmf::Vec3<int> ijk(iloc[0], jloc, kloc);
                cmf::Vec3<dtype> xyz(
                    info.blockBounds[0]+((dtype)ijk[0]+0.5)*info.dx[0],
                    info.blockBounds[2]+((dtype)ijk[1]+0.5)*info.dx[1],
                    info.blockBounds[4]+((dtype)ijk[2]+0.5)*info.dx[2]);
                prim_t<dtype> dataPrimLoc;
                for (auto n: range(0,5)) dataPrimLoc[n[0]] = data(n[0], ijk[0], ijk[1], ijk[2]);
                val_grad<3, prim_t<dtype>> qgrad;
                for (auto n: range(0,5))
                {
                    qgrad[0] = (data(n[0], ijk[0]+1, ijk[1],   ijk[2])   - data(n[0], ijk[0]-1, ijk[1],   ijk[2]))  /(2.0*info.dx[0]);
                    qgrad[1] = (data(n[0], ijk[0],   ijk[1]+1, ijk[2])   - data(n[0], ijk[0],   ijk[1]-1, ijk[2]))  /(2.0*info.dx[1]);
                    qgrad[2] = (data(n[0], ijk[0],   ijk[1],   ijk[2]+1) - data(n[0], ijk[0],   ijk[1],   ijk[2])-1)/(2.0*info.dx[2]);
                }
                int idx = iloc[0] + idx_block[0]*nxb[0];
                dtype dataNew = detail::forward_prim_callable_reduce_args(func, xyz, ijk, dataPrimLoc, qgrad, idx);
                output[idx] = dataNew;
            }
        }
    }
}