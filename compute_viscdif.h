#pragma once

#include "cmf.h"
#include "m9.h"

template <typename vlaw_t, typename dtype=double> void compute_viscdif(cmf::CartesianMeshArray& ar,
    const vlaw_t& visc,
    std::vector<dtype>& tau_00_0,
    std::vector<dtype>& tau_01_1,
    std::vector<dtype>& tau_02_2,
    std::vector<dtype>& tau_10_0,
    std::vector<dtype>& tau_11_1,
    std::vector<dtype>& tau_12_2,
    std::vector<dtype>& tau_20_0,
    std::vector<dtype>& tau_21_1,
    std::vector<dtype>& tau_22_2)
{
    {
        for (auto& d:tau_00_0) d = 0.0;
        for (auto& d:tau_01_1) d = 0.0;
        for (auto& d:tau_02_2) d = 0.0;
        for (auto& d:tau_10_0) d = 0.0;
        for (auto& d:tau_11_1) d = 0.0;
        for (auto& d:tau_12_2) d = 0.0;
        for (auto& d:tau_20_0) d = 0.0;
        for (auto& d:tau_21_1) d = 0.0;
        for (auto& d:tau_22_2) d = 0.0;
        int numBlocksY = ar.Mesh()->GetBlockDim()[1];
        auto vb = ar.Mesh()->Blocks()->GetBlockBounds();
        double dxDom = (vb[3]-vb[2])/numBlocksY;
        std::vector<std::size_t> counts;
        counts.resize(tau_20_0.size(), 0);
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
                
                auto calc_tau = [&](int ii, int jj, int kk) -> m9<double>
                {
                    val_grad<3, prim_t<double>> qgrad;
                    prim_t<double> q;
                    for (auto n: range(0,5))
                    {
                        q[n[0]]        = data(n[0],ii,jj,kk);
                        qgrad[0][n[0]] = (data(n[0], ii+1, jj,   kk)   - data(n[0], ii-1, jj,   kk))  /(2.0*info.dx[0]);
                        qgrad[1][n[0]] = (data(n[0], ii,   jj+1, kk)   - data(n[0], ii,   jj-1, kk))  /(2.0*info.dx[1]);
                        qgrad[2][n[0]] = (data(n[0], ii,   jj,   kk+1) - data(n[0], ii,   jj,   kk)-1)/(2.0*info.dx[2]);
                    }
                    double mu   = visc.calc_visc(q);
                    double beta = visc.calc_visc(q);
                    m9<double> output;
                    output(0,0) = visc.calc_visc(q)*(qgrad[0].u()+qgrad[0].u()) + visc.calc_beta(q)*(qgrad[0].u() + qgrad[1].v() + qgrad[2].w());
                    output(0,1) = visc.calc_visc(q)*(qgrad[1].u()+qgrad[0].v());
                    output(0,2) = visc.calc_visc(q)*(qgrad[2].u()+qgrad[0].w());
                    output(1,0) = visc.calc_visc(q)*(qgrad[0].v()+qgrad[1].u());
                    output(1,1) = visc.calc_visc(q)*(qgrad[1].v()+qgrad[1].v()) + visc.calc_beta(q)*(qgrad[0].u() + qgrad[1].v() + qgrad[2].w());
                    output(1,2) = visc.calc_visc(q)*(qgrad[2].v()+qgrad[1].w());
                    output(2,1) = visc.calc_visc(q)*(qgrad[0].w()+qgrad[2].u());
                    output(2,0) = visc.calc_visc(q)*(qgrad[1].w()+qgrad[2].v());
                    output(2,2) = visc.calc_visc(q)*(qgrad[2].w()+qgrad[2].w()) + visc.calc_beta(q)*(qgrad[0].u() + qgrad[1].v() + qgrad[2].w());
                    return output;
                };
                
                auto c = lb->GetBlockCenter();
                std::size_t jblock = (int)((c[1]-vb[2]) / (dxDom));
                int idx = i[1] + (j1-j0)*jblock;
                int ii0 = i[0];
                int jj0 = i[1];
                int kk0 = i[2];
                auto tau_0p = calc_tau(ii0+1,jj0,kk0);
                auto tau_0m = calc_tau(ii0-1,jj0,kk0);
                auto tau_1p = calc_tau(ii0,jj0+1,kk0);
                auto tau_1m = calc_tau(ii0,jj0-1,kk0);
                auto tau_2p = calc_tau(ii0,jj0,kk0+1);
                auto tau_2m = calc_tau(ii0,jj0,kk0-1);
                double tau_00_0_loc = (tau_0p(0,0)-tau_0m(0,0))/(2.0*info.dx[0]);
                double tau_01_1_loc = (tau_1p(0,1)-tau_1m(0,1))/(2.0*info.dx[1]);
                double tau_02_2_loc = (tau_2p(0,2)-tau_2m(0,2))/(2.0*info.dx[2]);
                double tau_10_0_loc = (tau_0p(1,0)-tau_0m(1,0))/(2.0*info.dx[0]);
                double tau_11_1_loc = (tau_1p(1,1)-tau_1m(1,1))/(2.0*info.dx[1]);
                double tau_12_2_loc = (tau_2p(1,2)-tau_2m(1,2))/(2.0*info.dx[2]);
                double tau_20_0_loc = (tau_0p(2,0)-tau_0m(2,0))/(2.0*info.dx[0]);
                double tau_21_1_loc = (tau_1p(2,1)-tau_1m(2,1))/(2.0*info.dx[1]);
                double tau_22_2_loc = (tau_2p(2,2)-tau_2m(2,2))/(2.0*info.dx[2]);
                
                tau_00_0[idx] += tau_00_0_loc;
                tau_01_1[idx] += tau_01_1_loc;
                tau_02_2[idx] += tau_02_2_loc;
                tau_10_0[idx] += tau_10_0_loc;
                tau_11_1[idx] += tau_11_1_loc;
                tau_12_2[idx] += tau_12_2_loc;
                tau_20_0[idx] += tau_20_0_loc;
                tau_21_1[idx] += tau_21_1_loc;
                tau_22_2[idx] += tau_22_2_loc;
                counts[idx] ++;
            }
        }
        for (std::size_t i = 0; i < counts.size(); i++)
        {
            tau_00_0[i] /= counts[i];
            tau_01_1[i] /= counts[i];
            tau_02_2[i] /= counts[i];
            tau_10_0[i] /= counts[i];
            tau_11_1[i] /= counts[i];
            tau_12_2[i] /= counts[i];
            tau_20_0[i] /= counts[i];
            tau_21_1[i] /= counts[i];
            tau_22_2[i] /= counts[i];
        }
    }
}