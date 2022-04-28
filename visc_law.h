#pragma once

#include "fluid_state.h"
#include <cmath>

template <typename dtype> struct visc_power_law
{
    visc_power_law(const dtype& mu_ref_in, const dtype& t_ref_in, const dtype& power_in)
    {
        mu_ref = mu_ref_in;
        power  = power_in;
        t_ref  = t_ref_in;
    }
    dtype mu_ref, t_ref, power;
    dtype calc_visc(const prim_t<dtype>& q) const
    {
        return mu_ref*pow(q.T()/t_ref, power);
    }
    dtype calc_beta(const prim_t<dtype>& q) const
    {
        return -0.66666666666667*calc_visc(q);
    }
};