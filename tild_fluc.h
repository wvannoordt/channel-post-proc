#pragma once
#include <vector>
static void comp_tild(std::vector<double>& phi_tilde, const std::vector<double>& rho_phi, const std::vector<double>& rho)
{
    for (int i = 0; i < phi_tilde.size(); i++)
    {
        phi_tilde[i] = rho_phi[i]/rho[i];
    }
}

static void comp_fluc(
    std::vector<double>& phi_fluc,
    const std::vector<double>& phi,
    const std::vector<double>& rho,
    const std::vector<double>& rho_phi,
    const std::vector<double>& rho2_phi2)
{
    for (int i = 0; i < phi_fluc.size(); i++)
    {
        phi_fluc[i] = phi[i]*phi[i] - 2.0*phi[i]*rho_phi[i]/rho[i] + rho2_phi2[i]/(rho[i]*rho[i]);
    }
}

static void comp_fluc_i(
    std::vector<double>& phi_fluc,
    const std::vector<double>& ab,
    const std::vector<double>& a,
    const std::vector<double>& b)
{
    for (int i = 0; i < phi_fluc.size(); i++)
    {
        phi_fluc[i] = ab[i] - a[i]*b[i];
    }
}
