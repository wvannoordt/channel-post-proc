#include "cmf.h"
#include "transform_inplace.h"
#include "range.h"
#include "compute_average.h"
#include "compute_viscdif.h"
#include "output.h"
#include "tild_fluc.h"
#include "timesteps.h"
#include "visc_law.h"
#include "avg_qty_t.h"

using cmf::print;
using cmf::strformat;
using cmf::ZFill;

int main(int argc, char** argv)
{
	int start = timesteps::startstep;
	int end   = timesteps::endstep;
	int skip  = timesteps::skipstep;
	
	std::string inFile = "input.ptl";
    cmf::ReadInput(inFile);
    cmf::globalSettings = cmf::GlobalSettings(cmf::mainInput["GlobalSettings"]);
	cmf::CreateParallelContext(&argc, &argv);
	
	std::string gid_file = cmf::strformat("data/gridInterpolationInfo_{}.bin", cmf::ZFill(start,8));
	std::string rba_file = cmf::strformat("data/restart_block_arrangement_nt_{}.dat", cmf::ZFill(start,8));
	
	cmf::LegacyRestartReader reader(gid_file, rba_file);
	cmf::CartesianMeshInputInfo inputInfo = reader.ReadMeshInfo();
	cmf::CartesianMesh domain(inputInfo);
	reader.ConformMesh(domain);
	// auto& primsAvg  = domain.DefineVariable("primsAvg",  cmf::CmfArrayType::CmfDouble, {5});
	auto& primsInst = domain.DefineVariable("primsInst", cmf::CmfArrayType::CmfDouble, {5});
	
	// transform_inplace(primsAvg, []() -> prim_t<double> {return prim_t<double>(0.0);});
	
	
	perfect_gas_t<double> air;
	air.R = 287.15;
	air.gamma = 1.4;
	
	visc_power_law<double> visc(3.0e-4, 100.0, 0.76);
	
	auto cons2prim = [=](const prim_t<double>& current_state) -> prim_t<double>
	{
		cons_t<double> actual;
		prim_t<double> output;
		for (auto n:range(0,5)) actual[n[0]] = current_state[n[0]];
		convert_state(actual, output, air);
		return output;
	};
	
	int nfiles = 1+(end-start)/skip;
	
	auto nxb = domain.GetMeshDataDim();
	auto nb  = domain.GetBlockDim();
	int ny = nxb[1]*nb[1];
	
	//To output:
	
	std::vector<avg_qty_t<double>*> registry;
	
	avg_qty_t<double>            y_bar(ny, registry);
	avg_qty_t<double>            U_bar(ny, registry);
	avg_qty_t<double>            V_bar(ny, registry);
	avg_qty_t<double>            W_bar(ny, registry);
	avg_qty_t<double>            P_bar(ny, registry);
	avg_qty_t<double>            T_bar(ny, registry);
	avg_qty_t<double>           mu_bar(ny, registry);
	avg_qty_t<double>          rho_bar(ny, registry);
	avg_qty_t<double>         rhoT_bar(ny, registry);
	avg_qty_t<double>         rhoU_bar(ny, registry);
	avg_qty_t<double>         rhoV_bar(ny, registry);
	avg_qty_t<double>         rhoW_bar(ny, registry);
	avg_qty_t<double>        rhoT2_bar(ny, registry);
	avg_qty_t<double>        rhoU2_bar(ny, registry);
	avg_qty_t<double>        rhoV2_bar(ny, registry);
	avg_qty_t<double>        rhoW2_bar(ny, registry);
	avg_qty_t<double>        rhoUT_bar(ny, registry);
	avg_qty_t<double>        rhoVT_bar(ny, registry);
	avg_qty_t<double>        rhoWT_bar(ny, registry);
	avg_qty_t<double>       rho2T2_bar(ny, registry);
	avg_qty_t<double>       rho2U2_bar(ny, registry);
	avg_qty_t<double>       rho2V2_bar(ny, registry);
	avg_qty_t<double>       rho2W2_bar(ny, registry);
	avg_qty_t<double>           U2_bar(ny, registry);
	avg_qty_t<double>           V2_bar(ny, registry);
	avg_qty_t<double>           W2_bar(ny, registry);
	avg_qty_t<double>           T2_bar(ny, registry);
	avg_qty_t<double>       tau_00_bar(ny, registry);
	avg_qty_t<double>       tau_01_bar(ny, registry);
	avg_qty_t<double>       tau_02_bar(ny, registry);
	avg_qty_t<double>       tau_10_bar(ny, registry);
	avg_qty_t<double>       tau_11_bar(ny, registry);
	avg_qty_t<double>       tau_12_bar(ny, registry);
	avg_qty_t<double>       tau_20_bar(ny, registry);
	avg_qty_t<double>       tau_21_bar(ny, registry);
	avg_qty_t<double>       tau_22_bar(ny, registry);
	avg_qty_t<double>       gru_00_bar(ny, registry);
	avg_qty_t<double>       gru_01_bar(ny, registry);
	avg_qty_t<double>       gru_02_bar(ny, registry);
	avg_qty_t<double>       gru_10_bar(ny, registry);
	avg_qty_t<double>       gru_11_bar(ny, registry);
	avg_qty_t<double>       gru_12_bar(ny, registry);
	avg_qty_t<double>       gru_20_bar(ny, registry);
	avg_qty_t<double>       gru_21_bar(ny, registry);
	avg_qty_t<double>       gru_22_bar(ny, registry);
	avg_qty_t<double>     tau_u_00_bar(ny, registry);
	avg_qty_t<double>     tau_u_01_bar(ny, registry);
	avg_qty_t<double>     tau_u_02_bar(ny, registry);
	avg_qty_t<double>     tau_u_10_bar(ny, registry);
	avg_qty_t<double>     tau_u_11_bar(ny, registry);
	avg_qty_t<double>     tau_u_12_bar(ny, registry);
	avg_qty_t<double>     tau_u_20_bar(ny, registry);
	avg_qty_t<double>     tau_u_21_bar(ny, registry);
	avg_qty_t<double>     tau_u_22_bar(ny, registry);
	avg_qty_t<double>         dpdy_bar(ny, registry);
	avg_qty_t<double>   tau_g_00_0_bar(ny, registry);
	avg_qty_t<double>   tau_g_01_1_bar(ny, registry);
	avg_qty_t<double>   tau_g_02_2_bar(ny, registry);
	avg_qty_t<double>   tau_g_10_0_bar(ny, registry);
	avg_qty_t<double>   tau_g_11_1_bar(ny, registry);
	avg_qty_t<double>   tau_g_12_2_bar(ny, registry);
	avg_qty_t<double>   tau_g_20_0_bar(ny, registry);
	avg_qty_t<double>   tau_g_21_1_bar(ny, registry);
	avg_qty_t<double>   tau_g_22_2_bar(ny, registry);
	avg_qty_t<double>         p_ux_bar(ny, registry);
	avg_qty_t<double>         p_vy_bar(ny, registry);
	avg_qty_t<double>         p_wz_bar(ny, registry);
	avg_qty_t<double>           ux_bar(ny, registry);
	avg_qty_t<double>           vy_bar(ny, registry);
	avg_qty_t<double>           wz_bar(ny, registry);
	
	
	auto accumulate = [](std::vector<double>& a, const std::vector<double>& b) -> void
	{
		for (int i = 0; i < a.size(); i++) a[i] += b[i];
	};
	
	for (auto i: range(0,nfiles))
	{
		std::string filename = strformat("data/restart_unk_nt_{}.dat", ZFill(start+skip*i[0], 8));
		print("Reading", filename);
		reader.LoadData(primsInst, filename);
		
		transform_inplace(primsInst, cons2prim);
		
		double Tref = 100.0;
		double mu_ref = 3e-4;
		double npow = 0.76;
		compute_average(primsInst,        y_bar.inst_data, [=](const cmf::Vec3<double>& xyz) -> double {return xyz[1];});
		compute_average(primsInst,        U_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.u();});
		compute_average(primsInst,        V_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.v();});
		compute_average(primsInst,        W_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.w();});
		compute_average(primsInst,        P_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.p();});
		compute_average(primsInst,        T_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.T();});
		compute_average(primsInst,       mu_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return mu_ref*pow(prim.T()/Tref, npow);});
		compute_average(primsInst,      rho_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return          prim.p()/(air.R*prim.T());});
		compute_average(primsInst,     rhoT_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,     rhoU_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,     rhoV_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,     rhoW_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoT2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.T()*prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoU2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.u()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoV2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.v()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoW2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.w()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoUT_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoVT_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,    rhoWT_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,   rho2T2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return (prim.T()*prim.p()/(air.R*prim.T()))*(prim.T()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst,   rho2U2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return (prim.u()*prim.p()/(air.R*prim.T()))*(prim.u()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst,   rho2V2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return (prim.v()*prim.p()/(air.R*prim.T()))*(prim.v()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst,   rho2W2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return (prim.w()*prim.p()/(air.R*prim.T()))*(prim.w()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst,       U2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.u();});
		compute_average(primsInst,       V2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.v();});
		compute_average(primsInst,       W2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.w();});
		compute_average(primsInst,       T2_bar.inst_data, [=](const prim_t<double>& prim)   -> double {return prim.T()*prim.T();});
		compute_average(primsInst,   tau_00_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[0].u()+prim_grad[0].u()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst,   tau_01_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[1].u()+prim_grad[0].v());});
		compute_average(primsInst,   tau_02_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[2].u()+prim_grad[0].w());});
		compute_average(primsInst,   tau_10_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[0].v()+prim_grad[1].u());});
		compute_average(primsInst,   tau_11_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[1].v()+prim_grad[1].v()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst,   tau_12_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[2].v()+prim_grad[1].w());});
		compute_average(primsInst,   tau_20_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[0].w()+prim_grad[2].u());});
		compute_average(primsInst,   tau_21_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[1].w()+prim_grad[2].v());});
		compute_average(primsInst,   tau_22_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return visc.calc_visc(prim)*(prim_grad[2].w()+prim_grad[2].w()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst, tau_u_00_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].u()*visc.calc_visc(prim)*(prim_grad[0].u()+prim_grad[0].u()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst, tau_u_01_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].u()*visc.calc_visc(prim)*(prim_grad[1].u()+prim_grad[0].v());});
		compute_average(primsInst, tau_u_02_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].u()*visc.calc_visc(prim)*(prim_grad[2].u()+prim_grad[0].w());});
		compute_average(primsInst, tau_u_10_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].v()*visc.calc_visc(prim)*(prim_grad[0].v()+prim_grad[1].u());});
		compute_average(primsInst, tau_u_11_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].v()*visc.calc_visc(prim)*(prim_grad[1].v()+prim_grad[1].v()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst, tau_u_12_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].v()*visc.calc_visc(prim)*(prim_grad[2].v()+prim_grad[1].w());});
		compute_average(primsInst, tau_u_20_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].w()*visc.calc_visc(prim)*(prim_grad[0].w()+prim_grad[2].u());});
		compute_average(primsInst, tau_u_21_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].w()*visc.calc_visc(prim)*(prim_grad[1].w()+prim_grad[2].v());});
		compute_average(primsInst, tau_u_22_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].w()*visc.calc_visc(prim)*(prim_grad[2].w()+prim_grad[2].w()) + visc.calc_beta(prim)*(prim_grad[0].u() + prim_grad[1].v() + prim_grad[2].w());});
		compute_average(primsInst,   gru_00_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].u();});
		compute_average(primsInst,   gru_01_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].u();});
		compute_average(primsInst,   gru_02_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].u();});
		compute_average(primsInst,   gru_10_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].v();});
		compute_average(primsInst,   gru_11_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].v();});
		compute_average(primsInst,   gru_12_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].v();});
		compute_average(primsInst,   gru_20_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].w();});
		compute_average(primsInst,   gru_21_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].w();});
		compute_average(primsInst,   gru_22_bar.inst_data, [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].w();});
		compute_average(primsInst,   dpdy_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].p();});
		compute_average(primsInst,   p_ux_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim.p()*prim_grad[0].u();});
		compute_average(primsInst,   p_vy_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim.p()*prim_grad[1].v();});
		compute_average(primsInst,   p_wz_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim.p()*prim_grad[2].w();});
		compute_average(primsInst,     ux_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[0].u();});
		compute_average(primsInst,     vy_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[1].v();});
		compute_average(primsInst,     wz_bar.inst_data,   [=](const prim_t<double>& prim, const val_grad<3, prim_t<double>>& prim_grad) -> double {return prim_grad[2].w();});
		
		compute_viscdif(primsInst, visc,
			tau_g_00_0_bar.inst_data,
			tau_g_01_1_bar.inst_data,
			tau_g_02_2_bar.inst_data,
			tau_g_10_0_bar.inst_data,
			tau_g_11_1_bar.inst_data,
			tau_g_12_2_bar.inst_data,
			tau_g_20_0_bar.inst_data,
			tau_g_21_1_bar.inst_data,
			tau_g_22_2_bar.inst_data);
		for (auto a:registry) a->accum();
	}
	
	std::vector<double> upp_upp(ny, 0.0);
	std::vector<double> vpp_vpp(ny, 0.0);
	std::vector<double> wpp_wpp(ny, 0.0);
	std::vector<double> Tpp_Tpp(ny, 0.0);
	
	std::vector<double> up_up(ny, 0.0);
	std::vector<double> vp_vp(ny, 0.0);
	std::vector<double> wp_wp(ny, 0.0);
	std::vector<double> Tp_Tp(ny, 0.0);
	
	std::vector<double> u_tilde(ny, 0.0);
	std::vector<double> v_tilde(ny, 0.0);
	std::vector<double> w_tilde(ny, 0.0);
	std::vector<double> T_tilde(ny, 0.0);

	std::vector<double> u2_tilde(ny, 0.0);
	std::vector<double> v2_tilde(ny, 0.0);
	std::vector<double> w2_tilde(ny, 0.0);
	std::vector<double> T2_tilde(ny, 0.0);

	std::vector<double> uT_tilde(ny, 0.0);
	std::vector<double> vT_tilde(ny, 0.0);
	std::vector<double> wT_tilde(ny, 0.0);

	std::vector<double> uT_pp(ny, 0.0);
	std::vector<double> vT_pp(ny, 0.0);
	std::vector<double> wT_pp(ny, 0.0);
	
	std::vector<double> A00(ny, 0.0);
	std::vector<double> A01(ny, 0.0);
	std::vector<double> A02(ny, 0.0);
	std::vector<double> A10(ny, 0.0);
	std::vector<double> A11(ny, 0.0);
	std::vector<double> A12(ny, 0.0);
	std::vector<double> A20(ny, 0.0);
	std::vector<double> A21(ny, 0.0);
	std::vector<double> A22(ny, 0.0);
	std::vector<double> C10(ny, 0.0);
	
	std::vector<double> D0 (ny, 0.0);
	std::vector<double> D1 (ny, 0.0);
	std::vector<double> D2 (ny, 0.0);
	std::vector<double> B00(ny, 0.0);
	std::vector<double> B01(ny, 0.0);
	std::vector<double> B02(ny, 0.0);
	
	comp_fluc_i(A00, tau_u_00_bar.accumulated_data, tau_00_bar.accumulated_data, gru_00_bar.accumulated_data);
	comp_fluc_i(A01, tau_u_01_bar.accumulated_data, tau_01_bar.accumulated_data, gru_01_bar.accumulated_data);
	comp_fluc_i(A02, tau_u_02_bar.accumulated_data, tau_02_bar.accumulated_data, gru_02_bar.accumulated_data);
	comp_fluc_i(A10, tau_u_10_bar.accumulated_data, tau_10_bar.accumulated_data, gru_10_bar.accumulated_data);
	comp_fluc_i(A11, tau_u_11_bar.accumulated_data, tau_11_bar.accumulated_data, gru_11_bar.accumulated_data);
	comp_fluc_i(A12, tau_u_12_bar.accumulated_data, tau_12_bar.accumulated_data, gru_12_bar.accumulated_data);
	comp_fluc_i(A20, tau_u_20_bar.accumulated_data, tau_20_bar.accumulated_data, gru_20_bar.accumulated_data);
	comp_fluc_i(A21, tau_u_21_bar.accumulated_data, tau_21_bar.accumulated_data, gru_21_bar.accumulated_data);
	comp_fluc_i(A22, tau_u_22_bar.accumulated_data, tau_22_bar.accumulated_data, gru_22_bar.accumulated_data);
	
	comp_tild(u_tilde, rhoU_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(v_tilde, rhoV_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(w_tilde, rhoW_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(T_tilde, rhoT_bar.accumulated_data, rho_bar.accumulated_data);

	comp_tild(u2_tilde, rhoU2_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(v2_tilde, rhoV2_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(w2_tilde, rhoW2_bar.accumulated_data, rho_bar.accumulated_data);

	comp_tild(uT_tilde, rhoUT_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(vT_tilde, rhoVT_bar.accumulated_data, rho_bar.accumulated_data);
	comp_tild(wT_tilde, rhoWT_bar.accumulated_data, rho_bar.accumulated_data);
	
	comp_tild(T2_tilde, rhoT2_bar.accumulated_data, rho_bar.accumulated_data);
	
	comp_fluc(upp_upp, U_bar.accumulated_data, rho_bar.accumulated_data, rhoU_bar.accumulated_data, rho2U2_bar.accumulated_data);
	comp_fluc(vpp_vpp, V_bar.accumulated_data, rho_bar.accumulated_data, rhoV_bar.accumulated_data, rho2V2_bar.accumulated_data);
	comp_fluc(wpp_wpp, W_bar.accumulated_data, rho_bar.accumulated_data, rhoW_bar.accumulated_data, rho2W2_bar.accumulated_data);
	comp_fluc(Tpp_Tpp, T_bar.accumulated_data, rho_bar.accumulated_data, rhoT_bar.accumulated_data, rho2T2_bar.accumulated_data);
	
	comp_fluc_i(up_up, U2_bar.accumulated_data, U_bar.accumulated_data, U_bar.accumulated_data);
	comp_fluc_i(vp_vp, V2_bar.accumulated_data, V_bar.accumulated_data, V_bar.accumulated_data);
	comp_fluc_i(wp_wp, W2_bar.accumulated_data, W_bar.accumulated_data, W_bar.accumulated_data);
	comp_fluc_i(Tp_Tp, T2_bar.accumulated_data, T_bar.accumulated_data, T_bar.accumulated_data);

	comp_fluc_i(upp_upp, u2_tilde, u_tilde, u_tilde);
	comp_fluc_i(vpp_vpp, v2_tilde, v_tilde, v_tilde);
	comp_fluc_i(wpp_wpp, w2_tilde, w_tilde, w_tilde);
	comp_fluc_i(Tpp_Tpp, T2_tilde, T_tilde, T_tilde);

	comp_fluc_i(uT_pp, uT_tilde, u_tilde, T_tilde);
	comp_fluc_i(vT_pp, vT_tilde, v_tilde, T_tilde);
	comp_fluc_i(wT_pp, wT_tilde, w_tilde, T_tilde);

	for (std::size_t i = 0; i < C10.size(); i++)
	{
		C10[i] = (U_bar[i] - u_tilde[i]) * dpdy_bar[i];
		D0[i]  = P_bar[i]*ux_bar[i]-p_ux_bar[i];
		D1[i]  = P_bar[i]*vy_bar[i]-p_vy_bar[i];
		D2[i]  = P_bar[i]*wz_bar[i]-p_wz_bar[i];
		B00[i] = (u_tilde[i] - U_bar[i])*tau_g_00_0_bar[i];
		B01[i] = (u_tilde[i] - U_bar[i])*tau_g_01_1_bar[i];
		B02[i] = (u_tilde[i] - U_bar[i])*tau_g_02_2_bar[i];
	}
	

	std::vector<avg_qty_t<double>*> registry2;
	avg_qty_t<double> q1_turb(ny, registry2);
	avg_qty_t<double> q2_turb(ny, registry2);
	avg_qty_t<double> q3_turb(ny, registry2);
	
	
	for (auto i: range(0,nfiles))
	{
		std::string filename = strformat("data/restart_unk_nt_{}.dat", ZFill(start+skip*i[0], 8));
		print("Reading", filename);
		reader.LoadData(primsInst, filename);
		transform_inplace(primsInst, cons2prim);
		
		compute_average(primsInst, q1_turb.inst_data, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.u()-U_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		
		compute_average(primsInst, q2_turb.inst_data, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.v()-V_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		
		compute_average(primsInst, q3_turb.inst_data, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.w()-W_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		for (auto a:registry2) a->accum();
	}
	
	std::vector<std::string> names {"y", "mu_bar", "rho_bar", "U_bar", "V_bar",
									"W_bar", "T_bar", "P_bar", "u_tilde", "v_tilde",
									"w_tilde", "T_tilde", "up_up", "vp_vp", "wp_wp",
									"Tp_Tp", "upp_upp", "vpp_vpp", "wpp_wpp", "Tpp_Tpp",
									"P_bar", "u_tilde", "v_tilde", "w_tilde", "T_tilde",
									"q1_turb", "q2_turb", "q3_turb", "uT_pp", "vT_pp", "wT_pp",
									"A00", "A01", "A02", "A10", "A11", "A12", "A20", "A21", "A22",
									"C10", "dpdy_bar",
									"dtau_00_dx", "dtau_01_dy", "dtau_02_dz",
									"dtau_10_dx", "dtau_11_dy", "dtau_12_dz",
									"dtau_20_dx", "dtau_21_dy", "dtau_22_dz", "D0", "D1", "D2", "B00", "B01", "B02"};
									
	save_csv("output/data.csv", y_bar, mu_bar, rho_bar, U_bar,
								V_bar, W_bar, T_bar, P_bar, u_tilde, v_tilde, w_tilde,
								T_tilde, up_up, vp_vp, wp_wp, Tp_Tp, upp_upp, vpp_vpp,
								wpp_wpp, Tpp_Tpp, P_bar, u_tilde, v_tilde, w_tilde, T_tilde,
								q1_turb, q2_turb, q3_turb, uT_pp, vT_pp, wT_pp,
								A00, A01, A02, A10, A11, A12, A20, A21, A22, C10, dpdy_bar,
								tau_g_00_0_bar, tau_g_01_1_bar, tau_g_02_2_bar, tau_g_10_0_bar,
								tau_g_11_1_bar, tau_g_12_2_bar, tau_g_20_0_bar, tau_g_21_1_bar,
								tau_g_22_2_bar,
								D0, D1, D2, B00, B01, B02);
	save_names("output/names.dat", names);
		
    return 0;
}
