#include "cmf.h"
#include "transform_inplace.h"
#include "range.h"
#include "compute_average.h"
#include "output.h"
#include "tild_fluc.h"

using cmf::print;
using cmf::strformat;
using cmf::ZFill;

int main(int argc, char** argv)
{
	std::string inFile = "input.ptl";
    cmf::ReadInput(inFile);
    cmf::globalSettings = cmf::GlobalSettings(cmf::mainInput["GlobalSettings"]);
	cmf::CreateParallelContext(&argc, &argv);
	
	std::string gid_file = "data/gridInterpolationInfo_01057000.bin";
	std::string rba_file = "data/restart_block_arrangement_nt_01057000.dat";
	
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
	
	auto cons2prim = [=](const prim_t<double>& current_state) -> prim_t<double>
	{
		cons_t<double> actual;
		prim_t<double> output;
		for (auto n:range(0,5)) actual[n[0]] = current_state[n[0]];
		convert_state(actual, output, air);
		return output;
	};
	
	int start = 980000;
	int end   = 1057000;
	int skip = 7000;
	int nfiles = 1+(end-start)/skip;
	
	auto nxb = domain.GetMeshDataDim();
	auto nb  = domain.GetBlockDim();
	int ny = nxb[1]*nb[1];
	
	//To output:
	
	std::vector<double>          y_bar(ny, 0.0);
	std::vector<double>          U_bar(ny, 0.0);
	std::vector<double>          V_bar(ny, 0.0);
	std::vector<double>          W_bar(ny, 0.0);
	std::vector<double>          P_bar(ny, 0.0);
	std::vector<double>          T_bar(ny, 0.0);
	std::vector<double>         mu_bar(ny, 0.0);
	std::vector<double>        rho_bar(ny, 0.0);
	std::vector<double>       rhoT_bar(ny, 0.0);
	std::vector<double>       rhoU_bar(ny, 0.0);
	std::vector<double>       rhoV_bar(ny, 0.0);
	std::vector<double>       rhoW_bar(ny, 0.0);
	std::vector<double>     rho2T2_bar(ny, 0.0);
	std::vector<double>     rho2U2_bar(ny, 0.0);
	std::vector<double>     rho2V2_bar(ny, 0.0);
	std::vector<double>     rho2W2_bar(ny, 0.0);
	std::vector<double>         U2_bar(ny, 0.0);
	std::vector<double>         V2_bar(ny, 0.0);
	std::vector<double>         W2_bar(ny, 0.0);
	std::vector<double>         T2_bar(ny, 0.0);
	
	std::vector<double>      y_bar_loc(ny, 0.0);
	std::vector<double>      U_bar_loc(ny, 0.0);
	std::vector<double>      V_bar_loc(ny, 0.0);
	std::vector<double>      W_bar_loc(ny, 0.0);
	std::vector<double>      P_bar_loc(ny, 0.0);
	std::vector<double>      T_bar_loc(ny, 0.0);
	std::vector<double>     mu_bar_loc(ny, 0.0);
	std::vector<double>    rho_bar_loc(ny, 0.0);
	std::vector<double>   rhoT_bar_loc(ny, 0.0);
	std::vector<double>   rhoU_bar_loc(ny, 0.0);
	std::vector<double>   rhoV_bar_loc(ny, 0.0);
	std::vector<double>   rhoW_bar_loc(ny, 0.0);
	std::vector<double> rho2T2_bar_loc(ny, 0.0);
	std::vector<double> rho2U2_bar_loc(ny, 0.0);
	std::vector<double> rho2V2_bar_loc(ny, 0.0);
	std::vector<double> rho2W2_bar_loc(ny, 0.0);
	std::vector<double>     U2_bar_loc(ny, 0.0);
	std::vector<double>     V2_bar_loc(ny, 0.0);
	std::vector<double>     W2_bar_loc(ny, 0.0);
	std::vector<double>     T2_bar_loc(ny, 0.0);
	
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
		compute_average(primsInst,      y_bar_loc, [=](const cmf::Vec3<double>& xyz) -> double {return xyz[1];});
		compute_average(primsInst,      U_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.u();});
		compute_average(primsInst,      V_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.v();});
		compute_average(primsInst,      W_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.w();});
		compute_average(primsInst,      P_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.p();});
		compute_average(primsInst,      T_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.T();});
		compute_average(primsInst,     mu_bar_loc, [=](const prim_t<double>& prim)   -> double {return mu_ref*pow(prim.T()/Tref, npow);});
		compute_average(primsInst,    rho_bar_loc, [=](const prim_t<double>& prim)   -> double {return          prim.p()/(air.R*prim.T());});
		compute_average(primsInst,   rhoT_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.T()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,   rhoU_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,   rhoV_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst,   rhoW_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.p()/(air.R*prim.T());});
		compute_average(primsInst, rho2T2_bar_loc, [=](const prim_t<double>& prim)   -> double {return (prim.T()*prim.p()/(air.R*prim.T()))*(prim.T()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst, rho2U2_bar_loc, [=](const prim_t<double>& prim)   -> double {return (prim.u()*prim.p()/(air.R*prim.T()))*(prim.u()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst, rho2V2_bar_loc, [=](const prim_t<double>& prim)   -> double {return (prim.v()*prim.p()/(air.R*prim.T()))*(prim.v()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst, rho2W2_bar_loc, [=](const prim_t<double>& prim)   -> double {return (prim.w()*prim.p()/(air.R*prim.T()))*(prim.w()*prim.p()/(air.R*prim.T()));});
		compute_average(primsInst,     U2_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.u()*prim.u();});
		compute_average(primsInst,     V2_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.v()*prim.v();});
		compute_average(primsInst,     W2_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.w()*prim.w();});
		compute_average(primsInst,     T2_bar_loc, [=](const prim_t<double>& prim)   -> double {return prim.T()*prim.T();});
		
		accumulate(     y_bar,      y_bar_loc);
		accumulate(     U_bar,      U_bar_loc);
		accumulate(     V_bar,      V_bar_loc);
		accumulate(     W_bar,      W_bar_loc);
		accumulate(     P_bar,      P_bar_loc);
		accumulate(     T_bar,      T_bar_loc);
		accumulate(    mu_bar,     mu_bar_loc);
		accumulate(   rho_bar,    rho_bar_loc);
		accumulate(  rhoT_bar,   rhoT_bar_loc);
		accumulate(  rhoU_bar,   rhoU_bar_loc);
		accumulate(  rhoV_bar,   rhoV_bar_loc);
		accumulate(  rhoW_bar,   rhoW_bar_loc);
		accumulate(rho2T2_bar, rho2T2_bar_loc);
		accumulate(rho2U2_bar, rho2U2_bar_loc);
		accumulate(rho2V2_bar, rho2V2_bar_loc);
		accumulate(rho2W2_bar, rho2W2_bar_loc);
		accumulate(    U2_bar,     U2_bar_loc);
		accumulate(    V2_bar,     V2_bar_loc);
		accumulate(    W2_bar,     W2_bar_loc);
		accumulate(    T2_bar,     T2_bar_loc);
	}
	
	auto scl = [](double a, std::vector<double>& d) -> void {for (auto& p: d) p*= a;};
	scl(1.0/nfiles,      y_bar);
	scl(1.0/nfiles,      U_bar);
	scl(1.0/nfiles,      V_bar);
	scl(1.0/nfiles,      W_bar);
	scl(1.0/nfiles,      P_bar);
	scl(1.0/nfiles,      T_bar);
	scl(1.0/nfiles,     mu_bar);
	scl(1.0/nfiles,    rho_bar);
	scl(1.0/nfiles,   rhoT_bar);
	scl(1.0/nfiles,   rhoU_bar);
	scl(1.0/nfiles,   rhoV_bar);
	scl(1.0/nfiles,   rhoW_bar);
	scl(1.0/nfiles, rho2T2_bar);
	scl(1.0/nfiles, rho2U2_bar);
	scl(1.0/nfiles, rho2V2_bar);
	scl(1.0/nfiles, rho2W2_bar);
	scl(1.0/nfiles,     U2_bar);
	scl(1.0/nfiles,     V2_bar);
	scl(1.0/nfiles,     W2_bar);
	scl(1.0/nfiles,     T2_bar);
	
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
	
	comp_tild(u_tilde, rhoU_bar, rho_bar);
	comp_tild(v_tilde, rhoV_bar, rho_bar);
	comp_tild(w_tilde, rhoW_bar, rho_bar);
	comp_tild(T_tilde, rhoT_bar, rho_bar);
	
	comp_fluc(upp_upp, U_bar, rho_bar, rhoU_bar, rho2U2_bar);
	comp_fluc(vpp_vpp, V_bar, rho_bar, rhoV_bar, rho2V2_bar);
	comp_fluc(wpp_wpp, W_bar, rho_bar, rhoW_bar, rho2W2_bar);
	comp_fluc(Tpp_Tpp, T_bar, rho_bar, rhoT_bar, rho2T2_bar);
	
	comp_fluc_i(up_up, U_bar, U2_bar);
	comp_fluc_i(vp_vp, V_bar, V2_bar);
	comp_fluc_i(wp_wp, W_bar, W2_bar);
	comp_fluc_i(Tp_Tp, T_bar, T2_bar);
	
	std::vector<std::string> names;
	names.push_back("y");
	names.push_back("mu_bar");
	names.push_back("rho_bar");
	names.push_back("u_bar");
	names.push_back("v_bar");
	names.push_back("w_bar");
	names.push_back("T_bar");
	names.push_back("P_bar");
	names.push_back("u_tilde");
	names.push_back("v_tilde");
	names.push_back("w_tilde");
	names.push_back("T_tilde");
	names.push_back("u\'");
	names.push_back("v\'");
	names.push_back("w\'");
	names.push_back("T\'");
	names.push_back("u\'\'");
	names.push_back("v\'\'");
	names.push_back("w\'\'");
	names.push_back("T\'\'");

	std::vector<double>     q1_turb_loc(ny, 0.0);
	std::vector<double>     q1_turb    (ny, 0.0);
	std::vector<double>     q2_turb_loc(ny, 0.0);
	std::vector<double>     q2_turb    (ny, 0.0);
	std::vector<double>     q3_turb_loc(ny, 0.0);
	std::vector<double>     q3_turb    (ny, 0.0);
	
	for (auto i: range(0,nfiles))
	{
		std::string filename = strformat("data/restart_unk_nt_{}.dat", ZFill(start+skip*i[0], 8));
		print("Reading", filename);
		reader.LoadData(primsInst, filename);
		transform_inplace(primsInst, cons2prim);
		
		compute_average(primsInst, q1_turb_loc, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.u()-U_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		
		compute_average(primsInst, q2_turb_loc, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.v()-V_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		
		compute_average(primsInst, q3_turb_loc, [&](const prim_t<double>& prim, const std::size_t& y_index) -> double
		{
			double dens = prim.p()/(air.R*prim.T());
			return dens*(prim.w()-W_bar[y_index])*(prim.T()-T_bar[y_index]);
		});
		
		accumulate(q1_turb, q1_turb_loc);
		accumulate(q2_turb, q2_turb_loc);
		accumulate(q3_turb, q3_turb_loc);
	}
	
	scl(1.0/nfiles, q1_turb);
	scl(1.0/nfiles, q2_turb);
	scl(1.0/nfiles, q3_turb);
	
	save_csv("output/data.csv",// names,
		y_bar, mu_bar, rho_bar, U_bar, V_bar, W_bar, T_bar, P_bar, 
		u_tilde, v_tilde, w_tilde, T_tilde,
		up_up, vp_vp, wp_wp, Tp_Tp,
		upp_upp, vpp_vpp, wpp_wpp, Tpp_Tpp);
		
	save_csv("output/flowmeandata.csv",// names,
		y_bar, P_bar, u_tilde, v_tilde, w_tilde, T_tilde);
	
	save_csv("output/qt.csv", y_bar, q1_turb, q2_turb, q3_turb);
		
    return 0;
}
