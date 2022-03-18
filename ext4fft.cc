#include "cmf.h"
#include "transform_inplace.h"
#include "range.h"
#include "compute_average.h"
#include "output.h"
#include "tild_fluc.h"
#include "timing.h"
#include "extract_profile.h"

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
	int nx = nxb[0]*nb[0];
    int ny = nxb[1]*nb[1];
    int nz = nxb[2]*nb[2];
    
	for (auto i: range(0,nfiles))
	{
		std::string filename = strformat("data/restart_unk_nt_{}.dat", ZFill(start+skip*i[0], 8));
		print("Reading", filename);
		reader.LoadData(primsInst, filename);
        transform_inplace(primsInst, cons2prim);
        for (auto k: range(0, nz))
        {
            std::vector<double> profile_X(nx, 0.0);
            std::vector<double> profile_P(nx, 0.0);
            std::vector<double> profile_U(nx, 0.0);
            std::vector<double> profile_V(nx, 0.0);
            std::vector<double> profile_W(nx, 0.0);
            std::vector<double> profile_T(nx, 0.0);
            std::vector<double> profile_R(nx, 0.0);
            
            for (auto j: range(0, 9))
            {
                {using cmf::strformat; using cmf::ZFill;
                    filename = strformat("xpuvwtr_nt{}_j{}_k{}.csv", ZFill(i[0], 3), ZFill(j[0], 3), ZFill(k[0], 3));
                }
                extract_profile(primsInst, j[0], k[0], profile_X, [=](const cmf::Vec3<double>& xyz) -> double {return xyz[0];});
                extract_profile(primsInst, j[0], k[0], profile_P, [=](const prim_t<double>& prim  ) -> double {return prim.p();});
                extract_profile(primsInst, j[0], k[0], profile_U, [=](const prim_t<double>& prim  ) -> double {return prim.u();});
                extract_profile(primsInst, j[0], k[0], profile_V, [=](const prim_t<double>& prim  ) -> double {return prim.v();});
                extract_profile(primsInst, j[0], k[0], profile_W, [=](const prim_t<double>& prim  ) -> double {return prim.w();});
                extract_profile(primsInst, j[0], k[0], profile_T, [=](const prim_t<double>& prim  ) -> double {return prim.T();});
                extract_profile(primsInst, j[0], k[0], profile_R, [=](const prim_t<double>& prim  ) -> double {return prim.p()/(air.R*prim.T());});
                print("output", filename);
                save_csv("output_fftp/" + filename, profile_X, profile_P, profile_U, profile_V, profile_W, profile_T, profile_R);
            }
        }
	}	
    return 0;
}
