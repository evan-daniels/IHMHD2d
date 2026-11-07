// #include "centpack_2d_FD2.h"
#include "centpack_2d_SD2.h"

using namespace std;

// ORIGINAL (pre-2026-11-06):
// void CENTPACK::run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p)
// NEW (2026-11-06): + parameters, to record the scheme settings
void CENTPACK::run_info_2d(double& dt, double& sum_t, long& J, long& K, double& cfl, const int& id, const int& p, const doublearray1d& parameters)
{
	ofstream OutFile;
	
	OutFile.open("run_info.txt");
	OutFile<<"dt = "<< dt <<endl;
	OutFile<<"total time = "<< sum_t <<endl;
	OutFile<<"grid size = "<< p*J <<" x "<< K << endl;
	OutFile<<"CFL = "<< cfl <<endl;
	OutFile<<"Number of processors: "<< p <<endl;
	// NEW (2026-11-06): scheme settings (IHMHD2d)
	{
		const long n = parameters.getIndex1Size();
		const double de = (n > 8) ? parameters(8) : 0.0;
		OutFile<<"d_e = "<< de << ((de > 0.0) ? "  (Formulation B, electron inertia)" : "  (MHD2d step)") <<endl;
		OutFile<<"eta_H = "<< ((n > 6) ? parameters(6) : 0.0) <<endl;
		if (de > 0.0)
			OutFile<<"c_whistler_inertia = "<< ((n > 9) ? parameters(9) : 0.5) << ((n > 9) ? "" : "  (default; par_count < 10)") <<endl;
		else
			OutFile<<"c_whistler_inertia = not used (d_e = 0)"<<endl;
		if (n > 7 && parameters(7) != 0.0)
			OutFile<<"c_whistler (parameters(7)) = "<< parameters(7) <<"  IGNORED: the uncapped c_whistler option was removed"<<endl;
	}
	OutFile<<"Processor "<< p - 1 <<" finished"<<endl;
	OutFile.close();
}
