/*****************************************************************************************[Main.cc]
 *
 * GPMC -- Copyright (c) 2017-2022, Kenji Hashimoto (Nagoya University, Japan)
 *

GPMC sources are based on Glucose 3.0 and SharpSAT 12.08.1.
We will follow the permissions and copyrights of them.

-----------------------------------------------------------

 * SharpSAT https://sites.google.com/site/marcthurley/sharpsat

MIT License

Copyright (c) 2019 marcthurley

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

-----------------------------------------------------------

 Glucose -- Copyright (c) 2009, Gilles Audemard, Laurent Simon
				CRIL - Univ. Artois, France
				LRI  - Univ. Paris Sud, France

Glucose sources are based on MiniSat (see below
 copyrights). Permissions and copyrights of
Glucose are exactly the same as Minisat on which it is based on. (see below).

---------------

Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
Copyright (c) 2007-2010, Niklas Sorensson

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 **************************************************************************************************/

#include <errno.h>

#include <signal.h>
#include <zlib.h>

#include "utils/System.h"
#include "utils/ParseUtils.h"
#include "utils/Options.h"
#include "core/Dimacs.h"
#include "preprocessor/Preprocessor.h"
// #include "core/Solver.h"
#include "core/Counter.h"

using namespace Glucose;
using namespace GPMC;

//=================================================================================================
static Solver *counter;

// Terminate by notifying the solver and back out gracefully. This is mainly to have a test-case
// for this feature of the Solver as it may take longer than an immediate call to '_exit()'.
static void SIGINT_interrupt(int signum) {
	printf("c o\n");
	printf("c o *** INTERRUPTED by signal %d ***\n", signum);
	printf("c o Elapsed time %.2lf s\n", cpuTime());
	fflush(stdout);
	counter->interrupt();
}

// Note that '_exit()' rather than 'exit()' has to be used. The reason is that 'exit()' calls
// destructors and may cause deadlocks if a malloc/free function happens to be running (these
// functions are guarded by locks for multithreaded use).
static void SIGINT_exit(int signum) {
	printf("c o\n");
	printf("c o *** INTERRUPTED by signal %d ***\n", signum);
	printf("c o Elapsed time %.2lf s\n", cpuTime());
	fflush(stdout);
	_exit(1); }

static void SetSigAct(Solver *s) {
	counter = s;
	signal(SIGINT,  SIGINT_interrupt);  //  2, SIGINT
	signal(SIGABRT, SIGINT_exit);       //  6, SIGABRT
	signal(SIGSEGV, SIGINT_exit);       // 11, SIGSEGV
	signal(SIGTERM, SIGINT_interrupt);  // 15, SIGTERM
	signal(SIGXCPU, SIGINT_interrupt);
}
//=================================================================================================
static PPMC::Instance::Mode getMode(int opt_mode)
{
	PPMC::Instance::Mode mode;
	switch(opt_mode) {
		case 0:	mode = PPMC::Instance::MC; 	break;
		case 1:	mode = PPMC::Instance::WMC;	break;
		case 2:	mode = PPMC::Instance::PMC; break;
		case 3:	mode = PPMC::Instance::WPMC; break;
		default:	mode = PPMC::Instance::MC; 	break;
	}
	return mode;
}
//=================================================================================================
// Main:


int main(int argc, char** argv)
{

	try {
		double start = realTime();

		setUsageHelp("c USAGE: %s [options] <input-file>\n");
		// setUsageHelp("c USAGE: %s [options] <input-file> <result-output-file>\n\n  where input may be either in plain or gzipped DIMACS.\n");
		// printf("This is MiniSat 2.0 beta\n");

#if defined(__linux__)
		fpu_control_t oldcw, newcw;
		_FPU_GETCW(oldcw); newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE; _FPU_SETCW(newcw);
		printf("c o WARNING: for repeatability, setting FPU to use double precision\n");
#endif
		// Extra options:
		//
		IntOption    verb("MAIN", "verb",   "Verbosity level (0=silent, 1=some).", 1, IntRange(0, 1));
		IntOption    cpu_lim("MAIN", "cpu-lim","Limit on CPU time allowed in seconds.\n", INT32_MAX, IntRange(0, INT32_MAX));
		IntOption    mem_lim("MAIN", "mem-lim","Limit on memory usage in megabytes.\n", INT32_MAX, IntRange(0, INT32_MAX));

		IntOption opt_mode("GPMC -- MAIN", "mode", "Counting mode (0=mc, 1=wmc, 2=pmc, 3=wpmc).", 0, IntRange(0, 3));

		IntOption opt_varlimit("GPMC -- PP", "varlim", "limit on #Vars in Preprocessing.", 150000, IntRange(0, INT32_MAX));
		DoubleOption opt_pptimelimit("GPMC -- PP", "pptimelim", "Time shreshold for deciding if it performs the last step of preprocessing.", 120, DoubleRange(0, true, DBL_MAX, true));
		IntOption opt_ppverb("GPMC -- PP", "ppverb", "Preprocessing verbosity level (0=some, 1=more).", 0, IntRange(0, 1));

		BoolOption opt_td("GPMC -- MAIN", "td", "Tree Decomposition", true);
		IntOption  opt_td_varlim("GPMC -- MAIN", "tdvarlim", "Limit on #Vars in Tree Decomposition", 150000, IntRange(0, INT32_MAX));
		DoubleOption  opt_td_dlim("GPMC -- MAIN", "tddenlim", "Limit on density of graph in Tree Decomposition", 0.10, DoubleRange(0, true, 1, true));
		DoubleOption  opt_td_rlim("GPMC -- MAIN", "tdratiolim", "Limit on ratio (edges/vars) of graph in Tree Decomposition", 30.0, DoubleRange(0, true, 1000, true));
		DoubleOption  opt_td_twvar("GPMC -- MAIN", "twvarlim", "Limit on tw/vars", 0.25, DoubleRange(0, true, 1, true));
		DoubleOption opt_td_to("GPMC -- MAIN", "tdtime", "Time Limit on Tree Decomposition", 0, DoubleRange(0, true, INT32_MAX, true));
		DoubleOption opt_coef("GPMC -- MAIN", "coef", "TDscore coefficient", 100, DoubleRange(0, true, 10000000, true));
		StringOption opt_tdout("GPMC -- MAIN", "tdout", "Outfile for Tree Decomposition", "NULL");

		parseOptions(argc, argv, true);


		PPMC::Preprocessor PP(false, opt_varlimit, opt_pptimelimit, opt_ppverb);
		PPMC::Instance::Mode mode = getMode(opt_mode);

		if(verb) {
			printf("c o "GPMC_VERSION"\n");
			switch(mode) {
			  case PPMC::Instance::MC:		printf("c o Mode: Model Counting\n"); break;
			  case PPMC::Instance::WMC:		printf("c o Mode: Weighted Model Counting\n"); break;
			  case PPMC::Instance::PMC:		printf("c o Mode: Projected Model Counting\n"); break;
			  case PPMC::Instance::WPMC:	printf("c o Mode: Weighted Projected Model Counting\n"); break;
			}
			fflush(stdout);
		}

		// Use signal handlers that forcibly quit until the solver will be able to respond to
		// interrupts:
		signal(SIGINT, SIGINT_exit);
		signal(SIGTERM, SIGINT_exit);  // 15, SIGTERM
		//  signal(SIGXCPU,SIGINT_exit);

		// Set limit on CPU-time:
		if (cpu_lim != INT32_MAX){
			rlimit rl;
			getrlimit(RLIMIT_CPU, &rl);
			if (rl.rlim_max == RLIM_INFINITY || (rlim_t)cpu_lim < rl.rlim_max){
				rl.rlim_cur = cpu_lim;
				if (setrlimit(RLIMIT_CPU, &rl) == -1)
					printf("WARNING! Could not set resource limit: CPU-time.\n");
			} }

		// Set limit on virtual memory:
		if (mem_lim != INT32_MAX){
			rlim_t new_mem_lim = (rlim_t)mem_lim * 1024*1024;
			rlimit rl;
			getrlimit(RLIMIT_AS, &rl);
			if (rl.rlim_max == RLIM_INFINITY || new_mem_lim < rl.rlim_max){
				rl.rlim_cur = new_mem_lim;
				if (setrlimit(RLIMIT_AS, &rl) == -1)
					printf("c o WARNING! Could not set resource limit: Virtual memory.\n");
			} }

		// Loading input
		if (argc == 1) {
			printf("c o Reading from standard input... Use '--help' for help.\n");
			PP.loadFromStdin(mode);
		}
		else {
			printf("c o Reading from the file %s...\n", argv[1]);
			PP.loadFromFile(argv[1], mode);
		}
		printf("c o Reading finished.\n");
		printf("c o Elapsed time %.2lf s\nc o\n", cpuTime());
		fflush(stdout);

		// Simplifying
		printf("c o Simplification starts...\n");
		PP.Simplify();		// Simplify includes SAT solving.
		printf("c o Simplification finished.\n");
		printf("c o Elapsed time %.2lf s\nc o\n", cpuTime());
		fflush(stdout);

		std::vector<int> dists;
		if(opt_td && !PP.ins.unsat && PP.ins.npvars > 0 && PP.ins.vars > 20) {
			printf("c o Tree Decomposition starts... \n");
			TreeDecomposition td = PP.getTD(opt_td_varlim, opt_td_dlim, opt_td_rlim, opt_td_to);
			if(td.numNodes() > 0) {
				int centroid = td.centroid(PP.ins.npvars);
				if(strcmp(opt_tdout, "NULL") != 0) {
					ofstream out(opt_tdout);
					td.toDimacs(out, centroid+1, PP.ins.npvars);
				}
				if((double)td.width()/PP.ins.npvars < opt_td_twvar)
					dists = td.distanceFromCentroid(PP.ins.npvars);
				else
					printf("c o ignore td\n");
			}
			printf("c o Tree Decomposition finished.\n");
			printf("c o Elapsed time %.2lf s\nc o\n", cpuTime());
		}
		fflush(stdout);

		// Counting
		printf("c o GPMC Counting starts... \n");
		if(mode == PPMC::Instance::MC || mode == PPMC::Instance::PMC) {
			// non-weighted
			Counter<mpz_class> S;
			SetSigAct((Solver *)&S);
			S.import(PP.ins);
			S.computeTDScore(dists, opt_coef);
			S.countModels();
			S.printStats();
		} else {
			// weighted
			Counter<mpfr::mpreal> S;
			SetSigAct((Solver *)&S);
			S.import(PP.ins);
			S.computeTDScore(dists);
			S.countModels();
			S.printStats();
		}
		printf("c o GPMC Counting finished.\nc o\n");
		printf("c o CPU time    = %.2lf s\n", cpuTime());
		printf("c o Real time   = %.2lf s\n", realTime()-start);
		double mem_used = memUsedPeak();
		if (mem_used != 0)
			printf("c o Memory used = %.2f MB\n", mem_used);

		fflush(stdout);
		exit(0);

	} catch (const std::invalid_argument&) {
		printf("c o Invalid argument\n");
		exit(1);
	} catch (OutOfMemoryException&){
		printf("c o INDETERMINATE (out of memory)\n");
		exit(1);
	}
}
