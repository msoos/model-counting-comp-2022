#ifndef Preprocessor_h
#define Preprocessor_h

#include <vector>
#include <iostream>
#include <fstream>
#include "mpfr/mpreal.h"

#include "preprocessor/Instance.h"
#include "preprocessor/IFlowCutter.h"
#include "preprocessor/TreeDecomposition.h"

#include "mtl/Vec.h"
#include "utils/System.h"

namespace PPMC {

class Identifier {
public:
	Identifier(int vars) { cidx.resize(2*vars, -1); }

	void identify(Glucose::Lit l1, Glucose::Lit l2);

	bool hasEquiv() { return eqc.size() > 0; }
	std::vector<std::vector<Glucose::Lit>>& getEquivClasses() { return eqc; }
	Glucose::Lit delegateLit(Glucose::Lit l) { return (cidx[toInt(l)] == -1) ? l : eqc[cidx[toInt(l)]][0]; }

private:
	void MergeEquivClasses(int c1, int c2);

	std::vector<std::vector<Glucose::Lit>> eqc;
	std::vector<int> cidx;
};

class Preprocessor {
	friend class Solver;
public:
	Preprocessor(bool outfile, int varlimit, double timelimit, int verbose) : var_limit(varlimit), time_limit(timelimit), outputDimacs(outfile), verbose(verbose) { }

	void loadFromStdin(Instance::Mode mode) { ins.load(std::cin, mode); }
	void loadFromFile(std::string filename, Instance::Mode mode);

	bool Simplify();
	TreeDecomposition getTD(int vlim, double dlim, double rlim, double timeout);

	Instance ins;

private:
	bool SAT_FLE();
	bool Strengthen();
	bool MergeAdjEquivs();
	bool VariableEliminate();
	bool DefVariableEliminate();
	void pickVars(std::vector<int>& vars);
	void pickDefVars(std::vector<int>& vars);
	int ElimVars(const std::vector<Glucose::Var>& del);

	void Compact(const Glucose::vec<Glucose::lbool>& assigns);
	void CompactClauses(const Glucose::vec<Glucose::lbool>& assigns, std::vector<std::vector<Glucose::Lit>>& cls, std::vector<bool>& occurred, int& varnum);
	void RewriteClauses(std::vector<std::vector<Glucose::Lit>>& cls, const std::vector<Glucose::Var>& map);
	void RewriteClauses(std::vector<std::vector<Glucose::Lit>>& cls, const std::vector<Glucose::Lit>& map);
	void Subsume();

	void printCNFInfo(const char* ppname = "", bool initial = false);

	int var_limit;				// if #vars is over var_limit, it does not try EquivEl/VE/DefVE.
	double time_limit;			// if it does not reach the last step (SAT_FLE) withn timelimit, the last step is skipped.
	bool outputDimacs;

	int verbose;
};

inline void Preprocessor::printCNFInfo(const char* ppname, bool initial)
{
	if(initial)
		printf("c o [%-7s] %d vars (%d pvars), %d cls\n", ppname, ins.vars, ins.npvars, ins.clauses.size());
	else
		printf("c o [%-7s] %d vars (%d pvars), %d cls, %d lrnts, %d fvars, elap. %.2lf s\n", ppname, ins.vars, ins.npvars, ins.clauses.size(), ins.learnts.size(), ins.freevars, Glucose::cpuTime());

	if(ins.weighted)
		std::cout << "c o gweight " << ins.gweight << std::endl;

	fflush(stdout);
}

}

#endif /* Preprocessor_h */
