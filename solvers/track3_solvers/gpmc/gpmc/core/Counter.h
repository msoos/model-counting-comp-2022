#ifndef Counter_h
#define Counter_h

#include "core/Solver.h"
#include "core/Component.h"
#include "preprocessor/Instance.h"

using namespace Glucose;

namespace GPMC {

#define GPMC_VERSION "gpmc-1.1"

enum btStateT {
	EXIT, RESOLVED, GO_TO_NEXT_COMP
};

template <class T_data>
class Counter: public Solver {
public:
	// Constructor/Destructor
	//
	Counter(bool weighted=false);
	~Counter();

	// Methods
	//
//	bool presimplify();			// Simplification before model counting
	bool simplify();
	bool countModels();		    // Main count method

//	void registerAsPVar(Var v, bool b);	// Register a projection variable (if b is true)
//	void registerAllVarsAsPVar(int nvars);
//	void registerPVars(vec<Var>& vars);
	int nPVars() const;			// The number of non-isolated projection variables
	int nIsoPVars() const;			// The number of isolated projection variables

	void import(const PPMC::Instance &ins);
	void computeTDScore(const vector<int>& dists, double coef=100.0);

	void initWeights(int nlits);
	void setWeight(Lit l, double weight);

	void printProblemStats(double parsed_time, const char* msg) const;
	void printStats() const;
	void toDimacsRaw(const char *file);	// For debug

	// Data members
	//
	bool sat;
	T_data npmodels;

	// options
	bool on_bj;
	double bjthd;
	bool on_simp;

	int verbosity_c;
	bool mc;
	bool wc;

	/// statistics
	uint64_t conflicts_pre, decisions_pre, propagations_pre;
	uint64_t conflicts_sg, decisions_sg, propagations_sg;
	uint64_t sats, nbackjumps, nbackjumps_sp;
	uint64_t reduce_dbs_pre, simp_dbs;
	double real_stime;
	double simplify_time;

protected:
	// Methods
	//
#if 0
	/// Preprocessing
	bool FailedLiterals();
	void Compact();
	void CompactClauses(vec<CRef>& cs, vec<bool>& occurred, int& varnum);
	void RewriteClauses(const vec<CRef>& cs, const vec<Var>& map);
	void removeClauseNoDetach(CRef cr);
	void reinit(const int varnum);
#endif
	/// Load weight
	void loadWeight(const PPMC::Instance& ins);

	/// Count main method
	void count_main();

	/// Backtrack
	btStateT backtrack(int backtrack_level = -1, Lit lit = lit_Undef, CRef cr = CRef_Undef);
	void cancelCurDL();

	/// SAT solving for the current component
	lbool solveSAT();
	lbool searchBelow(int start_level);

	bool analyzeMC(CRef confl, vec<Lit>& out_learnt, vec<Lit> & selectors,
			int& out_btlevel, unsigned int &nblevels,
			unsigned int &szWithoutSelectors);
	unsigned int computeLBDMC(const vec<Lit> & lits, int end);
	unsigned int computeLBDMC(const Clause &c);

	/// Print stats of the component manager.
	void printStatsOfCM() const;

	// Data members
	//
	vec<char> ispvar;		// This is used only when preprocessing.
	int npvars;				// boundary number between projection vars and non-projection vars
	int npvars_isolated;	// Number of isolated projection vars

	ComponentManager<T_data> cmpmgr;	// The manger for processing components
	vec<vec<int>> occ_lists;

	int limlevel;
	vec<Lit> unitcls;			// The list of learnt unit clauses

	bool last_suc;
	int  last_bklevel;
	CRef last_cr;
	Lit  last_lit;

	vec<T_data> lit_weight;
	T_data gweight;

	vec<double> tdscore;
};

// Inline methods
#if 0
template <typename T_data>
inline void Counter<T_data>::registerAsPVar(Var v, bool b) {
	if (v >= ispvar.size())
		ispvar.growTo(v + 1, 0);
	ispvar[v] = b;
	npvars += b;
}
template <typename T_data>
inline void Counter<T_data>::registerAllVarsAsPVar(int nvars) {
	ispvar.clear();
	ispvar.growTo(nvars+1, true);
	ispvar[nvars] = 0;
	npvars = nvars;
}
template <typename T_data>
inline void Counter<T_data>::registerPVars(vec<Var>& vars) {
	ispvar.clear();
	ispvar.growTo(nVars()+1, false);
	for(int i=0; i<vars.size(); i++) {
		assert(vars[i] < nVars());
		registerAsPVar(vars[i], true);
	}
}
#endif
template <typename T_data>
inline int Counter<T_data>::nPVars() const {
	return npvars;
}
template <typename T_data>
inline int Counter<T_data>::nIsoPVars() const {
	return npvars_isolated;
}
#if 0
template <typename T_data>
inline void Counter<T_data>::initWeights(int nlits) {
	lit_weight.clear();
	lit_weight.growTo(nlits, 1);
}
template <typename T_data>
inline void Counter<T_data>::setWeight(Lit l, double weight) {
	lit_weight[toInt(l)] = weight;
}
#endif
template <typename T_data>
inline void Counter<T_data>::printProblemStats(double time, const char *msg) const {
	printf("c o Problem: %d vars (%d pvars), %d clauses, %"PRIu64" literals, %s %.2f s\n", nVars(), nPVars(), nClauses(), clauses_literals, msg, time);
	fflush(stdout);
}
template <typename T_data>
inline void Counter<T_data>::printStatsOfCM() const {
	cmpmgr.printStats();
}

template <>
void Counter<mpz_class>::loadWeight(const PPMC::Instance &ins);
template <>
void Counter<mpfr::mpreal>::loadWeight(const PPMC::Instance &ins);

}

#endif
