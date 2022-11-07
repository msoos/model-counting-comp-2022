#ifndef Instance_h
#define Instance_h

#include <vector>
#include "core/SolverTypes.h"

#include <gmpxx.h>
#include "mpfr/mpreal.h"

namespace PPMC {

class Instance {
public:
	enum Mode { MC, WMC, PMC, WPMC };

	Instance();

	void load			(std::istream& in, Mode mode=MC);
	bool addClause	(std::vector<Glucose::Lit>& lits, bool learnt=false);
	void reflectMode	(Mode mode);

	Glucose::lbool value (Glucose::Var x) const;
	Glucose::lbool value (Glucose::Lit p) const;

	// For Debug
	void toDimacs();

	// CNF formula
	int vars;
	std::vector<std::vector<Glucose::Lit>> clauses;
	std::vector<std::vector<Glucose::Lit>> learnts;

	// Counter Mode
	Mode mode;
	bool weighted;
	bool projected;

	// For WMC/WPMC
	std::vector<mpfr::mpreal> lit_weights;

	// For PMC/WPMC
	int npvars;
	std::vector<bool> ispvars;

	// Instance keeps temporal fixed literals. Preprocessor will eliminate fixed variables.
	std::vector<Glucose::lbool>	assigns;
	std::vector<Glucose::Lit>	assignedLits;

	// additional information
	int freevars;
	mpfr::mpreal gweight;

	// State
	bool isSteady;		// true when no fixed literals, no free variables, all clauses are sorted.
	bool unsat;

};

inline Glucose::lbool Instance::value (Glucose::Var x) const { return assigns[x]; }
inline Glucose::lbool Instance::value (Glucose::Lit p) const { return assigns[var(p)] ^ sign(p); }

}

#endif
