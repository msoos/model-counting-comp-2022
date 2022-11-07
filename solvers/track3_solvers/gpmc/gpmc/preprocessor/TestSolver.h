#ifndef TestSolver_h
#define TestSolver_h

#include <vector>

#include "core/Solver.h"

namespace PPMC {

class TestSolver: public Glucose::Solver {
public:
	TestSolver(int nvars) { newVars(nvars); }
	TestSolver(int nvars,
			std::vector<std::vector<Glucose::Lit>> clauses,
			std::vector<std::vector<Glucose::Lit>> learnts,
			std::vector<Glucose::Lit> assignedLits);

	void addClauseWith(const std::vector<Glucose::Lit>& ps, bool learnt=false);
	void resetClauses(std::vector<std::vector<Glucose::Lit>>& clauses);

	Glucose::lbool Solve() { return solve_(); }
	Glucose::lbool Solve(const Glucose::vec<Glucose::Lit>& assumptions) { budgetOff(); setConfBudget(clauses.size()*10); return solveLimited(assumptions);}

	bool falsifiedBy(Glucose::Lit l);
	bool falsifiedBy(Glucose::Lit l1, Glucose::Lit l2);
	bool falsifiedBy(Glucose::vec<Glucose::Lit>& assump);
	bool FailedLiterals();

	void exportLearnts(std::vector<std::vector<Glucose::Lit>>& learnts);

	void assign(Glucose::Lit l) { enqueue(l); }
	bool bcp() { return propagate() == Glucose::CRef_Undef; }

	Glucose::vec<Glucose::lbool>& getAssigns() { return assigns; }

private:
	void newVars(int nvars);
	void backTo(int pos);

};
}

#endif
