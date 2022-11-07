#include "preprocessor/Instance.h"
#include "mtl/Sort.h"
//#include "core/SolverTypes.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

using namespace PPMC;
using namespace Glucose;
using namespace std;

Instance::Instance() :
vars(0),
mode(MC),
weighted(false),
projected(false),
npvars(0),
freevars(0),
gweight(1),
isSteady(false),
unsat(false)
{}

static inline Lit SignedIntToLit(int signed_int) {
	assert(signed_int != 0);
	int v = abs(signed_int)-1;
	return (signed_int > 0) ? mkLit(v) : ~mkLit(v);
}

static void Tokens(string buf, vector<string>& ret) {
	ret.push_back("");
	for (char c : buf) {
		if (std::isspace(c)) {
			if (!ret.back().empty()) {
				ret.push_back("");
			}
		} else {
			ret.back() += c;
		}
	}
	if (ret.back().empty()) ret.pop_back();
}

void Instance::load(istream& in, Mode mode) {
	this->mode = mode;
	reflectMode(mode);

	string buf;
	vector<Lit> ps;
	vector<string> tokens;
	bool parse_error = false;
	bool typeSpecified = false;
	npvars = 0;
	ispvars.clear();

	while (getline(in, buf)) {
		if (buf.empty()) continue;
		tokens.clear();
		Tokens(buf, tokens);

		if(tokens.size() == 0) continue;
		else if (tokens[0] == "p") {
			if (tokens.size() == 4 && tokens[1] == "cnf") {
				vars = stoi(tokens[2]);
				int ncls = stoi(tokens[3]);
				assigns.resize(vars, l_Undef);
				assignedLits.clear();

				if(projected) {
					if(ispvars.size() < vars) {
						ispvars.resize(vars, false);
					}
				} else {
					ispvars.resize(vars, true);
					npvars = vars;
				}

				if(weighted)
					lit_weights.resize(2*vars, 1);
			}
			else
				cerr << "c c Header Error!" << endl;
		}
		else if (tokens[0][0] == 'c') {
			if(tokens[0] == "c" && tokens[1] == "p") {
				if(tokens[2] == "weight") {	// Read the weight of a literal
					if(weighted && tokens.size() == 6 && tokens.back() == "0") {
						int lit = stoi(tokens[3]);
						Lit l = SignedIntToLit(lit);
						lit_weights[toInt(l)] = tokens[4];
					}
				}
				else if(tokens[2] == "show") { // Read the list of projected vars
					if(projected && tokens.back() == "0") {
						for(int i=3; i<tokens.size()-1; i++) {
							int v = stoi(tokens[i]);
							if(v-1 >= ispvars.size())
								ispvars.resize(v, false);
							ispvars[v-1] = true;
							npvars++;
						}
					}
				}
			}
			else if(tokens[0] == "cr") {
				if(projected && tokens.back() == "0") {
					for(int i=1; i<tokens.size()-1; i++) {
						int v = stoi(tokens[i]);
						if(v-1 >= ispvars.size())
							ispvars.resize(v, false);
						ispvars[v-1] = true;
						npvars++;
					}
				}
			}
		}
		else { // Read clause
			int lit = 0;
			assert(tokens.size() > 0);
			for(auto token : tokens) {
				lit = stoi(token);
				if(lit == 0) break;
				ps.push_back(SignedIntToLit(lit));
			}
			if(lit == 0) {
				if(!addClause(ps)) {
					return;
				}
				ps.clear();
			}
		}
	}

	learnts.clear();
	isSteady = false;
}

bool Instance::addClause(vector<Lit>& ps, bool learnt) {
	sort(ps.begin(), ps.end());

	vec<Lit> oc;
	oc.clear();

	Lit p; int i, j, flag = 0;

	for (i = j = 0, p = lit_Undef; i < ps.size(); i++)
		if (value(ps[i]) == l_True || ps[i] == ~p)
			return true;
		else if (value(ps[i]) != l_False && ps[i] != p)
			ps[j++] = p = ps[i];
	ps.resize(j);

	if (ps.size() == 0) {
		unsat = true;
		return false;
	}
	else if (ps.size() == 1) {
		assigns[var(ps[0])] = lbool(!sign(p));
		assignedLits.push_back(ps[0]);
		// No check by UP here. Preprocessor will do later.
	}
	else {
		clauses.push_back({});
		copy(ps.begin(), ps.end(), back_inserter(clauses.back()));
	}
	return true;
}

void Instance::reflectMode(Mode mode) {
	switch(mode) {
	case MC:
		weighted = false, projected = false; break;
	case WMC:
		weighted = true, projected = false; break;
	case PMC:
		weighted = false, projected = true; break;
	case WPMC:
		weighted = true, projected = true; break;
	}
}

// For Debug
 void Instance::toDimacs()
 {
	using namespace std;

	cout << "p cnf " << vars << " " << clauses.size() << endl;

	if(mode == Instance::MC)
		cout << "c t mc" << endl;
	else if(mode == Instance::WMC)
		cout << "c t wmc" << endl;
	else if(mode == Instance::PMC)
		cout << "c t pmc" << endl;
	else if(mode == Instance::WPMC)
		cout << "c t wpmc" << endl;

	if(projected) {
		cout << "c p show ";
		for (int i=0; i<npvars; i++)
			cout << (i+1) << " ";
		cout << "0" << endl;
	}
	if(weighted) {
		cout.precision(10);
		for (int i=0; i<npvars; i++) {
			cout << "c p weight " << (i+1) << " " << lit_weights[toInt(mkLit(i))] << " 0" << endl;
			cout << "c p weight -" << (i+1) << " " << lit_weights[toInt(~mkLit(i))] << " 0" << endl;
		}
	}

	for(auto c : clauses) {
		for (auto l : c)
			cout << (sign(l) ? "-":"") << (var(l)+1) << " ";
		cout << "0" << endl;
	}
}
