#pragma once
#include <ilcplex/cplex.h>
#include <stdexcept>

/*--------------------------------------------------------------*
 |  RAII wrappers for                                              |
 |  • CPXENVptr  (environment)                                   |
 |  • CPXLPptr   (problem)                                       |
 *--------------------------------------------------------------*/

class CplexEnv {
    CPXENVptr env_{nullptr};

public:
    CplexEnv() {
        int st = 0;
        env_ = CPXopenCPLEX(&st);
        if (!env_) {
            char err[1024] = {};
            CPXgeterrorstring(nullptr, st, err);
            throw std::runtime_error(
                std::string("CPXopenCPLEX failed: ") + err);
        }
    }
    // non-movable, non-copyable: one owner only
    CplexEnv(const CplexEnv&)            = delete;
    CplexEnv& operator=(const CplexEnv&) = delete;

    ~CplexEnv() {
        if (env_) CPXcloseCPLEX(&env_);
    }

    // implicit conversion when a raw pointer is needed
    operator CPXENVptr() const noexcept { return env_; }
};

class CplexProb {
    CPXLPptr  lp_{nullptr};
    CPXENVptr env_;

public:
    CplexProb(CPXENVptr env, const char* name) : env_(env) {
        int st = 0;
        lp_ = CPXcreateprob(env_, &st, name);
        if (!lp_) {
            char err[1024] = {};
            CPXgeterrorstring(env_, st, err);
            throw std::runtime_error(
                std::string("CPXcreateprob failed: ") + err);
        }
    }
    // non-movable, non-copyable
    CplexProb(const CplexProb&)            = delete;
    CplexProb& operator=(const CplexProb&) = delete;

    ~CplexProb() {
        if (lp_) CPXfreeprob(env_, &lp_);
    }

    operator CPXLPptr() const noexcept { return lp_; }
};
