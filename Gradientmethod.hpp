#ifndef GRADIENTMETHOD_H_
#define GRADIENTMETHOD_H_

#include <functional>
#include <iostream>
#include <limits>
#include <vector>
#include "Point.hpp"
#include <cmath>
#include "MuparserFun.hpp"


class Gradientmethod {
enum class Choice {
    Exponential_Decay,
    Inverse_Decay,
    Approx_line_search
};

private:
    //that are the methods to calculate the parameter alpha
    double Arm_rule(Point m_x);

    double Exp_Decay(Point m_x);

    double Inv_Decay(Point m_x);

    template<Choice choice>
        double calculateAlpha(Point x);

    //in this fields we have the variables of the gradient method
    MuparserFun m_fun;

    std::vector<MuparserFun> m_dfun;
    const unsigned int m_n_max_it;

    const double m_tol_fun;

    const double m_tol_x;

    Point m_init;

    unsigned int m_dim;

    // Variables employed by the solver

    double m_alpha;

    double m_res;

    unsigned int m_iter;

public:
    //constructor of the Gradientmethod
    Gradientmethod(unsigned int dim, MuparserFun fun,
    std::vector<MuparserFun> dfun, const unsigned int n_max_iter, double tol_fun, const double tol_x, const Point & init);

    //methods
    Point solve(); //this method gives to us the solution of the minimization problem

    double get_residuals() const {
        return m_res;
    };

    unsigned int get_iter() const {
        return m_iter;
    };

    std::size_t get_dim() const {
        return m_init.get_n_dimensions();
    };
};

#endif 