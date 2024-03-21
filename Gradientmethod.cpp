#include "Gradientmethod.hpp"

//I create a template to show the different possible choices for the computation of the parameter alpha 
template<Gradientmethod::Choice choice>
    double Gradientmethod::calculateAlpha(Point x, Point grad) {
    if constexpr (choice == Choice::Exponential_Decay) {
        return Gradientmethod::Exp_Decay(x);
    } else if constexpr (choice == Choice::Inverse_Decay) {
        return Gradientmethod::Inv_Decay(x);
    } else if constexpr (choice == Choice::Approx_line_search) {
        return Gradientmethod::Arm_rule(x, grad);
    } else {
        throw std::invalid_argument("Invalid choice");
    }
}   

//definition of the constructor
Gradientmethod::Gradientmethod(unsigned int dim, MuparserFun fun, std::vector<MuparserFun> dfun, 
                                const unsigned int n_max_it,const double tol_fun, const double tol_x, const Point & init):
                                m_fun(fun), m_n_max_it(n_max_it), m_tol_fun(tol_fun), m_tol_x(tol_x), m_init(init.get_coords()),m_dim(dim), m_alpha(0.1), m_res(0.0),m_iter(0) {
                                    m_dfun.reserve(m_dim);
                                    for(std::size_t i=0; i < m_dim; ++i){
                                        m_dfun.emplace_back(dfun[i]);
                                    }
                                }; 


//implementation of the method solve with the analitical gradient
Point Gradientmethod::solve(){
    //we initialize the variables
    Point m_x_old(m_init.get_coords());
    Point m_x(m_init.get_coords());
    Point grad(m_init.get_coords());
    double m_dis = 0.0;
    bool flag = true;
    
    //I solve the problem using the Approximation line search
    constexpr Choice choice = Choice::Approx_line_search;

    while(m_iter < m_n_max_it && flag){
        for(std::size_t i=0; i < m_dim; ++i){
            grad.set_coord(i, m_dfun[i](m_x.get_coords()));
        }
        m_alpha = Gradientmethod::calculateAlpha<choice>(m_x, grad);

        //update of the solution
        for(std::size_t i=0; i<m_dim; ++i){
            m_x.set_coord(i, m_x_old.get_coord(i) - m_alpha * grad.get_coord(i));
        }
        m_dis = m_x.distance(m_x_old);
        m_res = grad.euclidean_norm();
        flag = m_dis >= m_tol_x &&  m_res >= m_tol_fun; //check if the condition are satisfied or not

        //update of m_x_old
        for(std::size_t i=0; i<m_dim; ++i){
            m_x_old.set_coord(i, m_x.get_coord(i));
        }

        m_iter++;
    }

    return m_x;
}

//implementation of the method solve with the numerical gradient
Point Gradientmethod::solve_numerical(){
    //we initialize the variables
    Point m_x_old(m_init.get_coords());
    Point m_x(m_init.get_coords());
    Point grad(m_init.get_coords());
    double m_dis = 0.0;
    bool flag = true;
    double h = 0.05;
    std::vector<double> val(m_dim);
    
    //I solve the problem using the Approximation line search
    constexpr Choice choice = Choice::Approx_line_search;

    while(m_iter < m_n_max_it && flag){
        for(std::size_t i=0; i < m_dim; ++i){
            Point m_plus(m_x.get_coords());
            Point m_minus(m_x.get_coords());
            m_plus.set_coord(i,m_plus.get_coord(i)+h);
            m_minus.set_coord(i,m_minus.get_coord(i)-h);
            grad.set_coord(i, (m_fun(m_plus.get_coords())-m_fun(m_minus.get_coords()))/(2*h));
        }
        m_alpha = Gradientmethod::calculateAlpha<choice>(m_x, grad);

        //update of the solution
        for(std::size_t i=0; i<m_dim; ++i){
            m_x.set_coord(i, m_x_old.get_coord(i) - m_alpha * grad.get_coord(i));
        }
        m_dis = m_x.distance(m_x_old);
        m_res = grad.euclidean_norm();
        flag = m_dis >= m_tol_x &&  m_res >= m_tol_fun; //check if the condition are satisfied or not

        //update of m_x_old
        for(std::size_t i=0; i<m_dim; ++i){
            m_x_old.set_coord(i, m_x.get_coord(i));
        }

        m_iter++;
    }

    return m_x;
}

double Gradientmethod::Arm_rule(Point m_x,Point grad){
    double m_sigma = 0.25;
    unsigned int cont=0; //to check that the maximum number of iterations is not reached
    double cond;
    double fval;
    fval = m_fun(m_x.get_coords());

    Point diff = m_x.difference(grad*m_alpha);
    cond = fval - m_fun(diff.get_coords());

    while(cond < m_sigma*m_alpha*grad.euclidean_norm()*grad.euclidean_norm() && cont < m_n_max_it){
        m_alpha = m_alpha / 2;
        cont++;
    }
    return m_alpha;
}

double Gradientmethod::Exp_Decay(Point m_x){
    double mu = 0.2;
    m_alpha = m_alpha * exp(-mu*m_iter);
    return m_alpha;
}

double Gradientmethod::Inv_Decay(Point m_x){
    double mu = 0.2;
    m_alpha = m_alpha / (1+ mu * m_iter);
    return m_alpha;
}

