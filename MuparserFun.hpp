#ifndef MUPARSERFUN_H_
#define MUPARSERFUN_H_

#include <muParser.h>
#include "Point.hpp"
#include <memory>
#include <string>


class MuparserFun
{
private:
  mu::Parser m_parser;
  unsigned int dim;
  std::vector<double> m_var;
  std::string m_expr;

public:
  MuparserFun(const MuparserFun &m)
    : m_parser(m.m_parser), dim(m.dim)
  {
    m_expr = m.m_expr;
    m_parser.SetExpr(m_expr);
    m_var.resize(dim);
    for(std::size_t i=0; i<dim; ++i){
      m_parser.DefineVar("x"+std::to_string(i+1), &m_var[i]);
    }    
  };

  MuparserFun(const std::string &s, unsigned int n) : dim(n), m_expr(s)
  {
    try
      {
        m_var.resize(n);
        m_parser.SetExpr(s);
        for(std::size_t i=0; i<dim; ++i){
          m_parser.DefineVar("x"+std::to_string(i+1), &m_var[i]);
        }
      }
    catch (mu::Parser::exception_type &e)
      {
        std::cerr << e.GetMsg() << std::endl;
      }
  };

  double
  operator()(const std::vector<double> &val)
  {
    if(val.size() != dim){
        return std::numeric_limits<double>::quiet_NaN();
    }
    for(std::size_t i=0; i < dim; ++i){
        m_var[i] = val[i];
    }
    return m_parser.Eval();
  };

  std::string get_expression(){
    return m_expr;
  };

};

#endif 