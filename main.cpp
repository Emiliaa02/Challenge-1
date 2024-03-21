#include <iostream>
#include <string>
#include <vector>
#include <json.hpp>
#include <fstream>
#include "Point.hpp"
#include "MuparserFun.hpp"
#include "Gradientmethod.hpp"

using json = nlohmann::json;


int main(){
    //I read all the parameters from data.json
    std::ifstream file ("data.json");
    json input = json::parse(file);
    
    std::string function_s = input.value("f" , " ");
    const unsigned int dim = input.value("dim", 1);
    std::vector<std::string> dfunction_s(dim, "");
    for(std::size_t i=0; i < dim; ++i){
        dfunction_s[i]=input.value("df"+std::to_string(i+1), "");
    }      
    const unsigned int max_iter = input.value("max_it", 100);
    const double tol_res = input.value("tol_res",1e-7);
    const double tol_x = input.value("tol_x",1e-7);
    std::vector<double> coords(dim, 0.0);
    Point initial_guess(coords);
    for(std::size_t i=0; i < dim; ++i){
        initial_guess.set_coord(i, input.value("init"+std::to_string(i+1), 2.0));
    }  

    //I create the function and the gradient of the function using the strings just read from the json file
    MuparserFun function(function_s,dim);
    std::vector<MuparserFun> dfunction;
    dfunction.reserve(dim);
    for(std::size_t i=0; i<dim; ++i){
        std::string var = dfunction_s[i];
    	dfunction.emplace_back(var, dim);
    }
    
    std::cout<<"The function to minimize is :"<<function.get_expression()<<std::endl;
    std::cout<<"It's exact gradient is: ";
    std::cout<<dfunction[0].get_expression()<<", ";
    std::cout<<dfunction[1].get_expression()<<std::endl;
    
    //I create a Gradientmethod object in order to solve the problem
    Gradientmethod gm(dim, function, dfunction, max_iter, tol_res, tol_x, initial_guess);

    int flag;
    Point result(initial_guess.get_coords());
    std::cout << "Do you want to calculate the gradient analitical or numerical? (0 for analitical, 1 for numerical)"<<std::endl;
    std::cin >> flag;
    if(flag==0){
        for(std::size_t i=0; i<dim; ++i){
            result.set_coord(i,gm.solve().get_coord(i));
        }
    }
    if(flag==1){
        for(std::size_t i=0; i<dim; ++i){
            result.set_coord(i,gm.solve_numerical().get_coord(i));
        }   
    }

    std::cout<<"The point that minimize our function is :"<<std::endl;
    result.print();

    std::cout<<"The value of the function in that point is :"<< function(result.get_coords()) <<std::endl;
    
    return 0;
}
