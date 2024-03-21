
The purpose of the code is to solve a minimization problem using the gradient method. To accomplish this, I've designed a class called GradientMethod, within which I have implemented the method for solving the problem.

Additionally, within this class, there are three different approaches to compute the step size, denoted as alpha, to usein the solve method: the Armijo rule, Exponential Decay, or Inverse Decay.

To avoid the need for an if statement inside the loop, I've created a function template with an enumerator as a parameter. At compile time, I select the choice using if constexpr.

There are also two other classes in the code:

1) The Point.hpp class enables the creation of Point objects, where coordinates are defined as vectors of doubles.
1) The MuParserFun.hpp class allows for the parsing of mathematical expressions written as strings.

In the main function, after reading all the parameters from a JSON file and creating the MuParserFun from the obtained strings, an object of the GradientMethod class is created. Through this object, the solution is found.

TO RUN THE CODE:
After changing the path in the Makefile, you have to execute in the terminal the following commands:

LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/host/folder/pacs-examples/Examples/lib
make
./main 
