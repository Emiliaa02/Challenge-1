#ifndef POINT_H_
#define POINT_H_

#include <vector>

class Point
{
  typedef std::vector<double> coords_type;

protected:
  coords_type x;

public:
  explicit Point (coords_type const & coords): x (coords) {};

  //compute distance to Point p
  double distance (const Point & p) const;

  void print (void) const;

  //all the getters and setters
  std::size_t get_n_dimensions (void) const;
  double get_coord (std::size_t i) const;
  void set_coord (std::size_t i, double val);
  coords_type get_coords (void) const;

  //we define the methods that we need to use in the class Gradientmethod
  double euclidean_norm (void) const;

  Point difference(Point y);

  Point operator*(const double val);
};

#endif /* POINT_H_ */
