//-*-C++-*-
#ifndef _VOOM_MATH_
#define _VOOM_MATH_
 
#include "voom.h"

namespace voom
{

  inline double sqr(double a) {return (a*a);}
  
  double norm2(const Array1D  & v);

  void invert(const Tensor3D & a, Tensor3D & b);

  void invert(const Tensor2D & a, Tensor2D & b);

  void invert(const tvmet::Matrix<Real,4,4>& M, tvmet::Matrix<Real,4,4>& invM);

  double determinant(const Tensor2D & F);

  double determinant(const Tensor3D & F);

  void tensorProduct(const Vector3D u, const Vector3D v, Tensor3D & T);

  void tensorProduct(const Vector2D u, const Vector2D v, Tensor2D & T);

};

#endif // _VOOM_MATH_
