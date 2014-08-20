#include "TetMesh.h"

namespace voom{
  /*
    Routine to compute the values of the shapefunctions and its
    derivatives with respect to the position (not isoparametric values) at 
    each gauss point in the mesh.
  */
  void TetMesh::Compute() {
    const int numElems = _connectivity.size();
    const int dimension = 3, nodePerElem = 4;
    Real nds[nodePerElem][dimension]; // Nodal position
    TetQuadrature quad(_quadOrder); // Quad Object

    for( int i = 0; i < numElems; i++) { 
      for( int j = 0 ; j < nodePerElem; j++){
	int id = _indexMap[ _connectivity[i][j] - 1 ];
	for( int k = 0; k < dimension; k++) 
	  if ( id > 0 ) nds[j][k] = _nodes[ id - 1][k];
	  else          nds[j][k] = _ghNodes[ -(id + 1) ][k];
      }
	  
      //We will loop over the Quad Points now 
      for(Quadrature<3>::ConstPointIterator p = quad.begin(); p!= quad.end();
	  p++){
	ShapeTet4 shp(p->coords);
	const Shape<3>::FunctionContainer &  N = shp.functions();
	const Shape<3>::DerivativeContainer & DN = shp.derivatives();
	
	// Computing the Spatial Derivatives of the shape functions
	/*! 
	  We need \f$\frac{\partial N}{\partial x} \f$. Using the chain rule
	  we can compute this as
	  \f$
	  \frac{\partial N}{\partial \xi} \frac{\partial \xi}{\partial x}
	  \f$
	  The second term is the inverse of the jacobian. The jacobian is 
	  computed as
	  \f$
	  \frac{\partial x_a}{\partial \xi_b} = \sum_{m=1}^{3} 
	  \frac{\partial N_m}{\partial \xi_j} x_m^i
	  \f$
	*/
	Tensor3D dxds(0.), invJac(0.);
	for( int a = 0; a < dimension; a++)
	  for( int b = 0; b < dimension; b++) 
	    for( int c = 0; c < nodePerElem; c++)
	      dxds(a,b) += DN[c](b) * nds[c][a];

	// Inverse of the Jacobian
	invert( dxds, invJac );
	const Real J = determinant( dxds );
	std::vector< std::vector< Real > > 
	  dndx(nodePerElem, std::vector<Real> (dimension, 0.) );
	for(int a = 0; a < nodePerElem; a++) 
	  for(int b = 0; b < dimension; b++) 
	    for(int alpha = 0; alpha < dimension; alpha++) 
	      dndx[a][b] += DN[a](alpha)*invJac(alpha,b);
	  
	// Create a QuadPoint Struct
	struct QuadPointStruct qp = { N, dndx, p->weight*J };

	_points.push_back(qp);	
      }// Quad point loop
    } // i loop

    if (_isReferenceConfiguration){
	_isReferenceConfiguration = false;
	_initialPoints = _points;
    }
  } // Compute function
  
}// namespace
  
