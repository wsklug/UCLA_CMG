#include "MixedMesh3D.h"

namespace voom{
  /*
    Routine to compute values of shapefunctions and its derivatives with 
    respect to the position (not isoparametric values) at each point in the 
    mesh.
  */

  void MixedMesh3D::Compute() {
    const int numElems = _connectivity.size();
    int nBeams = 0, nHexs = 0, nTets = 0;
    Quadrature<3> *quad;
    bool isPmjNode = false;
    std::vector<Real> pmjData(4, 0.);
    Real J; // Jacobian
    _numQPoints.clear();
    for(int i = 0; i < numElems; i++) {
      const int nodePerElem = _connectivity[i].size();
      Real nds[ nodePerElem][_dimension];
      for(int j = 0; j < nodePerElem; j++) {
	int id = _indexMap[ _connectivity[i][j] - 1 ];
        for( int k = 0; k < _dimension; k++) {
          if ( id > 0 )
            nds[j][k] = _nodes[ id - 1][k];
          else
            nds[j][k] = _ghNodes[ -(id + 1) ][k];
        }
      }
      // See if first node is in the orient list
      if (_orient.find( _connectivity[i][0] - 1 ) != _orient.end()) {
	isPmjNode = true;
	pmjData = _orient.find(_connectivity[i][0]-1)->second;
      }
      // Hex Element      
      if (nodePerElem == 8) { 
	quad = new HexQuadrature(_quadOrder);	nHexs++;
	switch(_quadOrder) {
	case 1: _numQuadPoints = 1; break;
	case 2: _numQuadPoints = 8; break;
	default:
	  {
	    std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	    exit(0);
	  }
	}
	_numQPoints.push_back( _numQuadPoints );
      }
      // Tet Element
      if (nodePerElem == 4) {
	quad = new TetQuadrature(_quadOrder);	nTets++;
	switch(_quadOrder) {
	case 1: _numQuadPoints = 1; break;
	case 2: _numQuadPoints = 4; break;
	case 3: _numQuadPoints = 5; break;
	default:
	  {
	    std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	    exit(0);
	  }
	}
	_numQPoints.push_back( _numQuadPoints );
      }
      // Beam Element
      if (nodePerElem == 2) {
	int qord = 2*_quadOrder - 1;
	LineQuadrature myQuad(qord);
	switch(qord) {
	case 1: _numQuadPoints = 1; break;
	case 3: _numQuadPoints = 2; break;
	case 5: _numQuadPoints = 3; break;
	default:
	  {
	    std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	    exit(0);
	  }
	}
	_numQPoints.push_back( _numQuadPoints );	  
	Real mag = 0.;
	tvmet::Vector<Real, 3> ang(0.);
	nBeams++;
	for(int j = 0; j < 3; j++) ang(j) = (nds[1][j] - nds[0][j]);
	mag = tvmet::norm2( ang );
	ang /= mag;
	if (isPmjNode) {
	  ang(0) = pmjData[0]; ang(1) = pmjData[1]; ang(2) = pmjData[2];
	  mag = pmjData[3];
	}

	//We will loop over the Quad Points now 
	for(Quadrature<1>::ConstPointIterator p = myQuad.begin(); 
	    p!= myQuad.end(); p++){
	  ShapeBar shp(p->coords);
	  const Shape<1>::FunctionContainer &  Np = shp.functions();
	  const Shape<1>::DerivativeContainer & DNp = shp.derivatives();  
	  Shape<3>::FunctionContainer N;
	  const Real fac = 2./mag;
	  N.resize(2);
	  N[0] = Np[0];	  N[1] = Np[1];
	  std::vector< std::vector< Real > >
	    dndx(nodePerElem, std::vector<Real> (_dimension) );
	  dndx[0][0] = DNp[0](0)*ang(0)*fac; dndx[0][1] = DNp[0](0)*ang(1)*fac; 
	  dndx[0][2] = DNp[0](0)*ang(2)*fac;
	  dndx[1][0] = DNp[1](0)*ang(0)*fac; dndx[1][1] = DNp[1](0)*ang(1)*fac; 
	  dndx[1][2] = DNp[1](0)*ang(2)*fac;
	  // Create a QuadPoint Struct
	  struct QuadPointStruct qp = {N,dndx,p->weight*M_PI*sqr(_radius)*mag};
	  _points.push_back(qp);
	}// Quad point loop
      } else if( nodePerElem == 8 || nodePerElem == 4) {
	//We will loop over the Quad Points now 
	for(Quadrature<3>::ConstPointIterator p = quad->begin(); 
	    p!= quad->end(); p++){
	  Shape<3>* shp;
	  if (nodePerElem == 8) shp = new ShapeHex8(p->coords);
	  else shp = new ShapeTet4(p->coords);
	  const Shape<3>::FunctionContainer &  N = shp->functions();
	  const Shape<3>::DerivativeContainer & DN = shp->derivatives(); 
	  Tensor3D dxds(0.), invJac(0.);
	  for( int a = 0; a < _dimension; a++)
	    for( int b = 0; b < _dimension; b++)
	      for( int c = 0; c < nodePerElem; c++)
		dxds(a,b) += DN[c](b) * nds[c][a];
	  
	  // Inverse of the Jacobian
	  invert( dxds, invJac );
	  const Real J = determinant( dxds );
	  
	  // Spatial Derivative
	  std::vector< std::vector< Real > >
	    dndx(nodePerElem, std::vector<Real> (_dimension, 0.) );
	  for(int a = 0; a < nodePerElem; a++)
	    for(int b = 0; b < _dimension; b++) 
	      for(int alpha = 0; alpha < _dimension; alpha++)
		dndx[a][b] += DN[a](alpha)*invJac(alpha,b);      
	  struct QuadPointStruct qp = {N, dndx, p->weight*J};
	  delete shp;
	  _points.push_back(qp);
	}// Quad point loop	
	delete quad;
      } // nodeperelem =8/4
      // Unknown Element
      else {
	std::cerr << "Unknown Element Type. Mixed mesh currently supports\n";
	std::cerr << "Hex, tets and beam elements\n";
	exit(1);
      }
    } // i loop
    if (_isReferenceConfiguration) {
      _isReferenceConfiguration = false;
      _initialPoints = _points;
      
      std::cout << "Mixed Mesh Statistics \n";
      std::cout << "Number of Beam Elements:  " << nBeams << "\n";
      std::cout << "Number of Hex Elements :  " << nHexs << "\n";
      std::cout << "Number of Tet Elements :  " << nTets << "\n";
    }
  } // Compute routine
  
  /*
    GetQuadPointList overloaded here in mixed mesh derived class. The number 
    of quadpoint is not a constant since this is a mixed mesh. It will vary
    with each element. Hence we now have a vector of quad point
  */
  
  void MixedMesh3D::getIndices(const int qpoint, IDList& indices) {
    int ctr = 0;
    for(int i = 0 ; i < qpoint; i++)
      ctr += _numQPoints[i];
    indices[0] = ctr;
    ctr += _numQPoints[qpoint];
    indices[1] = ctr;
  }
  
  // Computing Deformation Gradient
  void MixedMesh3D::computeDeformationGradient() {
    _F.clear();
    Tensor3D I( tvmet::identity<Tensor3D>() );
    int numElems = _connectivity.size();
    IDList index(2,0.);
    for(int i = 0; i < numElems; i++) {
      const int nodePerElem = _connectivity[i].size();
      Real nds[nodePerElem][_dimension];
      for(int j = 0; j < nodePerElem; j++) {
	int id = _indexMap[ _connectivity[i][j] - 1 ];
        for( int k = 0; k < _dimension; k++) {
          if ( id > 0 ) nds[j][k] = _nodes[ id - 1][k];
          else nds[j][k] = _ghNodes[ -(id + 1) ][k];
        }
      } // j loop
      
      this->getIndices( i, index );
      // Looping over the Quad Points
      for(int j = index[0]; j < index[1]; j++) {
	Tensor3D F(0.);
	if ( (nodePerElem > 2) && ( nodePerElem <=8 ) ) 
	  for(int k = 0; k < nodePerElem; k++)
	    for(int m = 0; m < _dimension; m++) 
	      for(int n = 0; n < _dimension; n++)
		F(m,n) += nds[k][m]*_initialPoints[j].shapeDerivatives[k][n];
	else F = I; // Identity
	_F.push_back( F );
      } // Quad point loop - j loop
    }// Number of elements loop
  } // end of routine
}
