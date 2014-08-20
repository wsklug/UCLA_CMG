// -*- C++ -*-

#ifndef _voom_Tri3Mesh_h_
#define _voom_Tri3Mesh_h_
#include "Mesh.h"

#include "TriangleQuadrature.h"
#include "ShapeTri3.h"

namespace voom {
  /*!
    A Mesh with Linear Triangles
  */
  class Tri3Mesh: public Mesh {
  public:
    //! Constructor
    Tri3Mesh(const Position& nodes, const Position& ghNodes, 
	     const DiffusionVectors& diffusion,
	     const IDTable& elementConnectivity, const int nGlobalnodes, 
	     const int nGhostnodes, int *globalNodeID, int *ghostNodeID, 
	     int *indexMap, int *dlocalNodeID, int *dghostNodeID, 
	     int *dindexMap, const int quadOrder, const int nNodes,
	     std::set<int> & PurkinjeNodes, std::vector<int>& TissueType, 
	     BC& essentialBC, BC& naturalBC)
      :Mesh(nodes, ghNodes, diffusion, elementConnectivity, nGlobalnodes, 
	    nGhostnodes, globalNodeID, ghostNodeID, indexMap, dlocalNodeID,
	    dghostNodeID, dindexMap, quadOrder, nNodes, PurkinjeNodes,
	    TissueType, essentialBC, naturalBC){
      _dimension = 2;
      _nDof = 2;
      switch(quadOrder) {
      case 1: _numQuadPoints = 1; break;
      case 2: _numQuadPoints = 3; break;
      case 3: _numQuadPoints = 4; break;
      case 4: _numQuadPoints = 7; break;
      default: 
	{
	  std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	  exit(0);
	}
      }
    }

    //! Destructor
    ~Tri3Mesh() {;}
    
    //! Compute routine
    void Compute();

  };
}
#endif
