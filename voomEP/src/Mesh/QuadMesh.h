// -*- C++ -*-

#ifndef _voom_QuadMesh_h_
#define _voom_QuadMesh_h_
#include "Mesh.h"

#include "QuadQuadrature.h"
#include "ShapeQ4.h"

namespace voom {
  /*!
    A Mesh with Linear Quads
  */
  class QuadMesh: public Mesh {
  public:
    //! Constructor
    QuadMesh(const Position& nodes, const Position& ghNodes, 
	     const DiffusionVectors& diffusion,
	     const IDTable& elementConnectivity, const int nGlobalnodes, 
	     const int nGhostnodes, int *globalNodeID, int *ghostNodeID, 
	     int *indexMap, int *dlocalNodeID, int *dghostNodeID, 
	     int *dindexMap, const int quadOrder, const int nNodes,
	     std::set<int>& PurkinjeNodes, std::vector<int>& TissueType, 
	     BC& essentialBC, BC& naturalBC)
      :Mesh(nodes, ghNodes, diffusion, elementConnectivity, nGlobalnodes, 
	    nGhostnodes, globalNodeID, ghostNodeID, indexMap, dlocalNodeID,
	    dghostNodeID, dindexMap, quadOrder, nNodes, PurkinjeNodes,
	    TissueType, essentialBC, naturalBC){
      _dimension = 2;
      _nDof = 2;
      switch(quadOrder) {
      case 1: _numQuadPoints = 1; break;
      case 2: _numQuadPoints = 4; break;
      case 3: _numQuadPoints = 9; break;
      default:
	{
	  std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	  exit(0);
	}
      }
    }

    //! Destructor
    ~QuadMesh() {;}
    
    //! Compute routine
    void Compute();

  };
}
#endif
