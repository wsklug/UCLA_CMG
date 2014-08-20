// -*- C++ -*-

#ifndef _voom_TetMesh_h_
#define _voom_TetMesh_h_
#include "Mesh.h"

#include "TetQuadrature.h"
#include "ShapeTet4CP.h"

namespace voom {
  /*!
    Only the compute routine has been added on to this. This will compute 
    the necessary items in the structor Points
  */
  class TetMesh: public Mesh {
  public:
    //! Constructor
    TetMesh(const Position& nodes, const Mesh::Position& ghNodes, 
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
	    TissueType, essentialBC, naturalBC) 
    {
      _dimension = 3;
      _nDof = 3;
      switch(quadOrder) {
      case 1: _numQuadPoints = 1; break;
      case 2: _numQuadPoints = 4; break;
      case 3: _numQuadPoints = 5; break;
      default:
	{
	  std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	  exit(0);
	}
      }
    }

    //! Destructor
    ~TetMesh() {;}
    
    //! Compute routine
    void Compute();

  };
}
#endif
