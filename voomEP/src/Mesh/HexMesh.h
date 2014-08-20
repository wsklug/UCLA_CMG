// -*- C++ -*-

#ifndef _voom_HexMesh_h_
#define _voom_HexMesh_h_
#include "Mesh.h"

#include "HexQuadrature.h"
#include "ShapeHex8.h"

namespace voom {
  /*!
    Only the compute routine has been added on to this. This will compute 
    the necessary items in the structor Points
  */
  class HexMesh: public Mesh {
  public:
    //! Constructor
    HexMesh(const Position& nodes, const Mesh::Position& ghNodes, 
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
      case 2: _numQuadPoints = 8; break;
      default:
	{
	  std::cerr << "ERROR: Higher Quad Order Specified" << std::endl;
	  exit(0);
	}
      }
    }

    //! Destructor
    ~HexMesh() {;}
    
    //! Compute routine
    void Compute();

  };
}
#endif
