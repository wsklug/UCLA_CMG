// -*-C++-*-
//----------------------------------------------------------------------------
// Shankarjee Krishnamoorthi
// 3/25/2011
//
//----------------------------------------------------------------------------
/*!
  A Mixed mesh class for 3D applications. Currently it can handle Tets, Hex 
  and beam Elements. Based on number of nodes per elements the code determines 
  which elements it is

*/
#ifndef _voom_MixedMesh3D_h_
#define _voom_MixedMesh3D_h_

#include "Mesh.h"
#include "HexQuadrature.h"
#include "TetQuadrature.h"
#include "LineQuadrature.h"
#include "ShapeHex8.h"
#include "ShapeTet4CP.h"
#include "ShapeBar.h"

namespace voom {
  typedef struct QuadPointStruct QP;
  class MixedMesh3D: public Mesh {
  public:
    //! Constructor
    MixedMesh3D(const Position& nodes, const Mesh::Position& ghNodes,
		const DiffusionVectors &diffusion,
		const IDTable& elementConnectivity, const int nGlobalnodes,
		const int nGhostnodes, int *globalNodeID, int *ghostNodeID,
		int *indexMap, int *dlocalNodeID, int *dghostNodeID, 
		int *dindexMap, const int quadOrder, const int nNodes,
		std::set<int> &PurkinjeNodes, std::vector<int>& TissueType, 
		BC& essentialBC, BC& naturalBC,
		std::map<int, std::vector<Real> >& orient)
      :Mesh(nodes, ghNodes, diffusion, elementConnectivity, nGlobalnodes,
	    nGhostnodes, globalNodeID, ghostNodeID, indexMap, dlocalNodeID, 
	    dghostNodeID, dindexMap, quadOrder, nNodes, PurkinjeNodes,
	    TissueType, essentialBC, naturalBC), _orient(orient){
      _dimension = 3;
      _nDof      = 3;
      _radius    = 1.;
    }

    //! Destructor
    ~MixedMesh3D() {;}

    //! Compute routine
    void Compute();

    //! getIndices routine overloaded here
    void getIndices(const int qpoint, IDList& indices);

    //! Compute deformation gradient overloaded here due to bar elements
    void computeDeformationGradient();

  private:
    //! For PMJ master node read in angular cosines and length data
    std::map< int, std::vector<Real> > _orient;

    //! A vector of numquadpoints. This varies with each element
    std::vector<int> _numQPoints;
  };
}

#endif
