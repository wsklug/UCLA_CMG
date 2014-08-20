// -*- C++ -*-
// A Mesh Class

#ifndef _voom_Mesh_h_
#define _voom_Mesh_h_
#include "voom.h"
#include <fstream>
#include <set>
#include <sys/stat.h>
#include "VoomMath.h"

#include "Epetra_ConfigDefs.h"
#include "Epetra_Vector.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include <vector>
#include <map>

namespace voom {
  //! Meshtype definition
  enum MeshType {HEXMESH, TETMESH, QUADMESH, TRI3MESH, MESHFREE2D, MESHFREE3D,
		 MIXEDMESH3D};

  /*! 
    \brief
    A base mesh class. This containes nodal positions, connectivity table in 
    the case of finite elements or neighbor list in the case of mesh free.
    Derived classes will perform the activity of computing the quadpoint 
    structure
    
    In the main routine the user needs to call the static function onl
    Example
    
    #include "Mesh.h"
    ...
    ...
    Mesh* myMesh = Mesh::SetMesh(nodes, connectivity, myQuadOrder, HEXMESH);
    myMesh->Compute(); // Computing N and dN/dx
    Mesh::QuadPoints Quadvals = myMesh->getPoints();
    
  */
  class Mesh {
  public:
    struct QuadPointStruct {
      /*! 
	\brief
	The Quad Point Struct holds, shapfunction values evaluated at the 
	quadrature points and the shape function derivatives (not with respect 
	to isoparametric coordinates) also evaluated at the quadrature points 
	and the weight associated with each quadrature point.
	
      */
      std::vector < Real > shapeFunctions;
      std::vector < std::vector < Real > > shapeDerivatives;
      Real weight;
    };

    struct BoundaryCondition {
      /*!
	\brief
	This struct is used to defines values for essential and natural 
	boundary conditions. The first vector defines the nodal ID where 
	these conditions are enforced, second vector defines the DOF and the 
	last vector defines the value of the constraint
      */
      std::vector<int> nodalID;
      std::vector<int> nodalDOF;
      std::vector<Real> values;
    };

    //! To store list of ints
    typedef std::vector< int > IDList;

    //! To form ID Table
    typedef std::vector< IDList > IDTable;

    //! QuadValues. Stored as a 1D array for fast access
    typedef std::vector< struct QuadPointStruct > QuadPoints;

    //! List of boundary conditions
    typedef struct BoundaryCondition BC;

    //! Nodal Postion 2D/3D
    typedef std::vector< std::vector< Real > > Position;

    //! Diffusion Vector
    typedef std::vector< std::vector< Real > > DiffusionVectors;

    //! Vector of Deformation Gradient
    typedef std::vector< Tensor3D > DefGradient;

    //! Constructor
    Mesh(const Position& nodes, const Position& ghNodes, 
    	 const DiffusionVectors& diffusion, const IDTable& elementConnectivity,
	 const int nLocalnodes, const int nGhostnodes, int *localNodeID, 
	 int *ghostNodeID, int *indexMap, int *dlocalNodeID,
	 int *dghostNodeID, int *dindexMap, const int quadOrder, 
	 const int nNodes, std::set<int>& PurkinjeNodes, 
	 std::vector<int>& TissueType, BC& essentialBC, BC& naturalBC);
    
    //! Destructor
    ~Mesh() {
      delete _localNodeID; delete _ghostNodeID; delete _indexMap;
      delete _dlocalNodeID; delete _dghostNodeID; delete _dindexMap;
    }

    //! Compute Virtual function
    virtual void Compute() = 0;

    //! Update Shape function based on displacement
    void Update(){
      _points.clear();
      this->Compute();
    }

    //! Set Bar Element radius
    void setRadius(const Real radius) {
      _radius = radius;
    }

    //! Set Constraint element length
    void setLength(const Real length) {
      _length = length;
    }

    //! Dimensionality of the Mesh 1D/2D/3D
    int getDimension() const { return _dimension; }

    //! Get All QuadPoint Struct
    const QuadPoints& getPoints() { return _points; }

    //! Get Reference Configuration QuadPoint Struct
    const QuadPoints& getReferenceConfigPoints() { return _initialPoints; }

    //! Get indices from the 1D array relevant to qpoint
    virtual void getIndices(const int qpoint, IDList& indices);

    //! Get Nodes
    const Position& getNodes() { return _nodes; }

    //! Get ghostnode position
    const Position& getGhostNodes() { return _ghNodes; }

    //! Get size of struct
    int QuadSize() { return _points.size(); }

    //! Number of Nodes in the model
    int numberOfNodes() { return _nNodes; }

    // Get ID List
    const IDList& getIDList(const int qpoint) const { 
      return _connectivity[qpoint];
    }

    //! Get entire connectivity table
    const IDTable& getIDTable() const { return _connectivity; }

    //! Get diffusion vectors Vector of Vector
    const DiffusionVectors& getDiffusionVectors() const { return _diffusion; }

    //! Get IDTable size
    int getIDTableSize() { return _connectivity.size(); }

    //! Is model isotropic check
    bool isIsotropic() const { return _isIsotropic; }

    //! Get number of Locall nodes
    int getNumLocalNodes() { return _nLocalNodes; }

    //! Get number of Ghost nodes
    int getNumGhostNodes() { return _nGhostNodes; }

    //! Get local node IDs
    int* getLocalNodeID() const { return _localNodeID; }

    //! Get ghost node IDs
    int* getGhostNodeID() const { return _ghostNodeID; }

    //! Get index Map
    int* getIndexMap() const { return _indexMap; }

    //! Get Displacement LocalNodeID
    int* getDisplacementLocalNodeID() const { return _dlocalNodeID; }

    //! Get Displacement GhostNodeID
    int *getDisplacementGhostNodeID() const { return _dghostNodeID; }

    //! Get displacement IndexMap
    int *getDisplacementIndexMap() const { return _dindexMap; }

    //! Get Number of Restraint
    int getNumberOfRestraints() const { return _essentialBC.values.size(); }

    //! Get Purkinje Node Set
    const std::set<int>& getPurkinjeNodes() { return _PurkinjeNodes; }

    //! Get TissueType definition
    const std::vector<int>& getTissueType() { return _TissueType; }

    //! Smooth mesh. Virtual function. Overwritten in derived class
    virtual void smooth(Mesh* smoothMesh) {;}

    //! compute deformation gradient at all Quad Points
    virtual void computeDeformationGradient();

    //! Get deformation gradient
    const DefGradient& getDeformationGradient() {
      return _F;
    }

    //! Update Nodal position. Input is displacement
    void updateNodalPosition(Epetra_Vector* Displacement,
			     Epetra_Vector* ghostDisplacement);

    //! Get DOF of the mesh object
    int getDOF() const { return _nDof; }

    //! Get restraint nodal ID
    const std::vector<int>& getRestraintNodes() { 
      return _essentialBC.nodalID; 
    }

    //! Get restraint DOF
    const std::vector<int>& getRestraintDOF() { 
      return _essentialBC.nodalDOF;
    }

    //! Get restraint values
    const std::vector<Real>& getRestraintValues() { 
      return _essentialBC.values;
    }

    //! Get natural BC nodal ID
    const std::vector<int>&  getForceNodes() {
      return _naturalBC.nodalID;
    }

    //! Get natural BC DOF
    const std::vector<int>& getForceDOF() {
      return _naturalBC.nodalDOF;
    }

    //! Get natural BC values
    const std::vector<Real>& getForceValues() {
      return _naturalBC.values;
    }

    //! Interpolate Scalar to integration point
    std::vector<Real> interpolateScalar(std::vector<Real>& gamma,
					std::map<int, int>& nodeID);
    

    //! Get list of PMJ elements
    const IDList& getPMJElements() { return _pmjList; }

    //! Gateway to the Mesh class. User needs to call this function
    static Mesh* New(Epetra_MpiComm* mpicomm, std::string modelName, 
		     const int quadOrder, MeshType Type);

  protected:
    //! Nodal Positions
    Position     _nodes;

    //! Ghost Nodal Position
    Position     _ghNodes;

    //! Connectivity Table
    IDTable      _connectivity;

    //! Vector of QuadPoint Struct
    QuadPoints   _points;

    //! Vector of Reference configuration Shape function and Derivatives
    QuadPoints   _initialPoints;

    //! Quadrature Order
    int          _quadOrder;

    //! Total number of nodes in the entire model
    int          _nNodes;

    //! Isotropic or Anisoropic Modes
    bool         _isIsotropic;

    //! Number of QuadPoints for a given quadorder. To be set be derived class
    int          _numQuadPoints;

    //! Dimension of the Mesh
    int          _dimension;

    //! This data from the Mesh class will be used Ep_Solver class
    //! Number of local nodes
    int          _nLocalNodes;
    
    //! Number of ghost nodes
    int          _nGhostNodes;

    //! Local nodal ID
    int         *_localNodeID;

    //! Ghost nodal ID
    int         *_ghostNodeID;

    //! Index map. Map of local nodal ID's
    int         *_indexMap;

    /*! 
      Data from mesh class used in EPSolver class for displacement data
      For displacement we need to store nDof*numNodes sized array. Hence 
      we need to create localNodeID, ghostNodeID and indexMap
    */
    //! Local nodal ID for displacement
    int         *_dlocalNodeID;

    //! Ghost nodal ID for displacement
    int         *_dghostNodeID;

    //! Index map. Map of local nodal ID's for displacement
    int         *_dindexMap;

    //! Diffusion Vector
    DiffusionVectors _diffusion;

    //! Radius of bar element
    Real         _radius;

    //! Length of Constraint element
    Real         _length;

    //! Set of nodes which are in Purkinje Structure
    std::set<int> _PurkinjeNodes;

    /*! 
      For cell we have 9 choices. Vector holds this choice. 
            Apex Center Base
      =======================
      Epi    0     1     2
      Myo    3     4     5
      Endo   6     7     8
    */
    std::vector<int> _TissueType;

    //! Essential Boundary Conditions
    BC            _essentialBC;
    
    //! Natural Boundary Conditions
    BC            _naturalBC;

    //! Degrees of freedom for the Model
    int           _nDof;

    //! List of Total Deformation Gradient
    DefGradient  _F;

    //! Reference Configuration Check
    bool         _isReferenceConfiguration;

    //! list of PMJ Elements. Holds elemental ID
    IDList       _pmjList;
  };

}

#endif
