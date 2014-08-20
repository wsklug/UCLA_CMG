#include "Mesh.h"
#include "HexMesh.h"
#include "TetMesh.h"
#include "QuadMesh.h"
#include "Tri3Mesh.h"
#include "MeshFree2D.h"
#include "MeshFree3D.h"
#include "MixedMesh3D.h"


void OpenFile(std::ifstream& inp, const std::string filename) {
  inp.open(filename.c_str());
  if (!inp.is_open()) {
    std::cerr << "Error opening file " << filename << std::endl;
    exit(0);
  }
}

namespace voom {
  //! The constructor for the base class
  Mesh::Mesh(const Mesh::Position& nodes, const Mesh::Position& ghNodes,
	     const DiffusionVectors& diffusion,
	     const Mesh::IDTable& elementConnectivity, const int nLocalnodes, 
	     const int nGhostnodes, int *localNodeID, 
	     int *ghostNodeID, int* indexMap, int *dlocalNodeID, 
	     int* dghostNodeID, int *dindexMap, const int quadOrder, 
	     const int nNodes, std::set<int> &PurkinjeNodes,
	     std::vector<int> &TissueType, BC& essentialBC, BC& naturalBC):
    _quadOrder(quadOrder), _nNodes(nNodes), _PurkinjeNodes(PurkinjeNodes),
    _TissueType(TissueType), _essentialBC(essentialBC), _naturalBC(naturalBC) {
    _nodes        = nodes;
    _ghNodes      = ghNodes;
    _diffusion    = diffusion;
    _connectivity = elementConnectivity;
    _nLocalNodes  = nLocalnodes;
    _nGhostNodes  = nGhostnodes;
    _ghostNodeID  = ghostNodeID;
    _localNodeID  = localNodeID;
    _indexMap     = indexMap;
    _dlocalNodeID = dlocalNodeID;
    _dghostNodeID = dghostNodeID;
    _dindexMap    = dindexMap;
    _isIsotropic = true;
    // Stable Element length for constraint element. Default value
    _length       = 0.01;
    // Default Radius value
    _radius       = 0.01;
    // Set Reference Configuration
    _isReferenceConfiguration = true;
    // See if model is isotropic or not
    if ( _diffusion.size() != 0) {
      if ( _diffusion.size() != _connectivity.size() ) {
	std::cerr << "Diffusion Vector and Connectivity size mismatch \n";
	exit(0);
      }
      _isIsotropic = false;
    }
  }

  // Get QuadPointList. This is the Quadpoint struct at an element/node
  void Mesh::getIndices(const int qpoint, Mesh::IDList& indices) {
    indices[0] = _numQuadPoints*qpoint;
    indices[1] = indices[0] + _numQuadPoints;
  }

  //! Update Nodal Position based on displacement
  void Mesh::updateNodalPosition(Epetra_Vector* Displacement, 
				 Epetra_Vector* ghostDisplacement) {
    Real *disp, *ghostDisp;
    Displacement->ExtractView( &disp );
    ghostDisplacement->ExtractView( &ghostDisp );
    for(int i = 0; i < _nLocalNodes; i++) {
      int id = _indexMap[ _localNodeID[i] ] - 1;
      for(int j = 0; j < _nDof; j++) 
	_nodes[i][j] += disp[_nDof*id + j];
    }
    for(int i = 0; i < _nGhostNodes; i++) {
      int id = -_indexMap[ _ghostNodeID[i] ] -1;
	for(int j = 0; j < _nDof; j++) 
	  _ghNodes[i][j] += ghostDisp[_nDof*id + j];
    }
  }

  //! Compute deformation gradients
  void Mesh::computeDeformationGradient(){
    _F.clear();
    int numElems = _connectivity.size();
    IDList index(2,0);
    for(int i = 0; i < numElems; i++) {
      int nodePerElem = _connectivity[i].size();
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
	if (_dimension == 2) F(2,2) = 1.;
	for(int k = 0; k < nodePerElem; k++)
	  for(int m = 0; m < _dimension; m++) 
	    for(int n = 0; n < _dimension; n++)
	      F(m,n) += nds[k][m]*_initialPoints[j].shapeDerivatives[k][n];
	_F.push_back( F );
	//	std::cout << F << "\n";
      } // Quad point loop
    }// Number of elements loop
  }

  //! Interpolate Scalar to integration points
  std::vector<Real> Mesh::interpolateScalar(std::vector<Real>& gamma,
					    std::map<int, int>& nodeID) {
    std::vector<Real> ipValues;
    int numElems = _connectivity.size();
    IDList index(2,0);
    for(int i = 0; i < numElems; i++) {
      int nodePerElem = _connectivity[i].size();
      Real nds[nodePerElem];
      for(int j = 0; j < nodePerElem; j++) {
	int id = nodeID[ _connectivity[i][j] - 1 ];
	nds[j] = gamma[ id ];
      } // j loop

      this->getIndices( i, index );
      // Looping over the Quad Points
      for(int j = index[0]; j < index[1]; j++) {
	Real value = 0.;
	for(int k = 0; k < nodePerElem; k++)
	  value += nds[k]*_points[j].shapeFunctions[k];
	ipValues.push_back( value );
      } // Quad point loop
    }// Number of elements loop 
    return ipValues;
  }

  //! Sets the mesh to the desired meshtype and returns a const pointer to it
  Mesh* Mesh::New(Epetra_MpiComm *mpicomm, std::string modelName, 
		  const int quadOrder, MeshType type) {
    Mesh* myMesh;
    char nProcessors[10];
    Mesh::Position nodes, ghNodes;
    Mesh::IDTable elementConnectivity;
    int nLocalNodes, nGhostNodes, *localNodeID, *ghostNodeID, *indexMap;
    int *dlocalNodeID, *dghostNodeID, *dindexMap, nDof;
    std::ifstream inp;
    int nNodes, dim, nElems, temp, procNumber, id, dataPerLine;
    std::map<MeshType,int> numDiff;
    DiffusionVectors diffusion;
    std::vector<Real> metric, support, support_Hat; // Meshfree relevant data
    std::vector<int> nodalID; // Meshfree relevant data
    std::set<int> PurkinjeNodes;
    struct stat fileInfo;
    BC essentialBC, naturalBC;
    std::map<int, std::vector<Real> > orientation; // Angular Orientation

    numDiff[HEXMESH]  = numDiff[TETMESH]  = numDiff[MESHFREE3D] = 9;
    numDiff[QUADMESH] = numDiff[TRI3MESH] = numDiff[MESHFREE2D] = 4;
    numDiff[MIXEDMESH3D] = 9;

    sprintf(nProcessors, "%d", mpicomm->NumProc());
    std::string nodefile   = modelName + ".node";
    std::string elemfile   = modelName + ".ele";
    std::string fiberfile  = modelName + ".fibers";
    std::string npartsfile = modelName + ".metis.npart." + nProcessors;
    std::string epartsfile = modelName + ".metis.epart." + nProcessors;
    std::string mffile     = modelName + ".mf";
    std::string purkfile   = modelName + ".purk";
    std::string restraints = modelName + ".rest";
    std::string forces     = modelName + ".force";
    std::string orient     = modelName + ".angle";
    std::string cells      = modelName + ".cell"; // Epi Endo Myo Ap base Cen
    std::vector<int> TissueType;
    /*!
      For tissue cell we have 9 choices. 
            Apex Center Base
      =======================
      Epi    0     1     2
      Myo    3     4     5
      Endo   6     7     8

      Be Default all tissue cell type is assigned to 0

      For Purkinje nodes we use the following choice
               LV        RV
      =======================
      Purkinje  11        12
     
      By default all Purkinje nodes assigned to 11
    !*/

    if ( mpicomm->MyPID() == 0 ) {
      std::cout << "Modelname: " << modelName << std::endl;
      std::cout << "Number of Processors: " << nProcessors << std::endl;
    }
    OpenFile(inp, nodefile);

    inp >> nNodes >> dim >> temp >> temp;
    nLocalNodes = 0;
    inp.close();

    // See if there is a file of restraints
    if(stat(restraints.c_str(), &fileInfo) == 0) {    
      int N;
      OpenFile(inp, restraints);
      inp >> N;
      essentialBC.nodalID.resize(N);
      essentialBC.nodalDOF.resize(N);
      essentialBC.values.resize(N);
      for(int i = 0; i < N; i++) 
	inp >> essentialBC.nodalID[i] >> essentialBC.nodalDOF[i]
	    >> essentialBC.values[i];
      inp.close();
    }
    // See if there is a file of forces
    if(stat(forces.c_str(), &fileInfo) == 0) {    
      int N;
      OpenFile(inp, forces);
      inp >> N;
      naturalBC.nodalID.resize(N);
      naturalBC.nodalDOF.resize(N);
      naturalBC.values.resize(N);
      for(int i = 0; i < N; i++) 
	inp >> naturalBC.nodalID[i] >> naturalBC.nodalDOF[i]
	    >> naturalBC.values[i];
      inp.close();
    }

    // Find number of nodes per processor
    OpenFile(inp, npartsfile);
    for(int i = 0; i < nNodes; i++) {
      inp >> procNumber;
      if ( mpicomm->MyPID() == procNumber ) nLocalNodes++;
    }
    inp.close();
    
    // List which points to the node number for a processor
    localNodeID = new int[nLocalNodes];
    if (numDiff[type] == 9) nDof = 3; else nDof = 2;
    dlocalNodeID = new int[nDof*nLocalNodes];

    OpenFile(inp, npartsfile);
    id = 0;
    for(int i = 0; i < nNodes; i++){
      inp >> procNumber;
      if( mpicomm->MyPID() == procNumber ) { 
	localNodeID[id] = i;
	for(int p = 0; p < nDof; p++) 
	  dlocalNodeID[nDof*id +p] = nDof*i + p;
	id++;
      }
    }
    inp.close();
    
    // Used to create epetra map
    indexMap = new int[nNodes];
    dindexMap = new int[nDof*nNodes];
    for(int i = 0; i < nNodes; i++) {
      indexMap[i] = 0;
      for(int p = 0; p < nDof; p++)
	dindexMap[nDof*i + p] = 0;
    }

    for(int i = 0; i < nLocalNodes; i++) {
      indexMap[ localNodeID[i] ] = i+1;
      for(int p = 0; p < nDof; p++) 
	dindexMap[ dlocalNodeID[nDof*i + p] ] = nDof*i + p + 1;
    }

    // Read Elemental mapping information
    OpenFile(inp, elemfile);
    std::ifstream eparts;
    OpenFile(eparts, epartsfile);

    inp >> nElems >> temp >> temp;
    nGhostNodes = 0;
    for(int i = 0; i < nElems; i++){
      eparts >> procNumber;;
      inp >> id >> dataPerLine;
      for(int j = 0; j < dataPerLine; j++){
	inp >> id;
	if (mpicomm->MyPID() == procNumber) {
	  if ( indexMap[id - 1] == 0 ) {
	    indexMap[id - 1] = -2;
	    nGhostNodes++;
	  }
	}
      } // j loop
    } // i loop

    inp.close();
    eparts.close();

    ghostNodeID = new int[nGhostNodes];
    dghostNodeID = new int[nDof*nGhostNodes];
    id = 0;
    for(int i = 0; i < nNodes; i++) {
      if ( indexMap[i] == -2 ) {
	ghostNodeID[id] = i;
	indexMap[i] = -(id + 1);
	for(int p = 0; p < nDof; p++) 
	  dghostNodeID[ nDof*id + p ] = nDof*i + p;
	id++;
      }
    }

    for(int i = 0; i < nGhostNodes; i++) {
      for(int p = 0; p < nDof; p++) 
	dindexMap[ dghostNodeID[nDof*i + p] ] = -(nDof*i + p + 1);
    }

    // Read position info for local and ghost nodes
    OpenFile(inp, nodefile);
    inp >> nNodes >> dim >> temp >> temp;
    for(int i = 0; i < nNodes; i++) {
      std::vector<Real> points(dim, 0.);
      inp >> id;
      for(int j = 0; j < dim; j++) inp >> points[j];
      if (indexMap[i] > 0 ) nodes.push_back(points);
      else if (indexMap[i] < 0) ghNodes.push_back(points);
    }
    inp.close();

    // ID Table parsing
    OpenFile(inp, elemfile);
    OpenFile(eparts, epartsfile);
    std::ifstream fiber(fiberfile.c_str()); // Fiber data file
    inp >> nElems >> temp >> temp;
    for(int i = 0; i < nElems; i++) {
      int eId;
      inp >> eId >> dataPerLine;
      eparts >> procNumber;
      std::vector<int> conn(dataPerLine);
      for(int j = 0; j < dataPerLine; j++) inp >> conn[j];
      // Read fiber data
      if ( fiber.is_open() ) {
	std::vector< Real >  datas(numDiff[type]);
	for(int j = 0 ; j < numDiff[type]; j++) fiber >> datas[j];
	if ( mpicomm->MyPID() == procNumber) diffusion.push_back(datas);
      }
      if ( mpicomm->MyPID() == procNumber) {
	elementConnectivity.push_back(conn);
	nodalID.push_back(eId);
      }
    }
    inp.close();
    eparts.close();
    if ( fiber.is_open() ) fiber.close();

    // If Meshfree input read the MF File
    if ( type == MESHFREE2D || type == MESHFREE3D ){
      Real metricValue, supportValue, supportHat;
      OpenFile(inp, mffile);

      metric.resize( nLocalNodes + nGhostNodes );
      support.resize( nLocalNodes + nGhostNodes );
      support_Hat.resize( nLocalNodes + nGhostNodes );
      for(int i = 0; i < nNodes; i++) {
	inp >> metricValue >> supportValue >> supportHat;
	int id = indexMap[i];
	if (id == 0) continue;
	if (id > 0 ) id = id - 1;
	else id = nLocalNodes - id - 1;
	support[id ] = supportValue;
	metric[id] = metricValue;
	support_Hat[ id ] = supportHat;
      }
      inp.close();
    }

    TissueType.resize( nNodes, 0);
    // Read nodes which define Purkinje file if present
    if(stat(purkfile.c_str(), &fileInfo) == 0) {
      int nPurkNodes, nodeID;
      OpenFile(inp, purkfile);
      inp >> nPurkNodes;
      for(int i = 0; i < nPurkNodes; i++) {
	inp >> nodeID;
	PurkinjeNodes.insert( nodeID );
	// Set PurkinjeNode Type to 11 (Default is LV)
	TissueType[ nodeID - 1 ] = 11;
      }
      inp.close();
    }

    // See if cell file is present
    if (stat(cells.c_str(), &fileInfo) == 0) {
      if ( mpicomm->MyPID() == 0 ) 
	std::cout << "** TISSUETYPE DEFINITION EXISTS\n";
      std::ifstream cell; int ID, cellID;
      OpenFile( cell, cells);
      while ( !cell.eof() ) {
	cell >> ID >> cellID;
	TissueType[ ID - 1] = cellID;
      }
      cell.close();
    }
    
    // Read Angular Orientation file
    if(stat(orient.c_str(), &fileInfo) == 0) {
      if ( mpicomm->MyPID() == 0 ) 
	std::cout << "** PMJ INTRINSIC  DEFINITION EXISTS\n";
      std::ifstream ang; int ID; 
      OpenFile( ang, orient);
      while ( !ang.eof() ) {
	std::vector<Real> data(4, 0.);
	ang >> ID >> data[0] >> data[1] >> data[2] >> data[3];
	if (ang.eof() ) break;
	orientation[ID-1] = data;
      }
      ang.close();
    }
    
    switch(type)
      {
      case 0:
        myMesh = new HexMesh(nodes, ghNodes, diffusion, elementConnectivity, 
			     nLocalNodes, nGhostNodes, localNodeID, 
			     ghostNodeID, indexMap, dlocalNodeID, dghostNodeID,
			     dindexMap, quadOrder, nNodes, PurkinjeNodes,
			     TissueType, essentialBC, naturalBC);
        break;
      case 1:
        myMesh = new TetMesh(nodes, ghNodes, diffusion, elementConnectivity, 
			     nLocalNodes, nGhostNodes, localNodeID, 
			     ghostNodeID, indexMap, dlocalNodeID, dghostNodeID,
			     dindexMap, quadOrder, nNodes, PurkinjeNodes,
			     TissueType, essentialBC, naturalBC);
        break;
      case 2:
        myMesh = new QuadMesh(nodes, ghNodes, diffusion, elementConnectivity, 
			      nLocalNodes, nGhostNodes, localNodeID, 
			      ghostNodeID, indexMap, dlocalNodeID,dghostNodeID,
			      dindexMap, quadOrder, nNodes, PurkinjeNodes,
			      TissueType, essentialBC, naturalBC);
	break;
      case 3:
        myMesh = new Tri3Mesh(nodes, ghNodes, diffusion, elementConnectivity, 
			      nLocalNodes, nGhostNodes, localNodeID, 
			      ghostNodeID, indexMap, dlocalNodeID,dghostNodeID,
			      dindexMap, quadOrder, nNodes, PurkinjeNodes,
			      TissueType, essentialBC, naturalBC);
	break;
      case 4:
	myMesh = new MeshFree2D(nodes, ghNodes, diffusion, elementConnectivity,
				nLocalNodes, nGhostNodes, localNodeID, 
				ghostNodeID,indexMap,dlocalNodeID,dghostNodeID,
				dindexMap, quadOrder, nNodes, 
				nodalID, metric, support, support_Hat,
				PurkinjeNodes, TissueType, essentialBC, 
				naturalBC);
	break;
      case 5:
	myMesh = new MeshFree3D(nodes, ghNodes, diffusion, elementConnectivity,
				nLocalNodes, nGhostNodes, localNodeID, 
				ghostNodeID,indexMap,dlocalNodeID,dghostNodeID,
				dindexMap, quadOrder, nNodes,
				nodalID, metric, support, support_Hat,
				PurkinjeNodes, TissueType, essentialBC, 
				naturalBC);
	break;
      case 6:
	myMesh = new MixedMesh3D(nodes, ghNodes, diffusion,elementConnectivity,
				 nLocalNodes, nGhostNodes, localNodeID, 
				 ghostNodeID, indexMap, dlocalNodeID,
				 dghostNodeID, dindexMap, quadOrder, nNodes, 
				 PurkinjeNodes, TissueType, essentialBC, 
				 naturalBC, orientation);
	break;
      default:
        {
          std::cout << "ERROR: Unknown Mesh Type\n";
          exit(0);
	}
      }
    return myMesh;
  }
}
