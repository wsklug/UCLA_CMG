#include "CardiacBody.h"
// Set up parameters for Ionic Models
#include "SetUpIonicConstants.h"
// #include "Stewart.h"

namespace voom {
  // Constructor
  CardiacBody::CardiacBody(Mesh* const mesh, const IonicModels IonicModel){
    _myMesh                = mesh;
    _ionicModel            = IonicModel;
    _useDiffusiveSmoothing = false;
    _alpha                 = 0.1;
    const int numNodes = _myMesh->getNumLocalNodes() + 
      _myMesh->getNumGhostNodes();
    // stimulus current size = numlocal + numghost
    _stimulusCurrent.resize( numNodes );
    _myMesh->Compute(); // No need to call Compute later
    // UserSpecifiedXi
    _userXiSet = false;
    _C_m = 1.; // Nodal capacitance 1.0 uF/cm^2
    for(int i = 0; i < numNodes; i++) _stimulusCurrent[i] = 0.;
    // Initialize diffusion to default values 
    _D[0] = 0.001;    _D[1] = 0.001;    _D[2] = 0.001;
    // Initialize Purkinje diffusion values
    _Dp[0] = 0.0016;    _Dp[1] = 0.0016;    _Dp[2] = 0.0016;
    // Initialize LV Purkinje diffusion values
    _DpLeft[0] = 0.0016;    _DpLeft[1] = 0.0016;    _DpLeft[2] = 0.0016;
  }


  // Lumped Capacitance computation
  void CardiacBody::computeLumpedCapacitance(Epetra_Vector *EVcapacitance,
					     Epetra_Vector *EVghostCap) {
    // Elements in FE or Nodes in Mesh Free
    int numElems = _myMesh->getIDTableSize();
    const int *indexMap = _myMesh->getIndexMap();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    const int nGhostNodes = _myMesh->getNumGhostNodes();
    Real *capacitance, *ghostCapacitance, Xi;
    EVcapacitance->ExtractView(&capacitance);
    EVghostCap->ExtractView(&ghostCapacitance);

    // Shape Function and Deriatives
    const Mesh::QuadPoints& QuadVals = _myMesh->getPoints();
    for(int i = 0; i < nLocalNodes; i++) capacitance[i] = 0.;
    for(int i = 0; i < nGhostNodes; i++) ghostCapacitance[i] = 0.;

    // Idlist
    Mesh::IDList connectivity, indices(2);

    for(int i = 0; i < _myMesh->getIDTableSize(); i++ ) {
      // Get indices in the quadstruct which hold the necessary values
      _myMesh->getIndices(i, indices);
      connectivity = _myMesh->getIDList(i);
      // Loop over the number of Quad Points for the element/node
      for(int j = indices[0]; j < indices[1]; j++) {
	for(int k = 0; k < connectivity.size(); k++) {
	  Real Value = 
	    QuadVals[j].shapeFunctions[k] * _C_m * QuadVals[j].weight;
	  // Get correct ID position in capacitance array
	  int id = indexMap[ connectivity[k] - 1];
	  if ( id > 0 ) capacitance[id - 1] += Value;
	  else ghostCapacitance[ -id - 1]   += Value;
	}// k loop
      } // j loop      
    } // i loop
  }

  /*! \brief
    Conductivity Tensor Calculation for each element. This computation
    needs to be performed once and this will be used by the explicit and 
    implicit solver as needed. Also The full capacitance matrix is 
    computer
    \f[
    C_{ab}=\int_{\Omega} c_m N_a N_b\,dV
    \f]
  */
  void CardiacBody::Compute() {
    // Elements in FE or Nodes in Mesh Free
    const int numElems = _myMesh->getIDTableSize();
    const int dim_t = _myMesh->getDimension(); 
    const std::vector<  std::vector< Real > >& 
      diffVec = _myMesh->getDiffusionVectors();
    const bool isIsotropic = _myMesh->isIsotropic();
    const std::vector<int>& tissueType = _myMesh->getTissueType();
    std::vector< Vector3D > Identity;
    Identity.push_back( Vector3D(1.,0.,0.) );
    Identity.push_back( Vector3D(0.,1.,0.) );
    Identity.push_back( Vector3D(0.,0.,1.) );

    // Shape Function and Deriatives
    const Mesh::QuadPoints& QuadVals = _myMesh->getPoints();
    // Idlist
    Mesh::IDList connectivity, indices(2);

    // Set diffusion based on if model is Isotropic or Anisotropic
    Real D[3];

    for(int i = 0; i < numElems; i++ ) {
      // Get indices in the quadstruct which hold the necessary values
      _myMesh->getIndices(i, indices);
      // Stiffness Matrix
      ElementMatrix diffusion, Cmat;
      connectivity = _myMesh->getIDList(i);
      const int numCon = connectivity.size();
      /* 
	 Decide whether Purkinje is in LV or RV. If both nodes are 11
	 or atleast one is 11 we are in LV. Else RV
      */
      bool isInLV = false;
      if ( numCon == 2) {
	if ( (tissueType[connectivity[0]-1] == 11) || 
	     (tissueType[connectivity[1]-1] == 11) ) isInLV = true;
      }
      if ( isIsotropic) {
	// For purkinje(2 noded elements) use Dp values
	if (numCon == 2 ) 
	  D[0] = D[1] = D[2] = _Dp[0];
	// All other elements
	else D[0] = D[1] = D[2] = _D[0];
      }
      else {
	// 2 noded beam elements
	if (numCon == 2) {
	  if ( isInLV) {
	    D[0] = _DpLeft[0]; D[1] = _DpLeft[1]; D[2] = _DpLeft[2]; 
	  } else {
	    D[0] = _Dp[0]; D[1] = _Dp[1]; D[2] = _Dp[2]; 
	  }
	}
	// all other elements
	else {
	  D[0] = _D[0];	D[1] = _D[1]; D[2] = _D[2]; 
	}
      }

      diffusion.resize( numCon, numCon);
      Cmat.resize(numCon, numCon);
      diffusion = 0.; Cmat = 0.;

      // Loop over the number of Quad Points for the element/node
      for(int j = indices[0]; j < indices[1]; j++) {
	// Diffusion Matrix
	Real conductivity_matrix[dim_t][dim_t];
	std::vector< Vector3D > Ibar;
	for(int m = 0; m < 3; m++)
	  Ibar.push_back( Vector3D(0.,0.,0.) );
	  
	for( int m = 0; m < dim_t; m++)
	  for( int n = 0; n < dim_t; n++) conductivity_matrix[m][n] = 0.;
	
	// Calculate the conductivity tensor using anisotropic fiber angles
	if (isIsotropic) {
	  Ibar = Identity;
	  for(int p = 0; p < dim_t; p++) conductivity_matrix[p][p] = D[p];
	} 
	else {
	  Real f[dim_t];
	  for(int p = 0; p < dim_t; p++){
	    for(int l = 0; l < dim_t; l++) {
	      f[l] = diffVec[i][p*dim_t + l]; 
	      Ibar[p](l) = f[l];
	    }
	    for(int m = 0; m < dim_t; m++)
	      for(int n = 0; n < dim_t; n++) 
		conductivity_matrix[m][n] += D[p]*f[m]*f[n];
	  } // p loop
	} // else loop


	for(int a = 0; a< numCon; a++)
	  for(int b = 0; b < numCon; b++) {
	    Cmat(a,b) += QuadVals[j].shapeFunctions[a] *
	      QuadVals[j].shapeFunctions[b] * QuadVals[j].weight;
	    for(int m = 0; m < dim_t; m++)
	      for(int n = 0; n < dim_t; n++) {
		diffusion(a,b) += 
		  conductivity_matrix[m][n]*QuadVals[j].shapeDerivatives[a][m]
		  *QuadVals[j].shapeDerivatives[b][n]*QuadVals[j].weight*_C_m;
	      } // n loop
	  } // b loop
	_Ibar.push_back( Ibar );
      } // j loop i.e Quadpoint loop index    
      _diffusion.push_back(diffusion);
      _fullMass.push_back( Cmat );
    } // i loop        
  } // End of Function call

  /*!
    Explicit Function Call
    This computes Ionic Current from Diffusion. This part performs
    \f$ I_i^t = D_{ij} V_j^t \f$. Used Voltage and time "t" to compute the
    diffusive current at \f$t^{n+1}\f$. Current is also a Epetra_Vector
    \deprecated This function call is probably faster for small models. The 
    routine loops over each element performs \f$ D\times V\f$ and sums up 
    the current. Hence for a large model run time will probably be bad. The 
    Explicit and Crank-Nicolson solver Assemble \f$D\f$ and perform a Matrix 
    Vector multiplication which is faster.
  */
  void CardiacBody::getDiffusionCurrent(Epetra_Vector *DiffusiveCurrent,
					Epetra_Vector *GhostDiffusiveCurrent,
					Epetra_Vector *Voltage,
					Epetra_Vector *GhostVoltage){
    int numEntity = _diffusion.size();
    Real *diffusiveCurrent, *ghostDiffusiveCurrent;
    Real *voltage, *ghostVoltage;
    const int *indexMap = _myMesh->getIndexMap();
    Voltage->ExtractView(&voltage);
    GhostVoltage->ExtractView(&ghostVoltage);
    DiffusiveCurrent->ExtractView(&diffusiveCurrent);
    GhostDiffusiveCurrent->ExtractView(&ghostDiffusiveCurrent);

    Mesh::IDList connectivity;
    Real sum;
    for(int i = 0; i < numEntity; i++) {
      connectivity = _myMesh->getIDList(i);
      // Get voltage at relevant nodes for element in consideration	
      int size = connectivity.size();

      Real nodalVoltage[size];
      for(int j = 0; j < size; j++) {
	int id = indexMap[ connectivity[j] - 1 ];
	if ( id == 0 ) continue;
	if ( id > 0 )
	  nodalVoltage[j] = voltage[ id - 1 ];
	else
	  nodalVoltage[j] = ghostVoltage[ -id - 1 ];
      }

      // Now Perform Matrix multiplication of stiffness * nodalVoltage
      for(int a = 0 ; a < size; a++) {
	int id = indexMap[ connectivity[a] - 1 ];
	if (id == 0) continue;			
	sum = 0.;
	for(int b = 0; b < size; b++) sum -= _diffusion[i](a,b)*nodalVoltage[b];
	if ( id > 0 ) diffusiveCurrent[id - 1] += sum;
	else ghostDiffusiveCurrent[-id - 1] +=sum;
      } // a loop
    } // i loop
  } // End of Function

  // Set stimulus Current for a region
  void CardiacBody::setStimulusCurrent(const int Axis, const Real limit, 
				       const Real stim) {
    const int numNodes = _myMesh->numberOfNodes();
    const int *indexMap = _myMesh->getIndexMap();
    const Mesh::Position& nodes = _myMesh->getNodes();
    const Mesh::Position& ghNodes = _myMesh->getGhostNodes();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    Real value;
    for(int i =0; i < numNodes; i++) {
      int id = indexMap[i];
      if ( id > 0 ) 
	value = nodes[id-1][Axis];
      else if (id < 0 ) 
	value = ghNodes[-id - 1][Axis];
      else continue;
      if (value <= limit) {
	id = ( id > 0 ) ? id - 1 : nLocalNodes - id - 1;
	_stimulusCurrent[id] = stim;
      }
    }
  }

  // Set stimulus for nodes below a plane
  void CardiacBody::setStimulusCurrent(Real xo[], Real xn[], const Real stim) {
    const int *indexMap = _myMesh->getIndexMap();
    const Mesh::Position& nodes = _myMesh->getNodes();
    const Mesh::Position& ghNodes = _myMesh->getGhostNodes();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    const std::set<int> &purkinjeNodes = _myMesh->getPurkinjeNodes();
    
    tvmet::Vector<Real, 3> x(0);
    tvmet::Vector<Real, 3> normal(xn[0],xn[1],xn[2]); tvmet::normalize(normal);
    tvmet::Vector<Real, 3> origin(xo[0],xo[1],xo[2]);
    for(unsigned int i = 0; i < _myMesh->numberOfNodes(); i++) {
      // If node is a Purkinje Node then do not stimulate it.
      if ( purkinjeNodes.find( i + 1) != purkinjeNodes.end() ) continue;      
      int id = indexMap[i];
      if (id == 0) continue;
      for(int j = 0; j < _myMesh->getDimension(); j++)
	if ( id > 0 )   x(j) = nodes[id-1][j];
	else  x(j) = ghNodes[-id - 1][j];
      id = ( id > 0 ) ? id - 1: nLocalNodes - id - 1;
      if ( tvmet::dot( x-origin, normal) < 0.) _stimulusCurrent[id] = stim;
      else _stimulusCurrent[id] = 0.;
    }
  }

  // Stimulus for a circular region
  void CardiacBody::setStimulusCurrent(Real xc[], Real rad, const Real stim) {
    const int *indexMap = _myMesh->getIndexMap();
    const Mesh::Position& nodes = _myMesh->getNodes();
    const Mesh::Position& ghNodes = _myMesh->getGhostNodes();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    const std::set<int> &purkinjeNodes = _myMesh->getPurkinjeNodes();
    tvmet::Vector<Real, 3> x(0.), origin(xc[0], xc[1], xc[2]);
    for(int i = 0; i < _myMesh->numberOfNodes(); i++) {
      // If node is a Purkinje Node then do not stimulate it.
      if ( purkinjeNodes.find( i + 1) != purkinjeNodes.end() ) continue;
      int id = indexMap[i];
      if (id == 0) continue;
      for(int j = 0; j < _myMesh->getDimension(); j++)
	if ( id > 0 ) x(j) = nodes[id-1][j];
	else  x(j) = ghNodes[-id - 1][j];
      // Find the correct id value
      id = ( id > 0 ) ? id - 1: nLocalNodes - id - 1;
      if ( tvmet::norm2(origin - x )<= rad) _stimulusCurrent[id] = stim;
      else _stimulusCurrent[id] = 0.;
    }
  }

  // Set stimulus for a nodeset
  void CardiacBody::setStimulusCurrent(std::vector<int>& nodeSet, 
				       const Real istim) {
    const int* indexMap = _myMesh->getIndexMap();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    for(std::vector<int>::iterator it = nodeSet.begin(); it != nodeSet.end();
	++it) {
      int id = indexMap[ *it - 1];
      if (id > 0) _stimulusCurrent[id - 1] = istim;
      else _stimulusCurrent[nLocalNodes - id - 1] = istim;
    }
  }


  // Set stimulus voltage value
  void CardiacBody::setStimulusVoltage(Epetra_Vector *Voltage,
				       Epetra_Vector* GhostVoltage,
				       std::vector<int>& nodeSet,
				       const Real value) {
    Real *voltage;
    Voltage->ExtractView(&voltage);
    const int *indexMap   = _myMesh->getIndexMap();
    for(std::vector<int>::iterator it = nodeSet.begin(); it != nodeSet.end();
	++it) {
      // Get ID based on indexMap and i. Used in stimcurrent
      int id = indexMap[ *it - 1 ];
      if (id > 0 ) voltage[id - 1] = value;
    } // i loop
  }

  /*! 
    Implicit Function Call. The local quadrature point diffusion matrix
    is assembled into the global Diffusion matrix. Compute Routine in Body
    needs to be called before this routine is called.
  */ 
  void CardiacBody::getDiffusionMatrix(Epetra_FECrsMatrix *CondMat, 
				       Real dt) {
    const int numEntity = _diffusion.size();
    Real Value;
    int gRow, gCol; // Global Row and Col position

    for(int i = 0; i < numEntity; i++ ) {
      Mesh::IDList connectivity = _myMesh->getIDList(i); 
      const int numCon = connectivity.size();
      for(int row = 0; row < numCon; row++) {
	for(int col = 0; col < numCon; col++) {
	  Value = _diffusion[i](row, col) * dt;
	  gRow = connectivity[row] - 1; gCol = connectivity[col] - 1;
	  CondMat->InsertGlobalValues(1, &gRow, 1, &gCol, &Value);
	} // col loop
      } // row loop
    }// i loop
  } // End of Routine

  /*!
    Fill in the entries of Capacitance in a matrix. 
    @param C Relevant entries will be placed in this Sparse Matrix yielding 
    a full capacitance matrix
  */
  void CardiacBody::getCapacitanceMatrix(Epetra_FECrsMatrix *C, 
					 const bool isCap){
    const int numEntity = _fullMass.size();
    Real Value;
    const Real fac = (isCap) ? _C_m : 1.;
    int gRow, gCol; // Global Row and Col position

    for(int i = 0; i < numEntity; i++ ) {
      Mesh::IDList connectivity = _myMesh->getIDList(i); 
      const int numCon = connectivity.size();
      for(int row = 0; row < numCon; row++) {
	for(int col = 0; col < numCon; col++) {
	  Value = _fullMass[i](row, col) * fac;
	  gRow = connectivity[row] - 1; gCol = connectivity[col] - 1;
	  C->InsertGlobalValues(1, &gRow, 1, &gCol, &Value);
	} // col loop
      } // row loop
    }// i loop

  }

  /*! \brief 
    Get Ionic current from the nodes. Voltage does not change. Essentially 
    we are doing
    \f[
    I = I_{ion}(\mathbf{u}^{n+1}, V^n)
    \f]
    The state variables get updated with time and not Voltage. This function
    needs to be called by Solver which uses NCI. We also call compute at 
    the ghost node locations. We want the internal variables to updated 
    at the ghost nodes also. This is because we need updated internal 
    variables in Mechanics solve. 
  */
  void CardiacBody::getIonicCurrent(Epetra_Vector *IonicCurrent, 
				    Epetra_Vector *Voltage,
				    Epetra_Vector *GhostVoltage,
				    const Real dt_Ion) {
    Real *ic;
    IonicCurrent->ExtractView(&ic);
    Real *voltage, *ghostVoltage;
    Voltage->ExtractView(&voltage);
    GhostVoltage->ExtractView(&ghostVoltage);

    const int *localNodes = _myMesh->getLocalNodeID();
    const int *ghostNodes = _myMesh->getGhostNodeID();
    const int *indexMap   = _myMesh->getIndexMap();
    const int numNodes = _myMesh->getNumLocalNodes();

    for(int i = 0; i < numNodes; i++) {
      int id = indexMap[ localNodes[i] ] - 1;      
      ic[id] = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,dt_Ion, 
					    voltage[id],
					    _stimulusCurrent[id]);
    } // for loop

    // Update the Ghost nodes also
    for(int i = 0 ; i < _myMesh->getNumGhostNodes(); i++) {
      int id = numNodes - indexMap[ ghostNodes[i]] - 1;
      int lid = - indexMap[ ghostNodes[i] ] - 1;
      _ionicModel[id]->Compute_Ion(_userXi, _userXiSet, _C_m, dt_Ion,
				   ghostVoltage[lid],
				   _stimulusCurrent[id]);
    }
  }

  // Get cell model internal data
  void CardiacBody::getStateVariables( std::vector<int>&numVar,
				       std::vector< std::vector<Real> >& data) {
    const int *indexMap   = _myMesh->getIndexMap();
    const int *localNodes = _myMesh->getLocalNodeID();
    numVar.resize( _myMesh->getNumLocalNodes() );
    data.resize( _myMesh->getNumLocalNodes() );
    int nData = 0, id;
    for(int i = 0; i < _myMesh->getNumLocalNodes(); i++) {
      id = indexMap[ localNodes[i] ] - 1;
      const std::vector<Real>& intData = 
	_ionicModel[id]->getInternalParameters( nData );
      numVar[i] = nData;
      data[i]   = intData;
    }
  }

  // Set cell model internal data
  void CardiacBody::setStateVariables(const std::vector< std::vector<Real> >& 
				      data) {
    const int *indexMap   = _myMesh->getIndexMap();
    const int *localNodes = _myMesh->getLocalNodeID();
    for(int i = 0; i < _myMesh->getNumLocalNodes(); i++) {
      int id = indexMap[ localNodes[i] ] - 1;
      _ionicModel[id]->setInternalParameters( data[i] );
    }
  }  

  // Get Internal Force
  void CardiacBody::getInternalForce(Material *material,
				     Epetra_Vector *Force,
				     Epetra_Vector *GForce,
				     std::vector<Tensor3D>& FaInv) {
    Real *force, *gForce;
    Force->ExtractView( &force );
    GForce->ExtractView( &gForce );
    
    // Get deformation gradient at all Quad points
    const Mesh::DefGradient& F = _myMesh->getDeformationGradient();
    // Elements in FE or Nodes in Mesh Free
    const int numElems = _myMesh->getIDTableSize();
    const int dim_t = _myMesh->getDimension(); 
    const int* indexMap = _myMesh->getIndexMap();
    const std::vector<int>& nodalID = _myMesh->getRestraintNodes();
    const std::vector<int>& nodalDOF = _myMesh->getRestraintDOF();
    const std::vector<Real>& values = _myMesh->getRestraintValues();
    // Convert ID*DOF to a set for easy search
    std::set<int> boundarySet;
    for(int i = 0; i < nodalID.size(); i++) {
      int id = (nodalID[i]-1)*dim_t + (nodalDOF[i]-1);
      boundarySet.insert( id );
    }
    // Shape Function and Deriatives in Reference Configuration
    const Mesh::QuadPoints& QuadVals = _myMesh->getReferenceConfigPoints();
    const int dof = _myMesh->getDimension();
    // Idlist
    Mesh::IDList connectivity, indices(2);

    for(int i = 0; i < numElems; i++ ) {
      // Get indices in the quadstruct which hold the necessary values
      _myMesh->getIndices(i, indices);
      connectivity = _myMesh->getIDList(i);
      const int numCon = connectivity.size();
      // Ignore Beam and constraint elements
      if ( numCon <= 2 ) continue;

      // Loop over the number of Quad Points for the element/node
      for(int j = indices[0]; j < indices[1]; j++) {
	// Set the deformation gradient
	Tensor3D Fe(F[j] * FaInv[j]);
	material->setDeformationGradient( Fe );
	material->updateState(false, true, false, _Ibar[j][0], _Ibar[j][1]);
	Tensor3D P( material->piolaStress() );

	P = P * FaInv[j];
	for(int a = 0; a< numCon; a++) {
	  int conn = connectivity[a] - 1;
	  int id = indexMap[ connectivity[a] - 1 ], lid;
	  if ( id == 0 ) continue;
	  for(int b = 0; b < dim_t; b++) 
	    for(int p = 0; p < dim_t; p++) {
	      if ( id > 0 ) {
		lid = id - 1;
		if ( boundarySet.find(conn*dof+b) == boundarySet.end() ) 
		  force[lid*dof+b] += P(b,p)*QuadVals[j].shapeDerivatives[a][p]
		    *QuadVals[j].weight;
	      } else {
		lid = -id - 1;
		if (boundarySet.find(conn*dof+b) == boundarySet.end() ) 
		  gForce[lid*dof+b] += P(b,p)*QuadVals[j].shapeDerivatives[a][p]
		    *QuadVals[j].weight;
	      }
	    } // p loop
	} // a loop
      } // j loop i.e Quadpoint loop index    
    } // i loop
  }
  
  // Get Stiffness Matrix
  void CardiacBody::getStiffnessMatrix(Epetra_FECrsMatrix *K, 
				       Material* material,
				       std::vector<Tensor3D>& FaInv) {
    // Get deformation gradient at all Quad points
    const Mesh::DefGradient& F = _myMesh->getDeformationGradient();
    // Elements in FE or Nodes in Mesh Free
    const int numElems = _myMesh->getIDTableSize();
    const int dim_t = _myMesh->getDimension(); 
    const std::vector<int>& nodalID = _myMesh->getRestraintNodes();
    const std::vector<int>& nodalDOF = _myMesh->getRestraintDOF();
    const std::vector<Real>& values = _myMesh->getRestraintValues();
    // Convert ID*DOF to a set for easy search
    std::set<int> boundarySet;
    for(int i = 0; i < nodalID.size(); i++) 
      boundarySet.insert( (nodalID[i]-1)*dim_t + (nodalDOF[i]-1) );

    // Shape Function and Deriatives in Reference Configuration
    const Mesh::QuadPoints& QuadVals = _myMesh->getReferenceConfigPoints();

    // Idlist
    Mesh::IDList connectivity, indices(2);

    for(int i = 0; i < numElems; i++ ) {
      // Get indices in the quadstruct which hold the necessary values
      _myMesh->getIndices(i, indices);    
      // Stiffness Matrix
      ElementMatrix stiffness;
      connectivity = _myMesh->getIDList(i);
      const int numCon = connectivity.size();
      stiffness.resize( dim_t*numCon, dim_t*numCon);
      stiffness = 0.; 

      // Loop over the number of Quad Points for the element/node
      for(int j = indices[0]; j < indices[1]; j++) {
	// Set the deformation gradient
	Tensor3D Fe(F[j] * FaInv[j]);
	// Compute Lagrangian Moduli
	Array4D Cm; Cm.resize(3,3,3,3);Cm = 0.;
	if (numCon > 2) {
	  material->setDeformationGradient( Fe );
	  material->updateState(false, false, true, _Ibar[j][0], _Ibar[j][1]);
	  const Array4D& C = material->getLagrangianModulus();
	  // Compute C froom C_e
	  // C[iJkL] = C[kNiR]*FaInv(JR)*FaInv(LN)
	  for(int I = 0; I < 3; I++)
	    for(int J = 0; J <3; J++)
	      for(int K = 0; K < 3; K++ ) 
		for(int L = 0; L < 3; L++ ) {
		  Real sum = 0;
		  for(int N = 0; N < 3; N++)
		    for(int R = 0; R < 3; R++) 
		      sum += C(K,N,I,R)*FaInv[j](J,R)*FaInv[j](L,N);
		  Cm(I,J,K,L) = sum;
		}

	  // k[I,A,K,B] = C(I,J,K,L)*N(A,J)*N(B,L)
	  for(int I = 0; I < dim_t; I++)
	    for(int K = 0; K < dim_t; K++)
	      for(int A = 0; A < numCon; A++)
		for(int B = 0; B < numCon; B++) {
		  Real sum = 0.;
		  for(int J = 0; J < dim_t; J++)
		    for(int L = 0; L <dim_t; L++)
		      sum += Cm(I,J,K,L)*QuadVals[j].shapeDerivatives[A][J]
			*QuadVals[j].shapeDerivatives[B][L]*QuadVals[j].weight;
		  stiffness(A*dim_t + I, B*dim_t + K) += sum;
		}
	}// if statement
      } // j loop i.e Quadpoint loop index    
      // Make stiffness diagonal for 1D and constraint elements
      if ( numCon == 2 ) {
	stiffness = 0.;
	for(int a = 0; a < numCon*dim_t; a++) stiffness(a,a) = 1.;
      }	
      // Enforce elemental symmetry
      for (int a = 0; a < numCon*dim_t; a++)
      	for(int b = a+1; b < numCon*dim_t; b++)
	  stiffness(b,a) = stiffness(a,b);

      // Insert stiffness to Epetra matrix
      for(int a = 0; a< numCon; a++) {
	int Row = connectivity[a] - 1;
	for(int b = 0; b < numCon; b++) {
	  int Col = connectivity[b] - 1;
	  for(int p = 0; p < dim_t; p++)
	    for(int q = 0; q < dim_t; q++) {
	      int gRow = Row*dim_t + p, gCol = Col*dim_t + q;
	      if ( (boundarySet.find( gRow ) == boundarySet.end()) && 
		   (boundarySet.find( gCol ) == boundarySet.end()) ) {
		Real Value = stiffness(a*dim_t + p, b*dim_t + q);
		K->InsertGlobalValues(1, &gRow, 1, &gCol, &Value);
	      } // Check if not a boundary node. If so ignore
	    } // q loop
	} // b loop
      } // a loop      
    } // i loop
  }
  
  //! Calculate internal force
  void CardiacBody::getInternalForce(Material *material, 
				     const Mesh::QuadPoints &QuadVals,
				     const int start, const int end,
				     vector< vector<Real> >& x,
				     const int numCon, const int dim_t,
				     blitz::Array<Real, 1>& fanal) {
    for(int q = start; q < end; q++){
      Tensor3D Fe(0.); 
      if (dim_t == 2) Fe(2,2) = 1.;
      for (int a = 0; a < numCon; a++)
	for(int i = 0; i < dim_t; i++)
	  for(int j = 0 ; j < dim_t; j++)
	    Fe(i,j) += x[a][i]*QuadVals[q].shapeDerivatives[a][j];
      material->setDeformationGradient( Fe );
      material->updateState(false, true, false);
      Tensor3D P( material->piolaStress() );
      for(int a = 0; a < numCon; a++)
	for(int i = 0; i < dim_t; i++)
	  for(int j = 0; j < dim_t; j++)
	    fanal(a*dim_t + i) += P(i,j)*QuadVals[q].shapeDerivatives[a][j]*
	      QuadVals[q].weight;
    }// q loop
  }
  
  //! Consistency checks
  bool CardiacBody::consistencyCheck(Material* material) {
    bool success = false;
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n";
    std::cout << "** CONSISTENCY CHECKS \n";
    const int dim_t = _myMesh->getDimension(); 
    // Shape Function and Deriatives in Reference Configuration
    const Mesh::QuadPoints& QuadVals = _myMesh->getReferenceConfigPoints();
    const int* indexMap = _myMesh->getIndexMap();
    // Idlist
    Mesh::IDList connectivity, indices(2);
    // Get indices in the quadstruct which hold the necessary values
    _myMesh->getIndices( 0, indices );
    // Nodal Position
    const Mesh::Position& nodes = _myMesh->getNodes();
    const Mesh::Position& ghNodes = _myMesh->getGhostNodes();
    const Real eps = 1E-6;
    ElementMatrix Kanal, Knum;
    connectivity = _myMesh->getIDList( 0 );
    const int numCon = connectivity.size();
    blitz::Array<Real, 1> fanal, fnum;
    Kanal.resize( dim_t*numCon, dim_t*numCon); Kanal = 0.; 
    Knum.resize( dim_t*numCon, dim_t*numCon);  Knum = 0.; 
    fanal.resize(dim_t*numCon); fanal = 0.;
    fnum.resize(dim_t*numCon);  fnum = 0.;
    srand( time(NULL) );
    vector< vector<Real> > X, x; // Reference and final position
    for(int i = 0; i < numCon; i++) {
      vector<Real> pos;
      int id = indexMap[ connectivity[i] - 1];
      for(int j = 0; j < dim_t; j++) {
	if ( id > 0) pos.push_back( nodes[ id - 1 ][j]);
	else pos.push_back( ghNodes[-id - 1][j]);
      }
      x.push_back( pos ); X.push_back( pos );
    }
    // Perturb final position
    while(1) {
      for(int i = 0; i < numCon; i++)
	for(int j = 0; j < dim_t; j++)
	  x[i][j] = X[i][j] + 1e-3 * Real( rand() )/Real( RAND_MAX );
      // Get def Grad
      bool success = true;
      for(int q = indices[0]; q < indices[1]; q++){
	Tensor3D Fe(0.); 
	if (dim_t == 2) Fe(2,2) = 1.;
	for (int a = 0; a < numCon; a++)
	  for(int i = 0; i < dim_t; i++)
	    for(int j = 0 ; j < dim_t; j++)
	      Fe(i,j) += x[a][i]*QuadVals[q].shapeDerivatives[a][j];
	if ( determinant( Fe ) < 0. ) success = false;;
      }
      if (success) break;
    }
    // Compute internal force analytically
    this->getInternalForce(material,QuadVals, indices[0], indices[1], x, numCon,
			   dim_t, fanal);

    // Compute Internal force numerically
    for(int a = 0; a < numCon; a++)
      for(int i = 0; i < dim_t; i++) {
	x[a][i] += eps;
	Real Wplus = 0., Wminus = 0.;
	for(int q = indices[0]; q < indices[1]; q++) {
	  Tensor3D Fe(0.);
	  if (dim_t == 2) Fe(2,2) = 1.;
	  for (int p = 0; p < numCon; p++)
	    for(int r = 0; r < dim_t; r++)
	      for(int s = 0 ; s < dim_t; s++)
		Fe(r,s) += x[p][r]*QuadVals[q].shapeDerivatives[p][s];
	  material->setDeformationGradient( Fe );
	  material->updateState(true, false, false);
	  Wplus += material->strainEnergy()*QuadVals[q].weight;
	} // q loop
	x[a][i] -= 2*eps;
	for(int q = indices[0]; q < indices[1]; q++) {
	  Tensor3D Fe(0.);
	  if (dim_t == 2) Fe(2,2) = 1.;
	  for (int p = 0; p < numCon; p++)
	    for(int r = 0; r < dim_t; r++)
	      for(int s = 0 ; s < dim_t; s++)
		Fe(r,s) += x[p][r]*QuadVals[q].shapeDerivatives[p][s];
	  material->setDeformationGradient( Fe );
	  material->updateState(true, false, false);
	  Wminus += material->strainEnergy()*QuadVals[q].weight;
	} // q loop
	x[a][i] += eps;
	fnum(a*dim_t + i) = (Wplus - Wminus)/(2.*eps);
      } // i loop
    Real ferror = 0;
    for(int i = 0; i < numCon*dim_t; i++)
      ferror += ( fanal(i) - fnum(i) )*( fanal(i) - fnum(i) );
    ferror = sqrt( ferror );
    std::cout << "Internal force Error  : " << ferror << std::endl;

    // Get stiffness Analytically
    for(int q = indices[0]; q < indices[1]; q++) {
      Tensor3D Fe(0.); 
      if (dim_t == 2) Fe(2,2) = 1.;
      for (int a = 0; a < numCon; a++)
	for(int i = 0; i < dim_t; i++)
	  for(int j = 0 ; j < dim_t; j++)
	    Fe(i,j) += x[a][i]*QuadVals[q].shapeDerivatives[a][j];
      // Set the deformation gradient
      material->setDeformationGradient( Fe );
      material->updateState(false, false, true);
      
      const Array4D& C = material->getLagrangianModulus();
      for(int a = 0; a< numCon; a++) 
	for(int b = 0; b < numCon; b++) 
	  for(int i = 0; i < dim_t; i++)
	    for(int j = 0; j < dim_t; j++) 
	      for(int k = 0; k < dim_t; k++)
		for(int l = 0; l < dim_t; l++) 
		  Kanal(a*dim_t+i, b*dim_t+k) += C(i,j,k,l)*
		    QuadVals[q].shapeDerivatives[a][j]
		    *QuadVals[q].shapeDerivatives[b][l]*QuadVals[q].weight;
    } // q loop i.e Quadpoint loop index    
    
    // Compute Stiffness Numerically
    for(int a = 0; a < numCon; a++) 
      for(int i = 0; i < dim_t; i++) 
	for(int b = 0; b < numCon; b++)
	  for( int j = 0; j < dim_t; j++) {
	    x[b][j] += eps;
	    blitz::Array<Real, 1> fplus; fplus.resize(numCon*dim_t);
	    fplus = 0.;
	    getInternalForce( material, QuadVals, indices[0], indices[1], x, 
			      numCon, dim_t, fplus);
	    x[b][j] -= 2.*eps;
	    blitz::Array<Real, 1> fminus; fminus.resize(numCon*dim_t);
	    fminus = 0.;
	    getInternalForce( material, QuadVals, indices[0], indices[1], x, 
			      numCon, dim_t, fminus);
	    x[b][j] += eps;
	    Knum(a*dim_t+i, b*dim_t+j) = (fplus(a*dim_t+i)-fminus(a*dim_t+i))/
	      (2.*eps);
	  } // b loop
    Real Kerror = 0.;
    for(int i = 0; i < numCon*dim_t; i++)
      for(int j = 0; j < numCon*dim_t; j++)
	Kerror += pow( (Kanal(i,j) - Knum(i,j)), 2);
    Kerror = sqrt(Kerror);

    std::cout << "Stiffness Matrix Error: " << Kerror << std::endl;
    std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"; 
    const Real toler = 1000.*eps;
    if ( (ferror < eps) && (Kerror < toler) ) success = true;
    return success;
  }

  /*! \brief 
    Get Ionic current from the nodes. Then solve ODE to update Voltage
    \f$
    \xi C\frac{\partial V}{\partial t} = \xi I_{ion} + I_{stim} \\
    \approx V^{n+1} = V^n + \Delta t*\frac {I}{C}
    \f$
    The ionic model should return current as uA/cc. This is should be
    dealt in the Ionic model class. The body does NOT do this.
  */
  void CardiacBody::getVoltageFromIonicCurrent(Epetra_Vector *Voltage, 
					       Epetra_Vector *GhostVoltage,
					       const Real dT, 
					       const Real dt_Ion,
					       const int numNodes) {
    Real *voltage, *ghostVoltage;
    Voltage->ExtractView(&voltage);
    GhostVoltage->ExtractView(&ghostVoltage);
    Real iion, dt;
    int id, lid;
    const int *localNodes = _myMesh->getLocalNodeID();
    const int *ghostNodes = _myMesh->getGhostNodeID();
    const int *indexMap   = _myMesh->getIndexMap();
    const int nLocalNodes = _myMesh->getNumLocalNodes();

    for(int i = 0; i < numNodes; i++) {
      // Get ID based on indexMap and i. Used in stimcurrent
      if ( i < nLocalNodes) 
	lid = id = indexMap[ localNodes[i] ] - 1;
      else { 
	id = nLocalNodes - indexMap[ ghostNodes[i-nLocalNodes]] - 1;
	lid = - indexMap[ ghostNodes[i-nLocalNodes] ] - 1;
      }

      // We need to time stepped based on dT and dt_ion
      dt = dt_Ion; // dt gets changed inside while loop. Reset it
      if( dT/dt > 1. ){
	Real totT=0, delT=0;
	int slope;
	if ( i < nLocalNodes) {
	  iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,dt,
					      voltage[lid],
					      _stimulusCurrent[id]);
	  voltage[lid] += dt * iion;
	}
	else {
	  iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,dt,
					      ghostVoltage[lid], 
					      _stimulusCurrent[id]);
	  ghostVoltage[lid] += dt * iion;
	}
	totT = dt;
	slope = iion;
	// for negative slope we can increase stride by a factor of 5
	if(!slope) dt *=5;
	bool flag = true;
	delT = dt;
	while(flag){
	  if( ( totT + dt ) > dT){
	    delT = dT - totT;
	    flag = false;
	  }else{
	    delT = dt;
	  }
	  
	  if ( i < nLocalNodes) {
	    iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,delT, 
						voltage[lid],
						_stimulusCurrent[id]);
	    voltage[lid] += delT * iion;
	  }
	  else {
	    iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,delT, 
						ghostVoltage[lid],
						_stimulusCurrent[id]);
	    ghostVoltage[lid] += delT * iion;
	  }	      
	  totT += delT;
	} // while loop
      }else{
	// just one stride using the critical timestep
	if ( i < nLocalNodes) {
	  iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,dt, 
					      voltage[lid],
					      _stimulusCurrent[id]);
	  voltage[lid] += dt * iion;
	}
	else {
	  iion = _ionicModel[id]->Compute_Ion(_userXi,_userXiSet,_C_m,dt, 
					      ghostVoltage[lid],
					      _stimulusCurrent[id]);
	  ghostVoltage[lid] += dt * iion;
	}
      } // dT/dt loop
    } // for loop
  }


  //! Get active Deformation
  void CardiacBody::getActiveContraction(std::vector<Real>& Gamma,
					 std::map<int, int>& nodeID,
					 const int numNodes) {
    const int *localNodes = _myMesh->getLocalNodeID();
    const int *ghostNodes = _myMesh->getGhostNodeID();
    const int *indexMap   = _myMesh->getIndexMap();
    const int nLocalNodes = _myMesh->getNumLocalNodes();
    int id, lid;

    for(int i = 0; i < nLocalNodes; i++) {
      int id = indexMap[ localNodes[i] ] - 1;
      Gamma.push_back( _ionicModel[id]->getGamma() );
      nodeID[ localNodes[i] ] = i;
    }

    for(int i = 0 ; i < numNodes - nLocalNodes; i++) {
      int id = nLocalNodes - indexMap[ ghostNodes[i] ] - 1;
      Gamma.push_back( _ionicModel[id]->getGamma() );
      nodeID[ ghostNodes[i-nLocalNodes] ] = i+nLocalNodes;
    }
  }

  /*! 
    Compute active deformation matrix and its inverse
    \f[
    F_a = \gamma \mathbf{I}
    \f]
    Return \f[ F_a^{-1} \f]
  */
  std::vector<Tensor3D> 
  CardiacBody::computeActiveDeformationInverse(std::vector<Real>& gamma) {
    std::vector<Tensor3D> FaInvValues;
    const int numElems = _myMesh->getIDTableSize();
    // Shape Function and Deriatives in Reference Configuration
    const Mesh::QuadPoints& QuadVals = _myMesh->getReferenceConfigPoints();
    Mesh::IDList indices(2);
    const int dim_t = _myMesh->getDimension(); 
    for(int i = 0; i < numElems; i++ ) {
      // Get indices in the quadstruct which hold the necessary values
      _myMesh->getIndices(i, indices);
      for(int j = indices[0]; j < indices[1]; j++) {
	Tensor3D Fa(0.), FaInv(0.);
	// I bar stores e1 in first row e2 in second row and e3 in third row
	/*!
	  We compute
	  \f[
	  F_a = \gamma e_1 \otimes e_1 + \frac{1}{\sqrt{\gamma}} (
	  e_2 \otimes e_2 + e_3 \otimes e_3 )
	  \f]
	*/
	const Real fac0 = gamma[j], fac1 = 1. - 0.35*log(gamma[j]);
	const Real fac2 = 1./(fac0*fac1);
	for(int k = 0; k < dim_t; k++) 
	  for(int l = 0; l < dim_t; l++) 
	    // Fa = gamma e1xe1 + gamma e2xe2 + 1/gamma^2 e3xe3
	    Fa(k,l) = fac2 *_Ibar[j][2](k)*_Ibar[j][2](l) +
	      fac1* _Ibar[j][1](k)*_Ibar[j][1](l) +
	      fac0* _Ibar[j][0](k)*_Ibar[j][0](l) ;
	if ( dim_t == 2 ) Fa(2,2) = 1.;
	invert( Fa, FaInv);
	FaInvValues.push_back( FaInv );
      }   
    }
    return FaInvValues;
  }

  //! Creates Cardiac body with the needed CellModelType
  CardiacBody* CardiacBody::New(Mesh* const mesh, CellModelType Type) {
    CardiacBody *myBody;
    const int nLocalNodes = mesh->getNumLocalNodes();
    const int nGhostNodes = mesh->getNumGhostNodes();
    const int numNodes = mesh->getNumLocalNodes() + mesh->getNumGhostNodes(); 
    const int* localNodeID = mesh->getLocalNodeID();
    const int *ghostNodeID = mesh->getGhostNodeID();
    const std::set<int> &purkinjeNodes = mesh->getPurkinjeNodes();
    const std::vector<int>& tissueType = mesh->getTissueType();
    IonicModels Models;
    Real *Constants;
    SetUpPurkinjeParameters(&Constants);

    switch(Type){
    case(0):
      {
	// LuoRudy model without Gate Model
	LuoRudyGateTable *gate = NULL;
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1 ) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new LuoRudy(false, gate);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;
      }
    case(1):
      {
	// LuoRudy model with Gate Model
	LuoRudyGateTable *gate = new LuoRudyGateTable();
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new LuoRudy(true, gate);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;
      }
    case(2):
      {
	// Tusscher Cell model
	Real *TConstants;
	SetUpTusscherParameters(&TConstants);
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new Tusscher(TConstants);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;	
      }
    case(3):
    {
      // UCLA Cell model
      Real *MConstants[9];
      for(int i = 0; i < 9; i++)
	SetUpMahajanParameters( &MConstants[i], i );

      
      for(int i = 0; i < numNodes; i++) {
	int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	  ghostNodeID[i - nLocalNodes];
	if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	  /*
	    int typeID = 0; 
	    IonicModel *ptr = new Mahajan( MConstants[typeID], typeID );
	    Models.push_back(ptr);
	  */
	  IonicModel *ptr = new Purkinje(Constants);
	  Models.push_back(ptr);
	} else {
	  int typeID = tissueType[ id ]; 
	  IonicModel *ptr = new Mahajan( MConstants[typeID], typeID );
	  Models.push_back(ptr);
	}

      } // i loop
      myBody = new CardiacBody(mesh, Models);
      break;
    }
    case(4):
    {
      // UCLA HF Cell model
      Real *MConstants[9];
      for(int i = 0; i < 9; i++)
	SetUpMahajan_failParameters( &MConstants[i], i );
      for(int i = 0; i < numNodes; i++) {
	int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	  ghostNodeID[i - nLocalNodes];
	if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	  IonicModel *ptr = new Purkinje(Constants);
	  Models.push_back(ptr);
	} else {
	  int typeID = tissueType[ id ]; 
	  IonicModel *ptr = new Mahajan_fail( MConstants[typeID], typeID );
	  Models.push_back(ptr);
	}
      } // i loop
      myBody = new CardiacBody(mesh, Models);
      break;
    }
    default:
      {
	std::cerr << "ERROR: Unknown Cell Model Type\n";
	exit(0);
      }
    }
    return myBody;
  }

  //! Creates Cardiac body with the needed CellModelType
  CardiacBody* CardiacBody::New(Mesh* const mesh, CellModelType Type, const std::vector<Real>& purkinjeStates, const std::vector<Real>& myocardiumStates) {
        CardiacBody *myBody;
    const int nLocalNodes = mesh->getNumLocalNodes();
    const int nGhostNodes = mesh->getNumGhostNodes();
    const int numNodes = mesh->getNumLocalNodes() + mesh->getNumGhostNodes(); 
    const int* localNodeID = mesh->getLocalNodeID();
    const int *ghostNodeID = mesh->getGhostNodeID();
    const std::set<int> &purkinjeNodes = mesh->getPurkinjeNodes();
    const std::vector<int>& tissueType = mesh->getTissueType();
    IonicModels Models;
    Real *Constants;
    SetUpPurkinjeParameters(&Constants);

    switch(Type){
    case(0):
      {
	// LuoRudy model without Gate Model
	LuoRudyGateTable *gate = NULL;
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1 ) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    ptr->setInternalParameters(purkinjeStates);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new LuoRudy(false, gate);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;
      }
    case(1):
      {
	// LuoRudy model with Gate Model
	LuoRudyGateTable *gate = new LuoRudyGateTable();
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    ptr->setInternalParameters(purkinjeStates);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new LuoRudy(true, gate);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;
      }
    case(2):
      {
	// Tusscher Cell model
	Real *TConstants;
	SetUpTusscherParameters(&TConstants);
	for(int i = 0; i < numNodes; i++) {
	  int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	    ghostNodeID[i - nLocalNodes];
	  if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	    IonicModel *ptr = new Purkinje(Constants);
	    ptr->setInternalParameters(purkinjeStates);
	    Models.push_back(ptr);
	  } else {
	    IonicModel *ptr = new Tusscher(TConstants);
	    Models.push_back(ptr);
	  }
	}
	myBody = new CardiacBody(mesh, Models);
	break;	
      }
    case(3):
    {
      // UCLA Cell model
      Real *MConstants[9];
      for(int i = 0; i < 9; i++)
	SetUpMahajanParameters( &MConstants[i], i );

      for(int i = 0; i < numNodes; i++) {
	int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	  ghostNodeID[i - nLocalNodes];
	if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	  /*
	    int typeID = 0; 
	    IonicModel *ptr = new Mahajan( MConstants[typeID], typeID );
	    Models.push_back(ptr);
	  */
	  IonicModel *ptr = new Purkinje(Constants);
	  ptr->setInternalParameters(purkinjeStates);
          Models.push_back(ptr);
	} else {
	  int typeID = tissueType[ id ]; 
	  IonicModel *ptr = new Mahajan( MConstants[typeID], typeID );

	  vector<Real> myocardiumStatesTypeID;
	  if (myocardiumStates.size() == 26){
	    // Only default inital values specified
	    myocardiumStatesTypeID = myocardiumStates;
	  }
	  else {
	    for (int iter = typeID * 26; iter < (typeID + 1) * 26; iter++)
	      myocardiumStatesTypeID.push_back(myocardiumStates[iter]);
	  }
	  ptr->setInternalParameters(myocardiumStatesTypeID);
	  
	  Models.push_back(ptr);
	}

      } // i loop
      myBody = new CardiacBody(mesh, Models);
      break;
    }
    case(4):
    {
      // UCLA HF Cell model
      Real *MConstants[9];
      for(int i = 0; i < 9; i++)
	SetUpMahajan_failParameters( &MConstants[i], i );
      
      for(int i = 0; i < numNodes; i++) {
	int id = (i < nLocalNodes) ? localNodeID[ i ] : 
	  ghostNodeID[i - nLocalNodes];
	if ( purkinjeNodes.find( id + 1) != purkinjeNodes.end() ) {
	  /*
	    int typeID = 0; 
	    IonicModel *ptr = new Mahajan( MConstants[typeID], typeID );
	    Models.push_back(ptr);
	  */
	  IonicModel *ptr = new Purkinje(Constants);
	  ptr->setInternalParameters(purkinjeStates);
	  Models.push_back(ptr);
	} else {
	  int typeID = tissueType[ id ]; 
	  IonicModel *ptr = new Mahajan_fail( MConstants[typeID], typeID );

	  vector<Real> myocardiumStatesTypeID;
	  if (myocardiumStates.size() == 26){
	    // Only default inital values specified
	    myocardiumStatesTypeID = myocardiumStates;
	  }
	  else {
	    for (int iter = typeID * 26; iter < (typeID + 1) * 26; iter++)
	      myocardiumStatesTypeID.push_back(myocardiumStates[iter]);
	  }
	  ptr->setInternalParameters(myocardiumStatesTypeID);
	  
	  Models.push_back(ptr);
	}

      } // i loop

      myBody = new CardiacBody(mesh, Models);
      break;
    }
    default:
      {
	std::cerr << "ERROR: Unknown Cell Model Type\n";
	exit(0);
      }
    }
    return myBody;

 }

}// namespace 
