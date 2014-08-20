//-*-C++-*-
// Cardiac Body Class
/*!
  Body class performs the tasks of
  - EP Calculations
      - Computing Conductance matrix
      - Computing Mass matrix for Ionic Current Interpolation

  -  Mechanics Calculations
      - Internal force computation
      - Stiffness Matrix computation

   Pointer to Mesh class is provided as input and Body performs the 
   necessary computations. 
!*/


#if !defined(__CardiacBody_h__)
#define __CardiacBody_h__
#include "voom.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Export.h"
#include "Epetra_Import.h"
#include "Epetra_FECrsMatrix.h"
#include "Epetra_CrsMatrix.h"

#include <vector>

// Mesh Class
#include "Mesh.h"
// Ionic Model class
#include "IonicModel.h"
// Sub Ionic Model classes
#include "LR.h"
// Purkinje Model
#include "Purkinje.h"
// Tusscher Model
#include "Tusscher.h"
// UCLA Cell Model
#include "Mahajan.h"
//UCLA HF model
#include "Mahajan_fail.h"
// Material Class
#include "Material.h"


namespace voom {
  //! Cell Model Type definition
  enum CellModelType {LUORUDY, LUORUDYWITHGATE, TUSSCHER, MAHAJAN, MAHAJAN_FAIL};
  typedef blitz::Array<Real, 2> ElementMatrix;
  typedef std::vector< ElementMatrix > DiffusionStiffness;
  typedef std::vector<IonicModel*> IonicModels;
  /*! \brief CardaicBody:
    A Cardiac Body class. 
    The class takes the output from the Mesh object and performs necessary 
    operation. The solver class provides the body class with pertinent 
    objects and the body class fills it with required data.
    
    The class does not depend on the type of the mesh/solver. This can be
    used for Meshfree/Finite Element/Implicit or Explicit Solver.

    We are solving the PDE
    \f[
    \chi \left( C\frac{dV}{dt} + I_{ion} \right ) = \nabla \cdot (\boldmath{
    \sigma} \nabla V) + I_{stim}
    \f]					       
    where \f$\chi\f$ is the Surface Area to Volume Ratio in cm\f$^{-1}\f$, 
    \f$C\f$ is the membrane capacitance, \f$I_{ion}\f$ is the Ionic current 
    (\f$\mu A/cm^2\f$), \f$\sigma \f$ is the conductivity in \f$mS/cm\f$ and 
    \f$I_{stim}\f$ is the stimulus current in \f$\mu A/cc \f$. 
    \note The \b units are important. Ensure the input is consistent with 
    these units 
  */
  class CardiacBody {
  public:    
    //! Constructor definition
    CardiacBody(Mesh* const mesh, const IonicModels IonicModel);

    //! Destructor
    ~CardiacBody() {
      for(int i = 0; i < _ionicModel.size(); i++)
	delete _ionicModel[i];
    }

    /*! Lumped Capacitance Computation \f$ C_L \f$
      @param capacitance Epetra Vector of local Nodes capacitance
      @param ghostCapacitance Epetra Vector of ghost nodes 
      capacitance
    */
    void computeLumpedCapacitance(Epetra_Vector *capacitance, 
				  Epetra_Vector *ghostCapacitance);


    /*! 
      Computes each element/nodal diffusion matrix. Also each element/node
      Capacitance matrix is computed. Compute call should be done by the 
      solver class
     */
    void Compute();

    /*! 
      Get Diffusion Matrix. Will be called by Implicit Solver
      \f[
      [D_{ij}] = \int_{\Omega} N_{i,a}N_{j,b} D_{ab}\, dV
      \f]
      The values of \f$D_{ij}\f$ are done in compute routine. This routine 
      places the value of \f$D\f$ in the right location
      \param CondMatrix Epetra_FECrsMatrix which values of \f$D_{ij}\f$ 
      are placed.
      \param dt This is a scalar by which each value of \f$D_{ij}\f$ is 
      scaled. This is used by Crank-Nicholson solver.
    */
    void getDiffusionMatrix(Epetra_FECrsMatrix *CondMatrix, Real dt);

    /*! Get Capacitance Matrix.
      \f[
      M_{jk} = k\int_{\Omega} N_j N_k \, dV
      \f]

      @param M Epetra_FECrsMatrix. The entries of \f$M_{ij}\f$ are placed 
      in this matrix. This Matrix needs to be created outside the class
      @param isCap If isCap is true \f$ k = C_m \f$ else \f$ k = 1\f$. User 
      can get Capacitance or Mass Matrix based on this. Default true
    */
    void getCapacitanceMatrix(Epetra_FECrsMatrix *M, const bool isCap = true);

    /*!
      Get internal force
      \f[
      f^{int}_{ia} = \int_V P_{iJ} N_{a,J}\, dV
      \f]
    */
    void getInternalForce(Material *material, Epetra_Vector* force, 
			  Epetra_Vector* gForce, std::vector<Tensor3D>& FaInv);

    /*! Get Stiffness Matrix
      \f[
      K_{iakb}=\int C_{iJkL}N_{a,J}N_{b,L}\, dV
      \f]
      where \f$C\f$ is the Lagrangian Moduli. Material Class return \f$C\f$
      computed only from the elastic part. This needs to be multiplied by 
      the active Deformation Gradient. \f$N\f$ is the shape function

      @param K Epetra_FECrsMatrix. The entires in this matrix are filled by 
      this function.
      @param Material Pointer to material class. This will be used to 
      perform the Lagrangian Modulus Computation
    */
    void getStiffnessMatrix(Epetra_FECrsMatrix *K, Material* material, 
			    std::vector<Tensor3D>& FaInv);

    /*!
      Perform consistency checks. This will perform force and stiffness 
      consistency check for one element. We will take this to be the first 
      element in the model. The internal force is calculated as
      \f[
      f^{int}_{ia} = \int_V P_{i,J} N_{a,J}\, dV
      \f]
      This is the analytical value. Numerically we can compute this as
      \f[
      f_{num}_{ia} = \dfrac{\partial \Pi}{\partial x_{ia}}
      \f]
      The difference between these two should be within some tolerance. The 
      stiffness matrix similarly is
      \f[
      K_{iakb} = \int_V C_{iJkL}N_{a,J}N_{b,L}\, dV
      \f]
      and the numerical value of stiffness is
      \f[
      K^{num}_{iakb} = \dfrac{\partial f_{ia}}{\partial x_{kb}}
      \f]
      Again the difference between these should be within the tolerance. Note
      we are not using the active deformation gradient in this function.
    */
    bool consistencyCheck
( Material* material);

    /*!
      Get Diffusion Current Vector. Will be called by Explicit Solver
      \deprecated The solver class performs a matrix Vector multiplication 
      to get this value. Using this function call on a big model increases 
      computation time.
    */
    void getDiffusionCurrent(Epetra_Vector *DiffusiveCurrent,
			     Epetra_Vector *GhostDiffusiveCurrent,
			     Epetra_Vector *Voltage,
			     Epetra_Vector *GhostVoltage);

    //! Get Ionic Current. Makes call to Ionic Model object to get current
    void getVoltageFromIonicCurrent(Epetra_Vector *Voltage, 
				    Epetra_Vector *GhostVoltage, const Real dT,
				    const Real dt_Ion, const int numNodes);

    /*! 
      Get Ionic Current Only. Voltage is not updated. Routine returns
      \f$\dfrac{dV}{dt}\f$. Voltage updates needs to be performed
      outside the class
      @param IC EpetraVector where dV/dt is stored
      @param Voltage Voltage value a \f$n^{th}\f$ time step
      @param dt_ion Time step for which Ionic Update needs to be 
      performed      
    */
    void getIonicCurrent(Epetra_Vector *IC, Epetra_Vector *Voltage, 
			 Epetra_Vector* GhostVoltage, const Real dt_ion);

    //! Apply stimulus Current for a region along an axis
    void setStimulusCurrent(const int Axis, const Real limit, 
			    const Real stim = 0.);

    /*! Apply stimulus current in circular region. Purkinje Nodes not 
      stimulated. Mostly this is used for S2 stimulus and hence this 
      condition. For purkinje stimulation use nodesets.
    */
    void setStimulusCurrent(Real xc[], const Real radius, const Real stim=0.);

    /*! Apply stimulus to nodes below a plane. Input is the normal to a plane
      and a point on the plane. All nodes whch lie below the plane will be
      stimulated
     */
    void setStimulusCurrent(Real xp[], Real xn[], const Real stim = 0.);

    //! Apply stimulus current for a node set
    void setStimulusCurrent(std::vector<int>& nodeSet, const Real istim=0.);

    /*!
      Set stimulus voltage for a nodeset. 
      \warning This needs to be used for Scroll breakup testing purposes 
      only. For regular S1-S2 protocol use setStimulusCurrent member function
      @param Voltage EpetraVector which has voltage values. Need to sent by 
      the solver class
      @param nodeSet Nodeset containing node ID's where voltage needs to be 
      reset
      @param value Voltage value to be used
    */
    void setStimulusVoltage(Epetra_Vector* Voltage, 
			    Epetra_Vector* GhostVoltage,
			    std::vector<int>&
			    nodeSet, const Real value = 0.);

    //! Get Nodal Capacitance
    Real getNodalCapacitance() const { return _C_m; }

    //! Set Nodal Capacitance
    void setNodalCapacitance(const Real C_m) {
      _C_m = C_m;
    }

    //! User provided values of diffusion
    void setDiffusion(Real D[3]) {
	for(int i = 0; i < 3; i++) _D[i] = D[i];
    }

    //! User provided values of Purkinje diffusion
    void setPurkinjeDiffusion( Real Dp[3] ) {
      for(int i = 0; i < 3; i++) _Dp[i] = Dp[i];
    }

    /*! 
      User provided values of Purkinje Diffusion in LV. We want to make 
      Purkinje conduction faster in LV compared to RV. By default the 
      diffusion for Purkinje fiber in LV and RV is the same. Mesh class
      provides input of whether a purkinje node is in LV or RV.
    */
    void setPurkinjeDiffusionInLV( Real DpLeft[3] ) {
      for(int i = 0; i < 3; i++) _DpLeft[i] = DpLeft[i];
    }

    //! Set User Specified Xi value
    void setUserXi(const Real Xi) {
      _userXi     = Xi;
      _userXiSet  = true;
    }

    //! get active contraction
    void getActiveContraction(std::vector<Real>& gamma, 
			      std::map<int, int>& nodeID,
			      const int numNodes);

    //! Compute active deformation tensor inverse
    std::vector<Tensor3D> computeActiveDeformationInverse(std::vector<Real>&
							  gamma);

    //! Get fiber direction as GaussPoints
    const std::vector< std::vector<Vector3D> >& getFiberDataAtGaussPoints() {
      return _Ibar;
    }



    // Debug
    std::vector<double > getInitialVoltage() {
      int nLocalNodes = _myMesh->getNumLocalNodes();
      std::vector<double > InitialVoltage(nLocalNodes, 0.0);
      for (uint i = 0; i < nLocalNodes; i++) {
	InitialVoltage[i] = _ionicModel[i]->getVoltage();
      }
      return InitialVoltage;
    }

    std::vector<double > getInitialGhostVoltage() {
      int nLocalNodes = _myMesh->getNumLocalNodes();
      int nGhostNodes = _myMesh->getNumGhostNodes();
      std::vector<double > InitialGhostVoltage(nGhostNodes, 0.0);
      for (uint i = 0; i < nGhostNodes; i++) {
	InitialGhostVoltage[i] = _ionicModel[i+nLocalNodes]->getVoltage();
      }
      return InitialGhostVoltage;
    }
    //End debug



    /*! 
      Get internal variables for all cell models at each node.
      @param numVar is a list which says number of internal variables at 
      each node.
      @param data is a 2D list which stores all the internal variables at
      each node.
    */
    void getStateVariables( std::vector<int>& numVar,
			    std::vector< std::vector<Real> >& data);

    /*!
      Set internal variables for each cell model. Used as part of restart
      capability.
    */
    void setStateVariables(const std::vector< std::vector<Real> >& data);

    //! Set diffusive smoothing on/off
    void setDiffusiveSmoothing(bool smoothing, const Real alpha = 0.1) {
      _useDiffusiveSmoothing = smoothing;
      _alpha = alpha;
    }

    bool useDiffusiveSmoothing() {
      return _useDiffusiveSmoothing;
    }

    Real getSmoothingCoeff() { return _alpha; }

    //! Gateway to the Body Class. 
    static CardiacBody* New(Mesh* const mesh, CellModelType Type);

    //! Gateway to the Body Class with Specified Internal Variables
    static CardiacBody* New(Mesh* const mesh, CellModelType Type, const std::vector<Real>& purkinjeStates, const std::vector<Real>& myocardiumStates);

  private:
    //! A Pointer to the Mesh Object
    Mesh*                _myMesh;
 
    //! Diffusion Stiffness. Needs a better name
    DiffusionStiffness   _diffusion;

    //! Full Mass Matrix
    DiffusionStiffness   _fullMass;

    //! Ionic Models at the nodes
    IonicModels          _ionicModel;

    //! Stimulus Current at Nodes
    std::vector<Real>    _stimulusCurrent;

    //! Capacitance value
    Real                 _C_m;

    //! Diffusion Values for the Tissue
   Real                 _D[3];

    //! Diffusion Values for the Purkinje Fiber
    Real                 _Dp[3];

    //! Diffusion Values for Purkinje in LV
    Real                 _DpLeft[3];

    //! userSpecified Xi
    Real                 _userXi;

    //! Is user Specified Xi set
    bool                 _userXiSet;

    //! Tensor defining active fiber componentes at Guass Points
    std::vector< std::vector<Vector3D> > _Ibar;

    //! Calculate internal force
    void getInternalForce(Material* material, const Mesh::QuadPoints &QuadVals,
			  const int start, const int end, 
			  vector< vector<Real> >& x, const int numCon, 
			  const int dim_t, blitz::Array<Real, 1>& fanal) ;

    //! Smoothing Capacitance matrix
    bool               _useDiffusiveSmoothing;
    Real               _alpha;
  };
}// namespace voom
#endif

