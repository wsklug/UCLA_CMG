//-*-C++-*-
#ifndef _EPSolver_h_
#define _EPSolver_h_
#include "voom.h"

#include "Epetra_ConfigDefs.h"
#ifdef HAVE_MPI
#include "mpi.h"
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

// Aztec gives us Iterative Solvers
#include "AztecOO_config.h"
#include "AztecOO.h"
#include "AztecOO_Operator.h"


#include "Mesh.h"
#include "CardiacBody.h"
#include "Material.h"


namespace voom{
  //! Specify the Solution technique to be used
  enum Method {EXPLICIT, BACKWARDEULER, CRANKNICHOLSON};

  //! Typedefs
  typedef std::vector<Mesh *> MeshContainer;
  typedef std::vector<CardiacBody *> BodyContainer;
  typedef std::vector<Real > SmoothCoeff;
  typedef Epetra_MultiVector                MV;
  typedef Epetra_Operator                   OP;
  

  /*! \brief
    A solver class implemented using Epetra functionalities. This base 
    class implements the basis functions of the Solver
    - Create typical Epetra Objects to be used by sub classes
    - Output of Result
    - Computing Lumped Capacitance. This is also placed in an Epetra
    Vector for use
  */
  class EPSolver {
  private:
    //! Tolerance
    Real           _tolerance;
  protected:
    //! Mesh Pointer
    MeshContainer   _myMesh;
    
    //! Cardiac Body pointer
    BodyContainer   _myBody;

    //! Index Map
    int*           _indexMap;

    //! Number of Local Nodes
    int             _nLocalNodes;

    //! Number of ghost Nodes
    int             _nGhostNodes;

    //! Localnodes ID
    int*            _localNodesID;
    int*            _dlocalNodeID;

    //! Ghostnodes ID
    int*            _ghostNodesID;
    int*            _dghostNodeID;

    //! Epetra Objects
    Epetra_MpiComm* _comm;
    Epetra_Map*     _sourceMap;
    Epetra_Map*     _targetMap;
    Epetra_Export*  _exporter;

    //! Relevant to displacement MPI
    Epetra_Map*     _dispSourceMap;
    Epetra_Map*     _dispTargetMap;
    Epetra_Export*  _dispExporter;

    //! Diffusive Current
    Epetra_Vector* _diffusiveCurrent;

    //! Voltage at Nodes
    Epetra_Vector* _Voltage; 

    //! Lumped Capacitance at Nodes
    Epetra_Vector* _capacitance;

    //! Ghost Diffusive Current
    Epetra_Vector* _ghostDiffusiveCurrent;

    //! Nodal Ionic Currnet
    Epetra_Vector* _ionicCurrent;

    //! Ghost Voltage
    Epetra_Vector* _ghostVoltage;

    //! Ghost Capacitance
    Epetra_Vector* _ghostCapacitance;

    //! Ionic Time Step
    Real           _ionicTimeStep;

    //! Solve Time
    Real           _dT;

    //! Smoothing Coefficients
    SmoothCoeff    _coeff;

    //! Verbose Output
    bool           _verbose;

    //! Displacement Vector (All Components rolled into one vector)
    Epetra_Vector* _displacement;
    Epetra_Vector* _U;

    //! Ghost displacement
    Epetra_Vector* _ghostDisplacement;
    Epetra_Vector* _ghostU;
    
    //! External Force Vector
    Epetra_Vector* _force;

    //! Internal Force Vector
    Epetra_Vector* _fInternal;
    Epetra_Vector* _ghostfInternal;

    //! Material properties to be used
    Material*          _material;
    
    //! Check if displacement values exist
    bool              _displacementData;

    //! Max NR iteration
    int               _maxIter;

    //! Max CG Iteration
    int               _maxCGIter;

    //! CG Solver Tolerance
    Real              _cgSolverTolerance;

  public:
    /*! Constructor
      @param comm MPI object. Should be created in the application
      @param mesh A vector of meshes. For all finite element applications
      we need only one mesh object. This will be used in future for 
      MeshFree applications
      @param cardiacBody A vector of cardiacbody. Similar to mesh finite 
      element applications use only one body. Meshfree or some other
      application might need more than one body
      @param coeff Vector of Real value. This is list of smoothing
      coefficients. Results from different bodies/meshes are weighted 
      using these values
      @param mechanicsSolve. Informs solver whether to solve 
      mechanics or not
    */
    EPSolver(Epetra_MpiComm *comm, const MeshContainer& mesh, 
	     const BodyContainer& cardiacBody, const SmoothCoeff& coeff);
    
    //! Destructor
    ~EPSolver() {
      delete _targetMap; delete _sourceMap; delete _exporter;
      delete _diffusiveCurrent; delete _Voltage;
      delete _ghostDiffusiveCurrent; delete _ghostVoltage;
      delete _capacitance; delete _ghostCapacitance;
      delete _ionicCurrent; delete _force;
      delete _dispTargetMap; delete _dispSourceMap; delete _U; 
      delete _dispExporter; delete _displacement; delete _ghostDisplacement;
      delete _ghostU; delete _fInternal; delete _ghostfInternal;
    }

    /*! 
      Pure virtual function. Will be implemented in derived classes
      @param dtDiff Diffusion time step. Default is 0.1 ms
     */
    virtual void Solve(const Real dtDiff = 0.1) = 0;

    //! Set Ionic Time Step
    void setIonicTimeStep(const Real dt) { _ionicTimeStep = dt; }

    //! Set Solve Time
    void setSolveTime(const Real dT) { _dT = dT; }

    //! Compute Lumped Capacitance
    void computeCapacitance();
    
    /*! 
      Write Output (Binary File)
      @param time Time at which output is written. 
      @param path Directory where output files are written
    */
    void WriteOutput(const Real time, std::string path = "");

    //! Initialize the Voltage Epetra vector to some value
    void InitializeVoltage(const Real value){
	_Voltage->PutScalar( value );
	_ghostVoltage->PutScalar( value );
    }



    void InitializeVoltage(const std::vector<Real > & localValues, 
			   const std::vector<Real > & ghostValues) {
      Real *_extractedVoltage, *_extractedGhostVoltage;
      _Voltage->ExtractView(&_extractedVoltage);
      _ghostVoltage->ExtractView(&_extractedGhostVoltage);

      for(uint i = 0; i < localValues.size(); i++)
	_extractedVoltage[i] = localValues[i];
      for(uint i = 0; i < ghostValues.size(); i++)
	_extractedGhostVoltage[i] = ghostValues[i];
    };



    /*!
      Set Voltage value at some nodes. \warning Use this for debugging 
      S2 protocal. For normal S1-S2 use body class member functions
    */
    void setVoltage(std::vector<int>& nodeSet, const Real value = 0.) {
      _myBody[0]->setStimulusVoltage(_Voltage, _ghostVoltage, nodeSet,
				     value);
    }

    //! Verbose output
    void setVerbose(const bool cond = true) {
      _verbose = cond;
    }

    //! Set Material
    void setMaterial(Material *material) {
      _material = material;
    }

    //! Consistency Check
    void consistencyCheck();

    //! Mechanics Solve. Returns number of NR iterations
    int solveMechanics();

    //! Set tolerance for mechanics solver
    void setTolerance(const Real tolerance){
      _tolerance = tolerance;
    }

    //! Set max iteration value
    void setMaxIteration(const int maxIter) {
      _maxIter = maxIter;
    }

    //! Set max iteration for linear solver
    void setMaxLinSolveIteration(const int maxIter) {
      _maxCGIter = maxIter;
    }

    //! Set linear solver tolerance
    void setLinearSolverTolerance(const Real toler){
      _cgSolverTolerance = toler;
    }

    /*! 
      Impose Essential Boundary Condition. For a given NodalID and 
      the degree of freedom we do not have any entries added in the
      stiffness matrix. We place a value of 1. in the diagonal and the
      corresponding entry for the force vector we place the constraint
      value
    */
    void imposeEssentialBoundaryConditions(Epetra_FECrsMatrix *K);

    /*!
      Impose natural boundary conditions. Relevant external force components
      are added to the external force Epetra_Vector
    */
    void imposeNaturalBoundaryConditions();

    //! Write Restart data. Write cell model internal data
    void writeRestart(const Real time, std::string path);

    /*! 
      Read Restart data - Restart capability
      Read's voltage data, cell model internal data also.
    */
    void readRestart(const Real time, std::string path);

    /*!
      Gateway to Solver Class. Use this function to create Solver objects.
      @param comm Epetra_MpiComm objects
      @param mesh Vector of meshes
      @param cardiacBody Vector of cardiacBody
      @param coeff Vector of Real
      @param SolverType Solver Type to be used as defined in Method
      
    */
    static EPSolver* New(Epetra_MpiComm *comm, const MeshContainer& mesh,
			 const BodyContainer& cardiacBody, 
			 const SmoothCoeff& coeff, Method SolverType);
  };
}
#endif
