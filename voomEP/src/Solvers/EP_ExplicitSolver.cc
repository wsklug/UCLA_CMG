#include "EP_ExplicitSolver.h"

namespace voom{
  void EP_ExplicitSolver::Initialize() {
    if (_comm->MyPID() == 0)
      std::cout << "** Begin Constructing System Matrices " << std::endl;
    // Get _D assembled globally. Getting DV is faster with Epetra Calls
    _myBody[0]->getDiffusionMatrix(_D, -1.);
    _D->GlobalAssemble();
    _check = false;
    if (_comm->MyPID() == 0)
      std::cout << "** Finished Constructing System Matrices" << std::endl;
  }
  
  /*!
    The steps performed are
   - Solve for diffusive current for a time dT/2
   - Solve for ionic current for time dT taking _ionicTimestep as 
   incremental time step
   - Solve for diffusive current for a time dT/2
   \note We use Operator Splitting
  */    
  void EP_ExplicitSolver::Solve(const Real dtDiff) {
    int numNodes = _myMesh[0]->getNumLocalNodes() + 
      _myMesh[0]->getNumGhostNodes();
    Real *voltage, *ghostVoltage, *diffCurrent, *ghostDiffCurrent, 
      *capacitance, *ghostCapacitance;
    const int* localNodes  = _myMesh[0]->getLocalNodeID();
    const int* ghostNodes  = _myMesh[0]->getGhostNodeID();
    const int* indexMap    = _myMesh[0]->getIndexMap();
    const int nLocalNodes  = _myMesh[0]->getNumLocalNodes();
    const int nGhostNodes   = _myMesh[0]->getNumGhostNodes(); 
    const Real TOLER = 1E-6;
    int id;
    // Choose lower value
    const Real dtDiffusion = ( dtDiff < 0.5 * _dT ) ? dtDiff : 0.5 * _dT; 
    Real solveTime = 0.5 * _dT;
    Real stepTime = dtDiffusion;
    
    // Check if D Matrix has been built
    if (_check) Initialize();

    // Extract view to do stuff faster
    _Voltage->ExtractView(&voltage);
    _ghostVoltage->ExtractView(&ghostVoltage);
    _capacitance->ExtractView(&capacitance);
    _ghostCapacitance->ExtractView(&ghostCapacitance);
    _diffusiveCurrent->ExtractView(&diffCurrent);
    _ghostDiffusiveCurrent->ExtractView(&ghostDiffCurrent);
    
    
    // Step 1: Compute the diffusive Current from body
    // Do this inside a loop
    while (solveTime > TOLER) {      
      // Step 0: Zero out diffusive current Vector for first body
      _diffusiveCurrent->PutScalar(0.0);
      _ghostDiffusiveCurrent->PutScalar(0.0);
      
      // Finish for first body
      _D->Multiply(false, *_Voltage, *_diffusiveCurrent);

      _ghostVoltage->Import(*_Voltage, *_exporter, Insert);

      // Import Export epetra stuff
      _diffusiveCurrent->Export(*_ghostDiffusiveCurrent, *_exporter, Add);
      _ghostDiffusiveCurrent->Import(*_diffusiveCurrent, *_exporter, Insert);
      

      // Update Voltage from first body
      for(int i = 0; i < numNodes; i++) {
	if ( i < nLocalNodes) id = indexMap[ localNodes[i] ] - 1;
	else id = - indexMap[ ghostNodes[i - nLocalNodes] ] - 1;
	if ( i < nLocalNodes ) 
	  voltage[id] += stepTime*diffCurrent[id]/capacitance[id];
	else 
	  ghostVoltage[id] += stepTime*ghostDiffCurrent[id]/
	    ghostCapacitance[id];
      }
      
      solveTime -= stepTime;
      stepTime = ( stepTime < solveTime) ? stepTime: solveTime;
    } // While Loop

    // Step 2: Get Voltage from Ionic current only first body only
    _myBody[0]->getVoltageFromIonicCurrent(_Voltage, _ghostVoltage, _dT, 
					   _ionicTimeStep, numNodes);
    solveTime = 0.5 * _dT;
    stepTime = dtDiffusion;
    while (solveTime > TOLER) {
      // Step 0: Zero out diffusive current Vector for first body
      _diffusiveCurrent->PutScalar(0.0);
      _ghostDiffusiveCurrent->PutScalar(0.0);
      
      // Finish for first body
      _D->Multiply(false, *_Voltage, *_diffusiveCurrent);

      _ghostVoltage->Import(*_Voltage, *_exporter, Insert);
      
      // Import Export epetra stuff
      _diffusiveCurrent->Export(*_ghostDiffusiveCurrent, *_exporter, Add);
      _ghostDiffusiveCurrent->Import(*_diffusiveCurrent, *_exporter, Insert);
      

      // Update Voltage from first body
      for(int i = 0; i < numNodes; i++) {
	if ( i < nLocalNodes) id = indexMap[ localNodes[i] ] - 1;
	else id = - indexMap[ ghostNodes[i - nLocalNodes] ] - 1;
	if ( i < nLocalNodes ) 
	  voltage[id] += stepTime*diffCurrent[id]/capacitance[id];
	else 
	  ghostVoltage[id] += stepTime*ghostDiffCurrent[id]/
	    ghostCapacitance[id];
      }
      
      solveTime -= stepTime;
      stepTime = ( stepTime < solveTime) ? stepTime: solveTime;
    }// while Loop
  } // End of Solve routine

} // namespace
