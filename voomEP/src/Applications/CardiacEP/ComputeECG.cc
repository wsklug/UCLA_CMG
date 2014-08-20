/*!
 Program to compute ECG from Voltage output fata. ECG is computed as
 \f[
 \text{ECG}=\int\int\intD_{ij}(\nabla V)_i(\nabla \frac{1}{R})_jdV
 \f]
 where R is distance of the lead to a point on the heart

 Shankarjee Krishnamoorthi 
 UCLA
*/
#include "ComputeECG.h"

//! A Small helper function
MeshType getMeshType(std::string meshType) {
    if (meshType =="MIX3D") return MIXEDMESH3D;
    if (meshType == "HEX") return HEXMESH;
    if (meshType == "TET") return TETMESH;
    if (meshType == "QUAD") return QUADMESH;
    if (meshType == "TRIA") return TRI3MESH;
    cerr << "Unknown Mesh Type provided. Allowed Types: ";
    cerr<< "QUAD, TRIA, HEX, TET (Case Sensitive)\n";
    exit(1);
}
 

int main(int argc, char** argv){
  string meshdir, outputPath, outDir;
  string modelName, pointsFile, meshType, nparts;
  int nPoints, quadOrder, dim;
  Real startStep, endStep, frequency;
  vector< vector<Real> > data;
  string fmt;
  
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << " parameters_file\n";
    exit(0);
  }
  
  ifstream inp(argv[1]);
  if (!inp.is_open()){
    cerr << "ERROR opening file: " << argv[1] << endl;
    exit(0);
  }
  
  inp >> modelName >> meshdir >> outputPath >> outDir >> quadOrder >> meshType;
  inp >> nPoints >> dim;
  for(int i = 0; i < nPoints ;i++){
    vector<Real> point;
    Real p;
    for(int j = 0; j < dim; j++) { inp >> p; point.push_back(p); }
    data.push_back(point);
  }
  inp >> nparts >> startStep >> endStep >> frequency >> fmt;
  inp.close();
  // Mpi Related Stuff 
  MPI_Init(&argc, &argv);
  Epetra_MpiComm* mpicomm = new Epetra_MpiComm(MPI_COMM_WORLD);

  if (mpicomm->MyPID() == 0 ) {
    cout << "***********************************************************\n";
    cout << "ModelName        : " << modelName << endl; 
    cout << "Voltage Output   : " << outputPath << endl;
    cout << "Output Location  : " << outDir << endl;
    cout << "Quadrature order : " << quadOrder << endl;
    cout << "Mesh Type        : " << meshType << endl;
    cout << "Number of Procs  : " << nparts << endl;
    cout << "Start Step       : " << startStep << endl;
    cout << "Stop Step        : " << endStep << endl;
    cout << "Frequency        : " << frequency << endl;
    cout << "***********************************************************\n";
  }
  
  
  // Create Mesh Object
  MeshType mType = getMeshType(meshType);
  Mesh *myMesh = Mesh::New( mpicomm, meshdir+"/"+modelName, quadOrder, mType );
  if (mpicomm->MyPID() == 0 ) cout << "Created Mesh Object\n";
  myMesh->setRadius ( 0.0275 );
  myMesh->Compute(); // Compute Shapefunction and derivatives
  
  // Compute Grad (1/R) for all leads
  vector< struct GradVals > gradOneByR;
  computeGradOneByR(myMesh, data, gradOneByR);
  if (mpicomm->MyPID() == 0) cout << "GradOneByR computed\n";

  
  // Compute Diffusion Tensor for all QuadPoints
  Real D[3] = {0.001,0.0005,0.00025};
  Real Dp[3] = { 0.0016, 0.0016, 0.0016};
  vector< Diffusion> diffusion_matrix;
  computeDiffusionTensor(myMesh, D, Dp, diffusion_matrix);
  if (mpicomm->MyPID() == 0 ) cout << "Diffusion Matrix Computed\n";
  
  const int numQuadPoints = diffusion_matrix.size();
  const int numNodes = myMesh->numberOfNodes();
  Real results[numNodes];
  int nparts_int = atoi(nparts.c_str());
  ifstream files[nparts_int];
  float nparts_float = nparts_int;
  string npartsfname = meshdir + "/" + modelName + ".metis.npart." + nparts;
  ifstream npartsfp;
  char outname[200];
  float myproc = mpicomm->MyPID();
  sprintf(outname, "%s/ECG_Data.%02.0f",outDir.c_str(), myproc );
  ofstream out(outname);
 
  // Computation of ECG
  for(float step = startStep; step < endStep; step = step + frequency) {
    if (mpicomm->MyPID() == 0)
      cout << "Processing Step: " << step << endl;
    npartsfp.open(npartsfname.c_str());
    for (int np=0; np<nparts_int; np++) {
      char fname[300];
      float step_float = step;
      float p_float=np;
      sprintf(fname, "%s/Voltage_time_%08.2f.%02.0f",outputPath.c_str(), 
	      step, p_float);
      
      files[np].open(fname);
      if ( !files[np].is_open() ) {
        cerr << "Error opening " << fname << endl;
      }
      if( (!files[np]) ) {
        cout<< mpicomm->MyPID() <<" upto "<<step<<endl;
        MPI_Finalize();
        return 0;
      }
    } // np loop

    // Real results from results file
    int p; double v;
    for(int i=0; i < numNodes;i++){
      npartsfp >> p;
      if (fmt=="ASCII")  files[p] >> results[i];
      else  {
        files[p].read((char*)(&v), sizeof(double));
        results[i] = v;
      }
    } // i loop
    // Compute Voltage at quadPoints
    vector< GradVoltage > voltage;
    computeVoltageAtQuadPoints(myMesh, results, voltage);

    // Compute ECG with all data
    Real* ECG = computeECG(diffusion_matrix, voltage, gradOneByR,
			  myMesh->getDimension());
    out << step << "\t";
    for(int p = 0; p < nPoints; p++)
      out << ECG[p] << "\t";
    out << endl;
    // Close all Files
    npartsfp.close();
    for (int p=0; p < nparts_int; p++) files[p].close();
    delete [] ECG;
    
  } // step loop
  
  out.close();
  // Finalize MPI
  MPI_Finalize();
  return 0;
}
