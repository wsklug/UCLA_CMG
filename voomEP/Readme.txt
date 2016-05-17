Steps for Compiling voomEP:
1. First run the autoconfigure tools command "autoreconf -vif"
2. Next, Run the configure command with the appropriate paths for the the TVMET and BLITZ libraries:
"./configure CXX=mpicxx F77=gfortran --with-blitz-root=BLITZ_LIBRARY_PATH --with-tvmet-root=TVMET_LIBRARY_PATH"
3. Run the "make" command to compile the voomEP libraries.
4. Run the same commands, i.e. autoreconf, configure, make, in the src/Applications/CardiacEP folder to create the EPSimulator and ComputeECG executables.
