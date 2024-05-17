. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
rm log.de*
rm log.re*
rm -rf proc*
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar