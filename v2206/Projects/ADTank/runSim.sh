. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar