. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
rm log.decomposePar
rm log.simpleFoam
rm log.reconstructParMesh
rm log.reconstructPar

runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar