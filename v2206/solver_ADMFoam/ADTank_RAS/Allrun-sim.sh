. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
rm -rf proc*
rm -rf post* 
rm log*
find . -type d -regextype posix-extended -regex '\./[0-9]+(\.[0-9]+)?' -exec rm -rf {} +

restore0Dir
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar
rm -rf proc*
# Post -------------------
touch foam.foam