
# #!/bin/sh
# # cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
restore0Dir

runApplication blockMesh
runApplication surfaceFeatureExtract
runApplication decomposePar
runParallel snappyHexMesh -overwrite
runApplication reconstructParMesh -constant
runApplication transformPoints -scale '(0.001 0.001 0.001)'
