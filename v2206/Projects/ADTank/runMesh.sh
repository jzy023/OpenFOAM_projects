

# blockMesh
# surfaceFeatureExtract
# decomposePar #-force
# snappyHexMesh -overwrite # > log.SHM 
# # transformPoints -scale '(0.001 0.001 0.001)'


# #!/bin/sh
# # cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
runApplication blockMesh
runApplication surfaceFeatureExtract
runApplication decomposePar
runParallel snappyHexMesh -overwrite
runApplication transformPoints -scale '(0.001 0.001 0.001)'


# runApplication reconstructParMesh -constant
# # runApplication reconstructPar