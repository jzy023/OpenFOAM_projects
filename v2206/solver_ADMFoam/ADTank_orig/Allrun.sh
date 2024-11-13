
# #!/bin/sh
# # cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
# ------------------------------------------------------------------------------

restore0Dir
# Mesh -------------------
runApplication blockMesh
runApplication surfaceFeatureExtract
runApplication decomposePar
runParallel snappyHexMesh -overwrite
runApplication reconstructParMesh -constant
runApplication transformPoints -scale '(0.001 0.001 0.001)'
# Sim --------------------
rm -rf proc*
rm log.de*
rm log.re*
runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar
rm -rf proc*
# Post -------------------
touch foam.foam


# =====================================================
# blockMesh
# surfaceFeatureExtract
# snappyHexMesh -overwrite # > log.SHM 

# transformPoints -scale '(0.001 0.001 0.001)'

# decomposePar
# mpirun -np 4 simpleFoam -parallel
# reconstructPar

