# blockMesh
# surfaceFeatureExtract
# snappyHexMesh -overwrite # > log.SHM 

# transformPoints -scale '(0.001 0.001 0.001)'

# decomposePar
# mpirun -np 4 simpleFoam -parallel
# reconstructPar

# =====================================================
# #!/bin/sh
# cd ${0%/*} || exit 1    # Run from this directory

# # Source tutorial run functions
# . $WM_PROJECT_DIR/bin/tools/RunFunctions

# application=$(getApplication)

# runApplication blockMesh
# runApplication surfaceFeatures
# runApplication snappyHexMesh -overwrite
# runApplication decomposePar
# runParallel $(getApplication)


# =====================================================
# #!/bin/sh
# # cd "${0%/*}" || exit                                # Run from this directory
. ${WM_PROJECT_DIR:?}/bin/tools/RunFunctions        # Tutorial run functions
#------------------------------------------------------------------------------
runApplication blockMesh
runApplication surfaceFeatureExtract
runApplication decomposePar
runParallel snappyHexMesh -overwrite
runApplication transformPoints -scale '(0.001 0.001 0.001)'
runApplication reconstructParMesh -constant
rm -r proc*

rm log.decomposePar
# rm log.simpleFoam
rm log.reconstructParMesh
# rm log.reconstructPar

runApplication decomposePar
runParallel $(getApplication)
runApplication reconstructParMesh -constant
runApplication reconstructPar