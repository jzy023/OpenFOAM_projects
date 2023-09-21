blockMesh
surfaceFeatures
snappyHexMesh -overwrite
simpleFoam
foamToVTK
# decomposePar
# mpirun -np 4 DPMFoam -parallel > log &

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
