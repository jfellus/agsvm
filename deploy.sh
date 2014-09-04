#!/bin/bash


scp gossip_pca jerofell@plaintcontrix:/users/jerofell/agpca
scp -r src jerofell@plaintcontrix:/users/jerofell/agpca
scp CMakeLists.txt jerofell@plaintcontrix:/users/jerofell/agpca

echo "deployed !!"
