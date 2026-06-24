#!/bin/bash

# Création du dossier Debug
echo "Configuring DEBUG version..."
echo "----------------------------"
mkdir -p build/debug
cd build/debug
cmake ../../ -DCMAKE_BUILD_TYPE=Debug
cd ../../
echo ""

# Création du dossier Release
echo "Configuring RELEASE version..."
echo "------------------------------"
mkdir -p build/release
cd build/release
cmake ../../ -DCMAKE_BUILD_TYPE=Release
cd ../../
echo ""