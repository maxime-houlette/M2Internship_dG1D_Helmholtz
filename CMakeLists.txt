cmake_minimum_required(VERSION 3.15)

# Quelques bonnes pratiques reprises de votre ancien fichier
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_VERBOSE_MAKEFILE OFF)

# Configuration C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

project(Helmholtz1D_dG LANGUAGES CXX)

# Recherche de Eigen (Assurez-vous qu'il est installé via Homebrew)
find_package(Eigen3 REQUIRED)

###############################################################################
# 1. LA BIBLIOTHÈQUE PRINCIPALE (Le coeur de votre solveur)
###############################################################################
add_library(dg_core
    src/dg/mesh.cpp
    # src/dg/solver.cpp
)

# Dossiers d'inclusion pour la bibliothèque
target_include_directories(dg_core PUBLIC 
    include
    include/dg
)

# Options de compilation strictes
target_compile_options(dg_core PRIVATE 
    -Wall -Wextra -Wpedantic
)

# Lier Eigen
target_link_libraries(dg_core PUBLIC Eigen3::Eigen)

###############################################################################
# 2. L'EXÉCUTABLE DE TEST (Le main)
###############################################################################
add_executable(main_exec
    src/bin/main.cpp
)

target_compile_options(main_exec PRIVATE 
    -Wall -Wextra -Wpedantic
)

# On lie l'exécutable à notre propre bibliothèque dg_core
target_link_libraries(main_exec PRIVATE dg_core)