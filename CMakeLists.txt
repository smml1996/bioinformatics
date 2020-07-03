cmake_minimum_required(VERSION 3.15)
project(BioInformatics)

set(CMAKE_CXX_STANDARD 14)

add_executable(BioInformatics main.cpp SequenceAligner.cpp SequenceAligner.h NeedlemanWunsch.cpp NeedlemanWunsch.h SmithWaterman.cpp SmithWaterman.h MSA.cpp MSA.h ProgressiveMSA.cpp ProgressiveMSA.h Clustering.cpp Clustering.h DissociativeClustering.cpp DissociativeClustering.h UPGMA.cpp UPGMA.h NeighbourJoin.cpp NeighbourJoin.h AdditiveTree.cpp AdditiveTree.h SecondaryStructurePredictor.cpp SecondaryStructurePredictor.h)