This folder contains scripts written to reproduce results obtained in DeepSite's article [1].
For initial implementation, I use a smaller subset of scPDB v.2017 database [2].

I've added data preprocessing, code for which is not presented in the original article.

As a starting point I use original neural network implementation (taken from Supplementary materials) and plan to experiment with parameters, network architecture, probably set of chemical descriptors used for feature extraction.

These and other scripts (written by my student as a part of Bioinformatics Institute research project - I've planned to do this by myself, but after a while I thought that the best way to have a hard deadline for this is to outsource) can be found at [bitbucket repo](https://bitbucket.org/one_project_acc/bideepsite/src).

The original idea also included part 2 - "to experiment with capsule networks". 
This part I still plan to implement myself next autumn and to add to this subfolder.

[1]: J. Jiménez, S. Doerr, G. Martínez-Rosell, A. S. Rose, G. De Fabritiis; DeepSite: protein-binding site predictor using 3D-convolutional neural networks, Bioinformatics, Volume 33, Issue 19, 1 October 2017, Pages 3036–3042, DOI: 10.1093/bioinformatics/btx350.

[2]: Jérémy Desaphy, Guillaume Bret, Didier Rognan, Esther Kellenberger; sc-PDB: a 3D-database of ligandable binding sites—10 years on, Nucleic Acids Research, Volume 43, Issue D1, 28 January 2015, Pages D399–D404, DOI 10.1093/nar/gku928.
