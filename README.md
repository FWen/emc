# Efficient Maximum Consensus Algorithms for Robust Fitting in Computer Vision

This code is used to reproduce the results in "F. Wen, R. Ying, Z. Gong, and P. Liu, Efficient Algorithms for Maximum Consensus
Robust Fitting, IEEE Transactions on Robotics, 2019". 

To use the code, firstly unzip 'src.zip'.

This code is modified from the code of Huu Le at https://www.researchgate.net/publication/320707327demo_pami, details please see the paper "H. Le, T. J. Chin, A. Eriksson, and D. Suter, “Deterministic approx-imate methods for maximum consensus robust fitting,” arXiv preprint, arXiv:1710.10003, 2017."

Note that, some codes of Huu Le are directly coped here from https://www.researchgate.net/publication/320707327demo_pami to facilitate the ease of use for interested readers who want to reproducing the results in our paper. The coped codes (contained in the '/src' folder) include the EP method, RANSAC method and their dependency. We copy them here only for academic use purpose to illustrate the comparison results of the algorithms in our paper.

SeDuMi is need in solving the LP problems, which is available at http://sedumi.ie.lehigh.edu/, we used SeDuMi 1.3.

# Application in SLAM
The “C++ for Application in ORB-SLAM2” folder contains the C++ implementation of the ADMM algorithm for homography and fundamental matrix estimation, which can be directly used in the popular open source ORB-SLAM2 system (https://github.com/raulmur/ORB_SLAM2). 

We have tested the code in ORB-SLAM2 in a monocular example using a seuence from the KITTI dadaset.
Users can choose RANSAC or ADMM via setting a parameter in "Initializer.h", e.g. USE_RANSAC or USE_ADMM.
