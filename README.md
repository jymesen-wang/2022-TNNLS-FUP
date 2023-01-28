# FUP

To use the function "FUP" or "OFUP", please follow the input/output format:

 [ W,M,P,objva,r ] = FUP( X,P,m,d1,K,iteration,a)
 [ W,M,P,objva,r ] = OFUP(X,P,m,d1,K,iteration,lamda,a)
 
X  ：  d*n, sample matrix, d is the dimension of the sample, n is the number of the number
P  ：  n*m, similarity matrix between samples and representitive points, m is the number of the representative points
m  :   number of the representative points
d1 :   dimension of subspace, d1<d 
K  ：  number of nearest neighbors
iteration : iteration times
a  :   a==1, isometric sampling initialization for M， else K-means initialization for M
W  ：  d*d1 projection matrix
objva: 1*iteration, the objective function values of all iterations
r  ：  optimized regularization parameter
lamda ： hyperparameter


Please make sure that the documents Eu2_distance.m, EProjSimplex_new.m and ClusteringMeasure.m are in the same folder as FUP.m and PFUP.m

Use the codes, please cite
Wang J, Wang L, Nie F, et al. Fast Unsupervised Projection for Large-Scale Data[J]. IEEE Transactions on Neural Networks and Learning Systems, 33(8), pp. 3634-3644, 2022.

If you have any questions, please connect wanglinjun@mail.nwpu.edu.cn
