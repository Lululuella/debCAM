context("Valid CAM")

K <- 3
S <- matrix(2 ^ rnorm(1000*K, 5, 2.5), K)
# A <- matrix(c(0.2,0.4,0.6,0.8,
#               0.8,0.6,0.4,0.2), ncol = K)
A <- matrix(c(0.1,0.3,0.5,0.7,0.9,
              0.2,0.4,0.5,0.1,0.0,
              0.7,0.3,0.0,0.2,0.1), ncol=K)
X <- log2(A%*%S)  + rnorm(1000*nrow(A), 0 , 0.2)
X <- t (2 ^ X)
rownames(X) <- seq_len(nrow(X))

K_candi <- 2:4


CAMResult <- CAM(X, K=K_candi, corner.strategy=2, dim.rdc=10,
                 thres.low=0.3, thres.high=0.9, cluster.method = 'K-Means',
                 cluster.num=20, MG.num.thres=5, lof.thres=0.02,
                 qhull.enable=T, cores=20, seed=111)



CAMResult1 <- CAM(X, K=K_candi, corner.strategy=2, dim.rdc=10,
                 thres.low=0.3, thres.high=0.9, cluster.method = 'apcluster',
                 cluster.num=20, MG.num.thres=5, lof.thres=0.02,
                 qhull.enable=T, cores=20, seed=111)


MDLres <- MDL(CAMResult)

plot(MDLres)
plot(MDLres, data.term = TRUE)


plot(MDL(CAMResult1))
plot(MDL(CAMResult1), data.term = TRUE)
