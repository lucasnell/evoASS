# setwd('C:/Users/Tobin Northfield/Documents/EcoEvo/Community coevolution/2018 code/Multiple communities')
#
# library(scatterplot3d)

# Rcpp::sourceCpp("tony_sim_code/quant_gen_cr.cpp")


sc <- function(x1, x2, precision = 10^-2) {
    flag <- FALSE
    if (dim(x1)[1] == dim(x2)[1]) {
        n <- dim(x1)[1]
        flag <- TRUE
        i <- 0
        while (i < n & flag) {
            flag <- FALSE
            i <- i + 1
            for (j in 1:n) {
                if (mean((x1[i, ] - x2[j, ])^2) < precision^2) flag <- TRUE
            }

        }
    }
    return(flag)
}


# collist <- c("black", "blue", "dark green", "red", "green")
# pchlist <- c(19, 22, 24, 23, 4, 8)

ntrials <- 50
n <- 20  # number of resources
p <- 20  # number of consumers
q <- 3   # number of traits
precision <- 0.1


density.threshold <- 1e-5

maxtime <- 2e+06
starttime <- 1000

traj.flag <- FALSE
int <- 0.01 * maxtime

# turns evolution "on" with TRUE\n
evoN <- TRUE
evoP <- TRUE

# additive genetic variances:
sig2N <- 1e-10
sig2P <- 1e-10

sdN <- 0
sdP <- 0
sdV <- 0
sdU <- 0

if (evoN) {
    sig2N <- 0.05
}
if (evoP) {
    sig2P <- 0.05
}

# For resource...
r <- 0.05  # growth rate
a <- 2     # baseline competition coefficient
f <- 0.1   # cost of defensive traits
b <- 0.05  # coefficient for attack rate
# For consumer...
cc <- 1    # growth rate
m <- 0.01  # mortality rate
g <- 0.01  # cost of offense

# Non-additive tradeoffs for increasing >1 traits
etaN <- 0.3
etaP <- 0.2
C <- matrix(etaN, nrow = q, ncol = q)
diag(C) <- 1
D <- matrix(etaP, nrow = q, ncol = q)
diag(D) <- 1

initN <- 0.4
initP <- 0.1
initV <- 0.1
initU <- 1

# # These keep track of which species are extant
# Nsp.ident.vec <- rep(0, 10)
# Nsp.ident.vec[1] <- 1
# Psp.ident.vec <- rep(0, 10)
# Psp.ident.vec[1] <- 1
# extant.N <- 1:n
# extant.P <- 1:p
#
# etalist <- c(-.05, 0, .05, .1, .15, .2)


deltavec <- c(0.001, 0.01, 0.1, 0.3)


# ntrials <- ntrials / 10
maxtime <- maxtime / 100


# start - delta loop ----
for (delta.step in 1:length(deltavec)){
    # delta.step = 1

    delta <- deltavec[delta.step]

    preylist <- NULL
    predlist <- NULL
    num.preylist <- NULL
    num.predlist <- NULL

    # start - trial loop ----
    t0 <- Sys.time()
    for (trial in 1:ntrials) {
        # trial = 1
        N <- matrix(initN, nrow = n, ncol = 1)
        P <- matrix(initP, nrow = p, ncol = 1)
        V <- matrix(initV, nrow = n, ncol = q)
        U <- matrix(initU, nrow = p, ncol = q)

        # Nlist <- array(0, c(n, maxtime))
        # Plist <- array(0, c(p, maxtime))
        # Vlist <- array(0, c(n, q, maxtime))
        # Ulist <- array(0, c(p, q, maxtime))

        for (time in 1:starttime) {
            # time = 1

            A <- a + f * diag(V %*% C %*% t(V))
            B <- b * exp(-(V^2) %*% t(U^2))
            M <- m + g / diag(U %*% D %*% t(U))

            Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
            Pt <- P * exp(cc * t(B) %*% N - M)
            Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V *
                                   exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
            Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*%
                                   (V^2 * array(N, c(n, q))) + 2 * g *
                                   array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))

            N <- Nt
            P <- Pt
            V <- abs(Vt)
            U <- abs(Ut)
        }

        # perturbation
        # N <- exp(rnorm(n = n)) * N
        # P <- exp(rnorm(n = n)) * P
        V <- matrix(exp(rnorm(n = n * q, mean = 0, sd = delta)), nrow = n, ncol = q) * V
        U <- matrix(exp(rnorm(n = p * q, mean = 0, sd = delta)), nrow = p, ncol = q) * U

        for (time in 1:maxtime) {

            if (time == maxtime/2) {
                V <- matrix(exp(rnorm(n = n * q, mean = 0, sd = .3)),
                            nrow = n, ncol = q) * V
                U <- matrix(exp(rnorm(n = p * q, mean = 0, sd = .3)),
                            nrow = p, ncol = q) * U
            }
            A <- a + f * diag(V %*% C %*% t(V))
            B <- b * exp(-(V^2) %*% t(U^2))
            M <- m + g/diag(U %*% D %*% t(U))

            Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
            Pt <- P * exp(cc * t(B) %*% N - M)
            Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V *
                                   exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
            Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*%
                                   (V^2 * array(N, c(n, q))) + 2 * g *
                                   array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))

            N <- Nt
            P <- Pt
            V <- abs(Vt)
            U <- abs(Ut)

            # catch low values of U
            # U[U < 10^-4] <- 10^-4

            # Nlist[, time] <- N
            # Plist[, time] <- P
            # Vlist[, , time] <- V
            # Ulist[, , time] <- U
        }

        # Calculated final fitnesses and selection pressure to see if we're at equilibrium
        WN <-  exp(r * (1 - A * sum(N) - B %*% P))
        WP <- exp(cc * t(B) %*% N - M)


        Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V *
                               exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
        SV <- (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V *
                   exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))


        Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*%
                               (V^2 * array(N, c(n, q))) + 2 * g *
                               array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))
        SU <- (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*% (V^2 * array(N, c(n, q))) +
                   2 * g * array(1/diag(U %*% D %*% t(U))^2, c(p, q)) * (U %*% D))


        Fitness.N <- prod(WN[N>density.threshold])
        Fitness.P <- prod(WP[P>density.threshold])

        Selection.V <- sum(SV[N>density.threshold,])
        Selection.U <- sum(SU[P>density.threshold,])


        # if (traj.flag){
        #
        #     pdf(paste("traj",trial,".pdf", sep=""), width=8, height=6)
        #     par(mfcol = c(2, 1 + q))
        #
        #     matplot(t(Nlist[, int * (1:(dim(Nlist)[2]/int))]), type = "l",
        #             ylab = "prey density")
        #     matplot(t(Plist[, int * (1:(dim(Plist)[2]/int))]), type = "l",
        #             ylab = "pred density")
        #
        #     for (i in 1:q) {
        #         matplot(t(abs(Vlist[, i, int * (1:(dim(Vlist)[3]/int))])), type = "l",
        #                 ylab = paste("prey trait", i))
        #         matplot(t(abs(Ulist[, i, int * (1:(dim(Ulist)[3]/int))])), type = "l",
        #                 ylab = paste("pred trait", i))
        #     }
        #     dev.off()
        # }

        ##############################*
        # truncate any values of V or U that are over 3 to 3 so that unique species
        # can be identified
        V[V>3] <- 3
        U[U>3] <- 3
        ##############################*

        # find unique species.
        unique.prey <- 1:n
        for (i in 1:(n - 1)) for(j in (i+1):n) {
            if (sc(matrix(V[i,], nrow=1), matrix(V[j,], nrow=1), precision)) {
                unique.prey[j] <- 0
            }
        }

        unique.prey <- unique.prey[unique.prey > 0 & N > density.threshold]
        U[P < density.threshold,] <- 100
        unique.pred <- 1:p
        for (i in 1:(p - 1)) for (j in (i+1):p) {
            if (sc(matrix(U[i,], nrow=1), matrix(U[j,], nrow=1), precision)) {
                unique.pred[j] <- 0
            }
        }
        unique.pred <- unique.pred[unique.pred > 0 & P > density.threshold]

        num.prey <- length(unique.prey)
        if(num.prey == 0) {
            unique.prey <- NA
        }
        num.pred <- length(unique.pred)

        num.prey.pred <- c(trial, num.prey, num.pred)
        show(num.prey.pred)

        if (num.prey > 0) {
            preylist <- rbind(preylist, matrix(c(array(trial, dim = num.prey),
                                                 V[unique.prey, ]),nrow=num.prey))
        } else{
            preylist <- rbind(preylist, matrix(c(trial, array(NA, dim = q)),nrow=1))
        }
        predlist <- rbind(predlist, matrix(c(array(trial, dim = num.pred),
                                             U[unique.pred, ]),nrow=num.pred))

        num.preylist <- rbind(num.preylist,num.prey)
        num.predlist <- rbind(num.predlist,num.prey)
    }
    Sys.time() - t0; rm(t0)
    # end - trial loop ----

    preylist <- matrix(preylist, ncol=4)
    predlist <- matrix(predlist, ncol=4)

    output.prey <- data.frame(prey.pred = "prey", delta=delta, etaN = etaN, etaP = etaP,
                              rep=preylist[,1], tr1=preylist[,2], tr2=preylist[,3],
                              tr3=preylist[,4], Fitness.N=Fitness.N, Fitness.P=Fitness.P,
                              Selection.V=Selection.V,Selection.U=Selection.U)
    output.pred <- data.frame(prey.pred = "pred", delta=delta, etaN = etaN, etaP = etaP,
                              rep=predlist[,1], tr1=predlist[,2], tr2=predlist[,3],
                              tr3=predlist[,4], Fitness.N=Fitness.N, Fitness.P=Fitness.P,
                              Selection.V=Selection.V,Selection.U=Selection.U)
    output <- rbind(output.prey, output.pred)
    # if(delta.step ==1) {
    #     write.table(file="output coev instant etaN3_etaP2_22Oct18.csv", row.names = F,
    #                 output, sep=',', col.names = T)
    # }else{
    #     write.table(file="output coev instant etaN3_etaP2_22Oct18.csv", row.names = F,
    #                 output, sep=',', col.names = F, append=T)
    # }

    preylist[is.na(preylist)] <- 1

    # find unique communities w/ precision <- 10^-2
    unique.prey.coms <- 1:ntrials
    for (i in 1:(ntrials - 1)) for (j in (i + 1):ntrials) {
        if (sc(matrix(preylist[preylist[, 1] == i, 2:4], ncol=3),
               matrix(preylist[preylist[, 1] == j, 2:4], ncol=3),
               precision)) {
            unique.prey.coms[j] <- 0
        }
    }
    unique.prey.coms <- unique.prey.coms[unique.prey.coms > 0]

    unique.pred.coms <- 1:ntrials
    for (i in 1:(ntrials - 1)) for (j in (i + 1):ntrials) {
        if (sc(matrix(predlist[predlist[, 1] == i, 2:4], ncol=3),
               matrix(predlist[predlist[, 1] == j, 2:4], ncol=3),
               precision)) {
            unique.pred.coms[j] <- 0
        }
    }
    unique.pred.coms <- unique.pred.coms[unique.pred.coms > 0]

    unique.prey.coms
    unique.pred.coms

    num.unique.prey.com <- length(unique.prey.coms)
    num.unique.pred.com <- length(unique.pred.coms)

    unique.prey.coms <- union(unique.prey.coms, unique.pred.coms)
    unique.pred.coms <- unique.prey.coms

    preylist <- preylist[is.element(preylist[, 1], unique.prey.coms), ]
    predlist <- predlist[is.element(predlist[, 1], unique.pred.coms), ]

    preylist <- matrix(preylist, ncol=4)
    predlist <- matrix(predlist, ncol=4)

    num.unique.comms <- length(unique.prey.coms)
    mean.prey <- dim(preylist)[1]/length(unique.prey.coms)
    min.prey <- min(num.preylist)
    max.prey <- max(num.preylist)

    mean.pred <- dim(predlist)[1]/length(unique.pred.coms)
    min.pred <- min(num.predlist)
    max.pred <- max(num.predlist)

    mean.prey.pred <- c(mean.prey, mean.pred)

    output.prey <- data.frame(prey.pred = "prey",  delta=delta, etaN = etaN, etaP = etaP,
                              rep=preylist[,1], tr1=preylist[,2], tr2=preylist[,3],
                              tr3=preylist[,4],Fitness.N=Fitness.N, Fitness.P=Fitness.P,
                              Selection.V=Selection.V,Selection.U=Selection.U)
    output.pred <- data.frame(prey.pred = "pred",  delta =delta, etaN = etaN, etaP = etaP,
                              rep=predlist[,1], tr1=predlist[,2], tr2=predlist[,3],
                              tr3=predlist[,4],Fitness.N=Fitness.N, Fitness.P=Fitness.P,
                              Selection.V=Selection.V,Selection.U=Selection.U)
    output <- rbind(output.prey, output.pred)
    # if (delta.step == 1) {
    #     write.table(file="output unique coev instant  etaN3_etaP2_22Oct18_multdelta.csv",
    #                 row.names = F, output, sep=',', col.names = T)
    # } else{
    #     write.table(file="output unique coev instant  etaN3_etaP2_22Oct18_multdelta.csv",
    #                 row.names = F, output, sep=',', col.names = F, append=T)
    # }

}
# end - delta loop ----









# ######################################################################*
# ######################################################################*
# # data analysis
# ######################################################################*
# ######################################################################*
#
# d <- read.table(file="output unique coev instant  etaN3_etaP2_22Oct18_multdelta.csv",
#                 sep=',', header = T)
# summary(d)
#
# #etalist <- c(-.05, 0, .05, .1, .15, .2)
# #etalist <- c(.05, .1)
#
# com.no <- matrix(0, nrow=length(unique(d$etaP)), ncol=length(unique(d$etaN)))
# prey.mean.no <- matrix(0, nrow=length(unique(d$etaP)), ncol=length(unique(d$etaN)))
# pred.mean.no <- matrix(0, nrow=length(unique(d$etaP)), ncol=length(unique(d$etaN)))
# prey.max.no <- matrix(0, nrow=length(unique(d$etaP)), ncol=length(unique(d$etaN)))
# pred.max.no <- matrix(0, nrow=length(unique(d$etaP)), ncol=length(unique(d$etaN)))
#
# rownames(com.no) <- etalist
# rownames(prey.mean.no) <- etalist
# rownames(pred.mean.no) <- etalist
# rownames(prey.max.no) <- etalist
# rownames(pred.max.no) <- etalist
# colnames(com.no) <- etalist
# colnames(prey.mean.no) <- etalist
# colnames(pred.mean.no) <- etalist
# colnames(prey.max.no) <- etalist
# colnames(pred.max.no) <- etalist
#
# for(i in 1:length(etalist)) for(j in 1:length(etalist)){
#     dd <- d[d$etaN==etalist[i] & d$etaP==etalist[j],]
#     com.no[i, j] <- length(unique(dd$rep[dd$prey.pred == "prey"]))
#
#     prey.mean.no[i, j] <- length(dd$rep[dd$prey.pred == "prey"])/com.no[i, j]
#     pred.mean.no[i, j] <- length(dd$rep[dd$prey.pred == "pred"])/com.no[i, j]
#
#     prey.max <- 0
#     for(k in unique(dd$rep[dd$prey.pred == "prey"]))
#         prey.max <- max(prey.max, length(dd$rep[dd$prey.pred == "prey" & dd$rep == k]))
#     prey.max.no[i, j] <- prey.max
#
#     pred.max <- 0
#     for(k in unique(dd$rep[dd$prey.pred == "pred"]))
#         pred.max <- max(pred.max, length(dd$rep[dd$prey.pred == "pred" & dd$rep == k]))
#     pred.max.no[i, j] <- pred.max
# }
# com.no
# prey.mean.no
# pred.mean.no
# prey.max.no
# pred.max.no
#
# # d[d$etaN == .05 & d$etaP == 0,]
#
# # pdf("eta figs 22Jan16.pdf", width=4, height=9)
#
# # par(mfrow=c(2,1))
#
# # x <- prey.max.no
# # image(matrix(x, nrow=length(etalist), ncol=length(etalist)), ylab=expression(eta[P]),
# #       xlab=expression(eta[N]), main = "prey", col=gray(level=(1:max(x))/max(x)),
# #       xaxt="n", yaxt="n", cex.lab=1.4)
# # axis(side=1,at=c(0,.2,.4,.6,.8,1), labels=etalist)
# # axis(side=2,at=c(0,.2,.4,.6,.8,1), labels=etalist)
#
# # x <- pred.max.no
# # image(matrix(x, nrow=length(etalist), ncol=length(etalist)), ylab=expression(eta[P]),
# #       xlab=expression(eta[N]), main = "pred", col=gray(level=(1:max(x))/max(x)),
# #       xaxt="n", yaxt="n", cex.lab=1.4)
# # axis(side=1,at=c(0,.2,.4,.6,.8,1), labels=etalist)
# # axis(side=2,at=c(0,.2,.4,.6,.8,1), labels=etalist)
#
# # dev.off()
