# setwd('/Users/jc269188/Documents/EcoEvo/Community coevolution/2016 code/2 traits/2018 code')

# -----------------
# Quantitative genetics version
# -----------------

# sc <- function(x1, x2, precision = 10^-2) {
# 	flag <- F
# 	if (dim(x1)[1] == dim(x2)[1]) {
# 		n <- dim(x1)[1]
# 		flag <- T
# 		i <- 0
# 		while (i < n & flag == T) {
# 			flag <- F
# 			i <- i + 1
# 			for (j in 1:n) {
# 				if (mean((x1[i, ] - x2[j, ])^2) < precision^2)
# 					flag <- T
# 			}
#
# 		}
# 	}
# 	return(flag)
# }
#
#
# collist <- c("black", "blue", "dark green", "red", "green", "dark orange", "gray",
#              "deepskyblue2", "firebrick4", "darkorchid")
# pchlist <- c(19, 22, 24, 23, 4, 8)
#
# traj.flag <- T


q <- 2  # number of traits
# precision <- 0.02
# Nthreshold <- 1e-05
# Pthreshold <- 1e-05

# Ninit <- 0.001
# Pinit <- 0.001

# delta.flag <- F
# delta <- 0.001
nreps <- 1

# turns evolution "on" with TRUE\n
evoN <- T
evoP <- T

# Non-additive tradeoffs for increasing >1 traits
etaN <- 0.1
etaP <- 0.1


# at timeDis values of N, P, V, and U are all perturbed\n
timeAdd.long <- 120000
timeAdd.long <- timeAdd.long / 1000  # <-- added by LAN (Nov 2018)

# additive genetic variances:
sig2N <- 10^-10
sig2P <- 10^-10

# sdN <- 0
# sdP <- 0
# sdV <- 0
# sdU <- 0

if (evoN) {
	#	sig2N <- 1
	sig2N <- 0.05
	# sdN <- 0.1
	# sdV <- 0.1
}
if (evoP) {
	#	sig2P <- 1
	sig2P <- 0.05
	# sdP <- 0.1
	# sdU <- 0.1
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


###################################################
ntrial <- 4
ntrial <- 1  # <-- added by LAN (Nov 2018)


# These keep track of which species are extant
Nsp.ident.vec <- rep(0, 10)
Nsp.ident.vec[1] <- 1
Psp.ident.vec <- rep(0, 10)
Psp.ident.vec[1] <- 1

C <- matrix(etaN, nrow = q, ncol = q)
diag(C) <- 1
D <- matrix(etaP, nrow = q, ncol = q)
diag(D) <- 1

commSizelist <- data.frame(rep = array(0, dim = nreps), prey = 0, pred = 0)

graphstep <- 1

timevec <- c(1:(timeAdd.long * ntrial))
timecount <- 1

commV <- list(0)
commU <- list(0)

commSize <- matrix(0, nrow = ntrial, ncol = 2)

addcounter <- 0

Nlist.full <- array(0, c(20, length(timevec)))
Plist.full <- array(0, c(20, length(timevec)))
Vlist.full <- array(0, c(20, q, length(timevec)))
Ulist.full <- array(0, c(20, q, length(timevec)))

n.list <- c(1,2,3,4)
p.list <- c(1,2,2,2)

startV <- list(matrix(c(0.4840364, 0.4840364), ncol=q, byrow = T),
	matrix(c(0.5093952, 0.472734, 1.446557e-05, 0.6329789), ncol=q, byrow = T),
	matrix(c(1.457756e-13, 0.6443325, 0.6427733, 0.0522391, 0.5194171, 0.4975425), ncol = 2, byrow = T),
	matrix(c(9.274664e-09, 0.7338789, 9.274664e-09, 0.7338789, 0.7382802, 1.216864e-08, 0.7382802, 1.640028e-05), ncol=q, byrow = T))

startU <- list(matrix(c(0.8851652, 0.8851652), ncol=q, byrow = T),
	matrix(c(1.149305, 0.541385, 0.5761137, 1.1261267), ncol=q, byrow = T),
	matrix(c(0.0524799, 1.6050365, 1.60053494, 0.05322435), ncol=q, byrow = T),
	matrix(c(0.4458685, 1.2101178, 1.2112638, 0.4426656), ncol = 2, byrow = T))

output <- data.frame(x = 0, y = 0, fitness = 0, trial = 0, time = 0, sp = 0)
write.table(output, file = "V.csv", sep = ",", append = F, col.names = T, row.names = F)
write.table(output, file = "U.csv", sep = ",", append = F, col.names = T, row.names = F)

for (trial in 1:ntrial) {

	n <- n.list[trial]
	p <- p.list[trial]
	extant.N <- 1:n
	extant.P <- 1:p

	initN <- 0.1
	initP <- 0.1

	N <- matrix(initN, nrow = n, ncol = 1)
	P <- matrix(initP, nrow = p, ncol = 1)

	V <- startV[[trial]]
	U <- startU[[trial]]

	timeAdd <- timeAdd.long

	Nlist <- array(0, c(length(Nsp.ident.vec), timeAdd))
	Plist <- array(0, c(length(Psp.ident.vec), timeAdd))
	Vlist <- array(0, c(length(Nsp.ident.vec), q, timeAdd))
	Ulist <- array(0, c(length(Psp.ident.vec), q, timeAdd))

	for (time in 1:timeAdd) {

		A <- a + f * diag(V %*% C %*% t(V))
		B <- b * exp(-(V^2) %*% t(U^2))
		M <- m + g/diag(U %*% D %*% t(U))

		Nt <- N * exp(r * (1 - A * sum(N) - B %*% P))
		Pt <- P * exp(cc * t(B) %*% N - M)
		Vt <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V * exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
		Ut <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*% (V^2 * array(N, c(n, q))) + 2 * g * array(1/diag(U %*% D %*% t(U))^2,
			c(p, q)) * (U %*% D))

		N <- Nt
		P <- Pt
		V <- abs(Vt)
		U <- abs(Ut)

		#########################
		# catch low values of U
		U[U < 10^-4] <- 10^-4
		#########################

		Nlist[extant.N, time] <- N
		Plist[extant.P, time] <- P
		Vlist[extant.N, , time] <- V
		Ulist[extant.P, , time] <- U

		if (timecount%%graphstep == 0) {
			Nlist.full[extant.N, (timecount)] <- N
			Plist.full[extant.P, (timecount)] <- P
			Vlist.full[extant.N, , (timecount)] <- V
			Ulist.full[extant.P, , (timecount)] <- U
		}

		#		if (time == timeAdd.long | time == 1) {
		if (time == timeAdd.long) {

			###### Jacobian
			eps <- -10^-6
			nn <- (n + p) * q
			Jacobian <- matrix(NA, nrow = nn, ncol = nn)
			VU <- matrix(c(t(V), t(U)), ncol = nn)
			FV <- V + sig2N * (-2 * r * f * sum(N) * V %*% C + 2 * r * b * V * exp(-(V^2) %*% t(U^2)) %*% (U^2 * array(P, c(p, q))))
			FU <- U + sig2P * (-2 * b * cc * U * exp(-(U^2) %*% t(V^2)) %*% (V^2 * array(N, c(n, q))) + 2 * g * array(1/diag(U %*% D %*% t(U))^2,
				c(p, q)) * (U %*% D))
			FVU <- matrix(c(t(FV), t(FU)), ncol = nn)
			for (i in 1:nn) {
				VU1 <- VU
				VU1[i] <- VU[i] + eps
				V1 <- matrix(VU1[1:(n * q)], ncol = q, byrow = T)
				U1 <- matrix(VU1[(n * q + 1):nn], ncol = q, byrow = T)
				FV1 <- V1 + sig2N * (-2 * r * f * sum(N) * V1 %*% C + 2 * r * b * V1 * exp(-(V1^2) %*% t(U1^2)) %*% (U1^2 * array(P, c(p, q))))
				FU1 <- U1 + sig2P * (-2 * b * cc * U1 * exp(-(U1^2) %*% t(V1^2)) %*% (V1^2 * array(N, c(n, q))) + 2 * g * array(1/diag(U1 %*%
					D %*% t(U1))^2, c(p, q)) * (U1 %*% D))
				FVU1 <- matrix(c(t(FV1), t(FU1)), ncol = nn)
				Jacobian[i, ] <- (FVU1 - FVU)/eps
			}
			HessianN <- array(0, dim = c(q, q, n))
			for (i in 1:n) HessianN[, , i] <- Jacobian[(q * (i - 1) + 1):(q * i), (q * (i - 1) + 1):(q * i)] - diag(1, q, q)
			HessianP <- array(0, dim = c(q, q, p))
			for (i in 1:p) HessianP[, , i] <- Jacobian[(n * q + q * (i - 1) + 1):(n * q + q * i), (n * q + q * (i - 1) + 1):(n * q + q * i)] -
				diag(1, q, q)

			print("")
			print("")
			print(paste("NEW TRIAL, Trial number", trial, " Last Time = ", time != 1))

			print(Jacobian)
			print(paste("Jacobian eigenvalue = ", max(abs(eigen(Jacobian)$value))))

			for (sp in 1:n) {
				print(paste("Resource species", sp))
				print("Density")
				print(N[sp])
				print("Traits")
				print(V[sp, ])
				print("Hessian")
				print(HessianN[, , sp])
				print(paste("eigenvalues", eigen(HessianN[, , sp])$value))
			}
			print(paste("V Jacobian eigenvalue = ", max(abs(eigen(Jacobian[1:(n * q), 1:(n * q)])$value))))

			for (sp in 1:p) {
				print(paste("Consumer species", sp))
				print("Density")
				print(P[sp])
				print("Traits")
				print(U[sp, ])
				print("Hessian")
				print(HessianP[, , sp])
				print(paste("eigenvalues", eigen(HessianP[, , sp])$value))
				print(paste("det", det(HessianP[, , sp])))
			}
			print(paste("U Jacobian eigenvalue = ", max(abs(eigen(Jacobian[(1 + n * q):nn, (1 + n * q):nn])$value))))
		}
		################################################################
		# save files to plot fitness
		################################################################
		if (time == timeAdd.long) {
			# this code plots fitness only for trait 1
			for(spN in 1:n){
				output$trial[1] <- trial
				output$time[1] <- time
				output$sp[1] <- spN
				VV <- V
				for (ii in 0:100) for(jj in 0:100){
					VV[spN, 1] <- 1 * ii/100
					VV[spN, 2] <- 1 * jj/100
					AA <- a + f * diag(VV %*% C %*% t(VV))
					BB <- b * exp(-(VV^2) %*% t(U^2))

					output$fitness[1] <- exp(r * (1 - AA * sum(N) - BB %*% P))[spN]
					output$x[1] <- 1 * ii/100
					output$y[1] <- 1 * jj/100
					write.table(output, file = "V.csv", sep = ",", append = T, col.names = F, row.names = F)
				}
			}
			for(spP in 1:p){
				output$trial[1] <- trial
				output$time[1] <- time
				output$sp[1] <- spP
				UU <- U
				for (ii in 0:100) for(jj in 0:100){
					UU[spP, 1] <- 2 * ii/100
					UU[spP, 2] <- 2 * jj/100
					MM <- m + g/diag(UU %*% D %*% t(UU))
					BB <- b * exp(-(V^2) %*% t(UU^2))

					output$fitness[1] <- exp(cc * t(BB) %*% N - MM)[spP]
					output$x[1] <- 2 * ii/100
					output$y[1] <- 2 * jj/100
					write.table(output, file = "U.csv", sep = ",", append = T, col.names = F, row.names = F)
				}
			}
		}

		timecount <- timecount + 1
	}
}

timevec.short <- timevec/(10^3)



############################
# plotting figure 2
############################
ntrial <- 4
ntrial <- 1  # added by LAN (Nov 2018)


colvec <- c("black", "blue", "darkgreen", "red")
pchvec <- 15:18
collist <- colvec
cex <- 1.5
cex.lab <- 2.5
# Next two lines changed by LAN (Nov 2018)
# png("Northfield et al. Fig 2 10Oct18.png", height = 1200, width = 1100)
source(".Rprofile")
layout(matrix(c(1, 1, 3, 4, 1, 1, 6, 5, 2, 2, 7, 8, 2, 2, 10, 9), ncol = 4))
par(mai = c(0.8, 0.8, 0.1, 0.1), oma = c(5, 0, 0, 0))

start.vec <- (0:(ntrial - 1)) * timeAdd.long

# resources
plot(Vlist.full[1, 1, ], Vlist.full[1, 2, ], type = "n", ylab = "Resource trait 2", xlim = c(0, 1), ylim = c(0, 1), xlab = "Resource trait 1",
	cex = cex, cex.lab = cex.lab)
text(x = 0.005, y = 0.95, labels = "A. Resources", xpd = TRUE, adj = c(0, 0), cex = cex.lab)
for (t in 1:ntrial) for (i in 1:n.list[t]) lines(Vlist.full[i, 1, (start.vec[t] + 1):(start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] +
	1):(start.vec[t] + timeAdd.long)], col = colvec[t], lty = t)
for (t in 1:ntrial) for (i in 1:n.list[t]) {
	points(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t],
		cex = 5)
	show(c(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)]))
}
# consumersstart.vec <- c(1, 2 * timeAdd.long + 1)
plot(Ulist.full[1, 1, ], Ulist.full[1, 2, ], type = "n", ylab = "Consumer trait 2", xlim = c(0, 2), ylim = c(0, 2), xlab = "Consumer trait 1",
	cex = cex, cex.lab = cex.lab)
text(x = 0.005, y = 2 * 0.95, labels = "B. Consumers", xpd = TRUE, adj = c(0, 0), cex = cex.lab)
for (t in 1:ntrial) for (i in 1:p.list[t]) lines(Ulist.full[i, 1, (start.vec[t] + 1):(start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] +
	1):(start.vec[t] + timeAdd.long)], col = colvec[t], lty = t)
for (t in 1:ntrial) for (i in 1:p.list[t]) {
	points(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t],
		cex = 5)
	show(c(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)]))
}

##################
# fitness panels
par(mai = c(0.2, 0.8, .7, 0.1))
fitnessV <- read.csv(file = "V.csv", header = T)
for (trial in 1:ntrial) {
	ff <- fitnessV[fitnessV$trial == trial & fitnessV$sp == 1, ]
	x <- unique(ff$x)
	y <- unique(ff$y)
	z <- matrix(ff$fitness, nrow=length(x), byrow=T)^20
	if (trial == 1) {
		minf <- -.01
		contour(x, y, z, nlevels=200, drawlabels=F)
		t <- trial
 		for (i in 1:n.list[trial]) 	points(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("1 saddle node", side=3, cex = 1.5, col = collist[trial])
	}
	if (trial == 2) {
		minf <- -0.03
		contour(x, y, z, nlevels=200, drawlabels=F)
		t <- trial
 		for (i in 1:n.list[trial]) 	points(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("2 stable points", side=3, cex = 1.5, col = collist[trial])
		mtext("Resource trait 2", side = 2, cex = cex, padj = -2.5, adj=1.7)
	}
	if (trial == 3) {
		minf <- -0.2
		contour(x, y, z, nlevels=200, drawlabels=F)
		t <- trial
 		for (i in 1:n.list[trial]) 	points(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("3 stable points", side=3, cex = 1.5, col = collist[trial])
		mtext("Resource trait 1", side = 1, cex = cex, padj = 2.5, adj=-1.4)
	}
	if (trial == 4) {
		minf <- -0.03
		contour(x, y, z, nlevels=200, drawlabels=F)
		t <- trial
 		for (i in 1:n.list[trial]) 	points(Vlist.full[i, 1, (start.vec[t] + timeAdd.long)], Vlist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("2 stable points", side=3, cex = 1.5, col = collist[trial])
	}
}

fitnessU <- read.csv(file = "U.csv", header = T)
for (trial in 1:ntrial) {
	ff <- fitnessU[fitnessU$trial == trial & fitnessU$sp == 1, ]
	x <- unique(ff$x)
	y <- unique(ff$y)
	z <- matrix(ff$fitness, nrow=length(x), byrow=T)^20
	if (trial == 1) {
		minf <- -.01
		contour(x, y, z, nlevels=500, drawlabels=F)
		t <- trial
 		for (i in 1:p.list[trial]) 	points(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("1 stable point", side=3, cex = 1.5, col = collist[trial])
	}
	if (trial == 2) {
		minf <- -0.03
		contour(x, y, z, nlevels=500, drawlabels=F)
		t <- trial
 		for (i in 1:p.list[trial]) 	points(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("2 stable points", side=3, cex = 1.5, col = collist[trial])
		mtext("Consumer trait 2", side = 2, cex = cex, padj = -2.5, adj=1.7)
	}
	if (trial == 3) {
		minf <- -0.2
		contour(x, y, z, nlevels=500, drawlabels=F)
		t <- trial
 		for (i in 1:p.list[trial]) 	points(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("2 stable points", side=3, cex = 1.5, col = collist[trial])
		mtext("Consumer trait 1", side = 1, cex = cex, padj = 2.5, adj=-1.4)
	}
	if (trial == 4) {
		minf <- -0.03
		contour(x, y, z, nlevels=500, drawlabels=F)
		t <- trial
 		for (i in 1:p.list[trial]) 	points(Ulist.full[i, 1, (start.vec[t] + timeAdd.long)], Ulist.full[i, 2, (start.vec[t] + timeAdd.long)], col = colvec[t], pch = pchvec[t], cex = 3)
		mtext("2 stable points", side=3, cex = 1.5, col = collist[trial])
	}
}
dev.off()



# [1] "NEW TRIAL, Trial number 1  Last Time =  TRUE"
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]          [,9]         [,10]
# [1,]  9.999307e-01 -1.149341e-04  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -3.684880e-04  8.285328e-05 -3.684880e-04  8.285328e-05
# [2,] -1.149341e-04  9.999307e-01  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  8.285328e-05 -3.684880e-04  8.285328e-05 -3.684880e-04
# [3,]  0.000000e+00  0.000000e+00  9.999307e-01 -1.149341e-04  0.000000e+00  0.000000e+00 -3.684880e-04  8.285328e-05 -3.684880e-04  8.285328e-05
# [4,]  0.000000e+00  0.000000e+00 -1.149341e-04  9.999307e-01  0.000000e+00  0.000000e+00  8.285328e-05 -3.684880e-04  8.285328e-05 -3.684880e-04
# [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  9.999307e-01 -1.149341e-04 -3.684880e-04  8.285328e-05 -3.684880e-04  8.285328e-05
# [6,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.149341e-04  9.999307e-01  8.285328e-05 -3.684880e-04  8.285328e-05 -3.684880e-04
# [7,]  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  9.993618e-01 -5.708481e-04  0.000000e+00  0.000000e+00
# [8,] -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04 -5.708481e-04  9.993618e-01  0.000000e+00  0.000000e+00
# [9,]  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  0.000000e+00  0.000000e+00  9.993618e-01 -5.708481e-04
# [10,] -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04 -2.518791e-05  1.120229e-04  0.000000e+00  0.000000e+00 -5.708481e-04  9.993618e-01
# [1] "Jacobian eigenvalue =  1.00004562192391"
# [1] "Resource species 1"
# [1] "Density"
# [1] 0.1520725
# [1] "Traits"
# [1] 0.4840364 0.4840364
# [1] "Hessian"
# [,1]          [,2]
# [1,] -6.931214e-05 -1.149341e-04
# [2,] -1.149341e-04 -6.931214e-05
# [1] "eigenvalues 4.5621923906225e-05"   "eigenvalues -0.000184246200628735"
# [1] "Resource species 2"
# [1] "Density"
# [1] 0.1520725
# [1] "Traits"
# [1] 0.4840364 0.4840364
# [1] "Hessian"
# [,1]          [,2]
# [1,] -6.931214e-05 -1.149341e-04
# [2,] -1.149341e-04 -6.931214e-05
# [1] "eigenvalues 4.5621923906225e-05"   "eigenvalues -0.000184246200628735"
# [1] "Resource species 3"
# [1] "Density"
# [1] 0.1520725
# [1] "Traits"
# [1] 0.4840364 0.4840364
# [1] "Hessian"
# [,1]          [,2]
# [1,] -6.931214e-05 -1.149341e-04
# [2,] -1.149341e-04 -6.931214e-05
# [1] "eigenvalues 4.5621923906225e-05"   "eigenvalues -0.000184246200628735"
# [1] "V Jacobian eigenvalue =  1.00004562192391"
# [1] "Consumer species 1"
# [1] "Density"
# [1] 0.9246216
# [1] "Traits"
# [1] 0.8851652 0.8851652
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0006381588 -0.0005708481
# [2,] -0.0005708481 -0.0006381588
# [1] "eigenvalues -6.73106838036119e-05" "eigenvalues -0.00120900697686466"
# [1] "det 8.13790863360977e-08"
# [1] "Consumer species 2"
# [1] "Density"
# [1] 0.9246216
# [1] "Traits"
# [1] 0.8851652 0.8851652
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0006381588 -0.0005708481
# [2,] -0.0005708481 -0.0006381588
# [1] "eigenvalues -6.73106838036119e-05" "eigenvalues -0.00120900697686466"
# [1] "det 8.13790863360977e-08"
# [1] "U Jacobian eigenvalue =  0.999932689316196"
# [1] ""
# [1] ""
# [1] "NEW TRIAL, Trial number 2  Last Time =  TRUE"
# [,1]          [,2]         [,3]          [,4]          [,5]          [,6]          [,7]          [,8]          [,9]         [,10]
# [1,]  9.998471e-01 -8.175288e-05  0.00000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -3.895555e-04  6.825918e-05 -3.672774e-04  4.644973e-05
# [2,] -8.175294e-05  9.999357e-01  0.00000e+00  0.000000e+00  0.000000e+00  0.000000e+00  4.040057e-05 -2.067430e-04  1.166853e-04 -4.397338e-04
# [3,]  0.000000e+00  0.000000e+00  1.00008e+00 -2.283751e-05  0.000000e+00  0.000000e+00  4.440892e-10 -1.110223e-10  1.110223e-10  0.000000e+00
# [4,]  0.000000e+00  0.000000e+00 -2.28376e-05  9.997917e-01  0.000000e+00  0.000000e+00  0.000000e+00 -2.958850e-04  0.000000e+00 -2.157952e-04
# [5,]  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00  9.998471e-01 -8.175288e-05 -4.109368e-04  7.200562e-05 -3.874359e-04  4.899925e-05
# [6,]  0.000000e+00  0.000000e+00  0.00000e+00  0.000000e+00 -8.175294e-05  9.999357e-01  4.261813e-05 -2.180903e-04  1.230898e-04 -4.638692e-04
# [7,]  9.548973e-05 -9.903245e-06  0.00000e+00  0.000000e+00  9.548973e-05 -9.903245e-06  9.991737e-01 -4.632605e-04  0.000000e+00  0.000000e+00
# [8,] -1.673206e-05  5.067791e-05  0.00000e+00  1.100963e-04 -1.673206e-05  5.067791e-05 -4.632601e-04  9.997041e-01  0.000000e+00  0.000000e+00
# [9,]  1.108043e-04 -3.520295e-05  0.00000e+00  0.000000e+00  1.108043e-04 -3.520295e-05  0.000000e+00  0.000000e+00  9.996343e-01 -5.238177e-04
# [10,] -1.401346e-05  1.326639e-04  0.00000e+00  9.882495e-05 -1.401346e-05  1.326638e-04  0.000000e+00  0.000000e+00 -5.238184e-04  9.991704e-01
# [1] "Jacobian eigenvalue =  1.00008178816169"
# [1] "Resource species 1"
# [1] "Density"
# [1] 0.1683157
# [1] "Traits"
# [1] 0.5601175 0.4077974
# [1] "Hessian"
# [,1]          [,2]
# [1,] -1.528625e-04 -8.175288e-05
# [2,] -8.175294e-05 -6.427345e-05
# [1] "eigenvalues -0.000201549366624508" "eigenvalues -1.55865517304193e-05"
# [1] "Resource species 2"
# [1] "Density"
# [1] 0.1108824
# [1] "Traits"
# [1] 2.371692e-20 6.657516e-01
# [1] "Hessian"
# [,1]          [,2]
# [1,]  8.04784e-05 -2.283751e-05
# [2,] -2.28376e-05 -2.082709e-04
# [1] "eigenvalues -0.00021006601948875" "eigenvalues 8.22734895546378e-05"
# [1] "Resource species 3"
# [1] "Density"
# [1] 0.177554
# [1] "Traits"
# [1] 0.5601175 0.4077974
# [1] "Hessian"
# [,1]          [,2]
# [1,] -1.528625e-04 -8.175288e-05
# [2,] -8.175294e-05 -6.427345e-05
# [1] "eigenvalues -0.000201549366624508" "eigenvalues -1.55865517304193e-05"
# [1] "V Jacobian eigenvalue =  1.00008227348955"
# [1] "Consumer species 1"
# [1] "Density"
# [1] 0.8251679
# [1] "Traits"
# [1] 1.1763063 0.5068916
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0008263306 -0.0004632605
# [2,] -0.0004632601 -0.0002959472
# [1] "eigenvalues -0.00109493362130374" "eigenvalues -2.7344258433767e-05"
# [1] "det 2.99401479087499e-08"
# [1] "Consumer species 2"
# [1] "Density"
# [1] 1.015587
# [1] "Traits"
# [1] 0.606490 1.109233
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0003656904 -0.0005238177
# [2,] -0.0005238184 -0.0008296449
# [1] "eigenvalues -0.00117055401619254" "eigenvalues -2.4781312524312e-05"
# [1] "det 2.90078649018558e-08"
# [1] "U Jacobian eigenvalue =  0.999975218687476"
# [1] ""
# [1] ""
# [1] "NEW TRIAL, Trial number 3  Last Time =  TRUE"
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]          [,9]         [,10]
# [1,]  1.000667e+00 -2.223499e-05  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.387779e-11  0.000000e+00  1.776357e-09 -7.632783e-11
# [2,] -2.223498e-05  9.993835e-01  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  2.557767e-04  0.000000e+00 -8.628951e-05
# [3,]  0.000000e+00  0.000000e+00  9.993835e-01 -2.223498e-05  0.000000e+00  0.000000e+00 -4.314486e-05  0.000000e+00  1.278888e-04  0.000000e+00
# [4,]  0.000000e+00  0.000000e+00 -2.223499e-05  1.000667e+00  0.000000e+00  0.000000e+00 -4.163336e-11  8.881784e-10  0.000000e+00  6.938894e-12
# [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  9.993835e-01 -2.223506e-05 -4.314486e-05  0.000000e+00  1.278888e-04  0.000000e+00
# [6,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -2.223510e-05  1.000667e+00  1.297573e-09 -2.819966e-08  0.000000e+00 -2.289835e-10
# [7,]  1.016440e-13  0.000000e+00  2.691869e-05 -5.556536e-13  2.691869e-05 -8.346154e-10  9.995389e-01 -6.412604e-05  0.000000e+00  0.000000e+00
# [8,] -3.388132e-15 -7.979239e-05  0.000000e+00  1.211257e-11  0.000000e+00  1.817114e-08 -6.412620e-05  9.998201e-01  0.000000e+00  0.000000e+00
# [9,]  1.211596e-11  0.000000e+00 -7.979251e-05 -3.388132e-15 -7.979251e-05 -6.854191e-12  0.000000e+00  0.000000e+00  9.998201e-01 -6.412619e-05
# [10,] -5.556536e-13  2.691869e-05  0.000000e+00  9.825582e-14  0.000000e+00  1.491862e-10  0.000000e+00  0.000000e+00 -6.412582e-05  9.995389e-01
# [1] "Jacobian eigenvalue =  1.00066779963391"
# [1] "Resource species 1"
# [1] "Density"
# [1] 0.2223496
# [1] "Traits"
# [1] 1.091100e-08 7.362297e-01
# [1] "Hessian"
# [,1]          [,2]
# [1,]  6.674159e-04 -2.223499e-05
# [2,] -2.223498e-05 -6.164681e-04
# [1] "eigenvalues 0.000667800817017003" "eigenvalues -0.00061685303202019"
# [1] "Resource species 2"
# [1] "Density"
# [1] 0.111175
# [1] "Traits"
# [1] 7.362298e-01 1.091098e-08
# [1] "Hessian"
# [,1]          [,2]
# [1,] -6.164686e-04 -2.223498e-05
# [2,] -2.223499e-05  6.674147e-04
# [1] "eigenvalues 0.000667799633145129"  "eigenvalues -0.000616853587320233"
# [1] "Resource species 3"
# [1] "Density"
# [1] 0.111175
# [1] "Traits"
# [1] 7.362298e-01 1.637004e-05
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0006164686 -2.223506e-05
# [2,] -0.0000222351  6.674147e-04
# [1] "eigenvalues 0.000667799634633999"  "eigenvalues -0.000616853590537055"
# [1] "V Jacobian eigenvalue =  1.00066780081702"
# [1] "Consumer species 1"
# [1] "Density"
# [1] 1.387285
# [1] "Traits"
# [1] 0.05287194 1.60270985
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0004610942 -6.412604e-05
# [2,] -0.0000641262 -1.798853e-04
# [1] "eigenvalues -0.000475027024796279" "eigenvalues -0.000165952471848599"
# [1] "det 7.8831908959828e-08"
# [1] "Consumer species 2"
# [1] "Density"
# [1] 1.387287
# [1] "Traits"
# [1] 1.60270988 0.05287193
# [1] "Hessian"
# [,1]          [,2]
# [1,] -1.798846e-04 -6.412619e-05
# [2,] -6.412582e-05 -4.610928e-04
# [1] "eigenvalues -0.000475025595535203" "eigenvalues -0.000165951819441504"
# [1] "det 7.88313618603508e-08"
# [1] "U Jacobian eigenvalue =  0.999834048180559"
# [1] ""
# [1] ""
# [1] "NEW TRIAL, Trial number 4  Last Time =  TRUE"
# [,1]          [,2]          [,3]          [,4]          [,5]          [,6]          [,7]          [,8]          [,9]         [,10]
# [1,]  1.000088e+00 -2.284029e-05  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  1.110223e-10  0.000000e+00  4.440892e-10 -1.110223e-10
# [2,] -2.284039e-05  9.997803e-01  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -1.426865e-04  0.000000e+00 -2.189975e-04
# [3,]  0.000000e+00  0.000000e+00  9.997803e-01 -2.284273e-05  0.000000e+00  0.000000e+00 -2.189976e-04  0.000000e+00 -1.426865e-04  0.000000e+00
# [4,]  0.000000e+00  0.000000e+00 -2.284262e-05  1.000088e+00  0.000000e+00  0.000000e+00  3.219647e-09 -1.265654e-08  8.881784e-10 -3.330669e-09
# [5,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00  9.998738e-01 -7.707734e-05 -4.616365e-04  7.288325e-05 -7.495631e-04  1.679018e-04
# [6,]  0.000000e+00  0.000000e+00  0.000000e+00  0.000000e+00 -7.707740e-05  9.998738e-01  1.679017e-04 -7.495631e-04  7.288325e-05 -4.616364e-04
# [7,]  0.000000e+00  0.000000e+00  1.223910e-04 -1.867481e-09  7.288636e-05 -2.650946e-05  9.997005e-01 -4.648066e-04  0.000000e+00  0.000000e+00
# [8,]  0.000000e+00  7.974310e-05  0.000000e+00  7.249992e-09 -1.150735e-05  1.183462e-04 -4.648074e-04  9.991559e-01  0.000000e+00  0.000000e+00
# [9,]  0.000000e+00  0.000000e+00  7.974299e-05 -5.008167e-10  1.183462e-04 -1.150730e-05  0.000000e+00  0.000000e+00  9.991559e-01 -4.648073e-04
# [10,]  0.000000e+00  1.223911e-04  0.000000e+00  1.944286e-09 -2.650946e-05  7.288636e-05  0.000000e+00  0.000000e+00 -4.648066e-04  9.997005e-01
# [1] "Jacobian eigenvalue =  1.00008938626761"
# [1] "Resource species 1"
# [1] "Density"
# [1] 0.08246102
# [1] "Traits"
# [1] 2.316582e-13 6.569550e-01
# [1] "Hessian"
# [,1]          [,2]
# [1,]  8.787380e-05 -2.284029e-05
# [2,] -2.284039e-05 -2.196530e-04
# [1] "eigenvalues -0.000221340164287069" "eigenvalues 8.95609233367556e-05"
# [1] "Resource species 2"
# [1] "Density"
# [1] 0.08246102
# [1] "Traits"
# [1] 6.569550e-01 1.500379e-05
# [1] "Hessian"
# [,1]          [,2]
# [1,] -2.196532e-04 -2.284273e-05
# [2,] -2.284262e-05  8.787380e-05
# [1] "eigenvalues -0.000221340617970369" "eigenvalues 8.95612657299673e-05"
# [1] "Resource species 3"
# [1] "Density"
# [1] 0.2918857
# [1] "Traits"
# [1] 0.4975606 0.4975606
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0001261889 -7.707734e-05
# [2,] -0.0000770774 -1.261888e-04
# [1] "eigenvalues -0.000203266208419706" "eigenvalues -4.91114638947465e-05"
# [1] "V Jacobian eigenvalue =  1.00008956126573"
# [1] "Consumer species 1"
# [1] "Density"
# [1] 0.9216988
# [1] "Traits"
# [1] 0.5089913 1.1725664
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0002994614 -0.0004648066
# [2,] -0.0004648074 -0.0008441266
# [1] "eigenvalues -0.00111050602610542"  "eigenvalues -3.30820294782614e-05"
# [1] "det 3.67377930914064e-08"
# [1] "Consumer species 2"
# [1] "Density"
# [1] 0.9216988
# [1] "Traits"
# [1] 1.1725664 0.5089913
# [1] "Hessian"
# [,1]          [,2]
# [1,] -0.0008441266 -0.0004648073
# [2,] -0.0004648066 -0.0002994613
# [1] "eigenvalues -0.00111050595076093"  "eigenvalues -3.30819938004423e-05"
# [1] "det 3.67377509784275e-08"
# [1] "U Jacobian eigenvalue =  0.9999669180062"



png("Northfield et al. fig 2 Trajectory 23Sep18.png", height = 1000, width = 1000)

par(mfcol = c(2, 1 + q), cex = 1.5, cex.lab = 1.5)

plot(timevec.short, Nlist.full[1, ], type = "n", ylab = "Resource density", ylim = c(0, 0.8), xlab = expression(paste("Time (10"^5, ")", sep = "")))
text(x = 0.5, y = 0.95 * 0.8, labels = "A", xpd = TRUE, adj = c(0, 0), cex = 1.2)
for (i in 1:length(Nlist.full[, 1])) lines(timevec.short, Nlist.full[i, ], col = collist[i])
plot(timevec.short, Plist.full[1, ], type = "n", ylab = "Consumer density", xlab = expression(paste("Time (10"^3, ")")), ylim = c(0, 5))
text(x = 0.5, y = 0.95 * 5, labels = "D", xpd = TRUE, adj = c(0, 0), cex = 1.2)
for (i in 1:length(Plist.full[, 1])) lines(timevec.short, Plist.full[i, ], col = collist[i])
lettervec = c("B", "C", "E", "F", "G", "H")
for (j in 1:q) {
	plot(timevec.short, Vlist.full[1, j, ], type = "n", ylab = paste("Resource trait", j), ylim = c(0, 1), xlab = expression(paste("Time (10"^3,
		")")))
	text(x = 0.5, y = 0.95, labels = paste(lettervec[j]), xpd = TRUE, adj = c(0, 0), cex = 1.2)
	for (i in 1:length(Nlist.full[, 1])) lines(timevec.short, Vlist.full[i, j, ], col = collist[i])
	plot(timevec.short, Ulist.full[1, j, ], type = "n", ylab = paste("Consumer trait", j), ylim = c(0, 2), xlab = expression(paste("Time (10"^3,
		")")))
	text(x = 0.5, y = 2 * 0.95, labels = paste(lettervec[j + 2]), xpd = TRUE, adj = c(0, 0), cex = 1.2)
	for (i in 1:length(Plist.full[, 1])) lines(timevec.short, Ulist.full[i, j, ], col = collist[i])
}
dev.off()

