require(foreach)
require(doParallel)
require(KFAS)
require(lme4)
require(vars)

# generate simulation data ------------------------------------------------
# uncomment the settings you want to generate
# n is time length (variable T in the article)
# N is number of subjects
# p is expectation of a

cl <- parallel::makeCluster(20)
registerDoParallel(cl)
iterN <- 200
for (n in c(15,30,60)){
  N_set <- if (n==30){
    c(20,40,60)
  }else{
    40
  }
  for (N in N_set){
    p_set <- if (n==30 & N==40){
      seq(0,0.9,0.15)
    }else{
      0.3
    }
    for (p in p_set){
      all_res_var <- foreach(iter = 1:iterN, .packages = c('lme4')) %dopar% {
        set.seed(iter)
        sa <- function(start, c, int, sd, n){
          y <- vector('numeric',n)
          y[1] <- start
          for (i in 2:n){
            y[i] <- y[i-1]*c + int + rnorm(1, 0, sd)
          }
          return(y)
        }
        all_dat <- vector('list',N)
        all_x <- data.frame(X1=numeric(),X2=numeric())
        for (i in 1:N){
          a1s = sa(start = 10, c = 0.6, int = 1.2, sd = 1, n = n)
          a2s = sa(start = 10, c = 0.8, int = 2, sd = sqrt(1.5), n = n)
          a <- rbinom(n, 1, p)
          lam1 <- c(-0.5,0.2,-1)
          lam2 <- c(0.1,0.2,1)
          mu1 <- t(t(lam1)) %*% a1s
          mu1t <- sweep(mu1, 2, a, '*')
          mu2 <- t(t(lam2)) %*% a2s
          
          X1 <- rnorm(1,-5,sqrt(2))
          X2 <- rbinom(1,4,0.5)
          
          mu <- mu1t+mu2+X1+2*X2+0.03*matrix(rep(1:n,3), byrow = T, nrow=3)
          y1 <- rbinom(n,1,exp(mu[1,])/(exp(mu[1,])+1))
          y2 <- rpois(n, exp(mu[2,]))
          y3 <- rnorm(n, mu[3,])
          dat <- data.frame(y1,y2,y3,a,time=1:n)
          all_dat[[i]] <- dat
          all_x[i,] <- c(X1,X2)
        }
        list(y=all_dat, x=all_x)
      }
      save.image(paste0('./simu_sample_data/','N',N,'T',n,'p',p,'.Rdata')) # save simulated data
      rm(list = setdiff(ls(),c('n','N','p','N_set','p_set','iterN','cl')))
    }
  }
}

stopCluster(cl)


# mixed effect model -------------------------------------------------------------------
rm(list=ls())
input <- list.files(path = './simu_sample_data/')
input <- input[grepl('N',input)]
glmm <- vector('list',0)
for (fl in input){
  load(fl)
  print(fl)
  coeff <- NULL
  for (i in 1:iterN){
    x <- all_res_var[[i]]$x
    y <- all_res_var[[i]]$y
    A_dat <- do.call(rbind, y)
    A_dat$id <- rep(1:length(y), each=NROW(y[[1]]))
    x$id <- 1:length(y)
    A_dat <- merge(A_dat, x)
    suppressWarnings(m1 <- glmer(y1~X1+X2+time+(a+1|id), data=A_dat, family = binomial()))
    suppressWarnings(m2 <- glmer(y2~X1+X2+time+(a+1|id), data=A_dat, family = poisson()))
    m3 <- lmer(y3~X1+X2+time+(a+1|id), data=A_dat)
    coeff <- rbind(coeff,
                   c(summary(m1)$coefficients[,1], summary(m2)$coefficients[,1], summary(m3)$coefficients[,1]))
  }
  glmm[[fl]] <- coeff
}



# proposed method ---------------------------------------------------------
# rm(list=ls())
source("../../mrss.R")
input <- list.files(path = './simu_sample_data/')
input <- input[grepl('N',input)]
ss <- vector('list',0)

cl <- parallel::makeCluster(8)
registerDoParallel(cl)
for (fl in input){
  load(fl)
  print(fl)
  
  state <- data.frame(ind = c('alpha_1','alpha_2'),stringsAsFactors = F)
  state <- cbind(state, matrix(c(1,0,0,1.5),ncol=2))
  response <- data.frame(var = c('y1','y2','y3'), stringsAsFactors = F)
  response$scale <- c(0,0,0)
  response$dis <- c('binomial','poisson','gaussian')
  response$u <- c(1,1,1)
  response <- cbind(response, data.frame(alpha_1=rep('a',3),alpha_2=rep(1,3),stringsAsFactors = F))
  coeff <-  foreach(j = 1:iterN, .combine = 'rbind', .packages = c('KFAS'), .errorhandling = 'remove') %dopar% {
    x <- all_res_var[[j]]$x
    y <- all_res_var[[j]]$y
    
    x$ID <- as.character(1:nrow(x))
    names(y) <- as.character(1:nrow(x))
    ss_model <- mrss(ID='ID', state=state, covariate=c('X1','X2'), covariate_scale = F, 
                     time_covariate = 'time', response=response, demo=x, data=y, day='time')
    ss_model <- list(ss_model)
    all_para <- get_para(ss_model)
    all_para_vec <- para_vec(all_para)
    all_para_vec$x <- c(0.6, 0.8, 1.2, 2, 1, 0, 1.5, 1, -0.5, 0.2, -1, 0.1, 0.2, 1, 1, 1, 1, 2, 2, 2, 0.03, 0.03, 0.03)
    call_ss <- function(pars,trace){
      all_para_vec$x[5:23] <- pars
      all_para <- vec_para(all_para, all_para_vec)
      ss_model <- plug_para(ss_model, all_para)
      callike(ss_model,trace)
    }
    res <- optim(c(1,0,1.5,1,-0.5,0.2,-1,0.1,0.2,1,1,1,1, 2,2,2, 0.03,0.03,0.03), fn=call_ss, method = 'BFGS',
                 control = list(fnscale=-1), trace=NULL)
    res$par
  }
  ss[[fl]] <- coeff
}
stopCluster(cl)


# VAR_individual method ---------------------------------------------------
# rm(list=ls())
input <- list.files(path = './simu_sample_data/')
input <- input[grepl('N',input)]
var_ind <- vector('list',0)

for (fl in input){
  load(fl)
  print(fl)
  model <- vector('list',0)
  for (j in 1:iterN){
    y <- all_res_var[[j]]$y
    N <- length(y)
    for (kk in 1:N){
      y[[kk]]$y2 <- log(0.5+y[[kk]]$y2)
    }
    n <- NROW(y[[1]])
    all_sub <- vector('list',N)
    
    for (i in 1:N){
      m <- VAR(y[[i]][,1:3], p=1, exogen = y[[i]][,4,drop=F], type='const')
      possible <- tryCatch(sm <- summary(m), error = function(e)e)
      if (inherits(possible, "error")) next
      Tt <- rbind(matrix(c(m$varresult$y1$coefficients[1:4], m$varresult$y2$coefficients[1:4], m$varresult$y3$coefficients[1:4]), nrow=3,byrow = T),
                  c(0,0,0,1))
      Tt <- array(Tt, dim=c(dim(Tt), NROW(y[[i]])))
      Tt[is.na(Tt)] <- 0
      Tt[1:3,4,y[[i]][,4]==1] <- Tt[1:3,4,y[[i]][,4]==1] + c(m$varresult$y1$coefficients[5], m$varresult$y2$coefficients[5], m$varresult$y3$coefficients[5])
      Tt[is.na(Tt)] <- 0
      Qt <- sm$covres
      Rt <- rbind(diag(c(1,1,1)), matrix(0,ncol=3,nrow=1))
      P1 <- 0
      P1inf <- diag(c(1,1,1,0))
      Zt <- cbind(diag(1,3),c(0,0,0))
      
      Ht <- matrix(0,nrow=3,ncol=3)
      a1c <- c(as.numeric(y[[i]][1,1:3]),1)
      all_sub[[i]] <- SSModel(cbind(y1,y2,y3) ~  -1+SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt,
                                                              a1=a1c, P1inf=P1inf, P1=P1,
                                                              state_names = c('a1','a2','a3','int')),
                              data = y[[i]], distribution = 'gaussian', H = Ht)
    }
    model[[j]] <- all_sub
  }
  var_ind[[fl]] <- model
}


# VAR_pool method ---------------------------------------------------------
# rm(list=ls())
input <- list.files(path = './simu_sample_data/')
input <- input[grepl('N',input)]
var_pool <- vector('list',0)

for (fl in input){
  load(fl)
  print(fl)
  model <- vector('list',0)
  for (j in 1:iterN){
    y <- all_res_var[[j]]$y
    N <- length(y)
    for (kk in 1:N){
      y[[kk]]$y2 <- log(0.5+y[[kk]]$y2)
    }
    n <- NROW(y[[1]])
    all_sub <- vector('list',N)
    allT_s <- array(NA, dim=c(3,5,N))
    allQ <- array(NA, dim=c(3,3,N))
    
    for (i in 1:N){
      m <- VAR(y[[i]][,1:3], p=1, exogen = y[[i]][,4,drop=F], type='const')
      possible <- tryCatch(sm <- summary(m), error = function(e)e)
      if (inherits(possible, "error")) next
      allT_s[,,i] <- rbind(matrix(c(m$varresult$y1$coefficients, m$varresult$y2$coefficients, m$varresult$y3$coefficients), nrow=3,byrow = T))
      allQ[,,i] <- sm$covres
    }
    allTt <- apply(allT_s, c(1,2), function(x)mean(x,na.rm=T))
    Qt <- apply(allQ, c(1,2), function(x)mean(x,na.rm=T))
    for (i in 1:N){
      Tt <- rbind(allTt[,1:4], c(0,0,0,1))
      Tt <- array(Tt, dim=c(dim(Tt), NROW(y[[i]])))
      Tt[1:3,4,y[[i]][,4]==1] <- Tt[1:3,4,y[[i]][,4]==1] + allTt[,5]
      Rt <- rbind(diag(c(1,1,1)), matrix(0,ncol=3,nrow=1))
      P1 <- 0
      P1inf <- diag(c(1,1,1,0))
      Zt <- cbind(diag(1,3),c(0,0,0))
      
      Ht <- matrix(0,nrow=3,ncol=3)
      a1c <- c(as.numeric(y[[i]][1,1:3]),1)
      all_sub[[i]] <- SSModel(cbind(y1,y2,y3) ~  -1+SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt,
                                                              a1=a1c, P1inf=P1inf, P1=P1,
                                                              state_names = c('a1','a2','a3','int')),
                              data = y[[i]], distribution = 'gaussian', H = Ht)
    }
    model[[j]] <- all_sub
  }
  var_pool[[fl]] <- model
}

# Figure 3 (a1)(a2), Supp Figure 3(a) -------------------------------------
rm(list=ls())
# load fitting results of glmm and ss
load("./simu_results/glmm for f3.RData")
load("./simu_results/ss for f3.RData")
iterN <- 100
dat <- data.frame(rbind(glmm$N20n30.RData, glmm$N40n30.RData, glmm$N60n30.RData))[,c(2,3,5,6,8,9)]
dat$N <- rep(c(20,40,60),each=iterN)

tmp <- data.frame(rbind(ss$N20n30, ss$N40n30.RData, ss$N60n30.RData))[,c(6,9,7,10,8,11)]
tmp$N <- rep(c(20,40,60),each=iterN)
colnames(tmp) <- colnames(dat)
dat <- rbind(dat,tmp)
dat$method <- rep(c('glmm','ss'), each=3*iterN)
boxplot(X2~method+N, data=dat,at=c(1,2,4,5,7,8))

aggregate(dat[,1:6], list(dat$N, dat$method), function(x)mean((x-1)^2))
aggregate(dat[,1:6], list(dat$N, dat$method), function(x)sd((x-1)^2))

Y1 <- dat[,c(1,7,8)]
Y1$var <- 'X1'
tmp <- dat[c(2,7,8)]
tmp$var <- 'X2'
tmp$X2 <- tmp$X2-1
colnames(tmp) <- colnames(Y1)
Y1 <- rbind(Y1,tmp)

Y2 <- dat[,c(3,7,8)]
Y2$var <- 'X1'
tmp <- dat[c(4,7,8)]
tmp$var <- 'X2'
tmp$X2.1 <- tmp$X2.1-1
colnames(tmp) <- colnames(Y2)
Y2 <- rbind(Y2,tmp)

Y3 <- dat[,c(5,7,8)]
Y3$var <- 'X1'
tmp <- dat[c(6,7,8)]
tmp$var <- 'X2'
tmp$X2.2 <- tmp$X2.2-1
colnames(tmp) <- colnames(Y3)
Y3 <- rbind(Y3,tmp)

boxplot(X1~method+N+var, data=Y1,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1~method+N+var, data=Y1,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('20','40','60','20','40','60'))
axis(4, at = seq(0.6,1.6,0.2), labels=as.character(format(seq(0.6,1.6,0.2)+1,digits=2)))
text(0.5,1.55,expression('X'[1]), cex=1.5)
text(10.25,1.55,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Number of subjects (N)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 

boxplot(X1.1~method+N+var, data=Y2,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),ylim=c(0.8,1.15), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.1~method+N+var, data=Y2,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('20','40','60','20','40','60'))
axis(4, at = seq(0.8,1.15,0.05), labels=c('1.80','','1.90','','2.00','','2.10',''))
text(0.5,1.14,expression('X'[1]), cex=1.5)
text(10.25,1.14,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Number of subjects (N)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 

boxplot(X1.2~method+N+var, data=Y3,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.2~method+N+var, data=Y3,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('20','40','60','20','40','60'))
axis(4, at = seq(0.5,1.5,0.5), labels=as.character(format(seq(0.5,1.5,0.5)+1,digits=2)))
text(0.5,1.62,expression('X'[1]), cex=1.5)
text(10.25,1.62,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Number of subjects (N)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 


# Figure 3 (b1)(b2), Supp Figure 3(b) -------------------------------------
rm(list=ls())
# load fitting results of glmm and ss
load("./simu_results/glmm for f3.RData")
load("./simu_results/ss for f3.RData")
iterN <- 100
dat <- data.frame(rbind(glmm$N40n15.RData, glmm$N40n30.RData, glmm$N40n60.RData))[,c(2,3,5,6,8,9)]
dat$n <- rep(c(15,30,60),each=iterN)

tmp <- data.frame(rbind(ss$N40n15.RData, ss$N40n30.RData, ss$N40n60.RData))[,c(6,9,7,10,8,11)]
tmp$n <- rep(c(15,30,60),each=iterN)
colnames(tmp) <- colnames(dat)
dat <- rbind(dat,tmp)
dat$method <- rep(c('glmm','ss'), each=iterN*3)

aggregate(dat[,1:6], list(dat$n, dat$method), function(x)mean((x-1)^2))
aggregate(dat[,1:6], list(dat$n, dat$method), function(x)sd((x-1)^2))

Y1 <- dat[,c(1,7,8)]
Y1$var <- 'X1'
tmp <- dat[c(2,7,8)]
tmp$var <- 'X2'
tmp$X2 <- tmp$X2-1
colnames(tmp) <- colnames(Y1)
Y1 <- rbind(Y1,tmp)

Y2 <- dat[,c(3,7,8)]
Y2$var <- 'X1'
tmp <- dat[c(4,7,8)]
tmp$var <- 'X2'
tmp$X2.1 <- tmp$X2.1-1
colnames(tmp) <- colnames(Y2)
Y2 <- rbind(Y2,tmp)

Y3 <- dat[,c(5,7,8)]
Y3$var <- 'X1'
tmp <- dat[c(6,7,8)]
tmp$var <- 'X2'
tmp$X2.2 <- tmp$X2.2-1
colnames(tmp) <- colnames(Y3)
Y3 <- rbind(Y3,tmp)


boxplot(X1~method+n+var, data=Y1,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),ylim=c(0.5,1.4), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1~method+n+var, data=Y1,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('15','30','60','15','30','60'))
axis(4, at = seq(0.6,1.4,0.2), labels=as.character(format(seq(0.6,1.4,0.2)+1,digits=2)))
text(0.5,1.37,expression('X'[1]), cex=1.5)
text(10.25,1.37,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Time series length (T)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0))    

boxplot(X1.1~method+n+var, data=Y2,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), ylim=c(0.85,1.12), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.1~method+n+var, data=Y2,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('15','30','60','15','30','60'))
axis(4, at = seq(0.85,1.1,0.05), labels=as.character(format(seq(0.85,1.1,0.05)+1),digits=2))
text(0.5,1.11,expression('X'[1]), cex=1.5)
text(10.25,1.11,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Time series length (T)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 


boxplot(X1.2~method+n+var, data=Y3,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), ylim=c(0.25,1.5), xlab='', ylab='')
lines(c(-1,21),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.2~method+n+var, data=Y3,at=c(1,2,4,5,7,8, 11,12,14,15,17,18),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=9.5)
axis(1, at = c(1.5, 4.5, 7.5, 11.5,14.5, 17.5), labels=c('15','30','60','15','30','60'))
axis(4, at = seq(0.2,1.4,0.2), labels=as.character(format(seq(0.2,1.4,0.2)+1,digits=2)))
text(0.5,1.45,expression('X'[1]), cex=1.5)
text(10.25,1.45,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = "Time series length (T)", mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 


# Figure 3 (c1)(c2), Supp Figure (3c) -------------------------------------
rm(list=ls())
# load fitting results of glmm and ss
load("./simu_results/glmm for f3.RData")
load("./simu_results/ss for f3.RData")
iterN <- 100
dat <- data.frame(rbind(glmm$N40n30_a0.RData, glmm$N40n30_a0.15.RData, glmm$N40n30.RData, glmm$N40n30_a0.45.RData))[,c(2,3,5,6,8,9)]
dat$a <- rep(c(0,15,30,45),each=iterN)

tmp <- data.frame(rbind(ss$N40n30_a0.RData, ss$N40n30_a0.15.RData, ss$N40n30.RData, ss$N40n30_a0.45.RData))[,c(6,9,7,10,8,11)]
tmp$a <- rep(c(0,15,30, 45),each=iterN)
colnames(tmp) <- colnames(dat)
dat <- rbind(dat,tmp)
dat$method <- rep(c('glmm','ss'), each=400)
boxplot(X2~method+a, data=dat)

Y1 <- dat[,c(1,7,8)]
Y1$var <- 'X1'
tmp <- dat[c(2,7,8)]
tmp$var <- 'X2'
tmp$X2 <- tmp$X2-1
colnames(tmp) <- colnames(Y1)
Y1 <- rbind(Y1,tmp)

Y2 <- dat[,c(3,7,8)]
Y2$var <- 'X1'
tmp <- dat[c(4,7,8)]
tmp$var <- 'X2'
tmp$X2.1 <- tmp$X2.1-1
colnames(tmp) <- colnames(Y2)
Y2 <- rbind(Y2,tmp)

Y3 <- dat[,c(5,7,8)]
Y3$var <- 'X1'
tmp <- dat[c(6,7,8)]
tmp$var <- 'X2'
tmp$X2.2 <- tmp$X2.2-1
colnames(tmp) <- colnames(Y3)
Y3 <- rbind(Y3,tmp)

boxplot(X1~method+a+var, data=Y1,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6))
lines(c(-1,26),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1~method+a+var, data=Y1,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=12.5)
axis(1, at = c(1.5, 4.5, 7.5, 10.5, 14.5, 17.5, 20.5, 23.5), labels=c('0','0.15','0.3','0.45','0','0.15','0.3','0.45'))
axis(4, at = seq(0.4,1.4,0.2), labels=as.character(format(seq(0.4,1.4,0.2)+1,digits=2)))
text(0.5,1.42,expression('X'[1]), cex=1.5)
text(13.4,1.42,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = expression(paste("Expectation of ",a[t],' (p)')), mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 

boxplot(X1.1~method+a+var, data=Y2,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6))
lines(c(-1,26),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.1~method+a+var, data=Y2,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6),add=T)
abline(v=12.5)
axis(1, at = c(1.5, 4.5, 7.5, 10.5, 14.5, 17.5, 20.5, 23.5), labels=c('0','0.15','0.3','0.45','0','0.15','0.3','0.45'))
axis(4, at = seq(0.9,1.1,0.05), labels=as.character(format(seq(0.9,1.1,0.05)+1,digits=3)))
text(0.5,1.11,expression('X'[1]), cex=1.5)
text(13.4,1.11,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = expression(paste("Expectation of ",a[t],' (p)')), mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 

boxplot(X1.2~method+a+var, data=Y3,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), ylim=c(0.50,1.5))
lines(c(-1,26),c(1,1),lty=2, col='grey60',lwd=2)
boxplot(X1.2~method+a+var, data=Y3,at=c(1,2,4,5,7,8,10,11, 14,15,17,18,20,21,23,24),xaxt='n',col=rep(c('#FBDE8D','#4A83B7'),6), add=T)
abline(v=12.5)
axis(1, at = c(1.5, 4.5, 7.5, 10.5, 14.5, 17.5, 20.5, 23.5), labels=c('0','0.15','0.3','0.45','0','0.15','0.3','0.45'))
axis(4, at = seq(0.6,1.4,0.2), labels=as.character(format(seq(0.6,1.4,0.2)+1,digits=2)))
text(0.5,1.46,expression('X'[1]), cex=1.5)
text(13.4,1.46,expression('X'[2]), cex=1.5)
legend('bottomleft',legend = c('glmm','MRSS','true value'),fill=c('#FBDE8D','#4A83B7',NA), 
       lty=c(NA,NA,2), col=c('#FBDE8D','#4A83B7','grey60'),lwd=2,border =NA, merge=T)
title(xlab = expression(paste("Expectation of ",a[t],' (p)')), mgp = c(2.2, 0, 0))    # Add x-axis text
title(ylab = "Estimated coefficient", mgp = c(2.2, 1, 0)) 



# Figure 4 preprocess -----------------------------------------------------
rm(list=ls())
# load the fitting result by using all 30 time points
load("./simu_results/var_ind.RData")
load("./simu_results/var_pool.RData")
load("./simu_results/ss for f3.RData")
all_varm <- var_ind$N40n30.RData
all_varmpool <- var_pool$N40n30.RData
all_ss <- ss$N40n30.RData
# load the fitting result by using first 25 time points
load("./simu_results/25var_ind.RData")
load("./simu_results/25var_pool.RData")
load("./simu_results/25ss.RData")
varm <- var_ind$N40n30.RData
varmpool <- var_pool$N40n30.RData
ss <- ss$N40n30.RData
# load all predictors and all natural parameters
load("./simu_results/N40n30.RData")
load("./simu_results/N40n30_mu.RData")
rm(list=c('var_ind','var_pool','cl','iterN'))

get_ss <- function(ss_rep, x, y){
  n <- NROW(y)
  Tt <- matrix(0, ncol = 5, nrow = 5)
  diag(Tt) <- c(0.6, 0.8, 1,1,1)
  Tt[1,3] <- 1.2
  Tt[2,3] <- 2
  Qt <- matrix(0, nrow = 2, ncol = 2)
  diag(Qt) <- c(1, 1.5)
  P1 <- 0
  P1inf <- matrix(0, ncol = 5, nrow = 5)
  Zt <- array(c(-0.5,0.2,-1,0.1,0.2,1,0,0,0,1,1,1,2,2,2), dim = c(3,5,n))
  Zt[,1,y$a==0] <- 0
  Rt <- rbind(diag(c(1,1)), matrix(0,ncol=2,nrow=3))
  u <- matrix(1,ncol=3,nrow=n)
  a1c <- c(10,10,1,as.numeric(x))
  model <- SSModel(cbind(y1,y2,y3) ~  -1+SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt,
                                                   a1=a1c, P1inf=P1inf, P1=P1,
                                                   state_names = c('a1','a2','int','X1','X2')),
                   data = y, distribution = c('binomial','poisson','gaussian'), u = u)
  
  
  diag(model["Q", etas = "custom"]) <- ss_rep[1:2]
  model["Z", etas = "custom"][,2,] <- ss_rep[3:5]
  model["Z", etas = "custom"][,4,] <- ss_rep[6:8]
  model["Z", etas = "custom"][,5,] <- ss_rep[9:11]
  model
}

topre_last5 <- function(obj,ssmodel=F){
  # predict last 5 observations
  if (!ssmodel){
    obj$y <-obj$y[26:30,]
    obj$T <- obj$T[,,26:30]
    obj$y[,2] <- NA
    obj$y[,1] <- NA
    obj$y[,3] <- NA
    attributes(obj)$n <- 5L
  }else{
    obj$y <- obj$y[26:30,]
    obj$Z <- obj$Z[,,26:30]
    obj$u <- obj$u[26:30,]
    obj$y[,2] <- NA
    obj$y[,1] <- NA
    obj$y[,3] <- NA
    attributes(obj)$n <- 5L
  }
  return(obj)
}


expit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  x <- log(x/(1-x))
  return(x)
}


iterN <- 100
all_in_err <-  matrix(NA, ncol=9, nrow=100)
all_out_err <- matrix(NA, ncol=9, nrow=100)
for (j in 1:iterN){
  print(j)
  # extract model
  varm_rep <- varm[[j]]
  varmpool_rep <- varmpool[[j]]
  ss_rep <- ss[j,]
  all_varm_rep <- all_varm[[j]]
  all_varmpool_rep <- all_varmpool[[j]]
  all_ss_rep <- all_ss[j,]
  
  in_err_rep <- matrix(NA, ncol=9, nrow=40)
  out_err_rep <- matrix(NA, ncol=9, nrow=40)
  for (i in 1:40){
    varm_s <- varm_rep[[i]]
    varmpool_s <- varmpool_rep[[i]]
    ss_s <- get_ss(ss_rep, all_res_var[[j]]$x[i,], all_res_var[[j]]$y[[i]][1:25,])
    all_varm_s <- all_varm_rep[[i]]
    all_varmpool_s <- all_varmpool_rep[[i]]
    all_ss_s <- get_ss(all_ss_rep, all_res_var[[j]]$x[i,], all_res_var[[j]]$y[[i]])
    
    time <- 0.03*(1:30)
    a <- 1-all_mu[[j]]$y[[i]]$a
    #### for response y2
    # var ind
    mu2 <- all_mu[[j]]$y[[i]]$y2+time
    if (!is.null(varm_s)){
      varm_p <- predict(varm_s, filtered = T, interval = 'confidence')$y2
      in_err_rep[i, 4] <- mean(((varm_p[,1]+time[1:25]-mu2[1:25])[2:25]^2)[a[2:25]==1])
      varm_pp <- predict(varm_s, newdata = topre_last5(all_varm_s), interval = 'confidence')$y2
      out_err_rep[i, 4] <- mean(((varm_pp[,1]+time[26:30]-mu2[26:30])^2)[a[26:30]==1])
    }
    # var pool
    varmpool_p <- predict(varmpool_s, filtered = T, interval = 'confidence')$y2
    in_err_rep[i, 5] <- mean(((varmpool_p[,1]+time[1:25]-mu2[1:25])[2:25]^2)[a[2:25]==1])
    varmpool_pp <- predict(varmpool_s, newdata = topre_last5(all_varmpool_s), interval = 'confidence')$y2
    out_err_rep[i, 5] <- mean(((varmpool_pp[,1]+time[26:30]-mu2[26:30])^2)[a[26:30]==1])
    # ss
    ss_p <- predict(ss_s, filtered = T, type = 'link', interval = 'confidence')$y2
    in_err_rep[i, 6] <- mean(((ss_p[,1]+time[1:25]-mu2[1:25])[2:25]^2)[a[2:25]==1])
    ss_pp <- predict(ss_s, newdata = topre_last5(all_ss_s,ssmodel = T), type = 'link', interval = 'confidence')$y2
    out_err_rep[i, 6] <- mean(((ss_pp[,1]+time[26:30]-mu2[26:30])^2)[a[26:30]==1])
    #### for response y3
    # var ind
    mu3 <- all_mu[[j]]$y[[i]]$y3+time
    if (!is.null(varm_s)){
      varm_p <- predict(varm_s, filtered = T, interval = 'confidence')$y3
      in_err_rep[i, 7] <- mean(((varm_p[,1]+time[1:25]-mu3[1:25])[2:25]^2)[a[2:25]==1])
      varm_pp <- predict(varm_s, newdata = topre_last5(all_varm_s), interval = 'confidence')$y3
      out_err_rep[i, 7] <- mean(((varm_pp[,1]+time[26:30]-mu3[26:30])^2)[a[26:30]==1])
    }
    # var pool
    varmpool_p <- predict(varmpool_s, filtered = T, interval = 'confidence')$y3
    in_err_rep[i, 8] <- mean(((varmpool_p[,1]+time[1:25]-mu3[1:25])[2:25]^2)[a[2:25]==1])
    varmpool_pp <- predict(varmpool_s, newdata = topre_last5(all_varmpool_s), interval = 'confidence')$y3
    out_err_rep[i, 8] <- mean(((varmpool_pp[,1]+time[26:30]-mu3[26:30])^2)[a[26:30]==1])
    # ss
    ss_p <- predict(ss_s, filtered = T, type = 'link', interval = 'confidence')$y3
    in_err_rep[i, 9] <- mean(((ss_p[,1]+time[1:25]-mu3[1:25])[2:25]^2)[a[2:25]==1])
    ss_pp <- predict(ss_s, newdata = topre_last5(all_ss_s,ssmodel = T), type = 'link', interval = 'confidence')$y3
    out_err_rep[i, 9] <- mean(((ss_pp[,1]+time[26:30]-mu3[26:30])^2)[a[26:30]==1])
    #### for response y1
    # var ind
    mu1 <- all_mu[[j]]$y[[i]]$y1+time
    if (!is.null(varm_s)){
      varm_p <- predict(varm_s, filtered = T, interval = 'confidence')$y1
      in_err_rep[i, 1] <- mean(((logit(varm_p[,1])+time[1:25]-mu1[1:25])[2:25]^2)[a[2:25]==1],na.rm = T)
      varm_pp <- predict(varm_s, newdata = topre_last5(all_varm_s), interval = 'confidence')$y1
      out_err_rep[i, 1] <- mean(((logit(varm_pp[,1])+time[26:30]-mu1[26:30])^2)[a[26:30]==1],na.rm = T)
    }
    # var pool
    varmpool_p <- predict(varmpool_s, filtered = T, interval = 'confidence')$y1
    in_err_rep[i, 2] <- mean(((logit(varmpool_p[,1])+time[1:25]-mu1[1:25])[2:25]^2)[a[2:25]==1],na.rm = T)
    varmpool_pp <- predict(varmpool_s, newdata = topre_last5(all_varmpool_s), interval = 'confidence')$y1
    out_err_rep[i, 2] <- mean(((logit(varmpool_pp[,1])+time[26:30]-mu1[26:30])^2)[a[26:30]==1],na.rm = T)
    # ss
    ss_p <- predict(ss_s, filtered = T, type = 'link', interval = 'confidence')$y1
    in_err_rep[i, 3] <- mean(((ss_p[,1]+time[1:25]-mu1[1:25])[2:25]^2)[a[2:25]==1],na.rm = T)
    ss_pp <- predict(ss_s, newdata = topre_last5(all_ss_s,ssmodel = T), type = 'link', interval = 'confidence')$y1
    out_err_rep[i, 3] <- mean(((ss_pp[,1]+time[26:30]-mu1[26:30])^2)[a[26:30]==1],na.rm = T)
  }
  in_err_rep[is.infinite(in_err_rep)] <- NA
  all_in_err[j,] <- colMeans(in_err_rep, na.rm = T)
  out_err_rep[is.infinite(out_err_rep)] <- NA
  all_out_err[j,] <- colMeans(out_err_rep, na.rm = T)
}

save.image('./simu_results/total_prediction.RData')

# Figure 4 plot -----------------------------------------------------------
rm(list=ls())
# load results from last section
load("./simu_results/total_prediction.RData")
# columns: var ind (mu1), var pool (mu1), ss(mu1), var ind (mu2), var pool (mu2), ss(mu2), 
# var ind (mu3), var pool (mu3), ss(mu3)

boxplot(log10(all_in_err[,c(2,1,3,5,4,6,8,7,9)]),at=c(1,2,3,5,6,7,9,10,11), 
        col=c('#00B8A9','#F6416C','#FFDE7D'), medlwd = 1.1, yaxt='n', xaxt='n', xlab='Response', 
        ylab=expression(paste(L[2]," loss in natural parameter space")))
axis(2, at = seq(-1,1,0.5), labels=expression(10^-1,10^-0.5,10^-0,10^-0.5,10^1))
legend('bottomright',legend = c('VAR (pool)','VAR (ind)','MRSS'),fill=c('#00B8A9','#F6416C','#FFDE7D'),border = NA)
axis(1, at = c(2,6,10), labels=expression(bold(Y)^(1),bold(Y)^(2),bold(Y)^(3)))
abline(v=4)
abline(v=8)


boxplot(log10(all_out_err[,c(2,1,3,5,4,6,8,7,9)]),at=c(1,2,3,5,6,7,9,10,11), 
        col=c('#00B8A9','#F6416C','#FFDE7D'),  medlwd = 1.1, yaxt='n', xaxt='n', xlab='Response', 
        ylab=expression(paste(L[2]," loss in natural parameter space")))
axis(2, at = seq(-1,1,0.5), labels=expression(10^-1,10^-0.5,10^-0,10^-0.5,10^1))
legend('bottomright',legend = c('VAR (pool)','VAR (ind)','MRSS'),fill=c('#00B8A9','#F6416C','#FFDE7D'),border = NA)
axis(1, at = c(2,6,10), labels=expression(bold(Y)^(1),bold(Y)^(2),bold(Y)^(3)))
abline(v=4)
abline(v=8)

colMeans(all_in_err)
apply(all_in_err, 2, sd)

colMeans(all_out_err)
apply(all_out_err, 2, sd)


# Figure 5 ----------------------------------------------------------------
rm(list=ls())
# load the fitting result by using all 30 time points
load("./simu_results/var_ind.RData")
load("./simu_results/var_pool.RData")
load("./simu_results/ss for f3.RData")
all_varm <- var_ind$N40n30.RData
all_varmpool <- var_pool$N40n30.RData
all_ss <- ss$N40n30.RData
# load the fitting result by using first 25 time points
load("./simu_results/25var_ind.RData")
load("./simu_results/25var_pool.RData")
load("./simu_results/25ss.RData")
varm <- var_ind$N40n30.RData
varmpool <- var_pool$N40n30.RData
ss <- ss$N40n30.RData
# load all predictors and all natural parameters
load("./simu_results/N40n30.RData")
load("./simu_results/N40n30_mu.RData")
rm(list=c('var_ind','var_pool','cl','iterN'))


# use first repetition
j=1
varm_rep <- varm[[j]]
varmpool_rep <- varmpool[[j]]
ss_rep <- ss[j,]
all_varm_rep <- all_varm[[j]]
all_varmpool_rep <- all_varmpool[[j]]
all_ss_rep <- all_ss[j,]

getci <- function(p,pp){
  tmp1 <- c(p[,2]+time[1:25],pp[,2]+time[26:30])
  tmp2 <- c(p[,3]+time[seq(25,1)],pp[,3]+time[seq(30,26)])
  tmp1 <- c(tmp1, tmp2[seq(30,1)])
  tmp1 <- tmp1[c(-1,-60)]
  return(tmp1)
}

getcip <- function(p,pp){
  tmp1 <- trantime(c(p[,2],pp[,2]), time)
  tmp2 <- trantime(c(p[,3],pp[,3]), time[seq(30,1)])
  tmp1 <- c(tmp1, tmp2[seq(30,1)])
  tmp1 <- tmp1[c(-1,-60)]
  return(tmp1)
}

trantime <- function(x, time){
  ind <- x>0 & x<1
  x[ind] <- expit(logit(x[ind])+time[ind])
  return(x)
}

get_ss <- function(ss_rep, x, y){
  n <- NROW(y)
  Tt <- matrix(0, ncol = 5, nrow = 5)
  diag(Tt) <- c(0.6, 0.8, 1,1,1)
  Tt[1,3] <- 1.2
  Tt[2,3] <- 2
  Qt <- matrix(0, nrow = 2, ncol = 2)
  diag(Qt) <- c(1, 1.5)
  P1 <- 0
  P1inf <- matrix(0, ncol = 5, nrow = 5)
  Zt <- array(c(-0.5,0.2,-1,0.1,0.2,1,0,0,0,1,1,1,2,2,2), dim = c(3,5,n))
  Zt[,1,y$a==0] <- 0
  Rt <- rbind(diag(c(1,1)), matrix(0,ncol=2,nrow=3))
  u <- matrix(1,ncol=3,nrow=n)
  a1c <- c(10,10,1,as.numeric(x))
  model <- SSModel(cbind(y1,y2,y3) ~  -1+SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt,
                                                   a1=a1c, P1inf=P1inf, P1=P1,
                                                   state_names = c('a1','a2','int','X1','X2')),
                   data = y, distribution = c('binomial','poisson','gaussian'), u = u)
  
  
  diag(model["Q", etas = "custom"]) <- ss_rep[1:2]
  model["Z", etas = "custom"][,2,] <- ss_rep[3:5]
  model["Z", etas = "custom"][,4,] <- ss_rep[6:8]
  model["Z", etas = "custom"][,5,] <- ss_rep[9:11]
  model
}

topre_last5 <- function(obj,ssmodel=F){
  # predict last 5 observations
  if (!ssmodel){
    obj$y <-obj$y[26:30,]
    obj$T <- obj$T[,,26:30]
    obj$y[,2] <- NA
    obj$y[,1] <- NA
    obj$y[,3] <- NA
    attributes(obj)$n <- 5L
  }else{
    obj$y <- obj$y[26:30,]
    obj$Z <- obj$Z[,,26:30]
    obj$u <- obj$u[26:30,]
    obj$y[,2] <- NA
    obj$y[,1] <- NA
    obj$y[,3] <- NA
    attributes(obj)$n <- 5L
  }
  return(obj)
}


expit <- function(x){
  exp(x)/(1+exp(x))
}

logit <- function(x){
  x <- log(x/(1-x))
  return(x)
}



# plot the 9th subject as an example
i <- 9
varm_s <- varm_rep[[i]]
varmpool_s <- varmpool_rep[[i]]
ss_s <- get_ss(ss_rep, all_res_var[[j]]$x[i,], all_res_var[[j]]$y[[i]][1:25,])
all_varm_s <- all_varm_rep[[i]]
all_varmpool_s <- all_varmpool_rep[[i]]
all_ss_s <- get_ss(all_ss_rep, all_res_var[[j]]$x[i,], all_res_var[[j]]$y[[i]])
time <- 0.03*(1:30)
a <- 1-all_mu[[j]]$y[[i]]$a

#### figure 5(a)
# var ind
mu2 <- all_mu[[j]]$y[[i]]$y2+time
print(i)
varm_p <- predict(varm_s, filtered = T, interval = 'confidence')$y2
varm_pp <- predict(varm_s, newdata = topre_last5(all_varm_s), interval = 'confidence')$y2

# var pool
varmpool_p <- predict(varmpool_s, filtered = T, interval = 'confidence')$y2
varmpool_pp <- predict(varmpool_s, newdata = topre_last5(all_varmpool_s), interval = 'confidence')$y2

# ss
ss_p <- predict(ss_s, filtered = T, type = 'link', interval = 'confidence')$y2
ss_pp <- predict(ss_s, newdata = topre_last5(all_ss_s,ssmodel = T), type = 'link', interval = 'confidence')$y2

lim <- c(getci(varm_p, varm_pp), getci(varmpool_p, varmpool_pp), getci(ss_p, ss_pp))
plot(mu2,type='p', cex=0.8, xlim=c(0,30),ylim=c(min(lim),max(lim)), ylab=expression(paste("Predicted ",mu^(2))),
     xlab='Time')
abline(v=25,lty=2,col='red')
lines(c(varm_p[,1]+time[1:25],varm_pp[,1]+time[26:30]), col='#F6416C',lwd=2)
polygon(c(2:30,seq(30,2)),getci(varm_p, varm_pp), border = NA, col=adjustcolor('#F6416C',alpha.f = 0.15))

lines(c(varmpool_p[,1]+time[1:25], varmpool_pp[,1]+time[26:30]), col='#00B8A9',lwd=2)
polygon(c(2:30,seq(30,2)),getci(varmpool_p, varmpool_pp), border = NA, col=adjustcolor('#00B8A9',alpha.f = 0.15))

lines(c(ss_p[,1]+time[1:25],ss_pp[,1]+time[26:30]), col='#1f7bcf',lwd=2)
polygon(c(2:30,seq(30,2)),getci(ss_p, ss_pp), border = NA, col=adjustcolor('#1f7bcf',alpha.f = 0.2))
points(mu2,type='p', ylim=c(-2,4),cex=0.8, pch=16,col=1+all_mu[[j]]$y[[i]]$a)


legend('topleft',c('VAR (pool)','VAR (ind)','MRSS'), lty=1, col=c('#00B8A9','#F6416C','#1f7bcf'),lwd=2)

#### figure 5(b)
# var ind
varm_p <- predict(varm_s, filtered = T, interval = 'confidence')$y3
varm_pp <- predict(varm_s, newdata = topre_last5(all_varm_s), interval = 'confidence')$y3
varmpool_p <- predict(varmpool_s, filtered = T, interval = 'confidence')$y3
varmpool_pp <- predict(varmpool_s, newdata = topre_last5(all_varmpool_s), interval = 'confidence')$y3
ss_p <- predict(ss_s, filtered = T, type = 'link', interval = 'confidence')$y3
ss_pp <- predict(ss_s, newdata = topre_last5(all_ss_s,ssmodel = T), type = 'link', interval = 'confidence')$y3

lim <- c(getci(varm_p, varm_pp), getci(varmpool_p, varmpool_pp), getci(ss_p, ss_pp))
plot(all_mu[[j]]$y[[i]]$y3+time,type='p', cex=0.8, ylab=expression(paste("Predicted ",mu^(3))), xlim=c(0,30),ylim=c(min(lim),max(lim)),
     xlab='Time')
abline(v=25,lty=2,col='red')
lines(c(varm_p[,1]+time[1:25],varm_pp[,1]+time[26:30]), col='#F6416C',lwd=2)
polygon(c(2:30,seq(30,2)),getci(varm_p, varm_pp), border = NA, col=adjustcolor('#F6416C',alpha.f = 0.15))

lines(c(varmpool_p[,1]+time[1:25], varmpool_pp[,1]+time[26:30]), col='#00B8A9',lwd=2)
polygon(c(2:30,seq(30,2)),getci(varmpool_p, varmpool_pp), border = NA, col=adjustcolor('#00B8A9',alpha.f = 0.15))

lines(c(ss_p[,1]+time[1:25],ss_pp[,1]+time[26:30]), col='#1f7bcf',lwd=2)
polygon(c(2:30,seq(30,2)),getci(ss_p, ss_pp), border = NA, col=adjustcolor('#1f7bcf',alpha.f = 0.2))
points(all_mu[[j]]$y[[i]]$y3+time,type='p', ylim=c(-2,4),cex=0.8, pch=16,col=1+all_mu[[j]]$y[[i]]$a)

legend('topleft',c('VAR (pool)','VAR (ind)','MRSS'), lty=1, col=c('#00B8A9','#F6416C','#1f7bcf'),lwd=2)


