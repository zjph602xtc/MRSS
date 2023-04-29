library(lubridate)
library('KFAS')
library('foreach')

# clean demo --------------------------------------------------------------
demo <- read.csv('./mpower_data/other/demographics.csv',
                 colClasses = c(createdOn='POSIXct',recordId='character',healthCode='character'))
demo$professional.diagnosis <- as.numeric(demo$professional.diagnosis)
demo$gender <- as.numeric(demo$gender)-1 # 1 is male, 0 is female
demo$smoked <- as.numeric(demo$smoked)
demo$deep.brain.stimulation <- as.numeric(demo$deep.brain.stimulation)

demo <- demo[!is.na(demo$age),]
demo$smoked[!is.na(demo$packs.per.day)] <- demo$packs.per.day[!is.na(demo$packs.per.day)]
demo <- demo[!is.na(demo$smoked),]
demo$med <- ifelse(demo$medication.start.year > 100, 1, ifelse(!is.na(demo$medication.start.year),0, NA))
usedid <- ((demo$med==1) & (demo$professional.diagnosis==1)) |
  ((demo$med==0) & (demo$professional.diagnosis==1)) |
  ((demo$med==0) & (demo$professional.diagnosis==0))
usedid[is.na(usedid)] <- FALSE
demo <- demo[usedid, ]
usedid <- ((demo$professional.diagnosis==1) & !is.na(demo$deep.brain.stimulation)) |
  ((demo$deep.brain.stimulation==0) & (demo$professional.diagnosis==0))
usedid[is.na(usedid)] <- FALSE
demo <- demo[usedid, ]
demo$duration <- 2021-demo$onset.year
usedid <- ((demo$professional.diagnosis==1) & !is.na(demo$onset.year)) |
  ((demo$professional.diagnosis==0))
demo <- demo[usedid, ]
demo$duration[demo$professional.diagnosis==0] <- 0

series_l <- 10

# clean tapping -----------------------------------------------------------
tap <- read.csv('./mpower_data/other/tapping.csv', 
                colClasses = c(createdOn='POSIXct',healthCode='character'))
tapF <- read.csv('./mpower_data/extracted features/tapFeatures.csv')
colnames(tapF)[1] <- 'tapping_results.json.TappingSamples'
tap <- merge(tap, tapF, by = 'tapping_results.json.TappingSamples', all = F)
tap <- tap[order(tap$healthCode,tap$createdOn),]

t_sub_all <- split(tap[,6,drop=F], tap$healthCode)
t_sub_all <- t_sub_all[sapply(t_sub_all, function(x)NROW(x)>=series_l)]

usedid <- sapply(t_sub_all, function(x)x$healthCode[1])
demo <- demo[demo$healthCode %in% usedid,]

tap <- tap[tap$healthCode %in% demo$healthCode,]
t_sub_all <- split(tap, tap$healthCode)
rm(list=c('tap','tapF','usedid'))


# clean walking -----------------------------------------------------------
walk <- read.csv('./mpower_data/other/walking.csv', 
                 colClasses = c(createdOn='POSIXct', healthCode='character'))

walking <- read.csv('./mpower_data/extracted features/walkFeatures.csv')
walking_r <- read.csv('./mpower_data/extracted features/walk_rest_Features.csv')

walk <- merge(walk, walking, by = 'accel_walking_outbound.json.items', all = F)
colnames(walking_r) <- paste0(colnames(walking_r),'_r')
walk <- merge(walk, walking_r, by.x = 'accel_walking_rest.json.items', by.y = 'accel_walking_rest.json.items_r',all = F)
rm(list=c('walking','walking_r'))
walk <- walk[order(walk$healthCode,walk$createdOn),]

w_sub_all <- split(walk[,7,drop=F], walk$healthCode)
w_sub_all <- w_sub_all[sapply(w_sub_all, function(x)NROW(x)>=series_l)]

usedid <- sapply(w_sub_all, function(x)x$healthCode[1])
demo <- demo[demo$healthCode %in% usedid,]

walk <- walk[walk$healthCode %in% demo$healthCode,]
w_sub_all <- split(walk, walk$healthCode)
rm(list=c('walk','usedid'))

t_sub_all <- t_sub_all[names(t_sub_all) %in% demo$healthCode]

# pool tapping and walking ------------------------------------------------
data <- list()
for (i in 1:NROW(demo)){
  print(i)
  id <- demo$healthCode[i]
  ### walking
  sub <- w_sub_all[[id]]
  # inresp <- resp[resp %in% colnames(sub)]
  # sub <- sub[,c('healthCode', 'createdOn', inresp, 'medTimepoint')]
  # sub[,inresp] <- t(apply(sub[,inresp], 1, function(x)(x-colMeans(scale_w))/apply(scale_w, 2, sd)))
  sub$afternoon <- as.numeric(hour(sub$createdOn)>9)
  sub$createdOn <- round(sub$createdOn, units='day')
  v_sub <- data.frame(date = rep(seq(min(sub[,'createdOn']), max(sub[,'createdOn']), by = 'day'),each=4))
  sub$mtrt <- as.numeric(sub$medTimepoint=='Just after Parkinson medication (at your best)')
  sub <- aggregate(sub, by = list(as.factor(sub$createdOn), sub$mtrt, sub$afternoon), FUN = function(x)mean(x,na.rm=T))
  
  v_sub$medTimepoint <- c(0,1)
  v_sub$afternoon <- c(0,0,1,1)
  v_sub <- merge(v_sub, sub[,1:(NCOL(sub)-2)], by.x = c('date','medTimepoint','afternoon'), by.y = c('Group.1','Group.2','Group.3'), all.x = T, all.y = F)
  v_sub$mtrt_w <- v_sub$medTimepoint
  v_sub$mtrt_w[is.na(v_sub$medianAA_r)] <- NA
  v_sub$after_w <- v_sub$afternoon
  v_sub$after_w[is.na(v_sub$medianAA_r)] <- NA
  ### tapping
  sub <- t_sub_all[[id]]
  # inresp <- resp[resp %in% colnames(sub)]
  # sub <- sub[,c('healthCode', 'createdOn', inresp, 'medTimepoint')]
  sub$rangeTapInter <- as.numeric(sub$rangeTapInter)
  # sub[,inresp] <- t(apply(sub[,inresp], 1, function(x)(x-colMeans(scale_t))/apply(scale_t, 2, sd)))
  sub$afternoon <- as.numeric(hour(sub$createdOn)>9)
  sub$createdOn <- round(sub$createdOn, units='day')
  tmp <- data.frame(date = rep(seq(min(sub[,'createdOn']), max(sub[,'createdOn']), by = 'day'),each=4))
  sub$mtrt <- as.numeric(sub$medTimepoint=='Just after Parkinson medication (at your best)')
  sub <- aggregate(x = sub, by = list(as.factor(sub$createdOn), sub$mtrt, sub$afternoon), FUN = function(x)mean(x,na.rm=T))
  
  tmp$medTimepoint <- c(0,1)
  tmp$afternoon <- c(0,0,1,1)
  tmp <- merge(tmp, sub[,1:(NCOL(sub)-2)], by.x = c('date','medTimepoint','afternoon'), by.y = c('Group.1','Group.2','Group.3'), all.x = T, all.y = F)
  tmp$mtrt_t <- tmp$medTimepoint
  tmp$mtrt_t[is.na(tmp$rangeTapInter)] <- NA
  tmp$after_t <- tmp$afternoon
  tmp$after_t[is.na(tmp$rangeTapInter)] <- NA
  
  v_sub <- merge(v_sub, tmp, by = c('date','afternoon','medTimepoint'), all = T)
  v_sub <- v_sub[order(v_sub$date, v_sub$afternoon, v_sub$medTimepoint),]
  v_sub$time <- rep(1:(NROW(v_sub)/4),each=4)
  
  data[[id]] <- v_sub
}

rm(list = c('id', 'sub', 'tmp', 'v_sub', 'series_l'))
rm(list = c('t_sub_all', 'w_sub_all'))
## x,y not scaled yet


# build model -------------------------------------------------------------
state <- data.frame(ind = c('b1_trt','b2_morning','v1','v2'),stringsAsFactors = F)
state <- cbind(state, matrix(c(rep(1,11),0,1,1,0,1),ncol=4))
covariate <- c('age','gender','smoked','duration','deep.brain.stimulation')
time_covariate <- c('time')
response <- data.frame(var = c('mtrt_t','rangeTapInter','numberTaps',
                               'sdTapInter','mtrt_w','medianAA_r','sdAJ','dfaY'), stringsAsFactors = F)
response$scale <- c(0,1,1,1,0,1,1,1)
response$dis <- c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian')
response$u <- c(1,1,1,1,1,1,1,1)
response <- cbind(response, data.frame(b1_trt=c(0,'mtrt_t','mtrt_t','mtrt_t',0,'mtrt_w','mtrt_w','mtrt_w'),
                                       b2_morning=c('after_t','after_t','after_t','after_t',
                                                    'after_w','after_w','after_w','after_w'),
                                       v1=rep(1,8),v2=rep(1,8),stringsAsFactors = F))
source('../mrss.R')
md1 <- mrss(ID = 'healthCode', state = state, covariate = covariate, time_covariate = 'time', 
     response = response, demo = demo, subset = 'professional.diagnosis==1 & med==1',
     data = data, day = 'time')

state <- data.frame(ind = c('b2_morning','v1','v2'),stringsAsFactors = F)
state <- cbind(state, matrix(c(rep(1,5),0,1,0,1),ncol=3))
response <- data.frame(var = c('mtrt_t','rangeTapInter','numberTaps',
                               'sdTapInter','mtrt_w','medianAA_r','sdAJ','dfaY'), stringsAsFactors = F)
response$scale <- c(0,1,1,1,0,1,1,1)
response$dis <- c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian')
response$u <- c(1,1,1,1,1,1,1,1)
response <- cbind(response, data.frame(b2_morning=c('after_t','after_t','after_t','after_t',
                                                    'after_w','after_w','after_w','after_w'),
                                       v1=rep(1,8),v2=rep(1,8),stringsAsFactors = F))
md2 <- mrss(ID = 'healthCode', state = state, covariate = covariate, time_covariate = 'time', 
            response = response, demo = demo, subset = 'professional.diagnosis==1 & med==0',
            data = data, day = 'time')

state <- data.frame(ind = c('b2_morning','v1'),stringsAsFactors = F)
state <- cbind(state, matrix(c(rep(1,4)),ncol=2))
covariate <- c('age','gender','smoked')
response <- data.frame(var = c('mtrt_t','rangeTapInter','numberTaps',
                               'sdTapInter','mtrt_w','medianAA_r','sdAJ','dfaY'), stringsAsFactors = F)
response$scale <- c(0,1,1,1,0,1,1,1)
response$dis <- c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian')
response$u <- c(1,1,1,1,1,1,1,1)
response <- cbind(response, data.frame(b2_morning=c('after_t','after_t','after_t','after_t',
                                                    'after_w','after_w','after_w','after_w'),
                                       v1=rep(1,8),stringsAsFactors = F))
md3 <- mrss(ID = 'healthCode', state = state, covariate = covariate, time_covariate = 'time', 
            response = response, demo = demo, subset = 'professional.diagnosis==0 & med==0',
            data = data, day = 'time')

new_md <- list(md1, md2, md3)
all_para <- get_para(new_md)
all_para$mode[all_para$name=='time'] <- 2
# do not reorder all_para !!!!


# fit para ----------------------------------------------------------------
all_para_vec <- para_vec(all_para)
rownames(all_para_vec) <- 1:NROW(all_para_vec)
require(beepr)


call_ss <- function(pars,trace){
  all_para_vec$x <- pars
  all_para <- vec_para(all_para, all_para_vec)
  ss_model <- plug_para(new_md, all_para)
  callike(new_md,trace)
}
res <- optim(all_para_vec$x, fn=call_ss, method = 'BFGS',
             control = list(fnscale=-1), trace=NULL)


like <- NULL
range <- seq(-1,1,length.out=50)
# select a variable for optimize
# i.e., k <- 1
for (i in range){
  print(i)
  all_para_vec$x[k] <- i
  all_para <- vec_para(all_para, all_para_vec)
  new_md <- plug_para(new_md, all_para)
  like <- c(like, callike(new_md))
}
plot(range, log(-like))


# Figure 6 ----------------------------------------------------------------
#### glmm
library(lubridate)
library('KFAS')
library('foreach')
demo <- read.csv('./mpower_data/other/demographics.csv',colClasses = c(createdOn='POSIXct',recordId='character',healthCode='character'), stringsAsFactors = T)
demo$professional.diagnosis <- as.numeric(demo$professional.diagnosis)
demo$gender <- as.numeric(demo$gender)-1
demo$smoked <- as.numeric(demo$smoked)
demo$deep.brain.stimulation <- as.numeric(demo$deep.brain.stimulation)

demo <- demo[!is.na(demo$age),]
demo$smoked[!is.na(demo$packs.per.day)] <- demo$packs.per.day[!is.na(demo$packs.per.day)]
demo <- demo[!is.na(demo$smoked),]
demo$med <- ifelse(demo$medication.start.year > 100, 1, ifelse(!is.na(demo$medication.start.year),0, NA))
usedid <- ((demo$med==1) & (demo$professional.diagnosis==1)) |
  ((demo$med==0) & (demo$professional.diagnosis==1)) |
  ((demo$med==0) & (demo$professional.diagnosis==0))
usedid[is.na(usedid)] <- FALSE
demo <- demo[usedid, ]
demo$deep.brain.stimulation
usedid <- ((demo$professional.diagnosis==1) & !is.na(demo$deep.brain.stimulation)) |
  ((demo$deep.brain.stimulation==0) & (demo$professional.diagnosis==0))
usedid[is.na(usedid)] <- FALSE
demo <- demo[usedid, ]
demo$duration <- 2021-demo$onset.year
usedid <- ((demo$professional.diagnosis==1) & !is.na(demo$onset.year)) |
  ((demo$professional.diagnosis==0))
demo <- demo[usedid, ]
demo$duration[demo$professional.diagnosis==0] <- 0

# tapping
tap <- read.csv('./mpower_data/other/tapping.csv', colClasses = c(createdOn='POSIXct'))
tap$healthCode <- as.character(tap$healthCode)
tapF <- read.csv('./mpower_data/extracted features/tapFeatures.csv')
colnames(tapF)[1] <- 'tapping_results.json.TappingSamples'
tap <- merge(tap, tapF, by = 'tapping_results.json.TappingSamples', all = F)
tap <- tap[order(tap$healthCode,tap$createdOn),]

t_sub_all <- split(tap[,6,drop=F], tap$healthCode)
t_sub_all <- t_sub_all[sapply(t_sub_all, function(x)NROW(x)>=20)]

usedid <- sapply(t_sub_all, function(x)x$healthCode[1])
demo <- demo[demo$healthCode %in% usedid,]

tap <- tap[tap$healthCode %in% demo$healthCode,]

demo$age <- scale(demo$age)
demo$duration <- scale(demo$duration, center = F)

tap <- merge(tap, demo, by.x = 'healthCode', by.y = 'healthCode', all.x=T, all.y = F)

resp <- colnames(tap)[c(17:59)]
tap[,resp] <- sapply(tap[,resp], as.numeric)
for (i in resp){
  tap <- tap[!is.na(tap[,i]),]
}
tap[,resp] <- scale(tap[,resp])

summary(lmer(as.formula(paste0('numberTaps','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = tap))
summary(lmer(as.formula(paste0('rangeTapInter','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = tap))
summary(lmer(as.formula(paste0('sdTapInter','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = tap))

# walking 
walk <- read.csv('./mpower_data/other/walking.csv', colClasses = c(createdOn='POSIXct'))
walk$healthCode <- as.character(walk$healthCode)

walking <- read.csv('./mpower_data/extracted features/walkFeatures.csv')
walking_r <- read.csv('./mpower_data/extracted features/walk_rest_Features.csv')

walk <- merge(walk, walking, by = 'accel_walking_outbound.json.items', all = F)
colnames(walking_r) <- paste0(colnames(walking_r),'_r')
walk <- merge(walk, walking_r, by.x = 'accel_walking_rest.json.items', by.y = 'accel_walking_rest.json.items_r',all = F)
rm(list=c('walking','walking_r'))

walk <- walk[order(walk$healthCode,walk$createdOn),]

w_sub_all <- split(walk[,7,drop=F], walk$healthCode)
w_sub_all <- w_sub_all[sapply(w_sub_all, function(x)NROW(x)>=20)]

usedid <- sapply(w_sub_all, function(x)x$healthCode[1])
demo <- demo[demo$healthCode %in% usedid,]

walk <- walk[walk$healthCode %in% demo$healthCode,]

demo$age <- scale(demo$age)
demo$duration <- scale(demo$duration, center = F)

walk <- merge(walk, demo, by.x = 'healthCode', by.y = 'healthCode', all.x=T, all.y = F)

resp <- colnames(walk)[c(18:130,132:150)]
walk[,resp] <- sapply(walk[,resp], as.numeric)
for (i in resp){
  walk <- walk[!is.na(walk[,i]),]
}
walk[,resp] <- scale(walk[,resp])

summary(lmer(as.formula(paste0('saAJ','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = walk))
summary(lmer(as.formula(paste0('dfaY','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = walk))
summary(lmer(as.formula(paste0('medianAA_r','~mtrt+professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode)')), data = walk))

walk$mtrt <- as.numeric(walk$medTimepoint=='Just after Parkinson medication (at your best)')
glmer(mtrt~professional.diagnosis+age+gender+med+smoked+duration+deep.brain.stimulation+(1 | healthCode), data = walk, family = binomial(), control = glmerControl(optimizer='bobyqa'))

# use similar codes in Figure 4 plot in simulation codes to extract residuals
boxplot(all_in_pearson_err,at=c(1,2,3,5,6,7,9,10,11,13,14,15), 
        col=c('#ffde7d','#00b8a9','#4a83b7'), medlwd = 1.1, xaxt='n', xlab='Response', 
        ylab='Mean squares of Pearson residuals')
legend('topright',legend = c('MRSS','GLMM','VAR'),fill=c('#ffde7d','#00b8a9','#4a83b7'),border = NA)
axis(1, at = c(2,6,10,14), labels=c('A_walking','median C1','sd C2','dfa Y'))
abline(v=4)
abline(v=8)
abline(v=12)

boxplot(all_out_pearso_err,at=c(1,2,3,5,6,7,9,10,11,13,14,15),
        col=c('#ffde7d','#00b8a9','#4a83b7'),  medlwd = 1.1, xaxt='n', xlab='Response', 
        ylab='Mean squares of Pearson residuals')
legend('topright',legend = c('MRSS','GLMM','VAR'),fill=c('#ffde7d','#00b8a9','#4a83b7'),border = NA)
axis(1, at = c(2,6,10,14), labels=c('A_tapping','range Tap Inter','number Taps','sd Tap Inter'))
abline(v=4)
abline(v=8)
abline(v=12)

# Figure 7 ----------------------------------------------------------------
#### sdc2
ind <- 229 # plot number 229 subject
md_n <- new_md[[1]][[ind]]
md_n$y <- md_n$y[(140-1),]
md_n$Z <- md_n$Z[,,4]
md_n$R <- md_n$R[,,(140-1)]
md_n$u <- md_n$u[(140-1),]
newdat <- as.data.frame(t(md_n$y))
newdat <- newdat[c(1,1,1,1,1,1,1,1),]
newdat[1:8,] <- NA
md_n <- SSModel(cbind(mtrt_w,rangeTapInter,numberTaps,sdTapInter,mtrt_w,medianAA_r,sdAJ,dfaY) ~  -1+SSMcustom(Z = md_n$Z, T = md_n$T, R = md_n$R, Q = md_n$Q, a1=md_n$a1, P1inf=md_n$P1inf, P1=md_n$P1), data=newdat, distribution = c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian'), u=md_n$u, tol = 1e-8)

kfs <- predict(md, newdata = md_n, filtered = F, interval = 'confidence')
kfs_org <- predict(md, filtered = F, interval = 'confidence')
y_org <- kfs_org$sdAJ
y_org <- y_org[1:(nrow(y_org)-1),]
y_pre <- kfs$sdAJ

y_true <- m_all[[ind]]$y[-1,]
y_true <- y_true[,4]
y_ind <- m_all[[ind]]$y[-1,1]
y_ind[is.na(y_ind)] <- 2

plot((2:140)[y_ind==0],-y_true[y_ind==0], ylab=TeX('$sd(C_2)$'), xlab = 'day', pch=20, col = 'grey60', xlim = c(0,150))
points((2:140)[y_ind==1],-y_true[y_ind==1],pch=1,  cex=0.8, lwd=1.5, col = 'grey40')
lines((2:140),-y_org[,1], col = 'red')
lines((2:140),-y_org[,2], lty=2, col = '#FF6969')
lines((2:140),-y_org[,3], lty=2, col = '#FF6969')
polygon(c(140:146,rev(140:146)), c(-c(y_org[139,2],y_pre[,2]), rev(-c(y_org[139,3],y_pre[,3]))), col='#FF9191', border=F)
lines(140:146,-c(y_org[139,1],y_pre[,1]), col = 'red')
lines(140:146,-c(y_org[139,2],y_pre[,2]), lty=2, col = 'red')
lines(140:146,-c(y_org[139,3],y_pre[,3]), lty=2, col = 'red')

legend('topleft', legend = c('Sample points','Sample points (just after trt)','Prediction curve'), col = c('grey60','grey40','red'),lty=c(NA,NA,1), pch=c(20,1,NA))


#### number of taps
ind <- 229 # plot number 229 subject
md_n <- new_md[[1]][[ind]]
md_n$y <- md_n$y[(140-1),]
md_n$Z <- md_n$Z[,,4]
md_n$R <- md_n$R[,,(140-1)]
md_n$u <- md_n$u[(140-1),]
newdat <- as.data.frame(t(md_n$y))
newdat <- newdat[c(1,1,1,1,1,1,1,1),]
newdat[1:8,] <- NA
md_n <- SSModel(cbind(mtrt_w,rangeTapInter,numberTaps,sdTapInter,mtrt_w,medianAA_r,sdAJ,dfaY) ~  -1+SSMcustom(Z = md_n$Z, T = md_n$T, R = md_n$R, Q = md_n$Q, a1=md_n$a1, P1inf=md_n$P1inf, P1=md_n$P1), data=newdat, distribution = c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian'), u=md_n$u, tol = 1e-8)

kfs <- predict(md, newdata = md_n, filtered = F, interval = 'confidence')
kfs_org <- predict(md, filtered = F, interval = 'confidence')
y_org <- kfs_org$numberTaps
y_org <- y_org[1:(nrow(y_org)-1),]
y_pre <- kfs$numberTaps

y_true <- m_all[[ind]]$y[-1,]
y_true <- y_true[,7]
y_ind <- m_all[[ind]]$y[-1,1]
y_ind[is.na(y_ind)] <- 2
plot((2:140)[y_ind==0],-y_true[y_ind==0], ylab='Number of Taps', xlab = 'day', pch=20, col = 'grey60', xlim = c(0,150), ylim=c(-3,2))
points((2:140)[y_ind==1],-y_true[y_ind==1],pch=1,  cex=0.8, lwd=1.5, col = 'grey40')
lines((2:140),-y_org[,1], col = 'red')
lines((2:140),-y_org[,2], lty=2, col = '#FF9191')
lines((2:140),-y_org[,3], lty=2, col = '#FF9191')
polygon(c(140:146,rev(140:146)), c(-c(y_org[139,2],y_pre[,2]), rev(-c(y_org[139,3],y_pre[,3]))), col='#FF9191', border=F)
lines(140:146,-c(y_org[139,1],y_pre[,1]), col = 'red')
lines(140:146,-c(y_org[139,2],y_pre[,2]), lty=2, col = 'red')
lines(140:146,-c(y_org[139,3],y_pre[,3]), lty=2, col = 'red')

legend('bottomleft', legend = c('Sample points','Sample points (just after trt)','Prediction curve'), col = c('grey60','grey40','red'),lty=c(NA,NA,1), pch=c(20,1,NA))


# Figure 8 ----------------------------------------------------------------
year <- 1998:2015
res_year <- NULL
for (yr in year){
  id <- demo$healthCode[demo$medication.start.year==yr]
  model_use <- new_md[[1]][[id]]
  for (i in length(model_use)){
    md_n <- model_use[i]
    nT <- nrow(md_n$y)
    md_n$Z <- md_n$Z[,,4]
    md_n$Z[,,4] <- 0
    md_n$R <- md_n$R[,,(nT-1)]
    md_n$u <- md_n$u[(nT-1),]
    newdat <- as.data.frame(t(md_n$y))
    newdat <- newdat[c(1,1,1,1,1,1,1,1),]
    newdat[1:8,] <- NA
    md_n <- SSModel(cbind(mtrt_w,rangeTapInter,numberTaps,sdTapInter,mtrt_w,medianAA_r,sdAJ,dfaY) ~  -1+SSMcustom(Z = md_n$Z, T = md_n$T, R = md_n$R, Q = md_n$Q, a1=md_n$a1, P1inf=md_n$P1inf, P1=md_n$P1), data=newdat, distribution = c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian'), u=md_n$u, tol = 1e-8)
    kfs <- predict(md, newdata = md_n, filtered = F, interval = 'confidence')
    y_pre <- kfs$sdAJ
    
    md_n <- model_use[i]
    nT <- nrow(md_n$y)
    md_n$Z <- md_n$Z[,,4]
    md_n$Z[,,4] <- 1
    md_n$R <- md_n$R[,,(nT-1)]
    md_n$u <- md_n$u[(nT-1),]
    newdat <- as.data.frame(t(md_n$y))
    newdat <- newdat[c(1,1,1,1,1,1,1,1),]
    newdat[1:8,] <- NA
    md_n <- SSModel(cbind(mtrt_w,rangeTapInter,numberTaps,sdTapInter,mtrt_w,medianAA_r,sdAJ,dfaY) ~  -1+SSMcustom(Z = md_n$Z, T = md_n$T, R = md_n$R, Q = md_n$Q, a1=md_n$a1, P1inf=md_n$P1inf, P1=md_n$P1), data=newdat, distribution = c('binomial','gaussian','gaussian','gaussian','binomial','gaussian','gaussian','gaussian'), u=md_n$u, tol = 1e-8)
    kfs <- predict(md, newdata = md_n, filtered = F, interval = 'confidence')
    y_pre1 <- kfs$sdAJ
    
    res_year <- rbind(res_year, c(yr, y_pre1-y_pre))
  }
}
colnames(res_year) <- c('yr','trt')
boxplot(trt~yr, data=res_year, xlab='First year taking Levodopa', ylab='Treatment effect in Sd(C2)', ylim=c(-3,3))
abline(h=0, lty=2, lwd=1)