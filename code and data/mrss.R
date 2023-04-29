
# mrss function is the main function of the model -------------------------
mrss <- function(ID, state, covariate, covariate_scale=TRUE, time_covariate, response, demo, subset=TRUE, data, day){
  # ID: the variable stands for subject ID; demo[,ID] must be equal to names(data)!
  # state: a data frame: state$ind, the name of state; state[,2:NCOL(state)], covariance matrix
  # covariate: time invariant covariate
  # covariate_scale: if true, all covariates will be scaled; if false, none will be scaled
  # time_covariate: time varying covariate; will NOT be scaled
  # day: variable indicates day, can be null (each row is one day); day must be 1,2,3.. cannot be 1,3,4... adding NA rows in this case
  # data: the dataset including y and time varying covariates; each element name is subject ID
  # demo: the dataset including time invariant covariates; one column should be subject ID
  # subset: a set of constraints. subset=TRUE means no constraint.
  # response: response$var, response variable; response$scale, whether can be scaled (1 yes, 0 no)
  # response$dis, response distribution; response$u, u vector (see help of KFAS package)
  # response$state_name, loading on variables (1 is constant; or a variable name in data, e.g., treatment indicator)
  
  # for example:
  # ID <- 'healthCode'
  # day <- 'time'
  # state <- data.frame(ind = c('b1_trt','b2_morning','v1','v2'),stringsAsFactors = F)
  # state <- cbind(state, matrix(c(rep(1,11),0,1,1,0,1),ncol=4))
  # 
  # subset <- 'professional.diagnosis==1 & med==1'
  # 
  # covariate <- c('age','gender','smoked','duration','deep.brain.stimulation')
  # time_covariate <- c('time')
  # 
  # response <- data.frame(var = c('mtrt_t','rangeTapInter','numberTaps',
  #                                'sdTapInter','mtrt_w','medianAA_r','sdAJ','dfaY'), stringsAsFactors = F)
  # response$scale <- c(0,1,1,1,0,1,1,1)
  # response$dis <- c('binomial','poisson','gaussian','gaussian','gaussian','gaussian','gaussian','gaussian')
  # response$u <- c(1,1,1,1,1,1,1,1)
  # response <- cbind(response, data.frame(b1_trt=c(0,'mtrt_t','mtrt_t','mtrt_t',0,'mtrt_w','mtrt_w','mtrt_w'),
  #                                        b2_morning=c('after_t','after_t','after_t','after_t',
  #                                                     'after_w','after_w','after_w','after_w'),
  #                                        v1=rep(1,8),v2=rep(1,8),stringsAsFactors = F))
  
  
  
  # scale -------------------------------------------------------------------
  if (covariate_scale) {
    demo[,covariate] <- scale(demo[,covariate])
  }
  
  can_scale_y <- response$var[response$dis=='gaussian' & response$scale==1]
  tmp <- do.call(rbind,lapply(data, function(x)(x[,can_scale_y])))
  mean_y <- colMeans(tmp, na.rm=T)
  sd_y <- apply(tmp, 2, sd, na.rm=TRUE)
  
  # matrix ------------------------------------------------------------------
  Tt <- matrix(0, ncol = NROW(state)+1, nrow = NROW(state)+1)
  diag(Tt) <- c(rep(1, NROW(state)), 1)
  Tt[1:NROW(state), NROW(state)+1] <- 0.5
  
  Qt <- as.matrix(state[,2:NCOL(state)])
  colnames(Qt) <- state$ind
  
  P1 <- 0
  P1inf <- matrix(0, ncol = NROW(state)+1, nrow = NROW(state)+1)
  diag(P1inf)[1:NROW(state)] <- 1
  
  if (all(response$dis=='gaussian')){
    Ht <- diag(0, NROW(response))
    diag(Ht) <- 1
  }else{
    u <- response$u
  }
  
  demo_sub <- eval(parse(text=paste('subset(demo,', subset, ')')))
  loopid <- sapply(demo_sub[,ID], function(x){which(x==demo[,ID])})
  m_all <- list() # all models
  dat_all <- list() # used data in the model
  for (i in loopid){
    id <- demo[i, ID]
    sub <- data[[id]]
    # scale Y -----------------------------------------------------------------
    sub[,can_scale_y] <- sweep(sub[,can_scale_y],2,mean_y)
    sub[,can_scale_y] <- sweep(sub[,can_scale_y],2,sd_y,'/')
    
    dat_all[[id]] <- sub[,c(response$var, time_covariate)]
    
    Zt <- array(0, dim = c(NROW(response), NROW(state)+1, NROW(sub)))
    for (j in unique(unlist(response[,5:NCOL(response)]))) {
      if (j=='0' || j==0){
        next
      }else if (j=='1' || j==1) {
        Zt[,1:NROW(state),][response[,5:NCOL(response)]==1 | response[,5:NCOL(response)]=='1'] <- 1
      }else {
        Zt[,1:NROW(state),][response[,5:NCOL(response)]==j] <- rep(sub[,j], each=sum(response[,5:NCOL(response)]==j))
      }
    }
    # beta <- rep(1, length(c(covariate,time_covariate)))
    if (is.null(time_covariate)){
      Zt[,NROW(state)+1,] <- sum(demo[i,covariate])
    }else{
      Zt[,NROW(state)+1,] <- sum(demo[i,covariate])+rep(rowSums(as.matrix(sub[,time_covariate,drop=F])), each = NROW(response))
    }
    
    Rt <- array(0, dim = c(NROW(state)+1,NROW(state),NROW(sub)))
    if (!is.null(day)){
      change <- c(diff(sub[,day]),0)
    }else{
      change <- rep(1, NROW(sub))
    }
    for (j in 1:NROW(sub)){
      if (change[j]==1){
        Rt[,,j] <- rbind(diag(1, nrow = NROW(state), ncol = NROW(state)),matrix(0, nrow = 1, ncol = NROW(state)))
      }
    }
    
    a1 <- matrix(c(rep(0, NROW(state)), 1), ncol = 1)
    if (all(response$dis=='gaussian')){
      m_all[[id]] <- SSModel(as.formula(paste0('cbind(',paste(response$var, collapse = ','),') ~ -1+SSMcustom(Z=Zt,T=Tt,R=Rt,Q=Qt,a1=a1,P1inf=P1inf,P1=P1, state_names = c(\"', paste(state$ind, collapse = '","'),'\",\"int\"))')), 
                             data = dat_all[[id]], distribution = 'gaussian', H=Ht)
    }else{
      m_all[[id]] <- SSModel(as.formula(paste0('cbind(',paste(response$var, collapse = ','),') ~ -1+SSMcustom(Z=Zt,T=Tt,R=Rt,Q=Qt,a1=a1,P1inf=P1inf,P1=P1, state_names = c(\"', paste(state$ind, collapse = '","'),'\",\"int\"))')), 
                             data = dat_all[[id]], distribution = response$dis, u=matrix(u, byrow=T, nrow = NROW(sub), ncol = NROW(response)))
    }
  }
  Z_para <- as.matrix(response[,5:NCOL(response)])
  Z_para[Z_para=='0' | Z_para==0] <- NA
  Z_para[!is.na(Z_para)] <- 1
  Z_para[is.na(Z_para)] <- 0
  Z_para <- apply(Z_para, 2, as.numeric)
  rownames(Z_para) <- response$var
  tmp <- setdiff(unique(c(apply(as.matrix(response[,5:NCOL(response)]),2,as.character))), c('0','1'))
  attr(m_all, 'response') <- response
  attr(m_all, 'state') <- state
  attr(m_all, 'covariate') <- covariate
  attr(m_all, 'time_covariate') <- time_covariate
  attr(m_all, 'demo') <- demo
  # attr(m_all, 'beta') <- rep(1, length(c(covariate,time_covariate)))
  # names(attr(m_all, 'beta')) <- c(covariate,time_covariate)
  attr(m_all, 'beta') <- matrix(1, ncol=length(c(covariate,time_covariate)), nrow=NROW(response))
  colnames(attr(m_all, 'beta')) <- c(covariate,time_covariate)
  rownames(attr(m_all, 'beta')) <- rownames(Z_para)
  attr(m_all, 'Z_para') <- Z_para
  attr(m_all, 'data') <- lapply(data[loopid], '[', unique(c(response$var, time_covariate, tmp, day)))
  attr(m_all, 'ID') <- ID
  attr(m_all, 'date') <- day
  return(m_all)
}

get_para <- function(model){
  # model: list(model1, model2, ....)
  #      model1, model2,... are different TYPE of models. In our real data example, modol1 is `PD patients taking Levodopa' model,
  #      modol2 is `PD patients not taking Levodopa' model, model3 is `healthy controls' model
  # get all parameters from a list of mrss models
  # you should NEVER reorder the rows of returned data frame!!!
  # see the help of KFAS to learn T,Q,u,H... under column 'source'
  # mode = 0, fixed values; mode = 1, same values across models; mode = 2, different values in difference TYPE models
  num_group <- length(model)
  all_para <- list()
  Zfix_name <- NULL
  ufix_name <- NULL
  for (i in 1:num_group){
    para <- list()
    md <- model[[i]][[1]]
    
    num_r <- NROW(md$Z[,,1]) # number of responses
    num_s <- NROW(md$T)-1 # number of states
    name_t <- colnames(md$T)[1:num_s]
    para$T <- data.frame(source='T',name=c(name_t, paste0(name_t,'_int')),
                         value=c(diag(md$T[,,1])[1:num_s], md$T[1:num_s,num_s+1,1]),stringsAsFactors = F)
    
    tmp <- md$Q[,,1]
    for (l1 in 1:num_s){
      for (l2 in 1:num_s){
        tmp[l1,l2] <- paste(name_t[l1], name_t[l2], sep = ':')
      }
    }
    para$Q <- data.frame(source='Q',name=tmp[upper.tri(tmp, diag=T)],
                         value=md$Q[,,1][upper.tri(md$Q[,,1], diag = T)],stringsAsFactors = F)
    
    if (all(md$u=='Omitted')){
      para$H <- data.frame(source='H',name=attr(model[[i]], 'response')$var,value=diag(md$H[,,1]),stringsAsFactors = F)
    }else{
      para$u <- data.frame(source='u',name=attr(model[[i]], 'response')$var,value=md$u[1,],stringsAsFactors = F)
      ufix_name <- c(ufix_name,response$var[response$dis=='gaussian'])
    }
    
    tmp <- data.frame(name=character(), value=numeric(),stringsAsFactors = F)
    tmp1 <- attr(model[[i]],'Z_para')
    Zfix <- suppressWarnings(as.numeric(as.matrix(attr(model[[i]], 'response')[,-(1:4)])))
    for (l1 in 1:num_s){
      for (l2 in 1:num_r){
        tmp <- rbind(tmp, data.frame(name=paste(rownames(tmp1)[l2], colnames(tmp1)[l1], sep = ':'), 
                                     value=tmp1[l2,l1],stringsAsFactors = F))
      }
    }
    para$Z <- tmp
    para$Z <- cbind(source='Z',para$Z)
    Zfix_name <- c(Zfix_name, na.omit(para$Z$name[Zfix==0]))
    
    # para$beta <- data.frame(source='Beta',name=names(attr(model[[i]],'beta')),
    # value=attr(model[[i]],'beta'),row.names = NULL,stringsAsFactors = F)
    
    tmp <- attr(model[[i]],'beta')
    para$beta <- data.frame(source='Beta',name=paste(rownames(tmp),rep(colnames(tmp),each=NROW(tmp)),sep=':'),
                            value=as.numeric(tmp),row.names = NULL,stringsAsFactors = F)
    
    all_para[[i]] <- do.call(rbind, para)
    position_t <- all_para[[i]]
    position_t$value <- 1:NROW(position_t)
    if (i == 1){
      final <- all_para[[i]]
      position <- position_t
    }else{
      final <- merge(final, all_para[[i]], by=c('source','name'), all = T, sort = T)
      position <- merge(position, position_t, by=c('source','name'), all = T, sort = T)
    }
    colnames(final)[colnames(final)=='value'] <- paste0('value', i)
    colnames(position)[colnames(position)=='value'] <- paste0('value', i)
  }
  # mode = 0, fixed values; mode = 1, same values across models; mode = 2, different values
  final$mode <- 0
  final$mode[final$source %in% c('H','Q','Z')] <- 1
  final$mode[final$source %in% c('u') & final$name %in% ufix_name] <- 1
  final$mode[final$source %in% c('T','Beta')] <- 2
  final$mode[final$source %in% c('Z') & final$name %in% Zfix_name] <- 0
  attr(final, 'position_info') <- position
  final$ind <- seq(1, (NCOL(final)-3)*NROW(final), NCOL(final)-3)
  attr(final, 'position_info')$ind <- seq(1, (NCOL(final)-4)*NROW(final), NCOL(final)-4)
  return(final)
}



para_vec <- function(para){
  tmp <- split(para, para$mode)
  
  final1 <- data.frame(name = paste(tmp$`1`$source, tmp$`1`$name, sep = '_'),
                       x = rowMeans(tmp$`1`[,3:(NCOL(para)-2),drop=F],na.rm = T), 
                       ind = tmp$`1`$ind,stringsAsFactors = F)
  if (!is.null(tmp$`2`)){
    final2 <- data.frame(name = paste0(rep(paste(tmp$`2`$source, tmp$`2`$name, sep = '_'), each = NCOL(para)-4),'(',1:(NCOL(para)-4),')'),
                         x = c(t(tmp$`2`[,3:(NCOL(para)-2),drop=F])),
                         ind = c(sapply(tmp$`2`$ind, function(x)seq(x, x+NCOL(para)-5))),stringsAsFactors = F)
    final2 <- final2[!is.na(final2$x),,drop=F]
    final <- rbind(final1, final2)
    final <- final[order(final$ind),]
    return(final)
  }else{
    return(final1)
  }
}

vec_para <- function(para, para_vec){
  tmp <- para
  
  tmp1 <- data.frame(ind = 1:(NROW(para)*(NCOL(para)-4)))
  para_vec <- merge(tmp1, para_vec,  by = 'ind', all = T)
  ok <- which(!is.na(para_vec$x))
  gaps <- diff(c(ok, NROW(para_vec) + 1L))
  para_vec$x <- rep(para_vec$x[ok], gaps)
  para[,3:(NCOL(para)-2)] <- matrix(para_vec$x, byrow = T, nrow = NROW(para))
  para[,3:(NCOL(para)-2)][is.na(tmp[,3:(NCOL(para)-2)])] <- NA
  para[tmp$mode==0,3:(NCOL(para)-2)] <- tmp[tmp$mode==0,3:(NCOL(para)-2)]
  return(para)
}

plug_para <- function(model, paralist){
  
  num_group <- length(model)
  
  for (i in 1:num_group){
    para <- paralist[,c(1:2, i+2)]
    pos <- attr(paralist, 'position_info')[,c(1:2, i+2)]
    para <- para[!is.na(pos[,3]),]
    pos <- pos[!is.na(pos[,3]),]
    para <- para[order(pos[,3]),]
    
    num_r <- NROW(model[[i]][[1]]$Z[,,1]) # number of responses
    num_s <- NROW(model[[i]][[1]]$T)-1 # number of states
    ustat <- all(model[[i]][[1]]$u=='Omitted')
    for (j in 1:length(model[[i]])){
      
      diag(model[[i]][[j]]$T[,,1])[1:num_s] <- para[para$source=='T',3][1:num_s]
      model[[i]][[j]]$T[1:num_s,num_s+1,1] <- para[para$source=='T',3][(num_s+1):(2*num_s)]
      
      model[[i]][[j]]$Q[,,1][upper.tri(model[[i]][[j]]$Q[,,1], diag = T)] <- para[para$source=='Q',3]
      model[[i]][[j]]$Q[,,1][lower.tri(model[[i]][[j]]$Q[,,1])] <- t(model[[i]][[j]]$Q[,,1])[lower.tri(model[[i]][[j]]$Q[,,1])]
      
      if (ustat){
        diag(model[[i]][[j]]$H[,,1]) <- para[para$source=='H',3]
      }else{
        model[[i]][[j]]$u <- matrix(para[para$source=='u',3], byrow = T, nrow=NROW(model[[i]][[j]]$u), ncol = num_r)
      }
      
      # Z -----------------------------------------------------------------------
      model[[i]][[j]]$Z[,1:num_s,] <- para[para$source=='Z',3]
      tmp <- attr(model[[i]],'response')[,-c(1:4)]
      for (jj in unique(unlist(tmp))) {
        if (jj!='0' && jj!=0 && jj!='1' && jj!=1){
          model[[i]][[j]]$Z[,1:num_s,][tmp==jj] <- rep(attr(model[[i]],'data')[[names(model[[i]][j])]][,jj], each=sum(tmp==jj))*model[[i]][[j]]$Z[,1:num_s,][tmp==jj]
        }
      }
      
      # beta <- rep(1, length(c(covariate,time_covariate)))
      demo <- attr(model[[i]],'demo')
      sub <- attr(model[[i]],'data')[[names(model[[i]][j])]]
      int_no_time <- matrix(para[para$source=='Beta',3],nrow=num_r)[,1:length(attr(model[[i]],'covariate')),drop=F] %*% t(demo[demo[,attr(model[[i]],'ID')]==names(model[[i]][j]),attr(model[[i]],'covariate'),drop=F])
      if (is.null(attr(model[[i]],'time_covariate'))){
        model[[i]][[j]]$Z[,num_s+1,] <- int_no_time
      }else{
        model[[i]][[j]]$Z[,num_s+1,] <- as.numeric(int_no_time)+matrix(para[para$source=='Beta',3],nrow=num_r)[,-c(1:length(attr(model[[i]],'covariate'))),drop=F] %*% t(as.matrix(sub[,attr(model[[i]],'time_covariate'),drop=F]))
      }
    }
    
    if (!ustat){
      attr(model[[i]],'response')$u <- para[para$source=='u',3]
    }
    attr(model[[i]],'Z_para')[1:num_r,1:num_s] <- para[para$source=='Z',3]
    attr(model[[i]],'beta')[1:length(attr(model[[i]],'beta'))] <- para[para$source=='Beta',3]
  }
  return(model)
}

callike <- function(m_all, trace=NULL){
  m_all <- unlist(m_all, recursive = F)
  setl <- sapply(m_all, function(x)sum(!is.na(rowMeans(x$y))))
  l <- NULL
  for (i in 1:length(m_all)){
    l <- c(l, logLik(m_all[[i]], maxiter=100, H_tol=1e10))
  }
  l <- l/setl
  
  l <- l[!is.infinite(l)]
  l <- l[!is.na(l)]
  l <- l[l>quantile(l,0.05) & l<quantile(l,0.95)]
  if (!is.null(trace)){
    write.table(mean(l), paste0(trace,'.txt'), append = T, quote = F, row.names = F, col.names = F)
  }
  return(mean(l))
}