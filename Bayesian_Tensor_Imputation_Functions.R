
##########################################################################################
####      
####             Bayesian multiple tensor imputation method with i.i.d error 
####      Generate the posterior samples of the underlying structure X1 X2 X3... of the decomposition with standard error \sigma
####      Users can reconstruct the imputated tensor with [[X1,X2,X3...]]+N(0,\sigma^2)
####      ts: observed tensor with missing elements
####      num_com: rank of the tensor decomposition (number of rank 1 tensor components)
####      max: max MCMC iteration 
####      burn.in: burnin value for the posterior samples
####      initial: particular initial value for the Bayes sampling
####      cpem_as_initial: whether use frequentist cp decomposition result as initial value
####      thin: thin value for the posterior samples, thin = k means only keep the 1/k posterior samples as result
####      
##########################################################################################

cpbayeimp_Jef_check_effi=function(ts,num_com,max =1000,burn.in=500,initial=NA,cpem_as_initial=FALSE,thin=1){
  
  
  
  # Missing Value imputed as 0 in initial
  dims=dim(ts)
  N=length(dims)
  M=prod(dims)
  tts=1:M
  
  for(e in 1:M){
    if(ts[e]=="NA"){
      #if(is.na(ts[e])){
      tts[e]=0
      ts[e]="0"
    }else{
      tts[e]=1
    }
    
  }
  
  ts=as.numeric(ts)
  
  ts=array(ts,dim = dims)
  
  u_list=vector("list",N)
  unfolded=vector("list",N)
  
  ts_norm=fnormimp(ts)
  
  for(m in 1:N){
    #unfolded[[m]]=unfoldimp(ts,mode = m)
    u_list[[m]]=matrix(rnorm(dims[m]*num_com), nrow = dims[m],ncol = num_com)
  }
  
  if(!is.na(initial)){
    u_list=initial
  }
  est=ts
  now=ts
  
  if(cpem_as_initial){
    a1=cpimp(ts,num_com)
    u_list=a1$U
    est=a1$est
  }
  
  inxxx=1
  
  
  #print(a1$norm_percent)
  
  w=sum(((now*tts-est*tts))^2)
  mm=sum(tts)
  draw=rgamma(1,shape = mm/2,scale = 2/w)
  sigma_est=1/draw
  
  
  
  for(ii in 1:M){
    if(tts[ii]==0){
      now[ii]=rnorm(1,mean=est[ii],sd=sqrt(sigma_est))
      
    }
    else{
      now[ii]=ts[ii]
      
    }
  }
  
  #print(solve(matrix(data = c(1,1,1,1),nrow = 2)))
  
  
  curr_it=1
  
  total=round((max-burn.in)/thin)
  #samp=array(dim = c(max,dims))
  
  #resu=array(dim = c(total,dims))
  resu1=array(0,dim=c(total,dims[1],num_com))
  resu2=array(0,dim=c(total,dims[2],num_com))
  resu3=array(0,dim=c(total,dims[3],num_com))
  sigma_err=1:total
  
  
  #resu_err=resu
  norm_vec=function(x){
    norm(as.matrix(x))
  }
  sigma_est<<-1
  while (curr_it<=max) {
    
    #print(curr_it)
    
    
    #print(TSigma_est)
    for(m in 1:N){
      Amat=outerimp_mat(u_list[-m])
      solv=solve(outerimp_updated(u_list[-m]))
      now_mat=cbind(apply(now,m,as.vector))
      beta_est=solv%*%t(Amat)%*%now_mat
      sigtest=(sigma_est)*solv
      rand_part=rbind(mvrnorm(n=dims[m],mu=rep(0,num_com),Sigma = sigtest))
      
      u_list[[m]]=t(beta_est)+rand_part
      #print(Amat)
      
    }
    
    totalprod=1
    partprod=rep(1,num_com)
    
    
    for(rr in 1:num_com){
      totalprod=1
      #print(num_com)
      for(kk in 1:N){
        totalprod=totalprod*((mean((u_list[[kk]][,rr])^2))^(1/(N*2)))
        #print((sum((u_list[[m]][,kk])^2))^(1/(num_com*2)))
      }
      #print(totalprod)
      for(kk in 1:N){
        u_list[[kk]][,rr]=u_list[[kk]][,rr]*(totalprod)/((mean((u_list[[kk]][,rr])^2))^(1/2))
        #print(u_list[[kk]][,rr])
      }
    }
    
    #U_check<<-u_list
    
    
    
    #print(1)
    
    
    # Imputation
    est=outerimp(u_list)
    
    
    
    
    w=sum((now-est)^2)
    draw=rchisq(1,M)
    sigma_est=(w)/draw
    
    if(curr_it>burn.in){
      if((curr_it-burn.in)%%thin==0){
        resu1[inxxx,,]=u_list[[1]] 
        resu2[inxxx,,]=u_list[[2]] 
        resu3[inxxx,,]=u_list[[3]] 
        sigma_err[inxxx]=sigma_est
        inxxx=inxxx+1
      }
      
    }
    
    
    
    
    
    for(ii in 1:M){
      if(tts[ii]==0){
        now[ii]=rnorm(1,mean=est[ii],sd=sqrt(sigma_est))
        
      }
      else{
        now[ii]=ts[ii]
        
      }
    }
    
    curr_it=curr_it+1
  }
  
  invisible(list(resu1=resu1,resu2=resu2,resu3=resu3,M=max,sigma_err=sigma_err))
}




##########################################################################################
####      
####             Bayesian multiple tensor imputation method with seperable covariance structure 
####      Generate the posterior samples of the tensor with missing value imputed, with covariance structure \Sigma_1, \Sigma_2, \Sigma_3,...
####      ts: observed tensor with missing elements
####      num_com: rank of the tensor decomposition (number of rank 1 tensor components)
####      max: max MCMC iteration 
####      burn.in: burnin value for the posterior samples
####      initial: particular initial value for the Bayes sampling
####      cpem_as_initial: whether use frequentist cp decomposition result as initial value
####      thin: thin value for the posterior samples, thin = k means only keep the 1/k posterior samples as result
####      identi: indicator of which dimension is assumed to be independent
####      fiber_mis: indicator of which dimension is assumed to be fiber-wise missing; 0 indicate no fiberwise missing and the elements are missing randomly
####      show_progress: whether show the progress bar
####      
##########################################################################################

cpbaye_covariance_appli_mix2_fiber_5=function(ts,num_com,max =1000,burn.in=500,initial=NA,cpem_as_initial=FALSE,thin=1,identi=2,fiber_mis=0,show_progress=FALSE){
  
  # Missing Value imputed as 0 in initial
  time_start=proc.time()
  
  
  
  dims=dim(ts)
  N=length(dims)
  M=prod(dims)
  tts=1:M
  sigma_identi=identi
  tausq=1
  
  norm_vec=function(x){
    norm(as.matrix(x))
  }
  
  
  for(e in 1:M){
    if(ts[e]=="NA"){
      tts[e]=0
      ts[e]="0"
    }else{
      tts[e]=1
    }
    
  }
  
  ts=as.numeric(ts)
  
  ts=array(ts,dim = dims)
  
  u_list=vector("list",N)
  unfolded=vector("list",N)
  sigma_est=vector("list",N)
  
  for(m in 1:N){
    sigma_est[[m]]=diag(dims[m])
  }
  
  
  ts_norm=fnormimp(ts)
  
  for(m in 1:N){
    u_list[[m]]=matrix(rnorm(dims[m]*num_com), nrow = dims[m],ncol = num_com)
  }
  
  if(!is.na(initial)){
    u_list=initial
  }
  est=ts
  now=ts
  
  if(cpem_as_initial){
    a1=cpimp(ts,num_com)
    u_list=a1$U
    est=a1$est
  }
  
  inxxx=1
  
  curr_it=1
  
  total=round((max-burn.in)/thin)
  
  
  resu=array(dim = c(total,dims))
  resu_under=array(dim = c(total,dims))
  sigma1=array(dim = c(total,dims[1],dims[1]))
  sigma2=array(dim = c(total,dims[2],dims[2]))
  sigma3=array(dim = c(total,dims[3],dims[3]))
  est_tau=1:total
  log_likelihood=1:total
  tausq=1
  
  if(show_progress){
    pb = txtProgressBar(min = 0, max = max, initial = 0,style = 3) 
  }
  timep1=proc.time()
  
  
  m_inx=array((tts==0),dim = dims)
  o_inx=array((tts==1),dim = dims)
  m_inx_fiber=array()
  comp_dim=0
  
  # identify the missing index for fiber_wise missing
  
  if((fiber_mis!=0)&(fiber_mis!=identi)){
    comp_dim=c(1:3)[-c(fiber_mis,identi)]
    
    m_inx_fiber=apply(m_inx,c(identi,comp_dim),FUN=function(x){
      x[1]
    })
    
  }
  
  
  timep2=proc.time()
  
  totaltime1=timep2-timep1
  totaltime2=0
  totaltime3=0
  totaltime4=0
  totaltime5=0
  totaltime6=0
  totaltime7=0
  totaltime8=0
  
  while (curr_it<=max) {
    
    timep3=proc.time()
    
    if(show_progress){
      setTxtProgressBar(pb,curr_it)
    }
    
    for(m in 1:N){
      
      timep4=proc.time()
      
      
      if((m %in% sigma_identi)){
        
        sigma_est_inv=lapply(sigma_est,FUN=solve)
        sigma_est_inv_sqrt=lapply(sigma_est_inv,FUN=sqrtm)
        
        
        #sqrt_inverse_sigma_minus_k=kroneckerlist(rev(sigma_est_inv_sqrt[-m]))
        
        u_minus_k=krlistimp(rev(u_list[-m]))
        
        #y_tilde=unfoldimp(now,mode = m) %*% sqrt_inverse_sigma_minus_k
        y_tilde=kronecker_prod_effi(unfoldimp(now,mode = m),rev(sigma_est_inv_sqrt[-m]))
        
        
        
        #x_tilde=t(u_minus_k)  %*% sqrt_inverse_sigma_minus_k
        x_tilde=kronecker_prod_effi(t(u_minus_k),rev(sigma_est_inv_sqrt[-m]))
        
        x_tilde=t(x_tilde)
        
        solv=solve(t(x_tilde)%*%x_tilde)
        sigtest=(tausq)*solv
        
        for(ind in 1:dims[m]){
          
          Yi=array(y_tilde[ind,],dim = c(length(y_tilde[ind,]),1))
          
          beta_est=solv%*%t(x_tilde)%*%Yi
          
          u_ir=mvrnorm(n=1,mu=beta_est,Sigma = sigtest)
          u_list[[m]][ind,]=u_ir
          
        }
        #alpha_p=prod(dims)/2+0.001
        #beta_p=sum((u_list[[m]]%*%t(x_tilde)-y_tilde)^2)/2+0.001
        #tausq=1/(rgamma(1,shape = alpha_p,rate = beta_p))
        
        tausq=1
        sigma_est[[m]]=tausq*diag(dims[m])
        
        
        timep5=proc.time()
        
        totaltime2=totaltime2+timep5-timep4
        
      }else{
        
        sigma_est_inv=lapply(sigma_est,FUN=solve)
        sigma_est_inv_sqrt=lapply(sigma_est_inv,FUN=sqrtm)
        
        #sqrt_inverse_sigma_minus_k=kroneckerlist(rev(sigma_est_inv_sqrt[-m]))
        
        
        u_minus_k=krlistimp(rev(u_list[-m]))
        
        #y_tilde=unfoldimp(now,mode = m) %*% sqrt_inverse_sigma_minus_k
        y_tilde=kronecker_prod_effi(unfoldimp(now,mode = m),rev(sigma_est_inv_sqrt[-m]))
        
        y_tilde=t(y_tilde)
        
        #x_tilde=t(u_minus_k)  %*% sqrt_inverse_sigma_minus_k
        x_tilde=kronecker_prod_effi(t(u_minus_k),rev(sigma_est_inv_sqrt[-m]))
        
        x_tilde=t(x_tilde)
        
        
        lambda0=diag(num_com)*0.0001
        
        betan=solve(t(x_tilde)%*%(x_tilde)+lambda0)%*%t(x_tilde)%*%(y_tilde)
        
        Sn= diag(dims[m])+t(y_tilde-x_tilde%*%betan)%*%(y_tilde-x_tilde%*%betan)+t(betan)%*%lambda0%*%betan
        
        
        sigma_est[[m]]=riwish(v=(dims[m]+2+prod(dims[-m])),S=Sn)
        
        u_list[[m]]=t(matrix_normal(M=betan,U=solve(t(x_tilde)%*%(x_tilde)+lambda0),V=sigma_est[[m]]))
        
        timep6=proc.time()
        
        totaltime3=totaltime3+timep6-timep4
      }
      
      
      
      
      
      
    }
    
    timep7=proc.time()
    
    
    #sigma_matrix=kronecker_list(rev(sigma_est[-sigma_identi]))
    
    # whether the data is fiber-wise missing
    
    if((fiber_mis==0)){
      
      sigma_matrix=kronecker_list(rev(sigma_est[-sigma_identi]))
      
      
      for(j in 1:dims[sigma_identi]){
        
        txt=paste("m_j=m_inx[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
        eval(parse(text = txt))
        
        m_j=array(m_j,dim = length(m_j))
        
        txt=paste("o_j=o_inx[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
        eval(parse(text = txt))
        
        #o_j=array(o_j,dim = length(o_j))
        o_j=array(o_j,dim = length(o_j))
        
        # whether the elements are all observed (if so, do not impute)
        
        if((sum(o_j)!=prod(dims[-sigma_identi]))){
          
          if(sum(o_j)==0){
            
            txt=paste("full_mu=outerimp(u_list)[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            condi_mu=full_mu
            
            
            #imputed=mvtnorm::rmvnorm(n=1,mean = condi_mu,sigma =condi_sigma,method = "chol")
            
            chol_matrix=rev(lapply(sigma_est[-sigma_identi],chol))
            
            #B_mat=kronecker_list(chol_matrix)
            
            U_normal=matrix(rnorm(prod(dims[-sigma_identi]),mean = 0,sd=1),nrow = 1)
            #imputed=condi_mu+array(U_normal%*%B_mat,dim = dim(condi_mu))
            imputed=condi_mu+array(kronecker_prod_effi(U_normal,chol_matrix),dim = dim(condi_mu))
            
            
            txt=paste("now[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]=imputed",sep = "")
            eval(parse(text = txt))
            
          }else{
            
            txt=paste("full_mu=outerimp(u_list)[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            mu1=full_mu[m_j]
            mu2=full_mu[o_j]
            
            txt=paste("ts_slide=ts[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            
            observed=ts_slide[o_j]
            
            #sigma11=sigma_matrix[m_j,m_j]
            #sigma21=sigma_matrix[o_j,m_j]
            #sigma12=sigma_matrix[m_j,o_j]
            #sigma22=sigma_matrix[o_j,o_j]
            
            
            
            sigma11=sigma_matrix[m_j,m_j]
            sigma21=sigma_matrix[o_j,m_j]
            
            #sigma22=solve(kronecker_list(rev(sigma_est[-sigma_identi]))[o_j,o_j])
            
            invsigma22=solve(sigma_matrix[o_j,o_j])
            
            condi_mu=(mu1)+array(t(sigma21)%*%invsigma22%*%(observed-mu2),dim = length(mu1))
            condi_sigma=sigma11-t(sigma21)%*%invsigma22%*%sigma21
            
            condi_sigma=(condi_sigma+t(condi_sigma))/2
            condi_sigma_chol=chol(condi_sigma)
            u_norm=rnorm(dim(condi_sigma_chol)[1],mean = 0,sd=1)
            #imputed=mvtnorm::rmvnorm(n=1,mean = condi_mu,sigma =condi_sigma,method = "chol")
            
            imputed=array(condi_mu,dim = dim(condi_sigma_chol)[1])+array(condi_sigma_chol%*%u_norm,dim = dim(condi_sigma_chol)[1])
            
            
            ts_array=array(ts_slide,dim = length(ts_slide))
            ts_array[m_j]=imputed
            
            ts_slide=array(ts_array,dim = dims[-sigma_identi])
            
            
            txt=paste("now[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]=ts_slide",sep = "")
            eval(parse(text = txt))
          }
        }
      }
      
      gc() 
    }else{
      for(j in 1:dims[sigma_identi]){
        
        txt=paste("m_j=m_inx[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
        eval(parse(text = txt))
        
        m_j=array(m_j,dim = length(m_j))
        
        txt=paste("o_j=o_inx[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
        eval(parse(text = txt))
        
        o_j=array(o_j,dim = length(o_j))
        
        if(sum(m_inx_fiber[j,])!=0){
          
          
          if(sum(m_inx_fiber[j,])==dims[comp_dim]){
            
            
            timep8=proc.time()
            
            txt=paste("full_mu=outerimp(u_list)[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            condi_mu=full_mu
            
            
            chol_matrix=rev(lapply(sigma_est[-sigma_identi],chol))
            
            #B_mat=kronecker_list(chol_matrix)
            
            
            U_normal=matrix(rnorm(prod(dims[-sigma_identi]),mean = 0,sd=1),nrow = 1)
            #imputed=condi_mu+array(U_normal%*%B_mat,dim = dim(condi_mu))
            imputed=condi_mu+array(kronecker_prod_effi(U_normal,chol_matrix),dim = dim(condi_mu))
            
            txt=paste("now[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]=imputed",sep = "")
            eval(parse(text = txt))
            
            timep9=proc.time()
            
            totaltime4=totaltime4+timep8-timep9
            
          }else{
            
            
            ##########################################################################
            
            timep10=proc.time()
            
            
            mm_j=ifelse(m_inx_fiber[j,]==1,TRUE,FALSE)
            oo_j=ifelse(m_inx_fiber[j,]==0,TRUE,FALSE)
            sig11=matrix(sigma_est[[comp_dim]][mm_j,mm_j],nrow = sum(mm_j))
            sig12=matrix(sigma_est[[comp_dim]][mm_j,oo_j],nrow = sum(mm_j))
            sig22=matrix(sigma_est[[comp_dim]][oo_j,oo_j],nrow = sum(oo_j))
            
            
            
            txt=paste("full_mu=outerimp(u_list)[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            mu1=full_mu[,mm_j]
            if(sum(mm_j)==1){
              mu1=matrix(mu1, ncol = sum(mm_j))
            }
            
            mu2=full_mu[,oo_j]
            if(sum(oo_j)==1){
              mu2=matrix(mu2, ncol = sum(oo_j))
            }
            
            
            txt=paste("ts_slide=ts[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]",sep = "")
            eval(parse(text = txt))
            observed=ts_slide[,oo_j]
            
            if(sum(oo_j)>0){
              observed=matrix(observed, ncol = sum(oo_j))
              
            }
            
            
            a12a22inv=sig12%*%solve(sig22)
            
            condi_mean=mu1
            
            
            
            timep11=proc.time()
            
            
            
            for(jj in 1:dims[fiber_mis]){
              condi_mean[jj,]=mu1[jj,]+a12a22inv%*%(observed[jj,]-mu2[jj,])
            }
            
            condi_sigma_chol=kronecker(chol(sig11-sig12%*%solve(sig22)%*%t(sig12)),chol(sigma_est[[fiber_mis]]))
            
            
            timep12=proc.time()
            
            
            u_norm=rnorm(dim(condi_sigma_chol)[1],mean = 0,sd=1)
            
            imputed=array(condi_mean,dim = dim(condi_sigma_chol)[1])+array(condi_sigma_chol%*%u_norm,dim = dim(condi_sigma_chol)[1])
            
            
            timep13=proc.time()
            
            ts_slide[,mm_j]=imputed
            
            
            
            txt=paste("now[",paste(rep(",",(sigma_identi-1)),collapse = ""),j,paste(rep(",",(N-sigma_identi)),collapse = ""),"]=ts_slide",sep = "")
            eval(parse(text = txt))
            
            
            timep14=proc.time()
            
            totaltime5=totaltime5+timep11-timep10
            totaltime6=totaltime6+timep12-timep11
            
            totaltime7=totaltime7+timep13-timep12
            
            totaltime8=totaltime8+timep14-timep13
            
            ##########################################################################
            
            
          }
          
        }
        
      }
    }
    
    
    
    
    
    
    
    if(curr_it>burn.in){
      if((curr_it-burn.in)%%thin==0){
        resu[inxxx,,,]=now
        
        resu_under[inxxx,,,]=outerimp(u_list)
        sigma1[inxxx,,]=sigma_est[[1]]
        sigma2[inxxx,,]=sigma_est[[2]]
        sigma3[inxxx,,]=sigma_est[[3]]
        est_tau[inxxx]=tausq
        inxxx=inxxx+1
        
      }
    }
    
    curr_it=curr_it+1
  }
  time_end=proc.time()
  
  
  invisible(list(resu=resu,resu_under=resu_under,M=max,
                 sigma1=sigma1,sigma2=sigma2,sigma3=sigma3,
                 tau=est_tau,
                 time=time_start-time_end))
  
}






##########################################################################################
####      
####             The following are some base functions that will be called in our Bayesian imputation methods
####  
##########################################################################################


unfoldimp=function(ts,mode=1,dims=NULL){
  dims=dim(ts)
  M=prod(dims)
  N=length(dims)
  rs=mode
  cs=(1:N)[-mode]
  perm=c(rs,cs)
  newm=c(prod(dims[rs]),prod(dims[cs]))
  mat=ts
  mat=aperm(mat,perm)
  dim(mat)=newm
  mat
}

foldimp=function(ts,mode=1,dims=NULL){
  M=prod(dims)
  N=length(dims)
  rs=mode
  cs=(1:N)[-mode]
  matm=dim(ts)
  iperm=match(1:N,c(rs,cs))
  aperm(array(ts,dim = c(dims[rs],dims[cs])),iperm)
}


fnormimp=function(ts){
  sum(ts^2)
}

klistimp=function(x){
  res=x[[1]]
  for(i in 2:length(x)){
    res=kronecker(res,x[[i]])
  }
  res
}
krprodimp=function(x,y){
  res=matrix(0,nrow = dim(x)[1]*dim(y)[1],ncol = dim(x)[2])
  for(j in 1:dim(x)[2]){
    res[,j]=kronecker(x[,j],y[,j])
  }
  res
}
krlistimp=function(x){
  ncols=unlist(lapply(x, ncol))
  ncols=ncols[1]
  nrows=unlist(lapply(x, nrow))
  res=matrix(0,nrow = prod(nrows),ncol = ncols)
  for(j in 1:ncols){
    Lj=lapply(x,function(x) x[,j])
    res[,j]=klistimp(Lj)
    
  }
  res
}
hdlistimp=function(x){
  res=x[[1]]
  for(i in 2:length(x)){
    res=res*x[[i]]
  }
  res
}

superdiagtsimp=function(num_modes,len,elements=1L){
  modes=rep(len,num_modes)
  arr=array(0,dim = modes)
  if(length(elements)==1){
    elements=rep(elements,len)
  }
  for(i in 1:len){
    txt=paste("arr[",paste(rep("i",num_modes),collapse = ","),"] <-",elements[i],sep = "")
    eval(parse(text = txt))
  }
  arr
  
}

ttlimp=function(x,listma,ms=NULL){
  num_ma=length(listma)
  mat_nrow=vector("list",num_ma)
  mat_ncol=vector("list",num_ma)
  for(i in 1:num_ma){
    mat=listma[[i]]
    m=ms[i]
    mat_dims=dim(mat)
    modein=dim(x)
    stopifnot(modein[m]==mat_dims[2])
    modeout=modein
    modeout[m]=mat_dims[1]
    x_m=unfoldimp(x,mode = m)
    retarr=mat%*%x_m
    x=foldimp(retarr,mode = m,dims = modeout)
    
  }
  x
}

outerimp=function(ulist){
  num_com=dim(ulist[[1]])[2]
  N=length(ulist)
  est=0
  for(r in 1:num_com){
    oute=outer(ulist[[1]][,r],ulist[[2]][,r])
    for(k in 3:N){
      oute=outer(oute,ulist[[k]][,r])
    }
    est=est+oute
    
  }
  est
}

outerimp_mat=function(ulist){
  num_com=dim(ulist[[1]])[2]
  N=length(ulist)
  len=1:N
  for(i in 1:N){
    len[i]=dim(ulist[[i]])[1]
  }
  M=prod(len)
  est=matrix(0,nrow = M,ncol = num_com)
  
  for(r in 1:num_com){
    oute=outer(ulist[[1]][,r],ulist[[2]][,r])
    if(N>=3){
      for(k in 3:N){
        oute=outer(oute,ulist[[k]][,r])
      }
    }
    #print(as.vector(oute))
    est[,r]=t(as.vector(oute))
    
  }
  est
}
outerimp_updated=function(ulist){
  
  est=(t(ulist[[1]])%*%ulist[[1]])*(t(ulist[[2]])%*%ulist[[2]])
  est
}

CCofBayes_update=function(u11,u12,u13,u21,u22,u23,sd1,sd2,x,orx,mis){
  res_mean1=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  res_mean2=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  res_cover11=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  res_cover12=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  res_cover21=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  res_cover22=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  srf=array(0,dim = c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  n=dim(u11)[1]
  store1=1:n
  store2=1:n
  for(i in 1:dim(u11)[2]){
    for(j in 1:dim(u12)[2]){
      for(k in 1:dim(u13)[2]){
        for(ii in 1:n){
          store1[ii]=0
          store2[ii]=0
          for(r in 1:dim(u11)[3]){
            store1[ii]=store1[ii]+u11[ii,i,r]*u12[ii,j,r]*u13[ii,k,r]
            store2[ii]=store2[ii]+u21[ii,i,r]*u22[ii,j,r]*u23[ii,k,r]
          }
        }
        res_mean1[i,j,k]=mean(store1[ii])
        res_mean2[i,j,k]=mean(store2[ii])
        store_obs1=rnorm(n,mean = store1,sd=sd1)
        store_obs2=rnorm(n,mean = store2,sd=sd2)
        
        
        res_cover11[i,j,k]=ifelse((quantile(store1,prob=0.025)<orx[i,j,k])&(quantile(store1,prob=0.975)>orx[i,j,k]),1,0)
        res_cover12[i,j,k]=ifelse((quantile(store_obs1,prob=0.025)<x[i,j,k])&(quantile(store_obs1,prob=0.975)>x[i,j,k]),1,0)
        
        res_cover21[i,j,k]=ifelse((quantile(store2,prob=0.025)<orx[i,j,k])&(quantile(store2,prob=0.975)>orx[i,j,k]),1,0)
        res_cover22[i,j,k]=ifelse((quantile(store_obs2,prob=0.025)<x[i,j,k])&(quantile(store_obs2,prob=0.975)>x[i,j,k]),1,0)
        
        srf[i,j,k]=sqrt((2*sd(c(store1,store2))^2)/(sd(store1)^2+sd(store2)^2))
        
      }
    }
  }
  coverage11=sum(res_cover11*mis)/sum(mis)
  coverage12=sum(res_cover12*mis)/sum(mis)
  coverage21=sum(res_cover21*mis)/sum(mis)
  coverage22=sum(res_cover22*mis)/sum(mis)
  invisible(list(resu1=res_mean1,
                 resu2=res_mean2,
                 cover_underly1=coverage11,
                 cover_underly2=coverage21,
                 cover_obs1=coverage12,
                 cover_obs2=coverage22,
                 srf=mean(srf)))
}

cpimp=function(ts,num_com,max_iter = 25, tol = 1e-05){
  
  # Missing Value imputed as 0 in initial
  dims=dim(ts)
  N=length(dims)
  M=prod(dims)
  tts=1:M
  
  for(e in 1:M){
    if(ts[e]=="NA"){
      tts[e]=0
      ts[e]="0"
    }else{
      tts[e]=1
    }
    
  }
  
  ts=as.numeric(ts)
  
  ts=array(ts,dim = dims)
  
  u_list=vector("list",N)
  unfolded=vector("list",N)
  
  ts_norm=fnormimp(ts)
  
  
  for(m in 1:N){
    unfolded[[m]]=unfoldimp(ts,mode = m)
    u_list[[m]]=matrix(rnorm(dims[m]*num_com), nrow = dims[m],ncol = num_com)
  }
  est=ts
  curr_it=1
  converg=FALSE
  fresi=rep(0,max_iter)
  checkcv=function(est){
    curr_res=fnormimp(est-ts)
    fresi[curr_it]<<-curr_res
    if(curr_it==1) return(FALSE)
    if(abs(curr_res-fresi[curr_it-1])/ts_norm<tol) return(TRUE)
    else {return(FALSE)}
  }
  
  
  
  
  
  norm_vec=function(x){
    norm(as.matrix(x))
  }
  
  while ((curr_it<max_iter)&&(!converg)) {
    
    for(m in 1:N){
      v=hdlistimp(lapply(u_list[-m],function(x){
        t(x)%*%x
      }))
      vinv=solve(v)
      tmp=unfolded[[m]]%*%krlistimp(rev(u_list[-m]))%*%vinv
      lambda=apply(tmp,2,norm_vec)
      u_list[[m]]=sweep(tmp,2,lambda,"/")
      z=superdiagtsimp(num_modes = N,len = num_com,elements = lambda)
      est=ttlimp(z,u_list,ms=1:N)
    }
    
    # Imputation
    for(e in 1:M){
      if(tts[e]==0){
        ts[e]=est[e]
      }
    }
    for(m in 1:N){
      unfolded[[m]]=unfoldimp(ts,mode = m)
    }
    
    if(checkcv(est)){
      converg=TRUE
      
    }
    else{
      curr_it=curr_it+1
    }
  }
  #  mis=array(as.numeric(tts),dim = dims)
  
  #  norm_imp=(1-sum((unfoldimp(ts)*mis-unfoldimp(est)*mis)^2)/sum((unfoldimp(ts)*mis)^2))
  fresi=fresi[fresi!=0]
  norm_percent=(1-(tail(fresi,1)/ts_norm))
  invisible(list(U=u_list,est=est,fnormr=tail(fresi,1),norm_percent=norm_percent,fresi=fresi))
}

SRF=function(ch1,ch2,dim){
  ch=abind(ch1,ch2,along = 1)
  srfre=array(0,dim=dim)
  srfre=sqrt((apply(ch,c(2,3,4),var))/((apply(ch1,c(2,3,4),var)+apply(ch2,c(2,3,4),var))/2))
  mean(srfre)
}

cp_crossv_check_effi=function(ts,x,orx,frac=0.5,num_com=3,max=1000,burn.in=100,p=5){
  dims=dim(ts)
  N=length(dims)
  M=prod(dims)
  tts=1:M
  valid=1:M
  hold_m=ts
  
  for(e in 1:M){
    if(hold_m[e]=="NA"){
      hold_m[e]="0"
    }
  }
  
  
  
  
  hold_m=as.numeric(hold_m)
  
  hold_m=array(hold_m,dim = dims)
  
  
  
  
  for(e in 1:M){
    if(ts[e]=="NA"){
      tts[e]=1
      valid[e]=0
    }else{
      roll=runif(1,min = 0,max = 1)
      if(roll<=frac){
        ts[e]="NA"
        tts[e]=0
        valid[e]=1
      }else{
        tts[e]=0
        valid[e]=0
      }
    }
  }
  valid=array(valid,dim = dims)
  
  
  #mse_res=array(dim = c(p,4,2))
  mse_res=matrix(nrow = p,ncol = 11)
  
  for(ii in 1:p){
    comp_error<<-FALSE
    #comp_error2=FALSE
    tryCatch({
      a1=cpbayeimp_Jef_check_effi(ts,num_com = num_com-(p+1)/2+ii,max = max,burn.in = burn.in,thin = 10)
      a2=cpbayeimp_Jef_check_effi(ts,num_com = num_com-(p+1)/2+ii,max = max,burn.in = burn.in,cpem_as_initial = TRUE,thin = 10)
    },
    error=function(e){
      #print(e)
      comp_error<<-TRUE
      #mse_res[ii,]=rep(0,11)
    })
    #a1=cpbayeimp_Jef_check_updated(ts,num_com = num_com-3+ii,max = max,burn.in = burn.in)
    
    if(comp_error){
      mse_res[ii,]=rep(0,11)
    }else{
      result_check=CCofBayes_update(u11=a1$resu1,u12=a1$resu2,u13=a1$resu3,u21=a2$resu1,u22=a2$resu2,u23=a2$resu3,orx = orx,x=x,mis = tts,sd1 = sqrt(a1$sigma_err),sd2 = sqrt(a2$sigma_err))
      
      a1est=result_check$resu1
      a2est=result_check$resu2
      
      
      srf=result_check$srf
      
      
      mse_res[ii,1]=sum((valid*hold_m-valid*a1est)^2)/sum((valid*hold_m)^2)
      mse_res[ii,2]=sum((tts*x-tts*a1est)^2)/sum((tts*x)^2)
      mse_res[ii,3]=result_check$cover_underly1
      mse_res[ii,4]=result_check$cover_obs1
      
      
      
      mse_res[ii,5]=sum((valid*hold_m-valid*a2est)^2)/sum((valid*hold_m)^2)
      mse_res[ii,6]=sum((tts*x-tts*a2est)^2)/sum((tts*x)^2)
      mse_res[ii,7]=result_check$cover_underly2
      mse_res[ii,8]=result_check$cover_obs2
      
      mse_res[ii,9]= srf
      
      a3=cpimp(ts,num_com = 3)
      mse_res[ii,10]=sum((tts*x-tts*a3$est)^2)/sum((tts*x)^2)
      
      a4=cpimp(ts,num_com = 3)
      mse_res[ii,11]=sum((tts*x-tts*a4$est)^2)/sum((tts*x)^2)
    }
    
    
    #print(ii)
  }
  
  
  mse_res
  
  
  
}

kronecker_prod_effi=function(x,ulist){
  res=matrix(nrow = dim(x)[1],ncol = dim(x)[2])
  n=length(ulist)
  if(n==2){
    for(ii in 1:dim(x)[1]){
      a=matrix(x[ii,],nrow = dim(ulist[[2]])[1],ncol = dim(ulist[[1]])[1])
      res[ii,]=c(ulist[[2]]%*%a%*%ulist[[1]])
    }
  }else{
    a=matrix(x[ii,],nrow = dim(ulist[[2]])[1],ncol = dim(ulist[[1]])[1])
    u1=ulist[-n]
    u2=ulist[[n]]
    res[ii,]=c(kronecker_prod_effi(u2%*%a,u1))
  }
  res
}

outerimp_comp=function(a1,a2,a3){
  if(!is.null(dim(a1))){
    num_com=dim(a1)[2]
    dim1=dim(a1)[1]
    dim2=dim(a2)[1]
    dim3=dim(a3)[1]
    resu=array(0,dim = c(dim1,dim2,dim3))
    for(i in 1:num_com){
      resu=resu+outer(outer(a1[,i],a2[,i]),a3[,i])
    }
  }else{
    num_com=1
    dim1=length(a1)
    dim2=length(a2)
    dim3=length(a3)
    resu=outer(outer(a1,a2),a3)
  }
  
  
  resu
}

CCofBayes_tensor=function(u11,u12,u13,u21,u22,u23,sd1,sd2,x,mis){
  
  n=dim(u11)[1]
  tens1=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  tens2=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  tens1e=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  tens2e=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  
  for(ii in 1:n){
    tens1[ii,,,]=outerimp_comp(u11[ii,,],u12[ii,,],u13[ii,,])
    tens2[ii,,,]=outerimp_comp(u21[ii,,],u22[ii,,],u23[ii,,])
    tens1e[ii,,,]=outerimp_comp(u11[ii,,],u12[ii,,],u13[ii,,])+rnorm(1,mean = 0,sd=sd1[ii])
    tens2e[ii,,,]=outerimp_comp(u21[ii,,],u22[ii,,],u23[ii,,])+rnorm(1,mean = 0,sd=sd2[ii])
  }
  res_mean1=apply(tens1,c(2,3,4),mean)
  res_mean2=apply(tens2,c(2,3,4),mean)
  srf=SRF(tens1,tens2,c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  CCofBayes(ts=x,est = tens1,mis = mis)
  
  invisible(list(resu1=res_mean1,
                 resu2=res_mean2,
                 cover1=CCofBayes(ts=x,est = tens1e,mis = mis),
                 cover2=CCofBayes(ts=x,est = tens2e,mis = mis),
                 srf=(srf)))
}

CCofBayes=function(ts,est,alpha=0.95,mis){
  lower=(1-alpha)/2
  upper=(1+alpha)/2
  dims=dim(ts)
  N=length(dims)
  M=prod(dims)
  dims2=dim(est)
  #max=dims2[1]
  N2=length(dims2)
  checkfl=function(x){
    quantile(x,prob=lower)
  }
  checkfu=function(x){
    quantile(x,prob=upper)
  }
  lo=apply(est,2:N2,checkfl)
  up=apply(est,2:N2,checkfu)
  count=0
  for(i in 1:M){
    if((ts[i]>=lo[i])&&(ts[i]<=up[i])&&(mis[i]==1)){
      #if((ts[i]<=up[i])&&(mis[i]==1)){
      count=count+1
      #}
    }
  }
  #dim(lo)
  (count/sum(mis))
  
}

kroneckerlist=function(x){
  N=length(x)
  res=x[[1]]
  for(i in 2:N){
    res=kronecker(res,x[[i]])
  }
  res
}



CCofBayes_tensor_3=function(u11,u12,u13,u21,u22,u23,sd1,sd2,x,mis,obs){
  
  n=dim(u11)[1]
  tens1=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  tens2=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  totaln=prod(c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  mis=ifelse(mis==1,TRUE,FALSE)
  obs=ifelse(obs==1,TRUE,FALSE)
  #r1=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  #r2==array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  var_obs1=1:n
  var_obs2=1:n
  var_mis1=1:n
  var_mis2=1:n
  
  
  cor1_obs1=1:n
  cor1_obs2=1:n
  cor2_obs1=1:n
  cor2_obs2=1:n
  cor1_mis1=1:n
  cor1_mis2=1:n
  cor2_mis1=1:n
  cor2_mis2=1:n
  cor1_sim1=1:n
  cor1_sim2=1:n
  cor2_sim1=1:n
  cor2_sim2=1:n
  
  fcor1_obs1=1:n
  fcor1_obs2=1:n
  fcor2_obs1=1:n
  fcor2_obs2=1:n
  fcor1_mis1=1:n
  fcor1_mis2=1:n
  fcor2_mis1=1:n
  fcor2_mis2=1:n
  fcor1_sim1=1:n
  fcor1_sim2=1:n
  fcor2_sim1=1:n
  fcor2_sim2=1:n
  
  
  #tens1e=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  #tens2e=array(0,dim = c(n,dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  
  for(ii in 1:n){
    
    s1=outerimp_comp(u11[ii,,],u12[ii,,],u13[ii,,])
    s2=outerimp_comp(u21[ii,,],u22[ii,,],u23[ii,,])
    
    tens1[ii,,,]=s1+array(rnorm(totaln,mean = 0,sd=sd1[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))*(!obs)
    tens2[ii,,,]=s2+array(rnorm(totaln,mean = 0,sd=sd2[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))*(!obs)
    
    ####### Chain 1
    
    
    
    r1=x-tens1[ii,,,]
    
    r1obs=r1
    r1obs[!obs]=NA
    r1obs1=data.frame(apply(r1obs,2,c))
    r1obs2=data.frame(apply(r1obs,3,c))
    
    
    # The simulated residuals for the held-out structure
    
    r1mis=array(rnorm(totaln,mean=0,sd=sd1[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
    r1mis[!mis]=NA
    r1mis1=data.frame(apply(r1mis,2,c))
    r1mis2=data.frame(apply(r1mis,3,c))
    
    # The simulated residuals for the all simulated structure
    
    r1sim=array(rnorm(totaln,mean=0,sd=sd1[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
    r1sim[obs]=NA
    r1sim1=data.frame(apply(r1sim,2,c))
    r1sim2=data.frame(apply(r1sim,3,c))
    
    
    # Clculate the Variance and Contrivance
    
    var_obs1[ii]=sd(r1[obs])^2
    var_mis1[ii]=sd1[ii]^2
    cor1_obs1[ii]=det(cor(r1obs1,use = "pairwise.complete.obs"))
    cor1_obs2[ii]=det(cor(r1obs2,use = "pairwise.complete.obs"))
    fcor1_obs1[ii]=sum((cor(r1obs1,use = "pairwise.complete.obs"))^2)
    fcor1_obs2[ii]=sum((cor(r1obs2,use = "pairwise.complete.obs"))^2)
    
    cor1_mis1[ii]=det(cor(r1mis1,use = "pairwise.complete.obs"))
    cor1_mis2[ii]=det(cor(r1mis2,use = "pairwise.complete.obs"))
    fcor1_mis1[ii]=sum((cor(r1mis1,use = "pairwise.complete.obs"))^2)
    fcor1_mis2[ii]=sum((cor(r1mis2,use = "pairwise.complete.obs"))^2)
    
    cor1_sim1[ii]=det(cor(r1sim1,use = "pairwise.complete.obs"))
    cor1_sim2[ii]=det(cor(r1sim2,use = "pairwise.complete.obs"))
    fcor1_sim1[ii]=sum((cor(r1sim1,use = "pairwise.complete.obs"))^2)
    fcor1_sim2[ii]=sum((cor(r1sim2,use = "pairwise.complete.obs"))^2)
    
    
    
    ####### Chain 2
    
    
    r2=x-tens2[ii,,,]
    r2obs=r2
    r2obs[!obs]=NA
    r2obs1=data.frame(apply(r2obs,2,c))
    r2obs2=data.frame(apply(r2obs,3,c))
    
    r2mis=array(rnorm(totaln,mean=0,sd=sd2[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
    r2mis[!mis]=NA
    r2mis1=data.frame(apply(r2mis,2,c))
    r2mis2=data.frame(apply(r2mis,3,c))
    
    r2sim=array(rnorm(totaln,mean=0,sd=sd2[ii]),dim=c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
    r2sim[obs]=NA
    r2sim1=data.frame(apply(r2sim,2,c))
    r2sim2=data.frame(apply(r2sim,3,c))
    
    
    # Clculate the Variance and Contrivance
    
    var_obs2[ii]=sd(r2[obs])^2
    var_mis2[ii]=sd2[ii]^2
    cor2_obs1[ii]=det(cor(r2obs1,use = "pairwise.complete.obs"))
    cor2_obs2[ii]=det(cor(r2obs2,use = "pairwise.complete.obs"))
    fcor2_obs1[ii]=sum((cor(r2obs1,use = "pairwise.complete.obs"))^2)
    fcor2_obs2[ii]=sum((cor(r2obs2,use = "pairwise.complete.obs"))^2)
    
    cor2_mis1[ii]=det(cor(r2mis1,use = "pairwise.complete.obs"))
    cor2_mis2[ii]=det(cor(r2mis2,use = "pairwise.complete.obs"))
    fcor2_mis1[ii]=sum((cor(r2mis1,use = "pairwise.complete.obs"))^2)
    fcor2_mis2[ii]=sum((cor(r2mis2,use = "pairwise.complete.obs"))^2)
    
    cor2_sim1[ii]=det(cor(r2sim1,use = "pairwise.complete.obs"))
    cor2_sim2[ii]=det(cor(r2sim2,use = "pairwise.complete.obs"))
    fcor2_sim1[ii]=sum((cor(r2sim1,use = "pairwise.complete.obs"))^2)
    fcor2_sim2[ii]=sum((cor(r2sim2,use = "pairwise.complete.obs"))^2)
    
  }
  
  
  res_mean1=apply(tens1,c(2,3,4),mean)
  res_mean2=apply(tens2,c(2,3,4),mean)
  srf=SRF(tens1,tens2,c(dim(u11)[2],dim(u12)[2],dim(u13)[2]))
  
  #CCofBayes(ts=x,est = tens1,mis = mis)
  
  invisible(list(resu1=res_mean1,
                 resu2=res_mean2,
                 cover1=CCofBayes(ts=x,est = tens1,mis = mis),
                 cover2=CCofBayes(ts=x,est = tens2,mis = mis),
                 cor1_obs1=cor1_obs1,
                 cor1_obs2=cor1_obs2,
                 cor1_mis1=cor1_mis1,
                 cor1_mis2=cor1_mis2,
                 cor1_sim1=cor1_sim1,
                 cor1_sim2=cor1_sim2,
                 
                 cor2_obs1=cor2_obs1,
                 cor2_obs2=cor2_obs2,
                 cor2_mis1=cor2_mis1,
                 cor2_mis2=cor2_mis2,
                 cor2_sim1=cor2_sim1,
                 cor2_sim2=cor2_sim2,
                 
                 
                 
                 fcor1_obs1=fcor1_obs1,
                 fcor1_obs2=fcor1_obs2,
                 fcor1_mis1=fcor1_mis1,
                 fcor1_mis2=fcor1_mis2,
                 fcor1_sim1=fcor1_sim1,
                 fcor1_sim2=fcor1_sim2,
                 
                 fcor2_obs1=fcor2_obs1,
                 fcor2_obs2=fcor2_obs2,
                 fcor2_mis1=fcor2_mis1,
                 fcor2_mis2=fcor2_mis2,
                 fcor2_sim1=fcor2_sim1,
                 fcor2_sim2=fcor2_sim2,
                 
                 var_obs1=var_obs1,
                 var_obs2=var_obs2,
                 var_mis1=var_mis1,
                 var_mis2=var_mis2,
                 srf=(srf)))
}

CCofBayes_tensor_mix=function(resu1,resu2,x,mis,obs){
  n=dim(resu1)[1]
  
  mis=ifelse(mis==1,TRUE,FALSE)
  obs=ifelse(obs==1,TRUE,FALSE)
  
  var_obs1=1:n
  var_obs2=1:n
  var_mis1=1:n
  var_mis2=1:n
  
  cor1_obs1=1:n
  cor1_obs2=1:n
  cor2_obs1=1:n
  cor2_obs2=1:n
  cor1_mis1=1:n
  cor1_mis2=1:n
  cor2_mis1=1:n
  cor2_mis2=1:n
  cor1_sim1=1:n
  cor1_sim2=1:n
  cor2_sim1=1:n
  cor2_sim2=1:n
  
  fcor1_obs1=1:n
  fcor1_obs2=1:n
  fcor2_obs1=1:n
  fcor2_obs2=1:n
  fcor1_mis1=1:n
  fcor1_mis2=1:n
  fcor2_mis1=1:n
  fcor2_mis2=1:n
  fcor1_sim1=1:n
  fcor1_sim2=1:n
  fcor2_sim1=1:n
  fcor2_sim2=1:n
  
  for(ii in 1:n){
    
    
    
  }
  
  
}


indximp=function(dims,inx){
  M=prod(dims)
  N=length(dims)
  res=c(1:N)
  rem=inx
  J=c(1:N)
  J[1]=1
  for(i in 2:N){
    J[i]=J[i-1]*dims[i-1]
  }
  
  for(i in N:1){
    if(inx==0){
      res[i]=dims[i]
    }else{
      res[i]=floor((inx-1)/J[i])
      inx=inx-res[i]*J[i]
    }
    
  }
  #res[1]=res[1]-1
  res+1
  
}

oldunfoldimp=function(ts,mode=1){
  dims=dim(ts)
  M=prod(dims)
  N=length(dims)
  resu=matrix(data=1:M,nrow = dims[mode])
  J=c(1:N)
  J[1]=1
  for(i in 2:N){
    if(i<=mode){
      J[i]=J[i-1]*dims[i-1]
    }
    if(i==mode+1){
      J[i]=J[i-1]
    }
    if(i>mode+1){
      J[i]=J[i-1]*dims[i-1]
    }
  }
  
  
  
  
  for(i in 1:M){
    inx=indximp(dims,i)
    ij=1
    for(k in 1:N){
      if(k!=mode){
        ij=ij+(inx[k]-1)*J[k]
      }
    }
    resu[inx[mode],ij]=ts[i]
  }
  return(resu)
}

unfoldimp=function(ts,mode=1,dims=NULL){
  dims=dim(ts)
  M=prod(dims)
  N=length(dims)
  rs=mode
  cs=(1:N)[-mode]
  perm=c(rs,cs)
  newm=c(prod(dims[rs]),prod(dims[cs]))
  mat=ts
  mat=aperm(mat,perm)
  dim(mat)=newm
  mat
}

oldfoldimp=function(ts,mode=1,dims=NULL){
  
  M=prod(dims)
  N=length(dims)
  res0=1:M
  J=c(1:N)
  J[1]=1
  for(i in 2:N){
    if(i<=mode){
      J[i]=J[i-1]*dims[i-1]
    }
    if(i==mode+1){
      J[i]=J[i-1]
    }
    if(i>mode+1){
      J[i]=J[i-1]*dims[i-1]
    }
  }
  for(i in 1:M){
    inx=indximp(dims,i)
    ij=1
    for(k in 1:N){
      if(k!=mode){
        ij=ij+(inx[k]-1)*J[k]
      }
    }
    res0[i]=ts[inx[mode],ij]
  }
  res=array(res0,dim = dims)
  res
}


