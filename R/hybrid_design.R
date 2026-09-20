# Design-construction functions used in the paper.

library(lpSolve)

#=====Functions for constructing hybrid design==================================
# function to get the target concurrence matrix 'NNPo' required for hybrid design
NNPo_for_alg = function(N,num_obs, gsize, w.C2){
  target1 = (N*num_obs)/gsize
  target2 = ((N*(num_obs - 1))/(gsize - 1)) * (w.C2+((1-w.C2)*num_obs/gsize))
  area = sum(seq(1,gsize - 2))
  target3 = (choose(num_obs, 2)*N - target2*(gsize-1))/area
  
  NNPo = matrix(NA, nrow = gsize, ncol = gsize)
  diag(NNPo) = target1
  diag(NNPo[-1,]) = diag(NNPo[,-1]) = target2
  NNPo[which(is.na(NNPo))] = target3
  return(NNPo)
}

# function to find solution for one row of the incidence matrix
LIP_dispersed=function(gsize,N,num_obs,NNPo,colPos,N1,T_mat,rownum,relaxed,dist = 5){
  #   gsize: grid size
  #       N: number of subjects
  # num_obs: number of observations per subjects
  #    NNPo: target concurrence matrix
  #  colPos: target incidence matrix only for snippet portion
  #      N1: incidence matrix
  #   T_mat: matrix to store old solutions
  #  rownum: current row 
  # relaxed: if maximal number reached and no solution found, relax constraint.
  #    dist: for snippet portion how many grid points away from the cluster
  
  total_col = gsize+choose(gsize, 2)
  combs = combn(gsize,2)
  total_combs = choose(gsize, 2)
  gap1 = which(diff(combs) == 1)
  rvec_obt=t(N1)%*%matrix(1,nrow(N1),1)
  w=matrix(0,1,gsize)
  for (j in 1:gsize){
    if (rvec_obt[j,]==0) w[,j]=1
    else w[,j]=1/rvec_obt[j,]
  }		
  obj=cbind(w, t(rep(0, total_combs))) 
  constr1=matrix(rep(c(1,0), times = c(gsize, total_col - gsize)),1,total_col) 
  constr2 = matrix(0, total_combs*4, total_col) 
  for(i in 1:total_combs){ 
    constr2[i+3*(i-1), combs[1,i]] = 1
    constr2[i+3*(i-1), gsize+i] = -1
    constr2[i+3*(i-1)+1, combs[2,i]] = 1
    constr2[i+3*(i-1)+1, gsize+i] = -1
    constr2[i+3*(i-1)+2, combs[,i]] = 1
    constr2[i+3*(i-1)+2, gsize+i] = -1
    constr2[i+3*(i-1)+3, gsize+i] = 1
  }
  
  #make-up process T
  constr3 = cbind(T_mat, matrix(0, nrow(T_mat), total_combs)) 
  for(i in 1:total_combs){
    constr3[,i+gsize] = T_mat[,combs[1,i]]*T_mat[,combs[2,i]]
  }
  
  # constr4: decide the cluster points for snippet subjects
  constr4 = c(colPos[rownum,], rep(0, total_combs)) 
  
  # constr5: jump points for snippet designs
  if(rownum <= max(which(rowSums(colPos)>0))){
    from = which(colPos[rownum,]!=0)[1]
    end = which(colPos[rownum,]!=0)[2]
    grid = 1:gsize
    cand_indx = c(rep(0, gsize), rep(0, total_combs))
    jump_left <- from - dist
    jump_right <- end+dist
    pts_exclude <- seq(jump_left,jump_right,by=1)
    nonjump_cand <- grid[grid %in% pts_exclude & !grid %in% which(colPos[rownum,]!=0)]
    cand_indx[nonjump_cand]=1
    constr5 = cand_indx
  }else{ #BIBD subjects
    constr5 = c(rep(0, gsize), rep(0, total_combs))
  }
  
  constr=rbind(constr1,constr2,constr3,constr4,constr5)
  if (relaxed>0){
    constr=rbind(constr1,constr3,constr4,constr5)
  }
  
  dir1=rep("==", times=(1))
  dir2=rep(c(">=", ">=", "<=", "<="),times=(total_combs))  
  dim(dir2)=c(4*total_combs,1)	
  
  dir3=rep("<",times=(nrow(constr3)))
  dim(dir3)=c(nrow(constr3),1)
  dir4=rep("==", times=(1))
  dir5=rep("<=", times=(1))
  dir=rbind(dir1,dir2,dir3,dir4,dir5)	
  if (relaxed>0){
    dir=rbind(dir1,dir3,dir4,dir5)
  }
  
  rhs1=num_obs
  rhs2=c()
  if(rownum > max(which(rowSums(colPos)>0))){
    for (i in 1:total_combs){
      if(i %in% gap1){
        rhs2 = c(rhs2, c(0, 0, 1, 0))
      }else{
        rhs2 = c(rhs2, c(0, 0, 1, NNPo[combs[1,i], combs[2,i]] +5 - t(N1[,combs[1,i]])%*%N1[,combs[2,i]]))
      }
    }		
  }else{
    for (i in 1:total_combs){
      rhs2 = c(rhs2, c(0, 0, 1, NNPo[combs[1,i], combs[2,i]] +5 - t(N1[,combs[1,i]])%*%N1[,combs[2,i]]))
    }
  }
  rhs2 = matrix(rhs2, total_combs*4, 1)
  rhs3=matrix((num_obs-0.5),nrow(constr3),1)
  rhs4=sum(colPos[rownum,])
  rhs5 = 0
  
  rhs=rbind(rhs1,rhs2,rhs3,rhs4,rhs5)
  if (relaxed>0){
    rhs=rbind(rhs1,rhs3,rhs4,rhs5)
  }
  
  sol=lp (direction = "max", obj, constr, dir, rhs,transpose.constraints = TRUE, all.bin=TRUE, use.rw=TRUE, timeout = 15)
  if (sol[[28]]==0) {		
    row=sol[[12]][1:gsize]
    dim(row)=c(1,gsize)
    if (rownum>nrow(N1)) N1=rbind(N1,row) else N1[rownum,]=row		
  } 	
  return(N1)	
}

detect=function(gsize,N,num_obs,NNPo,colPos,N1,T_mat,relaxed){
  row_detected=0	
  result=0
  k0=1   
  while (k0<=min(4,nrow(N1)) & row_detected==0){
    row_indices=combn(nrow(N1),k0)
    nr=ncol(row_indices)
    j=1
    while(j<=nr & row_detected==0){
      rows=row_indices[,j]
      T_temp=rbind(T_mat,N1[rows,])			
      N1_temp=N1
      N1_temp[rows,]=matrix(0,1,gsize)
      cnt=0
      for (m in 1:k0){
        rownum=rows[m]				
        N1_temp=LIP_dispersed(gsize,N,num_obs,NNPo,colPos,N1_temp,T_temp,rownum,relaxed)
        if (sum(N1_temp[rownum,])>0) cnt=cnt+1
      }
      if (cnt==k0) {
        row_detected=1
        result=list(rows,N1_temp)
      }
      j=j+1	
    }
    k0=k0+1				
  }
  return(result)
}

matSnippet = function(gsize,N,w){
  snippet_rows = floor(N*w)
  full = floor(snippet_rows/(gsize-1))
  left = snippet_rows%%(gsize-1)
  full_row_indx = rep(1:(gsize-1), full)
  left_row_indx = round(seq(1,(gsize-1), length = left))
  row_indx = c(full_row_indx,  left_row_indx)
  mat = matrix(0, nrow = N, ncol = gsize)
  for(i in 1:snippet_rows){
    mat[i,c(row_indx[i], row_indx[i]+1)] = 1
  }
  return(mat)
}

distant_pts = function(cluster_pos, dist, gsize){
  from = cluster_pos[1]
  end = cluster_pos[2]
  grid = 1:gsize
  cand_indx = rep(0, gsize)
  jump_left <- from - dist
  jump_right <- end + dist
  pts_exclude <- seq(jump_left,jump_right,by=1)
  cand <- grid[!grid %in% pts_exclude & !grid %in% cluster_pos]
  return(cand)
}

# the main function to build design that utilizes 'LIP_dispersed' and 'detect'
hybrid_design=function(sim,gsize,N,num_obs,NNPo,colPos,ntrial=5,dist = 5){
  ntrial = 5
  trial=1
  success=0				
  N1=matrix(0,1,gsize)
  col=c(which(colPos[1,]==1), sample(distant_pts(which(colPos[1,]==1), dist=5, gsize), num_obs-2))
  N1[1,col]=1
  T_mat=matrix(0,1,gsize)
  i=2	
  decision=0
  relaxed=0			
  while (i<=N & decision==0){
    nt=nrow(T_mat)
    if(nt>2*N){
      T_mat=matrix(0,1,gsize)
      decision=1 						
    }
    N1=LIP_dispersed(gsize,N,num_obs,NNPo,colPos,N1,T_mat,i,relaxed)
    if(nrow(N1)<i){
      temp=detect(gsize,N,num_obs,NNPo,colPos,N1,T_mat,relaxed)
      rows=temp[[1]]
      if (all(rows>0)) {
        T_mat=rbind(T_mat,N1[rows,])
        N1=temp[[2]]												
      } else {
        decision=1 								
      }			
    } 					
    if(nrow(N1)<i & trial==ntrial){
      relaxed=1							
    } 
    Sys.sleep(0.1)
    i=nrow(N1)+1	
  }
  
  if (nrow(N1)==N){
    success=1
    result=N1
  } else{
    design="Design not found"	
    result=list(gsize=gsize,N=N,num_obs=num_obs,design=design)
  }
  return(result)			
}

#=====Wrapper===================================================================
# Construct one Hybrid design from the study-level design parameters.
#
# This wrapper prepares the target snippet positions and target concurrence
# matrix required by hybrid_design(), so users do not need to construct those
# internal objects themselves.
generate_hybrid_design <- function(N, gsize, num_obs, W = 0.5,
                                   seed = NULL, simulation_id = 1,
                                   ntrial = 5, dist = 5) {
  if (length(N) != 1 || N < 1 || N != as.integer(N)) {
    stop("N must be a positive integer.")
  }
  if (length(gsize) != 1 || gsize < 2 || gsize != as.integer(gsize)) {
    stop("gsize must be an integer of at least 2.")
  }
  if (length(num_obs) != 1 || num_obs < 2 || num_obs > gsize ||
      num_obs != as.integer(num_obs)) {
    stop("num_obs must be an integer between 2 and gsize.")
  }
  if (length(W) != 1 || !is.finite(W) || W <= 0 || W > 1) {
    stop("W must be in (0, 1].")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  colPos <- matSnippet(gsize = gsize, N = N, w = W)
  NNPo <- NNPo_for_alg(N = N, num_obs = num_obs, gsize = gsize, w.C2 = W)

  hybrid_design(
    sim = simulation_id,
    gsize = gsize,
    N = N,
    num_obs = num_obs,
    NNPo = NNPo,
    colPos = colPos,
    ntrial = ntrial,
    dist = dist
  )
}
#=====Wrapper===================================================================

