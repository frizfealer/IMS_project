function A=spdiag(vec) 
A = spdiags(vec,0,length(vec),length(vec));