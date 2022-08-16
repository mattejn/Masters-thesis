pinv <- function (A, tol = .Machine$double.eps ^ (2 / 3))
{
  #A is the matrix we want to calculate pseudoinverse for
  
  #check if the matrix is valid
  stopifnot(is.numeric(A) || is.complex(A), is.matrix(A))
  #perform singular value decomposition
  s <- svd(A)
  
  #Check if input matrix contains complex numbers
  if (is.complex(A))
    #if so, return conjugates
    s$u <- Conj(s$u)
  
  #if a calculated eigenvalue is smaller than first
  #eigenvalue*tolerance(an extremely small number)
  #assume it is zero
  p <- (s$d > max(tol * s$d[1], 0))
  
  #if all eigenvalues are non-zero
  if (all(p)) {
    #use the complete eigenvector/eigenvalue set 
    #to calculate the inverse matrix
    mp <- s$v %*% (1 / s$d * t(s$u))
  }
  
  #if 0 eigenvalues are present
  else if (any(p)) {
    #use only the non-0 eigenvectors/values 
    #to calculate the pseudoinverse
    #and fill the remainder of the matrix with zeros
    mp <- s$v[, p, drop = FALSE] %*% (1 / s$d[p] * t(s$u[,
                                        p, drop = FALSE]))
  }
  
  #if all eigenvalues are 0
  else {
    #return full zero matrix as pseudoinverse
    mp <- matrix(0, nrow = ncol(A), ncol = nrow(A))
  }
  
  #return the resulting matrix
  return(mp)
}