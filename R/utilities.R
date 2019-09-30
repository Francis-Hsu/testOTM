#' Min-max scaling
#'
#' \code{scaling.min.max} rescales the range of a data set to [min, max].
#' @param data a numeric matrix for the input data.
#' @param min desired minimum of the transformed data.
#' @param max desired maximum of the transformed data.
#' @return a matrix containing the scaled data.
#' @keywords utilities
#' @export
scaling.min.max = function(data, min = 0, max = 1) {
  data = scale(data,
               center = apply(as.matrix(data), 2, min, na.rm = T),
               scale = diff(apply(as.matrix(data), 2, range, na.rm = T)))
  data = data * (max - min) + min
  
  return(data)
}

#' Undo min-max scaling
#' 
#' @param data a numeric matrix for the input data.
#' @param center a numeric vector containing the locations for transforming the data.
#' @param scale a numeric vector containing the scales for transforming the data.
#' @return a matrix containing the unscaled data.
#' @keywords utilities
#' @export
unscaling.min.max = function(data, center, scale) {
  data = (data - min(data)) / diff(range(data))
  data = t(t(data) * scale + center)
  
  return(data)
}

#' Generate \eqn{m} Distinct Permutations of \eqn{n} Elements
#' 
#' \code{pc.perm} generates the first \eqn{m} distinct permutations of \eqn{n} elements with the Plain Changes, 
#' a.k.a. Steinhaus-Johnson-Trotter, algorithm.
#' @param n number of elements to permute.
#' @param m number of permutations to generate.
#' @return a \eqn{m} by \eqn{n} matrix of indices where each row represents a permutation. 
#' Rows beyond the \eqn{n!}-th one will be filled with 0s.
#' @keywords utilities
#' @references Donald E. Knuth. \emph{The Art of Computer Programming: Combinatorial Algorithms, Part 1}. 1st. 
#' Vol. 4A. Art of Computer Programming. Addison-Wesley Professional, 2011. Sec. 7.2.1.2. 
#' isbn: 978-0-201-03804-0.
#' @export
pc.perm = function(n, m) {
  if (n < 1 || m < 1 || (n + m) %% 1 != 0) {
    stop("n and m must be integers greater than or equal to 1.")
  }
  
  if (m > factorial(n)) {
    warning(paste0(m, " exceeds the number of permutations of ", n, " elements."))
  }
  
  # Storage for result
  perm.indices = matrix(0, m, n)
  perm.indices[1, ] = 1:n
  
  # P1: Initialize.
  c = rep(0, n)
  o = rep(1, n)
  
  # P3: Prepare for change.
  i = 1
  j = n
  s = 0
  
  repeat {
    if (i >= m) {
      break
    }
    i = i + 1
    
    # P4: Ready to change?
    q = c[j] + o[j]
    
    if (q < 0) {
      # P7: Switch direction.
      o[j] = -o[j]
      j = j - 1
      i = i - 1
      next
    } else if (q == j) {
      # P6: Increase s.
      if (j == 1) {
        break
      } else {
        s = s + 1
        
        # P7: Switch direction.
        o[j] = -o[j]
        j = j - 1
        i = i - 1
        next
      }
    }
    
    # P5: Change.
    perm.indices[i, ] = perm.indices[i - 1, ]
    perm.indices[i, j - q + s] = perm.indices[i - 1, j - c[j] + s]
    perm.indices[i, j - c[j] + s] = perm.indices[i - 1, j - q + s]
    c[j] = q
    
    # Repeat P3
    j = n
    s = 0
  }
  
  return(perm.indices)
}

#' Generate \eqn{m} Distinct Permutations of \eqn{n} Elements Randomly.
#' 
#' @keywords internal
rand.perm = function(n, m) {
  if (n < 1 || m < 1 || (n + m) %% 1 != 0) {
    stop("n and m must be integers greater than or equal to 1.")
  }
  
  if (m > factorial(n)) {
    warning(paste0(m, " exceeds the number of permutations of ", n, " elements."))
  }
  
  # Storage for result
  perm.indices = matrix(0, m, n)
  uniqueness = 1
  for (i in 1:m) {
    while(uniqueness != i + 1) {
      perm.indices[i, ] = sample(1:n)
      uniqueness = nrow(unique(perm.indices))
      if (uniqueness == m) {
        break
      }
    }
  }
  
  return(perm.indices)
}

#' Helper for choosing vertices with the min/max l2 norm.
#' 
#' @keywords internal
choose.vert = function(V, type = 1) {
  if (type == 1) {
    as.numeric(V[which.max(rowSums(V^2)), , drop = F])
  } else if (type == 2) {
    as.numeric(V[which.min(rowSums(V^2)), , drop = F])
  } else {
    stop("Unknown method of rank mapping!")
  }
}

#' Helper for uniform sampling over convex polygons
#' 
#' @keywords internal
uniform.rank = function(V) {
  n.cell = length(V)
  
  # random numbers for barycentric sampling
  r1 = sqrt(runif(n.cell))
  r2 = runif(n.cell)
  
  # sample 1 point from each cell using fan triangulation implicitly
  unif.rank = matrix(0, n.cell, 2)
  for (i in 1:n.cell) {
    n.vert = nrow(V[[i]])
    
    # area of the fan triangles
    tri.area = rep(0, n.vert - 2)
    A = V[[i]][1, ]
    for (j in 1:(n.vert - 2)) {
      B = V[[i]][j + 1, ]
      C = V[[i]][j + 2, ]
      tri.area[j] = 0.5 * abs(A[1] * (B[2] - C[2]) + B[1] * (C[2] - A[2]) + C[1] * (A[2] - B[2]))
    }
    
    # sample a triangle, then sample uniformly from that triangle
    tri.id = sample(1:(n.vert - 2), 1, prob = tri.area)
    B = V[[i]][tri.id + 1, ]
    C = V[[i]][tri.id + 2, ]
    unif.rank[i, 1] = (1 - r1[i]) * A[1] + r1[i] * (1 - r2[i]) * B[1] + r1[i] * r2[i] * C[1]
    unif.rank[i, 2] = (1 - r1[i]) * A[2] + r1[i] * (1 - r2[i]) * B[2] + r1[i] * r2[i] * C[2]
  }
  
  return(unif.rank)
}