#' Min-Max scaling
#'
#' \code{scaling.min.max} rescales the range of a data set to [min, max].
#' @param data a numeric matrix for the input data.
#' @param min desired minimum of the transformed data.
#' @param max desired maximum of the transformed data.
#' @return a matrix containing the scaled data.
#' @examples
#' # generate some data
#' X = c(1.5, 3.2, 2.6, 0.3, 5.1)
#' Y = c(3.3, 6.4, 1.8, 0.7, 2.2)
#' 
#' # scale (X, Y) into [0, 1] range
#' scaling.min.max(cbind(X, Y), 0, 1)
#' @keywords utilities
#' @export
scaling.min.max = function(data, min = 0, max = 1) {
  data = scale(data,
               center = apply(as.matrix(data), 2, min, na.rm = TRUE),
               scale = diff(apply(as.matrix(data), 2, range, na.rm = TRUE)))
  data = data * (max - min) + min
  
  return(data)
}

#' Undo min-max scaling
#' 
#' @param data a numeric matrix for the input data.
#' @param center a numeric vector containing the locations for transforming the data.
#' @param scale a numeric vector containing the scales for transforming the data.
#' @return a matrix containing the unscaled data.
#' @examples
#' # generate some data
#' X = c(1.5, 3.2, 2.6, 0.3, 5.1)
#' Y = c(3.3, 6.4, 1.8, 0.7, 2.2)
#' 
#' # scale (X, Y) into [0, 1] range
#' Z = scaling.min.max(cbind(X, Y), 0, 1)
#' 
#' # undo the min-max scaling
#' unscaling.min.max(Z, attr(Z, "scaled:center"), attr(Z, "scaled:scale"))
#' @keywords utilities
#' @export
unscaling.min.max = function(data, center, scale) {
  data = (data - min(data)) / diff(range(data))
  data = t(t(data) * c(scale) + center)
  
  return(data)
}

#' Generate \eqn{m} Distinct Permutations of \eqn{n} Elements
#' 
#' \code{pc.perm} generates the first \eqn{m} distinct permutations of \eqn{n} elements with the Plain Changes, 
#' a.k.a. Steinhaus-Johnson-Trotter, algorithm.
#' @param n number of elements to permute.
#' @param m number of permutations to generate.
#' @return an \eqn{m} by \eqn{n} matrix of indices where each row represents a permutation. 
#' Rows beyond the \eqn{n!}-th one will be filled with 0s.
#' @examples 
#' # generate 3 permutations of 5 elements
#' pc.perm(5, 3)
#' @keywords utilities
#' @references Donald E. Knuth (2011). \emph{The Art of Computer Programming: Combinatorial Algorithms, Part 1}. 
#' Vol. 4A. Art of Computer Programming. Addison-Wesley Professional. Sec. 7.2.1.2. 
#' ISBN: 978-0-201-03804-0.
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

#' Generate \eqn{m} Distinct Permutations of \eqn{n} Elements Uniformly.
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

#' Helper for Choosing Vertices with the Min/Max \eqn{l_2} Norm.
#' 
#' @keywords internal
choose.vert = function(V, type = 1) {
  if (type == 1) {
    as.numeric(V[which.max(rowSums(V^2)), , drop = FALSE])
  } else if (type == 2) {
    as.numeric(V[which.min(rowSums(V^2)), , drop = FALSE])
  } else {
    stop("Unknown method of rank mapping!")
  }
}

#' Helper for Uniform Sampling over Convex Polygons
#' 
#' @importFrom stats runif
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

#' Helper for Assigning Ranks for Goodness-of-fit Test
#' 
#' @keywords internal
gof.assign.rank = function(elem, rank.id, d) {
  if (rank.id == 0) {
    r = elem
  } else if (rank.id == 3) {
    r = lapply(split(elem[, -1], elem[, 1]), matrix, ncol = d)
    r = uniform.rank(r)
  } else {
    r = lapply(split(elem[, -1], elem[, 1]), matrix, ncol = d)
    r = t(sapply(r, choose.vert, type = rank.id))
  }
  
  return(r)
}

#' Helper for Assigning Ranks for Test of Independence
#' 
#' @importFrom stats runif
#' @keywords internal
dep.assign.rank = function(XY, elem, rank.id) {
  N = NROW(XY)
  rank.x = rank(XY[, 1])
  rank.y = rank(XY[, 2])
  
  # rh is the joint rank, rt is the usual rank
  if (rank.id == 0) {
    rh = elem
    rt = cbind((2 * rank.x - 1) / (2 * N), (2 * rank.y - 1) / (2 * N))
  } else if (rank.id == 3) {
    rh = uniform.rank(lapply(split(elem[, -1], elem[, 1]), matrix, ncol = 2))
    rt = cbind(runif(N, min = rank.x - 1, max = rank.x) / N, 
               runif(N, min = rank.y - 1, max = rank.y) / N)
  } else {
    rh = lapply(split(elem[, -1], elem[, 1]), matrix, ncol = 2)
    rh = t(sapply(rh, choose.vert, type = rank.id))
    rt = cbind((rank.x - (rank.id == 2)) / N, 
               (rank.y - (rank.id == 2)) / N)
  }
  
  return(list(h = rh, t = rt))
}

#' Helper for extracting line segments for visualization
#' 
#' @keywords internal
extract.segment = function(V, id_col, vert_col) {
  n.edge = NROW(V)
  v.edges = matrix(0, n.edge, 2 * length(vert_col))
  curr.row = 1
  
  # extract edges
  for (i in unique(V[, id_col])) {
    curr.cell = subset(V, V[, id_col] == i, select = vert_col)
    curr.cell.sum = rowSums(curr.cell)
    for (j in 1:NROW(curr.cell)) {
      p.id = j
      q.id = 1 + j %% NROW(curr.cell)
      # sort by l1 norm, everything positive!
      if (curr.cell.sum[p.id] <= curr.cell.sum[q.id]) {
        v.edges[curr.row, ] = c(curr.cell[p.id, ], curr.cell[q.id, ])
      } else {
        v.edges[curr.row, ] = c(curr.cell[q.id, ], curr.cell[p.id, ])
      }
      curr.row = curr.row + 1
    }
  }
  
  # remove duplicated edges
  v.edges = unique(round(v.edges, 8))
  
  return(v.edges)
}