# This file is a modified version of the gllim() function from the
# xLLiM R package.
#
# Original package:
#   xLLiM: High Dimensional Locally-Linear Mapping
#   Authors: Emeline Perthame, Florence Forbes, Antoine Deleforge,
#            Emilie Devijver, Melina Gallopin
#   License: GNU General Public License (GPL) version 2 or later
#
# Original source:
#   https://github.com/epertham/xLLiM
#   Version: 2.3
#   Commit: dba703acda2955c9462c0afab8b2b9e6ce45827f
#
# Modifications:
#   Umberto Picchini modified the EM stopping criterion, see lines 590-601
#
# This modified version is distributed under the same GPL-2-or-later
# license.

loggausspdf <- getFromNamespace("loggausspdf", "xLLiM")
logsumexp   <- getFromNamespace("logsumexp", "xLLiM")
ginv <- MASS::ginv

# the original gllim code from xLLiM modified to have a different termination criterion for EM
gllim_strict_em = function (tapp, yapp, in_K, in_r = NULL, maxiter = 100, Lw = 0, 
          cstr = NULL, verb = 0, in_theta = NULL, em_tol = 1e-6,
          min_iter = 30,
          patience = 5, ...) # "patience > 1" means that several consecutive small improvements are required, instead of just one. It avoids stopping because of a single flat step.
{
  
  if (!"ct" %in% names(cstr)) 
    cstr$ct = NULL
  if (!"cw" %in% names(cstr)) 
    cstr$cw = NULL
  if (!"Gammat" %in% names(cstr)) 
    cstr$Gammat = NULL
  if (!"Gammaw" %in% names(cstr)) 
    cstr$Gammaw = NULL
  if (!"pi" %in% names(cstr)) 
    cstr$pi = NULL
  if (!"A" %in% names(cstr)) 
    cstr$A = NULL
  if (!"b" %in% names(cstr)) 
    cstr$b = NULL
  if (!"Sigma" %in% names(cstr)) 
    cstr$Sigma = "i"
  if (ncol(tapp) != ncol(yapp)) {
    stop("Observations must be in columns and variables in rows")
  }
  covParBloc_EM = function(yCen, rk_bar, model) {
    covY = tcrossprod(yCen)/rk_bar
    Dim = dim(covY)[1]
    covEmpB = matrix(0, nrow = Dim, ncol = Dim)
    for (b in 1:max(model)) {
      ind = which(model == b)
      covEmpB[ind, ind] = covY[ind, ind]
    }
    return(covEst = covEmpB)
  }
  ExpectationZ = function(tapp, yapp, th, verb) {
    if (verb >= 1) 
      print("  EZ")
    if (verb >= 3) 
      print(" k=")
    D = nrow(yapp)
    N = ncol(yapp)
    K = length(th$pi)
    Lt = nrow(tapp)
    L = nrow(th$c)
    Lw = L - Lt
    logr = matrix(NaN, N, K)
    for (k in 1:K) {
      if (verb >= 3) 
        print(k)
      muyk = th$b[, k, drop = FALSE]
      covyk = th$Sigma[, , k]
      if (Lt > 0) {
        if (L == 1) {
          Atk = th$A[, 1:Lt, k, drop = FALSE]
        }
        else {
          Atk = th$A[, 1:Lt, k]
        }
        muyk = sweep(Atk %*% tapp, 1, muyk, "+")
      }
      if (Lw > 0) {
        Awk = matrix(th$A[, (Lt + 1):L, k, drop = FALSE], 
                     ncol = Lw, nrow = D)
        Gammawk = th$Gamma[(Lt + 1):L, (Lt + 1):L, k]
        cwk = th$c[(Lt + 1):L, k]
        covyk = covyk + Awk %*% Gammawk %*% t(Awk)
        muyk = sweep(muyk, 1, Awk %*% cwk, "+")
      }
      logr[, k] = log(th$pi[k])
      logr[, k] = logr[, k] + loggausspdf(yapp, muyk, covyk)
      if (Lt > 0) 
        logr[, k] = logr[, k] + loggausspdf(tapp, th$c[1:Lt, 
                                                       k, drop = FALSE], th$Gamma[1:Lt, 1:Lt, k])
    }
    lognormr = logsumexp(logr, 2)
    LL = sum(lognormr)
    r = exp(sweep(logr, 1, lognormr, "-"))
    ec = rep(TRUE, K)
    for (k in 1:K) {
      if (sum(r[, k]) == 0 | !is.finite(sum(r[, k]))) {
        ec[k] = FALSE
        if (verb >= 1) {
          print(paste("     WARNING: CLASS ", k, " HAS BEEN REMOVED"))
        }
      }
    }
    if (sum(ec) == 0) {
      print("REINIT! ")
      r = emgm(rbind(tapp, yapp), K, 2, verb)$R
      ec = rep(TRUE, ncol(r))
    }
    else {
      r = r[, ec, drop = FALSE]
    }
    return(list(r = r, LL = LL, ec = ec))
  }
  ExpectationW = function(tapp, yapp, th, verb) {
    if (verb >= 1) 
      print("  EW")
    if (verb >= 3) 
      print(" k=")
    D = nrow(yapp)
    N = ncol(yapp)
    K = length(th$pi)
    Lt = nrow(tapp)
    L = nrow(th$c)
    Lw = L - Lt
    if (Lw == 0) {
      muw = NULL
      Sw = NULL
    }
    Sw = array(0, dim = c(Lw, Lw, K))
    muw = array(0, dim = c(Lw, N, K))
    for (k in 1:K) {
      if (verb >= 3) 
        print(k)
      Atk = th$A[, 1:Lt, k]
      Sigmak = th$Sigma[, , k]
      if (Lw == 0) {
        Awk = NULL
        Gammawk = NULL
        cwk = NULL
        invSwk = NULL
      }
      else {
        Awk = th$A[, (Lt + 1):L, k]
        Gammawk = th$Gamma[(Lt + 1):L, (Lt + 1):L, k]
        cwk = th$c[(Lt + 1):L, k]
        invSwk = diag(Lw) + tcrossprod(Gammawk, Awk) %*% 
          solve(Sigmak) %*% Awk
      }
      if (!is.null(tapp)) {
        Atkt = Atk %*% tapp
      }
      else {
        Atkt = 0
      }
      if (Lw == 0) {
        muw = NULL
        Sw = NULL
      }
      else {
        muw[, , k] = solve(invSwk, sweep(Gammawk %*% 
                                           t(Awk) %*% solve(Sigmak) %*% sweep(yapp - Atkt, 
                                                                              1, th$b[, k], "-"), 1, cwk, "+"))
        Sw[, , k] = solve(invSwk, Gammawk)
      }
    }
    return(list(muw = muw, Sw = Sw))
  }
  Maximization = function(tapp, yapp, r, muw, Sw, cstr, verb, 
                          model) {
    if (verb >= 1) 
      print("  M")
    if (verb >= 3) 
      print(" k=")
    K = ncol(r)
    D = nrow(yapp)
    N = ncol(yapp)
    Lt = nrow(tapp)
    Lw = ifelse(is.null(muw), 0, nrow(muw))
    L = Lt + Lw
    th = list()
    th$c = matrix(NaN, nrow = L, ncol = K)
    th$Gamma = array(0, dim = c(L, L, K))
    if (Lw > 0) {
      th$c[(Lt + 1):L, ] = cstr$cw
      th$Gamma[(Lt + 1):L, (Lt + 1):L, ] = cstr$Gammaw
    }
    th$pi = rep(NaN, K)
    th$A = array(NaN, dim = c(D, L, K))
    th$b = matrix(NaN, nrow = D, ncol = K)
    th$Sigma = array(NaN, dim = c(D, D, K))
    rk_bar = rep(0, K)
    for (k in 1:K) {
      if (verb >= 3) 
        print(k)
      rk = r[, k]
      rk_bar[k] = sum(rk)
      if (Lt > 0) {
        if (verb >= 3) {
          print("c")
        }
        if (is.null(cstr$ct)) {
          th$c[1:Lt, k] = rowSums(sweep(tapp, 2, rk, 
                                        "*"))/rk_bar[k]
        }
        else {
          th$c[1:Lt, k] = cstr$ct[, k]
        }
        if (verb >= 3) {
          print("Gt")
        }
        diffGamma <- sweep(sweep(tapp, 1, th$c[1:Lt, 
                                               k], "-"), 2, sqrt(rk), "*")
        if (is.null(cstr$Gammat) || (length(cstr$Gammat) == 
                                     1 & cstr$Gammat == "*")) {
          th$Gamma[1:Lt, 1:Lt, k] = tcrossprod(diffGamma)/rk_bar[k]
        }
        else {
          if (!is.character(cstr$Gammat)) {
            th$Gamma[1:Lt, 1:Lt, k] = cstr$Gammat[, , 
                                                  k]
          }
          else {
            if (cstr$Gammat[1] == "d" | cstr$Gammat[1] == 
                "i") {
              gamma2 = rowSums(diffGamma^2)/rk_bar[k]
              if (cstr$Gammat[1] == "d") {
                th$Gamma[1:Lt, 1:Lt, k] = diag(gamma2)
              }
              else {
                th$Gamma[1:Lt, 1:Lt, k] = mean(gamma2) * 
                  diag(Lt)
              }
            }
            else {
              if (cstr$Gammat[1] == "v") {
                th$Gamma[1:Lt, 1:Lt, k] = tcrossprod(diffGamma)/rk_bar[k]
              }
              else {
                stop("  ERROR: invalid constraint on Gamma.")
              }
            }
          }
        }
      }
      th$pi[k] = rk_bar[k]/N
      if (Lw > 0) {
        x = rbind(tapp, muw[, , k])
        Skx = rbind(cbind(matrix(0, Lt, Lt), matrix(0, 
                                                    Lt, Lw)), cbind(matrix(0, Lw, Lt), Sw[, , k]))
      }
      else {
        x = tapp
        Skx = matrix(0, Lt, Lt)
      }
      if (verb >= 3) {
        print("A")
      }
      if (is.null(cstr$b)) {
        yk_bar = rowSums(sweep(yapp, 2, rk, "*"))/rk_bar[k]
        if (L > 0) {
          xk_bar = rowSums(sweep(x, 2, rk, "*"))/rk_bar[k]
        }
        else {
          xk_bar = NULL
        }
      }
      else {
        yk_bar = cstr$b[, k]
        xk_bar = rep(0, L)
        th$b[, k] = cstr$b[, k]
      }
      weights = sqrt(rk)
      y_stark = sweep(yapp, 1, yk_bar, "-")
      y_stark = sweep(y_stark, 2, weights, "*")
      if (L > 0) {
        x_stark = sweep(x, 1, xk_bar, "-")
        x_stark = sweep(x_stark, 2, weights, "*")
      }
      else {
        x_stark = NULL
      }
      if (!all(Skx == 0)) {
        if (N >= L & det(Skx + tcrossprod(x_stark)) > 
            10^(-8)) {
          th$A[, , k] = tcrossprod(y_stark, x_stark) %*% 
            qr.solve(Skx + tcrossprod(x_stark))
        }
        else {
          th$A[, , k] = tcrossprod(y_stark, x_stark) %*% 
            ginv(Skx + tcrossprod(x_stark))
        }
      }
      else {
        if (!all(x_stark == 0)) {
          if (N >= L & det(tcrossprod(x_stark)) > 10^(-8)) {
            th$A[, , k] = tcrossprod(y_stark, x_stark) %*% 
              qr.solve(tcrossprod(x_stark))
          }
          else {
            if (N < L && det(crossprod(x_stark)) > 10^(-8)) {
              th$A[, , k] = y_stark %*% solve(crossprod(x_stark)) %*% 
                t(x_stark)
            }
            else {
              if (verb >= 3) 
                print("p")
              th$A[, , k] = y_stark %*% ginv(x_stark)
            }
          }
        }
        else {
          if (verb >= 1 & L > 0) 
            print("null var\n")
          th$A[, , k] = 0
        }
      }
      if (verb >= 3) 
        print("b")
      if (L > 0) {
        wk = yapp - th$A[, , k] %*% x
      }
      else {
        wk = yapp
      }
      if (is.null(cstr$b)) 
        th$b[, k] = rowSums(sweep(wk, 2, rk, "*"))/rk_bar[k]
      if (verb >= 3) 
        print("S")
      if (Lw > 0) {
        Awk = th$A[, (Lt + 1):L, k]
        Swk = Sw[, , k]
        ASAwk = Awk %*% tcrossprod(Swk, Awk)
      }
      else ASAwk = 0
      diffSigma = sweep(sweep(wk, 1, th$b[, k], "-"), 2, 
                        sqrt(rk), "*")
      if (cstr$Sigma %in% c("", "*")) {
        th$Sigma[, , k] = tcrossprod(diffSigma)/rk_bar[k]
        th$Sigma[, , k] = th$Sigma[, , k] + ASAwk
      }
      else {
        if (!is.character(cstr$Sigma)) {
          th$Sigma = cstr$Sigma
        }
        else {
          if (cstr$Sigma[1] == "d" || cstr$Sigma[1] == 
              "i") {
            sigma2 = rowSums(diffSigma^2)/rk_bar[k]
            if (cstr$Sigma[1] == "d") {
              th$Sigma[, , k] = diag(sigma2, ncol = D, 
                                     nrow = D)
              if (is.null(dim(ASAwk))) {
                th$Sigma[, , k] = th$Sigma[, , k] + diag(ASAwk, 
                                                         ncol = D, nrow = D)
              }
              else {
                th$Sigma[, , k] = th$Sigma[, , k] + diag(diag(ASAwk))
              }
            }
            else {
              th$Sigma[, , k] = mean(sigma2) * diag(D)
              if (is.null(dim(ASAwk))) {
                th$Sigma[, , k] = th$Sigma[, , k] + sum(diag(ASAwk, 
                                                             ncol = D, nrow = D))/D * diag(D)
              }
              else {
                th$Sigma[, , k] = th$Sigma[, , k] + sum(diag(ASAwk))/D * 
                  diag(D)
              }
            }
          }
          else {
            if (cstr$Sigma == "bSHOCK") {
              th$Sigma[, , k] = covParBloc_EM(yCen = diffSigma, 
                                              rk_bar = rk_bar[k], model = model[[k]])
            }
            else {
              stop("  ERROR: invalid constraint on Sigma.")
            }
          }
        }
      }
      if (verb >= 3) 
        print("n")
      if (!is.finite(sum(th$Gamma[1:Lt, 1:Lt, k]))) {
        th$Gamma[1:Lt, 1:Lt, k] = 0
      }
      th$Gamma[1:Lt, 1:Lt, k] = th$Gamma[1:Lt, 1:Lt, k] + 
        1e-08 * diag(Lt)
      if (!is.finite(sum(th$Sigma[, , k]))) {
        th$Sigma[, , k] = 0
      }
      th$Sigma[, , k] = th$Sigma[, , k] + 1e-08 * diag(D)
      if (verb >= 3) 
        print(",")
    }
    if (verb >= 3) 
      print("end")
    if (cstr$Sigma == "*") {
      th$Sigma = sweep(th$Sigma, 3, rk_bar, "*")
      th$Sigma = array(apply(th$Sigma, c(1, 2), mean), 
                       dim = c(D, D, K))
    }
    if (!is.null(cstr$Gammat) && cstr$Gammat == "v") {
      detG = rep(0, K)
      for (k in 1:K) {
        if (D == 1) {
          detG[k] = th$Gamma[1:Lt, 1:Lt, k]
        }
        else {
          detG[k] = det(th$Gamma[1:Lt, 1:Lt, k])
        }
        th$Gamma[1:Lt, 1:Lt, k] = th$Gamma[1:Lt, 1:Lt, 
                                           k]/detG[k]
      }
      th$Gamma[1:Lt, 1:Lt, ] = sum(detG^(1/Lt) * th$pi) * 
        th$Gamma[1:Lt, 1:Lt, ]
    }
    if (is.character(cstr$Gammat) && !is.null(cstr$Gammat) && 
        cstr$Gammat[length(cstr$Gammat)] == "*") {
      for (k in 1:K) {
        th$Gamma[1:Lt, 1:Lt, k] = th$Gamma[1:Lt, 1:Lt, 
                                           k] %*% diag(rk_bar)
        th$Gamma[1:Lt, 1:Lt, k] = matrix(1, Lt, Lt) * 
          sum(th$Gamma[1:Lt, 1:Lt, k])/N
      }
    }
    if (!is.character(cstr$pi) || is.null(cstr$pi)) {
      if (!is.null(cstr$pi)) {
        th$pi = cstr$pi
      }
    }
    else {
      if (!is.null(cstr$pi) && cstr$pi[1] == "*") {
        th$pi = 1/K * rep(1, K)
      }
      else {
        stop("  ERROR: invalid constraint on pi.")
      }
    }
    return(th)
  }
  remove_empty_clusters = function(th, cstr, ec) {
    if (sum(ec) != length(ec)) {
      if (!is.null(cstr$ct) && !is.character(cstr$ct)) 
        cstr$ct = cstr$ct[, ec]
      if (!is.null(cstr$cw) && !is.character(cstr$cw)) 
        cstr$cw = cstr$cw[, ec]
      if (!is.null(cstr$Gammat) && !is.character(cstr$Gammat)) 
        cstr$Gammat = cstr$Gammat[, , ec]
      if (!is.null(cstr$Gammaw) && !is.character(cstr$Gammaw)) 
        cstr$Gammaw = cstr$Gammaw[, , ec]
      if (!is.null(cstr$pi) && !is.character(cstr$pi)) 
        cstr$pi = cstr$pi[, ec]
      if (!is.null(cstr$A) && !is.character(cstr$A)) 
        cstr$A = cstr$A[, , ec]
      if (!is.null(cstr$b) && !is.character(cstr$b)) 
        cstr$b = cstr$b[, ec]
      if (!is.null(cstr$Sigma) && !is.character(cstr$Sigma)) 
        cstr$Sigma = cstr$Sigma[, , ec]
      th$c = th$c[, ec]
      th$Gamma = th$Gamma[, , ec]
      th$pi = th$pi[ec]
      th$A = th$A[, , ec]
      th$b = th$b[, ec]
      th$Sigma = th$Sigma[, , ec]
    }
    return(list(th = th, cstr = cstr))
  }
  Lt = nrow(tapp)
  L = Lt + Lw
  D = nrow(yapp)
  N = ncol(yapp)
  if (verb >= 1) {
    print("EM Initializations")
  }
  if (!is.null(in_theta)) {
    theta = in_theta
    K = length(theta$pi)
    if (is.null(cstr$cw)) {
      cstr$cw = matrix(0, L, K)
    }
    if (is.null(cstr$Gammaw)) {
      cstr$Gammaw = array(diag(Lw), dim = c(Lw, Lw, K))
    }
    tmp = ExpectationZ(tapp, yapp, theta, verb)
    r = tmp$r
    ec = tmp$ec
    tmp = remove_empty_clusters(theta, cstr, ec)
    theta = tmp$th
    cstr = tmp$cstr
    tmp = ExpectationW(tapp, yapp, theta, verb)
    muw = tmp$muw
    Sw = tmp$Sw
    if (verb >= 1) 
      print("")
  }
  else {
    if (is.null(in_r)) {
      r = emgm(rbind(tapp, yapp), in_K, 1000, verb = verb)$R
    }
    else {
      r = in_r$R
    }
    if (Lw == 0) {
      Sw = NULL
      muw = NULL
    }
    else {
      theta = Maximization(tapp, yapp, r, NULL, NULL, cstr, 
                           verb)
      K = length(theta$pi)
      if (is.null(cstr$cw)) {
        cstr$cw = matrix(0, Lw, K)
      }
      theta$c = rbind(theta$c, cstr$cw[, 1:K])
      Gammaf = array(0, dim = c(L, L, K))
      Gammaf[1:Lt, 1:Lt, ] = theta$Gamma
      if (is.null(cstr$Gammaw)) {
        cstr$Gammaw = array(diag(Lw), dim = c(Lw, Lw, 
                                              K))
      }
      Gammaf[(Lt + 1):L, (Lt + 1):L, ] = cstr$Gammaw[, 
                                                     , 1:K]
      theta$Gamma = Gammaf
      Aw = array(0, dim = c(D, Lw, K))
      for (k in 1:K) {
        rk_bar = sum(r[, k])
        bk = theta$b[, k]
        w = sweep(yapp, 1, bk, "-")
        if (Lt > 0) {
          Ak = theta$A[, , k]
          w = w - Ak %*% tapp
        }
        w = sweep(w, 2, sqrt(r[, k]/rk_bar), "*")
        C = tcrossprod(w)
        tmp = eigen(C)
        U = tmp$vectors[, 1:Lw]
        Lambda = tmp$values[1:Lw]
        sigma2k = (sum(diag(C)) - sum(Lambda))/(D - Lw)
        theta$Sigma[, , k] = sigma2k * diag(D)
        Aw[, , k] = U %*% sqrt(diag(Lambda, ncol = length(Lambda), 
                                    nrow = length(Lambda)) - sigma2k * diag(Lw))
      }
      theta$A = abind(theta$A, Aw, along = 2)
      tmp = ExpectationZ(tapp, yapp, theta, verb)
      r = tmp$r
      ec = tmp$ec
      tmp = remove_empty_clusters(theta, cstr, ec)
      theta = tmp$th
      cstr = tmp$cstr
      tmp = ExpectationW(tapp, yapp, theta, verb)
      muw = tmp$muw
      Sw = tmp$Sw
      if (verb >= 1) 
        print("")
    }
  }
  if (verb >= 1) 
    print("       Running EM")
  LL = rep(-Inf, maxiter)
  iter = 0
  converged = FALSE
  while (!converged & iter < maxiter) {
    iter = iter + 1
    if (verb >= 1) 
      print(paste("       Iteration ", iter, sep = ""))
    theta = Maximization(tapp, yapp, r, muw, Sw, cstr, verb, 
                         ...)
    tmp = ExpectationZ(tapp, yapp, theta, verb)
    r = tmp$r
    LL[iter] = tmp$LL
    if (verb >= 1) {
      print(LL[iter])
    }
    ec = tmp$ec
    tmp = remove_empty_clusters(theta, cstr, ec)
    theta = tmp$th
    cstr = tmp$cstr
    tmp = ExpectationW(tapp, yapp, theta, verb)
    muw = tmp$muw
    Sw = tmp$Sw
    if (iter >= min_iter && iter > patience) {  # "patience" means that several consecutive small improvements are required, instead of just one. It avoids stopping because of a single flat step.
      recent_delta <- diff(LL[(iter - patience):iter])
      recent_rel <- recent_delta / (abs(LL[(iter - patience):(iter - 1)]) + 1)
      
      converged <- all(is.finite(recent_rel)) &&
        all(recent_delta >= -1e-8) &&
        all(recent_rel < em_tol)
    }
  # In the original gllim function the stopping criterion followed the approach below 
  #  deltaLL_total <- max(LL[1:iter]) - min(LL[1:iter])
  #  deltaLL <- LL[iter] - LL[iter - 1]
  #  converged <- deltaLL <= 0.001 * deltaLL_total
    
    if (verb >= 1) 
      print("")
  }
  LLf = LL[iter]
  if (verb >= 1) 
    print(paste("Converged in ", iter, " iterations", sep = ""))
  theta$r = r
  theta$LLf = LLf
  theta$LL = LL[1:iter]
  if (cstr$Sigma == "i") {
    nbparSigma = 1
  }
  if (cstr$Sigma == "d") {
    nbparSigma = D
  }
  if (cstr$Sigma == "") {
    nbparSigma = D * (D + 1)/2
  }
  if (cstr$Sigma == "*") {
    nbparSigma = D * (D + 1)/(2 * in_K)
  }
  if (cstr$Sigma == "bSHOCK") {
    nbparSigma = Inf
  }
  if (!is.null(cstr$Gammat)) {
    if (cstr$Gammat == "i") {
      nbparGamma = 1
    }
    if (cstr$Gammat == "d") {
      nbparGamma = Lt
    }
    if (cstr$Gammat == "") {
      nbparGamma = Lt * (Lt + 1)/2
    }
    if (cstr$Gammat == "*") {
      nbparGamma = Lt * (Lt + 1)/(2 * in_K)
    }
  }
  if (is.null(cstr$Gammat)) {
    nbparGamma = Lt * (Lt + 1)/2
  }
  theta$nbpar = (in_K - 1) + in_K * (D * L + D + Lt + nbparSigma + 
                                       nbparGamma)
  return(theta)
}