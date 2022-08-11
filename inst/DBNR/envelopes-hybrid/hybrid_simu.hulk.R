simu.hulk <- function(m, 
                      s = 10, 
                      K1 = floor(m/(s * 4)), 
                      d = 1, 
                      barmu = 4,
                      setting = "const",
                      methods_ = c("tree", "part", "Simes", "Oracle", "hyb1", "hyb2", "hyb3"),
                      rho = 0, 
                      grouped = FALSE, 
                      random.leaves = FALSE, 
                      alphas = 0.05,
                      gamma = 0.02,
                      verbose = FALSE) {
    methods_ <- match.arg(methods_, several.ok = TRUE)
    dd <- dyadic.from.window.size(m, s, method = 2)
    leaf_list <- dd$leaf_list
    mu <- gen.mu.leaves2(m, K1, d, grouped, setting, barmu, leaf_list)
    # m1 <- sum(mu > 0)
    m1 <-  floor(d*s)*K1 ## covers the case where barmu==0
    # xmax <- min(2 * m1, m)
    xmax <- min(4/3 * m1, m)
    #idxs1 <- c(1:max(xmax, s))
    idxs1 <- round(seq(from = 1, to = xmax, length = 10))
    idxs2 <- round(seq(from = xmax + 1, to = m, length = 10))
    idxs <- c(idxs1, idxs2)
    
    pvalues <- gen.p.values(m, mu, rho)
    #hf <- cherry::hommelFast(pvalues)
    C <- dd$C
    Cs <- list(tree = C,
               part = C[length(C)],
               hyb1 = C, 
               hyb2 = C,
               hyb3 = C)
    rm(C)
    
    # orders <- list(p.value = order(pvalues),
    #           mu = order(mu, decreasing = TRUE))
    orders <- list(p.value = order(pvalues))
    configs <- expand.grid(order = names(orders), method = methods_, alpha = alphas)
    resList <- list()
    for (cc in 1:nrow(configs)) {
        if (verbose) print(configs[cc, ])
        ord <- configs[cc, "order"]
        oo <- orders[[ord]]
        meth <- configs[cc, "method"]
        alpha <- configs[cc, "alpha"]
        #print(configs[cc, ])
        if (verbose) print(meth)
        if (meth == "Oracle") {
            H0 <- which(mu == 0)  ## formally not true when barmu is 0...
            V <- cumsum(oo %in% H0)
            V <- V[idxs]
        } else if (meth == "Simes") {
            thr <- alpha * 1:m/m
            FP_Simes <- sanssouci:::curveMaxFP(pvalues, thr)
            V <- FP_Simes[idxs]
            # V2 <- idxs - sapply(idxs, FUN = function(ii) {
            #     posthocBySimes(pvalues, oo[1:ii], alpha) ## slow!
            # })
            # stopifnot(identical(V, V2))
        } else if (meth %in% c("tree", "part")) {
            ZL <- zetas.tree(Cs[[meth]], leaf_list, zeta.DKWM, pvalues, alpha = alpha)
            V <- sapply(idxs, FUN = function(ii) {
                V.star(oo[1:ii], Cs[[meth]], ZL, leaf_list)
            })
        # dans les lignes de hyb1 et hyb2, on a remplacé les Cs[[meth]] par dd$C car erreur sinon
        # #ligne de test
        # } else if (meth == "hyb1") {  
        #   thr <- alpha * 1:m/m
        #   FP_Simes <- sanssouci:::curveMaxFP(pvalues, thr)
        #   V <- FP_Simes[idxs]
        } else if (meth == "hyb1") {
          thr <- (1 - gamma) * alpha * 1:m/m
          FP_Simes <- sanssouci:::curveMaxFP(pvalues, thr)
          V1 <- FP_Simes[idxs]
          ZL <- zetas.tree(Cs[["hyb1"]], leaf_list, zeta.DKWM, pvalues, alpha = gamma * alpha)
          V2 <- sapply(idxs, FUN = function(ii) {
            V.star(oo[1:ii], Cs[["hyb1"]], ZL, leaf_list)
          })
          V <- pmin(V1, V2)
          
        } else if (meth == "hyb2") {
          ZL_DKWM <- zetas.tree(Cs[["hyb2"]], leaf_list, zeta.DKWM, pvalues, alpha = gamma * alpha)
          ZL_HB <- zetas.tree(Cs[["hyb2"]], leaf_list, zeta.HB, pvalues, alpha = (1 - gamma) * alpha)
          ZL = list()
          for (i in 1:length(ZL_HB)){
            ZL[[i]] <- pmin(ZL_DKWM[[i]], ZL_HB[[i]])
          }
          V <- sapply(idxs, FUN = function(ii) {
            V.star(oo[1:ii], Cs[["hyb2"]], ZL, leaf_list)
          })
        } else if (meth == "hyb3") {
          ZL <- zetas.tree(Cs[["hyb3"]], leaf_list, zeta.DKWM, pvalues, alpha = gamma * alpha)
          V1 <- sapply(idxs, FUN = function(ii) {
            V.star(oo[1:ii], Cs[["hyb3"]], ZL, leaf_list)
          })
          ZL <- zetas.tree(Cs[["hyb3"]], leaf_list, zeta.HB, pvalues, alpha = (1 - gamma) * alpha)
          V2 <- sapply(idxs, FUN = function(ii) {
            V.star(oo[1:ii], Cs[["hyb3"]], ZL, leaf_list)
          })
          V <- pmin(V1, V2)
          
        }
         
        res <- data.frame(idxs, value = V, method = meth, order = ord, alpha = alpha)
        resList[[cc]] <- res
    }
    resList
}
