################################
# Author: Bulent Ozel
# e-mail: bulent.ozel@gmail.com
# Collaborators:
#	Mario Eboli
#	Andrea Teglio
#	Andrea Toto
#################################

#################################
# Importnat Note: In order to use these drivers ModelLibrary.R needs to be loaded.
#################################


# Simulation Templates :::::::::::::::


# b.sheet.equity can be one of:
# - "homo", denoting that each bank has the same and equal equity, which is E0/N.
# - "hetero", denoting that equity of the bank depends on phi, eps and network position.
# - "planA.spec", holds only when K = 1, for K > 1 and it switches to "homo" configuration.

# shock.type can be one of:
# - "assets", the set of external assets is permuted for a shock sequence.
#             Each bank owns a distinct subset of assets.
# - "banks", the set of banks in the system is permuted for a shock sequence.
#             Set of all assets of a bank which is picked is incrementally hit before moving to the
#             next bank in the sequence.
# - "full", the set of banks in the system is permuted for a shock sequence.
#             All assets of a bank which is picked all at once.

connectivity.and.centralization <- function(N,
                                            KMin = 1,
                                            KMax = N,
                                            rCentMin = 0,
                                            rCentMax = N,
                                            symmetric = 1,
                                            eps = 0.1,
                                            phi = 0.4,
                                            A0 = 100,
                                            E0 = 10,
                                            b.sheet.equity = "hetero",
                                            d.min = 1,
                                            nRuns = N,
                                            shock.type = "banks",
                                            unit.shock = 1.0,
                                            r.shock = 1.0,
                                            n.sources = N) {
  for (K in KMin:KMax) {
    for (c in rCentMin:rCentMax) {
      cent.level = c / N
      
      ### 1: Configure and create the network structure
      
      ##: First a K regular graph is cretaed. Then keeping the connectivity constant
      ##: the network is centralized to the desired level.
      res <- centralize(N, K, cent.level, symm = symmetric)
      net <- res$centralized
      Centralization <- net %n% "Centralization"
      Connectivity <- net %n% "Connectivity"
      
      ### 2:  Initialize balancesheets
      if (b.sheet.equity == "hetero") {
        results <- initialize.balanceSheets(net, A0, phi, eps)
      }
      else if (b.sheet.equity == "homo") {
        results <-
          initialize.balanceSheets.planA.v1(net, A0, E0, d.min)
      }
      else if (b.sheet.equity == "planA.spec") {
        if (K == 1) {
          results <- initialize.balanceSheets.planA.v2(net, A0, E0, d.min)
        }
        else {
          results <- initialize.balanceSheets.planA.v1(net, A0, E0, d.min)
        }
      }
      else {
        results <- initialize.balanceSheets(net, A0, phi, eps)
      }
      
      BS <- results[[1]]
      net <- results[[2]]
      
      defaults <- vector(mode = "numeric")
      
      
      for (run in 1:nRuns) {
        ### 3:  Initialize shock vector
        
        shocks <- initialize.shockVector(r.shock,
                                         n.sources,
                                         net,
                                         BS,
                                         shock.size = unit.shock,
                                         shock.type = shock.type)
        
        
        
        ### 4: Apply shocks and simulate the contagion
        
        results <- simulate.contagion(net, BS, shocks)
        
        ### 5:  Report results
        if (INSPECT_MODE){
          Contagion.State <- results[[1]]
        }
        Contagion.Process <- results[[2]]
        outputs <- getContagionResults(Contagion.Process, N)
        nPrimary <- outputs[[1]]
        nSecondary <- outputs[[2]]
        
        if (SCAN) {
          Tk <- outputs[[3]]
          Tf <- outputs[[4]]
          outputs <-
            as.vector(c(
              Centralization,
              Connectivity,
              nPrimary,
              nSecondary,
              Tf,
              Tk
            ))
        }
        else {
          T1 <- outputs[[3]]
          T2 <- outputs[[4]]
          dT <- outputs[[5]]
          outputs <-
            as.vector(c(
              dT,
              T1,
              T2,
              Centralization,
              Connectivity,
              nPrimary,
              nSecondary
            ))
        }
        defaults <- c(defaults, outputs)
      }
      
      ### 5:  Report results
      if (SCAN) {
        Out.m <- matrix(defaults, ncol = (5 + N), byrow = TRUE)
        Out <- as.data.frame(Out.m)
        ncolnames <-
          c("Centralization",
            "Connectivity",
            "nPrimary",
            "nSecondary",
            "Tfin")
        for (i in 1:N) {
          ncolnames <- c(ncolnames, paste("T", i, sep = ""))
        }
        colnames(Out) <- ncolnames
      }
      else {
        Out.m <- matrix(defaults, ncol = 7, byrow = TRUE)
        Out <- as.data.frame(Out.m)
        colnames(Out) <-
          c("dT",
            "T1",
            "T2",
            "Centralization",
            "Connectivity",
            "nPrimary",
            "nSecondary")
        #assign(paste("Out","_N",N,"_K",K,"_c",c, sep=""), Out)
      }
      fname = paste("Out", "_N", N, "_K", K, "_c", c, sep = "")
      cap.text = paste("Results.")
      create.Table(Out,
                   name = fname,
                   type = "csv",
                   cap = cap.text)
      
      if (INSPECT_MODE){
        fname = paste("BS", "_N", N, "_K", K, "_c", c, sep = "")
        cap.text = paste("Balance Sheets.")
        create.Table(BS,
                     name = fname,
                     type = "csv",
                     cap = cap.text)
        
        fname = paste("Net", "_N", N, "_K", K, "_c", c, sep = "")
        plot.Network.v2(net, f.name = fname)
        
        fname = paste("Diffusion", "_N", N, "_K", K, "_c", c, sep = "")
        cap.text = paste("Final diffusion state of the latest run.")
        create.Table(Contagion.State,
                     name = fname,
                     type = "csv",
                     cap = cap.text)
      }
      
    }
  }
}