################################
# Author: Bulent Ozel
# e-mail: bulent.ozel@gmail.com
# Collaborators:
#	Mario Eboli
#	Andrea Teglio
#	Andrea Toto
#################################

# Used R libraries:::::::::::::::
require(sna)
#require(network)
library("network")
require(xtable)

# Newly implemennted tools:::::::

initialize.balanceSheets <- function(net,
                                     A0 = 100,
                                     phi = 0.4,
                                     eps = 0.1) {
  # Get banks summary:
  Bank <- net %v% "vertex.names"
  n.banks <- net %n% "n"
  
  # Degree of a node is total in and out degrees.
  # The links are directional and the network is symmetric:
  degrees <- net %v% "deg"
  degrees <- degrees / 2
  
  # Get number of node level connections:
  #nlinks <- net$gal$mnext - 1
  nlinks <- sum(net[, ])
  # ... the links are directional.
  
  # Determine the link weight
  H0 <- A0 * (1 - eps) * (1 + eps * phi) ^ -1
  D0 <- H0 * phi
  C0 <- D0
  E0 <- eps * (A0 + C0)
  w0 <- D0 / nlinks
  
  # Rest of balancesheet values will be driven from connectivity:
  d0 <- w0
  a0 <- A0 / nlinks
  h0 <- H0 / nlinks
  e0 <- E0 / nlinks
  
  zeros <- rep(0, n.banks)
  A <- zeros
  H <- zeros
  E <- zeros
  C <- zeros
  D <- zeros
  Balance.Sheets <- data.frame(Bank, A, C, H, D, E)
  
  for (i in 1:n.banks) {
    deg <- degrees[i]
    # Initialize the isolate
    if (deg == 0) {
      Balance.Sheets$D[i] <- 0
      Balance.Sheets$C[i] <- 0
      a <- (1 / n.banks) * H0 / (1 - eps)
      Balance.Sheets$A[i] <- a
      Balance.Sheets$E[i] <- eps * a
      Balance.Sheets$H[i] <- (1 - eps) * a
      next
    }
    
    # Initialize D
    d <- d0 * deg
    Balance.Sheets$D[i] <- d
    
    # Initialize C
    c <- d
    Balance.Sheets$C[i] <- c
    
    # Fixed ratios case initialization:
    # initialize H
    
    h <- h0 * deg
    Balance.Sheets$H[i] <- h
    
    # initialize A
    #a <- (1 + phi * eps)/(phi * (1 - eps)) * d
    a <- a0 * deg
    Balance.Sheets$A[i] <- a
    
    # initialize E
    #e <- eps * (a + c)
    e <- e0 * deg
    Balance.Sheets$E[i] <- e
  }
  
  # Return balance sheet
  net %n% "link.weight" <- w0
  return(list(Balance.Sheets, net))
}

initialize.shockVector <- function(shock.rate,
                                   n.sources,
                                   net,
                                   Balance.Sheets = NULL,
                                   shock.size = 1.0,
                                   shock.type = "assets") {
  bank.ids <- net %v% "vertex.names"
  n.banks <- net %n% "n"
  
  # In following shock schema a rondomly selected n.sources of nodes are selected.
  # Each one of them will have an external shock rate of equal to shock.rate.
  # The shock rate is the ratio of total external rates of a bank that is hit.
  # A list data structure of shocked nodes and shock rates is returned.
  if (shock.type == "full") {
    sequence <- sample(bank.ids, n.sources)
    rates <- rep(shock.rate, length(sequence))
    return(list(sequence, rates))
  }
  
  
  # For the followig two shock schemas a unit amount of shock should be given.
  # The input parameter for the unit is shock.size.
  # The balance sheets of the banks should be given.
  
  
  # In the default choice a random sequence of shocks on the banks are selected.
  # A selected bank is hit repeatedly until all of its external assets are removed.
  A <- Balance.Sheets$A
  sequence <- sample(bank.ids, n.banks)
  sequence.org <- sequence
  remainders <- sapply(sequence, function(i) {
    A[i] %% shock.size
  })
  repitions <- sapply(sequence, function(i) {
    A[i] %/% shock.size
  })
  repitions <-
    repitions + sapply(remainders, function(i) {
      as.integer(i > 0)
    })
  sequence <- rep(sequence, times = repitions)
  rates <- sapply(sequence, function(i) {
    shock.size / A[i]
  })
  
  
  for (i in which(remainders > 0)) {
    ind <- sequence.org[i]
    rates[which(sequence == ind)[1]] = (remainders[i] / A[ind])
  }
  
  # When each asset is considered as a stand alone project, a random sequence of shocks
  # on the list of external assets is created.
  # It is assumed that assets are scattered randomly in the system.
  if (shock.type == "assets") {
    nsequence <- length(sequence)
    indices <- sample(1:nsequence, nsequence)
    
    rates <- sapply(indices, function(i) {
      rates[i]
    })
    sequence <- sapply(indices, function(i) {
      sequence[i]
    })
  }
  
  # A list of two vectors is returned. These parallel vectors contain the id of the bank
  # and the rate of shock received by the bank.
  return(list(sequence, rates))
}


generate.RandomNetwork <- function(N, p.link) {
  g <- rgraph(
    n = N,
    m = 1,
    tprob = p.link,
    mode = "graph",
    diag = F
  )
  g.net <- as.network(g)
  return(g.net)
}


generate.RegularGraph <- function(N, r = 1) {
  g <- rgraph(
    n = N,
    m = 1,
    tprob = 0,
    mode = "digraph",
    diag = F
  )
  net <- as.network(g)
  
  for (from in 1:N) {
    for (j in 1:r) {
      to <- (from - 1 + j) %% N
      to <- to + 1
      net[from, to] <- 1
    }
  }
  
  # Collect Information About Generated Network
  
  Degree.Distribution <- degree(net, gmode = "digraph")
  Centralization <- centralization(net, degree)
  Connectivity <- gden(net)
  
  net %v% "deg" <- Degree.Distribution
  net %n% "Centralization" <- Centralization
  net %n% "Connectivity" <- Connectivity
  
  return(net)
}

generate.CorePeriphery <- function(net, k, n, p.link) {
  if (k >= n) {
    return(net)
  }
  for (i in seq(from = k + 1, to = n, by = 1)) {
    if (runif(1, 0, 1) > p.link) {
      next
    }
    to <- i %% k + 1
    net[i, to] <- 1
    net[to, i] <- 1
  }
  return(net)
}


plot.Network.v2 <-
  function(a.net,
           n.sizes = NULL,
           n.labels = NULL,
           f.name = NULL,
           insolvency = NULL,
           shock.vect = NULL,
           isSecondary = NULL,
           v.coords = NULL) {
    n.colors <- rep("green", a.net$gal$n)
    v.shapes <- rep(8, a.net$gal$n)
    
    
    if (!is.null(insolvency)) {
      for (i in 1:length(n.colors)) {
        if (insolvency[i] < 0.66) {
          n.colors[i] <- "green"
        } else if (insolvency[i] < 1) {
          n.colors[i] <- "grey"
        } else {
          n.colors[i] <- "red"
        }
      }
    }
    
    
    if (!is.null(isSecondary)) {
      for (i in 1:length(n.colors)) {
        if (isSecondary[i] > 0) {
          n.colors[i] <- "black"
        }
      }
    }
    
    
    if (!is.null(shock.vect)) {
      # mark the node that is hit with ext. shock.
      for (i in 1:length(v.shapes)) {
        if (shock.vect[i] > 0) {
          v.shapes[i] <- 3
        }
      }
    }
    
    if (is.null(n.sizes)) {
      degrees <- degree(a.net, gmode = "graph")
      n.sizes <- (degrees + 1) ^ 0.5
    }
    
    if (is.null(n.labels)) {
      n.labels <- network.vertex.names(a.net)
    }
    
    if (is.null(f.name)) {
      dev.new()
      if (length(n.sizes) < 2) {
        if (is.null(v.coords)) {
          v.coords <-
            gplot(
              a.net,
              vertex.cex = n.sizes,
              label.cex = 0.7,
              vertex.col = n.colors,
              vertex.sides = v.shapes
            )
        } else {
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label.cex = 0.7,
            vertex.col = n.colors,
            vertex.sides = v.shapes,
            coord = v.coords
          )
        }
        return(v.coords)
      }
      if (is.null(v.coords)) {
        v.coords <-
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label = n.labels,
            label.cex = 0.7,
            vertex.col = n.colors,
            vertex.sides = v.shapes
          )
      } else {
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label = n.labels,
          label.cex = 0.7,
          vertex.col = n.colors,
          vertex.sides = v.shapes,
          coord = v.coords
        )
      }
      return(v.coords)
    }
    
    pngfile <- paste0(PLOTS_DIR, f.name, ".png")
    png(pngfile,
        width = 1000,
        height = 1000,
        pointsize = 18)
    if (length(n.sizes) < 2) {
      if (is.null(v.coords)) {
        v.coords <-
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label.cex = 0.7,
            vertex.col = n.colors,
            vertex.sides = v.shapes
          )
      } else {
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label.cex = 0.7,
          vertex.col = n.colors,
          vertex.sides = v.shapes,
          coord = v.coords
        )
      }
      dev.off()
      return(v.coords)
    }
    if (is.null(v.coords)) {
      v.coords <-
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label = n.labels,
          label.cex = 0.7,
          vertex.col = n.colors,
          vertex.sides = v.shapes
        )
    } else {
      gplot(
        a.net,
        vertex.cex = n.sizes,
        label = n.labels,
        label.cex = 0.7,
        vertex.col = n.colors,
        vertex.sides = v.shapes,
        coord = v.coords
      )
    }
    dev.off()
    return(v.coords)
  }

plot.Network <-
  function(a.net,
           n.sizes = NULL,
           n.labels = NULL,
           f.name = NULL,
           insolvency = NULL,
           shock.vect = NULL,
           v.coords = NULL) {
    n.colors <- rep("green", a.net$gal$n)
    
    if (!is.null(insolvency)) {
      for (i in 1:length(n.colors)) {
        if (insolvency[i] < 0.66) {
          n.colors[i] <- "green"
        } else if (insolvency[i] < 1) {
          n.colors[i] <- "grey"
        } else {
          n.colors[i] <- "black"
        }
      }
    }
    
    if (!is.null(shock.vect)) {
      # paint sources:
      for (i in 1:length(n.colors)) {
        if (shock.vect[i] > 0) {
          n.colors[i] <- "yellow"
        }
      }
      # paint primary defaults.
      if (!is.null(insolvency)) {
        m <- insolvency * shock.vect
        for (i in 1:length(n.colors)) {
          if (m[i] == 1) {
            n.colors[i] <- "red"
          }
        }
      }
    }
    
    if (is.null(n.sizes)) {
      degrees <- degree(a.net, gmode = "graph")
      n.sizes <- (degrees + 1) ^ 0.5
    }
    
    if (is.null(n.labels)) {
      n.labels <- network.vertex.names(a.net)
    }
    
    if (is.null(f.name)) {
      dev.new()
      if (length(n.sizes) < 2) {
        if (is.null(v.coords)) {
          v.coords <-
            gplot(
              a.net,
              vertex.cex = n.sizes,
              label.cex = 0.7,
              vertex.col = n.colors
            )
        } else {
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label.cex = 0.7,
            vertex.col = n.colors,
            coord = v.coords
          )
        }
        return(v.coords)
      }
      if (is.null(v.coords)) {
        v.coords <-
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label = n.labels,
            label.cex = 0.7,
            vertex.col = n.colors
          )
      } else {
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label = n.labels,
          label.cex = 0.7,
          vertex.col = n.colors,
          coord = v.coords
        )
      }
      return(v.coords)
    }
    
    pngfile <- paste0(PLOTS_DIR, f.name, ".png")
    png(pngfile,
        width = 1000,
        height = 1000,
        pointsize = 18)
    if (length(n.sizes) < 2) {
      if (is.null(v.coords)) {
        v.coords <-
          gplot(
            a.net,
            vertex.cex = n.sizes,
            label.cex = 0.7,
            vertex.col = n.colors
          )
      } else {
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label.cex = 0.7,
          vertex.col = n.colors,
          coord = v.coords
        )
      }
      dev.off()
      return(v.coords)
    }
    if (is.null(v.coords)) {
      v.coords <-
        gplot(
          a.net,
          vertex.cex = n.sizes,
          label = n.labels,
          label.cex = 0.7,
          vertex.col = n.colors
        )
    } else {
      gplot(
        a.net,
        vertex.cex = n.sizes,
        label = n.labels,
        label.cex = 0.7,
        vertex.col = n.colors,
        coord = v.coords
      )
    }
    dev.off()
    return(v.coords)
  }

plot.DegreeDistro <-
  function(datavector,
           xlabel,
           ylabel,
           title,
           f.name,
           isgrid = FALSE,
           ltype = "l",
           ispoint = F) {
    pngfile <- paste0(PLOTS_DIR, f.name, ".png")
    png(pngfile,
        width = 800,
        height = 600,
        pointsize = 18)
    times <- (1:length(datavector))
    plot(
      datavector ~ times,
      type = "n",
      xlab = xlabel,
      ylab = ylabel,
      main = title
    )
    if (isgrid) {
      abline(
        v = times,
        h = datavector,
        col = "lightgray",
        lty = "dotted",
        lwd = 2
      )
    }
    lines(datavector ~ times,
          type = ltype,
          col = "blue",
          lwd = 40)
    if (ispoint) {
      points(datavector ~ times,
             col = "red",
             pch = 16,
             cex = 2)
    }
    dev.off()
  }


create.Table <-
  function(x,
           name = NULL,
           type = NULL,
           cap = "",
           lab = "") {
    if (is.null(name)) {
      name = deparse(substitute(x))
    }
    csvfile <- paste0(TABLES_DIR, name, ".csv")
    write.table(
      x,
      file = csvfile,
      sep = ",",
      col.names = NA,
      qmethod = "double"
    )
    if (DATA_MODE) {
      return(TRUE)
    }
    pdf_option <- paste0(" -output-directory ", TABLES_DIR, " ")
    texfile <- paste0(TABLES_DIR, name, ".tex")
    tt <- print(xtable(x, caption = cap, label = lab), type = "latex")
    if (is.null(type)) {
      cat(tt, sep = "", file = texfile)
      return(TRUE)
    }
    cat(
      "\\documentclass[12pt]{report}\n  \\usepackage[landscape]{geometry}\n  \\date{}\n  \\begin{document}",
      tt,
      "\\end{document}",
      sep = "",
      file = texfile
    )
    system(paste0(LATEX_COMMAND, pdf_option, texfile))
  }


build.CorePeriphery <-
  function(N = 30,
           K = 3,
           p.periphery = 0.02,
           p.core = 0.8,
           p.to.core = 0.8) {
    # Generate Base
    net.base <- generate.RandomNetwork(N, p.periphery)
    
    # Generate Core
    net.core <- generate.RandomNetwork(K, p.core)
    
    # Merge Base and Core
    new.edges <- as.edgelist.sna(net.core)
    net.all <- network.edgelist(new.edges, net.base)
    
    # Connect Periphery to the Core
    net.cp <- generate.CorePeriphery(net.all, K, N, p.to.core)
    
    # Collect Information About Generated Network
    
    Degree.Distribution <- degree(net.cp, gmode = "digraph")
    Centralization <- centralization(net.cp, degree)
    Connectivity <- gden(net.cp)
    
    net.cp %v% "deg" <- Degree.Distribution
    net.cp %n% "Centralization" <- Centralization
    net.cp %n% "Connectivity" <- Connectivity
    
    return(net.cp)
  }


get.debtorLists <- function(net) {
  # It returns list of banks who had taken
  # loans from the bank i. They are the channels of receiving
  # interbank shocks.
  m <- net[, ]
  n <- net %n% "n"
  edge.lists <- list()
  for (i in 1:n) {
    edges <- which(m[, i] > 0)
    edge.lists[[i]] <- as.integer(names(edges))
  }
  return(edge.lists)
}



apply.incremental.shocks <- function(rate, bs, state) {
  # a simple fixed rate all nodes shock.
  shock <- rep(rate, length(bs$Bank))
  return(shock)
}


apply.initialShock <- function(rate, bs, sources = NULL) {
  # a simple fixed rate initial shock.
  shock <- rep(0, length(bs$Bank))
  if (is.null(sources)) {
    shock <- rate
    return(shock)
  }
  
  for (i in 1:length(sources)) {
    index <- sources[i]
    shock[index] <- rate
  }
  
  return(shock)
}


update.internalShocks <- function(Contagion.State, Debtors, BW) {
  # M is the matrix with link weights
  # Debtors is the edge list of bank ids that has taken loans from the key
  # of the list.
  # BW is either a scalar or matrix. When it is matrix it reprsents
  # that the links in the network are not equal.
  
  n <- length(Contagion.State$Bank)
  for (i in 1:n) {
    c.i <- Debtors[[i]]
    m <- length(c.i)
    acc <- 0
    if (m == 0) {
      Contagion.State$r.Shock.Int[i] <- acc
      next
    }
    
    if (is.matrix(BW)) {
      for (j in 1:m) {
        ind <- c.i[j]
        acc <- acc + Contagion.State$r.Flow.Out[ind] * BW[ind, i]
      }
      Contagion.State$r.Shock.Int[i] <- acc
    }
    else {
      for (j in 1:m) {
        ind <- c.i[j]
        acc <- acc + Contagion.State$r.Flow.Out[ind]
      }
      Contagion.State$r.Shock.Int[i] <- acc * BW
    }
  }
  return(Contagion.State)
}

update.totalShocks <- function(Balance.Sheets, Contagion.State) {
  n <- length(Balance.Sheets$Bank)
  for (i in 1:n) {
    r.ext <- Contagion.State$r.Shock.Ext[i]
    total.int <- Contagion.State$r.Shock.Int[i]
    total.ext <- r.ext * Balance.Sheets$A[i]
    Contagion.State$Shock[i] <- total.ext + total.int
  }
  return(Contagion.State)
}

update.totalAbsorbtion <-
  function(Balance.Sheets, Contagion.State) {
    n <- length(Balance.Sheets$Bank)
    rate <- 0
    for (i in 1:n) {
      e <- Balance.Sheets$E[i]
      Shock <- Contagion.State$Shock[i]
      rate <- Shock / e
      Contagion.State$r.Absorb.E[i] <- min(1, rate)
    }
    return(Contagion.State)
  }

update.shockDiffusion <- function(Balance.Sheets, Contagion.State) {
  n <- length(Balance.Sheets$Bank)
  rate <- 0
  for (i in 1:n) {
    E <- Balance.Sheets$E[i]
    H <- Balance.Sheets$H[i]
    D <- Balance.Sheets$D[i]
    Lambda <- Contagion.State$Shock[i]
    rate <- (Lambda - E) / (H + D)
    Contagion.State$r.Flow.Out[i] <- max(0, rate)
  }
  return(Contagion.State)
}

is.draining <- function(Balance.Sheets, Contagion.State) {
  n <- length(Balance.Sheets$Bank)
  acc.shock <- 0
  acc.absorb <- 0
  acc.drain <- 0
  for (i in 1:n) {
    acc.shock <-
      acc.shock + Balance.Sheets$A[i] * Contagion.State$r.Shock.Ext[i]
    acc.absorb <-
      acc.absorb + Balance.Sheets$E[i] * Contagion.State$r.Absorb.E[i]
    acc.drain <-
      acc.drain + Balance.Sheets$H[i] * Contagion.State$r.Flow.Out[i]
  }
  return(acc.shock - (acc.absorb + acc.drain))
}



initialize.balanceSheets.planA.v1 <-
  function(net, A0, E0, d.min, fact.cent = 0) {
    #In this set-up:
    # Total values A0,E0 and hence a,h,e are constant. d will cahange depending on the centralization.
    
    
    # Get banks summary:
    Bank <- net %v% "vertex.names"
    n.banks <- net %n% "n"
    degrees <- net %v% "deg"
    
    centralization <- net %n% "Centralization"
    
    # The exposure to interbank incerases with the centralization.
    
    d.min <- (1 + centralization * fact.cent) * d.min
    
    # Get number of node level connections:
    #nlinks <- net$gal$mnext - 1
    nlinks <- sum(net[, ])
    # ... the links are directional.
    
    zeros <- rep(0, n.banks)
    A <- zeros
    H <- zeros
    E <- zeros
    C <- zeros
    D <- zeros
    
    Balance.Sheets <- data.frame(Bank, A, C, H, D, E)
    
    H0 <- A0 - E0
    a <- A0 / n.banks
    e <- E0 / n.banks
    h <- H0 / n.banks
    
    for (i in 1:n.banks) {
      Balance.Sheets$A[i] <- a
      Balance.Sheets$E[i] <- e
      Balance.Sheets$H[i] <- h
      
      deg <- degrees[i]
      deg <- deg / 2
      
      
      # Initialize the isolate
      if (deg == 0) {
        Balance.Sheets$D[i] <- 0
        Balance.Sheets$C[i] <- 0
        next
      }
      
      interbank.debt <- deg * d.min
      Balance.Sheets$D[i] <- interbank.debt
      Balance.Sheets$C[i] <- interbank.debt
      
      # if (deg > 2 * k) {
      # interbank.debt <- deg * d.min
      # Balance.Sheets$D[i] <- interbank.debt
      # Balance.Sheets$C[i] <- interbank.debt
      # }
      # else {
      # Balance.Sheets$D[i] <- d.min
      # Balance.Sheets$C[i] <- d.min
      # }
    }
    
    # Return balance sheet
    net %n% "link.weight" <- d.min
    net %n% "min.link.weight" <- d.min
    
    return(list(Balance.Sheets, net))
  }

initialize.balanceSheets.planA.v2 <-
  function(net, A0, E0, d.min) {
    # In this set-up:
    # Total values A0,E0 and hence a,h,e are constant. d will cahange depending on the connectivity.
    # The d of a node on the cycle or a pendant node doesn't change.
    
    ### Internal link weight function:
    assign.link.weights.PlanA <- function(net, d.min, degrees) {
      m <- net[, ]
      m[, ] <- 0
      for (r in 1:(N - 1)) {
        for (c in (r + 1):N) {
          if (net[r, c] == 0) {
            next
          }
          if ((degrees[r] == 4) | (degrees[c] == 4)) {
            m[r, c] <- d.min / 2
            m[c, r] <- d.min / 2
          }
          else {
            m[r, c] <- d.min
            m[c, r] <- d.min
          }
        }
      }
      net %n% "link.weight" <- m
      return(net)
    }
    ###############
    
    # Get banks summary:
    Bank <- net %v% "vertex.names"
    n.banks <- net %n% "n"
    degrees <- net %v% "deg"
    
    centralization <- net %n% "Centralization"
    
    # The exposure to interbank incerases with the centralization.
    
    # Get number of node level connections:
    # nlinks <- sum(net[,])
    net  <- assign.link.weights.PlanA(net, d.min, degrees)
    IB.Exp <- net %n% "link.weight"
    
    zeros <- rep(0, n.banks)
    A <- zeros
    H <- zeros
    E <- zeros
    C <- zeros
    D <- zeros
    
    Balance.Sheets <- data.frame(Bank, A, C, H, D, E)
    
    H0 <- A0 - E0
    a <- A0 / n.banks
    e <- E0 / n.banks
    h <- H0 / n.banks
    
    for (i in 1:n.banks) {
      Balance.Sheets$C[i] <- sum(IB.Exp[, i])
      Balance.Sheets$D[i] <- sum(IB.Exp[i, ])
      Balance.Sheets$A[i] <- a
      Balance.Sheets$E[i] <- e
      Balance.Sheets$H[i] <- h
    }
    
    # Return balance sheet
    return(list(Balance.Sheets, net))
  }




simulate.contagion.old <- function(net, Balance.Sheets, shocks) {
  Debtors <- get.debtorLists(net)
  w0 <- net %n% "link.weight"
  n.banks <- net %n% "n"
  
  Contagion.State <- initialize.contagionState(net)
  Shock.Nodes <- which(shocks > 0)
  Shock.Nodes <- sample(Shock.Nodes, length(Shock.Nodes))
  Shock.Rates <- rep(0, length(Shock.Nodes))
  
  for (i in 1:length(Shock.Nodes)) {
    Shock.Rates[i] <- shocks[Shock.Nodes[i]]
  }
  
  Default.Size <- rep(0, length(Shock.Nodes))
  Draining.Time <- rep(0, length(Shock.Nodes))
  Contagion.Process <-
    data.frame(Shock.Nodes, Shock.Rates, Draining.Time, Default.Size)
  
  current.shock.rates <- rep(0, n.banks)
  CS <- Contagion.State
  CP <- Contagion.Process
  BS <- Balance.Sheets
  Absorbtion.Levels <- matrix(CS$r.Absorb.E, ncol = n.banks)
  AL <- Absorbtion.Levels
  Shock.Sequence <- rep(0, n.banks)
  ASS <- matrix(Shock.Sequence, ncol = n.banks)
  for (i in 1:length(Shock.Nodes)) {
    # Stop if all nodes have defaulted
    if (all(CS$r.Absorb.E == 1)) {
      break
    }
    
    # Apply one shock at a time
    ind <- Shock.Nodes[i]
    current.shock.rates[ind] <- shocks[ind]
    CS <- apply.externalShocks(CS, current.shock.rates)
    current.shock.rates[ind] <- 0
    Shock.Sequence[ind] <- 1
    ASS <- rbind(ASS, Shock.Sequence)
    
    # Let the shock dissipate
    impact <- diffuse.shock(BS, CS, Debtors, w0)
    
    # Update the state of contagion
    CP$Draining.Time[i] <- impact[[1]]
    CS <- impact[[2]]
    CP$Default.Size[i] <- length(which(CS$r.Absorb.E == 1))
    AL <- rbind(AL, CS$r.Absorb.E)
  }
  return(list(CS, CP, AL, ASS))
}

simulate.contagion <- function(net, Balance.Sheets, shocks) {
    Debtors <- get.debtorLists(net)
    w0 <- net %n% "link.weight"
    n.banks <- net %n% "n"
    Contagion.State <- initialize.contagionState(net)
    
    Shock.Nodes <- shocks[[1]]
    Shock.Rates <- shocks[[2]]
    Default.Size <- rep(0, length(Shock.Nodes))
    Draining.Time <- rep(0, length(Shock.Nodes))
    Cumulative.Shock <- rep(0, length(Shock.Nodes))
    Secondary.Defaults <- rep(0, length(Shock.Nodes))
    Contagion.Process <-
      data.frame(
        Shock.Nodes,
        Shock.Rates,
        Draining.Time,
        Default.Size,
        Cumulative.Shock,
        Secondary.Defaults
      )
    
    CS <- Contagion.State
    CP <- Contagion.Process
    BS <- Balance.Sheets
    
    
    Absorbtion.Levels <- matrix(CS$r.Absorb.E, ncol = n.banks)
    AL <- Absorbtion.Levels
    Secondaries <- matrix(CS$is.Secondary, ncol = n.banks)
    
    current.shock.rates <- rep(0, n.banks)
    ASS <- matrix(current.shock.rates, ncol = n.banks)
    
    pre_cumulative_shock <- 0
    for (i in 1:length(Shock.Nodes)) {
      # Stop if all nodes have defaulted
      if (all(CS$r.Absorb.E == 1)) {
        break
      }
      
      # Keep the state of insolvency before the new shock.
      pre.r.Absorb.E <- CS$r.Absorb.E
      
      # Apply one shock at a time
      ind <- Shock.Nodes[i]
      current.shock.rates[ind] <- Shock.Rates[i]
      
      CS <- apply.externalShocks(CS, current.shock.rates)
      ASS <- rbind(ASS, current.shock.rates)
      current.shock.rates[ind] <- 0
      
      
      # Let the shock dissipate
      impact <- diffuse.shock(BS, CS, Debtors, w0)
      
      # Update the state of contagion
      CP$Draining.Time[i] <- impact[[1]]
      CS <- impact[[2]]
      CP$Default.Size[i] <- length(which(CS$r.Absorb.E == 1))
      size_of_new_shock <- Shock.Rates[i] * BS$A[ind]
      CP$Cumulative.Shock[i] <- pre_cumulative_shock + size_of_new_shock
      pre_cumulative_shock <- CP$Cumulative.Shock[i]
      
      for (j in 1:length(CS$r.Absorb.E)) {
        if (CS$r.Absorb.E[j] < 1) {
          next
        }
        if (pre.r.Absorb.E[j] == 1) {
          next
        }
        
        extShock = CS$r.Shock.Ext[j] * BS$A[j]
        if (extShock < BS$E[j]) {
          CS$is.Secondary[j] <- 1
        }
      }
      
      CP$Secondary.Defaults[i] <- sum(CS$is.Secondary)
      
      AL <- rbind(AL, CS$r.Absorb.E)
      Secondaries <- rbind(Secondaries, CS$is.Secondary)
      
    }
    return(list(CS, CP, AL, ASS, Secondaries))
  }


diffuse.shock <-
  function(Balance.Sheets,
           Contagion.State,
           Debtors,
           BandWidths) {
    BS <- Balance.Sheets
    CS <- Contagion.State
    BW <- BandWidths
    
    epsilon <- 0.01 * min(BW[BW > 0])
    
    t <- 0
    repeat {
      # Keep track of number of iterations
      t <- t + 1
      
      # Update contagion from t-1
      CS <- update.internalShocks(CS, Debtors, BW)
      
      # Calculate Lambda
      CS <- update.totalShocks(BS, CS)
      
      # Update beta
      CS <- update.totalAbsorbtion(BS, CS)
      
      # Update outgoing flow
      CS <- update.shockDiffusion(BS, CS)
      
      # Check Drainage
      if (is.draining(BS, CS) < epsilon) {
        break
      }
    }
    return(list(t, CS))
  }







initialize.contagionState <- function(net) {
  Bank <- net %v% "vertex.names"
  n.banks <- net %n% "n"
  
  # Initialize state of contagion:
  r.Shock.Ext <- rep(0, n.banks)
  r.Shock.Int <- rep(0, n.banks)
  r.Absorb.E <- rep(0, n.banks)
  r.Flow.Out <- rep(0, n.banks)
  Shock <- rep(0, n.banks)
  is.Secondary <- rep(0, n.banks)
  Contagion.State <-
    data.frame(Bank,
               Shock,
               r.Shock.Ext,
               r.Shock.Int,
               r.Absorb.E,
               r.Flow.Out,
               is.Secondary)
  return(Contagion.State)
}

apply.externalShocks <- function(Contagion.State, shock) {
  CS <- Contagion.State
  for (i in 1:length(shock)) {
    rate <- shock[i]
    r.ext.pre <- CS$r.Shock.Ext[i]
    CS$r.Shock.Ext[i] <- r.ext.pre + rate
  }
  return(CS)
}


run.contagion <-
  function(N,
           r.shock,
           n.sources,
           K,
           regular = TRUE,
           phi = 0.4,
           eps = 0.1,
           H0 = 100,
           pPeri = 0,
           pCore = 1,
           pToCore = 1,
           latex = NULL) {
    ### 1: Configure and create the network structure
    
    if (regular) {
      net <- generate.RegularGraph(N, K)
    } else {
      net <- build.CorePeriphery(N, K, pPeri, pCore, pToCore)
    }
    
    ### 2:  Initialize balancesheets
    results <- initialize.balanceSheets(net, H0, phi, eps)
    BS <- results[[1]]
    net <- results[[2]]
    
    
    ### 3:  Initialize shock vector
    shock.vector <- initialize.shockVector(r.shock, n.sources, net)
    
    ### 4: Apply shocks and simulate the contagion
    results <- simulate.contagion(net, BS, shock.vector)
    
    ### 5:  Report results
    Input.Parameters <-
      data.frame(H0, N, K, n.sources, r.shock, phi, eps)
    Balance.Sheets <- BS
    Contagion.State <- results[[1]]
    Contagion.Process <- results[[2]]
    
    Centralization <- net %n% "Centralization"
    Connectivity <- net %n% "Connectivity"
    isRegular <- regular
    Network.Properties <-
      data.frame(N,
                 K,
                 Connectivity,
                 Centralization,
                 isRegular,
                 pPeri,
                 pCore,
                 pToCore)
    
    create.Table(Input.Parameters, type = latex, cap = "Inputs for the Simulation.")
    create.Table(Balance.Sheets, type = latex, cap = "Balance Sheets of Banks.")
    create.Table(Contagion.State, type = latex, cap = "Final State of Contagion.")
    create.Table(Contagion.Process, type = latex, cap = "Growth of Contagion.")
    create.Table(Network.Properties, type = latex, cap = "Features of the Network.")
    
    if (DATA_MODE) {
      return(TRUE)
    }
    
    Coords <-
      plot.Network(net, f.name = "Ex-Ante", shock.vect = shock.vector)
    plot.Network(net, shock.vect = shock.vector, v.coords = Coords)
    
    dominos <- results[[3]]
    shocks <- results[[4]]
    for (r in 1:dim(dominos)[1]) {
      absorbtion <- dominos[r,]
      shock <- shocks[r,]
      fname <- paste("shock-", r - 1, sep = "")
      plot.Network(
        net,
        f.name = fname,
        insolvency = absorbtion,
        shock.vect = shock,
        v.coords = Coords
      )
      plot.Network(
        net,
        insolvency = absorbtion,
        shock.vect = shock,
        v.coords = Coords
      )
    }
    
    absorbtion <- results[[1]]$r.Absorb.E
    plot.Network(
      net,
      f.name = "Ex-Post",
      insolvency = absorbtion,
      shock.vect = shock,
      v.coords = Coords
    )
    
    # Degree Distribution:
    Degree.Distribution <- net %v% "deg"
    f <- hist(Degree.Distribution)
    dev.off()
    plot.DegreeDistro(
      f$counts,
      "Degree Centrality",
      "Frequency",
      title = "Degree Distribution",
      f.name = "DegreeDistro",
      isgrid = T,
      ltype = "h",
      ispoint = T
    )
    
    return(TRUE)
  }

run.experiment <- function(N,
                           r.shock,
                           n.sources,
                           K,
                           regular = TRUE,
                           phi = 0.4,
                           eps = 0.1,
                           H0 = 100,
                           pPeri = 0,
                           pCore = 1,
                           pToCore = 1,
                           latex = NULL,
                           n.runs = 10) {
  ### 1: Configure and create the network structure
  
  if (regular) {
    net <- generate.RegularGraph(N, K)
  } else {
    net <- build.CorePeriphery(N, K, pPeri, pCore, pToCore)
  }
  
  ### 2:  Initialize balance sheets
  results <- initialize.balanceSheets(net, H0, phi, eps)
  BS <- results[[1]]
  net <- results[[2]]
  
  defaults <- vector(mode = "integer")
  for (i in 1:n.runs) {
    ### 3:  Initialize shock vector
    shock.vector <-
      initialize.shockVector(r.shock, n.sources, net)
    
    ### 4: Apply shocks and simulate the contagion
    out0 <- simulate.contagion(net, BS, shock.vector)
    out1 <- out0[[2]]
    out2 <- as.vector(out1[, 4])
    defaults <- c(defaults, out2)
  }
  
  ### 5:  Report results
  Outputs.m <- matrix(defaults, ncol = n.runs)
  Outputs <- as.data.frame(Outputs.m)
  create.Table(Outputs)
  
  return(TRUE)
}



centralize <-
  function(N = 20,
           r = 3,
           prop = 1.0,
           keep = TRUE,
           symm = 1) {
    net.r <- generate.RegularGraph(N, r)
    net <- net.r
    shifted <- 0
    total <- sum(net[, ])
    
    for (j in r:1) {
      for (i in N:1) {
        rate <- shifted / total
        if (rate >= prop) {
          break
        }
        to <- (i + j - 1) %% N
        to <- to + 1
        net[i, to] <- 0
        if (keep) {
          net[i, j] <- 1
        }
        else {
          cent <- i %% r
          cent <- cent + 1
          net[i, cent] <- 1
        }
        
        shifted <- shifted + 1
      }
    }
    
    
    for (i in 1:N) {
      net[i, i] <- 0
    }
    
    if (prop > 0) {
      net[, ] <- symmetrize(net)
    }
    
    if (symm > 0) {
      net.r[, ] <- symmetrize(net.r)
      net[, ] <- symmetrize(net)
    }
    
    
    Degree.Distribution <- degree(net, gmode = "digraph")
    Centralization <- centralization(net, betweenness)
    Connectivity <- gden(net)
    
    net %v% "deg" <- Degree.Distribution
    net %n% "Centralization" <- Centralization
    net %n% "Connectivity" <- Connectivity
    
    Degree.Distribution <- degree(net.r, gmode = "digraph")
    Centralization <- centralization(net.r, betweenness)
    Connectivity <- gden(net.r)
    
    net.r %v% "deg" <- Degree.Distribution
    net.r %n% "Centralization" <- Centralization
    net.r %n% "Connectivity" <- Connectivity
    
    return(list("original" = net.r, "centralized" = net))
  }



centralize.complete.net <- function(N = 20,
                                    prop = 20,
                                    symm = 1) {
  net.r <- generate.RegularGraph(N, N)
  net <- net.r
  shifted <- 0
  
  for (i in N:2) {
    if (shifted >= prop) {
      break
    }
    net[i, ] <- 0
    net[, i] <- 0
    net[i, 1] <- 1
    shifted <- shifted + 1
  }
  
  for (i in 1:N) {
    net[i, i] <- 0
  }
  
  if (prop > 0) {
    net[, ] <- symmetrize(net)
  }
  
  if (symm > 0) {
    net.r[, ] <- symmetrize(net.r)
    net[, ] <- symmetrize(net)
  }
  
  
  Degree.Distribution <- degree(net, gmode = "digraph")
  Centralization <- centralization(net, betweenness)
  Connectivity <- gden(net)
  
  net %v% "deg" <- Degree.Distribution
  net %n% "Centralization" <- Centralization
  net %n% "Connectivity" <- Connectivity
  
  Degree.Distribution <- degree(net.r, gmode = "digraph")
  Centralization <- centralization(net.r, betweenness)
  Connectivity <- gden(net.r)
  
  net.r %v% "deg" <- Degree.Distribution
  net.r %n% "Centralization" <- Centralization
  net.r %n% "Connectivity" <- Connectivity
  
  return(list("original" = net.r, "centralized" = net))
}

run.contagion.2 <-
  function(N,
           r.shock,
           n.sources,
           K,
           cent.level = 1.0,
           fixed = F,
           regular = TRUE,
           phi = 0.4,
           eps = 0.1,
           H0 = 100,
           pPeri = 0,
           pCore = 1,
           pToCore = 1,
           latex = NULL) {
    ### 1: Configure and create the network structure
    
    res <- centralize(N, K, cent.level)
    net <- res$centralized
    
    ### 2:  Initialize balancesheets
    results <- initialize.balanceSheets(net, H0, phi, eps, fixed)
    BS <- results[[1]]
    net <- results[[2]]
    
    
    ### 3:  Initialize shock vector
    shock.vector <- initialize.shockVector(r.shock, n.sources, net)
    
    ### 4: Apply shocks and simulate the contagion
    results <- simulate.contagion(net, BS, shock.vector)
    
    ### 5:  Report results
    Input.Parameters <-
      data.frame(H0, N, K, n.sources, r.shock, phi, eps)
    Balance.Sheets <- BS
    Contagion.State <- results[[1]]
    Contagion.Process <- results[[2]]
    
    Centralization <- net %n% "Centralization"
    Connectivity <- net %n% "Connectivity"
    isRegular <- regular
    Network.Properties <-
      data.frame(N,
                 K,
                 Connectivity,
                 Centralization,
                 isRegular,
                 pPeri,
                 pCore,
                 pToCore)
    
    create.Table(Input.Parameters, type = latex, cap = "Inputs for the Simulation.")
    create.Table(Balance.Sheets, type = latex, cap = "Balance Sheets of Banks.")
    create.Table(Contagion.State, type = latex, cap = "Final State of Contagion.")
    create.Table(Contagion.Process, type = latex, cap = "Growth of Contagion.")
    create.Table(Network.Properties, type = latex, cap = "Features of the Network.")
    
    if (DATA_MODE) {
      return(TRUE)
    }
    
    Coords <-
      plot.Network(net, f.name = "Ex-Ante", shock.vect = shock.vector)
    plot.Network(net, shock.vect = shock.vector, v.coords = Coords)
    
    dominos <- results[[3]]
    shocks <- results[[4]]
    for (r in 1:dim(dominos)[1]) {
      absorbtion <- dominos[r,]
      shock <- shocks[r,]
      fname <- paste("shock-", r - 1, sep = "")
      plot.Network(
        net,
        f.name = fname,
        insolvency = absorbtion,
        shock.vect = shock,
        v.coords = Coords
      )
      plot.Network(
        net,
        insolvency = absorbtion,
        shock.vect = shock,
        v.coords = Coords
      )
    }
    
    absorbtion <- results[[1]]$r.Absorb.E
    plot.Network(
      net,
      f.name = "Ex-Post",
      insolvency = absorbtion,
      shock.vect = shock,
      v.coords = Coords
    )
    
    # Degree Distribution:
    Degree.Distribution <- net %v% "deg"
    f <- hist(Degree.Distribution)
    dev.off()
    plot.DegreeDistro(
      f$counts,
      "Degree Centrality",
      "Frequency",
      title = "Degree Distribution",
      f.name = "DegreeDistro",
      isgrid = T,
      ltype = "h",
      ispoint = T
    )
    
    return(TRUE)
  }

getContagionResults <- function(Contagion.Process, N) {
  CP <- Contagion.Process
  nSecondary <- max(CP$Secondary.Defaults)
  nDefaults <- max(CP$Default.Size)
  nPrimary <- nDefaults - nSecondary
  
  for (k in 1:length(CP$Shock.Nodes)) {
    if (CP$Default.Size[k] == N) {
      Tfin <- CP$Cumulative.Shock[k]
      break
    }
  }
  
  if (SCAN) {
    Tk <- rep(0, N)
    ind_pre = 0.1
    for (k in 1:length(CP$Shock.Nodes)) {
      ind <- CP$Secondary.Defaults[k]
      if (ind > ind_pre) {
        Tk[ind] <- CP$Cumulative.Shock[k]
        ind_pre <- ind
      }
    }
    return(list(nPrimary, nSecondary, Tk, Tfin))
    
    # An alternative way of display:
    val = 0
    for (k in N:1) {
      if (Tk[k]) {
        val = Tk[k]
      }
      else{
        Tk[k] = val
      }
    }
    return(list(nPrimary, nSecondary, Tk, Tfin))
  }
  
  T1 = 0
  T2 = 0
  dT = 0
  
  for (k in 1:length(CP$Shock.Nodes)) {
    if (CP$Secondary.Defaults[k] > 0) {
      T1 <- CP$Cumulative.Shock[k]
      break
    }
  }
  
  T2 <- Tfin
  dT <- T2 - T1
  
  return(list(nPrimary, nSecondary, T1, T2, dT))
}



exp.cent.with.fixed.assets <-
  function(N,
           K,
           A0,
           E0,
           d.min,
           cent.level,
           n.runs,
           unit.shock = 1.0) {
    r.shock = 1.0
    n.sources = N
    
    ### 1: Configure and create the network structure
    
    res <- centralize(N, K, cent.level)
    net <- res$centralized
    
    
    ### 2:  Initialize balancesheets
    if (VERSION == 2) {
      results <- initialize.balanceSheets.planA.v2(net, K, A0, E0, d.min)
    } else if (VERSION == 1) {
      results <-
        initialize.balanceSheets.planA.v1(net, K, A0, E0, d.min, fact.cent = SPEED.EXPOSURE)
    } else {
      results <-
        initialize.balanceSheets.planA.v1(net, K, A0, E0, d.min, fact.cent = 0)
    }
    
    BS <- results[[1]]
    net <- results[[2]]
    
    
    defaults <- vector(mode = "integer")
    
    for (i in 1:n.runs) {
      ### 3:  Initialize shock vector
      shock.vector <-
        initialize.shockVector(r.shock, n.sources, net, BS, unit.shock)
      
      ### 4: Apply shocks and simulate the contagion
      out0 <- simulate.contagion.fixed.shock(net, BS, shock.vector)
      out1 <- out0[[2]]
      out2 <- as.vector(out1[, 4])
      defaults <- c(defaults, out2)
    }
    
    ### 5:  Report results
    Outputs.m <- matrix(defaults, ncol = n.runs)
    Outputs <- as.data.frame(Outputs.m)
    create.Table(Outputs)
    
    create.Table(BS, type = latex, cap = "The Balance Sheets of the Banks")
    
    Centralization <- net %n% "Centralization"
    Connectivity <- net %n% "Connectivity"
    Network.Properties <-
      data.frame(N, K, Connectivity, Centralization)
    
    create.Table(Network.Properties, type = latex, cap = "The Features of the Network")
    
    return(TRUE)
  }

run.contagion.with.fixedShockLevels <-
  function(N,
           K,
           cent.level = 1.0,
           eps = 0.1,
           phi = 0.4) {
    r.shock = 1.0
    n.sources = N
    A0 = 100
    
    ### 1: Configure and create the network structure
    
    res <- centralize(N, K, cent.level)
    net <- res$centralized
    
    ### 2:  Initialize balancesheets
    results <- initialize.balanceSheets.planB(net, A0, phi, eps)
    
    BS <- results[[1]]
    net <- results[[2]]
    
    
    ### 3:  Initialize shock vector
    #A0 = (1 + eps * phi) / (1 - eps) * H0
    
    #A0 <- sum(BS$A)
    
    #fixed.shock.level = min(BS$A)
    
    fixed.shock.level = min(BS$A) / 2
    shocks <-
      initialize.shockVector(r.shock, n.sources, net, BS, level = fixed.shock.level)
    
    ### 4: Apply shocks and simulate the contagion
    results <- simulate.contagion.fixed.shock(net, BS, shocks)
    
    ### 5:  Report results
    Input.Parameters <- data.frame(N, K, A0, phi, eps)
    Balance.Sheets <- BS
    
    Contagion.State <- results[[1]]
    Internal.Shock <- Contagion.State$r.Shock.Int
    colnames(Contagion.State)[4] <- "Internal.Shock"
    Contagion.State$Internal.Shock <-  Internal.Shock
    
    Contagion.Process <- results[[2]]
    nShocks <- length(Contagion.Process$Shock.Nodes)
    
    
    Centralization <- net %n% "Centralization"
    Connectivity <- net %n% "Connectivity"
    
    
    if (SCAN) {
      outputs <- getContagionResults(Contagion.Process, N)
    }
    else {
      outputs <- getContagionResults(Contagion.Process, N)
    }
    nPrimary <- outputs[[1]]
    nSecondary <- outputs[[2]]
    if (SCAN) {
      Tk <- outputs[[3]]
      Tfin <- outputs[[4]]
      Summary.Results <-
        data.frame(N, Connectivity, Centralization, nPrimary, nSecondary, Tfin)
      for (i in 1:N) {
        Summary.Results$x <- Tk[i]
        col.name <- paste("T", i, sep = "")
        colnames(Summary.Results)[i + 6] <- col.name
      }
    }
    else {
      T1 <- outputs[[3]]
      T2 <- outputs[[4]]
      dT <- outputs[[5]]
      Summary.Results <-
        data.frame(N,
                   Connectivity,
                   Centralization,
                   nPrimary,
                   nSecondary,
                   T1,
                   T2 ,
                   dT)
    }
    
    if (DATA_MODE == 2) {
      create.Table(Summary.Results, type = latex, cap = "The Summary of the Results.")
      return(TRUE)
    }
    
    
    create.Table(Summary.Results, type = latex, cap = "The Summary of the Results.")
    create.Table(Input.Parameters, type = latex, cap = "Inputs for the Simulation.")
    create.Table(Balance.Sheets, type = latex, cap = "The Balance Sheets of the Banks.")
    create.Table(Contagion.State, type = latex, cap = "The Final State of the Contagion.")
    create.Table(Contagion.Process, type = latex, cap = "The Growth of Contagion.")
    
    
    Absorbtion.Levels <- as.data.frame(results[[3]])
    Shock.Sequences <- as.data.frame(results[[4]])
    Secondaries <- as.data.frame(results[[5]])
    
    
    create.Table(Absorbtion.Levels, type = latex, cap = "The Evolution of Absorption by Shareholders.")
    create.Table(Shock.Sequences, type = latex, cap = "The Accumulation of External Shock.")
    create.Table(Secondaries, type = latex, cap = "The Occurrances of Secondary Defaults.")
    
    if (PLOT_MODE == 0) {
      return(TRUE)
    }
    
    Coords <- plot.Network.v2(net, f.name = "Ex-Ante")
    plot.Network.v2(net, v.coords = Coords)
    
    absorbtion <- results[[1]]$r.Absorb.E
    plot.Network.v2(
      net,
      f.name = "Ex-Post",
      insolvency = absorbtion,
      isSecondary = Contagion.State$is.Secondary,
      v.coords = Coords
    )
    plot.Network.v2(
      net,
      insolvency = absorbtion,
      isSecondary = Contagion.State$is.Secondary,
      v.coords = Coords
    )
    dev.new()
    
    # Degree Distribution:
    Degree.Distribution <- net %v% "deg"
    f <- hist(Degree.Distribution)
    dev.off()
    plot.DegreeDistro(
      f$counts,
      "Degree Centrality",
      "Frequency",
      title = "Degree Distribution",
      f.name = "DegreeDistro",
      isgrid = T,
      ltype = "h",
      ispoint = T
    )
    
    
    if (PLOT_MODE < 2) {
      return(TRUE)
    }
    
    
    dominos <- results[[3]]
    shocks <- results[[4]]
    for (r in 1:dim(dominos)[1]) {
      absorbtion <- dominos[r,]
      shock <- shocks[r,]
      secs <- Secondaries[r,]
      fname <- paste("shock-", r - 1, sep = "")
      plot.Network.v2(
        net,
        f.name = fname,
        insolvency = absorbtion,
        shock.vect = shock,
        isSecondary = secs,
        v.coords = Coords
      )
      plot.Network.v2(
        net,
        insolvency = absorbtion,
        shock.vect = shock,
        isSecondary = secs,
        v.coords = Coords
      )
    }
    
    return(TRUE)
  }






run.test.experiment.planB <-
  function(N,
           K,
           eps,
           phi,
           nRuns = 100,
           rCentMin = 0,
           rCentMax = N,
           size.shock = 1.0) {
    n.sources = N
    A0 = 100
    
    defaults <- vector(mode = "numeric")
    
    #permutations <- generate_all_permutations(1:N)
    
    for (c in rCentMin:rCentMax) {
      cent.level = c / N
      
      ### 1: Configure and create the network structure
      
      res <- centralize(N, K, cent.level)
      net <- res$centralized
      Centralization <- net %n% "Centralization"
      Connectivity <- net %n% "Connectivity"
      
      ### 2:  Initialize balancesheets
      results <- initialize.balanceSheets.planB(net, A0, phi, eps)
      BS <- results[[1]]
      net <- results[[2]]
      
      if (MULTIPLE_FILES) {
        defaults <- vector(mode = "numeric")
      }
      
      #sequence <- c(4,20,5,19,6,18,7,17,8,16,9,15,10,14,11,13,12,3,2,1)
      #sequence <- c(2:20,1)
      #sequence <- c(20:2,1)
      #sequence <- c(1,2:20)
      sequence <- c(1, 20:2)
      
      reps <- sapply(sequence, function(x) {
        BS$A[x] %/% size.shock
      })
      sequence <- rep(sequence, times = reps)
      rates <- sapply(sequence, function(x) {
        size.shock / BS$A[x]
      })
      shocks <- list(sequence, rates)
      
      for (run in 1:nRuns) {
        ### 3:  Initialize shock vector
        
        #A0 = (1 + eps * phi) / (1 - eps) * H0
        
        #A0 <- sum(BS$A)
        #fixed.shock.level = A0 / N
        
        #fixed.shock.level = min(BS$A) / 2
        #shocks <- initialize.shockVector(r.shock, n.sources, net, BS, level = fixed.shock.level)
        
        
        #sequence <- permutations[[run]]
        #sequence <- rep(1:N, A0/(N * r.shock))
        
        
        
        ### 4: Apply shocks and simulate the contagion
        
        results <- simulate.contagion.fixed.shock(net, BS, shocks)
        
        ### 5:  Report results
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
      if (MULTIPLE_FILES) {
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
                     type = "pdf",
                     cap = cap.text)
      }
    }
    
    if (!MULTIPLE_FILES) {
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
      fname = "Results"
      cap.text = paste("Results.")
      create.Table(Out,
                   name = fname,
                   type = "pdf",
                   cap = cap.text)
    }
  }
