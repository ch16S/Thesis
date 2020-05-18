calculate_centralities<-function (x, except = NULL, include = NULL, weights = NULL){
  #y <- as_edgelist(x)
  y <- asNetwork(x)
  if (is_directed(x) && is_weighted(x)) {
    centrality_funcs <- list(`Alpha Centrality` = function(x) alpha.centrality(x, 
                                                                               weights = NULL), `Burt's Constraint` = function(x) constraint(x, 
                                                                                                                                             weights = NULL), `Page Rank` = function(x) page_rank(x)$vector, 
                             `Average Distance` = function(x) averagedis(x, weights = NULL), 
                             `Barycenter Centrality` = function(x) barycenter(x, 
                                                                              weights = NULL), `BottleNeck Centrality` = function(x) bottleneck(x), 
                             `Centroid value` = function(x) centroid(x, weights = NULL), 
                             `Closeness Centrality (Freeman)` = function(x) closeness.freeman(x, 
                                                                                              weights = NULL), ClusterRank = function(x) clusterrank(x), 
                             `Decay Centrality` = function(x) decay(x, weights = NULL), 
                             `Degree Centrality` = function(x) centr_degree(x)$res, 
                             `Diffusion Degree` = function(x) diffusion.degree(x), 
                             `DMNC - Density of Maximum Neighborhood Component` = function(x) dmnc(x), 
                             `Eccentricity Centrality` = function(x) eccentricity(x), 
                             `eigenvector centralities` = function(x) eigen_centrality(x, 
                                                                                       weights = NULL)$vector, `K-core Decomposition` = function(x) coreness(x), 
                             `Geodesic K-Path Centrality` = function(x) geokpath(x, 
                                                                                 weights = NULL), `Katz Centrality (Katz Status Index)` = function(x) katzcent(x), 
                             `Kleinberg's authority centrality scores` = function(x) authority_score(x, 
                                                                                                     weights = NULL)$vector, `Kleinberg's hub centrality scores` = function(x) hub_score(x, 
                                                                                                                                                                                         weights = NULL)$vector, `clustering coefficient` = function(x) transitivity(x, 
                                                                                                                                                                                                                                                                     weights = NULL, type = "local"), `Lin Centrality` = function(x) lincent(x, 
                                                                                                                                                                                                                                                                                                                                             weights = NULL), `Lobby Index (Centrality)` = function(x) lobby(x), 
                             `Markov Centrality` = function(x) markovcent(x), 
                             `Radiality Centrality` = function(x) radiality(x, 
                                                                            weights = NULL), `Shortest-Paths Betweenness Centrality` = function(x) betweenness(x), 
                             `Current-Flow Closeness Centrality` = function(x) closeness.currentflow(x, 
                                                                                                     weights = NULL), `Closeness centrality (Latora)` = function(x) closeness.latora(x, 
                                                                                                                                                                                     weights = NULL), `Communicability Betweenness Centrality` = function(x) communibet(x), 
                             `Community Centrality` = function(x) communitycent(x), 
                             `Cross-Clique Connectivity` = function(x) crossclique(x), 
                             `Entropy Centrality` = function(x) entropy(x, weights = NULL), 
                             `EPC - Edge Percolated Component` = function(x) epc(x), 
                             `Laplacian Centrality` = function(x) laplacian(x), 
                             `Leverage Centrality` = function(x) leverage(x), 
                             `MNC - Maximum Neighborhood Component` = function(x) mnc(x), 
                             `Hubbell Index` = function(x) hubbell(x, weights = NULL), 
                             `Semi Local Centrality` = function(x) semilocal(x), 
                             `Closeness Vitality` = function(x) closeness.vitality(x, 
                                                                                   weights = NULL), `Residual Closeness Centrality` = function(x) closeness.residual(x, 
                                                                                                                                                                     weights = NULL), `Stress Centrality` = function(x) stresscent(y), 
                             `Load Centrality` = function(x) loadcent(y), `Flow Betweenness Centrality` = function(x) flowbet(y), 
                             `Information Centrality` = function(x) infocent(y), 
                             `Weighted Vertex Degree` = function(x) strength(x, 
                                                                             vids = V(x), mode = "all", weights = NULL), `Harary Centrality` = function(x) graphcent(y, 
                                                                                                                                                                     gmode = "graph", diag = T, cmode = "directed"), 
                             `Dangalchev Closeness Centrality` = function(x) dangalchev_closeness_centrality(x, 
                                                                                                             vids = V(x), mode = "all", weights = NULL), `Group Centrality` = function(x) group_centrality(x, 
                                                                                                                                                                                                           vids = V(x)), `Harmonic Centrality` = function(x) harmonic_centrality(x, 
                                                                                                                                                                                                                                                                                 vids = V(x), mode = "all", weights = NULL), `Local Bridging Centrality` = function(x) local_bridging_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                 vids = V(x)), `Wiener Index Centrality` = function(x) wiener_index_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               vids = V(x), mode = "all", weights = NULL))
    if (!is.null(include)) {
      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), 
                                                     include)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
    else {
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), 
                                                   except)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
  }
  if (!is_directed(x) && is_weighted(x)) {
    centrality_funcs <- list(`subgraph centrality scores` = function(x) subgraph.centrality(x), 
                             `Topological Coefficient` = function(x) topocoefficient(x), 
                             `Average Distance` = function(x) averagedis(x, weights = NULL), 
                             `Barycenter Centrality` = function(x) barycenter(x, 
                                                                              weights = NULL), `BottleNeck Centrality` = function(x) bottleneck(x), 
                             `Centroid value` = function(x) centroid(x, weights = NULL), 
                             `Closeness Centrality (Freeman)` = function(x) closeness.freeman(x, 
                                                                                              weights = NULL), ClusterRank = function(x) clusterrank(x), 
                             `Decay Centrality` = function(x) decay(x, weights = NULL), 
                             `Degree Centrality` = function(x) centr_degree(x)$res, 
                             `Diffusion Degree` = function(x) diffusion.degree(x), 
                             `DMNC - Density of Maximum Neighborhood Component` = function(x) dmnc(x), 
                             `Eccentricity Centrality` = function(x) eccentricity(x), 
                             `eigenvector centralities` = function(x) eigen_centrality(x, 
                                                                                       weights = NULL)$vector, `K-core Decomposition` = function(x) coreness(x), 
                             `Geodesic K-Path Centrality` = function(x) geokpath(x, 
                                                                                 weights = NULL), `Katz Centrality (Katz Status Index)` = function(x) katzcent(x), 
                             `Kleinberg's authority centrality scores` = function(x) authority_score(x, 
                                                                                                     weights = NULL)$vector, `Kleinberg's hub centrality scores` = function(x) hub_score(x, 
                                                                                                                                                                                         weights = NULL)$vector, `clustering coefficient` = function(x) transitivity(x, 
                                                                                                                                                                                                                                                                     weights = NULL, type = "local"), `Lin Centrality` = function(x) lincent(x, 
                                                                                                                                                                                                                                                                                                                                             weights = NULL), `Lobby Index (Centrality)` = function(x) lobby(x), 
                             `Markov Centrality` = function(x) markovcent(x), 
                             `Radiality Centrality` = function(x) radiality(x, 
                                                                            weights = NULL), `Shortest-Paths Betweenness Centrality` = function(x) betweenness(x), 
                             `Current-Flow Closeness Centrality` = function(x) closeness.currentflow(x, 
                                                                                                     weights = NULL), `Closeness centrality (Latora)` = function(x) closeness.latora(x, 
                                                                                                                                                                                     weights = NULL), `Communicability Betweenness Centrality` = function(x) communibet(x), 
                             `Community Centrality` = function(x) communitycent(x), 
                             `Cross-Clique Connectivity` = function(x) crossclique(x), 
                             `Entropy Centrality` = function(x) entropy(x, weights = NULL), 
                             `EPC - Edge Percolated Component` = function(x) epc(x), 
                             `Laplacian Centrality` = function(x) laplacian(x), 
                             `Leverage Centrality` = function(x) leverage(x), 
                             `MNC - Maximum Neighborhood Component` = function(x) mnc(x), 
                             `Hubbell Index` = function(x) hubbell(x, weights = NULL), 
                             `Semi Local Centrality` = function(x) semilocal(x), 
                             `Closeness Vitality` = function(x) closeness.vitality(x, 
                                                                                   weights = NULL), `Residual Closeness Centrality` = function(x) closeness.residual(x, 
                                                                                                                                                                     weights = NULL), `Stress Centrality` = function(x) stresscent(y), 
                             `Load Centrality` = function(x) loadcent(y), `Flow Betweenness Centrality` = function(x) flowbet(y), 
                             `Information Centrality` = function(x) infocent(y), 
                             `Weighted Vertex Degree` = function(x) strength(x, 
                                                                             vids = V(x), mode = "all", weights = NULL), `Harary Centrality` = function(x) graphcent(y, 
                                                                                                                                                                     gmode = "graph", diag = T, cmode = "undirected"), 
                             `Dangalchev Closeness Centrality` = function(x) dangalchev_closeness_centrality(x, 
                                                                                                             vids = V(x), mode = "all", weights = NULL), `Group Centrality` = function(x) group_centrality(x, 
                                                                                                                                                                                                           vids = V(x)), `Harmonic Centrality` = function(x) harmonic_centrality(x, 
                                                                                                                                                                                                                                                                                 vids = V(x), mode = "all", weights = NULL), `Local Bridging Centrality` = function(x) local_bridging_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                 vids = V(x)), `Wiener Index Centrality` = function(x) wiener_index_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                               vids = V(x), mode = "all", weights = NULL))
    if (!is.null(include)) {
      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), 
                                                     include)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
    else {
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), 
                                                   except)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
  }
  if (!is_directed(x) && !is_weighted(x)) {
    centrality_funcs <- list(`subgraph centrality scores` = function(x) subgraph.centrality(x), 
                             `Topological Coefficient` = function(x) topocoefficient(x), 
                             `Average Distance` = function(x) averagedis(x, weights = NULL), 
                             `Barycenter Centrality` = function(x) barycenter(x, 
                                                                              weights = NULL), `BottleNeck Centrality` = function(x) bottleneck(x), 
                             `Centroid value` = function(x) centroid(x, weights = NULL), 
                             `Closeness Centrality (Freeman)` = function(x) closeness.freeman(x, 
                                                                                              weights = NULL), ClusterRank = function(x) clusterrank(x), 
                             `Decay Centrality` = function(x) decay(x, weights = NULL), 
                             `Degree Centrality` = function(x) centr_degree(x)$res, 
                             `Diffusion Degree` = function(x) diffusion.degree(x), 
                             `DMNC - Density of Maximum Neighborhood Component` = function(x) dmnc(x), 
                             `Eccentricity Centrality` = function(x) eccentricity(x), 
                             `eigenvector centralities` = function(x) eigen_centrality(x, 
                                                                                       weights = NULL)$vector, `K-core Decomposition` = function(x) coreness(x), 
                             `Geodesic K-Path Centrality` = function(x) geokpath(x, 
                                                                                 weights = NULL), `Katz Centrality (Katz Status Index)` = function(x) katzcent(x), 
                             `Kleinberg's authority centrality scores` = function(x) authority_score(x, 
                                                                                                     weights = NULL)$vector, `Kleinberg's hub centrality scores` = function(x) hub_score(x, 
                                                                                                                                                                                         weights = NULL)$vector, `clustering coefficient` = function(x) transitivity(x, 
                                                                                                                                                                                                                                                                     weights = NULL, type = "local"), `Lin Centrality` = function(x) lincent(x, 
                                                                                                                                                                                                                                                                                                                                             weights = NULL), `Lobby Index (Centrality)` = function(x) lobby(x), 
                             `Markov Centrality` = function(x) markovcent(x), 
                             `Radiality Centrality` = function(x) radiality(x, 
                                                                            weights = NULL), `Shortest-Paths Betweenness Centrality` = function(x) betweenness(x), 
                             `Current-Flow Closeness Centrality` = function(x) closeness.currentflow(x, 
                                                                                                     weights = NULL), `Closeness centrality (Latora)` = function(x) closeness.latora(x, 
                                                                                                                                                                                     weights = NULL), `Communicability Betweenness Centrality` = function(x) communibet(x), 
                             `Community Centrality` = function(x) communitycent(x), 
                             `Cross-Clique Connectivity` = function(x) crossclique(x), 
                             `Entropy Centrality` = function(x) entropy(x, weights = NULL), 
                             `EPC - Edge Percolated Component` = function(x) epc(x), 
                             `Laplacian Centrality` = function(x) laplacian(x), 
                             `Leverage Centrality` = function(x) leverage(x), 
                             `MNC - Maximum Neighborhood Component` = function(x) mnc(x), 
                             `Hubbell Index` = function(x) hubbell(x, weights = NULL), 
                             `Semi Local Centrality` = function(x) semilocal(x), 
                             `Closeness Vitality` = function(x) closeness.vitality(x, 
                                                                                   weights = NULL), `Residual Closeness Centrality` = function(x) closeness.residual(x, 
                                                                                                                                                                     weights = NULL), `Stress Centrality` = function(x) stresscent(y), 
                             `Load Centrality` = function(x) loadcent(y), `Flow Betweenness Centrality` = function(x) flowbet(y), 
                             `Information Centrality` = function(x) infocent(y), 
                             `Harary Centrality` = function(x) graphcent(y, gmode = "graph", 
                                                                         diag = T, cmode = "undirected"), `Dangalchev Closeness Centrality` = function(x) dangalchev_closeness_centrality(x, 
                                                                                                                                                                                        vids = V(x), mode = "all", weights = NULL), `Group Centrality` = function(x) group_centrality(x, 
                                                                                                                                                                                                                                                                                      vids = V(x)), `Harmonic Centrality` = function(x) harmonic_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                            vids = V(x), mode = "all", weights = NULL), `Local Bridging Centrality` = function(x) local_bridging_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            vids = V(x)), `Wiener Index Centrality` = function(x) wiener_index_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          vids = V(x), mode = "all", weights = NULL))
    if (!is.null(include)) {
      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), 
                                                     include)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
    else {
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), 
                                                   except)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
  }
  if (is_directed(x) && !is_weighted(x)) {
    centrality_funcs <- list(`Alpha Centrality` = function(x) alpha.centrality(x, 
                                                                               weights = NULL), `Bonacich power centralities of positions` = function(x) bonpow(x), 
                             `Page Rank` = function(x) page_rank(x)$vector, `Average Distance` = function(x) averagedis(x, 
                                                                                                                        weights = NULL), `Barycenter Centrality` = function(x) barycenter(x, 
                                                                                                                                                                                          weights = NULL), `BottleNeck Centrality` = function(x) bottleneck(x), 
                             `Centroid value` = function(x) centroid(x, weights = NULL), 
                             `Closeness Centrality (Freeman)` = function(x) closeness.freeman(x, 
                                                                                              weights = NULL), ClusterRank = function(x) clusterrank(x), 
                             `Decay Centrality` = function(x) decay(x, weights = NULL), 
                             `Degree Centrality` = function(x) centr_degree(x)$res, 
                             `Diffusion Degree` = function(x) diffusion.degree(x), 
                             `DMNC - Density of Maximum Neighborhood Component` = function(x) dmnc(x), 
                             `Eccentricity Centrality` = function(x) eccentricity(x), 
                             `eigenvector centralities` = function(x) eigen_centrality(x, 
                                                                                       weights = NULL)$vector, `K-core Decomposition` = function(x) coreness(x), 
                             `Geodesic K-Path Centrality` = function(x) geokpath(x, 
                                                                                 weights = NULL), `Katz Centrality (Katz Status Index)` = function(x) katzcent(x), 
                             `Kleinberg's authority centrality scores` = function(x) authority_score(x, 
                                                                                                     weights = NULL)$vector, `Kleinberg's hub centrality scores` = function(x) hub_score(x, 
                                                                                                                                                                                         weights = NULL)$vector, `clustering coefficient` = function(x) transitivity(x, 
                                                                                                                                                                                                                                                                     weights = NULL, type = "local"), `Lin Centrality` = function(x) lincent(x, 
                                                                                                                                                                                                                                                                                                                                             weights = NULL), `Lobby Index (Centrality)` = function(x) lobby(x), 
                             `Markov Centrality` = function(x) markovcent(x), 
                             `Radiality Centrality` = function(x) radiality(x, 
                                                                            weights = NULL), `Shortest-Paths Betweenness Centrality` = function(x) betweenness(x), 
                             `Current-Flow Closeness Centrality` = function(x) closeness.currentflow(x, 
                                                                                                     weights = NULL), `Closeness centrality (Latora)` = function(x) closeness.latora(x, 
                                                                                                                                                                                     weights = NULL), `Communicability Betweenness Centrality` = function(x) communibet(x), 
                             `Community Centrality` = function(x) communitycent(x), 
                             `Cross-Clique Connectivity` = function(x) crossclique(x), 
                             `Entropy Centrality` = function(x) entropy(x, weights = NULL), 
                             `EPC - Edge Percolated Component` = function(x) epc(x), 
                             `Laplacian Centrality` = function(x) laplacian(x), 
                             `Leverage Centrality` = function(x) leverage(x), 
                             `MNC - Maximum Neighborhood Component` = function(x) mnc(x), 
                             `Hubbell Index` = function(x) hubbell(x, weights = NULL), 
                             `Semi Local Centrality` = function(x) semilocal(x), 
                             `Closeness Vitality` = function(x) closeness.vitality(x, 
                                                                                   weights = NULL), `Residual Closeness Centrality` = function(x) closeness.residual(x, 
                                                                                                                                                                     weights = NULL), `Stress Centrality` = function(x) stresscent(y), 
                             `Load Centrality` = function(x) loadcent(y), `Flow Betweenness Centrality` = function(x) flowbet(y), 
                             `Information Centrality` = function(x) infocent(y), 
                             `Harary Centrality` = function(x) graphcent(y, gmode = "graph", 
                                                                         diag = T, cmode = "directed"), `Dangalchev Closeness Centrality` = function(x) dangalchev_closeness_centrality(x, 
                                                                                                                                                                                        vids = V(x), mode = "all", weights = NULL), `Group Centrality` = function(x) group_centrality(x, 
                                                                                                                                                                                                                                                                                      vids = V(x)), `Harmonic Centrality` = function(x) harmonic_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                            vids = V(x), mode = "all", weights = NULL), `Local Bridging Centrality` = function(x) local_bridging_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                            vids = V(x)), `Wiener Index Centrality` = function(x) wiener_index_centrality(x, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          vids = V(x), mode = "all", weights = NULL))
    if (!is.null(include)) {
      centrality_funcs <- centrality_funcs[intersect(names(centrality_funcs), 
                                                     include)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
    else {
      centrality_funcs <- centrality_funcs[setdiff(names(centrality_funcs), 
                                                   except)]
      n <- names(centrality_funcs)
      warningsText <- ""
      result <- lapply(setNames(n, n), function(functionName, 
                                                x) {
        f <- centrality_funcs[[functionName]]
        tryCatch(f(x), error = function(e) {
          warningsText <- paste0(warningsText, "\nError in ", 
                                 functionName, ":\n", e$message)
          return(NULL)
        })
      }, x)
      if (nchar(warningsText) > 0) 
        warning(warningsText)
      return(result)
    }
  }
}
