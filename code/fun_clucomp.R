if (FALSE) { # rows for debugging are not executed
    x = dm
    k = K
    methods = c("PAM", "average", "centroid", "complete", "abc", "1wew")
    criteria = c("avg.silwidth", "dunn", "dunn2", "abc", "1wew")
    TMP = clucomp(dm, k = 2:5,  methods = c("average", "centroid", "complete", "PAM"), criteria = criteria)
}

# function to cluster and compare different methods
clucomp <- function(x, # distance matrix
                    k, # integer vector of number of clusters to explore
                    methods = NULL,
                    criteria = NULL,
                    plot = TRUE) {
    # checking
    K <- k
#    K <- sort(unique(round(k[k <= nrow(x) & k > 0])))
#    if (length(K) < 1) {stop("need to specify potential number of clusters, k, in the range from 1 to nrow(x).")}

    M <- c("average", "centroid", "complete", "median", "PAM", "single", "ward.D", "ward.D2")
    if (is.null(methods)) {
        methods <- M
    }
    tmp <- base::setdiff(methods, M)
    if (length(tmp) > 0) {
        warning(paste("Some of the specified methods do not match the function's methods and will be disregarded:",
                      paste0(tmp, collapse = ", ")))
    }
    methods <- base::intersect(M, methods)
    if (length(methods) < 1) {
        stop(paste0("need to specify methods as one or more of: '", paste0(M, collapse = "', '"), "'"))
    }

    Crit <- c("avg.silwidth", "ch", "dunn", "dunn2", "entropy", "wb.ratio")
    if (is.null(criteria)) {
        criteria <- Crit
    }
    # -1 if larger is better; 1 if smaller is better
    CritDirection <- data.frame(Criterion = Crit, Direction = c(-1, -1, -1, -1, 1, 1))
    tmp <- base::setdiff(criteria, Crit)
    if (length(tmp) > 0) {
        warning(paste("Some of the specified criteria do not match the function's criteria and will be disregarded:",
                      paste0(tmp, collapse = ", ")))
    }
    criteria <- base::intersect(Crit, criteria)
    if (length(criteria) < 1) {
        stop(paste0("need to specify criteria as one or more of: '", paste0(Crit, collapse = "', '"), "'"))
    }

    # clustering
    CL <- lapply(methods, function(m) {
        if (identical(m, "PAM")) { # k-medoids
            lapply(K, function(k) {
                cl <- cluster::pam(x, k = k, diss = TRUE)$clustering
                fpc::cluster.stats(x, cl)
            })
        } else { # hierarchical
            cltree <- hclust(x, method = m)
            lapply(K, function(k) {
                cl <- cutree(cltree, k = k)
                fpc::cluster.stats(x, cl)
            })
        }
    })
    names(CL) <- methods

    # evaluating (extract performance criteria from the results)
    EV <- lapply(methods, function(m) { # outer loop is over methods; m = "average"
        tmp <- sapply(CL[[m]], function(clk) unlist(clk[criteria])) # inner loop is over K; clk = CL[[m]][[1]]
        data.frame(t(tmp), k = K, Method = m)
    })
    EV <- do.call(rbind, EV)
    EVlong <- EV %>%
        tidyr::pivot_longer(cols = 1:length(criteria),
                            names_to = "Criterion", values_to = "Value")

    # ranking (find ranks such that smaller rank (#1) is better)
    EVlongr <- EVlong %>%
        left_join(CritDirection, by = "Criterion") %>%
        mutate(Value = Value * Direction)
    for (cr in criteria) { # cr = "dunn"
        EVlongr$Value[EVlongr$Criterion == cr] <- base::rank(EVlongr$Value[EVlongr$Criterion == cr])
    }
    # calculate average ranks and find a combination with the smallest average
    rank_mk <- with(EVlongr, tapply(Value, list(Method, k), mean))
    # rank_mk
    tmp <- which(rank_mk == min(rank_mk), arr.ind = TRUE)
    opt_method <- rownames(tmp)[1]
    opt_k = tmp[1, 2] + min(K) - 1

    # robustness
    # best ranked method overall
    rank_m <- with(EVlongr, tapply(Value, list(Method), mean))
    # best ranked k overall
    rank_k <- with(EVlongr, tapply(Value, list(k), mean))

    # plotting
    clplot <- EVlong %>%
        ggplot(aes(x = k, y = Value, group = Method, color = Method)) +
        geom_point(aes(shape = Method)) + geom_line() +
        facet_wrap(vars(Criterion), scales = "free_y") +
        scale_color_viridis_d(end = 0.85, direction = -1) +
        theme_minimal()
    if (plot) {
        print(clplot)
    }

    # output
    return(list(EV = EV,
                clplot = clplot,
                opt_method = opt_method, opt_k = opt_k,
                rank_mk = rank_mk, rank_m = rank_m, rank_k = rank_k
    ))
    # Note: ranks are from 1 to length(K)*length(methods); smaller the better.
}
