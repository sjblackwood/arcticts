#' Detect Rare Events Using DBSCAN
#'
#' This function applies the DBSCAN (Density-Based Spatial Clustering of Applications
#' with Noise) algorithm to identify rare events (anomalies) in a dataset using the
#' `dbscan` package. Points assigned to the noise cluster (label 0) are considered
#' rare events. The function computes cluster assignments and returns anomaly labels
#' and scores based on the DBSCAN model.
#'
#' @param x A numeric matrix or data frame with no missing values. Rows are
#'   observations, and columns are features.
#' @param eps Numeric, the maximum distance between two points for them to be
#'   considered in the same neighborhood (default: estimated using k-nearest neighbors).
#' @param minPts Integer, the minimum number of points required to form a dense
#'   region (default: 5).
#' @param scale Logical, whether to scale the data to have mean 0 and standard
#'   deviation 1 before clustering (default: TRUE).
#' @param columns Optional character vector, names of columns in `x` to use for
#'   clustering (default: NULL, uses all numeric columns).
#' @param ... Additional arguments passed to `dbscan::dbscan`.
#' @return A list containing:
#'   \itemize{
#'     \item \code{scores}: Numeric vector of approximate anomaly scores (inverse of
#'       k-nearest neighbor distances, normalized; higher values indicate more likely
#'       anomalies).
#'     \item \code{is_anomaly}: Logical vector indicating whether each observation is
#'       a rare event (TRUE for noise points, cluster label 0).
#'     \item \code{cluster}: Integer vector of cluster assignments (0 for noise, 1+ for clusters).
#'     \item \code{model}: The fitted DBSCAN model object.
#'   }
#' @details
#' DBSCAN identifies clusters based on density, labeling points that do not belong to
#' any dense cluster as noise (cluster 0), which are treated as rare events. If `eps`
#' is not specified, it is estimated as the 90th percentile of the k-nearest neighbor
#' distances (where k = `minPts`). Scaling is recommended to ensure features contribute
#' equally to distance calculations. The `scores` are approximate, based on inverse
#' k-nearest neighbor distances, and should be interpreted cautiously.
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' data <- matrix(rnorm(1000), nrow = 500)
#' data[1:5, ] <- data[1:5, ] + 10  # Add some outliers
#' result <- rare_dbscan(data)
#' table(result$is_anomaly)
#' plot(data, col = ifelse(result$is_anomaly, "red", "blue"), pch = 19)
#'
#' result <- rare_dbscan(iris[, 1:4], minPts = 20)
#' table(result$is_anomaly)
#'
#' result <- rare_dbscan(swiss)
#' table(result$is_anomaly)
#' }
#' @importFrom dbscan dbscan kNNdist
#' @importFrom stats scale
#' @export
#'
rare_dbscan <- function(x,
                        eps = NULL,
                        minPts = NULL,
                        scale = TRUE,
                        ...) {
    # Define minPts
    if (is.null(minPts)) {
        minPts <- 2 * ncol(x) + 1
    }

    # Input validation
    if (!is.matrix(x) && !is.data.frame(x)) {
        stop("`x` must be a matrix or data frame.")
    }

    # Check for non-finite values (NA, NaN, Inf)
    if (any(!is.finite(as.matrix(x)))) {
        stop("`x` contains NA, NaN, or Inf values. DBSCAN does not support missing values.")
    }

    # Validate other parameters
    if (!is.null(eps) && (!is.numeric(eps) || eps <= 0)) {
        stop("`eps` must be a positive numeric value or NULL.")
    }
    # if (!is.integer(minPts) && !is.numeric(minPts) || minPts < 1) {
    #     stop("`minPts` must be a positive integer.")
    # }
    minPts <- as.integer(minPts)
    if (!requireNamespace("dbscan", quietly = TRUE)) {
        stop("Package 'dbscan' is required. Please install it using install.packages('dbscan').")
    }

    # Scale data if requested
    if (scale) {
        x <- as.matrix(scale(x))
    } else {
        x <- as.matrix(x)
    }

    # Estimate eps if not provided
    if (is.null(eps)) {
        knn_dist <- dbscan::kNNdist(x, k = minPts)
        # eps <- quantile(knn_dist, 0.95)  # Use 95th percentile of kNN distances
        # message("`eps` not specified. Estimated as ", round(eps, 4), " based on kNN distances.")

        knee <- kneedle(1:length(knn_dist), sort(knn_dist))
        eps <- knee[2]
    }

    # Fit DBSCAN model
    model <- dbscan::dbscan(x, eps = eps, minPts = minPts, ...)

    # Extract cluster assignments (0 for noise, 1+ for clusters)
    cluster <- model$cluster

    # Identify anomalies (noise points, cluster 0)
    is_anomaly <- cluster == 0

    # Compute approximate anomaly scores (inverse kNN distances, normalized)
    knn_dist <- dbscan::kNNdist(x, k = minPts)
    scores <- 1 / (knn_dist + 1e-10)  # Avoid division by zero
    scores <- (scores - min(scores)) / (max(scores) - min(scores))  # Normalize to [0,1]

    # Return results
    return(list(
        scores = scores,
        is_anomaly = is_anomaly,
        cluster = cluster,
        model = model
    ))
}
