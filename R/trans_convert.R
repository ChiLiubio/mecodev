#' @title
#' Abundance table transformation.
#'
#' @description
#' Abundance table standardization, transformation or normalization for the otu_table in microtable object.
#' The input dataset must be a microtable object.
#'
#' @export
trans_convert <- R6Class(classname = "trans_convert",
	public = list(
		#' @description
		#' Get a transposed abundance table in the object. Rows are samples and columns are features. 
		#' This can make the further operations same with the traditional ecological methods.
		#' @param dataset the object of \code{\link{microtable}} Class.
		#' @return transposed abundance table stored in the object.
		#' @examples
		#' data(dataset)
		#' t1 <- trans_convert$new(dataset = dataset)
		initialize = function(dataset = NULL)
			{
			abund_table <- dataset$otu_table
			abund_table <- t(abund_table)
			self$abund_table <- abund_table
			self$dataset <- dataset
		},
		#' @description
		#' The function offers some methods based on the \code{\link{decostand}} function in vegan package and several famous methods in microbial field.
		#' MARGIN = 1 to represent samples, i.e. rows and MARGIN = 2 to represent features, i.e. columns.
		#' @param method default NULL; See the following details and available options. \cr 
		#' \cr 
		#' Methods from \code{\link{decostand}}:
		#' \itemize{
		#'   \item \code{total}: divide by margin total (default MARGIN = 1, i.e. rows - samples).
		#'   \item \code{max}: divide by margin maximum (default MARGIN = 2, i.e. columns - features).
		#'   \item \code{normalize}:  make margin sum of squares equal to one (default MARGIN = 1).
		#'   \item \code{range}: standardize values into range 0...1 (default MARGIN = 2). If all values are constant, they will be transformed to 0.
		#'   \item \code{standardize}: scale x to zero mean and unit variance (default MARGIN = 2).
		#'   \item \code{pa}: scale x to presence/absence scale (0/1).
		#'   \item \code{log}: logarithmic transformation as suggested by Anderson et al. (2006): log_b (x) + 1 for x > 0, where b is the base of the logarithm; zeros are left as zeros. Higher bases give less weight to quantities and more to presences, and logbase = Inf gives the presence/absence scaling. Please note this is not log(x+1). Anderson et al. (2006) suggested this for their (strongly) modified Gower distance (implemented as method = "altGower" in vegdist), but the standardization can be used independently of distance indices.
		#' }
		#' Methods useful in microbial field:
		#' \itemize{
		#'   \item \code{CLR}: Centered log-ratio normalization (default MARGIN = 2, i.e. features).
		#' }
		#' @param MARGIN default NULL; 1 = samples, and 2 = features of abundance table.
		#' @param logbase default exp(1); The logarithm base used in method = "log" or "CLR".
		#' @param ... parameters pass to \code{\link{decostand}}.
		#' 
		#' @return a new microtable object; rows are features.
		#' @examples
		#' newdataset <- t1$convert(method = "log")
		convert = function(method = NULL, MARGIN = NULL, logbase = exp(1), ...)
			{
			abund_table <- self$abund_table
			
			# use decostand function
			if(method %in% c("total", "max", "frequency", "normalize", "range", "rank", "standardize", "pa", "chi.square", "hellinger", "log")){
				if(is.null(MARGIN)){
					MARGIN <- switch(method, total = 1, max = 2, frequency = 2, normalize = 1, range = 2, rank = 1, standardize = 2, chi.square = 1, NULL)
				}
				res_table <- vegan::decostand(x = abund_table, method = method, MARGIN = MARGIN, logbase = logbase, ...)
			}
			# Centered log-ratio normalization
			if(method == "CLR"){
				if(!require(SpiecEasi)){
					stop("SpiecEasi package is not installed! See https://github.com/zdk123/SpiecEasi for installation.")
				}
				if(is.null(MARGIN)){
					MARGIN <- 2
				}
				res_table <- SpiecEasi::clr(abund_table, marg = MARGIN, base = logbase)
			}
			
			res_dataset <- clone(self$dataset)
			res_dataset$otu_table <- as.data.frame(t(res_table))
			res_dataset
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
