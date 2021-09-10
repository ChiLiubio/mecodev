#' @title
#' Biomass evaluation and biological interaction prediction for Time Series Data analysis.
#'
#' @description
#' Biomass evaluation and biological interaction prediction for Time Series Data analysis based on the method from BEEM package, Li et al. 2019. <doi: 10.1186/s40168-019-0729-z>.
#'
#' @export
trans_ts <- R6Class(classname = "trans_ts",
	public = list(
		#' @description
		#' Prepare data for the following analysis.
		#' 
		#' @param dataset the object of \code{\link{microtable}} class.
		#'   Two columns with exact names in sample_table are necessary; one is 'Time', which is the time point and should be the numeric class;
		#'   the other is 'Rep', which represents the biological replicates and is also numeric class. If no replicates, use 1 to represent 1 replicate.
		#' @param filter_thres default 0.001; the relative abundance threshold of taxa. 
		#' @return abund_table and sample_table stored in the object.
		#' @examples
		#' data(gut_microb_ts)
		#' t1 <- trans_ts$new(dataset = gut_microb_ts, filter_thres = 0)
		initialize = function(dataset = NULL, filter_thres = 0.001)
			{
			if(is.null(dataset)){
				stop("No dataset provided !")
			}
			dataset1 <- clone(dataset)
			dataset1$tidy_dataset()
			
			sample_table <- dataset1$sample_table
			# check sample_table
			if(! any(colnames(sample_table) %in% "Time")){
				stop("No 'Time' column found in the dataset$sample_table column names!")
			}
			if(! any(colnames(sample_table) %in% "Rep")){
				stop("No 'Rep' column found in the dataset$sample_table column names! Please see the example data !")
			}
			if(!("SampleID" %in% colnames(sample_table))){
				sample_table$SampleID <- rownames(sample_table)
			}
			if(!("isIncluded" %in% colnames(sample_table))){
				sample_table$isIncluded <- 1
			}
			if(!("purterbID" %in% colnames(sample_table))){
				sample_table$purterbID <- 0
			}
			sample_table <- sample_table[, c("SampleID", "isIncluded", "Rep", "Time", "purterbID")]
			colnames(sample_table)[3:4] <- c("subjectID", "measurementID")
			
			abund_table <- dataset1$otu_table
			# only use high abundant species
			abund_table <- abund_table[apply(abund_table, 1, sum)/sum(abund_table) > filter_thres, ]
			
			# Now abund_table is a matrix, rows: sample, cols: features
			self$abund_table <- abund_table
			self$sample_table <- sample_table
		},
		#' @description
		#' Predict the biomass.
		#' 
		#' @param min_iter minimal number of iterations for the EM algorithm (default: 30)
		#' @param max_iter maximal number of iterations for the EM algorithm (default: 100)
		#' @param verbose print more information (default: TRUE)
		#' @param scaling median total biomass to scale all biomass data (default:10000)
		#' @param seed random seed used in BLASSO (default:NULL)
		#' @param ... parameters passed to \code{\link{EM}} function in beem package.
		#' @return res_biomass and res_param stored in the object.
		#' @examples
		#' t1$cal_biomass(min_iter = 50)
		cal_biomass = function(min_iter = 30, max_iter = 100, verbose = TRUE, scaling = 10000, seed = 1, ...)
			{
			
			if(!require(beem)){
				stop("beem package not installed! See https://github.com/ChiLiubio/mecodev for the installation !")
			}
			message("Please also cite the original article: Li et al. <doi:10.1186/s40168-019-0729-z>")
			message("An expectation-maximization algorithm enables accurate ecological modeling using longitudinal microbiome sequencing data. Microbiome, (2019) 7:118.\n")

			abund_table <- self$abund_table
			sample_table <- self$sample_table
			
			res <- EM(
				dat = abund_table, 
				meta = sample_table, 
				min_iter = min_iter,
				max_iter = max_iter,
				verbose = verbose,
				scaling = scaling,
				seed = seed,
				...
				)
			
			biomass <- biomassFromEM(res)
			param <- paramFromEM(res, abund_table, sample_table)
			
			self$res_biomass <- biomass
			self$res_param <- param
		},
		#' @description
		#' Extract network result from EM results.
		#' 
		#' @param sig significance threshold for the interaction evaluation (default: 0.05)
		#' @return res_network stored in the object.
		#' @examples
		#' t1$t1$cal_network()
		cal_network = function(sig = 0.05){
			if(!require(igraph)){
				stop("igraph package not installed !")
			}
			if(is.null(self$res_param)){
				stop("Please first run cal_biomass function !")
			}else{
				parms <- self$res_param
			}
			counts_mean <- rowMeans(self$abund_table)
			int <- parms[parms$parameter_type == 'interaction' & parms$source_taxon != parms$target_taxon, ]
			int.f <- int[int$significance > sig, 2:4]
			g <- graph.data.frame(int.f[, 1:2])
			# assign colors
			# V(g)$color <- annote[V(g)$name, ]$V4
			V(g)$size <- log(counts_mean[V(g)$name]) + 4
			E(g)$color <- ifelse(int.f$value > 0, "red", "green")
			E(g)$lty <- ifelse(int.f$value > 0, 1, 2)
			E(g)$width <- private$minmax(abs(int.f$value)) * 2 + 0.5
			self$res_network <- g
		}
	),
	private = list(
		minmax = function(x){
			(x-min(x))/(max(x)-min(x))
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)

