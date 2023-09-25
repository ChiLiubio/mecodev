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
		#' @return microtable object stored in the object.
		#' @examples
		#' data(gut_microb_ts)
		#' t1 <- trans_ts$new(dataset = gut_microb_ts, filter_thres = 0)
		initialize = function(dataset = NULL, filter_thres = 0.001)
			{
			if(is.null(dataset)){
				stop("No dataset provided !")
			}
			use_dataset <- clone(dataset)
			use_dataset$tidy_dataset()
			
			sample_table <- use_dataset$sample_table
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
			
			abund_table <- use_dataset$otu_table
			# only use high abundant species
			abund_table <- abund_table[apply(abund_table, 1, sum)/sum(abund_table) > filter_thres, ]
			# Now abund_table is a matrix, rows: sample, cols: features
			use_dataset$otu_table <- abund_table
			use_dataset$sample_table <- sample_table
			use_dataset$tidy_dataset()
			self$dataset <- use_dataset
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
		#' \dontrun{
		#' t1$cal_biomass(min_iter = 50)
		#' }
		cal_biomass = function(min_iter = 30, max_iter = 100, verbose = TRUE, scaling = 10000, seed = 1, ...)
			{
			
			if(!require(beem)){
				stop("beem package not installed! See https://github.com/ChiLiubio/mecodev for the installation !")
			}
			message("Please also cite the original article: Li et al. <doi:10.1186/s40168-019-0729-z>")
			message("An expectation-maximization algorithm enables accurate ecological modeling using longitudinal microbiome sequencing data. Microbiome, (2019) 7:118.\n")

			abund_table <- self$dataset$otu_table
			sample_table <- self$dataset$sample_table
			
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
		#' @param sig default 0.05 significance threshold for the interaction evaluation (default: 0.05).
		#' @param add_taxa_name default "Phylum"; one or more taxonomic rank name; used to add taxonomic rank name to network node properties.
		#' @return trans_network object
		#' @examples
		#' \dontrun{
		#' # t2 is a trans_network object used for network analysis
		#' t2 <- t1$cal_network()
		#' }
		cal_network = function(sig = 0.05, add_taxa_name = "Phylum"){
			if(!require(igraph)){
				stop("igraph package not installed !")
			}
			if(is.null(self$res_param)){
				stop("Please first run cal_biomass function !")
			}else{
				parms <- self$res_param
			}
			counts_mean <- rowMeans(self$dataset$otu_table)
			int <- parms[parms$parameter_type == 'interaction' & parms$source_taxon != parms$target_taxon, ]
			int.f <- int[int$significance > sig, 2:4]
			network <- graph.data.frame(int.f[, 1:2])

			E(network)$label <- ifelse(int.f$value > 0, "+", "-")
			E(network)$weight <- abs(int.f$value)
			E(network)$width <- private$minmax(abs(int.f$value)) * 2 + 0.5
			if(!is.null(add_taxa_name)){
				if(!is.null(self$dataset$tax_table)){
					for(i in add_taxa_name){
						if(i %in% colnames(self$dataset$tax_table)){
							network <- set_vertex_attr(network, i, value = V(network)$name %>% 
								self$dataset$tax_table[., i] %>% 
								gsub("^.__", "", .))
						}else{
							message("Skip adding taxonomy: ", i, " to node as it is not in colnames of tax_table ...")
						}
					}
				}else{
					message('Skip adding taxonomy to node as tax_table is not found ...')
				}
			}
			t1 <- trans_network$new(dataset = self$dataset, cor_method = NULL)
			# from microeco v0.12.1
			V(network)$RelativeAbundance <- t1$data_relabund[V(network)$name]
			t1$res_network <- network
			t1
		}
	),
	private = list(
		minmax = function(x){
			(x - min(x))/(max(x) - min(x))
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)

