#' @title
#' Simulate gamma-beta diversity relationships under specific species distribution.
#'
#' @description
#' Simulate gamma-beta diversity relationships under specific species distribution, and plot the result
#' based on <doi:10.1038/s41467-020-19228-4>. 
#' The beta diversity is defined as the average distance-to-centroid value and
#' measured as the average distance (or compositional dissimilarity) from one sample to the centroid of the group.
#'
#' @export
trans_gamma <- R6Class(classname = "trans_gamma",
	public = list(
		#' @description
		#' This function is used to simulate gamma-beta diversity relationships under specific species distribution.
		#'
		#' @param dataset default NULL; microtable object; used for the observed pattern.
		#' @param group default NULL; a column name in sample_table; used as the different species pool for groups.
		#' @param method default "bray"; dissimilarity indices; see \code{\link{vegdist}} function and method parameter in vegan package;
		#' or "wei_unifrac" or "unwei_unifrac" in microtable$cal_betadiv().
		#' @param seed default 123; random seed used for the fixed random number generator for the repeatability.
		#' @return parameters in object.
		#' @examples
		#' \donttest{
		#' library(microeco)
		#' data(dataset)
		#' test1 <- trans_gamma$new(dataset = dataset, group = "Type", method = "bray")
		#' }		
		initialize = function(dataset = NULL, group = NULL, method = "bray", seed = 123
			){
			set.seed(seed)
			if(!is.null(group)){
				if(is.null(dataset)){
					stop("group provided, but no dataset provided !")
				}else{
					self$group <- group
				}
			}
			if(!is.null(dataset)){
				use_dataset <- clone(dataset)
				use_dataset$rep_fasta <- NULL
				use_dataset$taxa_abund <- NULL
				use_dataset$alpha_diversity <- NULL
				use_dataset$beta_diversity <- NULL
				self$dataset <- use_dataset
			}
			self$method <- method
		},
		#' @description
		#' Calculate observed gamma and beta diversity for each group.
		#' The beta diversity is defined as mean dispersion from the centroid based on the distance matrix.
		#' The gamma diversity is defined as the total species observed.
		#'
		#' @param sample_size default NULL; a numeric vector for the rarefied and uniform individual numbers in each sample; 
		#'   If null, use the observed data; If provided, use the rarefied data, e.g. c(500, 2000).
		#' @return res_observed in object.
		#' @examples
		#' \donttest{
		#' test1$cal_observed(sample_size = NULL)
		#' }
		cal_observed = function(sample_size = NULL){
			method <- self$method
			if(method %in% c("wei_unifrac", "unwei_unifrac")){
				unifrac <- TRUE
				beta_method <- NULL
			}else{
				unifrac <- FALSE
				beta_method <- method
			}
			if(is.null(self$dataset)){
				stop("Please provide the dataset when creates the object !")
			}else{
				use_data <- self$dataset
			}
			if(is.null(sample_size)){
				used_sample_size = "observed"
			}else{
				used_sample_size = sample_size
			}
			if(is.null(self$group)){
				stop("Please provide the group parameter when creates the object !")
			}else{
				all_groups <- unique(use_data$sample_table[, self$group])
			}
			
			res_beta <- matrix(0, nrow = length(used_sample_size), ncol = length(all_groups))
			res_gamma <- matrix(0, nrow = length(used_sample_size), ncol = length(all_groups))

			for(i in seq_along(used_sample_size)){
				beta_obs_mean <- c()
				gamma_obs <- c()
				for(j in all_groups){
					use_data_subset <- clone(use_data)
					use_data_subset$sample_table %<>% .[.[, self$group] == j, ]

					if(is.null(sample_size)){
						use_data_subset$tidy_dataset()
					}else{
						suppressMessages(use_data_subset$rarefy_samples(sample.size = used_sample_size[i]))
					}
					suppressMessages(use_data_subset$cal_betadiv(method = beta_method, unifrac = unifrac))
					beta_partition_bdisper <- vegan::betadisper(as.dist(use_data_subset$beta_diversity[[method]]), 
						group = rep(1, times = nrow(use_data_subset$sample_table)), type="centroid")
					beta_partition_distance <- mean(beta_partition_bdisper$distances)
					beta_obs_mean <- c(beta_obs_mean, beta_partition_distance)
					gamma_obs <- c(gamma_obs, nrow(use_data_subset$otu_table))
				}
				res_beta[i, ] <- beta_obs_mean
				res_gamma[i, ] <- gamma_obs
			}

			res_beta <- as.data.frame(res_beta)
			rownames(res_beta) <- used_sample_size
			colnames(res_beta) <- all_groups

			res_gamma <- as.data.frame(res_gamma)
			rownames(res_gamma) <- used_sample_size
			colnames(res_gamma) <- all_groups
			
			res_beta <- res_beta %>%
				cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
				reshape2::melt(id.vars = "Sample")

			colnames(res_beta)[2:3] <- c("Group", "beta")

			res_gamma <- res_gamma %>%
				cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
				reshape2::melt(id.vars = "Sample")

			res_data <- cbind.data.frame(res_beta, gamma = res_gamma[, 3])
			res_data <- dropallfactors(res_data, unfac2num = TRUE)
			res_data$Sample %<>% as.character
			self$res_observed <- res_data

			message('The result is stored in object$res_observed ...')			
		},
		#' @description
		#' Plot the observed result.
		#'
		#' @param x_axis_title default "Gamma diversity"; x axis title.
		#' @param y_axis_title default "Beta diversity"; y axis title.
		#' @param ... parameters pass to the function plot_scatterfit in \code{\link{trans_env}} Class.
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' test1$plot_observed(cor_method = "spearman")
		#' }
		plot_observed = function(
			x_axis_title = "Gamma diversity",
			y_axis_title = "Beta diversity",
			...
			){
			if(is.null(self$res_observed)){
				stop("Please first use cal_observed function to obtain result !")
			}else{
				res_data <- self$res_observed	
			}
			# create trans_env object to use scatter plot
			t1 <- trans_env$new(dataset = self$dataset, env_cols = 1)

			p <- t1$plot_scatterfit(
				x = res_data$gamma, 
				y = res_data$beta,
				x_axis_title = x_axis_title,
				y_axis_title = y_axis_title,
				...
				)
			p
		},
		#' @description
		#' Simulate gamma-beta diversity relationships under log-normal abundance distribution without any measured data.
		#'
		#' @param gamma_vect default seq(1, 10000, by = 200); a vector as gamma diversity.
		#' @param ind_vect default c(500, 1500, 3000, 5500, 10000); a vector as individuals per plot.
		#' @param ncom default 100; number of communities; how many communities or samples in the region or the studied species pool.
		#' @return res_simulation in object.
		#' @examples
		#' \donttest{
		#' test1$cal_simulation(ncom = 20, ind_vect = c(200, 1000, 2000))
		#' }
		cal_simulation = function(
			gamma_vect = seq(1, 10000, by = 200),
			ind_vect = c(500, 1500, 3000, 5500, 10000), 
			ncom = 100
			){

			#simulation study with a lognormal species abundance distribution
			#Generating a vector with 10000 elements on regional gamma diversity;
			#Gamma diversity range from 1 to 10000, with 200 increments between settings;
			self$gamma_vect <- gamma_vect

			# set the number of plots per region(ncom) and individuals per plot (ind_vect)
			self$ncom <- ncom
			self$ind_vect <- ind_vect
			method <- self$method
			if(method %in% c("wei_unifrac", "unwei_unifrac")){
				stop("For the simulation, unifrac index is not available! Please switch to others!")
			}

			# create a matrix to store all the simulated results
			beta_obs_matrix_mean <- matrix(0, nrow = length(ind_vect), ncol = length(gamma_vect))

			for (i in seq_along(ind_vect)) {
				beta_obs_mean <- NULL
				cat(paste0("## Runs: ind_vect at ", ind_vect[i], "\n"))
				for(j in seq_along(gamma_vect)){
					# cat(paste0("Runs: ind_vect at ", ind_vect[i], " gamma_vect at ", gamma_vect[j], "\n"))
					# obtain the simulated communities
					simulated_communities = NULL
					k = 0
					while(k < ncom){
						# use cv_abund = 1 temporarily. Little effect
						simulated_communities <- rbind(simulated_communities, private$sim_com(gamma_d = gamma_vect[j], alpha_d = ind_vect[i], cv_abund = 1))
						k <- k + 1
					}

					# calculate beta-diversity as the distance (or compositional dissimilarity) from an individual 
					# plot to the centroid of the group of all plots within a region (distance-to-centroid) using the multivariate method
					beta_partition <- vegdist(simulated_communities, method = method)
					if(mean(beta_partition) == 0){
						beta_partition.distance <- 0
					}else{
						beta_partition.bdisper <- betadisper(beta_partition, group = rep(1, times=nrow(simulated_communities)), type="centroid")
						beta_partition.distance <- mean(beta_partition.bdisper$distances)
					}
					beta_obs_mean <- c(beta_obs_mean, beta_partition.distance)
				}
				beta_obs_matrix_mean[i,] <- beta_obs_mean 
			}
			res_sim <- as.data.frame(beta_obs_matrix_mean)
			rownames(res_sim) <- ind_vect
			colnames(res_sim) <- gamma_vect
			
			self$res_simulation <- res_sim

			message('The result is stored in object$res_simulation ...')
		},
		#' @description
		#' Plot the simulation result.
		#'
		#' @param color_values colors used for presentation.
		#' @param color default "SampleID"; color mapping in the plot.
		#' @param show_point default TRUE; whether show the point.
		#' @param point_size default .3; point size value.
		#' @param point_alpha default .6; point alpha value.
		#' @param add_fitting default FALSE; whether add fitted line.
		#' @param ... parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).
		#' @return ggplot.
		#' @examples
		#' \donttest{
		#' test1$plot_simulation()
		#' }
		plot_simulation = function(
			color_values = RColorBrewer::brewer.pal(8, "Dark2"),
			show_point = TRUE,
			point_size = .3,
			point_alpha = .6,
			add_fitting = FALSE,
			...
			){
			if(is.null(self$res_simulation)){
				stop("Please first run cal_simulation function!")
			}else{
				res_sim <- self$res_simulation
			}

			plot_data <- res_sim %>%
				cbind.data.frame(Sample = rownames(.), ., stringsAsFactors = FALSE) %>%
				reshape2::melt(id.vars = "Sample")
			
			plot_data <- dropallfactors(plot_data, unfac2num = TRUE)
			plot_data$Sample %<>% as.character
			plot_data$Sample %<>% factor(., levels = self$ind_vect)
			
			p <- ggplot(plot_data, aes_string(x = "variable", y = "value", color = "Sample", group = "Sample")) + 
				scale_color_manual(values = color_values) +
				xlab("Gamma diversity") + 
				ylab("Beta diversity")
				
			if(show_point == T){
				p <- p + geom_point(alpha = point_alpha, size = point_size)
			}
			if(add_fitting == T){
				p <- p + geom_smooth(se = FALSE, ...)
			}else{
				p <- p + geom_line(...)
			}
			p <- p + theme_bw()
			
			p

		},
		#' @description
		#' Print the trans_gamma object.
		print = function() {
			cat("trans_gamma class:\n")
			if(!is.null(self$res_simulation)){
				cat(paste("Simulation result have", nrow(self$res_simulation), "rows and", ncol(self$res_simulation), "columns ...\n"))
			}
			if(!is.null(self$res_observed)){
				cat(paste("Observed results have been calculated ...\n"))
			}
			invisible(self)
		}
	),
	private = list(
		# This function is used to generate simulated community given gamma diversity and number of individuals per plot
		# cv_abund: Setting lognormal distribution parameters
		sim_com = function(gamma_d, alpha_d, cv_abund = 1)
		{
			# hypothesize lognormal species abundance distribution
			# abund_pool~LogNorm(mean = mean_abund, sd = sd_abund)

			mean_abund <- alpha_d/gamma_d
			# the variable coefficient setup 1
			sd_abund <- mean_abund * cv_abund #given cv_abund=1, sd_abund=mean_abund;

			# Given log(abund_pool)~Norm(mean=mu1, sd=sigma1), 
			# calculating the corresponding normal distribution parameters;
			sigma1 <- sqrt(log(sd_abund^2/mean_abund^2 + 1))
			mu1 <- log(mean_abund) - sigma1^2/2

			# simulating species abundance distribution from lognormal distribution;
			abund_pool <- rlnorm(n = gamma_d, meanlog = mu1, sdlog = sigma1)

			# calulating simulated relative abundance;
			# ranking and sorting from largest relative abundance to the smallest;
			rel_abund_pool <- abund_pool/sum(abund_pool)
			rel_abund_pool <- sort(rel_abund_pool, decreasing = T)

			names(rel_abund_pool) <- paste0("species",seq_along(rel_abund_pool))

			# generate the simulated community from the pool via random sampling with replacement;
			# relative abundance is used to weight the probability on sample selection; 

			sample_vec <- sample(x = names(rel_abund_pool), size = alpha_d, 
							   replace = TRUE, prob = rel_abund_pool)
			sample_vec <- factor(sample_vec, levels = names(rel_abund_pool))

			abund_local <- table(sample_vec)

			class(abund_local) <- c("integer")
			abund_local
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)
