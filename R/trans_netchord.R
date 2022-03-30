#' @title
#' Sum and plot the links number of the network.
#'
#' @description
#' Sum and plot the links number from one taxa to another or in the same taxa in the network.
#' The input dataset must be a trans_network object.
#'
#' @export
trans_netchord <- R6Class(classname = "trans_netchord",
	public = list(
		#' @description
		#' This function is used to sum the links number from one taxa to another or in the same taxa, for example, at Phylum level.
		#' This is very useful to fast see how many nodes are connected between different taxa or within the taxa.
		#'
		#' @param dataset default NULL; trans_network object.
		#' @param taxa_level default "Phylum"; taxonomic rank.
		#' @return res_sum_links_pos and/or res_sum_links_neg in object.
		#' @examples
		#' \dontrun{
		#' test1 <- trans_netchord$new(dataset = trans_network_object, taxa_level = "Phylum")
		#' }
		initialize = function(dataset = NULL, taxa_level = "Phylum")
			{
			if(is.null(dataset)){
				stop("The dataset should not be NULL!")
			}
			if(!require(igraph)){
				stop("igraph package not installed!")
			}
			taxa_table <- dataset$use_tax
			network <- dataset$res_network
			link_table <- data.frame(t(sapply(1:ecount(network), function(x) ends(network, x))), label = E(network)$label, stringsAsFactors = FALSE)
			# check the edge label
			if(! any(c("+", "-") %in% link_table[, 3])){
				stop("Please check the edge labels! The labels should be + or - !")
			}
			if("+" %in% link_table[, 3]){
				link_table_1 <- link_table[link_table[, 3] %in% "+", ]
				self$res_sum_links_pos <- private$sum_link(taxa_table = taxa_table, link_table = link_table_1, taxa_level = taxa_level)
				message('The positive results are stored in object$res_sum_links_pos ...')
			}else{
				message('No positive edges found ...')
			}
			if("-" %in% link_table[, 3]){
				link_table_1 <- link_table[link_table[, 3] %in% "-", ]
				self$res_sum_links_neg <- private$sum_link(taxa_table = taxa_table, link_table = link_table_1, taxa_level = taxa_level)
				message('The negative results are stored in object$res_sum_links_neg ...')
			}else{
				message('No negative edges found ...')
			}
		},
		#' @description
		#' Plot the summed linkages among taxa using chorddiag package <https://github.com/mattflor/chorddiag>.
		#'
		#' @param plot_pos default TRUE; If TRUE, plot the summed positive linkages; If FALSE, plot the summed negative linkages.
		#' @param plot_num default NULL; number of taxa presented in the plot.
		#' @param color_values default NULL; If not provided, use default.
		#' @return chorddiag plot
		#' @examples
		#' \donttest{
		#' test1$plot_sum_links(plot_pos = TRUE, plot_num = 10)
		#' }
		plot_sum_links = function(plot_pos = TRUE, plot_num = NULL, color_values = NULL){
			if(plot_pos == T){
				if(is.null(self$res_sum_links_pos)){
					stop("No res_sum_links_pos found!\n")
				}else{
					use_data <- self$res_sum_links_pos
				}
			}else{
				if(is.null(self$res_sum_links_neg)){
					stop("No res_sum_links_neg found!\n")
				}else{
					use_data <- self$res_sum_links_neg
				}
			}
			if(!is.null(plot_num)){
				if(plot_num > ncol(use_data)){
					message("The plot_num provided is larger than the total taxa number! Use the taxa number instead of it ...")
					plot_num <- ncol(use_data)
				}
				use_data %<>% .[1:plot_num, 1:plot_num]
			}
			if(is.null(color_values)){
				if(nrow(use_data) <= 14){
					groupColors <- c(RColorBrewer::brewer.pal(12, "Paired"), "#FDE0EF", "#C51B7D")
				}else{
					groupColors <- unname(randomcoloR::distinctColorPalette(nrow(use_data)))
				}
			}
			chorddiag::chorddiag(use_data, groupColors = groupColors)
		},
		#' @description
		#' Print the trans_netchord object.
		print = function() {
			cat("trans_netchord class:\n")
			if(!is.null(self$res_sum_links_pos)){
				cat("Positive summed links have been calculated! \n")
			}
			if(!is.null(self$res_sum_links_neg)){
				cat("Negative summed links have been calculated! \n")
			}
			invisible(self)
		}
	),
	private = list(
		sum_link = function(taxa_table, link_table, taxa_level){
			# first obtain the taxa names
			all_names <- taxa_table[rownames(taxa_table) %in% unique(c(link_table[,1], link_table[,2])), ] %>%
				{table(.[, taxa_level])} %>%
				sort(., decreasing = TRUE) %>% 
				rownames
			com_group <- expand.grid(all_names, all_names)
			colnames(com_group) <- c("C1", "C2")
			# assign rownames irrespective of the order
			rownames(com_group) <- apply(com_group, 1, function(x) paste0(x, collapse = "-"))
			# get the unifrom combined name without regard to the order
			com_group$uni_name <- apply(com_group, 1, function(x) paste0(sort(x), collapse = "-"))
			com_group1 <- com_group[, -c(1,2), drop = FALSE]
			res <- link_table
			# use taxa name to replace the species name
			res[, 1] <- taxa_table[res[, 1], taxa_level]
			res[, 2] <- taxa_table[res[, 2], taxa_level]
			res$pname <- paste(res[, 1], res[, 2], sep = "-")
			res %<>% dplyr::group_by(pname) %>% 
				dplyr::summarise(count = dplyr::n()) %>%
				as.data.frame(stringsAsFactors = FALSE)
			res <- dplyr::left_join(res, rownames_to_column(com_group1), by = c("pname" = "rowname")) %>%
				dplyr::group_by(uni_name) %>% 
				dplyr::summarise(sum_count = sum(count)) %>%
				as.data.frame(stringsAsFactors = FALSE)
			res <- dplyr::left_join(res, com_group, by = c("uni_name" = "uni_name"))
			res <- reshape2::dcast(res, C1~C2, value.var = "sum_count") %>%
				`row.names<-`(.[,1]) %>%
				.[, -1, drop = FALSE] %>%
				.[all_names, all_names] %>%
				as.matrix
			res[is.na(res)] <- 0			
			res
		}
	),
	lock_objects = FALSE,
	lock_class = FALSE
)

