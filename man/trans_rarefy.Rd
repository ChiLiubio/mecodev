% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_rarefy.R
\name{trans_rarefy}
\alias{trans_rarefy}
\title{Rarefaction and plotting.}
\description{
Rarefaction based on the microtable rarefy_samples function and plotting based on the ggplot2.
}
\examples{

## ------------------------------------------------
## Method `trans_rarefy$new`
## ------------------------------------------------

\donttest{
library(microeco)
data(dataset)
t1 <- trans_rarefy$new(dataset = dataset, depth = c(0, 10, 50, 400, 800))
}

## ------------------------------------------------
## Method `trans_rarefy$plot_rarefy`
## ------------------------------------------------

\donttest{
t1$plot_rarefy(color = "Group")
}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-trans_rarefy-new}{\code{trans_rarefy$new()}}
\item \href{#method-trans_rarefy-plot_rarefy}{\code{trans_rarefy$plot_rarefy()}}
\item \href{#method-trans_rarefy-print}{\code{trans_rarefy$print()}}
\item \href{#method-trans_rarefy-clone}{\code{trans_rarefy$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_rarefy-new"></a>}}
\if{latex}{\out{\hypertarget{method-trans_rarefy-new}{}}}
\subsection{Method \code{new()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_rarefy$new(dataset = NULL, alphadiv = "Shannon", depth = NULL, ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{dataset}}{the object of \code{\link{microtable}} Class.}

\item{\code{alphadiv}}{default "Shannon"; alpha diversity measurement used for the rarefaction; see microtable$cal_alphadiv for all the measurement.}

\item{\code{depth}}{default NULL; an integer vecotr used for the rarefying.}

\item{\code{...}}{parameters passed to \code{rarefy_samples} function of \code{microtable} class, except the sample.size parameter.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
res_rarefy stored in the object.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
library(microeco)
data(dataset)
t1 <- trans_rarefy$new(dataset = dataset, depth = c(0, 10, 50, 400, 800))
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_rarefy-plot_rarefy"></a>}}
\if{latex}{\out{\hypertarget{method-trans_rarefy-plot_rarefy}{}}}
\subsection{Method \code{plot_rarefy()}}{
Plotting the rarefied result.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_rarefy$plot_rarefy(
  color_values = RColorBrewer::brewer.pal(8, "Dark2"),
  color = "SampleID",
  show_point = TRUE,
  point_size = 0.3,
  point_alpha = 0.6,
  add_fitting = FALSE,
  x_axis_title = "Sequence number",
  y_axis_title = NULL,
  show_legend = TRUE,
  show_samplename = FALSE,
  samplename_size = 3,
  samplename_color = "grey30",
  ...
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{color_values}}{colors used for presentation.}

\item{\code{color}}{default "SampleID"; color mapping in the plot.}

\item{\code{show_point}}{default TRUE; whether show the point.}

\item{\code{point_size}}{default .3; point size value.}

\item{\code{point_alpha}}{default .6; point alpha value.}

\item{\code{add_fitting}}{default FALSE; whether add fitted line.}

\item{\code{x_axis_title}}{default "Sequence number"; x axis title.}

\item{\code{y_axis_title}}{default NULL; default NULL represents the measure used.}

\item{\code{show_legend}}{default TRUE;    whether show the legend in the plot.}

\item{\code{show_samplename}}{default FALSE; whether show the sample name in the plot.}

\item{\code{samplename_size}}{default 3; the sample name text size. Only available when show_samplename is TRUE.}

\item{\code{samplename_color}}{default "grey70"; sample name text color. Only available when show_samplename is TRUE.}

\item{\code{...}}{parameters pass to ggplot2::geom_line (when add_fitting = FALSE) or ggplot2::geom_smooth (when add_fitting = TRUE).}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
ggplot.
}
\subsection{Examples}{
\if{html}{\out{<div class="r example copy">}}
\preformatted{\donttest{
t1$plot_rarefy(color = "Group")
}
}
\if{html}{\out{</div>}}

}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_rarefy-print"></a>}}
\if{latex}{\out{\hypertarget{method-trans_rarefy-print}{}}}
\subsection{Method \code{print()}}{
Print the trans_rarefy object.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_rarefy$print()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-trans_rarefy-clone"></a>}}
\if{latex}{\out{\hypertarget{method-trans_rarefy-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{trans_rarefy$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
