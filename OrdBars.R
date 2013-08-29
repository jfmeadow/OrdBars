####################################
####################################
##
##  R code to create ordinations with error bars 
##   summarizing groups of points. 
##
##  James Meadow jfmeadow@gmail.com 
##
## Copy and paste the whole function into R console.
##	there are examples at the bottom of the script. 
##
##  Arguements: 
##	ord = ordination object of any class of nmds or any of the constrained ordinations.
##		If it is a constrained ordination (like capscale or cca), do this: 
##		example.plot <- plot(example.cca) # and then use example.plot as the ord.
##	clustering = has to be a factor. this is how points are grouped
##	labels = empty unless you put them in. You can put them in like this labels = c('group1', 'group2', ...)
##		or you can use the levels() function as in dune example below.
## 	x.lim, y.lim = you can leave this empty, it will automatically scale to the limits of the nmds points,
##		but that might make ugly clustering in the center, 
##		so you can treat this just as xlim and ylim in other plots 
##	lab.off.x, lab.off.y = this offsets the labels from the points. 
##		If you have labels, you will have to tweak this. 
##		Can be either one single value or can be as many values as there are levels.
##	bg = color of center of points.
##	col = color of error bars
## 		both of these color arguements can be either lists of colors or single colors for all
##  new = asking whether to build a new plot from scratch (if you are just adding points to an 
##		existing plot, this should be FALSE. TRUE is default
##	add = goes with 'new'. If you just want to add points to an existing plot, set these
##  	like this: new=FALSE, add=TRUE. Or if you want to build a plot without points (to put things 
##		under the points) like this: new=TRUE, add=FALSE. The default setting will
##		make a new plot and add points. 
##
##	Also, I made the axes numberless as a style issue. 
##	If you want numbers (ordination scores) on axis, 
##   	you would have to build a plot and use this: new=FALSE, add=TRUE
##
####################################

OrdBars <- function(ord, clustering, labels='', x.lim=NULL, y.lim=NULL, lab.off.x=0.2, lab.off.y=0.2, bg='gray70', col='gray30', new=TRUE, add=TRUE) {
		require(vegan)
		options(warn=-1)
		if(class(ord) == 'list') {class(ord) <- 'nmds'}
		if(class(ord) == 'nmds') {names(ord)[1] <- 'sites'}
		if('metaMDS' %in% class(ord)) {ord <- list('sites' = scores(ord))}
		if(any(table(clustering) == 0)) {
			stop('One or more empty levels in clustering')
			}
		caps.table <- data.frame(matrix(NA, nlevels(clustering), 9))
		row.names(caps.table) <- levels(clustering)
		names(caps.table) <- c('mean.x', 'mean.y', 'sd.x', 'sd.y', 
			'n', 'se.x', 'se.y', 'lab.x', 'lab.y')
		caps.table$mean.x <- tapply(ord$sites[, 1], clustering, mean)
		caps.table$mean.y <- tapply(ord$sites[, 2], clustering, mean)
		caps.table$sd.x <- tapply(ord$sites[, 1], clustering, sd)
		caps.table$sd.y <- tapply(ord$sites[, 2], clustering, sd)
		caps.table$n <- table(clustering)
		caps.table$se.x <- caps.table$sd.x/sqrt(caps.table$n)
		caps.table$se.y <- caps.table$sd.y/sqrt(caps.table$n)
		caps.table$lab.x <- caps.table$mean.x + lab.off.x
		caps.table$lab.y <- caps.table$mean.y + lab.off.y
		if (is.na(labels)) {
			labels <- row.names(caps.table)
		}
		options(warn=0)
		# plot
	if(new) {
		if(is.null(x.lim) & is.null(y.lim)) {
			plot(caps.table$mean.x, caps.table$mean.y, type='n', 
				xlim=c(range(ord$sites[, 1])), ylim=c(range(ord$sites[, 2])), 
				xaxt='n', yaxt='n', 
				xlab='', ylab='') }
			else {plot(caps.table$mean.x, caps.table$mean.y, type='n', 
				xlim=x.lim, ylim=y.lim, 
				xaxt='n', yaxt='n', 
				xlab='', ylab='') }
		if(class(ord) == 'nmds') {
			mtext('NMDS 1', side=1, adj=1)
			mtext('NMDS 2', side=2, adj=1) }
		else{	mtext(dimnames(ord$sites)[[2]][1], side=1, adj=1)
			mtext(dimnames(ord$sites)[[2]][2], side=2, adj=1)}
		}
	if(add) {
		arrows(caps.table$mean.x + caps.table$se.x, 
			caps.table$mean.y,
			caps.table$mean.x - caps.table$se.x, 
			caps.table$mean.y,
			code=3, angle=90, lwd=1, col=col, length=.05) 
		arrows(caps.table$mean.x, 
			caps.table$mean.y + caps.table$se.y,
			caps.table$mean.x, 
			caps.table$mean.y - caps.table$se.y,
			code=3, angle=90, lwd=1, col=col, length=.05) 
		points(caps.table$mean.x, caps.table$mean.y, 
			pch=21, cex=1, bg=bg)
		text(caps.table$mean.x + lab.off.x, 
			caps.table$mean.y + lab.off.y, 
			labels)
		}		
		invisible(caps.table)
}


###################################
##
## this is an example using the dune data embedded in vegan
##

# require(vegan)
# data(dune); data(dune.env)
# dune.nmds <- isoMDS(vegdist(dune, 'bray'))
# OrdBars(dune.nmds, dune.env$Use, labels=levels(dune.env$Use), lab.off.x=.12, lab.off.y=.04)

## put into object to get coordinates table for adding other ornaments

# dune.ordbars <- OrdBars(dune.nmds, dune.env$Use)
# dune.ordbars

####################################
