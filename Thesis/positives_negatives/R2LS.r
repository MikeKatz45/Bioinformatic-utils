#source("~/Dropbox/AllScriptsApp/Scripts/R2LS.r")

#Version Mayo 21 2013 : bug quantize

#Version Marzo 1 2012 : bug quantize

#Version Oct 27 2011

#Version Dec 3rd 2010

# Procesing PubMed Results
# Replace: ([^\n])\n([^\n]) with $1 $2
# Replace: \n\n with \t

#cluster.exemplars <- function(data, threshold = 0.8, similarity = function(x) cor(x, method="spearman"), 
	#merge = function(m) apply(m, 1, mean),reorder=NA, verbose=FALSE) 

## Saca un histograma gráfico en 2 dimensiones
# plot.hist2d(x, y, nbins=50, col = c("white",heat.colors(16)), ...) 


## ALinea 2 secuencias por prog dinamica
#align.score <- function(seq1, seq2, gap=0, igap=0, r=0, match=1, mismatch=0,ret.scores=FALSE) 
# obtiene el consenso de 2 secuencias
#best.align <- function(seq1, seq2, ret.scores=TRUE, ret.pos=FALSE, ...)  ## ... are parameters for align.score
# consensus <- function(s)

# calcula la frecuencia de oligos de long len de cada secuencia s y regresa una matriz
# oligosignature <- function(s, len=6)  

# ls.usage <- function(xls=ls(envir=.GlobalEnv))

# Calcula el p-value u otras valores en un modelo de cox
#p.cox <- function(data, time, status, method="breslow", value=c("logtest","sctest","waldtest","rsq","coefficient")[1], attr=c("pvalue", "rsq")[1])

# hace blast p al ncbi
# ncbi.blastp

## Read LARGE FILES slow but more securely (memory friendly)
# soft.read.delim <- function(xFile, row.estimate=NULL, chunk.size=1000, wc="wc -l", skip=0, header=TRUE, sep="\t", nrows=-1, comment.char="")

## Convert data (matrix, dataframe or vector) to values quantized in an interval
# quantize <- function(data, n=2, min=min(data, na.rm=TRUE), max=max(data, na.rm=TRUE)) 

# uniformize <- function(x)
# scale.mad <- function(x, center=TRUE, scale=TRUE)

# correct <- function(x, na=median, nan=na, inf=TRUE, neg.inf=min, pos.inf=max, null=x.na, outlier=TRUE, outlier.min=0.05, outlier.max=.95) 

# similarity <- function(x, y=NULL) 

# redgreenblue <- function(col=-3, nc=ifelse(length(col)==2,16,8), len=NULL, values=NULL) 

# p.overlap <- function(balls, white, drawn, observed.white, lower.tail=FALSE) 

# cumaverage <- function(x) 

# copa <- function(x, q=NULL) 

# outlier.sum <- function(x, group=NULL, ref=levels(group)[1], f=NULL, f2=NULL) 

# min.t <- function(x, k=length(x)*.1, parts=max(length(x)-k*2+1,1), index=FALSE, test=t.test, ...) 

# rlog <- function(v, func=log2, offset=2, cutoff=0) 

# rantilog <- function(v, func=function(v) 2^v, offset=2, cutoff=0) 

# glog <- function(y, lambda=1, func=log2) 

# gantilog <- function(x, lambda=1, func=function(x) 2^x) 

# apply.df <- function(x, margin, fun) 

# circularize <- function(x, len) 

# plot.densities <- function(m, m2, m3, m4, lty=rep(1:6, each=length(palette())), legends=TRUE, 
# 	legend.pos="topleft", legend.cols=1, legend.cex=1, horiz=FALSE, xlim=if(horiz) c(mny,mxy) else 	c(mnx,mxx),
# 	ylim=if(!horiz) c(mny,mxy) else c(mnx,mxx), pch=NULL, type="l", npch=20, col=1:length(d), normalize=FALSE, 
# 	q=NULL, overall=length(d) > 4, ...) 

# plot.ma <- function(m, m2, m3, m4, gap=0, pch=20, xaxt="n",yaxt="n", normalise=FALSE, ...) 

# plot.hist <- function(m, m2, m3, m4, col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), 
# 	separate=FALSE, beside=TRUE, lines=FALSE, legends=TRUE, legend.cex=1, legend.pos="topleft", legend.cols=1, 
# 	log="", hollow=FALSE, border=col, labels=FALSE, srt=90, cex=1, 
# 	ylim=c(0,(if(beside) max(bp) else max(apply(bp,1,sum)))*y.adj), density=FALSE, round.mids=NULL, 
# 	normalize=NULL, test=FALSE, test.func=dist.chisq.test, round.test=8, y.adj=1.25, names.arg=NULL, total=FALSE, 
# 	cumulative=FALSE, roundfunc=round, ...) 

# histq <- function(x, q=c(.5, .25, .75, .05, .95, 0, 1), col=c(1,rep(2:10,times=rep(2,9)))+1, 
# 	lty=c(1,rep(2:10,times=rep(2,9))), cex=.75, ...) 

# histm <- function(x, m=c(mean(x), mean(x)+sd(x), mean(x)-sd(x), mean(x)+sd(x)*2, mean(x)-sd(x)*2), 
# 	col=c(1,rep(2:10,times=rep(2,9)))+1, lty=c(1,rep(2:10,times=rep(2,9))), cex=.75, ...) 

# plot.hist2 <- function(m, m2, m3, m4, las=0, base=1, col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 
# 	1:length(m), lty=rep(1:6, each=length(palette())), xaxt="s", yaxt="s", xlab="", ylab=if(!is.null(normalize) && 
# 	any(normalize)) "Normalized Frequency" else if(density) "Density" else "Frequency", separate=FALSE, 
# 	legends=TRUE, legend.cex=1, legend.pos="topleft", legend.cols=1, log="", hollow=FALSE, border=col, 
# 	labels=FALSE, srt=90, cex=1, ylim=c(0,max(bp)*y.adj), density=FALSE, splice=0, shading=-1, round.mids=NULL, 
# 	normalize=NULL, test=FALSE, test.func=dist.chisq.test, round.test=8, y.adj=1.25, angle=45, total=FALSE, 
# 	cumulative=FALSE, roundfunc=round, ...) 

# plot.hist3D <- function(m, m2, m3, m4, col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), 
# 	separate=FALSE, beside=TRUE, lines=FALSE, legends=FALSE, legend.cex=1, legend.pos="topleft", legend.cols=1, 
# 	log="", hollow=FALSE, border=col, labels=FALSE, srt=90, cex=1, ylim=c(0,(if(beside) max(bp) else max(apply(bp,
# 	1,sum)))*y.adj), density=FALSE, round.mids=NULL, normalize=NULL, test=FALSE, test.func=dist.chisq.test, 
# 	round.test=8, y.adj=1, ysep=.1, xsep=0.5, bw=0.9, las=1, names.arg=TRUE, xaxis.adj=if (length(names.arg) && 
# 	names.arg[1]) 4 else 0, main="", xlab="", ylab="Frequency", sub="", roundfunc=round, ...) 

# dist.chisq.test <- function(a, b, count=NULL, breaks=seq(min(a,b), max(a,b), len=12), simulate.p.value=FALSE, 
# 	B=2000) 

# plot.qq <- function(m, m2, m3, m4, las=0, col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), 
# 	pch=1:length(col), lty=rep(1:6, each=length(palette())), legends=TRUE,  legend.pos="topleft", legend.cols=1, 
# 	legend.cex=1, len=NULL, type="o", ...) 

# na.impute <- function(m, n=max(3,ncol(m)/10), imputfunc=mean, distfunc = dist, method=c("static","dynamic"), 
# 	verbose=TRUE, ...) 

# plot.xy.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),breaks=20,
# 	breaks.h=NULL,breaks.v=NULL,widths=c(5,2),heights=c(2.5,5),smooth=FALSE,smooth.col="red",
# 	smooth.lty=1,smooth.type="l",smooth.pch=1,smooth.lwd=1,restore.layout=TRUE,log="",hook.func=function(x,y) 
# 	0, ...) 

# plot.box.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),breaks=20,breaks.h=NULL,breaks.v=NULL,widths=c(5,2),heights=c(2.5,5),smooth=FALSE,smooth.col="red",smooth.lty=1,smooth.type="l",smooth.pch=1,smooth.lwd=1,restore.layout=TRUE,overlap=TRUE,pch=20,pcol=rgb(.9,.9,.9),lwd=if (overlap) 2 else 1,hook.func=function(x,y) 0, ...) 

# matplot.xy.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),
# 	breaks=20,widths=c(5,2),heights=c(2.5,5),smooth=FALSE,smooth.col="red",smooth.lty=1,smooth.type="l",
# 	smooth.pch=1,smooth.lwd=1,...) 

# plot.2scales <- function(func=hist, ..., main=NULL, xlab=NULL, xaxt=NULL, ylim1=c(.9,1), ylim2=c(0,.25), 
# 	heights=c(0.25,0.75), axes=list(c(2),if(is.null(xaxt) || xaxt=="n") c(2) else c(1,2)), axes1=FALSE, 
# 	axes2=FALSE, panel1=NULL, panel2=NULL, restore.par=TRUE) 

# plot.2scales2 <- function(heights=c(0.25,0.75), panel1=NULL, panel2=NULL) 

# plot.bars <- function(xt, space=0.25, adj=c(0.5, -0.25), labels=paste(round(xt*100/sum(xt),labels.round),"%",sep=""), 
# 	cex.text=1, ylim=c(0, max(xt)*1.10), ...) 

# plot.two.bars <- function(h1, h2, weight=.8, space.center=1, space.left=mean(h1)/4, space.right=mean(h2)/4, col=1:length(h1), border=1:length(h1), labels=TRUE, cex=-1, font=1, mar=c(0,0,3,0), zeroes=TRUE, main="h1 vs h2", center=FALSE, maxHeight=max(c(h1,h2))) 
#

# plot.signatures <- function(mat, xaxis, type="l", pch=ifelse(type=="l",-1,20), splice=0, lwd=1, col=1:ncol(mat), 
# 	lty=rep(1:6, each=length(palette())), legends=FALSE, legend.pos="topleft", legend.cex=1, legend.cols=1, 
# 	labels=FALSE, labels.dig=NULL, add=FALSE, ...) 

# plot.profiles.msd <- function(mat, xaxis,  splice=0, ...) 

# plot.profile.range <- function(x, m, breaks=20, q=0.9, 

# smooth.window <- function(v, size=1, left=size, right=size, func=mean) 

# smooth.window.matrix <- function(mat, ...) apply(mat, 2, smooth.window, ...)

# smooth.window.data.frame <- function(mat, ...) apply(data.matrix(mat), 2, smooth.window, ...)

# smooth.y.by.x.window <- function(x, y, size=max(2,length(x) / 1000), func=mean) 

# smooth.y.by.x.hist <- function(x, y, breaks=100, func=mean) 

# predict.smooth.y.by.x.hist <- function(w, z, x) 

# colourise.dendrogram <- function(d, use.order=TRUE, ...) 

# plot.hclust.colours <- function(data, distfun=dist, hclustfun=hclust, label.list=list(lab.col=2, lab.font=1, 
# 	lab.cex=1), labels.sorted=FALSE, hang=0.1, ...) 

#plot.hclust.colored.clusters <- function(data, distfun=dist, hclustfun=hclust, h=NULL, k=NULL, 
#	hang=NULL, 
#	label.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1, lab.cex=0.5),
#	branch.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1),
#	upper.attr=list(col="#888888", pch=1, lty=1, lwd=1.5, cex=1), 
#	silhouettes=NULL, sil.labadj=FALSE, 
#	sil.off=ifelse(sil.labadj,ifelse(sil.image,0.05,0.01),0.3), 
#	sil.size=0.7, ylab="Height", sil.ylab="Silhouette", sil.col=redgreenblue(-1), 
#	sep=TRUE, sil.image=TRUE, 
#	side.factors=NULL, factors.size=1, factors.cex=1,...) {


# filterize <- function(x, q1, q2, func=NULL)  

# my.quantile <- function(v, probs, na.rm=FALSE) 

# quantile.normalization.old <- function(e.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, 
# 	m.func=mean) 

# quantile.normalization <- function(l.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, 
# 	m.func=mean, qmin=min(avg), qmax=max(avg), qmid=1, qshrink=NULL) 

# pad <- function(chv, char=if(is.numeric(chv)) "0" else " ", pos=-1, len=max(nchar(chv))) 

# trim <- function(ch) sub(" +$","",sub("^ +","",ch))
# trim2 <- function(ch,char=" ") sub(paste(char,"+$",sep=""),"",sub(paste("^",char,"+",sep=""),"",ch))
# rtrim <- function(ch) sub(" +$","",ch)

# ltrim <- function(ch) sub("^ +","",ch)

# plot.confront.one.pair <- function(l, x.adj=0, col.adj=0, test=t.test, dec=8, ...) 

# plot.confront.pair <- function(l1, l2, l3, l4, main="", xlab="Group", ylab="Value", 

# contingency <- function(x)  

# plot.confront <- function(dframe, cls, main="", xlab="Group", ylab="Value", 

# medline.coocurrance <- function(var.terms, fix.terms="") 

# medline.apply <- function(var.terms, func=medline.count, fix.terms="", pause.by=if (length(var.terms) > 99) 10 
# 	else 20, ...) 

# medline.simplify.geneName <- function(x) 

# medline.generalize.geneName <- function(x) 

# medline.query <- function (..., program="esearch.fcgi", db="pubmed", parameter="&rettype=count", collapse="+AND+", verbose=FALSE, baseUrl=getOption("serviceUrl.entrez"), pause=FALSE)

# medline.count <- function (...) 

# medline.abstracts <- function (..., abstracts=TRUE) 

# medline.pmids <- function (..., retmax=NULL)

# medline.pauseBetweenQueries <- function (peak=15, offpeak=3)

# names2index <- function(access.names, index.names, by=10000) 

# list.2.data.frame <- function(x, func=cbind, transfunc=function(x) x) 

# list.2.matrix <- function(x) 

# left <- function(x, n) 

# right <- function(x, n) 

# mid <- function(x, n) 

# moda <- function(x) 

# extract the first pattern in each position from a string vector
# extract <- function(pattern, stvector) 

# extract all patterns in each position from a string vector
# extract.all <- function(pattern, st, ret.num=FALSE) 

# vioplot <- function(x,range=1.5,h=NULL,ylim=NULL,names=NULL, col=1:length(datas), ...)

# gradient.colours <- function(vrange, ncolours = 16, breaks, c.breaks) 

# align.text.matches <- function(pattern, sts, nchars=80, left.chars=round((nchars-nchar(pattern))/2,0), marks="[]", 
# 	ends="...", position=FALSE, filter=TRUE, ret.all=FALSE) 
# align.score <- function(seq1, seq2, gap=0, igap=0, r=0, match=1, mismatch=0,ret.scores=FALSE) 

# text.similitude <- function(a, b, ignore.case=TRUE, sep=" ", ret.int=FALSE, uniq=TRUE, aprox=FALSE, ...) 

# catf <- function(...) 

# stabilizedICA <- function(data, nic=5, mic=10, times=100, q=0.95, average=TRUE, ret.all=FALSE, debug=FALSE) 

# p.information.content <- function(values, classes, wclass, sorted=FALSE, large=TRUE, V=length(values)) 

# p.yeohchi2 <- function(values, classes, wclass, sorted=FALSE, V=length(values)) 

# p.bimodal <- function(x, boxcoxify = TRUE, breaks="Sturges", constr=mixconstr(), emsteps=3, df=6, ret.all = FALSE, 
# 	modals=2, ...) 

# de.test <- function(x, classes, test=c("all", "ttest", "kolmogorov", "kruskal", "ftest", "welch", "cochran", 
# 	"wilcoxon", "copa", "s2n", "osum", "sam", "infcont", "yeohchi2", "pca", "affinity", "sliding","bimodal", "brown"), 
# 	nperms=100, ret.all = FALSE, debug=FALSE, except=c("affinity","bimodal"), ...) 

# affinity.propagation.clustering <- function(s = NULL, data = NULL, similarityfun = function(x) cor(x, 
# 	method="pearson"), max.iter=100, lam = 0.5, stable = 10, save.memory=nrow(s) > 5000, debug=FALSE) 

# list.extract <- function(xlist, xname, unlist=TRUE) 

# list.bind <- function(xlist, byrow=1) 

# venn <- function(id, category, cutoff=1, duplicates = FALSE, tab, func=function(x) sum(x) >= cutoff, main, 
# 	labels=FALSE,

# incidence.table <- function(id, category, names = NULL, cutoff = 1, 

# plot.chars <- function(cex=1,col=1:130) 

# plot.vertical.densities <- function(x.round, cls, pch=c(20,17,18,15,1:25), col=1:nl, 
# 	ylim=c(min(x.round),max(x.round)), xlim=c(0,nl), x.levels=levels(cls),xlab="Groups",nperline=10,cex=1,...) 

# plot.contour <- function(x, y, z, xbreaks="Sturges", ybreaks="Sturges", ...) 

# rad2deg <- function(rad) 

# deg2rad <- function(deg) 

# polar2xy <- function(z, rad) 

# plot.radar <- function(m, alfa.off=90, minpoint=0.5, col=1:s, pch=1:s, lty=rep(1,s), lwd=rep(1,s), 
# 	type=rep("l",s),

# filter.by.variation <- function(data.n, data.sd=NULL, data.m=NULL, draw=TRUE, pch=20, cex=0.5, use=c("cv","sd"), 
# 	filter.th = 0.1, minq=0.025, maxq=0.975, ...) 

# plot.pca <- function(data, classes=factor(1:ncol(data)), col=1:nlevels(classes), npc=4, gap=0, pch=1:nlevels(classes), labels=rownames(data), log="")

# remove parts of the string that are repeated in a "fraction" (default=100%=all) of names
#simplify.names <- function(x, pattern="[^A-Za-z0-9]", fraction=1, ...) 

#plot a string "rotated" within a rectangle estimating the font size
#plot.string.rot <- function(x,y,xx,yy, stv, col=rep(1,length(stv)), srt=0, cex.min=.1, cex.max=3, word.sep=" ", plot=TRUE, adj=c(0, 0), f.adj=1.05, space.x=1, space.y=1.3, rect=FALSE) 
#plot.string <- function(x0=NULL, y0=NULL, xx=NULL, yy=NULL, x, sep=ifelse(rotation,"  "," "), max.chars=40,  precision.factor=1, col=NULL, add=TRUE, wrap=TRUE, formatter=format, cex.min=.1, cex.max=10.1, autoscale=missing(cex), cex=cex.max, vfactor=1.25, srt=c(0,90,180,270)[1], ...) 


# minmax <- function(x, min=0, max=1) 
# qminmax <- function(x, min=0, max=1)

# plot.text <- function(x, space=5, max.chars=40, autoscale=missing(cex), cex=1, precision.factor=1, col=NULL, box.col=NULL, add=FALSE, x0=NULL, y0=NULL, xx=NULL, yy=NULL, formatter=format, hlines=8, vlines=8, hadj=NULL, ...)

# plot.xy.discrete <- function(x, y, col=1, pch=19, cex=1, jitters=TRUE, use.swarm="none", dec=2, cex.labels=0.75, df=3, lining=min(length(unique(x)),length(unique(y))) < 20, lining.col=8, xlim=NULL, ylim=NULL, xaxt=ifelse(is.factor(x),"n","s"), yaxt=ifelse(is.factor(y),"n","s"), xper=TRUE, yper=TRUE, allper=TRUE, alln=FALSE, top.cex.axis=1, right.cex.axis=1, in.circles=trunc, ...)

# plot.pairs <- function (x=1:1, y=1:1, func=function(...) plot(..., xlab="", ylab="", axes=FALSE), gap=1, ...) 

# those <- function(pattern, x, ...) 

# rgb2rgb <- function(xrgb)

# panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
# agrega un panel con la correlación y tamaño del texto proporcional, para usar en pairs

#pairs.3.panels <- function(data, upper.panel=panel.cor, diag.panel=panel.hist, lower.panel=points, pch=20, gap=0, ...) 
# ejecuta un pairs con panels tipicos

#order.by.ranks <- function(data, rank.func=mean, return.order=TRUE, cluster=1)
#ordena los renglones de data dado el rank promedio de las columnas de data considerando rank correlacionados y anticorrelacionados
# Order rows of data analyzing the average ranks of the columns, cluster is used to reverse order for odd clusters

#tolower(x)
#toupper(x)
#capture.output
#dnorm - gives the puntual value of the probability at any given point
#pnorm - gives the cumulative probability until the specified point
#qnorm - gives the point which the cumulative probability is the equal to quantile
#e.g.
#dnorm(0) = .399, mean=0, sd=1
#pnorm(0) = .5
#qnorm(.5) = 0

#estblish the default palette

#string.contain <- function(pattern, x, ...)
#dice que strings x contienen el patron pattern
#si se quiere la posición, usar regrexpr mejor que hace eso

#read.zip <- function(zipfile, internalfile=NA, read.function=read.delim, verbose=TRUE, ...)

#read.z <- function(zfile, internalfile=NA, read.function=read.delim,...)


#perm = function(n, x)
#comb = function(n, x)

#text.non.overlap <- function(x, y, labels, adj=c(0.5, 0.5), cex=1, col=1, vfont=NULL, font=NULL, lines=FALSE, precision=3, marginx=1.1, marginy=1.1, background=NULL, erase=NULL) 
#dibuja texto evitando que se translape

#finites <- function(x) x[is.finite(x)]

#change <- function(x, w, r) { x[w] <- r; x }

#colorize <- function(xcol, f=1, alpha=FALSE) 

#readLines.and.grep <- function(xfile, grepbegin, grepend=rep("",length(grepbegin)), postreadfunc=trim, collapse="", block=10000, maxLines=-1, ...) 

#list to matrix
#l2m <- function(xlist) 


library(lattice)

palette(c("black","red","green3","blue","cyan3","magenta","orange","gray",
"darkred","darkgreen","darkblue","darkcyan","purple4","brown3","pink1","gray50"))

palette(c("black","red","green3","blue","cyan3","magenta","orange","gray",
"darkred","darkgreen","darkblue","darkcyan","purple4","brown3","burlywood4", "cornflowerblue", 
"pink2", "darkolivegreen3","tomato", "goldenrod4", "indianred3", #cornsilk3 instead of pink2 ?
"chocolate4","mediumvioletred","orchid3","salmon3","yellow4","royalblue4","tan","seashell4","gold3"))

uniformize <- function(x) {
	if (is.matrix(x) || is.data.frame(x)) {
		xx <- x
		for (i in 1:ncol(x)) {
			xx[,i] <- uniformize(xx[,i])
		}
		#xx <- apply(x, 2, uniformize)
		#rownames(xx) <- rownames(x)
		xx
	}
	else if (is.list(x)) {
		lapply(x, uniformize)
	}
	else {
		e <- ecdf(x)
		e(x)
	}
}


scale.mad <- function(x, center=TRUE, scale=TRUE) {
	if (is.vector(x)) {
		if (center) x <- x - median(x)
		if (scale) x <- x / mad(x)
	} else {
		x <- apply(x, 2, scale.mad, center, scale)
	}
	return (x)
}

correct <- function(x, na=median, nan=na, inf=TRUE, neg.inf=min, pos.inf=max, null=na, outlier=FALSE, outlier.min=0.05, outlier.max=.95) {
	
	do <- function(x, fx, gx) { z <- fx(x); if (any(z)) { f <- is.finite(x); if(is.function(gx)) x[z] <- gx(x[f]) else x[z] <- gx }; x }
	is.posinf <- function(x) { x==Inf }
	is.neginf <- function(x) { x==-Inf }
	
	y <- x
	y <- do(y, is.null, null)
	y <- do(y, is.nan, nan)
	y <- do(y, is.na, na)
	if (inf) y <- do(y, is.neginf, neg.inf)
	if (inf) y <- do(y, is.posinf, pos.inf)
	if (outlier && !is.null(outlier.min) && is.finite(outlier.min)) {
		q <- quantile(y, outlier.min)
		y[y <= q] <- q
	}
	if (outlier && !is.null(outlier.max) &&is.finite(outlier.max)) {
		q <- quantile(y, outlier.max)
		y[y >= q] <- q
	}
	y
}

similarity <- function(x, y=NULL) {
	if (is.matrix(x)  || is.data.frame(x)) {
		if (is.null(y)) {
			m <- matrix(1, ncol=ncol(x), nrow=ncol(x))			
			colnames(m) <- colnames(x)
			rownames(m) <- colnames(x)
			for (i in 1:(ncol(x)-1)) {
				for (j in (i+1):ncol(x)) {
					m[i,j] <- m[j,i] <- similarity(x[,i], x[,j])
				}
			}
			m
		} else {
			m <- matrix(0, ncol=ncol(x), nrow=ncol(y))
			colnames(m) <- colnames(x)
			rownames(m) <- colnames(y)
			for (i in 1:ncol(x)) {
				for (j in 1:ncol(y)) {
					m[i,j] <- similarity(x[,i], y[,j])
				}
			}
			m
		}
	} else {
		sum(x %in% y) / min(length(x), length(y))
	}
}


redgreenblue <- function(col=-3, nc=ifelse(length(col)==2,16,8), len=NULL, values=NULL) {
	
	if (is.numeric(col)) {
		if (length(col) == 1  && col > 0) { len = col; col = -3;}
		if (all(col==-1)) col <- c("green3","black","red3")
		else if (all(col==-2)) col <- c("blue","yellow")
		else if (all(col==-3)) col <- c("darkblue","white","darkred") 
		else if (all(col==-4)) col <- c("cyan3","white","magenta4") 
		else if (all(col==-5)) col <- c("darkgreen","white","purple3")
		else if (all(col==-6)) col <- c("cyan","black","magenta")
		else if (all(col==-7)) col <- c("cyan3","white","purple")
		else if (all(col==-8)) col <- c("yellow3","white","purple3")
		else if (all(col==-9)) col <- c("darkgreen","white","darkred")
		else if (all(col==-10)) col <- c("cyan","yellow","magenta")
		else if (all(col==-11)) col <- c("white","black")
		else if (all(col==-12)) col <- c("green","red")
		else if (all(col==-13)) col <- c("blue","red")
		else if (all(col==-14)) col <- c("green","blue")
		else if (all(col==-15)) col <- c("cyan","magenta")
		else if (all(col==-16)) col <- c("yellow","purple")
		else if (all(col==-17)) col <- c("cyan","blue")
		else if (all(col==-18)) col <- c("yellow","red")
		else if (all(col==-19)) col <- c("green","purple")
		else if (all(col==-20)) col <- c("cyan","magenta","yellow")
		else if (all(col==-21)) col <- c("red","green","blue")
		else if (all(col==-22)) col <- c("blue","green","red")
		else if (all(col==-23)) col <- c("blue","magenta","green","cyan","red")
		else if (all(col==-24)) col <- c("yellow","black","purple")
		else if (all(col==-25)) col <- c("deepskyblue4","lightgoldenrodyellow","salmon")
		else if (all(col < 0)) col <- c("green","black","red")
	} else {
		if (length(col)==1 && col %in% brewer.pal(n=3,name="x",return.names=TRUE))
			col <- brewer.pal(n=7,name=col)
		else if (length(col)==1 && is.character(col)) {
			if (col=="modern1") col <- c("deepskyblue4","lightgoldenrodyellow","salmon")
			else if (col=="modern2") col <- c("steelblue","wheat","tomato")
			else if (col=="modern3") col <- c("steelblue4","wheat","tomato3")
			else if (col=="modern4") col <- c("lightskyblue4","ivory3","chocolate4")
			else if (col=="modern5") col <- c("lightskyblue4","mistyrose","firebrick4")
			else if (col=="modern6") col <- c("indianred4","lavender","darkolivegreen4")
			else if (col=="modern7") col <- c("plum4","seashell1","royalblue4")
			else if (col=="modern8") col <- c("ivory4","snow2","goldenrod")
			else if (col=="modern9") col <- c("azure3","burlywood4","pink2")
			else if (col=="modern0") col <- c("antiquewhite4","ivory3","maroon")
			else if (col=="modernA") col <- c("coral1","snow2","steelblue3")
			else if (col=="improved") {
				col <- c(rev(redgreenblue("Paired")),redgreenblue(c("#A6CEE3","#FFFFFF"),nc=12))
				if (!is.null(len)) nc <- len
				if (nc > length(col)) col <- col[seq(1,length(col),length.out=nc)]
				else col <- redgreenblue(col, len=nc)
				return (col)
			}
			else col <- c("deepskyblue4","lightgoldenrodyellow","firebrick3")
		}
	}
	x <- col2rgb(col)/255
	x <- cbind(x, c(0,0,0))
	xcol <- list()
	for (i in 1:(length(col)-1)) {
	   	if (length(nc) > 1) {
	   		L <- 255*((nc[i]/sum(nc)) * length(col))
	   	} else {
	   		L <- 255
	   	}
		dec <- L:0/L
	   	cre <- rev(dec)
		xcol[[i]] <- rgb(pmax(0,pmin(1,dec*x[1,i]+cre*x[1,i+1])), 
						 pmax(0,pmin(1,dec*x[2,i]+cre*x[2,i+1])),
						 pmax(0,pmin(1,dec*x[3,i]+cre*x[3,i+1])))
		if (i > 1) 
			xcol[[i]] <- xcol[[i]][-1]
	}
	zcol <- unlist(xcol)
	if (is.null(len))
		len <- if (length(nc) == 1) nc * (length(col)-1) - 1 else sum(nc)
	col <- zcol[seq(1,length(zcol),length.out=len)]
	if (is.null(values)) return (col)
	ev <- ecdf(values)
	return (col[round(ev(values) * (length(col)-1) + 1,0)])

	### this is old versions
	if (!is.null(len))
		nc <- len / (length(col)-1)
   	zcol <- matrix(0, ncol=round(nc*(length(col)-1),0), nrow=3)
   	for (i in 1:(length(col)-1)) {
   		ini <- round((i-1)*nc+1,0)
   		fin <- round(i*nc,0)
   		n <- fin-ini+1
   		z <- if (n==1) 1 else (n-1):0/(n-1)
   		s <- rev(z) #if (n==1) 0 else 0:(n-1)/n
   		zcol[1,ini:fin] <- x[1,i]*z+x[1,i+1]*s
   		zcol[2,ini:fin] <- x[2,i]*z+x[2,i+1]*s
   		zcol[3,ini:fin] <- x[3,i]*z+x[3,i+1]*s
   	}
   	
	if (!is.null(len))
		nc <- len / (length(col)-1)
   	dec <- nc:0/nc #if (nc < 2) 1 else nc:0/nc
   	cre <- rev(dec)
	for (i in 1:(length(col)-1)) {
		xcol[[i]] <- rgb(dec*x[1,i]+cre*x[1,i+1], 
						 dec*x[2,i]+cre*x[2,i+1],
						 dec*x[3,i]+cre*x[3,i+1])
		#xcol[[i]] <- rgb(nc:0/nc*x[1,i]+0:nc/nc*x[1,i+1], 
		#				 nc:0/nc*x[2,i]+0:nc/nc*x[2,i+1],
		#				 nc:0/nc*x[3,i]+0:nc/nc*x[3,i+1])
		if (i > 1) 
			xcol[[i]] <- xcol[[i]][-1]
	}
	col <- unlist(xcol)
	if (!is.null(len)) col <- c(col,col)[1:len]
	if (is.null(values)) return (col)
	ev <- ecdf(values)
	col[round(ev(values) * (length(col)-1) + 1,0)]
}


#overlap between two lists
p.overlap <- function(balls, white, drawn, observed.white, lower.tail=FALSE) {
	black <- balls - white
	phyper(observed.white, white, black, drawn, lower.tail=lower.tail)
}

cumaverage <- function(x) {
	(x[-1] + x[-length(x)]) / 2
}


#Cancer Outlier Profile Analysis - Tomlins - Science - 2005
copa <- function(x, q=NULL) {
  if (is.vector(x)) {
  	  x <- na.omit(x)
  	  md <- mad(x)
  	  if (md == 0) md <- mad(jitter(x))
	  kp <- (x - median(x)) / md
	  if (missing(q) || is.null(q)) return(kp)
	  quantile(kp, q, na.rm=TRUE)
  } else if (is.list(x)) {
	  lapply(x, copa, q)
  } else {
	  kp <- apply(t(x), 2, copa, q)
	  if (!is.null(q)) kp
	  else t(kp)
  }
}


#Outlier Sum - Tibshirani - Biostatistics - 2007 ..... COPA-GENERALIZATION
# group - factor of two levels, level 1 is used as reference and level 2 is used to compute the outlier.sum
outlier.sum <- function(x, group=NULL, ref=levels(group)[1], f=NULL, f2=NULL) {
	if (is.vector(x)) {
		if (is.null(f)) {
			if (!is.null(group)) {
				f <- which(group == ref)
				f2 <- which(group != ref)
			} else {
				f <- f2 <- 1:length(x)
			}
		}
  		md <- mad(x[f], na.rm=TRUE)
  		if (md == 0) md <- mad(jitter(x[f]), na.rm=TRUE)
		qs <- quantile(x[f2], c(0.25, 0.5, 0.75), names=FALSE, na.rm=TRUE)
		kp <- (x[f2] - qs[2]) / md
		iqr <- qs[3] - qs[1]
		w <- max(abs(sum(kp[kp < qs[1] - iqr], na.rm=TRUE)), abs(sum(kp[kp > qs[3] + iqr], na.rm=TRUE)))
		w
	} else if (is.list(x)) {
		lapply(x, outlier.sum, group, ref, f=which(group == ref), f2=which(group != ref))
	} else {
		apply(t(x), 2, outlier.sum, group, ref, f=which(group == ref), f2=which(group != ref))
	}
}


# min.t computes the minimum t-test in sortered values making all possible partitions (minus 2*k + 1 samples)
# k is the minimum number of samples
min.t <- function(x, k=length(x)*.1, parts=max(length(x)-k*2+1,1), index=FALSE, test=t.test, ...) {
	if (is.vector(x)) {
		s <- sort(x)
		l <- length(x)
		parts <- seq(k+1, l-k+1, length.out=parts)
		ps <- sapply(parts, function(i) test(x[1:(i-1)],x[i:l], ...)$p.value)
		if (index) parts[which.min(ps)] else min(ps)
	} else if (is.list(x)) {
		lapply(x, min.t, k, parts, index, test, ...)
	} else {
		t(apply(t(x), 2, min.t, k, parts, index, test, ...))
	}
}


#rlog is monotonical for integer numbers (positive and negative)
#decimal numbers between (-1~0~+1) create asymptote
#previous version used offset=1, this create a discontinuity and could create artifacts
rlog <- function(v, func=log2, offset=2, cutoff=0) {
	suppressWarnings(x <- func(v))
	suppressWarnings(x[v <= cutoff] <- func(1/(abs(v[v <= cutoff])+offset)))
	x
}

rantilog <- function(v, func=function(v) 2^v, offset=2, cutoff=0) {
	suppressWarnings(x <- func(v))
	suppressWarnings(x[v <= cutoff] <- -(func(-v[v <= cutoff]))+offset)
	x
}


glog <- function(y, lambda=1, func=log2) {
	return(func(y+sqrt(y^2+lambda)))
}

gantilog <- function(x, lambda=1, func=function(x) 2^x) {
	z <- func(x)
	#return(sqrt((z^2-lambda)/4)*sign(x))
	return(z/2)
}


apply.df <- function(x, margin, fun, ...) {
	xl <- list()
	if (margin==1) {
		for (i in 1:nrow(x)) xl[[i]] <- fun(x[i,], ...)
		names(xl) <- rownames(x)
	} else {
		for (i in 1:ncol(x)) xl[[i]] <- fun(x[,i], ...)
		names(xl) <- colnames(x)
	}
	if (all(unlist(lapply(xl,length))==length(xl[[1]]))) if (length(xl[[1]]) > 1)xl <- data.frame(xl) else xl <- unlist(xl)
	xl
}



circularize <- function(x, len) {
	if (length(x) < len) {
		rep(x, length.out=len)
	} else {
		x[1:len]
	}
	
}

#FDR, p.adjust in stats package (should be present by default)

#plot densities from columns of a matrix (or variables)
plot.densities <- function(m, m2, m3, m4, lty=rep(1:6, each=length(palette())), lwd=1, cex=1, legends=TRUE, legend.pos="topleft", legend.cols=1, legend.cex=1, horiz=FALSE, xlim=if(horiz) c(mny,mxy) else c(mnx,mxx),ylim=if(!horiz) c(mny,mxy) else c(mnx,mxx), pch=NULL, type="l", npch=20, col=1:length(d), normalize=FALSE, q=NULL, overall=length(d) > 4, ylab=ifelse(normalize,"Relative Density","Density"), xlab="", main="", ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- as.vector(m)
		if (!missing(m2)) { m <- cbind(m,as.vector(m2)); xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m <- cbind(m,as.vector(m3)); xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m <- cbind(m,as.vector(m4)); xn <- c(xn,deparse(substitute(m4))) }
		colnames(m) <- xn
	}
	else if (!is.list(m)){
		m <- as.matrix(m)
	}
	if (is.list(m)) d <- lapply(m, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
	else d <- apply(m, 2, function(x) { suppressWarnings(density(x[is.finite(x)],...)) })
	lty <- circularize(lty, length(d))
	lwd <- circularize(lwd, length(d))
	cex <- circularize(cex, length(d))
	if (!is.null(q)) {
		if (is.list(m)) qq <- lapply(m, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
		else qq <- apply(m, 2, function(x) { suppressWarnings(quantile(x[is.finite(x)], q)) })
	}
	if (overall) {
		o <- unlist(m)
		n <- names(d)
		d[[length(d)+1]] <- density(o[is.finite(o)], ...)
		names(d) <- c(n, "OVERALL")
		lty[length(d)] <- 1
		lwd[length(d)] <- 3
		col[length(d)] <- 1
		cex[length(d)] <- 1
		if (!is.null(q)) qq[length(d)] <- suppressWarnings(quantile(o[is.finite(o)], q))
	}
	if (normalize) {
		for (i in 1:length(d)) {
			d[[i]]$y <- d[[i]]$y/max(d[[i]]$y)
		}
	}
	mxy <- max(unlist(lapply(d, function(x) max(x$y))))
	mny <- min(unlist(lapply(d, function(x) min(x$y))))
	mxx <- max(unlist(lapply(d, function(x) max(x$x))))
	mnx <- min(unlist(lapply(d, function(x) min(x$x))))
	lty <- rep(lty, length(d), length.out=length(d))
	pchs <- if (type=="l" || is.null(pch)) rep(0,length(d)) else rep(pch, length(d), length.out=length(d))
	suppressWarnings(plot(0,0,type="n",ylim=ylim,xlim=xlim, main=main, xlab=xlab, ylab=ylab, ...))
	cols <- if (horiz) c("y","x") else c("x","y")
	for (i in length(d):1) {
		suppressWarnings(lines(d[[i]][cols], col=col[i], lty=lty[i], type="l", lwd=lwd[i], cex=cex[i]))
		if (!is.null(q)) {
			for (k in 1:length(q)) {
				if (horiz) abline(h=qq[[i]][k], lty=k, col=col[i])
				else abline(v=qq[[i]][k], lty=k, col=col[i])
			}
		}
	}
	if (type != "l") {
		pchs <- rep(pchs, length(d), length.out=length(d))
		if (horiz) for (i in 1:length(d)) suppressWarnings(points(d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], d[[i]]$x[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pchs[i], ...))
		else	   for (i in 1:length(d)) suppressWarnings(points(d[[i]]$x[seq(1,length(d[[i]]$y), length.out=npch)], d[[i]]$y[seq(1,length(d[[i]]$y), length.out=npch)], col=col[i], pch=pchs[i], ...))
	}
	if (length(legends) > 1 || (!is.character(legends) && legends)) {
		if ((length(legends) == length(d)) && (length(legends) > 0)) xl <- legends
		else {
			xl <- if (is.list(m)) names(m) else colnames(m)
			if (overall && !is.null(xl)) xl <- c(xl, "OVERALL")
		}
		if (!is.null(xl)) {
			plot.legend(legend.pos, xl, col=col[1:length(d)], lty=lty, ncol=legend.cols, cex=legend.cex, pch=if (is.null(pch)) NULL else pchs)
		}
	}
	invisible(d)
}


plot.ma <- function(m, m2, m3, m4, gap=0, pch=20, xaxt="n",yaxt="n", normalise=FALSE, ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- as.vector(m)
		if (!missing(m2)) { m <- cbind(m,as.vector(m2)); xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m <- cbind(m,as.vector(m3)); xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m <- cbind(m,as.vector(m4)); xn <- c(xn,deparse(substitute(m4))) }
		colnames(m) <- xn
	}
	else {
		m <- as.matrix(m)
	}
	xc <- c()
	for (i in 1:ncol(m)) for (j in 1:ncol(m)) if (i != j) xc <- c(xc, i, j)
	xc.i <- 0
	mapanel <- function(x,y,col,...) {
		 cnm <- colnames(m)
		 if (is.null(cnm)) cnm <- paste("var",1:ncol(m))
		 usr <- par("usr"); on.exit(par(usr))
		 a <- (x+y)/2
		 m <- (x-y)
		 text(mean(par("usr")[1]),par("usr")[4],paste(cnm[xc[xc.i+1]],"-",cnm[xc[xc.i+2]]),adj=c(0,1))
		 xc.i <<- xc.i + 2
		 if (normalise) {
			xl <- loess.smooth(a,m)
			m <- unlist(lapply(2:length(xl$x)-1, function(i) m[a > xl$x[i] & a <= xl$x[i+1]] - mean(xl$y[c(i,i+1)])))
			a <- unlist(lapply(2:length(xl$x)-1, function(i) a[a > xl$x[i] & a <= xl$x[i+1]]))
		 }
		 par(usr = c(min(a)*0.95, max(a)*1.05, min(m)*0.95, max(m)*1.05))
		 abline(h=(2:-2), col=c("red","darkred","gray","darkgreen","green"), lty=3)
		 mf <- pmax(pmin(m,2),-2)
		 points(a,m,col=if(missing(col)) rgb(pmax(mf/2,0),pmax(-mf/2,0),0) else col,...)
		 lines(loess.smooth(a,m), col="blue", lwd=2, lty=1)
		 abline(h=0, col="gray")

	}
	pairs(m, upper.panel=mapanel, gap=gap, lower.panel=mapanel, pch=pch, xaxt=xaxt,yaxt=yaxt, ...)
}



#plot comparative histograms from columns of a matrix (or variables)
plot.hist <- function(m, m2, m3, m4, 
	col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), 
	separate=FALSE, beside=TRUE, lines=FALSE, 
	legends=TRUE, legend.cex=1, legend.pos="topleft", legend.cols=1, legend.n=FALSE,
	log="", hollow=FALSE, border=col, 
	labels=FALSE, srt=90, cex=1, 
	ylim=c(0,(if(beside) max(bp) else max(apply(bp,1,sum)))*y.adj), 
	density=FALSE, round.mids=NULL, normalize=NULL, 
	test=FALSE, test.func=dist.chisq.test, round.test=8, 
	y.adj=1.25, names.arg=NULL, total=FALSE, cumulative=FALSE, roundfunc=round, 
	quantiles=NULL, round.quantiles=4,
	...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- list(m)
		if (!missing(m2)) { m[[length(m)+1]] <- m2; xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m[[length(m)+1]] <- m3; xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m[[length(m)+1]] <- m4; xn <- c(xn,deparse(substitute(m4))) }
		names(m) <- xn
		if (missing(col)) col <- 1:length(m)
	} else if (is.data.frame(m) || is.matrix(m)) {
		nm <- colnames(m)
		m <- lapply(1:ncol(m), function(i) m[,i])
		if (!is.null(nm)) names(m) <- nm
		if (missing(col)) col <- 1:length(m)
	} else if (is.vector(m) && !is.list(m)) {
		m <- list(m)
	}
	n <- unlist(lapply(m, function(x) sum(is.finite(x))))
	x <- list(...)
	h  <- suppressWarnings(hist.default(unlist(m), plot=FALSE, ...))
	b  <- h$breaks
	if (any(names(x)=="breaks")) x[names(x)=="breaks"] <- NULL
	xl <- unlist(lapply(m, length))
	hs <- suppressWarnings(lapply(m, hist, breaks=b, plot=F, ...=x))
	bp <- data.matrix(data.frame(lapply(hs, function(x) x$counts)))
	obp <- bp
	if (log=="y") bp <- log10(bp)
	bp[!is.finite(bp)] <- 0
	if (density) if (any(xl==xl[1])) bp <- t(t(bp)/apply(bp,2,sum))
	if (!is.null(normalize)) if (length(normalize) == 1) bp <- bp/ifelse(bp[,normalize]==0,1,bp[,normalize]) else if (length(normalize)==nrow(bp)) bp <- bp/normalize
	
	suppressWarnings(barplot(if(separate) bp else t(bp), beside=beside, names.arg=if (!is.null(names.arg)) names.arg else if (is.null(round.mids)) h$mids else if (round.mids < 0 & max(h$mids) < 10000) paste(roundfunc(b[-length(b)],-round.mids),"~",roundfunc(b[-1],-round.mids),sep="") else roundfunc(h$mids,round.mids), col=if(hollow) 0 else col, border=border, ylim=ylim, ...))

	if (lines & beside) {
		for (i in 1:ncol(bp)) {
			lines(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,bp[,i],col=col[i],lwd=6-ncol(bp))
		}
	}
	if (cumulative) {
		for (i in 1:ncol(bp)) {
			lines(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,max(bp)*cumsum(bp[,i])/sum(bp[,i]),col=col[i],lwd=1, lty=2)
		}
		axis(side=4, at=axTicks(2), labels=paste(round(100*axTicks(2)/max(bp),0),"%",sep=""))
	}
	if (labels & beside) {
		for (i in 1:ncol(bp)) {
			suppressWarnings(text(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,bp[,i],obp[,i],col=col[i],adj=c(0,0.5),srt=srt,cex=cex,...))
		}
	}
	leg <- FALSE
	if (length(legends) > 1 || (!is.character(legends) && legends) || (!is.null(quantiles) && any(quantiles))) {
		if (!(length(legends) > 1 || (!is.character(legends) && legends))) {
			legends <- 1:length(m)
		}
		if (!is.null(quantiles) || any(quantiles)) {
			if (is.logical(quantiles) && any(quantiles)) {
				quantiles <- c(0, .05, .5, .95, 1)
			}
			qs <- unlist(sapply(1:length(m), function(i) {
					x <- m[[i]]
					q <- quantile(x,quantiles,na.rm=TRUE)
					#rug(((q-min(h$mids))/(max(h$mids)-min(h$mids)))*length(h$mids), col=col[i])
					#rug(((q-min(b))/(max(b)-min(b)))*(par("usr")[2]-par("usr")[1])+par("usr")[1], col=col[i],line=0.5)
					paste(paste(quantiles,"=",format(round(q,round.quantiles)),sep=""),collapse=", ")
					}
				))			
		} else {
			qs <- rep("",length(m))
		}
		legtext <- if(length(legends) > 1) legends else if (!total) names(m) else paste(names(m), " (", lapply(m, function(x) sum(!is.na(x))), ")", sep="")
		if (any(nchar(qs) > 1)) 
			legtext <- paste(legtext, "(Quantiles: ",qs,")")
		if (legend.n) { 
			legtext <- paste0(legtext, " (n=",n,")")
		}
		plot.legend(legend.pos[1], legtext, fill=col, ncol=legend.cols, cex=legend.cex)
		leg <- TRUE
	}
	if (test && is.function(test.func)) {
		ok <- 1:ncol(bp)
		p <- matrix(NA, ncol=length(ok), nrow=length(ok))
		for (i in 1:length(ok)) {
			for (j in 1:length(ok)) {
				suppressWarnings(p[i,j] <- test.func(bp[,ok[i]], bp[,ok[j]])$p.value)
			}
		}
		p <- round(p, round.test)
		if (ncol(bp) == 2) {
			text(legend.pos[if (leg  && legend.pos > 1) 2 else 1], p[1,2])
		} else {
			mode(p) <- "character"
			for (i in 1:ncol(p)) for (j in i:ncol(p)) p[i,j] <- "-"
			nombres <- (if(length(legends) > 1) legends else names(m))[ok]
			p <- rbind(nombres,p)
			p <- cbind(c("",if(!leg && total) paste(nombres," (",lapply(m, function(x) sum(!is.na(x))),")",sep="") else nombres),p)
			z <- matrix(0, ncol=length(ok)+1, nrow=length(ok)+1)
			xcol <- z
			xcol[2:nrow(xcol),1] <- xcol[1,2:nrow(xcol)] <- col[ok]
			xfill <- z
			xfill[2:nrow(xfill),1] <- xfill[1,2:nrow(xfill)] <- col[ok]
			plot.legend(legend.pos[if (leg  && legend.pos > 1) 2 else 1], p[,-ncol(p)], ncol=length(ok), col=xcol[,-ncol(xcol)], fill=xfill[,-ncol(xfill)], border.col=xcol[,-ncol(xcol)], cex=legend.cex)
		}
	}

	invisible(bp)
}



histq <- function(x, q=c(.5, .25, .75, .05, .95, 0, 1), 
	col=c(1,rep(2:10,times=rep(2,9)))+1, 
	lty=c(1,rep(2:10,times=rep(2,9))), 
	cex=.75, main=paste(deparse(substitute(x), 500), collapse = "\n"), 
	...) {
	a <- hist(x, main=main, ...)
	iqs <- quantile(x, q, na.rm=TRUE)
	for (i in 1:length(q)) {
		iq <- iqs[i]
		abline(v=iq, col=col[i], lty=lty[i])
		text(iq, max(a$counts), format(iq), col=col[i], srt=90, adj=c(1.05,1.05), cex=cex)
		text(iq, par("usr")[3], format(q[i]), col=col[i], srt=90, adj=c(0,1), cex=cex)
	}
}

histm <- function(x, m=c(mean(x), mean(x)+sd(x), mean(x)-sd(x), mean(x)+sd(x)*2, mean(x)-sd(x)*2), col=c(1,rep(2:10,times=rep(2,9)))+1, lty=c(1,rep(2:10,times=rep(2,9))), cex=.75, ...) {
	a <- hist(x, ...)
	for (i in 1:length(m)) {
		iq <- m[i]
		abline(v=iq, col=col[i], lty=lty[i])
		text(iq, max(a$counts), format(iq), col=col[i], srt=90, adj=c(1.05,1.05), cex=cex)
	}
}



#plot comparative histograms in hollow mode from columns of a matrix (or variables)
plot.hist2 <- function(m, m2, m3, m4, las=0, 
	base=if (missing(cumulative) && missing(shading)) 1 else 0, 
	col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), 
	lty=rep(1:6, each=length(palette())), xaxt="s", yaxt="s", xlab="", 
	ylab=if(!is.null(normalize) && any(normalize)) "Normalized Frequency" else if(density) "Density" else "Frequency", 
	separate=FALSE, 
	legends=TRUE, legend.cex=1, legend.pos=c("topleft","topright"), legend.cols=1, log="", legend.n=FALSE,
	hollow=FALSE, border=col, labels=FALSE, 
	srt=90, cex=1, ylim=c(0,max(bp)*y.adj), 
	density=FALSE, splice=0, shading=-1, 
	round.mids=NULL, normalize=NULL, 
	test=FALSE, test.func=dist.chisq.test, round.test=8, y.adj=ifelse((is.character(labels) && any(tolower(labels) == "top")),1.5,1.25), angle=45, 
	total=FALSE, cumulative=FALSE, roundfunc=round, 
	quantiles=NULL, round.quantiles=4, ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- list(m)
		if (!missing(m2)) { m[[length(m)+1]] <- m2; xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m[[length(m)+1]] <- m3; xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m[[length(m)+1]] <- m4; xn <- c(xn,deparse(substitute(m4))) }
		names(m) <- xn
		if (missing(col)) col <- 1:length(m)
	} else if (is.data.frame(m) || is.matrix(m)) {
		nm <- colnames(m)
		m <- lapply(1:ncol(m), function(i) m[,i])
		if (!is.null(nm)) names(m) <- nm
		if (missing(col)) col <- 1:length(m)
	} else if (is.vector(m) && !is.list(m)) {
		m <- list(m)
	}
	n <- unlist(lapply(m, function(x) sum(is.finite(x))))
	x <- list(...)
	h  <- suppressWarnings(hist.default(unlist(m), plot=FALSE, ...))
	b  <- h$breaks
	if (any(names(x)=="breaks")) x[names(x)=="breaks"] <- NULL
	xl <- unlist(lapply(m, length))
	hs <- suppressWarnings(lapply(m, hist, breaks=b, plot=F, ...=x))
	bp <- data.matrix(data.frame(lapply(hs, function(x) x$counts)))
	obp <- bp
	antilog <- function(x) x
	if (log=="y") { bp <- log10(bp); antilog <- function(x) ifelse(is.finite(x),10^x,0); }
	bpx <- bp
	bp[!is.finite(bp)] <- 0
	if (density) if (any(xl==xl[1])) bp <- t(t(bp)/apply(bp,2,sum))
	if (!is.null(normalize)) if (length(normalize) == 1) { if (is.numeric(normalize)) bp <- bp/ifelse(bp[,normalize]==0,1,bp[,normalize]) else if (normalize) bp <- t(t(bp)/apply(bp,2,max)) } else if (length(normalize)==nrow(bp)) bp <- bp/normalize
	suppressWarnings(plot(0,0,type="n",ylim=ylim, xlim=c(0,length(b)-1), xaxt="n", xlab=xlab, ylab=ylab, las=las, ...))
	if (is.null(xaxt) || xaxt != "n") suppressWarnings(axis(1,at=1:length(b)-1,labels=if (is.null(round.mids)) pretty(b, n=length(b)) else roundfunc(b,round.mids),las=las,...)) 
	if (is.null(yaxt) || yaxt != "n") suppressWarnings(axis(2, las=las, ...))
	#if (cumulative) abline(v=axTicks(1), col="#DDDDDD")
	b.m <- (b[-length(b)]+b[-1])/2
	xp <- rep(0:(length(b)-1),each=2)
	#xp <- xp[-c(1,length(b))]
	ok <- unique(c(base,1:ncol(bp)))
	ok <- ok[ok > 0 & ok <= ncol(bp)]
	lty <- rep(lty, length(ok), length.out=length(ok))
	angle <- rep(angle, length.out=length(ok))
	if (length(shading) > 1) shading <- rep(shading, length.out=length(ok))
	leg.dens <- rep(NA, ncol(bp))
	for (i in ok) {
		yp <- c(0,rep(bp[,i],each=2),0)
		leg.dens[i] <- if(length(shading)==ncol(bp)) shading[i] else shading*(length(ok)-i+1)
		suppressWarnings(polygon(c(xp,0)+splice*i,c(yp,0),col=ifelse(i %in% base,col[i],ifelse(leg.dens[i] > 0, col[i],0)),border=col[i], density=leg.dens[i], lty=lty[i], angle=angle[i], ...))
		if ((is.logical(labels) && labels[min(i,length(labels))]) || (i %in% labels)) suppressWarnings(text(x=1:nrow(bp)-0.9*i/ncol(bp), y=bp[,i], labels=round(antilog(bpx[,i]),0), col=col[i], srt=srt, cex=cex, adj=c(0,0.5), ...))
		if ((is.character(labels) && labels[min(i,length(labels))] == "%")) suppressWarnings(text(x=1:nrow(bp)-0.5, y=bp[,i], labels=paste(round(bpx[,i]*100/sum(bpx[,i]),1),"%",sep=""), col=col[i], srt=srt, cex=cex, adj=c(0,0.5), ...))
		if ((is.character(labels) && tolower(labels[min(i,length(labels))]) == "top")) suppressWarnings(text(x=1:nrow(bp)-0.5, y=par("usr")[4]-i*strheight("X",cex=cex)*1.1, labels=round(antilog(bpx[,i]),0), col=col[i], cex=cex, adj=c(0.5,0), ...))
	}
	if (length(base) > 1 || any(shading != -1) || all(shading[base] == 0))
		for (i in ok) {
			yp <- c(0,rep(bp[,i],each=2),0)
			suppressWarnings(polygon(c(xp,0)+splice*i,c(yp,0),col=0,border=col[i], lty=lty[i], ...))
		}
	if (log=="y") axis(4,at=axTicks(2),labels=round(antilog(axTicks(2)),0),las=las)
	abline(h=0, col="grey")

	if (cumulative) {
		mx <- max(max(bp),par("usr")[4]*.95,rev(axTicks(2))[1])
		for (i in 1:ncol(bp)) {
			#lines(1:nrow(bp)-0.5,mx*cumsum(bp[,i])/sum(bp[,i]),col=col[i],lwd=1, lty=2)
			lines(0:nrow(bp),mx*(ecdf(m[[i]])(b)),col=col[i],lwd=1, lty=2)
		}
		#axis(side=4, at=axTicks(2), labels=paste(round(100*axTicks(2)/max(bp),0),"%",sep=""))
		xp <- if (is.numeric(cumulative)) cumulative else 8 #rev(which.min(abs(1:5*2-length(axTicks(2)))))[1]*2
		xp <- 0:xp/xp
		axis(side=4, at=mx*xp, labels=paste(round(100*xp,0),"%",sep=""))
		abline(h=mx*xp, col="#DDDDDD")
	}
	
	#suppressWarnings(barplot(if(separate) bp else t(bp), beside=beside, names.arg=h$mids, col=if(hollow) 0 else col, border=border, ylim=ylim, ...))


	#if (length(legends) > 1 || (!is.character(legends) && legends)) 
	#	plot.legend(legend.pos[1],(if(length(legends) > 1) legends else names(m)), col=col, border.col=col, fill=sapply(1:length(ok),function(i) ifelse(i %in% base,col[i],ifelse(leg.dens[i] > 0,col[i],-1))), density=leg.dens, angle=angle, cex=legend.cex, ncol=legend.cols)
	if (length(legends) > 1 || (!is.character(legends) && legends) || (!is.null(quantiles) && any(quantiles))) {
		if (!(length(legends) > 1 || (!is.character(legends) && legends))) {
			legends <- 1:length(m)
		}
		if (!is.null(quantiles) || any(quantiles)) {
			if (is.logical(quantiles) && any(quantiles)) {
				quantiles <- c(0, .05, .5, .95, 1)
			}
			qs <- unlist(sapply(1:length(m), function(i) {
					x <- m[[i]]
					q <- quantile(x,quantiles,na.rm=TRUE)
					rug(((q-min(b))/(max(b)-min(b)))*(length(b)-1), col=col[i])
					#rug(((q-min(h$mids))/(max(h$mids)-min(h$mids)))*length(h$mids), col=col[i])
					paste(paste(quantiles,"=",format(round(q,round.quantiles)),sep=""),collapse=", ")
					}
				))
		} else {
			qs <- rep("",length(m))
		}
		legtext <- if(length(legends) > 1) legends else if (!total) names(m) else paste(names(m), " (", lapply(m, function(x) sum(!is.na(x))), ")", sep="")
		if (any(nchar(qs) > 1)) 
			legtext <- paste(legtext, "(Quantiles: ",qs,")")
		if (legend.n) { 
			legtext <- paste0(legtext, " (n=",n,")")
		}
		plot.legend(legend.pos[1], legtext, fill=col, ncol=legend.cols, cex=legend.cex)
		leg <- TRUE
	}

	if (test && is.function(test.func)) {
		p <- matrix(NA, ncol=length(ok), nrow=length(ok))
		for (i in 1:length(ok)) {
			for (j in 1:length(ok)) {
				suppressWarnings(p[i,j] <- test.func(obp[,ok[i]], obp[,ok[j]])$p.value)
			}
		}
		if (!is.null(round.test)) p <- round(p, round.test)
		if ((length(legends) > 1 || (!is.character(legends) && legends))  && length(ok)==2) {
			plot.legend(legend.pos[2], paste("p.value =",p[2,1]), cex=legend.cex, ncol=legend.cols)
		} else {
			mode(p) <- "character"
			for (i in 1:ncol(p)) for (j in i:ncol(p)) p[i,j] <- "-"
			nombres <- (if(length(legends) > 1) legends else names(m))[ok]
			p <- rbind(nombres,p)
			p <- cbind(c("p.values",nombres),p)
			p <- cbind(c("n",paste(apply(obp, 2, sum)[ok],sep="")), p)
			z <- matrix(0, ncol=length(ok)+2, nrow=length(ok)+1)
			xcol <- z
			xcol[2:nrow(xcol),2] <- xcol[1,3:ncol(xcol)] <- col[ok]
			xfill <- z
			xfill[2:nrow(xfill),2] <- xfill[1,3:ncol(xfill)] <- sapply(ok,function(i) ifelse(i %in% base,col[i],-1))
			plot.legend(legend.pos[2], p[,-ncol(p)], ncol=length(ok)+1, col=xcol[,-ncol(xcol)], fill=xfill[,-ncol(xfill)], border.col=xcol[,-ncol(xcol)], cex=legend.cex)
		}
	}

	invisible(bp)
}

plot.dist <- plot.hist2


#plot comparative histograms from columns of a matrix (or variables)
plot.hist3D <- function(m, m2, m3, m4, col=colorize(if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), alpha=0.666), 
	separate=FALSE, beside=TRUE, lines=FALSE, 
	legends=FALSE, legend.cex=1, legend.pos="topleft", legend.cols=1, legend.n=FALSE,
	log="", hollow=FALSE, border=col, labels=FALSE, srt=90, cex=1, ylim=c(0,(if(beside) max(bp) else max(apply(bp,1,sum)))*y.adj), 
	density=FALSE, round.mids=NULL, normalize=NULL, test=FALSE, test.func=dist.chisq.test, 
	round.test=8, y.adj=1, ysep=.1, xsep=0.5, bw=0.9, las=1, names.arg=TRUE, xaxis.adj=if (length(names.arg) && names.arg[1]) 4 else 0, 
	main="", xlab="", ylab="Frequency", sub="", roundfunc=round, bar.density=NA, bar.angle=45, 
	guidelines=TRUE, ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- list(m)
		if (!missing(m2)) { m[[length(m)+1]] <- m2; xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m[[length(m)+1]] <- m3; xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m[[length(m)+1]] <- m4; xn <- c(xn,deparse(substitute(m4))) }
		names(m) <- xn
		if (missing(col)) col <- 1:length(m)
	} else if (is.data.frame(m) || is.matrix(m)) {
		nm <- colnames(m)
		m <- lapply(1:ncol(m), function(i) m[,i])
		if (!is.null(nm)) names(m) <- nm
		if (missing(col)) col <- 1:length(m)
	} else if (is.vector(m) && !is.list(m)) {
		m <- list(m)
	}
	n <- unlist(lapply(m, function(x) sum(is.finite(x))))
	x <- list(...)
	h  <- suppressWarnings(hist.default(unlist(m), plot=FALSE, ...))
	b  <- h$breaks
	if (any(names(x)=="breaks")) x[names(x)=="breaks"] <- NULL
	xl <- unlist(lapply(m, length))
	hs <- suppressWarnings(lapply(m, hist, breaks=b, plot=F, ...=x))
	bp <- data.matrix(data.frame(lapply(hs, function(x) x$counts)))
	obp <- bp
	if (log=="y") bp <- log10(bp)
	bp[!is.finite(bp)] <- 0
	if (density) if (any(xl==xl[1])) bp <- t(t(bp)/apply(bp,2,sum))
	if (!is.null(normalize)) if (length(normalize) == 1) bp <- bp/bp[,normalize] else if (length(normalize)==nrow(bp)) bp <- bp/normalize

	plotbarheights <- function(h, xadj=0, yadj=0, fg, bg, w=0.5, bd=0, ba=0, ...) {
		#symbols(1:length(h)+xadj, rep(0,length(h))+yadj, rectangles=data.matrix(data.frame(width=rep(w, length(h)), height=h)), fg=fg, bg=bg, add=TRUE, inches=FALSE, adj=0, ...)
		for (i in 1:length(h)) {
			rect(i+xadj, yadj, i+xadj+w, yadj+h[i], col=fg, border=bg, density=bd, angle=ba, ...)
		}
	}
	x.adjust <- function(h, xadj=0, ...) {
		(1:length(h))+xadj
	}
	y.adjust <- function(h, yadj=0, ...) {
		h+yadj
	}
	plot.new()
	plot.window(xlim=c(1,round(nrow(bp)+ncol(bp)*xsep+0.4999+xaxis.adj,0)), ylim=c(ylim[1], ylim[2]*(1+ncol(bp)*ysep)))
	axis(1, at=1:nrow(bp)+0.5, labels=if (is.null(round.mids)) h$mids else if (round.mids < 0 & max(h$mids) < 10000) paste(roundfunc(b[-length(b)],-round.mids),"~",roundfunc(b[-1],-round.mids),sep="") else roundfunc(h$mids,round.mids), las=las)
	tx <- axTicks(side=1)
	xt <- pretty(c(ylim[1], ylim[2]), n=5)
	#axis(2, at=xt, pos=c((ncol(bp)-1)*xsep, (ncol(bp)-1)*ysep*ylim[2]), labels=rep("", length(xt)))
	for (i in 1:length(xt)) lines(x=c(par("usr")[1], (ncol(bp)-1)*xsep), y=xt[i]+c(0, (ncol(bp)-1)*ysep*ylim[2]), col=8, lty=2)
	for (i in 1:length(xt)) lines(x=c((ncol(bp)-1)*xsep, nrow(bp)+ncol(bp)*xsep), y=rep(xt[i]+(ncol(bp)-1)*ysep*ylim[2],2), col=8, lty=2)
	if (guidelines) for (i in 1+0:nrow(bp)) lines(x=i+c(0,(ncol(bp)-1)*xsep), y=c(0,(ncol(bp)-1)*ysep*ylim[2]), col=ifelse(is.logical(guidelines),8,guidelines), lty=2)
	axis(2, at=xt)
	col <- rep(col,len=ncol(bp))
	border <- rep(border,len=ncol(bp))
	for (i in ncol(bp):1) plotbarheights(bp[,i], xadj=(i-1)*xsep, yadj=(i-1)*ysep*ylim[2], fg=ifelse(hollow,0,col[i]), bg=border[i], w=bw, bd=bar.density, ba=bar.angle, ...)
	title(main = main, sub = sub, xlab = xlab, ylab = ylab)

	#suppressWarnings(barplot(if(separate) bp else t(bp), beside=beside, names.arg=if (is.null(round.mids)) h$mids else if (round.mids < 0 & max(h$mids) < 10000) paste(round(b[-length(b)],-round.mids),"~",round(b[-1],-round.mids),sep="") else round(h$mids,round.mids), col=if(hollow) 0 else col, border=border, ylim=ylim, ...))
	if (lines & beside) {
		for (i in 1:ncol(bp)) {
			#lines(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,bp[,i],col=col[i],lwd=6-ncol(bp))
			lines((1:nrow(bp)+(i-1)*xsep+0.5),(bp[,i]+(i-1)*ysep*ylim[2]),col=col[i],lwd=6-ncol(bp))
		}
	}
	if (labels & beside) {
		for (i in 1:ncol(bp)) {
			#suppressWarnings(text(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,bp[,i],obp[,i],col=col[i],adj=c(0,0.5),srt=srt,cex=cex,...))
			xw <- which(obp[,i] != 0)
			suppressWarnings(text((1:nrow(bp)+(i-1)*xsep+0.5)[xw],(bp[,i]+(i-1)*ysep*ylim[2])[xw],obp[xw,i],col=col[i],adj=c(0,0.5),srt=srt,cex=cex,...))
		}
	}
	if (is.character(names.arg) || length(names.arg) > 1 || names.arg) {
		nm <- if(is.character(names.arg) || length(names.arg) > 1) names.arg else names(m)
		text(nrow(bp)+1+(1:ncol(bp))*xsep, (1:ncol(bp)-1)*ysep*ylim[2], nm, col=col, adj=c(0,0.5), cex=cex, ...)
	}
	if (length(legends) > 1  || (!is.character(legends) && legends)) plot.legend(legend.pos, rev(paste(if(length(legends) > 1) legends else names(m), if (legend.n) paste0("(n=",n,")") else "")), fill=rev(col), bg="white", ncol=legend.cols, cex=legend.cex)
	if (test && is.function(test.func)) {
		ok <- 1:ncol(bp)
		p <- matrix(NA, ncol=length(ok), nrow=length(ok))
		for (i in 1:length(ok)) {
			for (j in 1:length(ok)) {
				suppressWarnings(p[i,j] <- test.func(bp[,ok[i]], bp[,ok[j]])$p.value)
			}
		}
		p <- round(p, round.test)
		mode(p) <- "character"
		for (i in 1:ncol(p)) for (j in i:ncol(p)) p[i,j] <- "-"
		nombres <- (if(length(legends) > 1) legends else names(m))[ok]
		p <- rbind(nombres,p)
		p <- cbind(c("",nombres),p)
		z <- matrix(0, ncol=length(ok)+1, nrow=length(ok)+1)
		xcol <- z
		xcol[2:nrow(xcol),1] <- xcol[1,2:nrow(xcol)] <- col[ok]
		xfill <- z
		xfill[2:nrow(xfill),1] <- xfill[1,2:nrow(xfill)] <- col[ok]
		plot.legend(legend.pos, rev(p[,-ncol(p)]), ncol=length(ok), col=rev(xcol[,-ncol(xcol)]), fill=rev(xfill[,-ncol(xfill)]), border.col=rev(xcol[,-ncol(xcol)]), cex=legend.cex, bg="white")
	}

	invisible(bp)
}



plot.bubbles <- function(m, m2, m3, m4, col=colorize(if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), alpha=0.666), 
	border=col, labels=FALSE, srt=90, cex=1, log="", density=TRUE,
	round.mids=NULL, normalize=NULL, las=1, names.arg=TRUE, hollow=FALSE,
	main="", xlab="", ylab="", sub="", roundfunc=round,
	guidelines=TRUE, margin=10, factor=1, ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- list(m)
		if (!missing(m2)) { m[[length(m)+1]] <- m2; xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m[[length(m)+1]] <- m3; xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m[[length(m)+1]] <- m4; xn <- c(xn,deparse(substitute(m4))) }
		names(m) <- xn
		if (missing(col)) col <- 1:length(m)
	} else if (is.data.frame(m) || is.matrix(m)) {
		nm <- colnames(m)
		m <- lapply(1:ncol(m), function(i) m[,i])
		if (!is.null(nm)) names(m) <- nm
		if (missing(col)) col <- 1:length(m)
	} else if (is.vector(m) && !is.list(m)) {
		m <- list(m)
	}
	n <- unlist(lapply(m, function(x) sum(is.finite(x))))
	x <- list(...)
	h  <- suppressWarnings(hist.default(unlist(m), plot=FALSE, ...))
	b  <- h$breaks
	if (any(names(x)=="breaks")) x[names(x)=="breaks"] <- NULL
	xl <- unlist(lapply(m, length))
	hs <- suppressWarnings(lapply(m, hist, breaks=b, plot=F, ...=x))
	bp <- data.matrix(data.frame(lapply(hs, function(x) x$counts)))
	obp <- bp
	if (log=="y") bp <- log10(bp)
	bp[!is.finite(bp)] <- 0
	if (density) if (any(xl==xl[1])) bp <- t(t(bp)/apply(bp,2,sum))
	if (!is.null(normalize)) if (length(normalize) == 1) bp <- bp/bp[,normalize] else if (length(normalize)==nrow(bp)) bp <- bp/normalize

	bp <- apply(bp, 2, function(x) sqrt(x/pi))
	bp <- bp/(2*max(bp,na.rm=TRUE)) # scale to max
	bp <- bp*factor

	plotbarheights <- function(h, xadj=0, yadj=0, fg, bg, w=0.5, bd=0, ba=0, ...) {
		#symbols(1:length(h)+xadj, rep(0,length(h))+yadj, rectangles=data.matrix(data.frame(width=rep(w, length(h)), height=h)), fg=fg, bg=bg, add=TRUE, inches=FALSE, adj=0, ...)
		#for (i in 1:length(h)) {
		#	rect(i+xadj, yadj, i+xadj+w, yadj+h[i], col=fg, border=bg, density=bd, angle=ba, ...)
		#}
		#print(h)
		symbols(1:length(h)+xadj,rep(yadj,length(h)),circles=h+1e-20/2,fg=fg,bg=bg,density=bd, add=TRUE, inches=FALSE)
	}
	par(mar=c(par("mar")[1:3],margin))
	plot.new()
	plot.window(xlim=c(1,nrow(bp)+0.5), ylim=c(1,ncol(bp)))
	axis(1, at=1:nrow(bp)+0.5, labels=if (is.null(round.mids)) h$mids else if (round.mids < 0 & max(h$mids) < 10000) paste(roundfunc(b[-length(b)],-round.mids),"~",roundfunc(b[-1],-round.mids),sep="") else roundfunc(h$mids,round.mids), las=las)
	if (guidelines) for (i in 1:nrow(bp)) lines(x=c(i,i)+0.5, y=c(1,ncol(bp)), col=ifelse(is.logical(guidelines),8,guidelines), lty=1, lwd=0.25)
	if (guidelines) for (i in 1:ncol(bp)) lines(x=c(1,nrow(bp))+0.5, y=c(i,i), col=ifelse(is.logical(guidelines),8,guidelines), lty=1, lwd=0.25)
	tx <- axTicks(side=1)
	xt <- pretty(par("usr")[3:4], n=5)
	axis(2, at=xt)
	col <- rep(col,len=ncol(bp))
	border <- rep(border,len=ncol(bp))
	for (i in 1:ncol(bp)) plotbarheights(bp[,i], yadj=i, xadj=0.5, fg=ifelse(hollow,0,col[i]), bg=border[i], ...)
	title(main = main, sub = sub, xlab = xlab, ylab = ylab)

	if (labels) {
		for (i in 1:ncol(bp)) {
			#suppressWarnings(text(0:(nrow(bp)-1)*(ncol(bp)+1)+i+0.5,bp[,i],obp[,i],col=col[i],adj=c(0,0.5),srt=srt,cex=cex,...))
			xw <- which(obp[,i] != 0)
			suppressWarnings(text((1:nrow(bp)+(i-1)*xsep+0.5)[xw],(bp[,i]+(i-1)*ysep*ylim[2])[xw],obp[xw,i],col=col[i],adj=c(0,0.5),srt=srt,cex=cex,...))
		}
	}
	if (is.character(names.arg) || length(names.arg) > 1 || names.arg) {
		nm <- if(is.character(names.arg) || length(names.arg) > 1) names.arg else names(m)
		for (i in 1:ncol(bp)) axis(4,at=i,nm[i],col.axis=col[i], tick=FALSE, cex=cex, las=2, ...)
	}
	invisible(bp)
}




# a and b should be "counts"
dist.chisq.test <- function(a, b, count=NULL, breaks=seq(min(a,b), max(a,b), len=12), simulate.p.value=FALSE, B=2000, df=NULL) {
	if ((length(a) != length(b)) && !missing(count) && !count)  stop("a and b length's should be the same.\n")
	if ((!missing(count) && count) | any(a != as.integer(a)) | any(b != as.integer(b))) {
		if (length(breaks) == 1) breaks=seq(min(a,b), max(a,b), len=breaks)
		a <- hist(a, breaks=breaks, plot=FALSE)$count
		b <- hist(b, breaks=breaks, plot=FALSE)$count
	}
	a.s <- sum(a)
	b.s <- sum(b)
	p <- (a + b) / (a.s+b.s)
	ea <- p*a.s
	eb <- p*b.s
	obs <- c(a, b)
	e <- c(ea, eb)
	error <- any(e < 1)
	if (error) {
		w <- which(e < 1)
		e <- e[-w]
		obs <- obs[-w]
		warning(paste("Ommited", length(w), "expected values < 1."))
	}
	chi2 <- sum( (obs - e)^2 / e )
	names(chi2) <- "X-squared"
	dfree <- if (is.null(df)) (2-1)*(length(a)-1) else df
	names(dfree) <- "df"
	if (simulate.p.value) {
		almost.1 <- 1 - 64 * .Machine$double.eps
		ss <- 0
		n <- min(a.s, b.s)
		ex <- if(n==a.s) ea else eb
		nx <- length(a)
		for (.b in 1:B) {
			sm <- sum((table(factor(sample(1:nx, n, TRUE, prob = p), levels=1:nx)) - ex)^2/ex)
			ss <- ss + sum(sm >= almost.1 * chi2)
		}
		p.value <- (1 + ss)/(B + 1)
	} else {
		p.value <- pchisq(chi2, df=dfree, lower.tail=FALSE)
	}
	if (any(e < 5) || error) 
		warning("Chi-squared approximation may be incorrect, some expected values are less than 5.")
	structure(list(statistic = chi2, parameter = dfree, 
		p.value = p.value, method = paste("Chi Square Test",if (simulate.p.value) "(simulated)" else "", ", null hypothesis: there is no difference"), data.name = c(deparse(substitute(a)),deparse(substitute(b))), observed = obs, 
		expected = e, residuals = (obs - e)/sqrt(e)), class = "htest")
}



#from R legend, just added border.col
plot.legend <- function (x, y = NULL, legend, fill = NULL, col = "black", lty, 
	lwd, pch, angle = 45, density = NULL, bty = "o", bg = par("bg"), 
	pt.bg = NA, cex = 1, pt.cex = cex, pt.lwd = lwd, xjust = 0, 
	yjust = 1, x.intersp = 1, y.intersp = 1, adj = c(0, 0.5), 
	text.width = NULL, text.col = par("col"), merge = do.lines && 
		has.pch, trace = FALSE, plot = TRUE, ncol = 1, horiz = FALSE, 
	title = NULL, inset = 0, border.col="black") 
{
	if (missing(legend) && !missing(y) && (is.character(y) || 
		is.expression(y))) {
		legend <- y
		y <- NULL
	}
	mfill <- !missing(fill) || !missing(density)
	if (length(title) > 1) 
		stop("invalid title")
	n.leg <- if (is.call(legend)) 
		1
	else length(legend)
	if (n.leg == 0) 
		stop("'legend' is of length 0")
	auto <- if (is.character(x)) 
		match.arg(x, c("bottomright", "bottom", "bottomleft", 
			"left", "topleft", "top", "topright", "right", "center"))
	else NA
	if (is.na(auto)) {
		xy <- xy.coords(x, y)
		x <- xy$x
		y <- xy$y
		nx <- length(x)
		if (nx < 1 || nx > 2) 
			stop("invalid coordinate lengths")
	}
	else nx <- 0
	xlog <- par("xlog")
	ylog <- par("ylog")
	rect2 <- function(left, top, dx, dy, density = NULL, angle, 
		...) {
		r <- left + dx
		if (xlog) {
			left <- 10^left
			r <- 10^r
		}
		b <- top - dy
		if (ylog) {
			top <- 10^top
			b <- 10^b
		}
		rect(left, top, r, b, angle = angle, density = density, 
			...)
	}
	segments2 <- function(x1, y1, dx, dy, ...) {
		x2 <- x1 + dx
		if (xlog) {
			x1 <- 10^x1
			x2 <- 10^x2
		}
		y2 <- y1 + dy
		if (ylog) {
			y1 <- 10^y1
			y2 <- 10^y2
		}
		segments(x1, y1, x2, y2, ...)
	}
	points2 <- function(x, y, ...) {
		if (xlog) 
			x <- 10^x
		if (ylog) 
			y <- 10^y
		points(x, y, ...)
	}
	text2 <- function(x, y, ...) {
		if (xlog) 
			x <- 10^x
		if (ylog) 
			y <- 10^y
		text(x, y, ...)
	}
	if (trace) 
		catn <- function(...) do.call("cat", c(lapply(list(...), 
			formatC), list("\n")))
	cin <- par("cin")
	Cex <- cex * par("cex")
	if (is.null(text.width)) 
		text.width <- max(strwidth(legend, units = "user", cex = cex))
	else if (!is.numeric(text.width) || text.width < 0) 
		stop("'text.width' must be numeric, >= 0")
	xc <- Cex * xinch(cin[1], warn.log = FALSE)
	yc <- Cex * yinch(cin[2], warn.log = FALSE)
	xchar <- xc
	xextra <- 0
	yextra <- yc * (y.intersp - 1)
	ymax <- max(yc, strheight(legend, units = "user", cex = cex))
	ychar <- yextra + ymax
	if (trace) 
		catn("  xchar=", xchar, "; (yextra,ychar)=", c(yextra, 
			ychar))
	if (mfill) {
		xbox <- xc * 0.8
		ybox <- yc * 0.5
		dx.fill <- xbox
	}
	do.lines <- (!missing(lty) && (is.character(lty) || any(lty > 
		0))) || !missing(lwd)
	n.legpercol <- if (horiz) {
		if (ncol != 1) 
			warning("horizontal specification overrides: Number of columns := ", 
				n.leg)
		ncol <- n.leg
		1
	}
	else ceiling(n.leg/ncol)
	if (has.pch <- !missing(pch) && length(pch) > 0) {
		if (is.character(pch) && !is.na(pch[1]) && nchar(pch[1], 
			type = "c") > 1) {
			if (length(pch) > 1) 
				warning("not using pch[2..] since pch[1] has multiple chars")
			np <- nchar(pch[1], type = "c")
			pch <- substr(rep.int(pch[1], np), 1:np, 1:np)
		}
		if (!merge) 
			dx.pch <- x.intersp/2 * xchar
	}
	x.off <- if (merge) 
		-0.7
	else 0
	if (is.na(auto)) {
		if (xlog) 
			x <- log10(x)
		if (ylog) 
			y <- log10(y)
	}
	if (nx == 2) {
		x <- sort(x)
		y <- sort(y)
		left <- x[1]
		top <- y[2]
		w <- diff(x)
		h <- diff(y)
		w0 <- w/ncol
		x <- mean(x)
		y <- mean(y)
		if (missing(xjust)) 
			xjust <- 0.5
		if (missing(yjust)) 
			yjust <- 0.5
	}
	else {
		h <- (n.legpercol + (!is.null(title))) * ychar + yc
		w0 <- text.width + (x.intersp + 1) * xchar
		if (mfill) 
			w0 <- w0 + dx.fill
		if (has.pch && !merge) 
			w0 <- w0 + dx.pch
		if (do.lines) 
			w0 <- w0 + (2 + x.off) * xchar
		w <- ncol * w0 + 0.5 * xchar
		if (!is.null(title) && (tw <- strwidth(title, units = "user", 
			cex = cex) + 0.5 * xchar) > w) {
			xextra <- (tw - w)/2
			w <- tw
		}
		if (is.na(auto)) {
			left <- x - xjust * w
			top <- y + (1 - yjust) * h
		}
		else {
			usr <- par("usr")
			inset <- rep(inset, length.out = 2)
			insetx <- inset[1] * (usr[2] - usr[1])
			left <- switch(auto, bottomright = , topright = , 
				right = usr[2] - w - insetx, bottomleft = , left = , 
				topleft = usr[1] + insetx, bottom = , top = , 
				center = (usr[1] + usr[2] - w)/2)
			insety <- inset[2] * (usr[4] - usr[3])
			top <- switch(auto, bottomright = , bottom = , bottomleft = usr[3] + 
				h + insety, topleft = , top = , topright = usr[4] - 
				insety, left = , right = , center = (usr[3] + 
				usr[4] + h)/2)
		}
	}
	if (plot && bty != "n") {
		if (trace) 
			catn("  rect2(", left, ",", top, ", w=", w, ", h=", 
				h, ", ...)", sep = "")
		rect2(left, top, dx = w, dy = h, col = bg, density = NULL)
	}
	xt <- left + xchar + xextra + (w0 * rep.int(0:(ncol - 1), 
		rep.int(n.legpercol, ncol)))[1:n.leg]
	yt <- top - 0.5 * yextra - ymax - (rep.int(1:n.legpercol, 
		ncol)[1:n.leg] - 1 + (!is.null(title))) * ychar
	if (mfill) {
		if (plot) {
			fill <- rep(fill, length.out = n.leg)
			rect2(left = xt, top = yt + ybox/2, dx = xbox, dy = ybox, 
				col = fill, density = density, angle = angle, 
				border = border.col)
		}
		xt <- xt + dx.fill
	}
	if (plot && (has.pch || do.lines)) 
		col <- rep(col, length.out = n.leg)
	if (missing(lwd)) 
		lwd <- par("lwd")
	if (do.lines) {
		seg.len <- 2
		if (missing(lty)) 
			lty <- 1
		ok.l <- !is.na(lty) & (is.character(lty) | lty > 0)
		lty <- rep(lty, length.out = n.leg)
		lwd <- rep(lwd, length.out = n.leg)
		if (trace) 
			catn("  segments2(", xt[ok.l] + x.off * xchar, ",", 
				yt[ok.l], ", dx=", seg.len * xchar, ", dy=0, ...)")
		if (plot) 
			segments2(xt[ok.l] + x.off * xchar, yt[ok.l], dx = seg.len * 
				xchar, dy = 0, lty = lty[ok.l], lwd = lwd[ok.l], 
				col = col[ok.l])
		xt <- xt + (seg.len + x.off) * xchar
	}
	if (has.pch) {
		pch <- rep(pch, length.out = n.leg)
		pt.bg <- rep(pt.bg, length.out = n.leg)
		pt.cex <- rep(pt.cex, length.out = n.leg)
		pt.lwd <- rep(pt.lwd, length.out = n.leg)
		ok <- !is.na(pch) & (is.character(pch) | pch >= 0)
		x1 <- (if (merge) 
			xt - (seg.len/2) * xchar
		else xt)[ok]
		y1 <- yt[ok]
		if (trace) 
			catn("  points2(", x1, ",", y1, ", pch=", pch[ok], 
				", ...)")
		if (plot) 
			points2(x1, y1, pch = pch[ok], col = col[ok], cex = pt.cex[ok], 
				bg = pt.bg[ok], lwd = pt.lwd[ok])
		if (!merge) 
			xt <- xt + dx.pch
	}
	xt <- xt + x.intersp * xchar
	if (plot) {
		if (!is.null(title)) 
			text(left + w/2, top - ymax, labels = title, adj = c(0.5, 
				0), cex = cex, col = text.col)
		text2(xt, yt, labels = legend, adj = adj, cex = cex, 
			col = col) # text.col
	}
	invisible(list(rect = list(w = w, h = h, left = left, top = top), 
		text = list(x = xt, y = yt)))
}





plot.qq <- function(m, m2, m3, m4, las=0, col=if (is.data.frame(m) | is.matrix(m)) 1:ncol(m) else 1:length(m), pch=1:length(col), lty=rep(1:6, each=length(palette())), legends=TRUE,  legend.pos="topleft", legend.cols=1, legend.cex=1, len=NULL, type="o", ...) {
	if (!missing(m2) || !missing(m3) || !missing(m4)) {
		xn <- deparse(substitute(m))
		m <- list(m)
		if (!missing(m2)) { m[[length(m)+1]] <- m2; xn <- c(xn,deparse(substitute(m2))) }
		if (!missing(m3)) { m[[length(m)+1]] <- m3; xn <- c(xn,deparse(substitute(m3))) }
		if (!missing(m4)) { m[[length(m)+1]] <- m4; xn <- c(xn,deparse(substitute(m4))) }
		names(m) <- xn
		if (missing(col)) col <- 1:length(m)
	} else if (is.data.frame(m) || is.matrix(m)) {
		nm <- colnames(m)
		m <- lapply(1:ncol(m), function(i) m[,i])
		if (!is.null(nm)) names(m) <- nm
		if (missing(col)) col <- 1:length(m)
	} else if (is.vector(m) && !is.list(m)) {
		m <- list(m)
	}
	if (is.null(len)) len=length(m[[1]])
	m <- lapply(m, function(x) approx(1:length(x), sort(x), n=len)$y)
	plot(m[[1]], m[[1]], col=col[1], pch=pch[1], type=type, lty=lty[1], ...)
	for (i in 2:length(m)) {
		points(m[[1]], m[[i]], col=col[i], pch=pch[i], lty=lty[i], type=type)
	}
	if (length(legends) > 1 || (!is.character(legends) && legends)) plot.legend(legend.pos,if(length(legends) > 1) legends else names(m), lty=lty[1:length(m)], col=col[1:length(m)], pch=if (type!="l") pch[1:length(m)] else NULL, cex=legend.cex, ncol=legend.cols)
	invisible(m)
}




# impute missing values using knn
# m - matrix, rows - variables (genes), cols - samples (patients)
# n - number of neighbours samples to estimate the missing values
# distfunc - distance function
# verbose - progress info
# method - static=impute using always the same "n" neighbours, dynamic=impute using always the closest "n" neighbours that have no missing values for the imputing gene, dynamic is a bit slower
# ... - distfunc parameters
na.impute <- function(m, n=max(3,ncol(m)/10), imputfunc=mean, distfunc = dist, method=c("static","dynamic"), verbose=TRUE, min.values=2, ...) {

	method <- match.arg(method)
	if (verbose) { cat("Computing Distances...\n"); flush.console() }
	m.dist <- as.matrix(distfunc(t(m),...))
	if (verbose) { cat("Distance Matrix:", dim(m.dist),"\n") }

	if (verbose) { cat("Imputing..."); flush.console() }
	m.imp <- m
	for (i in 1:ncol(m)) {
		if (method=="static") {
			nni <- order(m.dist[,i],decreasing=FALSE)[1:n]
			if (verbose) { cat("Imputing Sample",i,", Neighbours:",nni,"\n"); flush.console() }
			wna <- which(is.na(m.imp[,i]))
			m.imp[wna,i] <- apply(m[wna,nni,drop=FALSE], 1, imputfunc, na.rm=TRUE)
		} else {
			nni <- order(m.dist[,i],decreasing=FALSE)
			if (verbose) { cat("Imputing Sample",i,"\n"); flush.console() }
			wna <- which(is.na(m.imp[,i]))
			for (j in wna) {
				nona <- na.omit(m[j,nni])
				if (length(nona) >= min.values) { 
					m.imp[j,i] <- imputfunc(nona[1:min(length(nona),n)])
				}
			}
		}
	}
	m.imp
}



# plot a scatter plot along with distributions of x and y
#, panel=function() {plot(0,0,type="n",xaxt="n",yaxt="n",main="",xlab="",ylab="")},
plot.xy.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),
	breaks=20,breaks.h=NULL,breaks.v=NULL,widths=c(5,2),heights=c(2.5,5),
	smooth=FALSE,smooth.col="red",smooth.lty=1,smooth.type="l",smooth.pch=1,smooth.lwd=1,
	restore.layout=TRUE,log="",hook.func=function(x,y) 0, nrpoints=50, nrcex=1, nrpch=".", nrcol=1, 
	labels=FALSE, ...) {
	fnx <- if (length(grep("x",log))) log10 else function(x) x
	fny <- if (length(grep("y",log))) log10 else function(y) y
	xp <- par(mar=c(1,4,1,1))
	layout(matrix(c(2,1,4,3),2,2),widths=widths,heights=heights,respect=TRUE)
	par(mar=c(5,4,0,0))
	plot(x,y,xlab=xlab,ylab=ylab,main="",log=log,...)
	ylim <- par("usr")[3:4]
	xlim <- par("usr")[1:2]
	if (smooth) {
		lines(loess.smooth(x,y,...),col=smooth.col,type=smooth.type,pch=smooth.pch,lty=smooth.lty,lwd=smooth.lwd)
	}
	hook.func(x,y)
	par(mar=c(1,3.5,4,0))
	cntx <- hist(fnx(x),breaks=if (missing(breaks.h) || is.null(breaks.h)) c(-Inf,seq(xlim[1],xlim[2],len=breaks),Inf) else breaks.h,plot=FALSE)$counts
	#print(cntx)
	cntx <- cntx[3:length(cntx)-1]
	bpv <- barplot(cntx,space=0,col="white",xlab="",main=main,xaxt="n")
	if (any(labels)) {
		text(bpv[,1], pmin(cntx,par("usr")[4]-strheight("X")), cntx, cex=ifelse(is.numeric(labels),labels,0.75), adj=c(0.5,0))
	}
	par(mar=c(4.5,1,0,1))
	#panel()
	cnty <- hist(fny(y),breaks=if (missing(breaks.v) || is.null(breaks.v)) c(-Inf,seq(ylim[1],ylim[2],len=breaks),Inf) else breaks.v,plot=FALSE)$counts
	#print(cnty)
	cnty <- cnty[3:length(cnty)-1]
	bph <- barplot(cnty,horiz=TRUE,space=0,col="white",xlab="Frequency")
	if (any(labels)) {
		text(pmin(cnty,par("usr")[2]-max(strwidth(bph))), bph[,1], cnty, cex=ifelse(is.numeric(labels),labels,0.75), adj=c(0,0.5))
	}
	par(mar=c(1,1,4,1))
	smoothScatter(fnx(x),fny(y),xaxt="n", yaxt="n", xlab="", ylab="", frame.plot=FALSE, postPlotHook=function(...) { }, nrpoints=nrpoints, cex=nrcex, pch=nrpch, col=nrcol, transformation=function(x) x^.4)
	if (restore.layout) {
		par(xp)
		layout(matrix(1))
	}
}

# plot a boxplot along with distributions of x and y
plot.box.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),breaks=20,breaks.h=NULL,breaks.v=NULL,widths=c(5,2),heights=c(2.5,5),smooth=FALSE,smooth.col="red",smooth.lty=1,smooth.type="l",smooth.pch=1,smooth.lwd=1,restore.layout=TRUE,overlap=TRUE,pch=20,pcol=rgb(.9,.9,.9),lwd=if (overlap) 2 else 1,hook.func=function(x,y) 0, ...) {
	xp <- par(mar=c(1,4,1,1))
	layout(matrix(c(1,2,0,3),2,2),widths=widths,heights=heights,respect=T)
	par(mar=c(1,4,4,0))
	xh <- hist(x, breaks=if (missing(breaks.h) || is.null(breaks.h)) breaks else breaks.h, xlab="",main=main,xaxt="n")
	par(mar=c(5,4,0,0))
	lab <- cut(x, breaks=xh$breaks) # if (missing(breaks.h) || is.null(breaks.h)) breaks else breaks.h)
	plot(x,y,xlab=xlab,ylab=ylab,main="",pch=pch,col=pcol,type=if (overlap) "p" else "n")
	xr <- max(x)-min(x)
	ur <- par("usr")[2]-par("usr")[1]
	xmas <- nlevels(lab)*(ur/xr-1)
	par(usr=c(0, nlevels(lab)+1, par("usr")[3], par("usr")[4]))
	boxplot(y~lab,xlab=xlab,ylab=ylab,main="",xaxt="n",yaxt="n",add=TRUE,lwd=lwd,...)
	if (smooth) {
		lines(loess.smooth(x,y,...),col=smooth.col,type=smooth.type,pch=smooth.pch,lty=smooth.lty,lwd=smooth.lwd)
	}
	hook.func(x,y)
	par(mar=c(5,1,0,1))
	#panel()
	barplot(hist(y,breaks=if (missing(breaks.v) || is.null(breaks.v)) breaks else breaks.v,plot=F)$counts,horiz=T,space=0,col="white",xlab="Frequency")
	if (restore.layout) {
		par(xp)
		layout(matrix(1))
	}
}

# plot a scatter plot along with distributions of x and y
#, panel=function() {plot(0,0,type="n",xaxt="n",yaxt="n",main="",xlab="",ylab="")},
matplot.xy.hist <- function(x,y,main="",xlab=deparse(substitute(x)),ylab=deparse(substitute(y)),breaks=20,widths=c(5,2),heights=c(2.5,5),smooth=FALSE,smooth.col="red",smooth.lty=1,smooth.type="l",smooth.pch=1,smooth.lwd=1,...) {
	xp <- par(mar=c(1,4,1,1))
	layout(matrix(c(1,2,0,3),2,2),widths=widths,heights=heights,respect=T)
	par(mar=c(1,4,4,0))
	plot.densities(x,breaks=breaks,xlab="",ylab="Frequency",main=main,xaxt="n")
	par(mar=c(5,4,0,0))
	matplot(x,y,xlab=xlab,ylab=ylab,main="",...)
	if (smooth) {
		lines(loess.smooth(x,y,...),col=smooth.col,type=smooth.type,pch=smooth.pch,lty=smooth.lty,lwd=smooth.lwd)
	}
	par(mar=c(5,1,0,1))
	#panel()
	plot.densities(y,breaks=breaks,xlab="Frequency",main=main,yaxt="n",horiz=TRUE)
	par(xp)
	layout(matrix(1))
}


plot.2scales <- function(func=hist, ..., main=NULL, xlab=NULL, xaxt=NULL, ylim1=c(.9,1), ylim2=c(0,.25), heights=c(0.25,0.75), axes=list(c(2),if(is.null(xaxt) || xaxt=="n") c(2) else c(1,2)), axes1=FALSE, axes2=FALSE, panel1=NULL, panel2=NULL, restore.par=TRUE) {
	xpar <- par()
	xmar <- par("mar")
	layout(matrix(1:2,nrow=2),widths=c(1,1),heights=heights)
	#par(mfrow=c(2,1))
	par(mar=xmar/c(5,1,1,1))
	suppressWarnings(func(..., xaxt="n", xlab="", ylim=ylim1, main=main, axes=axes1))
	if (length(axes) > 0  && is.list(axes)) suppressWarnings(sapply(axes[[1]], axis))
	panel1
	par(mar=xmar/c(1,1,10,1))
	suppressWarnings(func(..., ylim=ylim2, main="", xlab=xlab, xaxt=xaxt, axes=axes2))
	if (length(axes) > 1  && is.list(axes)) suppressWarnings(sapply(axes[[2]], axis))
	panel2
	if (restore.par) suppressWarnings(par(xpar))
}

#y <- list(A=rnorm(10000), B=runif(10000,min=-4,max=4))
#plot.2scales2(panel1=plot.hist2(y, ylim=c(1500,2000), axes=FALSE, xaxt="n", yaxt="s"), panel2=plot.hist2(y, ylim=c(0,1000), xaxt="s", yaxt="s", axes=FALSE, legends=FALSE))
plot.2scales2 <- function(heights=c(0.25,0.75), panel1=NULL, panel2=NULL) {
	xpar <- par()
	xmar <- par("mar")
	layout(matrix(1:2,nrow=2),widths=c(1,1),heights=heights)
	par(mar=xmar/c(5,1,1,1))
	panel1
	par(mar=xmar/c(1,1,10,1))
	panel2
	suppressWarnings(par(xpar))
}



#plot a barplot then labels up to the bars
plot.bars <- function(xt, space=0.25, adj=if (horiz) c(0, 0.5) else c(0.5, -0.25), 
	labels=paste(xt,"\n",round(xt*100/sum(xt),labels.round),"%",sep=""), 
	cex.text=1, blim=c(0, max(xt)*ifelse(horiz,1.6,1.25)), col=8, 
	labels.round=3, omit.zeros=TRUE, horiz=FALSE, labels.col=col, 
	cex.axis=1, mar=par("mar"), restore.par=TRUE, axes=FALSE, ...) {
	xpar <- par()
	par(mar=mar)
	if (restore.par) on.exit(par(xpar))
	if (horiz) {
		xlim <- blim
		ylim <- NULL
	} else {
		ylim <- blim
		xlim <- NULL
	}
	if (blim[1] > 0) {
		barplot(xt-ylim[1], space=space, ylim=ylim-ylim[1], xlim=xlim, axes=axes, col=col, horiz=horiz, ...)
		axis(2, at=axTicks(2), labels=axTicks(2) + ylim[1])
	} else  { #if (ylim[1] == 0)
		barplot(xt, space=space, ylim=ylim, col=col, horiz=horiz, xlim=xlim, axes=axes, axisnames=FALSE ,...)
		side <- ifelse(horiz,2,1)
		col <- circularize(col, length(xt))
		for (t in 1:length(xt)) {
			axis(side,at=t-0.5+t*space,labels=names(xt)[t],col.axis=col[t],cex.axis=cex.axis,...)
		}
	}
	labels.col[labels.col == 0 | labels.col == "white" | labels.col == "#FFFFFF"] <- 1
	if (horiz) {
		text(pmin(xt-xlim[1],par("usr")[2]*.9), 1:length(xt)-0.5+1:length(xt)*space, ifelse(xt==0 & omit.zeros, "", gsub("\n",", ",labels)), adj=adj, cex=cex.text, col=labels.col)
	} else {
		text(1:length(xt)-0.5+1:length(xt)*space, pmin(xt-ylim[1],par("usr")[4]*.9), ifelse(xt==0 & omit.zeros, "", labels), adj=adj, cex=cex.text, col=labels.col)
	}
}



plot.two.bars <- function(h1, h2, weight=.8, space.center=1, space.left=mean(h1)/4, space.right=mean(h2)/4, col=1:length(h1), border=1:length(h1), labels=TRUE, cex=-1, font=1, mar=c(0,0,3,0), zeroes=TRUE, main="h1 vs h2", center=FALSE, maxHeight=max(c(h1,h2))) {

	stopifnot(length(h1)==length(h2))
	xpar <- par()
	par(mar=mar)
	on.exit(par(xpar))
	limitedHeight <- function(h) {
		pmin(h, maxHeight)+ifelse(h > maxHeight,sqrt(h),0)
	}
	lh1 <- limitedHeight(h1)
	lh2 <- limitedHeight(h2)
	mx1 <- max(lh1)
	mx2 <- max(lh2)
	mx <- max(mx1,mx2)
	plot(0,0,xlim=c(-mx1,mx2), ylim=c(0.9,length(h1)+0.1), frame.plot=FALSE, axes=FALSE, xlab="", ylab="", main=ifelse(length(main) > 1, "", main), type="n")
	if (length(main) > 1) {
		title(main[1], adj=0)
		if (length(main) > 2) {
			title(main[3], adj=1)
			title(main[2], adj=c(0.5))
		} else {
			title(main[2], adj=1)
		}
	}
	if (cex < 0) {
		cex <- max(0.1, 1.1 * weight/strheight("[](){}GgJj", font=font))
		while (cex > 0.1 & strheight("[](){}GgJj", cex=cex, font=font) > weight) {
			cex <- cex - 0.01
		}
	}
	spacec <- strwidth(paste(rep(" ",space.center),collapse=""))
	mxw1 <- max(strwidth(names(h1),cex=cex,font=font,units="inches"))
	mxw2 <- max(strwidth(names(h2),cex=cex,font=font,units="inches"))
	mxrat1 <- (mxw1)/par("fin")[1]
	mxrat2 <- (mxw2)/par("fin")[1]
	par(usr=c(-(mx1+space.left),(mx2+space.right),par("usr")[3:4]))
	abline(v=-spacec)
	symbols(-lh1/2-spacec, 1:length(h1), inches=FALSE, rectangles=data.matrix(data.frame(width=as.numeric(lh1),height=weight)), bg=col, fg=border, add=TRUE)
	abline(v=spacec)
	symbols(lh2/2+spacec, 1:length(h2), inches=FALSE, rectangles=data.matrix(data.frame(width=as.numeric(lh2),height=weight)), bg=col, fg=border, add=TRUE)
	if (center) {
		text(0, 1:length(h1)+(1-weight)/2, names(h1), col=col, cex=cex, font=font, adj=c(0.5,0.5))
	} else {
		text(par("usr")[1], 1:length(h1)+(1-weight)/2, names(h1), col=col, cex=cex, font=font, adj=c(0,0.5))
		text(par("usr")[2], 1:length(h2)+(1-weight)/2, names(h2), col=col, cex=cex, font=font, adj=c(1,0.5))
	} 
	if (labels) {
		h1s <- as.character(h1)
		if (!zeroes && any(h1s=="0")) h1s[h1s=="0"] <- ""
		h2s <- as.character(h2)
		if (!zeroes && any(h2s=="0")) h2s[h2s=="0"] <- ""
		#text(-lh1-spacec, 1:length(h1)+(1-weight)/2, h1s, col=col, cex=cex, font=font, adj=c(1,0.5))
		#text( lh2+spacec, 1:length(h2)+(1-weight)/2, h2s, col=col, cex=cex, font=font, adj=c(0,0.5))
		text(-lh1-spacec, -0.45+1:length(h1)+(weight)/2, h1s, col=col, cex=cex, font=font, adj=c(1,0.5))
		text( lh2+spacec, -0.45+1:length(h2)+(weight)/2, h2s, col=col, cex=cex, font=font, adj=c(0,0.5))
	}
}



#x is optional for x coordinates
plot.signatures <- function(mat, xaxis, type="l", pch=ifelse(type=="l",-1,20), splice=0, lwd=1, col=1:ncol(mat), lty=rep(1:6, each=length(palette())), legends=FALSE, legend.pos="topleft", legend.cex=1, legend.cols=1, labels=FALSE, labels.dig=NULL, add=FALSE, nsplines=0, ...) {
	if (missing(xaxis)) xaxis <- 1:nrow(mat)
	x <- mat[,1]
	x[1:2] <- c(max(mat)+splice*(ncol(mat)-1),min(mat))
	if (!add) plot(xaxis, x, type="n", ...)
	lty <- rep(lty, ncol(mat), length.out=ncol(mat))
	col <- rep(col, ncol(mat), length.out=ncol(mat))
	for (i in 1:ncol(mat)) {
		x <- xaxis
		y <- mat[,i]
		typex <- type
		if (type %in% c("p","o")) {
			suppressWarnings(points(x, y+splice*(i-1),col=col[i],type="p",pch=pch[(i-1) %% length(pch)+1],lty=lty[i], lwd=lwd, ...))
			typex <- "l"
		}
		if (nsplines > 0) {
			library(stats)
			s <- spline(x,y,nsplines)
			x <- s$x
			y <- s$y
		}
		suppressWarnings(lines(x, y+splice*(i-1),col=col[i],type=typex,pch=pch[(i-1) %% length(pch)+1],lty=lty[i], lwd=lwd, ...))
	}
	if (labels) {
		for (i in 1:ncol(mat)) {
			x <- xaxis
			y <- mat[,i]
			if (nsplines > 0) {
				library(stats)
				s <- spline(x,y,nsplines)
				x <- s$x
				y <- s$y
			}
			suppressWarnings(text(x=x, y=y+splice*(i-1), labels=if(is.null(labels.dig)) y+splice*(i-1) else round(y+splice*(i-1), labels.dig), ...))
		}
	}
	if (length(legends) > 1 || (!is.character(legends) && legends)) suppressWarnings(plot.legend(legend.pos, if(length(legends) > 1) legends else colnames(mat), col=col[1:ncol(mat)], lty=lty, pch=pch, lwd=lwd, cex=legend.cex, ncol=legend.cols))
}

plot.profiles <- plot.signatures





plot.profiles.msd <- function(mat, xaxis,  splice=0, ...) {

	if (is.matrix(mat)) mat <- lapply(apply(mat, 2, list), unlist)
	if (missing(xaxis)) xaxis <- 1:length(mat)
	x <- mat[[1]]
	x <- c(max(unlist(lapply(mat,max)))+splice*(length(mat)-1),min(unlist(lapply(mat,min))))
   	plot(xaxis, rep(x, length(xaxis), length.out=length(xaxis)), type="n", ...)
	media <- unlist(lapply(mat, mean))
	desv <- unlist(lapply(mat, sd))
	lines(xaxis, media+desv, lty=2)
	lines(xaxis, media+desv*2, lty=3)
	lines(xaxis, media, lty=1)
	lines(xaxis, media-desv*2, lty=3)
	lines(xaxis, media-desv, lty=2)
}


plot.profile.range <- function(x, m, breaks=20, q=0.9, 
	col=1:length(m), pch=1:length(m), type=rep("o", length(m)+1), 
	lty=1:(length(q)+1), ylim=range(unlist(l)), lwd=(length(q)+1):1,
	xlim=range(breaks), labels=NULL, legend.pos="bottomright", legend.ncol=1, 
	q.show.pch=rep(TRUE, length(m)), q.show=rep(TRUE,length(m)), q0.show=rep(TRUE,length(m)), ...) {

	# x is the xaxis coordinates of the profiles
	# m is a list of profiles matrices
	
	#if (is.matrix(m) || is.data.frame(m)) {
	#	m <- list()
	#	for (i in 1:ncol(m)) m[[i]] <- m[,i]
	#}
	
	h <- hist(x, breaks=breaks, plot=FALSE)
	breaks <- h$breaks
	mids <- h$mids
	xidx <- sapply(h$mids, function(mid) which.min(abs(x-mid)))
	
	q <- ifelse(q < 0.5, 1-q, q)
	qs <- c(0.5, q, 1-q)
	rx <- matrix(0, ncol=length(qs), nrow=length(mids))
	l <- list()
	for (i in 1:length(m)) {
		r2 <- rx
		for (j in 1:length(mids)) {
			a <- apply(m[[i]], 2, function(y) y[xidx[j]])
			r2[j,] <- quantile(a, qs)
		}
		l[[i]] <- r2
	}
	plot(x,rep(0,length(x)), type="n", xlim=xlim, ylim=ylim, ...)
	idx <- c(1, 1:length(q) + 1, 1:length(q) + 1)
	for (i in 1:length(l)) {
		for (j in 1:ncol(l[[i]])) {
			if ((j == 1  && q0.show[i]) || (j > 1  &&  q.show[i])) {
				suppressWarnings(lines(mids, l[[i]][,j], col=col[i], pch=ifelse(q.show.pch[i]  || j == 1, pch[i], ""), lty=lty[ifelse(length(q)==0,i,idx[j])], type=type[i], lwd=lwd[ifelse(length(q)==0,i,idx[j])], ...))
			}
		}
	}
	if (all(is.null(labels))) {
		if (!is.null(names(m))) labels <- names(m) else labels <- 1:length(m)
	}
	if (!any(is.na(labels))) legend(legend.pos, legend=labels, lty=lty, col=col, pch=ifelse(type=="o",pch,""), lwd=lwd[1], ncol=legend.ncol)
	invisible(l)
}


## to smooth a vector of values from its neighbours
smooth.window <- function(v, size=1, left=size, right=size, func=mean, ...) {
	u <- 1:length(v)
	from <- pmax(u - left, 1)
	to   <- pmin(u + right, length(v))
	a <- v
	for (i in 1:length(v)) a[i] <- func(v[from[i]:to[i]], ...)
	a
}

smooth.window.matrix <- function(mat, ...) apply(mat, 2, smooth.window, ...)
smooth.window.data.frame <- function(mat, ...) apply(data.matrix(mat), 2, smooth.window, ...)

smooth.y.by.x.window <- function(x, y, size=max(2,length(x) / 1000), func=mean) {
	o <- order(x)
	s <- smooth.window(y[o], size=size, func=func)
	z <- y
	z[o] <- s
	z
}
 
smooth.y.by.x.hist <- function(x, y, breaks=100, func=mean) {
	if (length(breaks) == 1) breaks <- seq(min(x)-0.0001, max(x), length.out=breaks)
	xc <- cut(x, breaks=breaks, labels=FALSE)
	u <- sort(unique(xc))
	z <- numeric(length(u))
	w <- numeric(length(u))
	for (i in 1:length(u)) {
		f <- xc == u[i]
		w[i] <- mean(x[f])
		z[i] <- func(y[f])
	}
	list(x=w, y=z)
}


predict.smooth.y.by.x.hist <- function(w, z, x) {
	y <- numeric(length(x))
	xx <- pmin(x, max(w))
	for (i in 1:length(x)) {
		b <- max(2,which(w >= xx[i])[1])
		a <- max(1,b-1)
		dx <- w[b] - w[a]
		dy <- z[b] - z[a]
		y[i] <- z[a] + (x[i] - w[a]) * (dy / dx)
	}
	y
}

colourise.dendrogram <- function(d, use.order=TRUE, ...) {
	if (all(class(d) != "dendrogram")) d <- as.dendrogram(d)
	if (! missing(...)) {
		local(
		colLab <<- function(n, ...) {
		   leafs <- unlist(lapply(n,is.leaf))
		   if(is.leaf(n) || (length(n)==2 & any(leafs))) {
			 a <- attributes(n)
			 i <- a$index
			 if (is.null(i)) i <- attributes(n[[which(leafs)[1]]])$index
			 attr(n, "index") <- i
			 #lab.col=
			 #lab.font=
			 x <- list(...)[[1]]
			 #thePar <- ifelse(is.leaf(n), "nodePar", "edgePar")
			 attr(n, "nodePar") <- c(a$nodePar, lapply(x, function(z, i) z[((i-1) %% length(z))+1], o[i]))
			 attr(n, "edgePar") <- c(a$edgePar, lapply(x, function(z, i) z[((i-1) %% length(z))+1], o[i]))
		   }
		   n
		})
		local(
		ordLab <<- function(n, ...) {
		   if(is.leaf(n)) {
			 i <<- i+1
			 attr(n, "index") <- i
		   }
		   n
		})
		i <- 0
		o <- order.dendrogram(d)
		if (!use.order) o <- 1:length(o)
		dc <- dendrapply(d, ordLab, ...)
		dc <- dendrapply(dc, colLab, ...)
		dc
	} else {
		d
	}
}

#dend <- colourise.dendrogram(as.dendrogram(hcl, hang=0.1), label.list=list(lab.col=xcol3, lab.font=2), ...)
#variables in label.list: 
#plot(dend)
#hc <- hclust(dist(USArrests), "ave")
#den <- as.dendrogram(hc) 
#plot.hclust.colours(den, label.list=list(lab.col=rep(1:5, each=10), col=rep(1:5, each=10), pch=19, bg=2, lty=1:3, lwd=2), labels.sorted=TRUE)
plot.hclust.colours <- function(data, distfun=dist, hclustfun=hclust, label.list=list(lab.col=2, lab.font=1, lab.cex=1), labels.sorted=FALSE, hang=0.1, ...) {
	if (all(class(data)!="hclust") && all(class(data)!="dendrogram")) {
		if (any(class(data)!="dist")) xdist <- distfun(data) else xdist <- data
		hc <- hclustfun(xdist)
	} else {
		hc <- data
	}
	if (all(class(hc) != "dendrogram")) hc <- as.dendrogram(hc, hang=hang)
	plot(colourise.dendrogram(hc, !labels.sorted, label.list), ...)
}



## plot a hierarchical clustering with colours 
## cutting the tree at specified height (h) or cluster number (k), thus k or h must be specified
plot.hclust.colored.clusters <- function(data, distfun=dist, hclustfun=hclust, h=NULL, k=NULL, hang=NULL, 
	label.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1, lab.cex=1),
	branch.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1),
	upper.attr=list(col="#888888", pch=1, lty=1, lwd=1.5, cex=1), 
	silhouettes=NULL, sil.labadj=FALSE, sil.off=ifelse(sil.labadj,ifelse(sil.image,0.05,0.01),0.3), sil.size=0.7, 
	ylab="Height", sil.ylab="Silhouette", sil.col=redgreenblue(-1), sep=TRUE, numbers=1.1, ## if 0 or FALSE
	sil.image=TRUE, side.factors=NULL, factors.size=1, factors.cex=1,...) {

	
	
	if (!any(class(data)=="hclust") && !any(class(data)=="dendrogram")) {
		if (!any(class(data)=="dist")) xdist <- distfun(data) else xdist <- data
		if (is.null(k) && h < 0 && h >= -1) h <- quantile(xdist, abs(h))
		if (is.logical(silhouettes) && silhouettes) silhouettes <- xdist
		hc <- hclustfun(xdist)
	} else {
		hc <- data
	}
	hc2 <- hc
	cl <- cutree(hc, h=h, k=k)
	heights <- range(hc$height)
	n <- length(unique(cl))

	if (all(class(hc) != "dendrogram")) {
		hc <- if (is.null(hang)) as.dendrogram(hc) else as.dendrogram(hc, hang=hang)
	}
	
	## to specify in parameters only those that changes	
	.label.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1, lab.cex=1)
	for (i in names(label.attr)) .label.attr[[i]] <- label.attr[[i]]
	 
	.branch.attr=list(lab.col=1:n, col=1:n, pch="", lty=1, lwd=1, cex=0.5, p.col=1:n, p.lwd=1, p.lty=1)
	for (i in names(branch.attr)) .branch.attr[[i]] <- branch.attr[[i]]

	.upper.attr=list(col="#888888", pch=1, lty=1, lwd=1.5, cex=1)
	for (i in names(upper.attr)) .upper.attr[[i]] <- upper.attr[[i]]
	
	## update the formal variables
	label.attr <- .label.attr
	branch.attr <- .branch.attr
	upper.attr <- .upper.attr

	applyAttr <- function(values, k) {
		lapply(values, function(z, i) z[((i-1) %% length(z))+1], k)
	}
	local(
		colLeaf <<- function(n, ...) {
	   if(is.leaf(n)) {
		 a <- attributes(n)
		 i <<- i+1
		 k <- cl[o[i]]
		 attr(n, "index") <- i
		 attr(n, "cluster") <- k
		 #attr(n, "nodePar") <- c(a$nodePar, lapply(label.attr, function(z, i) z[((i-1) %% length(z))+1], k))
		 #attr(n, "edgePar") <- c(a$edgePar, lapply(label.attr, function(z, i) z[((i-1) %% length(z))+1], k))
		 attr(n, "nodePar") <- c(a$nodePar, applyAttr(label.attr, k))
		 attr(n, "edgePar") <- c(a$edgePar, applyAttr(label.attr, k))
	   }
	   n
	})
	
	local(
		colBranch <<- function(n, ...) {
		 if(!is.leaf(n)) {
		 	.leafs <- c()
		 	extLeaf <- function(x) {
		 		if (is.leaf(x)) .leafs[length(.leafs)+1] <<- attr(x, "cluster") 
		 		else lapply(x, extLeaf)
		 	}
		 	extLeaf(n)
		 	if (length(.leafs) > 0) {
			 	clusters <- .leafs #unlist(lapply(.leafs, function(x) attr(x, "cluster")))
			 	xt <- table(clusters)
			 	wm <- which.max(xt)
			 	k <- as.numeric(names(xt)[wm])
				a <- attributes(n)
				if (xt[wm] == length(.leafs)) {
					 attr(n, "nodePar") <- c(a$nodePar, applyAttr(branch.attr, k))
					 attr(n, "edgePar") <- c(a$edgePar, applyAttr(branch.attr, k))
				} else {
					 attr(n, "nodePar") <- c(a$nodePar, applyAttr(upper.attr, k))
					 attr(n, "edgePar") <- c(a$edgePar, applyAttr(upper.attr, k))
				}
		 	}
		 }
		 n
	  })
	
	o <- order.dendrogram(hc)
	i <- 0
	#dc <- hc
	dc <- dendrapply(hc, colLeaf)
	dc <- dendrapply(dc, colBranch)
	#dc
	if (!is.null(k)) {
		# find height
		mi <- heights[1]
		mf <- heights[2]
		while (mf-mi > 0.01) {
			m <- (mi+mf)/2
			L <- cutree(hc2, h=m)
			if (length(unique(L)) > k) {
				mi <- m
			} else {
				mf <- m
			}
		}
		h = m
	}
	hh <- 0
	if (!is.null(silhouettes) || !is.null(side.factors)) {
		opar <- par("mar")
		on.exit(par(mar=opar))
		mar <- par("mar")
		par(mar=c(0, mar[-1]))
		mn <- min(0,min(heights))
		d <- max(heights)-mn
		sf <- sft <- 0
		if (!is.null(side.factors)) {
				#plot(dc, ylim=c(mn-d*sil.size, max(heights)), yaxt="n", ...)
			sf <- if (!is.null(side.factors)) ncol(side.factors)*factors.size*d/20 else 0
			#sft <- if (!is.null(side.factors)) ncol(side.factors) * strheight("JjYygGpPQq", cex=factors.cex)*1.15 else 0
			sft <- if (!is.null(side.factors)) ncol(side.factors)*factors.cex*d/20 else 0
			sf <- sf + sft
		}
			plot(dc, ylim=c(mn-d*sil.size-sf, max(heights)), yaxt="n", ...)
			sft <- max(if (!is.null(side.factors)) ncol(side.factors) * strheight("[](){}JjYygGpPQq", cex=factors.cex)*1.15 else 0,sft)
			mtext(ylab, 2, at=max(heights)/2, line=2)
			axis(2, pretty(c(0,heights)))
			if (is.matrix(data)) {
				rn <- rownames(data)
				nr <- nrow(data)
			} else {
				rn <- attr(data, "Labels")
				nr <- max(1,length(attr(data, "Labels")))
			}
			inh <- sapply(paste(rn,"  "), strwidth, cex=label.attr$lab.cex, units="inches")
			rny <- (par("usr")[4]-par("usr")[3])*(mean(inh)+sd(inh))/(par("pin")[2])
			
			##### Factores #####
			if (!is.null(side.factors)) {
				sfyy <- -rny
				sfy <- -rny-sf
				syp <- seq(sfy, sfyy, by=sf/ncol(side.factors))
				for (i in 1:ncol(side.factors)) {
					f <- side.factors[o,i]
					if (!is.factor(f)) { # &&  (length(unique(f)) <= nr / 10 || is.character(f) || is.logical(f))
						f <- factor(f)
					}
					if (is.factor(f)) {
						lf <- levels(f)
						f <- unclass(f)
					} else {
						lf <- unique(f)
					}
					u <- length(lf)
					fcol <- if(u > length(palette())) rainbow(u*1.1)[1:u] else palette()[1:u]
					.y <- syp[0:1+i]+c(sft/ncol(side.factors),0)
					image(1:nr, .y, matrix(f,ncol=1), col=fcol, add=TRUE, breaks=c(0,1:u))
					lu <- paste(c(colnames(side.factors)[i],as.character(lf))," ")
					abline(h=.y[1], col="white",lwd=1)
					text(left(cumsum(c(par("usr")[1],sapply(lu,strwidth,cex=factors.cex))),-1), .y[1], lu, col=c(1,fcol), cex=factors.cex, adj=c(0,1))
					#text(par("usr")[1],.y[2],colnames(side.factors)[i], adj=c(0,1))
				}
			}
			
			
			
			if (!is.null(silhouettes)) {
			
				library(cluster)
				sil <- silhouette(cl, silhouettes)
				if ("silhouette" %in% class(sil)) {
					rs <- range(c(0,sil[,"sil_width"]))
					s1 <- mn-d*sil.off
					s2 <- mn-d*sil.size
					if (sil.image) {
						image( 1:(length(o)+1)-0.5, c(-d/25, -d/50, -d/100)-sf, matrix(rep(sil[o,"sil_width"],2),ncol=2), add=TRUE, breaks=c(seq(-1,0,length.out=length(sil.col)/2), seq(0,1,length.out=length(sil.col)/2+1))[1:(length(sil.col)+1)],col=sil.col)
						text(0,-d/60-sf,round(min(sil[,"sil_width"]),2),col=sil.col[1],adj=c(1,1),cex=0.5)
						text(0,-d/60-sf,round(max(sil[,"sil_width"]),2),col=sil.col[length(sil.col)],adj=c(1,0),cex=0.5)
					}
					if (sil.size > 0) {
						silcoor <- function(s) { s2+(s-rs[1])*(s1-s2)/(rs[2]-rs[1]) }
						hh <- silcoor(0)
						kol <- branch.attr$col
						for (i in 1:length(o)) {
							lines(c(i,i), c(hh, silcoor(sil[o[i],"sil_width"]))-sf, col=kol[cl[o[i]]])
						}
						abline(h=silcoor(0)-sf)
						axis(4, at=silcoor(round(c(0,rs),2))-sf, labels=round(c(0,rs),2), line=-2.5)
						mtext(sil.ylab, 4, at=silcoor(mean(rs)*0), line=-0.66, cex=0.66)
						hh <- s2
					}
				}
			}
	} else {
			plot(dc, ylab=ylab, ...)
	}
	if (sep || numbers) {
		j <- 1
		ia <- 1
		lcex <- 0.75
		if (!is.null(.label.attr$lab.cex)) 
			lcex <- .label.attr$lab.cex[1] * ifelse(is.logical(numbers), 1.15, numbers)
	   	for (i in 1:(length(o)+1)) {
	   		if (i == 1 || i==(length(o)+1) || cl[o[i]] != cl[o[i-1]]) { 
	   			if (sep) lines(c(i-0.5,i-0.5), c(hh, h), lty=2, col=8)
	   			if (i > 1) {
	   				if (numbers) text(ia+(i-ia-1)/2, hh, j, col=.label.attr$lab.col[cl[o[i-1]]], cex=lcex)
	   				j <- j + 1
	   			}
	   			ia <- i
	   		}
	   	}
	}
   	abline(h=h, col=8, lty=2)
   	text(par("usr")[1], h, round(h,3), adj=c(0,0), cex=0.75)
   	## renumber cluster numbers by appearance
   	xcl <- cl - cl
   	cc <- 1
   	for (i in 1:length(cl)) {
   		if (xcl[o[i]]==0) {
   			xcl[cl == cl[o[i]]] <- cc
   			cc <- cc + 1
   		}
   	}
   	cl <- xcl
	invisible(list(cluster=cl,order=o))
}



#filter values in a vector depending on its quantiles, 
filterize <- function(x, q1, q2, func=NULL) { 
	a <- sort(quantile(x,p=c(q1,q2)))
	x <- x[x >= a[1] & x <= a[2]]
	if (!missing(func) && length(x) > 0) x <- func(x)
	x
}


my.quantile <- function(v, probs, na.rm=FALSE) {
	s <- sort(v, na.last=ifelse(na.rm,NA,TRUE))
	l <- length(s)
	s[pmax(pmin(l,probs*l),1)]
}

#e.mx - matrix to normalize
#center - use median to center e.mx to the median of medians
#averages - returns averages
#q.func - quantile function
#m.func - in the case averages is true
quantile.normalization.old <- function(e.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, m.func=mean) {
	#q.func=quantile	# if you want to use the embed slower R quantile (slightly different results)
	m <- 0 #nrow(e.mx)
	if (center) {
		md <- apply(e.mx, 2, median, na.rm = TRUE)
		mmd <- median(md)
		e.mx <- sweep(data.matrix(e.mx), 2, md-mmd)
	}
	for (i in 1:ncol(e.mx)) m <- max(m,sum(is.finite(e.mx[,i])))
	if (verbose) cat("Max-Quantiles:",m,"\n")
	q.mx <- matrix(0,nrow=m, ncol=ncol(e.mx))
	for (i in 1:ncol(e.mx)) { q.mx[,i] <- q.func(e.mx[,i],probs=(1:m)/m,na.rm=TRUE); if (verbose) cat("Quantiled:",i,"\n"); }
	avg <- numeric(length=m)
	if (verbose) cat("Averaging...\n")
	for (i in 1:m) avg[i] <- m.func(q.mx[i,], na.rm=TRUE) 
	#require(stats)
	q.mx <- matrix(0,nrow=nrow(e.mx), ncol=ncol(e.mx), dimnames=dimnames(e.mx))
	if (verbose) cat("Normalizing...\n")
	for (i in 1:ncol(q.mx)) { 
		if (verbose) cat("Computing ecdf",i,"...\n");
		avg.cdf <- ecdf(e.mx[is.finite(e.mx[,i]),i])
		r <- pmin(pmax(1,avg.cdf(e.mx[,i])*m),m)
		q.mx[,i] <- avg[r]
	}
	if (verbose) cat("Done!\n");
	if (averages) avg else q.mx
}

#center - median centering before quantile (almost no effect, just create an offset)
#qmin - scale the distribution forcing the minimum to be qmin shriniking or stretching the distribution
#qmax - scale the distribution forcing the maximum to be qmax shriniking or stretching the distribution
#qmid - scale the centered qmid part of the distribution to be between qmin and qmax shriniking or stretching the distribution
#qshrink - scale every element to have the centered qshrink part of the distribution equal before estimating the averages (helps to smooth eventual wild distributions)
quantile.normalization <- function(l.mx, verbose=FALSE, averages=FALSE, center=FALSE, q.func=my.quantile, m.func=mean, qmin=min(avg), qmax=max(avg), qmid=1, qshrink=NULL) {
	#q.func=quantile	# if you want to use the embed slower R quantile (slightly different results)
	
	#convert to list
	retfn <- function(x) x
	if (is.data.frame(l.mx) || is.matrix(l.mx)) { 
		.rn <- rownames(l.mx)
		l.mx <- lapply(as.list(data.frame(l.mx)), function(x) { names(x) <- .rn; x; })
		if (is.data.frame(l.mx))
			retfn <- function(x) { x <- data.frame(x); rownames(x) <- .rn; x; }
		else
			retfn <- function(x) { x <- data.matrix(data.frame(x)); rownames(x) <- .rn; x; }
	}

	#scale values to min~max range, or to part of the median/center quantile (qmid)
	q.scale <- function(v, qmin=min(v), qmax=max(v), qmid=1) {
		qs <- quantile(v, p=c((1-qmid)/2,(qmid+1)/2))
		vmin <- qs[1]
		vmax <- qs[2]
		#3-rule eq: (v - vmin) / (vmax - vmin) = (q - qmin) / (qmax - qmin)
		((v - vmin) / (vmax - vmin)) * (qmax - qmin) + qmin
	}

	#median centering
	if (center) {
		if (verbose) cat("Centering...\n")
		# center to 
		md <- unlist(lapply(l.mx, median, na.rm = TRUE))
		mmd <- median(md)
		l.mx <- lapply(l.mx, function(l) l-(median(l,na.rm=TRUE) - mmd))
	}

	#shrinking distributions to median/center quantile range
	if (!is.null(qshrink)) {
		if (verbose) cat("Shrinking...\n")
		minq <- median(unlist(lapply(l.mx, quantile, p=(1-qshrink)/2, na.rm = TRUE)))
		maxq <- median(unlist(lapply(l.mx, quantile, p=(qshrink+1)/2, na.rm = TRUE)))
		l.mx <- lapply(l.mx, q.scale, qmin=minq, qmax=maxq, qmid=qshrink)
	}

	#maximum number of values
	m <- 0 
	for (i in 1:length(l.mx)) m <- max(m,sum(is.finite(l.mx[[i]])))
	if (verbose) cat("Max-Quantiles:",m,"\n")

	#quantiling
	q.mx <- matrix(NA,nrow=m, ncol=length(l.mx))
	for (i in 1:length(l.mx)) { q.mx[,i] <- q.func(l.mx[[i]],probs=(1:m)/m,na.rm=TRUE); if (verbose) cat("Quantiled:",i,"\n"); }

	#averaging
	avg <- numeric(length=m)
	if (verbose) cat("Averaging...\n")
	for (i in 1:m) avg[i] <- m.func(q.mx[i,], na.rm=TRUE) 

	#scale averages?
	if (min(avg) != qmin  ||  max(avg) != qmax  ||  qmid < 1) avg <- q.scale(avg, qmin=qmin, qmax=qmax, qmid=qmid)

	#normalizing
	q.mx <- list()
	if (verbose) cat("Normalizing...\n")
	for (i in 1:length(l.mx)) { 
		if (verbose) cat("Computing ecdf",i,"...\n");
		avg.cdf <- ecdf(l.mx[[i]][is.finite(l.mx[[i]])])
		r <- pmin(pmax(1,avg.cdf(l.mx[[i]])*m),m)
		xq <- avg[r]
		names(xq) <- names(l.mx[[i]])
		q.mx[[i]] <- xq
	}
	names(q.mx) <- names(l.mx)

	#return
	if (verbose) cat("Done!\n");
	if (averages) avg else retfn(q.mx)
}


pad <- function(chv, char=if(is.numeric(chv)) "0" else " ", pos=-1, len=max(nchar(chv))) {
	sp <- paste(rep(char,len),collapse="")
	chv <- as.character(chv)
	nc <- nchar(chv)
	tc <- len-nc
	if (pos < 0) paste(substr(rep(sp,length(chv)),1,tc*nchar(char)),chv,sep="")
	else if (pos == 0) paste(substr(rep(sp,length(chv)),1,(tc-trunc(tc/2))*nchar(char)), chv, substr(rep(sp,length(chv)),1,trunc(tc/2)*nchar(char)), sep="")
	else  paste(chv, substr(rep(sp,length(chv)),1,tc*nchar(char)), sep="")
}

#trim <- function(ch) sub(" +$","",sub("^ +","",ch))
#rtrim <- function(ch) sub(" +$","",ch)
#ltrim <- function(ch) sub("^ +","",ch)
trim <- function(ch,char=" ") sub(paste(char,"+$",sep=""),"",sub(paste("^",char,"+",sep=""),"",ch))
rtrim <- function(ch,char=" ") sub(paste(char,"+$",sep=""),"",ch)
ltrim <- function(ch,char=" ") sub(paste("^",char,"+",sep=""),"",ch)






plot.pairs <- function (x=1:1, y=1:1, func=function(...) plot(..., xlab="", ylab="", axes=FALSE), gap=1, ...) {
	dots <- list(...)
	nmdots <- names(dots)
	oma <- if ("oma" %in% nmdots) 
		dots$oma
	else NULL
	main <- if ("main" %in% nmdots) 
		dots$main
	else NULL
	if (is.null(oma)) {
		oma <- c(4, 4, 4, 4)
		if (!is.null(main)) 
			oma[3] <- 6
	}
	if (length(x) == 1) x <- 1:x
	if (length(y) == 1) y <- 1:y
	opar <- par(mfrow = c(length(y), length(x)), mar = rep.int(gap/2, 4), oma = oma)
	on.exit(par(opar))
	for (i in y) {
		for (j in x) {
			func(j, i, ...)
		}
	}
	if (!is.null(main)) {
		font.main <- if ("font.main" %in% nmdots) 
			dots$font.main
		else par("font.main")
		cex.main <- if ("cex.main" %in% nmdots) 
			dots$cex.main
		else par("cex.main")
		mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
	}
	invisible(NULL)
}








plot.confront.one.pair <- function(l, x.adj=0, col.adj=0, test=t.test, dec=8, ...) {

	for (i in 1:length(l)) {
		x <- l[[i]]
		points(rep(i, length(x))+x.adj, x, col=i+col.adj, ...)
		axis(side=1, at=i+x.adj, names(l)[i], col=i+col.adj, col.axis=i+col.adj, lwd=3)
	}
	for (i in 2:length(l)) {
		x <- l[[i-1]]
		y <- l[[i]]
		xcol <- ifelse(x > y, i-1, ifelse(x == y, 8, i))+col.adj
		xlty <- ifelse(x > y, i-1, ifelse(x == y, 1, i))
		for (j in 1:length(x)) lines(c(i-1, i)+x.adj, c(x[j], y[j]), col=xcol[j], lty=xlty[j],...)
		r <- test(x, y)
		text((i-1+i)/2+x.adj, par("usr")[4], paste(paste(names(r), " = ", round(unlist(r),dec), sep=""), collapse="\n"), adj=c(.5,1))
	}
}
#plot.confront.pair(list("Normal"=normal.q[slit2,], "Tumour"=tumour.q[slit2,]), pch=20, cex=2)
#plot.confront.pair(list("Normal"=normal.q[slit2,], "Tumour"=tumour.q[slit2,]), list("Normal"=normal.q[slit2,], "Tumour"=tumour.q[slit2,]), pch=20, cex=2)
#this draw points in vertical where in horizontal is each sample
plot.confront.pair <- function(l1, l2, l3, l4, main="", xlab="Group", ylab="Value", 
	test=function(a, b) { list(t.test=t.test(a,b,paired=TRUE)$p.value, wilcoxon.test=wilcox.test(a,b,paired=TRUE)$p.value) }, xlim=NULL, ylim=NULL,
	x.adj=0, col.adj=0, ...) {
	xl <- l1
	if (!missing(l2)) xl <- c(xl, l2)
	if (!missing(l3)) xl <- c(xl, l3)
	if (!missing(l4)) xl <- c(xl, l4)
	if (missing(ylim)) ylim <- range(unlist(xl))*c(0.9,1.1)
	if (missing(xlim)) xlim = c(0, length(xl)+1)
	plot(0,0,type="n", main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, xaxt="n", ...)
	plot.confront.one.pair(l1, test=test, ...)
	if (!missing(l2)) plot.confront.one.pair(l2, test=test, x.adj=x.adj+length(l1), col.adj=col.adj, ...)
	if (!missing(l3)) plot.confront.one.pair(l3, test=test, x.adj=x.adj+length(l1)+length(l2), col.adj=col.adj, ...)
	if (!missing(l4)) plot.confront.one.pair(l4, test=test, x.adj=x.adj+length(l1)+length(l2)+length(l3), col.adj=col.adj, ...)
}


# x is a list
contingency <- function(x) { 
	xl <- unique(unlist(lapply(x, function(z) levels(as.factor(z))))) # xl <- levels(factor(unlist(x)))
	xn <- lapply(x, function(y) {
		yt <- table(y)
		z <- numeric(length(xl))
		names(z) <- xl
		z[names(yt)] <- yt
		z
	})
	data.matrix(data.frame(xn))
}

#rows are variables, cols are samples, cls is class for samples
#cols are variables, rows are samples, cls is class for samples
plot.confront <- function(dframe, cls, main="", xlab="Group", ylab="Value", 
	test=function(x,v,k) {  if (is.factor(v) || is.character(v)) unlist(list(kruskal.test=fisher.test(table(v,k))$p.value)) else unlist(list(kruskal.test=kruskal.test(x)$p.value)) },
	x.adj=0, col.adj=0, test.label="Kruskal-Wallis\nOr Fisher-Exact", alfa=0.05, 
	sort=TRUE, las=1, ncol=max(dfapply(cls, 2, nlevels)), correct=function(p) p.adjust(p, "fdr"), 
	left.margin=8, 
	legends=TRUE,
	legend.pos="top",
	legend.cols=ncol, 
	legend.cex=1,
	same.range=FALSE,
	filter=FALSE,
	rows.range=c(0,1),
	pch=1,
	cex=1,
	cex.points=1,
	adjust.cex=missing(cex.points), jitter=FALSE,
	outliers=0,
	col=1:30,
	...) {
	ocls <- cls
	if (is.vector(cls) || is.factor(cls)) {
		if (length(cls) != nrow(dframe)) stop("Class vector should be the same than the number of samples (rows).")
		cls <- data.frame(as.factor(cls))
		while (ncol(cls) < ncol(dframe)) {
			cls <- cbind(cls, cls[,1])
		}
	}
	multifactor <- 1*(!is.factor(ocls) && !is.vector(ocls))
	nlvl <- max(dfapply(cls, 2, nlevels)) #nlevels(cls)
	if (nrow(cls) != nrow(dframe)) stop("Class vector/data.frame should be the same than the number of samples (rows).")
	if (ncol(cls) != ncol(dframe)) stop("Class data.frame should be the same dimension than input (or a vector/factor length=rows data).")
	#plot(0,0,type="n", main=main, xlab=xlab, ylab=ylab, xlim=c(0,(nlvl+2)*nrow(dframe)), ylim=c(0,1), yaxt="n", xaxt="n", ...)
	pis <- c()
	#per.class <- table(cls)

	okdata <- function(x) {
		wok <- !is.na(x)
		L <- length(x)
		if (outliers > 0) {
			library(outliers)
			nrem <- ifelse(outliers < 1, round(sum(wok)*outliers), outliers)
			while (L-sum(wok) < nrem) {
				wo <- outlier(x[wok], logical=TRUE)
				if (any(wo)) { wok[which(wok)[wo]] <- FALSE }
				else { break }
			}
		}
		which(wok)
	}
	ndata <- rep(nrow(dframe), ncol(dframe))
	for (i in 1:ncol(dframe)) {
		wok <- okdata(dframe[,i])
		x <- dframe[wok,i]
		kcls <- cls[wok,i]
		vl <- lapply(levels(kcls), function(x) dframe[kcls==x,i])
		names(vl) <- levels(kcls)
		nl <- nlevels(kcls)
		p <- NULL
		try(p <- test(vl, x, kcls))
		pis[i] <- if (nl == 1 || is.null(p) || is.na(p) || !is.finite(p)) 1 else p
		ndata[i] <- length(wok)
	}
	names(pis) <- colnames(dframe)
	xpis <- pis
	if (sort) {
		dframe <- dframe[,order(pis),drop=FALSE]
		cls <- cls[,order(pis),drop=FALSE]
		pis <- pis[order(pis)]
	}
	qus <- NULL
	if (is.function(correct)) qus <- correct(pis)
	if (filter) {
		f <- pis <= alfa ### (if(is.null(qus)) pis else qus) <= alfa
		if (sum(f) > 0) {
			dframe <- dframe[,f,drop=FALSE]
			pis <- pis[f]
			qus <- qus[f]
		} else {
			warning("Filter produced no results. Relax your conditions. No filtering took place.")
		}
	}
	xpar <- par("mar")
	on.exit(par(mar=xpar), add=TRUE)
	np <- xpar
	np[2] <- left.margin
	if (multifactor) {
		np[3] <- 2
		np[1] <- 2.5
	}
	par(mar=np)
	plot.new()
	xr <- trunc(pmin(ncol(dframe),pmax(1,ncol(dframe)*rows.range))+0.5)
	if (xr[1] != 1  || xr[2] != ncol(dframe)) {
		w <- xr[1]:xr[2]
		dframe <- dframe[,w,drop=FALSE]
		pis <- pis[w]
		qus <- qus[w]
	}
	ncols <- xr[2]-xr[1]+1
	plot.window(c(-.2,1.55+ifelse(is.null(qus),0,0.3)), c((nlvl+1)*ncols+2*0.5,-1*0-nlvl+1), "", NA, ...)
	if (adjust.cex) {
		cex.points <- adjustCEXuser(pch, height=1, maxCex=1.5)
	}
	colfunc <- function(p,a) { return (ifelse(p <= a, "#770000", ifelse(p <= a*2, "#335577", "#888888"))) }
	xcex <- adjustCEXuser("X", height=2.5/nlvl)
	for (i in 1:ncol(dframe)) {
		wok <- okdata(dframe[,i])
		x <- dframe[wok,i]
		kcls <- cls[wok,i]
		per.class <- table(kcls)
		rect(0, (i-1)*(nlvl+1), 1, i*(nlvl+1), border=8)
		next.pos <- 0
		xcent <- (i-0.5)*(nlvl+1)
		if (is.character(x)) x <- factor(ifelse(nchar(x)==0,".",x))
		if (is.factor(x) || length(unique(x)) < 5 && any(table(x) > 5)) {
			if (!is.factor(x)) x <- factor(x)
			xtab <- table(x, kcls)
			xlvl <- levels(x)
			x.nlvl <- length(xlvl)
			mxbs <- 1/x.nlvl
			for (j in 1:x.nlvl) text((j-1)*mxbs*1.222, (i-1)*(nlvl+1)+nlvl/2, xlvl[j])
			for (k in 1:nlvl) {
				xlabp <- (k-1)*mxbs
				xv <- x[kcls==levels(kcls)[k]]
				xav <- numeric(x.nlvl)
				names(xav) <- xlvl
				xt <- table(xv)
				xav[names(xt)] <- xt
				ypos <- (i-1)*(nlvl+1)+k
				for (j in 1:x.nlvl) {
					xpos <- (j-1)*mxbs
					if (k==1) { 
						lines(c(xpos, xpos), c(i-1,i)*(nlvl+1), col=8)
						#lines(c(xpos, xpos)+mxbs*0.75, c(i-1,i)*(nlvl+1), col=8, lty=2)
						#text(xpos, (i-1)*(nlvl+1), paste0(xlvl[j],"=",sum(xtab[j,])), adj=c(0,1))
					}
					#+(xav[j]/per.class[k])
					rect(xpos+mxbs*0.333, ypos-0.666/nlvl, xpos+mxbs*0.333+(xav[j]/max(xtab*1.1))*mxbs*0.666, ypos+0.666/nlvl, border=k, col=k)
					text(xpos+mxbs*0.333, ypos, xav[j], adj=c(1,0.5), cex=xcex, ...)
					#paste0("(",round(100*xav[j]/sum(xtab[j,]),0),"%) ",xav[j])
					xlabperc <- paste0(round(100*xav[j]/sum(xtab[j,]),0),"% ")
					
					#if (x.nlvl < 3) {
						text((j-1)/x.nlvl, ypos, xlabperc, col=col[k], adj=c(0,0.5), cex=cex*0.75*0+xcex)
					#} else {
					#	#text((j-1)/x.nlvl+strwidth("00% ",cex=cex*0.75)*(k-1)+0.01, (i)*(nlvl+1), xlabperc, col=k, adj=c(0,0), cex=cex*0.75)
					#	text((j-1)/x.nlvl+strwidth("00% ",cex=xcex)*(k-1)+0.01, (i)*(nlvl+1), xlabperc, col=k, adj=c(0,0), cex=xcex)
					#}
				}
			}
		} else {
			knlvl <- nlevels(kcls)
			knval <- length(unique(x))
			xo <- x
			jitterea <- (jitter && (knval < length(x)/10 || any(table(x) > 3)))
			##if (jitter && knlvl < length(x)/10) x <- jitter(x)
			if (!same.range || i==1) r <- range(if(same.range) dframe else xo, na.rm=TRUE)
			dr <- r[2]-r[1]
			rcalc <- r
			if (jitterea) { rcalc[1] <- rcalc[1]-dr*0.05; rcalc[2] <- rcalc[2]+dr*0.05}
			dr <- rcalc[2]-rcalc[1]
			for (k in 1:knlvl) {
				xv <- x[kcls==levels(kcls)[k]]
				v <- (xv-rcalc[1])/dr
				ypos <- (i-1)*(nlvl+1)+k*nlvl/knlvl
				vok <- v
				if (jitterea) vok <- jitter(vok, factor=0.5)
				yposv <- rep(ypos,length(v))
				if (jitterea) yposv <- jitter(yposv, factor=0.25)
				points(vok, yposv, col=col[as.numeric(k)], cex=cex.points, pch=pch, ...)
				if (multifactor) {
					text(1, ypos, paste0("  ",length(v), " ",levels(kcls)[k]), col=col[as.numeric(k)], adj=c(0,0.5), cex=cex.points, ...)
				}
				if (any(is.na(x))) {
					text(1, ypos, sum(!is.na(xv)), adj=c(-0.25,0.5), col=col[as.numeric(k)], cex=cex*0.75*0+xcex)
				}
			}
			text(-0.01, xcent, round(r[1],5-log10(dr)), adj=c(1,.5), cex=cex, ...)
			text(1.2+multifactor*.1,   xcent, paste0(round(r[2],5-log10(dr)),ifelse(jitterea," (j)","")), adj=c(1,.5), cex=cex, ...)
		}
		text(1.5+multifactor*.1,   xcent, round(pis[i], 7), adj=c(1,.5), col=colfunc(pis[i],alfa), cex=cex, ...)
		if (!is.null(qus)) text(1.8+multifactor*.1,   xcent, round(qus[i], 7), adj=c(1,.5), col=colfunc(qus[i],alfa), cex=cex, ...)
	}
	if (length(legends) > 1 || (!is.character(legends) && legends)) {
		if (! multifactor) {
			cls <- cls[,1]
			plot.legend(legend.pos, paste(levels(cls), " (",table(cls),")",sep=""), fill=col[1:nlvl], ncol=legend.cols, cex=legend.cex, bty="n", col=col)
		} else {
			axis(side=1, at=1.05, labels="Levels")
		}
	}
	axis(side=2, at=((1:ncols-0.5)*(nlvl+1)), labels=paste(xr[1]:xr[2],": ",colnames(dframe), ifelse(ndata != nrow(dframe), paste0(" (n=",ndata,")"), ""), sep="")[xr[1]:xr[2]], las=las, ...)
	#axis(side=1, at=na.omit(c(-.1, 1.1, 1.4,ifelse(is.null(qus),NA,1.7))), labels=na.omit(c("Min","Max", test.label,ifelse(is.null(qus),NA,"q.value"))), las=las, ...)
	axis(side=1, at=na.omit(c(-.1, multifactor*.1+1.15, multifactor*.1+1.4,ifelse(is.null(qus),NA,multifactor*.1+1.7))), labels=na.omit(c("Min","Max", "p.value", ifelse(is.null(qus),NA,"q.value"))), las=las, ...)
	axis(side=1, at=multifactor*.1+1.4, labels=test.label, las=las, padj=-1.25, ...)
	title(main=main)
	invisible(xpis)
}
#plot.confront(data, cls)


options("serviceUrl.entrez" = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")

#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# query NCBI Pubmed to get the pmid (PubMed identifier)
########################################################################
# Example:
#  medline.query("&term=cancer+OR+diabetes")
#  medline.query("&term=Trevino+V[Author]")
#  
options("serviceUrl.entrez" = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
medline.query <- function (..., program="esearch.fcgi", db="pubmed", parameter="&rettype=count", collapse="+AND+", verbose=FALSE, baseUrl=getOption("serviceUrl.entrez"), pause=FALSE) {
	require (XML) || cat("Need the XML package from CRAN!")
	all.terms <- unlist(lapply(unlist(list(...)),trim))
	all.terms <- all.terms[nchar(all.terms) > 0]
	# Get the query string ready. This string should be in the
	# Pubmed syntax. The Pubmed syntax is documented at
	# http://eutils.ncbi.nlm.nih.gov/entrez/query/static/esearch_help.html
	query <- paste (baseUrl,program,"?","db=",db,parameter,paste(all.terms,collapse=collapse),sep="")
	if (verbose) { cat("FETCHING URL [",query,"]"); flush.console(); }
	xtext <- try(readLines(query))
	xml <- try (xmlTreeParse(xtext, asText=TRUE))   # in XML
	if (verbose) { cat("..."); flush.console(); }
	if (pause) medline.pauseBetweenQueries();
	if (verbose) { cat("\n"); flush.console(); }
	return (xml)
}

#based on MedlineR_v30.txt (MedlineR "package")
#examples:
#medline.count("SSX1")
#medline.count("SSX*")
#medline.count("SSX*[Sym]", verbose=TRUE)
#medline.count("Trevino+V[Author]") #5
#medline.count("Trevino+V") #5
medline.count <- function (...) {
	xml <- medline.query(..., parameter="&rettype=count&term=")
	return (as.numeric(xmlValue(xmlRoot(xml)[["Count"]])))
}


#based on MedlineR_v30.txt (MedlineR "package")
#examples:
#medline.pmids("16510496")
#[[1]]
#An object of class 'pubMedAbst':
#Title: GALGO: an R package for multivariate variable selection using genetic algorithms.
#PMID: 16510496
#Authors: V Trevino, F Falciani
#Journal: Bioinformatics
#Date: May 2006

#medline.abstracts(medline.pmids("Trevino+V[Author]"))
#
#lapply(medline.pmids(medline.abstracts("Trevino+V[Author]")), articleTitle)
#[1] "Synovial fluid leukocyte apoptosis is inhibited in patients with very early rheumatoid arthritis."			   
#[2] "GALGO: an R package for multivariate variable selection using genetic algorithms."							   
#[3] "Analysis of host response to bacterial infection using error model based gene expression microarray experiments."
#[4] "Making sense of molecular signatures in the immune system."													  
#[5] "Transcription unit conservation in the three domains of life: a perspective from Escherichia coli."			  
#medline.abstracts(c(123456,765432))

medline.abstracts <- function (..., abstracts=TRUE) {
	xml <- medline.query(..., program="efetch.fcgi", parameter="&retmode=xml&id=", collapse=",")
	if (abstracts) {
		#based on pubMedAbst class example in annotate package
		library(annotate)
		a <- xmlRoot(xml)
		n <- length(xmlChildren(a)) 
		xabs <- list()
		for (i in 1:n) xabs[[i]] <- buildPubMedAbst(a[[i]])
		return (xabs)
	}
	return (xml)
}

#based on MedlineR_v30.txt (MedlineR "package")
#medline.abstracts(16510496)
#medline.abstracts(medline.pmids("Trevino+V[Author]"))
#		Id		 Id		 Id		 Id		 Id 
#"16859518" "16510496" "15800204" "15134529" "11275307"
medline.pmids <- function (..., retmax=NULL) {
	xml <- medline.query(..., parameter=paste(if (is.null(retmax)) "" else paste("&retmax=",retmax,sep=""),"&rettype=uilist&term=",sep=""))
	return (sapply(xmlChildren(xmlRoot(xml)[["IdList"]]), xmlValue))
	#return (xml)
}

#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# pause between queries, according to the NCBI rule
########################################################################
# According to the NCBI rule, query larger than 100 requests should be 'nicely' run off hours and weekends
# The NCBI rule is documented at
#  http://eutils.ncbi.nlm.nih.gov/entrez/query/static/eutils_help.html
#Examples:
#medline.pauseBetweenQueries()
#medline.pauseBetweenQueries(30, 5)
medline.pauseBetweenQueries <- function (peak=15, offpeak=3) { # pause (in seconds) during peak and off-peak hours
	# Date example:
	# "Thu"	  "Jan"	  "15"	   "16:46:11" "2004"
	result.date<- unlist (strsplit(date(), split=" "))
	hour <- as.numeric(unlist (strsplit (result.date[4], split=':'))[1])
	# off peak hours are Sat, Sun or anytime between 9 pm and 5 am
	Sys.sleep ( if ( (result.date[1] %in% c("Sat","Sun")) | !(hour %in% 5:21) ) peak else offpeak )
}


#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# return a co-occurrence matrix for a list of terms
########################################################################
# Example:
#  find the co-occurance matrix of these terms in yeast
#  termList<- c("STE18", "DIG1", "STE12", "SST2", "HOG1",
#			   "pheromone", "cell+cycle")
#  result.matrix<- getAmatrix(termList,
#	termAdditional= "+AND+(Saccharomyces+OR+yeast")
#  # after building the matrix you can get an image for visualization
#  image (result.matrix)
#  image (result.matrix^0.2, col=terrain.colors(16))
# Note:
#  must transform the raw data matrix for display,
#  otherwise, the dynamic range of the color is not
#  suitable for human recognition!
medline.coocurrance <- function(var.terms, fix.terms="") {
	q.terms <- unlist(lapply(unlist(var.terms),trim))
	f.terms <- unlist(lapply(unlist(fix.terms),trim))
	# initialize the co-occurance matrix
	n.t <- length(q.terms)
	r.m <- matrix (0, ncol=n.t, nrow=n.t)
	for (i in 1:n.t) {
		for (j in i:n.t) {
			n <- medline.count(q.terms[i], q.terms[j], f.terms)
			medline.pauseBetweenQueries()
			result.matrix[i,j] <- result.matrix[j,i] <- n
		}
	}
   return (r.m)
}


#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# output the co-occurrence matrix in Pajek format for visualization
########################################################################
# The output is a graph, where each node is a term connected to other
# terms by the co-occurance matrix.  The Pajek software should be used
# to visualize the output.
# documentation for the Pajek file format can be found at
# http://vlado.fmf.uni-lj.si/pub/networks/pajek/
# Note:
#   Pajek file name must end with .net It is a requirement for Pajeck to recognize the file type.
# Example:
#  termList<- c("STE18", "DIG1", "STE12", "SST2", "HOG1",
#			   "pheromone", "cell+cycle")
#  result.matrix<- getAmatrix(termList,
#	termAdditional= "+AND+(Saccharomyces+OR+yeast")
#  write.pajek (result.matrix, termList, fileName="test.net")
#
medline.write.pajek<- function (result.matrix,	   # co-occurance matrix
						termList=NULL,	   # node name
						nodeColorList=NULL,  # node color, default
											 #   is yellow
						fileName){		   # must end with .net

  # Pajek files must end with .net
  # fileName<- "genes.net"

  # setup color list
  n.vertices<- dim (result.matrix)[1]
  if (is.null(nodeColorList)){			   # default: yellow
	nodeColorList<- I(rep ("Yellow", n.vertices))
  }										  # see help (I)

  # setup term list, default is numerical value
  if (is.null(termList)){
	termList<- 1:n.vertices
  }

  # builds data frame
  df.vertices<- data.frame (
	id= 1:n.vertices,
	name= termList,
	color= I(nodeColorList)
  )

  zz <- file(fileName, "w")			 # open an output file connection

  # write vertices to Pajek-formatted file
  cat ("*Vertices ", n.vertices,		# nodes ("Vertices")
	   "\n", sep="", file=zz)
  for (i in 1: n.vertices) {
   cat (df.vertices$id[i], " \"",
		as.character(df.vertices$name[i]), "\"",		 # name
		" x_fact ", round (result.matrix[i,i]^0.2)+1,	# size
		" y_fact ", round (result.matrix[i,i]^0.2)+1,	#
	  # " s_size ", round (result.matrix[i,i]^0.2)+1,	#
		" ic ", as.character(df.vertices$color[i]),	  # color
		"\n",
		sep="", file=zz)
  }

  # write edges to Pajek-formatted file
  cat ("*Edges", "\n", sep="", file=zz) # conncetions ("Edges")

  for (i in 1:(n.vertices-1)){
   if (result.matrix [i,i]==0) {next}
   for (j in (i+1):n.vertices) {
	   n.counts <-result.matrix[i,j]
	   if (n.counts==0) {next}
	   thickness.line<- n.counts^0.5   # transformation for visual effect
		cat (
		 i, " ",					   # from
		 j, " ",					   # to
		 round (thickness.line)+1,	 # width
		 "\n", sep="", file=zz)
	}
  }
  close(zz)
}


#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# mapping a term to MeSH term
########################################################################
# Example:
#   mapToMeSH (term="nosebleed")
medline.mesh <- function (...) {
	xml <- medline.query(..., db="mesh", parameter="&term=")
	#query <- paste (baseUrl,"esearch.fcgi?","db=mesh&term=",term,sep="")
	#parse resulting XML into a tree to be returned to user
	#result.xml <- try(xmlRoot(xmlTreeParse(file=query, isURL=T)))
	ids <- sapply(xmlChildren (xmlRoot(xml)[["IdList"]]), xmlValue)
	#step 2: retrieve the info according to the IDs
	#query <- paste(baseUrl,"esummary.fcgi?","db=mesh&id=",paste (ids, collapse=","),sep="")
	xml <- medline.query(ids, program="esummary.fcgi", db="mesh", parameter="&id=", collapse=",")
	# print (query)
	#result.xml<- try (xmlTreeParse(file=query, isURL=T))
	return (xml)
}

#based on MedlineR_v30.txt (MedlineR "package")
########################################################################
# return a list of stop words
########################################################################
medline.stopwords<- function () {
 return (c(
"a", "about", "again", "all", "almost", "also", "although", "always", "among", "an", "and", "another", "any", "are", "as", "at", "be", "because", "been", "before", "being", "between", "both", "but", "by", "can", "could", "did", "do", "does", "done", "due", "during", "each", "either", "enough", "especially", "etc", "for", "found", "from", "further", "had", "has", "have", "having", "here", "how", "however", "i", "if", "in", "into", "is", "it", "its", "itself", "just", "kg", "km", "made", "mainly", "make", "may", "mg", "might", "ml", "mm", "most", "mostly", "must", "nearly", "neither", "no", "nor", "obtained", "of", "often", "on", "our", "overall", "perhaps", "quite",
 "rather", "really", "regarding", "seem", "seen", "several", "should", "show", "showed", "shown", "shows", "significantly", "since", "so", "some", "such", "than", "that", "the", "their", "theirs", "them", "then", "there", "therefore", "these", "they", "this", "those", "through", "thus", "to", "upon", "use", "used", "using", "various", "very", "was", "we", "were", "what", "when", "which", "while", "with", "within", "without", "would"
))
}


#
#x <- medline.apply(c("MAGE*","TKTL*","DAZL*","SSX*","SSX4","MAEL*"), fix.terms="methylation", verbose=TRUE)
#MAGE* TKTL* DAZL*  SSX*  SSX4 MAEL* 
#  113	 0	 1	 7	 1	 5
medline.apply <- function(var.terms, func=medline.count, fix.terms="", pause.by=if (length(var.terms) > 99) 10 else 20, ...) {
	n.t <- length(var.terms)
	L <- list()
	for (i in 1:n.t) {
		o <- func(var.terms[i], fix.terms, ...)
		L[[i]] <- o
		if (i %% pause.by == 0) medline.pauseBetweenQueries()
	}
	try(names(L) <- if(is.null(names(var.terms))) as.character(var.terms) else names(var.terms), silent=TRUE)
	medline.pauseBetweenQueries()
	return (L)
}


medline.simplify.geneName <- function(x) {
	y <- sub("[^A-Z0-9]+.*","*",x)
	yn <- nchar(y)
	z <- gsub("[^A-Z0-9]","",x)
	ifelse(yn > 2, y, z)
}

medline.generalize.geneName <- function(x) {
	y <- sub("[0-9].*","",x)
	y <- sub("[^A-Z].*","",y)
	x <- medline.simplify.geneName(x) #sub("[^A-Z0-9]+.*","*",x)
	yn <- nchar(y)
	ifelse(yn > 2, paste(y,"*",sep=""), x)
}

##medline.generalize.geneName(c("MAGEB2","TKTL1","DAZL","GAGE2","SSX4","SSX4","MAEL","CT45-1","SSX1","TKTL1","FLJ32942","CT45-3","GAGE7B","GAGE7","MAGEA4","FMR1NB","NA","CTAG1B","SPANXA1","HORMAD1","GAGE6","GAGE4","SSX1","MAGEB1","ECAT8","LOC148756","MAGEA12","MAGEA6"))
# [1] "MAGEB*"  "TKTL*"   "DAZL*"   "GAGE*"   "SSX*"	"SSX*"	"MAEL*"  
# [8] "CT45*"   "SSX*"	"TKTL*"   "FLJ*"	"CT45*"   "GAGE*"   "GAGE*"  
#[15] "MAGEA*"  "FMR*"	"NA"	  "CTAG*"   "SPANXA*" "HORMAD*" "GAGE*"  
#[22] "GAGE*"   "SSX*"	"MAGEB*"  "ECAT*"   "LOC*"	"MAGEA*"  "MAGEA*" 
##medline.simplify.geneName(c("MAGEB2","TKTL1","DAZL","GAGE2","SSX4","SSX4","MAEL","CT45-1","SSX1","TKTL1","FLJ32942","CT45-3","GAGE7B","GAGE7","MAGEA4","FMR1NB","NA","CTAG1B","SPANXA1","HORMAD1","GAGE6","GAGE4","SSX1","MAGEB1","ECAT8","LOC148756","MAGEA12","MAGEA6"))
# [1] "MAGEB2"	"TKTL1"	 "DAZL"	  "GAGE2"	 "SSX4"	  "SSX4"	 
# [7] "MAEL"	  "CT45*"	 "SSX1"	  "TKTL1"	 "FLJ32942"  "CT45*"	
#[13] "GAGE7B"	"GAGE7"	 "MAGEA4"	"FMR1NB"	"NA"		"CTAG1B"   
#[19] "SPANXA1"   "HORMAD1"   "GAGE6"	 "GAGE4"	 "SSX1"	  "MAGEB1"   
#[25] "ECAT8"	 "LOC148756" "MAGEA12"   "MAGEA6"   



### to convert names to indexes, equivalent to several "which"
### example: names2index(c("A","B","C"), rownames(data.frame)), will return row numbers
### equivalent to "match" ???
names2index <- function(access.names, index.names, by=10000) {
	if (any(index.names=="")) {
		index.names[index.names==""] <- "@~@"
		access.names[access.names==""] <- "@~@"
	}
	v <- 1:length(index.names)
	names(v) <- index.names
	r <- c()
	i <- 1
	while (i <= length(access.names)) {
	r <- c(r, v[access.names[i:min(length(access.names),(i+by-1))]])
	i <- i + by
	}
	r
}



## converts a list to a data frame adding elements of the list to columns or rows
list.2.data.frame <- function(x, func=cbind, transfunc=function(x) x) {
	for (i in 1:length(x)) {
		if (i==1) xdf <- data.frame(x[[i]])
		else xdf <- func(xdf, transfunc(x[[i]]))
	}
	if (length(names(x)) == ncol(xdf)) colnames(xdf) <- names(x)
	else if (length(names(x)) == nrow(xdf)) rownames(xdf) <- names(x)
	xdf
}

## converts the elements of a list to a matrix, assuming all elements in the list are same length vectors
list.2.matrix <- function(x, name=NULL) {
	m <- matrix(, ncol=length(x), nrow=if (is.null(name)) length(x[[1]]) else length(x[[1]][[name]]))
	for (i in 1:length(x)) 
		if (is.null(name)) 
			m[,i] <- x[[i]]
		else
			m[,i] <- x[[i]][[name]]
	if (length(names(x)) == ncol(m)) colnames(m) <- names(x)
	m
}


left <- function(x, n) {
	if (n == 0) x
	if (n < 0) x[-((n+1):0+length(x))]
	else left(x, n-length(x))
}

right <- function(x, n) {
	if (n == 0) x
	if (n < 0) x[-(1:-n)]
	else right(x, n-length(x))
}

mid <- function(x, n) {
	if (n == 0) x
	if (n > 0) n <- n - length(x)
	l <- -trunc(abs(n)/2)
	r <- n-l
	x <- left(x, l)
	x <- right(x, r)
	x
}

moda <- function(x) {
	xt <- table(x)
	xn <- names(xt[xt==max(xt)])
	unique(x[which(as.character(x) %in% xn)])
}

#extract the first pattern in each position from a string vector
extract <- function(pattern, stvector, sub=FALSE, ...) {
	re <- regexpr(pattern, stvector, ...)
	xc <- character(len=length(stvector))
	w <- which(re > 0)
	xc[w] <- substr(stvector[w], re[w], re[w]+attr(re,"match.length")[w]-1)
	if (sub) {
		#xc[w] <- sub(pattern, "\\1", xc[w])
		xc[w] <- sub(pattern, paste0("\\",ifelse(is.numeric(sub), sub, 1)), xc[w])
	}
	xc
}

#extract all patterns in each position from a string vector
extract.all <- function(pattern, st, ret.num=FALSE, collapse=NULL, ...) {
	if (length(st) == 0) return (NULL);
	o <- gregexpr(pattern, st, ...)
	if (ret.num) return (unlist(lapply(o, function(x) if (x[1] == -1) 0 else length(x))))
	l <- list()
	for (i in 1:length(o)) {
		p <- o[[i]]
		if (length(p) > 0) {
			l[[i]] <- substr(rep(st[i],length(p)), p, p+attr(p,"match.length")-1)
		} else {
			l[[i]] <- c()
		}
	} 
	if (length(l) == 1) return (l[[1]])
	if (!is.null(collapse)) {
		v <- unlist(lapply(l, paste, collapse=collapse))
		names(v) <- names(o)
		return (v)
	}
	names(l) <- names(o)
	l
}




#violin plot
#vioplot( x, range=1.5, h, ylim, names, horizontal=FALSE,
#   col="magenta", border="black", lty=1, lwd=1, rectCol="black",
#   colMed="white", pchMed=19, at, add=FALSE, wex=1,
#   drawRect=TRUE)
# from vioplot package, Daniel Adler, http://wsopuppenkiste.wiso.uni-goettingen.de/~dadler
# http://addictedtor.free.fr/graphiques/RGraphGallery.php?graph=102
vioplot <- function(x,range=1.5,h=NULL,ylim=NULL,names=NULL, col=1:length(datas), ...)
{
	require(sm)
	# process multiple datas
	datas <- x
	n <- length(datas)
	# pass 1
	#
	# - calculate base range
	# - estimate density
	#
	# setup parameters for density estimation
	upper  <- vector(mode="numeric",length=n)
	lower  <- vector(mode="numeric",length=n) 
	q1	 <- vector(mode="numeric",length=n)
	q3	 <- vector(mode="numeric",length=n)
	med	<- vector(mode="numeric",length=n)
	base   <- vector(mode="list",length=n)
	height <- vector(mode="list",length=n)
	baserange <- c(Inf,-Inf)
	# global args for sm.density function-call   
	args <- list(display="none")
	if (!(is.null(h)))
		args <- c(args, h=h)
			
	for(i in 1:n) {
		data<-datas[[i]]
		# calculate plot parameters
		#   1- and 3-quantile, median, IQR, upper- and lower-adjacent
		data.min <- min(data)
		data.max <- max(data)
		q1[i]<-quantile(data,0.25)
		q3[i]<-quantile(data,0.75)
		med[i]<-median(data)
		iqd <- q3[i]-q1[i]
		upper[i] <- min( q3[i] + range*iqd, data.max )
		lower[i] <- max( q1[i] - range*iqd, data.min )
	   
		#   strategy:
		#	   xmin = min(lower, data.min))
		#	   ymax = max(upper, data.max))
		#
		est.xlim <- c( min(lower[i], data.min), max(upper[i], data.max) ) 
		# estimate density curve
		smout <- do.call("sm.density", c( list(data, xlim=est.xlim), args ) )
		# calculate stretch factor
		#
		#  the plots density heights is defined in range 0.0 ... 0.5 
		#  we scale maximum estimated point to 0.4 per data
		#
		hscale <- 0.4/max(smout$estimate)
		
		# add density curve x,y pair to lists
		base[[i]]   <- smout$eval.points
		height[[i]] <- smout$estimate * hscale
		
		# calculate min,max base ranges
		t <- range(base[[i]])
		baserange[1] <- min(baserange[1],t[1])
		baserange[2] <- max(baserange[2],t[2])
	}
	# pass 2
	#
	# - plot graphics
	# setup parameters for plot
	xlim <- c(0.5, n+0.5)
	if (is.null(ylim)) {
		ylim <- baserange
	}
	if (is.null(names)) {
		label <- c("",1:n,"")
	} else {
		label <- c("",names,"")
	}
	boxwidth <- 0.05
		
	# setup plot
	plot.new()
	plot.window(xlim = xlim, ylim = ylim, ...)
	
	# setup axis
	axis(1,at = c(0:(n+1)), label=label )
	axis(2)
	box()
	# plot data
	for(i in 1:n) {
		# plot left/right density curve
		lines ( i-height[[i]], base[[i]], col=col[i])
		lines ( i+height[[i]], base[[i]], col=col[i] )
		# close density curves
		last <- length(height[[i]])
		lines ( c( i-height[[i]][1] , i+height[[i]][1] ), c( base[[i]][1] , base[[i]][1] ), col=col[i])
		lines ( c( i-height[[i]][last] , i+height[[i]][last] ), c( base[[i]][last] , base[[i]][last] ), col=col[i] )
		# plot 50% KI box
		rect( i-boxwidth/2, q1[i], i+boxwidth/2, q3[i],col=col[i])
		# plot IQR
		lines( c( i, i), c(lower[i], upper[i]), col=col[i] )
		# plot median point
		points( i, med[i], pch=1, col=col[i])
	}
	invisible (list( upper=upper, lower=lower, median=med, q1=q1, q3=q3))
}

gradient.colours <- function(vrange, ncolours = 16, breaks, c.breaks) {
	# c.breaks - list or colours, rgb format, each element contains 3 values = rgb
	# breaks - vector of breaks within vrange
	# example green-black-red in a scale -2 to 2
	# c.breaks = list(c(0,1,0),c(0,0,0),c(1,0,0))
	# breaks = 0 # which assumes that first "break" is -2 and last break is 2
	if (length(vrange) != 2) vrange = range(vrange)
	delta <- (max(vrange) - min(vrange)) / (ncolours-0)
	if (breaks[1] != min(vrange)) breaks <- c(min(vrange),breaks)
	if (breaks[length(breaks)] != max(vrange)) breaks <- c(breaks, max(vrange))
	xcol <- rgb(c.breaks[[1]][1],c.breaks[[1]][2],c.breaks[[1]][3])
	for (i in 2:length(breaks)) {
		d <- breaks[i] - breaks[1]
		n <- round(d / delta) - length(xcol)
		r <- c.breaks[[i]][1] - c.breaks[[i-1]][1]
		g <- c.breaks[[i]][2] - c.breaks[[i-1]][2]
		b <- c.breaks[[i]][3] - c.breaks[[i-1]][3]
		if (n > 1) {
			f <- 1:(n-1)/n
			xcol <- c(xcol, rgb(c.breaks[[i-1]][1] + f * r, c.breaks[[i-1]][2] + f * g, c.breaks[[i-1]][3] + f * b))
		}
		xcol <- c(xcol, rgb(c.breaks[[i]][1],c.breaks[[i]][2],c.breaks[[i]][3]))
	}
	xcol
}

#gradient.colours(c(-2,2), ncolours=16, c.breaks = list(c(0,1,0),c(0,0,0),c(1,0,0)), breaks = 0)

# pattern - pattern to search
# sts - vector of strings
# nchars - total number of rounding chars
# left.chars - total number of chars in left part of the pattern
# ends - chars to put at the end and begining of lines to indicate that there are more chars
# position - indicates if position should be included
# filter - TRUE : non matched string are filtered, FALSE : all strings are included and blank lines are used for those unmatched
align.text.matches <- function(pattern, sts, nchars=80, left.chars=round((nchars-nchar(pattern))/2,0), marks="[]", ends="...", position=FALSE, filter=TRUE, ret.all=FALSE, ignore.case = TRUE, names=NULL) {
	p <- gregexpr(pattern, sts, ignore.case=ignore.case)
	out <- c()
	notends <- gsub("."," ",ends)
	spaces <- "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
	outi <- c()
	outp <- c()
	outmatch <- c()
	if (length(p) > 0)
	for (i in 1:length(p)) {
		pos <- p[[i]]
		if (is.na(pos)  ||  length(pos) == 0 || pos[1] < 0) {
			if (!filter) {
				out <- c(out, substr(spaces, 1, nchars+nchar(ends)*2+ifelse(position, 10, 0)))
			}
		} else {
			len <- attr(pos,"match.length")
			for (j in 1:length(pos)) {
				n <- pos[j]-left.chars
				s <- paste(substring(sts[i], max(1, n), pos[j]-1), substring(marks,1,1), substring(sts[i], pos[j], pos[j]+len[j]-1), substring(marks,2,3), sep="")
				if (pos[j] <= left.chars) {
					s <- paste(substr(spaces, 1, left.chars-pos[j]+1+nchar(ends)), s, sep="")
				} else {
					s <- paste(ends, s, sep="")
				}
				rc <- nchars - nchar(s)
				s <- paste(s, substring(sts[i], pos[j]+len[j], pos[j]+len[j]+rc), sep="")
				if (pos[j]+len[j]+rc < nchar(sts[i])) s <- paste(s, ends, sep="")
				s <- paste(s, substr(spaces, 1, nchars + nchar(ends)*2 - nchar(s) - 2), sep="")
				if (position) s <- paste(pad(i,char=" ",len=4),",",pad(pos[j], char=" ", len=4), ":", s, sep="")
				if (!is.null(names))
					s <- paste(names[i],":",s,sep="")
				outmatch <- c(outmatch, substring(sts[i], pos[j], pos[j]+len[j]-1))
				out <- c(out, s)
				outi <- c(outi, i)
				outp <- c(outp, pos[j])
			}
		}
	}
	if (ret.all) {
		list(text=out, index=outi, pos=outp, match=outmatch)
	} else {
		#attr(out, "index") <- outi
		#attr(out, "pos") <- outp
		out
	}
}
# align.text.matches("factor", "fibroblast growth factor receptor 1 (fms-related tyrosine kinase 2, pfeiffer syndrome)hsa04010:mapk signaling pathway, hsa04520:adherens junction, hsa04810:regulation of actin cytoskeleton, 3d-structure, atp, alternative splicing, atp-binding, autophosphorylation, chromosomal translocation, direct protein sequencing, disease mutation, duplication, extracellular protein, glycoprotein, growth factor receptor, growth factor receptor 1")




# compare two string separating by words
# sep es el separador de palabras
# ret.int indica si se regresan las palabras que hicieron "match"
# uniq indica si se usan palabras no repetidas
# aprox indica si se usa "agrep" para buscar si palabras
text.similitude <- function(a, b, ignore.case=TRUE, sep=" ", ret.int=FALSE, uniq=TRUE, aprox=FALSE, ...) {
	if (ignore.case) { a <- tolower(a); b <- tolower(b) }
	a.words <- unlist(strsplit(a, sep), use.names=FALSE)
	b.words <- unlist(strsplit(b, sep), use.names=FALSE)
	if (a == b) {
		return (1)
	}
	if (uniq) { a.words <- unique(a.words); b.words <- unique(b.words) }
	if (aprox) {
		int <- NULL
		if (sum(nchar(a.words)) < sum(nchar(b.words))) {
			words <- a.words
			zords <- b.words
		} else {
			words <- b.words
			zords <- a.words
		}		
		for (i in 1:length(words)) {
			if (length(agrep(words[i], zords, ...)) > 0) {
				int <- c(int, words[i])
			}
		}
	} else {
		int <- intersect(a.words, b.words)
	}
	if (ret.int) return (int)
	if (length(int) < 1) return (0);
	sum(nchar(int)) / max(1,mean(c(sum(nchar(a.words)), sum(nchar(b.words)))))
}


catf <- function(...) {
	cat(...)
	invisible(flush.console())
}


#nic - number of desired icas or percentage 
#mic - number of estimated icas, twice as nic is recommended to ensure all runs obtain similar icas
#times - number of times fastICA is ran, larger better
#q - correlation cut-off to determine two icas in different runs are the same
#ret.all - return diagnostic values if true, returns only the stablized components otherwise
#Example to observe the correlated icas
#icas <- stabilizedICA(data, times=10, ret.all=TRUE)
#pairs(icas$icas[[1]], gap=0)
#Rx.Sy : x is the RUN and y is the ica
#pairs(icas$icas[[2]], gap=0)
#To estimate correlations of icas
#hist(abs(cor(icas$icas[[1]])),breaks=100)
#To see counts per ica
#table(icas$icas.idx)


stabilizedNMF <- function(k=3, ...) {
	library(NMF)
	nmfs <- function(data, n.comp, method) {
		nm <- nmf((data), rank=n.comp)
		list(S=t(nm@fit@H))
	}
	stabilizedICA(..., compFunc=nmfs)
}


stabilizedISOMAP <- function(k=3, ...) {
	library(vegan)
	.k.stabilized.isomap <<- k
	iso <- function(data, n.comp, method) {
		# +rnorm(nrow(data)*ncol(data),m=0,sd=0.01)
		im <- isomap(dist(data), ndim=n.comp, k=.k.stabilized.isomap, fragmentedOK=TRUE)
		list(S=(im$points))
	}
	stabilizedICA(..., compFunc=iso, q=.90)
}


stabilizedICA <- function(data, nic=5, mic=10, times=100, q=0.95, average=TRUE, ret.all=FALSE, debug=FALSE, compFunc="fastICA") {
	if (debug) cat("nic=",nic,", mic=",mic,", times=",times,", q=",q,", ret.all=",ret.all,"\n")
	if (!is.function(compFunc)) {
		if (compFunc == "fastICA") library(fastICA)
		compFunc <- eval(parse(text=compFunc))
	}
	s <- list()
	icas <- list()
	# Compute icas
	if (debug) cat("Estimating Runs=")
	for (i in 1:times) {
		if (debug) cat(".")
		icas[[i]] <- compFunc(data, mic, method="C")
		s[[i]] <- icas[[i]]$S
	}
	icax <- icas
	
	# Compute correlation matrices
	scm <- list()
	for (i in 1:(times-1)) {
		if (debug) cat("\nCorr. Matrix=",i,paste(dim(s[[i]]),collapse="x"),"...")
		scm[[i]] <- list()
		for (j in (i+1):times) {
			if (debug) cat(paste(dim(s[[j]]),collapse="x"),".")
			scm[[i]][[j]] <- cor(s[[i]], s[[j]])
		}
	}
	
	# Estimate cut-off
	if (debug) {
		cat("\nQuantiles corr matrix=")
		print(quantile(abs(unlist(scm)), 0:100/100))
	}
	xq <- quantile(abs(unlist(scm)), q)
	icas.idx <- list()
	icas <- list()

	# look for equivalent icas
	for (i in 1:(times-1)) {
		if (debug) cat("\nGetting Clusters=",i,"...")
		## i ica run
		sm <- scm[[i]]
		for (j in (i+1):length(sm)) {
			## j ica run comparison
			if (debug) cat(".")
			m <- sm[[j]]
			#cat("i=",i,",li=",length(sm),", j=", j,", is.null(m)=",is.null(m),"\n")
			if (!is.null(m)) {
				aw <- which(abs(m) >= xq, arr.ind=TRUE)
				if (nrow(aw) > 0) {
					# correspondings
					
					# remove possible overlaps
					while (TRUE) {
						tr <- table(aw[,1])
						tc <- table(aw[,2])
						wr <- which(tr > 1) # one row is correlated to more than one column, remove smaller
						wc <- which(tc > 1) # one column is correlated to more than one row, remove smaller
						if (length(wr) > 0) {
							xw <- aw[which(aw[,1] == as.numeric(names(tr)[wr[1]])), ]
							aw <- aw[-order(abs(m)[xw])[1],, drop=FALSE]
						} else if (length(wc) > 0) {
							xw <- aw[which(aw[,2] == as.numeric(names(tc)[wc[1]])), ]
							aw <- aw[-order(abs(m)[xw])[1],, drop=FALSE]
						} else {
							break()
						}
					}
					
					
					for (p in 1:nrow(aw)) {					
						w <- aw[p,2]
						k <- aw[p,1]
						ni <- paste("R",i,".S",k,sep="")
						nj <- paste("R",j,".S",w,sep="")
						idx <- -1
						if (is.null(icas.idx[[ni]])) {
							# save ica for the first time
							icas[[ni]] <- matrix(s[[i]][,k], ncol=1)
							colnames(icas[[ni]]) <- ni
							idx <- which(ni == names(icas))
							icas.idx[[ni]] <- idx
						} else {
							idx <- icas.idx[[ni]]
						}
						if (is.null(icas.idx[[nj]])) {
							# this ica has never been correlated, thus assign the most similar ica
							xn <- colnames(icas[[idx]])
							icas[[idx]] <- cbind(icas[[idx]], s[[j]][,w])
							colnames(icas[[idx]]) <- c(xn, nj)
							icas.idx[[nj]] <- idx
						}
						m[,w] <- 0
						m[k,] <- 0
					}
				}
			}
		}
	}
	
	# get top nic icas
	xt <- table(unlist(icas.idx))
	if (debug) {
		print(xt)
	}
	xo <- order(xt, decreasing=TRUE)
	if (nic < 1) {
		minfreq <- max(round(times * nic,0),2)
	} else {
		minfreq <- xt[xo[nic]]
	}
	w <- which(xt >= minfreq)
	if (length(w) < 1) {
		w <- which(xt >= xt[xo[1]])
	}
	s <- NULL
	for (i in 1:length(w)) {
		s <- cbind(s, if (average) quantile.normalization(icas[[w[i]]])[,1] else icas[[w[i]]][,1])
	}
	
	if (!ret.all) {
		return (list(S=s, frequency=xt[xo], threshold=xq))
	}
	
	return(list(S=s, frequency=xt[xo], IC=icas[xo], icas=icas, icas.idx=unlist(icas.idx), scm=scm, .icas=icax, threshold=xq))

}


# p.information.content divides test that the values encountered within a range for a particular class is not by random chance
# using a hypergeometric distribution
# 
p.information.content <- function(values, classes, wclass, sorted=FALSE, large=TRUE, V=length(values)) {
	if (is.matrix(values)) {
		vt <- t(values)
		V <- ncol(values)		
		return (apply(vt, 2, function(v) p.information.content(v, classes, wclass, FALSE, large, V)) )
	}
	if (! sorted) {
		o <- order(values)
		classes <- classes[o]
		values <- values[o]
	}
	wclass <- which(classes == wclass)
	if (length(wclass) < 2) return (1)
	L <- length(wclass)
	i <- rep(1:(L-1), (L-1):1)
	j <- unlist(sapply(1:(L-1), function(i) (i+1):L))
	p <- p.overlap(V, L, wclass[j]-wclass[i]+1, j-i) # j-i+1 is the lower.tail
	w <- which.min(p)
	p.min <- p[w[1]]
	if (! large) return (p.min)
	w.min <- c(i[w[1]], j[w[1]])
	return (list(p.value = p.min, n=w.min[2]-w.min[1]+1, m=L, positions = wclass[w.min[1]:w.min[2]], size=V, range = values[wclass[w.min]]))
}


#Classification, subtype discovery, and prediction of outcome in pediatric acute lymphoblastic leukemia by gene expression profilingEng-Juh Yeoh 
p.yeohchi2 <- function(values, classes, wclass, sorted=FALSE, V=length(values)) {
	if (is.matrix(values)) {
		vt <- t(values)
		V <- ncol(values)		
		return (apply(vt, 2, function(v) p.yeohchi2(v, classes, wclass, FALSE, V)) )
	}
	if (! sorted) {
		o <- order(values)
		classes <- classes[o]
		values <- values[o]
	}
	cs <- cumsum(classes == wclass)
	d <- diff(cs)
	w <- which(d != 0)+1
	if (cs[1] == 1) w <- c(1, w)
	p <- p.overlap(V, cs[length(cs)], w, cs[w]-1) 
	return (min(p))
}


# does x follow a bimodal distribution?
# x vector
# boxcoxify - transform x by boxcox (which is a monotonical transformation to normalize data)
# http://hosho.ees.hokudai.ac.jp/~kubo/Rdoc/library/car/html/box.cox.powers.html
# fitting the mixture model: http://www.math.mcmaster.ca/peter/mix/mix.html
# package mixdist
p.bimodal <- function(x, boxcoxify = TRUE, breaks="Sturges", constr=mixconstr(), emsteps=3, df=6, ret.all = FALSE, modals=2, ...) {
	if (is.matrix(x) || is.data.frame(x)) return (apply(x, 2, p.bimodal, boxcoxify))
	if (is.list(x)) return (lapply(x,p.bimodal, boxcoxify))
	if (boxcoxify) {
		if (any(x <= 0)) x <- x+abs(min(x))+1 # otherwise box.cox.powers will fail
		library(car)
		bcp <- box.cox.powers(x)
		x <- box.cox(x, bcp$lambda)
	}
	library(mixdist)
	hx <- hist(x, breaks=breaks, plot=FALSE)
	pje <- seq(0,1,length.out=modals*2+1)
	qs <- quantile(x, pje)
	.means <- unlist(sapply(1:modals, function(i) mean(x[x > qs[i*2-1] & x <= qs[i*2+1]])))
	.sds <-   unlist(sapply(1:modals, function(i)   sd(x[x > qs[i*2-1] & x <= qs[i*2+1]])))
	mxp <- data.frame(pi=rep(1/modals,modals), mu=.means, sigma=.sds)
	mx <- mix(data.frame(right=hx$breaks[-1], freq=hx$counts), mixpar=mxp, constr=constr, emsteps = emsteps, ...)
	h0 <- shapiro.test(x)$p.value # probability of being normal
	h1 <- mx$P # a significance level (P-value) for the goodness-of-fit test.
	if (is.na(h1)) {
		h1 <- 1-pchisq(mx$chisq, df=1)
	}
	lambda <- h0 / h1
	L <- -2 * log(lambda)
	if (ret.all) {
		list(mix=mx, h0=h0, h1=h1, lambda=lambda, LR=L, p.value=1 - pchisq(L, df=df))
	} else {
		1 - pchisq(L, df=df)
	}
}

#(1 - pchisq(x$LR0, df))


#differential expression testER - de.test
#ref: selecting genes by test statistics Chen Hua J Biomed Biotech 2005
#ref:Comparison and evaluation of methods for generating differentially expressed gene lists from microarray data - Jeffrey Culhane - BMC Bioinf 2006
#ref: Comparison of various statistical methods for identifying differential gene expression in replicated microarray data Kim Lee Stat Met Med Res 2006
#x - matrix or data frame
#classes
de.test <- function(x, classes, test=c("all", "ttest", "kolmogorov", "kruskal", "ftest", "welch", "cochran", "wilcoxon", "copa", "s2n", "osum", "sam", "infcont", "yeohchi2", "pca", "affinity", "sliding","bimodal","brown","perm1"), nperms=100, ret.all = FALSE, debug=FALSE, except=c("affinity","bimodal","brown"), ...) {
	
	if (!is.factor(classes)) classes <- factor(classes)	
	tx <- t(x)
	nl <- nlevels(classes)
	l <- levels(classes)
	m <- matrix(Inf, nrow = nrow(x), ncol=nl)
	colnames(m) <- l	
	p <- list()
	if (("ttest" %in% test  ||  test == "all")  &&  ! "ttest" %in% except) {
		if (debug) cat("ttest...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) t.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["ttest"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("bimodal" %in% test  ||  test == "all")  &&  ! "bimodal" %in% except) {
		if (debug) cat("bimodal...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) p.bimodal(xwi[,k]))
		}
		m <- cbind(m, apply(t(x), 2, p.bimodal))
		p[["bimodal"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("infcont" %in% test  ||  test == "all") &&  ! "infcont" %in% except) {
		if (debug) cat("infcont...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- p.information.content(x, classes, l[i], large=FALSE)
		}
		p[["infcont"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("sliding" %in% test  ||  test == "all")  &&  ! "sliding" %in% except) {
		if (debug) cat("sliding...\n")
		uclasses <- unclass(classes)
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			.s <- sum(uclasses == i)
			.v <- numeric(nrow(x))
			for (j in 1:nrow(x)) {
				a <- order(x[j,])
				ac <- uclasses[a]
				.v[j] <- min(sapply(1:(ncol(x)-.s+1), function(k) {
					p.overlap(ncol(x), .s, .s, sum(ac[1:.s+k-1] == i))
				}))
			}
			m[,i] <- .v
		}
		p[["sliding"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("yeohchi2" %in% test  ||  test == "all")  &&  ! "yeohchi2" %in% except) {
		if (debug) cat("yeohchi2...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- p.yeohchi2(x, classes, l[i])
		}
		p[["yeohchi2"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("kolmogorov" %in% test  ||  test == "all")  &&  ! "kolmogorov" %in% except) {
		if (debug) cat("kolmogorov...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) ks.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["kolmogorov"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("kruskal" %in% test  ||  test == "all")  &&  ! "kruskal" %in% except) {
		if (debug) cat("kruskal...\n")
		p[["kruskal"]] <- apply(tx, 2, function(x) kruskal.test(x, classes)$p.value)
	}
	if ("ftest" %in% test  ||  test == "all") {
		if (debug) cat("ftest...\n")
		p[["ftest"]] <- apply(tx, 2, function(x) {
			sm <- summary(lm(x~classes))$fstatistic
			pf(sm[1], sm[2], sm[3], lower.tail=FALSE)
		} )
	#	if (debug) cat("ftest...\n")
	#	p[["ftest"]] <- apply(tx, 2, function(x) { tdf <- data.frame(x=x, g=classes); oneway.test(x ~ g, tdf, var.equal=TRUE)
	}
	if (("wilcoxon" %in% test  ||  test == "all")  &&  ! "wilcoxon" %in% except) {
		if (debug) cat("wilcoxon...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) wilcox.test(xwi[,k], xwj[,k])$p.value)
		}
		p[["wilcoxon"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("s2n" %in% test  ||  test == "all")  &&  ! "s2n" %in% except) {
		if (debug) cat("s2n...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			wi <- which(classes == l[i])
			wj <- which(classes != l[i])
			xwi <- t(x[,wi])
			xwj <- t(x[,wj])
			m[,i] <- sapply(1:nrow(x), function(k) (mean(xwi[,k],na.rm=TRUE)-mean(xwj[,k],na.rm=TRUE)) / (sd(xwi[,k],na.rm=TRUE)+sd(xwj[,k],na.rm=TRUE)))
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x,na.rm=TRUE))  # transform to a p-value like
		p[["s2n"]] <- if (ret.all) m else apply(m,1,min,na.rm=TRUE)
	}
	if (("welch" %in% test  ||  "oneway" %in% test  ||  test == "all")   &&  ! "welch" %in% except) {
		if (debug) cat("welch...\n")
		#xdf <- as.data.frame(tx)
		#xdf[,"testcluster"] <- classes
		tdf <- data.frame(x=tx[,1], g=classes)
		p[["welch"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=FALSE)$p.value })
	}
	if (("brown" %in% test  ||  "brown-forsythe" %in% test  ||  
		  "levene"  %in% test  || test == "all")   
		  &&  ! "brown" %in% except) {
		if (debug) cat("brown-forsythe/levene...\n")
		#xdf <- as.data.frame(tx)
		#xdf[,"testcluster"] <- classes
		library(lawstat)
		p[["brown"]] <- apply(tx, 2, function(x) { levene.test(x, g, location="median")$p.value })
	}
	if (("osum" %in% test  ||  test == "all")  &&  ! "osum" %in% except) {
		if (debug) cat("osum...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			m[,i] <- outlier.sum(x, group=classes, ref=l[i])
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x))  # transform to a p-value like
		p[["osum"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("copa" %in% test  ||  test == "all")  &&  ! "copa" %in% except) {
		if (debug) cat("copa...\n")
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			#m[,i] <- pmax(copa(x[,classes == l[i]], c(.75,.90,.95)), copa(x[,classes == l[i]], .90), copa(x[,classes == l[i]], .95))
			m[,i] <- apply(copa(x[,classes == l[i]], c(.75,.90,.95)), 2, max)
		}
		m <- glog(abs(m))
		m <- apply(m, 2, function(x) 1 - x / max(x))  # transform to a p-value like
		p[["copa"]] <- if (ret.all) m else apply(m,1,min)
	}
	if (("sam" %in% test  ||  test == "all")  &&  ! "sam" %in% except) {
		if (debug) cat("sam...\n")
		library(samr)
		data <- list(x=x, y=unclass(classes), geneid=as.character(1:nrow(x)),genenames=paste("g",as.character(1:nrow(x)),sep=""), logged2=TRUE)		
		samr.obj <- samr(data,  resp.type=ifelse(nl == 2, "Two class unpaired", "Multiclass"), nperms=nperms)
		p[["sam"]] <- if(ret.all) samr.obj else samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar)
	}
	if (("pca" %in% test  ||  test == "all")  &&  ! "pca" %in% except) {
		if (debug) cat("pca...\n")
		if (any(is.na(x))) {
			x <- na.impute(x)
		}
		for (i in 1:nl) {
			if (debug) cat(l[i],"\n")
			xsd <- apply(t(x[,classes == l[i]]),2,sd,na.rm=TRUE)
			wsd <- which(xsd == 0 | xsd < 0.0000000001)
			if (length(wsd) > 0) {
				# all values are the same, sd = 0, problems
				pc <- prcomp(t(x[-wsd,classes == l[i]]), scale=TRUE)
				m[-wsd,i] <- apply(abs(t(pc$rotation[,1:4])), 2, max)
				m[wsd,i] <- -Inf
			} else {
				pc <- prcomp(t(x[,classes == l[i]]), scale=TRUE)
				m[,i] <- apply(abs(t(pc$rotation[,1:4])), 2, max)
			}
		}
		pc <- apply(t(m), 2, max)
		p[["pca"]] <- if (ret.all) m else 1 - pc / max(pc)  # transform to a p-value like
	}
	if (("affinity" %in% test  ||  test == "all")  &&  ! "affinity" %in% except) {
		### this is not run in "all" mode because it takes long
		if (debug) cat("affinity propagation...\n")
		m[] <- 0
		npart <- trunc(nrow(x) / 4000 + 1)
		ngenes <- round(nrow(x) / npart,0)
		for (j in 0:(npart-1)) {
			g <- 1:ngenes + j * ngenes
			if (debug) cat("partition ",j,", genes ",min(g),":",max(g),"...\n")
			ptx <- tx[,min(ncol(tx),g)]
			for (i in 1:nl) {
				if (debug) cat(l[i],"\n")
				kpk <- affinity.propagation.clustering(data=tx[classes==l[i],g], debug=debug, max.iter=1000)
				m[g[kpk$exemplars],i] <- 1
			}
		}
		pc <- apply(t(m), 2, sum)
		pc[pc > 0] <- 1 / pc[pc > 0]
		p[["affinity"]] <- if (ret.all) m else 1 - pc / max(pc)  # transform to a p-value like
	}
	
	if (length(grep("perm1:",test)) > 0  &&  ! "perm1" %in% except) {
		xtest <- sub("perm1:","",test[grep("perm1:",test)])		
		thepis <- de.test(x, classes, test=xtest, debug=debug)
		minperm = 30
		if (nperms < minperm) nperms = minperm
		xstats <- list()
		xcounts <- c()
		#BREAKS ARE for p-values
		xbrks <- c(-1e300, 0, 10^-(100:20), 10^-(399:1/20), 1, +1e300)
		for (i in 1:nperms) {
			if (debug  ||  TRUE) {
				if ((i-1) %% 20 == 0) cat("\nPermutation ", i,"") else cat(".")
			}
			xp <- de.test(x, sample(classes), test=xtest, debug=debug)
			if (i <= minperm) {
				xstats[[i]] <- xp
			}
			if (i == minperm) {
				#xr <- range(unlist(lapply(xstats, range)))
				#xd <- abs(xr[2]-xr[1])				
				#xbrks <- c(-Inf, seq(xr[1]-xd,xr[2]+xd,len=100000), +Inf)
				xcounts <- numeric(length(xbrks))
				for (j in 1:minperm) {
					xh <- hist(xstats[[j]], breaks=xbrks, plot=FALSE)$counts
					xcounts <- xcounts + xh
				}
				xstats <- NULL
			} else {
				xh <- hist(xp, breaks=xbrks, plot=FALSE)$counts
				xcounts <- xcounts + xh
			}
		}
		if (debug) cat("\n")
		xsums <- cumsum(xcounts)
		thesum <- sum(xcounts)
		#wisna <- which(is.na(thepis))
		#if (length(wisna) > 0) thepis[wisna] <- 1
		xc <- cut(thepis, breaks=xbrks, labels=FALSE)
		newp <- (xsums[xc]+1)/thesum
		#if (length(wisna) > 0) {
		#	newp[wisna] <- NA
		#	thepis[wisna] <- NA
		#}
		p[["perm1"]] <- if (ret.all) data.frame(p.raw=thepis, perm.counts=xcounts, estimated.p=newp) else newp
	}

	#if ("ftest2" %in% test  ||  test == "all") {
	#	if (debug) cat("ftest2...\n")
	#	tdf <- data.frame(x=tx[,1], g=classes)
	#	p[["ftest2"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=TRUE)$p.value })
	#}
	#if ("cochran" %in% test  ||  test == "all") {
	#	### this seems to be the same than ftest
	#	if (debug) cat("cochran...\n")
	#	tdf <- data.frame(x=tx[,1], g=classes)
	#	p[["cochran"]] <- apply(tx, 2, function(x) { tdf$x <- x; oneway.test(x ~ g, tdf, var.equal=TRUE)$p.value })
	#}
	# Brown-Forsythe seems to be f-test with unequal variance so the same than welch or very similar
	#ebayes - limma
	#rank product
	#template matching
	#bga - between groups analysis
	#maxT
	#ROC
	#logistic regression
	if (length(p) == 1) p[[1]] else p
}

# based on Fray and ... Clustering by passing messages between data points Science 
# http://www.psi.toronto.edu/affinitypropagation/SimpleMATLABImplementation.txt
# http://www.psi.toronto.edu/affinitypropagation/software/apcluster.m
# MatLab Code:
#N=size(S,1); A=zeros(N,N); R=zeros(N,N); % Initialize messages
#S=S+(eps*S+realmin*100).*rand(N,N); % Remove degeneracies
#lam=0.5; % Set damping factor
#for iter=1:100
#	% Compute responsibilities
#	Rold=R;
#	AS=A+S; [Y,I]=max(AS,[],2);
#	for i=1:N AS(i,I(i))=-realmax; end;
#	[Y2,I2]=max(AS,[],2);
#	R=S-repmat(Y,[1,N]);
#	for i=1:N R(i,I(i))=S(i,I(i))-Y2(i); end;
#	R=(1-lam)*R+lam*Rold; % Dampen responsibilities
#
#	% Compute availabilities
#	Aold=A;
#	Rp=max(R,0); for k=1:N Rp(k,k)=R(k,k); end;
#	A=repmat(sum(Rp,1),[N,1])-Rp;
#	dA=diag(A); A=min(A,0); for k=1:N A(k,k)=dA(k); end;
#	A=(1-lam)*A+lam*Aold; % Dampen availabilities
#end;
#E=R+A; % Pseudomarginals
#I=find(diag(E)>0); K=length(I); % Indices of exemplars
#[tmp c]=max(S(:,I),[],2); c(I)=1:K; idx=I(c); % Assignments
affinity.propagation.clustering <- function(s = NULL, data = NULL, 
	similarityfun = function(x) cor(x, method="pearson"), 
	max.iter=100, lam = 0.5, 
	stable = 10, save.memory=nrow(s) > 5000, debug=FALSE) {
	setsdiag <- FALSE
	if (is.null(s)) {
		if (is.null(data)) {
			stop("similarities 's' nor 'data' was provided")
		}
		if (debug) catf("Estimating similarities...\n")
		#s = as.matrix(dist(data, method=method))
		#if (save.memory) for (k in 1:nrow(s))  s[,k] <- -s[,k]
		#else s <- -s
		s <- similarityfun(data)
		setsdiag <- TRUE
	}
	
	if (debug) catf("Initialization...\n")

	# initialize
	n <- nrow(s)
	a <- matrix(0, nrow=n, ncol=n) # affynities
	r <- matrix(0, nrow=n, ncol=n) # responsabilities
	as <- matrix(0, nrow=n, ncol=n) # sums of a and s
	
	#lam <- 0.5 # damping factor
	wm <- list()
	wm[[n]] <- 0 # just to create the list
	wmax <- matrix(0, nrow=n, ncol=2)
	wmax[,1] <- 1:n
	mx <- numeric(length=n)
	rold <- numeric(length=n)
	dia <- matrix(0, nrow=n, ncol=2)
	dia[,1] <- dia[,2] <- 1:n
	if (setsdiag) {
		if (debug) catf("Estimating shared similarity...\n")
		m <- s[2,1]
		for (k in 1:n) {
			m <- min(m, s[-k,k])
		}
		s[dia] <- m
	}
	# remove degenerancies
	if (debug) catf("Removing degenerancies...\n")
	for (i in 1:n) {
		if (i %% 100 == 0  && debug) catf(".")
		#s[,i] <- s[,i] + (.Machine$double.eps * s[,i] + 100 * .Machine$double.xmin) * runif(nrow(s))
	}
	
	pex <- 0	
	wmx <- numeric(n)
	yes <- 0
	we <- numeric(n)
	for (iter in 1:max.iter) {
		if (debug) catf("Iteration ", iter, ", clusters = ", length(we), ", stable iterations =", yes, "...\n")
		# compute responsabilities
		if (save.memory) for (k in 1:n)  as[,k] <- s[,k] + a[,k]
		else as <- s + a
		
		for (i in 1:n) {
			wm[[i]] <- which.max(as[i,])			# obtain which are the maximum in order to correct after r updates
			wmax[i,2] <- wm[[i]][1]
		}
		mx <- as[wmax]								# get and save maximums
		rold <- r[wmax]								# save maximums from r matrix
		for (k in 1:n) r[,k] <- (1-lam) * (s[,k] - mx) + lam * r[,k]			# update responsabilities assuming maximum in all 
		for (i in 1:n) {
			# correct only if there is only one maximum, otherwise removing that index would have the same maximum
			if (length(wm[[i]]) == 1) {
				r[i,wmax[i,2]] <- (1-lam) * (s[i,wmax[i,2]] - max(as[i,-wmax[i,2]])) + lam * rold[i]
			}
		}
		
		# compute availabilities
		sumr <- apply(r, 2, function(r) sum(pmax(0, r))) ## crude sum of columns but not considering
		rdia <- r[dia] 
		adiaold <- a[dia]
		maxrdia <- pmax(0, rdia)
		cte <- rdia + sumr - maxrdia
		for (i in 1:n) {
			a[i,] <- (1-lam) * pmin(0, cte - pmax(0, r[i,]))  + lam  * a[i,] ## use the crude sum to speed up but then correct
		}
		#for (i in 1:n)
		#	for (k in 1:n)
		#		a[i,k] <- min(0, rdia[k] + sumr[k] - maxrdia[k] - r[i,k])
		a[dia] <- (1-lam) * (sumr - maxrdia) + lam * adiaold

		#Obtain exemplars
		#Pseudomarginals
		for (i in 1:n) {
			wmx[i] <- which.max(a[i,] + r[i,])[1]
		}
		we <- which(wmx == 1:n)
		if (length(we) == length(pex)  &&  all(we == pex)) {
			yes <- yes + 1
			if (yes >= stable) {
				break
			}
		} else {
			yes <- 0
			pex <- we
		}
	}
	names(we) <- if (!is.null(data)) colnames(data)[we] else colnames(s)[we]
	k.clusters <- length(we)
	cluster <- numeric(n)
	if (!is.null(data)) names(cluster) <- colnames(data)
	if (debug) {
		#print(we)
		#print(wmx)
	}
	for (k in 1:n) {
		w <- which(we == wmx[k])
		if (length(w) > 0) cluster[k] <- w[1]
	}
	return (list(exemplars = we, k = k.clusters, cluster = cluster, iterations=iter, members=table(cluster)))
}


# adapted from invalid from gplots package
is.invalid <- function (x) 
{
	if (missing(x) || is.null(x) || length(x) == 0) 
		return(TRUE)
	if (is.list(x)) 
		return(all(sapply(x, is.invalid)))
	else if (is.vector(x)) 
		return(all(is.na(x)))
	else return(FALSE)
}


# adapted from heatmap.2 from gplots package
plot.heatmap <- 
function (x, 
	Rowv = TRUE, 
	Colv = if (symm) "Rowv" else TRUE, 
	distfun = dist, 
	hclustfun = function(x) { 
		if (any(is.na(x))) { 
			x[is.na(x)] <- max(x, na.rm=TRUE)
		}
		hclust(x, method=hclustmethod)
		}, 
	dendrogram = c("both", "row", "column", "none"), 
   	symm = FALSE, 
   	scale = c("none", "row", "column"), 
   	na.rm = TRUE, 
   	revC = identical(Colv, "Rowv"), 
   	add.expr, 
   	breaks, 
   	col = "redgreenblue", 
   	colsep, 
	rowsep, 
	sepcolor = "white", 
	sepwidth = c(0.05, 0.05), 
	cellnote, 
	notecex = 0.85, 
	notecol = "black", 
	na.color = par("bg"), 
	trace = c("none", "column", "row", "both"),
	tracecol = "gray", 
	hline = median(breaks), 
	vline = median(breaks), 
	linecol = tracecol, 
	margins = c(5, 5), 
	ColSideColors=NULL, 
	RowSideColors=NULL, 
	cexRow = 0.2 + 1/log10(nr), 
	cexCol = 0.2 + 1/log10(nc), 
	labRow = NULL, 
	labCol = NULL, 
	key = !uniform, 
	keysize = 1, 
	density.info = c("histogram", "density", "none"), 
	density.title = ifelse(density.info=="density", "Color Key & Density Plot", ifelse(density.info=="histogram","Color Key & Histogram","Color Key")),
	density.cex = 1,
	denscol = tracecol, 
	symkey = min(x <  0, na.rm = TRUE), 
	densadj = 0.25, 
	main = NULL, 
	xlab = NULL, 
	ylab = NULL, 
	lmat = NULL, 
	lhei = c(keysize*2, 4), 
	lwid = c(keysize, 4), 
	hclustmethod="ward.D", 
	ColSideLabels=NULL,
	RowSideLabels=NULL, 
	uniform=FALSE, 
	RowSideLabels.cols=min(5,length(unique(RowSideLabels))), 
	ColSideLabels.cols=min(10,length(unique(ColSideLabels))), 
	ColSideLabels.cex=3, 
	RowSideLabels.cex=2, 
	autozero=FALSE,
	...) 
{
	if (uniform) x <- apply(x,2,rank) / nrow(x)
	scale01 <- function(x, low = min(x), high = max(x)) {
		x <- (x - low)/(high - low)
		x
	}
	scale <- if (symm && missing(scale)) 
		"none"
	else match.arg(scale)
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (!missing(breaks) && (scale != "none")) 
		warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
			"specified can produce unpredictable results.", "Please consider using only one or the other.")
	if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1) 
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2) 
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote)) 
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	else if (is.logical(cellnote)  && cellnote) cellnote <- x
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
			c("both", "row"))) {
			if (is.logical(Colv) && (Colv)) 
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
				dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
			c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv)) 
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
				dendrogram, "'. Omitting column dendogram.")
		}
	}
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.integer(Rowv)) {
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd)) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		Rowv <- rowMeans(x, na.rm = na.rm)
		hcr <- hclustfun(distfun(x))
		ddr <- as.dendrogram(hcr)
		ddr <- reorder(ddr, Rowv)
		rowInd <- order.dendrogram(ddr)
		if (nr != length(rowInd)) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (is.numeric(Rowv) && length(Rowv) == nr) {
		rowInd <- Rowv
	}
	else {
		rowInd <- nr:1
	}

	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc) 
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.integer(Colv)) {
		hcc <- hclustfun(distfun(if (symm) 
			x
		else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd)) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		Colv <- colMeans(x, na.rm = na.rm)
		hcc <- hclustfun(distfun(if (symm) 
			x
		else t(x)))
		ddc <- as.dendrogram(hcc)
		ddc <- reorder(ddc, Colv)
		colInd <- order.dendrogram(ddc)
		if (nc != length(colInd)) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (is.numeric(Colv) && length(Colv) == nc) {
		colInd <- Colv
	}
	else {
		colInd <- 1:nc
	}
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow)) 
		labRow <- if (is.null(rownames(x))) 
			(1:nr)[rowInd]
		else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol)) 
		labCol <- if (is.null(colnames(x))) 
			(1:nc)[colInd]
		else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
		sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
		sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (autozero) {
		if (0 > min(x)*.9  && 0 < max(x)*.9) {
			xm <- min(abs(c(range(x, na.rm=TRUE))))
			x <- pmax(pmin(x,xm,na.rm=TRUE),-xm, na.rm=TRUE)
		}
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) 
		if (missing(col)) 
			breaks <- 16
		else breaks <- length(col) + 1
	if (length(breaks) == 1) {
		qs <- quantile(x, c(.025,.975), na.rm=TRUE)
		breaks <- seq(qs[1], qs[2], length = breaks)
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function") 
		col <- col(ncol)
	else if (is.character(col) && length(col) == 1) 
		col <- do.call(col, list(ncol))
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[] <- ifelse(x < min.breaks, min.breaks, x)
	x[] <- ifelse(x > max.breaks, max.breaks, x)
	if (missing(lhei) || is.null(lhei)) 
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid)) 
		lwid <- c(keysize, 4)
	ColSideLabels <- if (is.factor(ColSideColors)  &&  is.null(ColSideLabels)) levels(ColSideColors)[ColSideColors] else ColSideLabels
	RowSideLabels <- if (is.factor(RowSideColors)  &&  is.null(RowSideLabels)) levels(RowSideColors)[RowSideColors] else RowSideLabels
	if (missing(lmat) || is.null(lmat)) {
		lmat <- rbind(4:3, 2:1)
		if (!is.null(ColSideColors)) {
			if (length(ColSideColors) != nc) # !is.character(ColSideColors) || 
				stop("'ColSideColors' must be a character vector of length ncol(x)")
			if (is.factor(ColSideColors)) ColSideColors <- as.character(unclass(ColSideColors))
			else if (!is.character(ColSideColors)) ColSideColors <- as.character(ColSideColors)
			lmat <- rbind(lmat[1, ] + 1, c(NA, 1), lmat[2, ] + 
				1)
			lhei <- c(lhei[1], 0.2, lhei[2])
		}
		if (!is.null(RowSideColors)) {
			if (length(RowSideColors) != nr)  # !is.character(RowSideColors) ||  
				stop("LENGTH PROBLEMS:'RowSideColors' must be a character vector of length nrow(x)")
			if (is.factor(RowSideColors)) RowSideColors <- as.character(unclass(RowSideColors))
			else if (!is.character(RowSideColors)) RowSideColors <- as.character(RowSideColors)
			lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 
				1), 1), lmat[, 2] + 1)
			lwid <- c(lwid[1], 0.2, lwid[2])
		}
		lmat[is.na(lmat)] <- 0
	}
	if (length(lhei) != nrow(lmat)) 
		stop("lhei must have length = nrow(lmat) = ", nrow(lmat))
	if (length(lwid) != ncol(lmat)) 
		stop("lwid must have length = ncol(lmat) =", ncol(lmat))
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	if (! is.null(ColSideLabels)  ||  ! is.null(RowSideLabels)) {
		lmat <- rbind(lmat, max(lmat)+1)
		if (! is.null(ColSideLabels)  &&  ! is.null(RowSideLabels) ) lmat[nrow(lmat),ncol(lmat)] <- max(lmat)+1
		#lwid <- c(lwid, 1)
		lhei <- c(lhei, 1)
   	}
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	if (!is.null(RowSideColors)) {
		par(mar = c(margins[1], 0, 0, 0.5))
		image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
	}
	if (!is.null(ColSideColors)) {
		par(mar = c(0.5, 0, 0, margins[2]))
		image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
	}
	par(mar = c(margins[1], 0, 0, margins[2]))
	if (!symm || scale != "none") {
		x <- t(x)
		cellnote <- t(cellnote)
	}
	if (revC) {
		iy <- nr:1
		ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
		c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
		breaks = breaks, ...)
	if (!is.invalid(na.color) & any(is.na(x))) {
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
			col = na.color, add = TRUE)
	}
	axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
		cex.axis = cexCol)
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
		cex.axis = cexRow)
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	if (!missing(colsep)  && !is.null(colsep)) 
		for (csep in colsep) rect(xleft = csep + 0.5, ybottom = rep(0, 
			length(csep)), xright = csep + 0.5 + sepwidth[1], 
			ytop = rep(ncol(x) + 1, csep), lty = 1, lwd = 1, 
			col = sepcolor, border = sepcolor)
	if (!missing(rowsep) && !is.null(rowsep)) 
		for (rsep in rowsep) rect(xleft = 0, ybottom = (ncol(x) + 
			1 - rsep) - 0.5, xright = nrow(x) + 1, ytop = (ncol(x) + 
			1 - rsep) - 0.5 - sepwidth[2], lty = 1, lwd = 1, 
			col = sepcolor, border = sepcolor)
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		for (i in colInd) {
			if (!is.null(vline)) {
				vline.vals <- scale01(vline, min.scale, max.scale)
				abline(v = i - 0.5 + vline.vals, col = linecol, 
				  lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		for (i in rowInd) {
			if (!is.null(hline)) {
				hline.vals <- scale01(hline, min.scale, max.scale)
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote)) 
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
			col = notecol, cex = notecex)
	par(mar = c(margins[1], 0, 0, 0))
	if (dendrogram %in% c("both", "row")) {
		plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
	}
	else plot.new()
	par(mar = c(0, 0, if (!is.null(main)) 5 else 0, margins[2]))
	if (dendrogram %in% c("both", "column")) {
		plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
	}
	else plot.new()
	if (!is.null(main)) 
		title(main, cex.main = 1.5 * op[["cex.main"]])
	if (key) {
		par(mar = c(3, 3, 1, 1), cex = 0.75)
		#par(mar = c(5, 4, 2, 1), cex = 0.75)
		if (symkey) {
			max.raw <- max(abs(x), na.rm = TRUE)
			min.raw <- -max.raw
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
			xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		axis(1, at = xv, labels = lv)
		if (scale == "row") 
			mtext(side = 1, "Row Z-Score", line = 2, cex=density.cex)
		else if (scale == "column") 
			mtext(side = 1, "Column Z-Score", line = 2, cex=density.cex)
		else mtext(side = 1, "Value", line = 2, cex=density.cex)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
				lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y))
			title(density.title, cex.main=density.cex)
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2, cex=density.cex)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
				col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy))
			title(density.title, cex.main=density.cex)
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2, cex=density.cex)
		}
		else title(density.title, cex.main=density.cex)
	}
	else plot.new()
	if (! is.null(RowSideLabels)) {
		plot.new()
		par(mar=c(0,0,0,0))
		legend("bottom", legend=names(table(RowSideLabels)), fill=1:length(unique(RowSideLabels)), col=NULL, ncol=RowSideLabels.cols, box.lwd=0, cex=RowSideLabels.cex)
	}
	if (! is.null(ColSideLabels)) {
		plot.new()
		par(mar=c(0,0,0,0))
		legend("bottom", legend=names(table(ColSideLabels)), fill=1:length(unique(ColSideLabels)), col=NULL, ncol=ColSideLabels.cols, box.lwd=0, cex=ColSideLabels.cex)
	}
	invisible(list(rowInd = rowInd, colInd = colInd, rowDend = ddr, colDend = ddc))
}




# adapted from heatmap.2 from gplots package
plot.heatmap.info <- 
function (x, 
	Rowv = TRUE, # if function, this is ran to re-order
	Colv = if (symm) "Rowv" else TRUE, 
	RowvFactor = NULL, # character (colname) or col number of the RowSideFactor
	ColvFactor = NULL, # character (colname) or col number of the ColSideFactor
	distfun = dist, 
	hclustfun = function(x) { 
		if (any(is.na(x))) { 
			x[is.na(x)] <- max(x, na.rm=TRUE)
		}
		hclust(x, method=hclustmethod)
		}, 
	dendrogram = c("both", "row", "column", "none"), 
   	symm = FALSE, 
   	scale = c("none", "row", "column"), 
   	na.rm = TRUE,
   	revC = identical(Colv, "Rowv"), 
   	add.expr, 
   	breaks, 
   	col = "redgreenblue", 
   	colsep, 
	rowsep, 
	seplwd = 1,
	sepcolor = "white", 
	sepwidth = c(0.05, 0.05), 
	cellnote, 
	notecex = 0.85, 
	notecol = "black", 
	na.color = par("bg"), 
	trace = c("none", "column", "row", "both"),
	tracecol = "gray", 
	hline = median(breaks), 
	vline = median(breaks), 
	linecol = tracecol, 
	margins = c(5, 5), # bottom, right
	ColSideFactors=NULL, # data.frame, factors and numeric are "colored" by default or .col list, characters are assumed specific color
	ColSideFactors.counts=TRUE,
	ColSideFactors.type=c("i","h","l","p","o","t")[1], ### image, hist, line, point, overlap, any other pch char,"?" = default
	ColSideFactors.labels=c("left","right","none")[1], # if numeric, ==> "right" and the number of offset lines
	ColSideFactors.col=list(), # list, named color list
	ColSideFactors.height=rep(0.2,ncol(ColSideFactors)),
	ColSideFactors.correction=-0.5-.200*(if(!is.null(RowSideFactors)) ncol(RowSideFactors) else 0)*5*(keysize/sum(lwid)),
	ColSideFactors.sep=FALSE, # color or true
	ColSideFactors.block=FALSE, ## show blocks in each "cluster" (related to ColvFactor)
	ColSideFactors.frame=FALSE, ## show the labels of each factor within a frame, if # or string means color, if vector ==> per factor
	ColSideFactors.axis=1, ## cex for the axis labels
	ColSideFactors.order.function=NULL, ## if this is a function, it uses this function instead of the hiearchical clustering, it assumes that the returned values are relative order.
	ColSideFactors.ylim=NULL,
	ColSideFactors.labelling=FALSE,
	RowSideFactors=NULL, # data.frame, factors and numeric are "colored" by default or .col list, characters are assumed specific color
	RowSideFactors.counts=TRUE, 
	RowSideFactors.type=c("i","h","l","p", "o","t")[1], ### image, hist, line, point
	RowSideFactors.labels=c("top","bottom","none")[1], # if numeric ==> "bottom" and the number of offset lines
	RowSideFactors.col=list(), # list, named color list
	RowSideFactors.width=rep(0.2,ncol(RowSideFactors)),
	RowSideFactors.correction=0.7,
	RowSideFactors.sep=FALSE, # color or true and vector
	RowSideFactors.block=FALSE, ## show blocks in each "cluster" (related to RowvFactor)
	RowSideFactors.frame=FALSE, ## show the labels of each factor within a frame, if # or string means color, if vector ==> per factor
	RowSideFactors.axis=1, ## cex for the axis labels
	RowSideFactors.order.function=NULL, ## if this is a function, it uses this function instead of the hiearchical clustering, it assumes that the returned values are relative order.
	RowSideFactors.xlim=NULL,
	SideFactorsPriority=c("row","col")[1], # in the case both sidefactors are used, which should use the top-left space ?
	cexRow = NA, #0.2 + 1/log10(nr), 
	cexCol = NA, #0.2 + 1/log10(nc), 
	labRow = NULL, 
	labRow.col=TRUE,
	labCol = NULL, 
	labCol.col=TRUE,
	key = !uniform, 
	keysize = 1, 
	density.info = c("histogram", "density", "none"), 
	density.title = ifelse(density.info=="density", "Color Key & Density Plot", ifelse(density.info=="histogram","Color Key & Histogram","Color Key")),
	density.cex = 1,
	denscol = tracecol, 
	symkey = min(x <  0, na.rm = TRUE), 
	densadj = 0.25, 
	main = NULL, 
	xlab = NULL, 
	ylab = NULL, 
	lmat = NULL, 
	lhei = c(keysize*2, 4), 
	lwid = c(keysize, 4), 
	hclustmethod="ward.D", 
	uniform=FALSE, 
	factorize=FALSE,
	autozero=FALSE,
	overlap=NULL, # function (receives rowInd and colInd) to run just after image perhaps to overlap data in heatmap.
	...) 
{
	setdef <- function(default, frame, valores, as.list=FALSE) {
		if (is.null(frame)) return (NULL)
		v <- if(as.list) as.list(rep(default, ncol(frame))) else rep(default, ncol(frame))
		names(v) <- colnames(frame)
		if (!is.null(names(valores))) {
			v[names(valores)] <- if(as.list) valores else unlist(valores)
		} else if (length(valores) > 0) {
			v[1:length(valores)] <- unlist(valores)
			if (length(valores) < length(v)) {
				v[1:ncol(frame)] <- rep(valores,  ncol(frame))[1:ncol(frame)]
			}
		}
		v
	}
	if (!is.null(RowSideFactors)  && missing(RowSideFactors.type) && ncol(RowSideFactors) > 0) {
		for (i in 1:ncol(RowSideFactors)) 
			if ((is.numeric(RowSideFactors[,i])  ||  all(!is.na(as.numeric(as.character(RowSideFactors[,i]))))) &&
				length(unique(RowSideFactors[,i])) > length(palette())) {
				RowSideFactors[,i] <- round(as.numeric(as.character(RowSideFactors[,i])),5)
				RowSideFactors.type[i] = "h"
			}
	}
	RowSideFactors.type <- setdef("?", RowSideFactors, RowSideFactors.type)
	RowSideFactors.labels <- setdef("top", RowSideFactors, RowSideFactors.labels) ### top, bottom, none
	RowSideFactors.counts <- setdef(TRUE, RowSideFactors, RowSideFactors.counts)
	RowSideFactors.width <- setdef(0.2, RowSideFactors, RowSideFactors.width)
	RowSideFactors.col <- setdef(NA, RowSideFactors, RowSideFactors.col, as.list=TRUE)
	RowSideFactors.frame <- setdef(0, RowSideFactors, RowSideFactors.frame)
	RowSideFactors.sep <- setdef(0, RowSideFactors, ifelse(is.logical(RowSideFactors.sep),ifelse(RowSideFactors.sep,"gray","0"), RowSideFactors.sep))
	RowSideFactors.axis <- setdef(1, RowSideFactors, RowSideFactors.axis)


	ColSideFactors.type <- setdef("?", ColSideFactors, ColSideFactors.type)
	ColSideFactors.labels <- setdef("left", ColSideFactors, ColSideFactors.labels)  ## left, right, none
	ColSideFactors.counts <- setdef(TRUE, ColSideFactors, ColSideFactors.counts)
	ColSideFactors.height <- setdef(0.2, ColSideFactors, ColSideFactors.height)
	ColSideFactors.col <- setdef(NA, ColSideFactors, ColSideFactors.col, as.list=TRUE)
	ColSideFactors.frame <- setdef(0, ColSideFactors, ColSideFactors.frame)
	ColSideFactors.sep <- setdef(0, ColSideFactors, ifelse(is.logical(ColSideFactors.sep),ifelse(ColSideFactors.sep,"gray","0"), ColSideFactors.sep))
	ColSideFactors.axis <- setdef(0, ColSideFactors, ColSideFactors.axis)
	ColSideFactors.labelling <- setdef(0, ColSideFactors, ColSideFactors.labelling)
	if (!is.null(ColSideFactors)  && missing(ColSideFactors.type) && ncol(ColSideFactors) > 0) {
		for (i in 1:ncol(ColSideFactors)) 
			if ((is.numeric(ColSideFactors[,i])  ||  all(!is.na(as.numeric(as.character(ColSideFactors[,i]))))) &&
			length(unique(ColSideFactors[,i])) > length(palette())) {
				ColSideFactors[,i] <- round(as.numeric(as.character(ColSideFactors[,i])),5)	
				ColSideFactors.type[i] = "h"
			}
	}
		
		
	if (missing(lhei) || is.null(lhei)) 
		lhei <- c(keysize, 4)
	if (missing(lwid) || is.null(lwid)) 
		lwid <- c(keysize, 4)

	if (uniform) x <- apply(x,2,rank) / nrow(x)
	scale01 <- function(x, low = min(x, na.rm=TRUE), high = max(x, na.rm=TRUE)) {
		x <- (x - low)/(high - low)
		x
	}
	scale <- if (symm && missing(scale)) 
		"none"
	else match.arg(scale)
	if (factorize) {
		x <- matrix(as.numeric(factor(x)), nrow=nrow(x))
	}
	dendrogram <- match.arg(dendrogram)
	trace <- match.arg(trace)
	density.info <- match.arg(density.info)
	if (!missing(breaks) && (scale != "none")) 
		warning("Using scale=\"row\" or scale=\"column\" when breaks are", 
			"specified can produce unpredictable results.", "Please consider using only one or the other.")
	par.Colv <- Colv	
	if (is.function(Colv)) {
		Colv <- Colv(t(x))
	}
	if (is.function(Rowv)) {
		if (is.numeric(Colv)) {
			Rowv <- Rowv(x[,Colv])
		} else {
			Rowv <- Rowv(x)
		}
	}
	if ((Colv == "Rowv") && (!isTRUE(Rowv) || is.null(Rowv))) 
		Colv <- FALSE
	if (length(di <- dim(x)) != 2 || !is.numeric(x)) 
		stop("`x' must be a numeric matrix")
	nr <- di[1]
	nc <- di[2]
	if (nr <= 1 || nc <= 1) 
		stop("`x' must have at least 2 rows and 2 columns")
	if (!is.numeric(margins) || length(margins) != 2) 
		stop("`margins' must be a numeric vector of length 2")
	if (missing(cellnote)) 
		cellnote <- matrix("", ncol = ncol(x), nrow = nrow(x))
	else if (is.logical(cellnote)  && cellnote) { cellnote <- trim(format(x)); cellnote[!is.finite(x) | x == 0] <- ""; }
	if (!inherits(Rowv, "dendrogram")) {
		if (((!isTRUE(Rowv)) || (is.null(Rowv))) && (dendrogram %in% 
			c("both", "row"))) {
			if (is.logical(Colv) && (Colv)) 
				dendrogram <- "column"
			else dedrogram <- "none"
			warning("Discrepancy: Rowv is FALSE, while dendrogram is `", 
				dendrogram, "'. Omitting row dendogram.")
		}
	}
	if (!inherits(Colv, "dendrogram")) {
		if (((!isTRUE(Colv)) || (is.null(Colv))) && (dendrogram %in% 
			c("both", "column"))) {
			if (is.logical(Rowv) && (Rowv)) 
				dendrogram <- "row"
			else dendrogram <- "none"
			warning("Discrepancy: Colv is FALSE, while dendrogram is `", 
				dendrogram, "'. Omitting column dendogram.")
		}
	}
	ddr <- NULL
	ddc <- NULL
	if (inherits(Rowv, "dendrogram")) {
		ddr <- Rowv
		rowInd <- order.dendrogram(ddr)
	}
	else if (is.numeric(Rowv)) {
		#hcr <- hclustfun(distfun(x))
		#ddr <- as.dendrogram(hcr)
		#ddr <- reorder(ddr, Rowv)
		#rowInd <- order.dendrogram(ddr)
		rowInd = Rowv
		if (nr != length(rowInd)) 
			stop("row dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Rowv)) {
		if (is.null(RowvFactor)) {
			Rowv <- rowMeans(x, na.rm = na.rm)
			hcr <- hclustfun(distfun(x))
			ddr <- as.dendrogram(hcr)
			#ddr <- reorder(ddr, Rowv)
			rowInd <- order.dendrogram(ddr)
			if (nr != length(rowInd)) 
				stop("row dendrogram ordering gave index of wrong length")
		} else {
			f <- RowSideFactors[,RowvFactor]
			if (is.factor(f)) f <- as.factor(as.character(f)) # to remove labels not used
			f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
			u <- levels(as.factor(f))
			nl <- length(u) # works even for numeric
			ddr <- list()
			rowInd <- NULL
			.l <- NULL
			for (i in 1:nl) {
				w <- which(f == u[i])
				if (length(w) > 1) {
					Rowv <- rowMeans(x[w,,drop=FALSE], na.rm = na.rm)
					if (!is.function(RowSideFactors.order.function)) {
						hcc <- hclustfun(distfun(x[w,,drop=FALSE]))
						ddr[[i]] <- as.dendrogram(hcc)
						#ddr[[i]] <- reorder(ddr[[i]], Colv)
						rowInd <- c(rowInd, w[order.dendrogram(ddr[[i]])])
					} else {
						rof <- RowSideFactors.order.function(x[w,], cluster=i)
						rowInd <- c(rowInd, w[rof])
						ddr[[i]] <- FALSE
					}
				} else {
					ddr[[i]] <- FALSE
					rowInd <- c(rowInd, w)
				}
				.l <- c(.l, length(w))
			}
			if (length(lhei) != nl+1) {
				par(mar=c(margins[1],0,0,margins[2]))
				avnm <- (par("cin")[2] * (margins[1]-RowSideFactors.correction*margins[1]/5)) / (par("fin")[2])
				.l2 <- lhei[2]*.l*(1-avnm)/nrow(RowSideFactors)
				.l2[1] <- lhei[2]-sum(.l2[-1])
				lhei <- c(lhei[1], rev(.l2))
			}
			attr(ddr, "multiplier") <- .l/nrow(RowSideFactors)
		}
	}
	else {
		rowInd <- nr:1
	}
	if (inherits(Colv, "dendrogram")) {
		ddc <- Colv
		colInd <- order.dendrogram(ddc)
	}
	else if (identical(Colv, "Rowv")) {
		if (nr != nc) 
			stop("Colv = \"Rowv\" but nrow(x) != ncol(x)")
		if (exists("ddr")) {
			ddc <- ddr
			colInd <- order.dendrogram(ddc)
		}
		else colInd <- rowInd
	}
	else if (is.numeric(Colv)) {
		#hcc <- hclustfun(distfun(if (symm) 
		#	x
		#else t(x)))
		#ddc <- as.dendrogram(hcc)
		#ddc <- reorder(ddc, Colv)
		#colInd <- order.dendrogram(ddc)
		colInd <- Colv
		if (nc != length(colInd)) 
			stop("column dendrogram ordering gave index of wrong length")
	}
	else if (isTRUE(Colv)) {
		if (is.null(ColvFactor)) {
			Colv <- colMeans(x, na.rm = na.rm)
			hcc <- hclustfun(distfun(if (symm) 
				x
			else t(x)))
			ddc <- as.dendrogram(hcc)
			#ddc <- reorder(ddc, Colv)
			colInd <- order.dendrogram(ddc)
			if (nc != length(colInd)) 
				stop("column dendrogram ordering gave index of wrong length")
		} else {
			f <- ColSideFactors[,ColvFactor]
			if (is.factor(f)) f <- as.factor(as.character(f)) # to remove labels not used
			f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
			u <- levels(as.factor(f))
			nl <- length(u) # works even for numeric
			ddc <- list()
			colInd <- NULL
			.l <- NULL
			for (i in 1:nl) {
				w <- which(f == u[i])
				if (length(w) > 1) {
					Colv <- colMeans(x[,w], na.rm = na.rm)
					if (!is.function(ColSideFactors.order.function)) {
						hcc <- hclustfun(distfun(t(x[,w])))
						ddc[[i]] <- as.dendrogram(hcc)
						#ddc[[i]] <- reorder(ddc[[i]], Colv)
						colInd <- c(colInd, w[order.dendrogram(ddc[[i]])])
					} else {
						cof <- ColSideFactors.order.function(t(x[,w]), cluster=i)
						colInd <- c(colInd, w[cof])
						ddc[[i]] <- FALSE
					}
				} else {
					ddc[[i]] <- FALSE
					colInd <- c(colInd, w)
				}
				.l <- c(.l, length(w))
			}
			
			if (length(lwid) != nl+1) {
				#avnm <- (par("cin")[2] * (margins[2]-(11.8-3.5*par("pin")[1]))) / (par("pin")[1])
				par(mar=c(margins[1],0,0,margins[2]))
				avnm <- (par("cin")[1] * (margins[2]-ColSideFactors.correction*margins[2]/5)) / (par("fin")[1])
				#avnm <- par("pin")[1]-margins[2]*strheight("JjGgPp",cex=cexRow, units="inches")
				#avnm <- (strheight("JjGgPp",par("cex.axis"), units="inches") * (margins[2])) / (par("pin")[1])
				.l2 <- lwid[2]*.l*(1-avnm)/nrow(ColSideFactors)
				#.l[length(.l)] <- .l[length(.l)]+margins[2]*(1+1.5*margins[2]/10)*sum(lwid)*par("pin")[1]/par("fin")[1]
				#.l2 <- lwid[2]*.l/nrow(ColSideFactors)
				.l2[length(.l2)] <- lwid[2]-sum(.l2[-length(.l2)])
				lwid <- c(lwid[1], .l2)
				#lwid <- c(lwid[1], .l*lwid[2]/sum(.l))
			}
			attr(ddc, "multiplier") <- .l/nrow(ColSideFactors)
		}
	}
	else {
		colInd <- 1:nc
	}
	x <- x[rowInd, colInd]
	x.unscaled <- x
	cellnote <- cellnote[rowInd, colInd]
	if (is.null(labRow)) 
		labRow <- if (is.null(rownames(x))) 
			(1:nr)[rowInd]
		else rownames(x)
	else labRow <- labRow[rowInd]
	if (is.null(labCol)) 
		labCol <- if (is.null(colnames(x))) 
			(1:nc)[colInd]
		else colnames(x)
	else labCol <- labCol[colInd]
	if (scale == "row") {
		x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
		sx <- apply(x, 1, sd, na.rm = na.rm)
		x <- sweep(x, 1, sx, "/")
	}
	else if (scale == "column") {
		x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
		sx <- apply(x, 2, sd, na.rm = na.rm)
		x <- sweep(x, 2, sx, "/")
	}
	if (autozero) {
		if (0 > min(x, na.rm=TRUE)*.9  && 0 < max(x, na.rm=TRUE)*.9) {
			xm <- min(abs(c(range(x, na.rm=TRUE))))
			x <- pmax(pmin(x,xm,na.rm=TRUE),-xm, na.rm=TRUE)
		}
	}
	if (missing(breaks) || is.null(breaks) || length(breaks) < 1) 
		if (missing(col)  || is.null(col)) 
			breaks <- 16
		else breaks <- length(col) + 1
	if (length(breaks) == 1) {
		qs <- quantile(x, c(.025,.975,.5), na.rm=TRUE)
		if (qs[1] == qs[2] || qs[1] == qs[3] || qs[2] == qs[3]) qs <- quantile(x, c(0,1), na.rm=TRUE)
		breaks <- seq(qs[1], qs[2], length = breaks)
	}
	nbr <- length(breaks)
	ncol <- length(breaks) - 1
	if (class(col) == "function") 
		col <- col(ncol)
	else if (is.character(col) && length(col) == 1) 
		col <- do.call(col, list(ncol))
	min.breaks <- min(breaks)
	max.breaks <- max(breaks)
	x[] <- ifelse(x < min.breaks, min.breaks, x)
	x[] <- ifelse(x > max.breaks, max.breaks, x)
   	.c <- ifelse(class(ddc)=="list",length(ddc),1)
   	.r <- ifelse(class(ddr)=="list",length(ddr),1)
	if (missing(lmat) || is.null(lmat)) {
		#lmat <- rbind(4:3, 2:1)
		lmat <- matrix(1, ncol=.c+1, nrow=.r+1)
		lmat[1,1] <- .c+.r+2
		lmat[1,2:(.c+1)] <- 1:.c+1
		lmat[(.r+1):2,1] <- 1:.r+.c+1
		
		if (!is.null(ColSideFactors)  && !is.null(ColSideFactors)) {
			# colsidefactors is a matrix/data.matrix of factors
			if (nrow(ColSideFactors) != ncol(x)) 
				stop("ColSideFactors must have rows = ncol(x).", " Rows=",nrow(ColSideFactors), ", ncol(x)=", ncol(x))
			n <- ncol(ColSideFactors)
			lmat <- lmat + 2*n
			.new <- t(matrix(1:(n*2), nrow=2))
			for (i in 1:n) {
				lmat <- rbind(lmat[1:i,], 0, lmat[-(1:i),])
				lmat[i+1,1] <- .new[i,1]
				lmat[i+1,-1] <- .new[i,-1]
			}
			lhei <- c(lhei[1], ColSideFactors.height, lhei[-1]*(sum(lhei[-1])-ColSideFactors.height)/sum(lhei[-1]))
		}
		if (!is.null(RowSideFactors)) {
			# rowsidefactors is a matrix/data.matrix of factors
			if (nrow(RowSideFactors) != nrow(x)) 
				stop("RolSideFactors must have rows = nrow(x).", " Rows=",nrow(RowSideFactors), ", nrow(x)=", nrow(x))

			n <- ncol(RowSideFactors)
			lmat <- lmat + 2*n
			.new <- matrix(1:(n*2), nrow=2)
			for (i in 1:n) lmat <- cbind(lmat[,1], lmat)
			for (i in 1:n) {
				.rows <- if (SideFactorsPriority == "row") {
						 nrow(lmat)-if (class(ddr)=="list") length(ddr) else 1
						} else 1+if(!is.null(ColSideFactors)) ncol(ColSideFactors) else 0
				lmat[1:.rows, i+1] <- .new[1,i]
				lmat[-(1:.rows),i+1] <- .new[-1,i]
			}

			lwid <- c(lwid[1], RowSideFactors.width, lwid[-1]*(sum(lwid[-1])-RowSideFactors.width)/sum(lwid[-1]))
		}
		lmat[is.na(lmat)] <- 0
	}
	if (length(lhei) != nrow(lmat)) {
		if (length(lhei) < nrow(lmat)) stop("lhei must have length = nrow(lmat) = ", nrow(lmat), " instead of ",length(lhei), ":=",paste(lhei,collapse=","))
		lhei <- lhei[1:nrow(lmat)]
	}
	if (length(lwid) != ncol(lmat)) {
		if (length(lwid) < ncol(lmat)) stop("lwid must have length = ncol(lmat) =", ncol(lmat), " instead of ",length(lwid), ":=",paste(lwid,collapse=","))
		lwid <- lwid[1:ncol(lmat)]
	}
	op <- par(no.readonly = TRUE)
	on.exit(par(op))
	layout(lmat, widths = lwid, heights = lhei, respect = FALSE)
	getquasifactor <- function(f) {
			an <- as.numeric(as.character(f))
			if (any(is.na(an))) {
				if (is.factor(f)  ||  is.logical(f) || is.character(f)) 
					f <- as.factor(as.character(f)) # to remove labels not used
			} else {
				f <- an
			}
			f
	}
	if (!is.null(RowSideFactors)) {
		for (i in 1:ncol(RowSideFactors)) {
			f <- suppressWarnings(getquasifactor(RowSideFactors[,i]))
			f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
			f.n.col <- f.n
			u <- levels(if (is.numeric(f)) factor(f) else as.factor(as.character(f))) #levels(as.factor(as.character(f)))
			nl <- length(u) # works even for numeric
			kol <- if (nl > length(palette())) rainbow(nl*1.1)[1:nl] else 1:nl
			if (all(!is.na(RowSideFactors.col[[colnames(RowSideFactors)[i]]]))) {
				kol <- RowSideFactors.col[[colnames(RowSideFactors)[i]]]
				if(length(kol) == length(f)) {
					f.n.col <- 1:length(kol)
				} else {
					#if (length(kol) < nl) 
					## quantile like
					nok <- round(cumsum(rep(nl/length(kol), length(kol))),0)
					kol <- rep(kol, times=c(nok[1],diff(nok)))					
				}
			}
			tb <- table(f)

			ofl <- 3
			rsfl <- RowSideFactors.labels[[i]]
			if (is.numeric(rsfl)) {
				ofl <- rsfl
				rsfl <- "bottom"
			}
			if (is.character(f) || rsfl == FALSE) rsfl <- "none"
			par(mar = c(0, 0, 0, 0))
			if (RowSideFactors.type[i] == "?" || !(substr(RowSideFactors.type[i],1,1) %in% c("i","h","p","l","o","t","b"))) {
				if ((is.numeric(RowSideFactors[,i])  ||  all(!is.na(as.numeric(as.character(RowSideFactors[,i]))))) &&
					length(unique(RowSideFactors[,i])) > 8) {
					RowSideFactors[,i] <- round(as.numeric(as.character(RowSideFactors[,i])),5)
					RowSideFactors.type[i] = "h"
				} else {
					##RowSideFactors.type[i] = ifelse(nrow(RowSideFactors) < 30, "t", "i")
					RowSideFactors.type[i] = ifelse(length(unique(RowSideFactors[,i])) > length(palette()), "t", "i")
				}			
			}
			if (length(grep("\\+$",as.character(f))) > 0 && length(grep("\\-$",as.character(f))) == 0) {
				### Posible censoring
				xf <- trim(as.character(f))
				yf <- as.numeric(gsub("\\+","",xf))
				if (sum(is.finite(yf)) > 4) {
					### yes, it is censored
					f <- yf
					f[!is.finite(f)] <- NA
					kol <- c(1,7)
					f.n <- ifelse(gsub("[^+]","",xf)=="+",2,1)
					f.n.col <- f.n
					RowSideFactors.type[i] <- "b"
					RowSideFactors.counts[i] <- 1
					xind <- ifelse(!is.finite(f),"NA",ifelse(f.n[is.finite(f)]==2,"Censored","Event"))
					tb <- table(xind)[c(2,1,3)[1:length(unique(xind))]]
				}
			}
			if (RowSideFactors.type[i] != "t" && rsfl != "none") {
				plot(0,0,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",frame.plot=FALSE,type="n")
				box(col=RowSideFactors.frame[i])
			} else {
				plot.new()
			}
			if (rsfl == "top") {
				if ((!is.null(RowvFactor)  &&  
					(RowvFactor == i   ||  RowvFactor == colnames(RowSideFactors)[i])) ||
					(is.null(RowvFactor) && i == 1)) {
					if (!is.null(ddr)) attr(ddr, "kolor") <- kol
					if (is.logical(labRow.col) && labRow.col) {
						labRow.col <- kol[f.n.col]
					}
				}
				## plot.string adaptations
				if (RowSideFactors.counts[i]) u <- names(tb)
				if (any(nchar(u)==0)) u[nchar(u)==0] <- "NA"
				many <- length(u) > 20
				val <- if(many) seq(1,length(u),len=10) else 1:length(u)
					#	paste(u,if (many) cumsum(tb) else tb,sep=if (many) "->" else ":") 
				s <- if (RowSideFactors.counts[i]) 
						paste(u,if (many) ""		 else round(tb,4),sep=if (many) ""   else ":") 
					 else u
				if (RowSideFactors.type[i] != "t" || RowSideFactors.counts[i]) 
					plot.string(0,0,1,1,s[val],col=kol[val],wrap=TRUE,srt=90, cex.max=3.5)
				## 
#				
#				a <- 0
#				s <- if (RowSideFactors.counts[i]) paste(names(tb),tb,sep=":") else u
#				vin <- par("fin")[2]
#				sw <- strwidth(paste(s,collapse="."),units="inches")/vin
#				xadj <- 0.5
#				.x <- 0.5
#				t.cex = 1
#				print.levels <- 1:nl
#				if (sw > 1) {
#					if (sw > 2) {
#						t.cex = .66
#						sw <- strwidth(paste(s,collapse="."), cex=t.cex, units="inches") / vin
#						xnl <- trunc(nl / 2)
#						while (sw > 3) {
#							print.levels <- seq(1,nl,length.out=xnl)
#							sw <- strwidth(paste(s[print.levels],collapse="."), cex=t.cex, units="inches") / vin
#							xnl <- xnl - 1
#						}
#					}
#					xadj <- 1
#					.x <- 0
#				}
#				for (j in print.levels) {
#					uj <-  ifelse(j==1 && nchar(u[j]) == 0, 1, u[j])
#					ujx <- ifelse(j==1 && nchar(u[j]) == 0, "NA", u[j]) 
#					s <- if (RowSideFactors.counts[i]) paste(ujx,tb[uj],sep=":") else ujx
#					stw <- strwidth(s, cex=t.cex, units="inches")/vin
#					if (a+stw > 1) {
#						a <- 0
#						.x <- .x + strheight(s, cex=t.cex, units="inches")*1.15/par("fin")[1]
#					}
#					text(.x,a,labels=s,col=kol[j], cex=t.cex, srt=90, adj=c(0,xadj))
#					a <- a + stw + strwidth(".", cex=t.cex, units="inches")*2/vin
#				}
			}
			#text(0, 0, labels=colnames(RowSideFactors)[i], adj=c(0, 1), srt=90)
			par(mar = c(margins[1], 0, 0, 0.25))
			if (RowSideFactors.type[i] == "i") {
				image(rbind(1:nr), col = if(is.character(f) && all(f[rowInd] %in% colours())) f[rowInd] else kol[f.n.col][rowInd], axes = FALSE)		  
			}
			if (RowSideFactors.type[i] != "i") {				
				x.ty <- substr(RowSideFactors.type[i],1,1)
				x.type <- x.ty
				x.pch <- 20 #ifelse(x.ty %in% c("p","h","l","o"),if (x.ty=="o") 20 else ".",x.ty)
				x.cex <- if(nchar(RowSideFactors.type[i]) < 2) 1 else as.numeric(substr(RowSideFactors.type[i],2,20))
				x.lwd <- if(nchar(RowSideFactors.type[i]) < 2) 1 else as.numeric(substr(RowSideFactors.type[i],2,20))
				x.ty <- ifelse(x.ty %in% c("p","h","l","o","b","t"),x.ty,"p")
				x.lim <- NULL
				if (length(RowSideFactors.xlim) && any(names(RowSideFactors.xlim)==colnames(RowSideFactors)[i])) {
					x.lim <- RowSideFactors.xlim[[colnames(RowSideFactors)[i]]]
				}
				plot((if (x.type == "t" || any(is.na(as.numeric(f)))) 1:length(f) else as.numeric(f))[rowInd],
					 1:nr, type="n",axes=FALSE,xlab="",ylab="",main="",
					 frame.plot=RowSideFactors.frame[i], xlim=x.lim)
				if (!class(f) %in% c("character")) {
					if (is.null(x.lim)) try(par(usr=c(max(finites(as.numeric(f)),na.rm=TRUE),min(finites(as.numeric(f)),na.rm=TRUE)-ifelse(x.ty %in% c("p","h","t"),1,0),0.5,0.5+nr)))
					else try(par(usr=c(par("usr")[1:2],0.5,nr+0.5)))
				}
				mid <- mean(par("usr")[1:2])
				if (x.type == "h") {
					p1 = min(finites(as.numeric(f)),na.rm=TRUE)-1
					fri = as.numeric(f[rowInd])
					kfri = kol[f.n.col][rowInd]
					lr = 1:nr
					#for (lr in 1:nr)
					#   lines(c(p1,fri[lr]), c(lr,lr), type="l", col=kfri[lr], pch=x.pch, lwd=x.lwd, cex=x.cex)
					segments(p1,lr,fri,lr,col=kfri[lr], lwd=x.lwd)
				} else if (x.type == "t") {
					box(col=RowSideFactors.frame[i])
					par(usr=c(0,1,par("usr")[3:4]))
					x.cex <- 1.5
					#xs <- which.max(strwidth(trim(as.character(f))))
					#cat(as.character(f),"\n")
					#cat("Max=",xs,"[",as.character(f[xs]),"]\n")
					#while (	(strheight("[](){}GgJj",units="user",cex=x.cex) > 0.95 || 
					#		strwidth(as.character(f[xs]),units="user",cex=x.cex) > 0.95) && x.cex > 0.3) # strwidth(as.character(f[xs]),units="user",cex=x.cex) > 1*RowSideFactors.width[i]/0.2
					#	x.cex = x.cex - 0.05
					pin <- par("pin")
					x.cex <- adjustCEXinches(as.character(f), width=pin[1], height=pin[2])
					if (any(!is.na(f.n)) && any(f.n[!is.na(f.n)] < 0)) {
						f.n[!is.na(f.n) & f.n < 0] = 1
					}
					#text(0.5, 1:nr, labels=f[rowInd], col=kol[rowInd], pch=x.pch, lwd=x.lwd, cex=x.cex)
					text(0.5, 1:nr, labels=f[rowInd], col=kol[f.n.col][rowInd], pch=x.pch, lwd=x.lwd, cex=x.cex)
					#for (lr in 1:nr)
					#	text(0.5, lr, labels=f[rowInd][lr], type=x.ty, col=kol[f.n][rowInd][lr], pch=x.pch, lwd=x.lwd, cex=x.cex)
				} else if (x.ty == "b") {
					grueso <- as.numeric(substr(RowSideFactors.type[i],2,100))
					if (!is.finite(grueso)) grueso <- 0.9
					symbols(f[rowInd]/2,1:nr+0.5*0,rectangles=data.matrix(data.frame(w=abs(f[rowInd]),h=abs(grueso))),bg=kol[f.n.col][rowInd],fg=kol[f.n.col][rowInd], add=TRUE, inches=FALSE)
				} else {
					points(f[rowInd],1:nr, type=x.ty, col=kol[f.n.col][rowInd], pch=x.pch, lwd=x.lwd, cex=x.cex)
				}
			}
			if (RowSideFactors.sep[i] != 0) {
				abline(h=seq(par("usr")[3],par("usr")[4],length.out=nr+1), col=RowSideFactors.sep[i], lwd=0.5)
			}
			axis(1, at=mean(par("usr")[1:2]), labels=colnames(RowSideFactors)[i], las=3, cex.axis=RowSideFactors.axis[i])
			if (rsfl =="bottom") {
				os <- ofl
				for (j in 1:nl) {
					uj <-  ifelse(j==1 && nchar(u[j]) == 0, 1, u[j])
					ujx <- ifelse(j==1 && nchar(u[j]) == 0, "NA", u[j]) 
					s <- trim(if (RowSideFactors.counts[i]) paste(ujx,tb[uj],sep=":") else ujx)
					axis(at=mean(par("usr")[1:2]), labels=s, side=1, line=os, col.axis=kol[j], las=3, tick=FALSE)
					os <- os + ifelse(nchar(s) < 2, .75, ifelse(nchar(s) > 2, nchar(s) / 3 + 0.75, 1))
				}
			}
		}
	}
	if (!is.null(ColSideFactors)) {
		for (i in 1:ncol(ColSideFactors)) {
			f <- suppressWarnings(getquasifactor(ColSideFactors[,i]))
			f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
			f.n.col <- f.n
			u <- levels(if (is.numeric(f)) factor(f) else as.factor(as.character(f))) #levels(as.factor(as.character(f)))
			nl <- length(u) # works even for numeric
			kol <- if (nl > length(palette())) rainbow(nl*1.1)[1:nl] else 1:nl
			if (all(!is.na(ColSideFactors.col[[colnames(ColSideFactors)[i]]]))) {
				kol <- ColSideFactors.col[[colnames(ColSideFactors)[i]]]
				if(length(kol) == length(f)) {
					f.n.col <- 1:length(kol)
				} else  {
					#if (length(kol) < nl)
					## quantile like
					nok <- round(cumsum(rep(nl/length(kol), length(kol))),0)
					kol <- rep(kol, times=c(nok[1],diff(nok)))
				} 
			}
			tb <- table(f)

			ofl <- 3
			csfl <- ColSideFactors.labels[[i]]
			if (is.numeric(csfl)) {
				ofl <- csfl
				csfl <- "right"
			}
			if (is.character(f) || csfl == FALSE) csfl <- "none"
			par(mar = c(0, 0, 0, 0))
			if (ColSideFactors.type[i] == "?" || !(substr(ColSideFactors.type[i],1,1) %in% c("i","h","p","l","o","t","b"))) {
				if ((is.numeric(ColSideFactors[,i])  ||  all(!is.na(as.numeric(as.character(ColSideFactors[,i]))))) &&
				length(unique(ColSideFactors[,i])) > 8) {
					ColSideFactors[,i] <- round(as.numeric(as.character(ColSideFactors[,i])),5)	
					ColSideFactors.type[i] = "h"
				} else {
					ColSideFactors.type[i] = ifelse(length(unique(ColSideFactors[,i])) > length(palette()), "t", "i")
				}
			}
			if (length(grep("\\+$",as.character(f))) > 0 && length(grep("\\-$",as.character(f))) == 0) {
				### Posible censoring
				xf <- trim(as.character(f))
				yf <- as.numeric(gsub("\\+","",xf))
				if (sum(is.finite(yf)) > 4) {
					### yes, it is censored
					f <- yf
					f[!is.finite(f)] <- NA
					kol <- c(1,7)
					f.n <- ifelse(gsub("[^+]","",xf)=="+",2,1)
					f.n.col <- f.n
					ColSideFactors.type[i] <- "b"
					ColSideFactors.counts[i] <- 1
					xind <- ifelse(!is.finite(f),"NA",ifelse(f.n[is.finite(f)]==2,"Censored","Event"))
					tb <- table(xind)[c(2,1,3)[1:length(unique(xind))]]
				}
			}
			#cat("Col=",i,colnames(ColSideFactors)[i],ColSideFactors.type[i],"\n")
			if (ColSideFactors.type[i] != "t" && csfl != "none") {
				plot(0,0,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",frame.plot=FALSE,type="n")
				box(col=ColSideFactors.frame[i])
			} else {
				plot.new()
			}
			if (csfl == "left") {
				if (!is.null(ColvFactor)  &&  
					(ColvFactor == i   ||  ColvFactor == colnames(ColSideFactors)[i]) || 
					(is.null(ColvFactor) && i == 1)) {
					if (!is.null(ddc)) attr(ddc, "kolor") <- kol
					if (is.logical(labCol.col) && labCol.col) {
						labCol.col <- kol[f.n.col]
					}
				}
				## plot.string adaptations
				if (ColSideFactors.counts[i]) u <- names(tb)
				if (any(nchar(u)==0)) u[nchar(u)==0] <- "NA"
				many <- length(u) > 20
				val <- if(many) seq(1,length(u),len=10) else 1:length(u)
					#	paste(u,if (many) cumsum(tb) else tb,sep=if (many) "->" else ":") 
				s <- if (ColSideFactors.counts[i]) 
						paste(u,if (many) ""		 else round(tb,4),sep=if (many) ""   else ":") 
					 else u
				if (ColSideFactors.type[i] != "t" || ColSideFactors.counts[i])
					plot.string(1,1,0,0,s[val],col=kol[val],wrap=TRUE,srt=0, cex.max=3.5)
				## 
#				a <- 1
#				s <- if (ColSideFactors.counts[i]) paste(names(tb),tb,sep=":") else u
#				t.cex = 2
#				sw <- strwidth(paste(s,collapse="."), cex=t.cex)
#				yadj <- 0.5
#				y <- 0.5
#				print.levels <- 1:nl
#				if (sw > 1.03) {
#					if (sw > 2) {
#						t.cex = .66
#						sw <- strwidth(paste(s,collapse="."), cex=t.cex)
#						xnl <- trunc(nl / 2)
#						while (sw > 3) {
#							print.levels <- seq(1,nl,length.out=xnl)
#							sw <- strwidth(paste(s[print.levels],collapse="."), cex=t.cex)
#							xnl <- xnl - 1
#						}
#					}
#					yadj <- 1
#					y <- 1
#				}
#				for (j in print.levels) {
#					uj <-  ifelse(j==1 && nchar(u[j]) == 0, 1, u[j])
#					ujx <- ifelse(j==1 && nchar(u[j]) == 0, "NA", u[j]) 
#					s <- if (ColSideFactors.counts[i]) paste(ujx,tb[uj],sep=":") else as.character(uj)
#					stw <- strwidth(s, cex=t.cex)
#					if (stw-0.03 > a) {
#						a <- 1
#						y <- y - strheight(s, cex=t.cex)
#					}
#					text(a, y, labels=s, adj=c(1, yadj), col=kol[j], cex=t.cex)
#					a <- a - stw - strwidth(".", cex=t.cex)*0.75
#				}
			}
			#text(1, 0.5, labels=colnames(ColSideFactors)[i], adj=c(1, 0.5))
			par(mar = c(0.25, 0, 0, margins[2]))
   			if (ColSideFactors.type[i] == "i") {
				image(cbind(1:nc), col = if(is.character(f) && all(f[colInd] %in% colours())) f[colInd] else kol[f.n.col][colInd], axes = FALSE)
				if (ColSideFactors.labelling[i] > 0) {
					# Problem: it is not garantee that cols are ordered
				}
   			}
   			if (ColSideFactors.type[i] != "i") {
				x.ty <- substr(ColSideFactors.type[i],1,1)
				x.pch <- 20 #ifelse(x.ty %in% c("p","h","l","o"),if (x.ty=="o") 20 else ".",x.ty)
				x.ty <- ifelse(x.ty %in% c("p","h","l","o","b","t"),x.ty,"p")
				x.cex <- if(nchar(ColSideFactors.type[i]) < 2) 1 else as.numeric(substr(ColSideFactors.type[i],2,20))
				x.lwd <- if(nchar(ColSideFactors.type[i]) < 2) 1 else as.numeric(substr(ColSideFactors.type[i],2,20))
				x.lim <- NULL
				if (length(ColSideFactors.ylim) && any(names(ColSideFactors.ylim)==colnames(ColSideFactors)[i])) {
					x.lim <- ColSideFactors.ylim[[colnames(ColSideFactors)[i]]]
				}
				plot(1:nc,
					(if (x.ty == "t" || any(is.na(as.numeric(f)))) 1:length(colInd) else as.numeric(f))[colInd], 
					type="n",frame.plot=ColSideFactors.frame[i],axes=FALSE,xlab="",ylab="",main="", ylim=x.lim)
				#par(usr=c(0.5,nc+0.5,-1,1))
				if (x.ty == "t" || is.factor(f)) {
					par(usr=c(0.5,nc+0.5,-1,1))
					x.cex <- 2.05
					#xs <- which.max(strwidth(trim(as.character(f))))
					#mxnc <- max(nchar(trim(as.character(f))))
					#while ((strheight("X",units="user",cex=x.cex)*strwidth("W",units="user",cex=x.cex)/strheight("W",units="user",cex=x.cex) > 1.2 
					#	|| max(strwidth(trim(f[xs]),units="user",cex=x.cex)) > 6*ColSideFactors.height[i]/0.2) && x.cex > 0.3) 

					#try(while (strheight("[](){}0",units="user",cex=x.cex-0.05)*mxnc*0.8 > 2 ||
					#	   strwidth("0",units="user",cex=x.cex-0.5)*1/0.8 > 1 &&
					#	x.cex > 0.3) { x.cex = x.cex - 0.05 })
					pin <- par("pin")
					x.cex <- adjustCEXinches(as.character(f), width=pin[2], height=pin[1])
					#try(while (strheight("[](){}0",units="inches",cex=x.cex-0.05)*length(f) > .95*pin[1] ||
					#	   strwidth(as.character(f[xs]),units="inches",cex=x.cex-0.5) > .95*pin[2] &&
					#	x.cex > 0.1) { x.cex = x.cex - 0.02 })
					#for (lr in 1:nc)
					#	text(lr, 0, labels=f[colInd][lr], type=x.ty, col=kol[f.n][colInd][lr], pch=x.pch, lwd=x.lwd, cex=x.cex, srt=90)
					#text(1:nc, 0, labels=f[colInd], col=kol[colInd], cex=x.cex, srt=90)
					text(1:nc, 0, labels=f[colInd], col=kol[f.n.col][colInd], cex=x.cex, srt=90)
				} else if (x.ty == "b") {
					grueso <- as.numeric(substr(ColSideFactors.type[i],2,100))
					if (!is.finite(grueso)) grueso <- 0.9
					if (is.null(x.lim)) try(par(usr=c(0.5,nc+0.5,min(f,na.rm=TRUE)-ifelse(x.ty %in% c("p","h","t"),min(f,na.rm=TRUE)*0.1,0),max(f,na.rm=TRUE))))
					else try(par(usr=c(0.5,nc+0.5,par("usr")[3:4])))
					symbols(1:nc,f[colInd]/2,rectangles=data.matrix(data.frame(w=abs(grueso),h=abs(f[colInd]))),bg=kol[f.n.col][colInd],fg=kol[f.n.col][colInd], add=TRUE, inches=FALSE)
				} else {
					if (is.null(x.lim)) try(par(usr=c(0.5,nc+0.5,min(f,na.rm=TRUE)-ifelse(x.ty %in% c("p","h","t"),min(f,na.rm=TRUE)*0.1,0),max(f,na.rm=TRUE))))
					else try(par(usr=c(0.5,nc+0.5,par("usr")[3:4])))
					points(1:nc, f[colInd], type=x.ty, col=kol[f.n.col][colInd], pch=x.pch, lwd=x.lwd, cex=x.cex)
				}
			}
			if (ColSideFactors.sep[i] != 0) {
				abline(v=seq(par("usr")[1],par("usr")[2],length.out=nc+1), col=ColSideFactors.sep[i], lwd=0.5)
			}
			axis(4, at=mean(par("usr")[3:4]), labels=colnames(ColSideFactors)[i], las=1, cex.axis=ColSideFactors.axis[i])
			if (csfl =="right") {
				os <- ofl
				for (j in 1:nl) {
					uj <-  ifelse(j==1 && nchar(u[j]) == 0, 1, u[j])
					ujx <- ifelse(j==1 && nchar(u[j]) == 0, "NA", u[j]) 
					s <- trim(if (ColSideFactors.counts[i]) paste(ujx,tb[uj],sep=":") else ujx)
					axis(at=mean(par("usr")[3:4]), labels=s, side=4, line=os, col.axis=kol[j], las=1, tick=FALSE)
					os <- os + ifelse(nchar(s) < 2, .75, ifelse(nchar(s) > 2, nchar(s) / 3 + 0.75, 1))
				}
			}
		}
	}
	par(mar = c(margins[1], 0, 0, margins[2]))
	if (!symm || scale != "none") {
		x <- t(x)
		cellnote <- t(cellnote)
	}
	if (revC) {
		iy <- nr:1
		ddr <- rev(ddr)
		x <- x[, iy]
		cellnote <- cellnote[, iy]
	}
	else iy <- 1:nr
	
	
	###################################
	try(
		image(1:nc, 1:nr, x, xlim = 0.5 + c(0, nc), ylim = 0.5 + 
			c(0, nr), axes = FALSE, xlab = "", ylab = "", col = col, 
			breaks = breaks, ...)
		)

		
	if (!is.invalid(na.color) & any(is.na(x))) {
		mmat <- ifelse(is.na(x), 1, NA)
		image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "", 
			col = na.color, add = TRUE)
	}

	if (is.function(overlap)) {
		overlap(rowInd, colInd)
	}
	###################################

	cexRowFixed <- !missing(cexRow)
	if (!is.finite(cexRow) || cexRow < 0) {
		cexRow <- if (!is.finite(cexRow)) 3 else -cexRow
		#while (cexRow > 0.01 & strheight("[](){}GgJj",units="user",cex=cexRow) > 0.85) 
		#	cexRow = cexRow - 0.01
		xin <- par("fin")-par("pin")
		#win <- par()$pin
		#ml <- which.max(nchar(labRow))[1]
		#while (cexRow > 0.02 & (strheight("[](){}GgJj",units="inches",cex=cexRow)*length(labRow) > win[2] || strwidth(labRow[ml],units="inches",cex=cexRow) > 0.925*xin[1]))
		#	cexRow = cexRow - 0.01
		cexRow <- adjustCEXinches(labRow, width=xin[1], height=par("pin")[2], hf=.95, font=3)*0.90 #0.90 porque es italicas
	}
	if (!is.finite(cexCol)  ||  cexCol < 0) {
		cexCol <- if (!is.finite(cexCol)) 3 else -cexCol
		#while (cexCol > 0.01 && (strheight("[](){}X",units="user",cex=cexCol)*strwidth("W",units="user",cex=cexCol)/strheight("[](){}W",units="user",cex=cexCol) > 1.2)) 
		#	cexCol = cexCol - 0.01
		xin <- par("fin")-par("pin")
		#win <- par()$pin
		#ml <- which.max(nchar(labCol))[1]
		#while (cexCol > 0.02 && (strwidth(labCol[ml],units="inches",cex=cexCol) > 0.9*xin[2] || strheight(labCol[ml],units="inches",cex=cexCol)*length(labCol) > win[1]) ) 
		#	cexCol = cexCol - 0.01
		cexCol <- adjustCEXinches(labCol, width=xin[2], height=par("pin")[1])
	}

	if (length(labCol.col) < 2 && (is.null(labCol.col)  ||  labCol.col == FALSE  ||  (is.null(ColvFactor)  && is.logical(labCol.col))))
		axis(1, 1:nc, labels = labCol, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexCol)
	else {
		if (is.logical(labCol.col)) {
			#labCol.col <- as.factor(ColSideFactors[,ifelse(is.null(ColvFactor), 1, ColvFactor)])
			f <- getquasifactor(ColSideFactors[,ifelse(is.null(ColvFactor), 1, ColvFactor)])
			labCol.col <- as.numeric(if (is.numeric(f)) factor(f) else f)
		}
		labCol.col <- labCol.col[colInd]
		u <- unique(labCol.col)
		for (j in 1:length(u)) {
			w <- which(labCol.col == u[j])
			axis(1, (1:nc)[w], labels = labCol[w], las = 2, line = -0.5, tick = 0, 
				cex.axis = cexCol, col.axis=labCol.col[w[1]])
		}
	}
	if (!is.null(xlab)) 
		mtext(xlab, side = 1, line = margins[1] - 1.25)
	maxwidth = abs(diff(par("usr")[1:2])) * par("mai")[4] * 0.9 / par("pin")[1]
	#if (!cexRowFixed) {
	#	jcex = cexRow
	#	while (quantile(unlist(sapply(labRow, strwidth, cex=jcex)),.95) > maxwidth && jcex > 0.3333) {
	#		jcex = jcex - 0.1
	#	}
	#	cexRow = jcex
	#}
	if (length(labRow.col) < 2 && (is.null(labRow.col)  || labRow.col == FALSE  ||  (is.null(RowSideFactors) && is.logical(labRow.col))))
		axis(4, iy, labels = labRow, las = 2, line = -0.5, tick = 0, 
			cex.axis = cexRow, font=3)
	else {
		if (is.logical(labRow.col)) {
			#labRow.col <- as.factor(RowSideFactors[,ifelse(is.null(RowvFactor), 1, RowvFactor)])
			f <- getquasifactor(RowSideFactors[,ifelse(is.null(RowvFactor), 1, RowvFactor)])
			labRow.col <- as.numeric(if (is.numeric(f)) factor(f) else f)
		}
		labRow.col <- labRow.col[rowInd]
		u <- unique(labRow.col)
		for (j in 1:length(iy)) {
			w <- which(labRow.col == u[j])
			axis(4, iy[w], labels = labRow[w], las = 2, line = -0.5, tick = 0, 
				cex.axis = cexRow, col.axis=labRow.col[w[1]], font=3)
		}
	}
	if (!is.null(ylab)) 
		mtext(ylab, side = 4, line = margins[2] - 1.25)
	if (!missing(add.expr)) 
		eval(substitute(add.expr))
	if (!missing(colsep)  && !is.null(colsep)) 
		for (csep in colsep) {
			xc <- sepcolor[min(length(sepcolor),csep)]
			if (is.numeric(xc) && xc < 0) xc <- labCol.col[csep]
			rect(xleft = csep + 0.5, ybottom = rep(0, 
			length(csep))+0.5, xright = csep + 0.5 + sepwidth[1]*0, 
			ytop = rep(ncol(x) + 1, length(csep))+0.5, lty = 1, lwd = seplwd, 
			col = xc, border = xc)
		}
	if (!missing(rowsep)  && !is.null(rowsep)) 
		for (rsep in rowsep) {
			xc <- sepcolor[min(length(sepcolor),rsep)]
			if (is.numeric(xc) && xc < 0) xc <- labRow.col[rsep]
			rect(xleft = 0, ybottom = rsep -1+((ncol(x) + 
			1 - rsep) - 0.5)*0 + 0.5, xright = nrow(x) + 1*0 +0.5, ytop = rsep-1+ ((ncol(x) + 
			1 - rsep) - 0.5 - sepwidth[2])*0 + 0.5, lty = 1, lwd = seplwd, 
			col = xc, border = xc)
		}
	min.scale <- min(breaks)
	max.scale <- max(breaks)
	x.scaled <- scale01(t(x), min.scale, max.scale)
	if (trace %in% c("both", "column")) {
		for (i in colInd) {
			if (!is.null(vline)) {
				vline.vals <- scale01(vline, min.scale, max.scale)
				abline(v = i - 0.5 + vline.vals, col = linecol, 
				  lty = 2)
			}
			xv <- rep(i, nrow(x.scaled)) + x.scaled[, i] - 0.5
			xv <- c(xv[1], xv)
			yv <- 1:length(xv) - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (trace %in% c("both", "row")) {
		for (i in rowInd) {
			if (!is.null(hline)) {
				hline.vals <- scale01(hline, min.scale, max.scale)
				abline(h = i + hline, col = linecol, lty = 2)
			}
			yv <- rep(i, ncol(x.scaled)) + x.scaled[i, ] - 0.5
			yv <- rev(c(yv[1], yv))
			xv <- length(yv):1 - 0.5
			lines(x = xv, y = yv, lwd = 1, col = tracecol, type = "s")
		}
	}
	if (!missing(cellnote)) 
		text(x = c(row(cellnote)), y = c(col(cellnote)), labels = c(cellnote), 
			col = if(is.matrix(notecol)) t(notecol[rowInd,colInd]) else if (is.function(notecol)) notecol(cellnote) else notecol, cex = if(is.matrix(notecex)) t(notecex[rowInd,colInd]) else notecex)
	
	if (any(ColSideFactors.block != 0) && !is.null(ColvFactor)) {
	   	f <- ColSideFactors[,ColvFactor]
		if (is.factor(f)) f <- as.factor(as.character(f)) # to remove labels not used
		f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
		u <- levels(as.factor(f))
		nl <- length(u) # works even for numeric
		prev <- 0
		for (i in 1:nl) {
			w <- which(f == u[i])
			if (length(w) > 0) {
				#abline(v=c(prev+length(w)+0.5), col=i)
				#abline(v=c(prev+0.6,prev+length(w)+0.5), col=i)
				rect(prev+0.5,0+0.5,prev+length(w)+0.4+0.1*(i==nl),nr+0.5,border=ifelse(is.logical(ColSideFactors.block),i,ColSideFactors.block),lwd=ifelse(length(ColSideFactors.block) > 1,ColSideFactors.block[2],1.5))
				#rect(prev+0.5,0+0.5,prev+length(w)+0.5,nr+0.5,border=i,lwd=1.5)
				prev <- prev + length(w)
			}
		}
	}
	if (RowSideFactors.block != 0 && !is.null(RowvFactor)) {
	   	f <- RowSideFactors[,RowvFactor]
		if (is.factor(f)) f <- as.factor(as.character(f)) # to remove labels not used
		f.n <- as.numeric(if (is.numeric(f)) factor(f) else f)
		u <- levels(as.factor(f))
		nl <- length(u) # works even for numeric
		prev <- 0
		for (i in 1:nl) {
			w <- which(f == u[i])
			if (length(w) > 0) {
				rect(0+0.5,prev+0.5,nc+0.475,prev+length(w)+0.475+0.1*(i==nl),border=ifelse(is.logical(RowSideFactors.block),i,RowSideFactors.block),lwd=ifelse(length(RowSideFactors.block) > 1,RowSideFactors.block[2],1.5))
				prev <- prev + length(w)
			}
		}
	}
	
	.dd <- ddc
   	if (class(ddc) != "list") {
   		.dd <- list(ddc)
   		attr(.dd, "kolor") <- 1
   	}
	if (dendrogram %in% c("both", "column")) {
		for (i in 1:length(.dd)) {
			
			if (is.logical(.dd[[i]])  || is.null(.dd[[i]])  ||  is.na(.dd[[i]])) {
				par(mar=c(0,0,0,0))
				plot.new()
			}
			else {
				par(mar = c(0, 0, if (!is.null(main)) 3 else 0, if (i==length(.dd)) margins[2] else 0))
				plot(.dd[[i]], axes = FALSE, xaxs = "i", leaflab = "none", edgePar=list(col=attr(.dd, "kolor")[i]))
			}
			if (i==1 &&  !is.null(main))
				title(main, cex.main = 1.5 * op[["cex.main"]])
		}
	}
	else {
		#if (length(.dd) != 1  ||  !is.null(.dd[[1]]))
			lapply(.dd,  function(x) { par(mar=c(0,0,0,0)); plot.new() })
	  #title(main, cex.main = 1.5 * op[["cex.main"]])
	}

	par(mar = c(margins[1], 0, 0, 0))
   	.dd <- ddr
   	if (class(ddr) != "list") {
   		.dd <- list(ddr)
   		attr(.dd, "kolor") <- 1
   	}
	if (dendrogram %in% c("both", "row")) {
		for (i in 1:length(.dd)) {
			#par(mar = c(0, 0, if (!is.null(main)) 3 else 0, if (i==length(.dd)) margins[2] else 0))
			par(mar=c(if(i==1) margins[1] else 0,0,0,0))
			if (is.logical(.dd[[i]])  || is.null(.dd[[i]])  ||  is.na(.dd[[i]])) plot.new()
			else plot(.dd[[i]], horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", edgePar=list(col=attr(.dd, "kolor")[i]))
		}
		#lapply(.dd, function(ddr) plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none"))		
	}
	else lapply(.dd,  function(x) { par(mar=c(0,0,0,0)); plot.new() } )
	
	if (key) {
		par(mar = c(3, 3, 1, 1), cex = 0.75)
		if (symkey) {
			max.raw <- max(abs(x), na.rm = TRUE)
			min.raw <- -max.raw
		}
		else {
			min.raw <- min(x, na.rm = TRUE)
			max.raw <- max(x, na.rm = TRUE)
		}
		z <- seq(min.raw, max.raw, length = length(col))
		image(z = matrix(z, ncol = 1), col = col, breaks = breaks, 
			xaxt = "n", yaxt = "n")
		par(usr = c(0, 1, 0, 1))
		lv <- pretty(breaks)
		if (length(lv) < 5) {
			lv <- sort(c(lv,(lv[-length(lv)]+lv[-1])/2))
		}
		xv <- scale01(as.numeric(lv), min.raw, max.raw)
		if (sum(!is.finite(xv))==0) axis(1, at = xv, labels = lv)
		if (scale == "row") 
			mtext(side = 1, "Row Z-Score", line = 2, cex=density.cex)
		else if (scale == "column") 
			mtext(side = 1, "Column Z-Score", line = 2, cex=density.cex)
		else mtext(side = 1, "Value", line = 2, cex=density.cex)
		if (density.info == "density") {
			dens <- density(x, adjust = densadj, na.rm = TRUE)
			omit <- dens$x < min(breaks) | dens$x > max(breaks)
			dens$x <- dens$x[-omit]
			dens$y <- dens$y[-omit]
			dens$x <- scale01(dens$x, min.raw, max.raw)
			lines(dens$x, dens$y/max(dens$y) * 0.95, col = denscol, 
				lwd = 1)
			axis(2, at = pretty(dens$y)/max(dens$y) * 0.95, pretty(dens$y), cex.axis=density.cex)
			title(density.title, cex.main=density.cex)
			par(cex = 0.5)
			mtext(side = 2, "Density", line = 2, cex=density.cex)
		}
		else if (density.info == "histogram") {
			h <- hist(x, plot = FALSE, breaks = breaks)
			hx <- scale01(breaks, min.raw, max.raw)
			hy <- c(h$counts, h$counts[length(h$counts)])
			lines(hx, hy/max(hy) * 0.95, lwd = 1, type = "s", 
				col = denscol)
			axis(2, at = pretty(hy)/max(hy) * 0.95, pretty(hy), cex.axis=density.cex)
			title(density.title, cex.main=density.cex)
			par(cex = 0.5)
			mtext(side = 2, "Count", line = 2, cex=density.cex)
		}
		else title(density.title, cex.main=density.cex)
	}
	else plot.new()
	invisible(list(rowInd = rowInd, colInd = colInd))
}



#string.contain <- function(pattern, x, ...)
#dice que strings x contienen el patron pattern
#si se quiere la posición, usar regrexpr mejor que hace eso
string.contain <- function(pattern, x, ..., ret.pos=FALSE) {
	xg <- grep(pattern, x, ...)
	ok <- logical(length(x))
	if (length(xg) > 0) ok[xg] <- TRUE
	ok
	
}


# list.extract the field xname from xlist
list.extract <- function(xlist, xname, unlist=TRUE, parse=NULL) {
	if ((is.null(parse) && regexpr("\\$|\\[",xname) > 0) || parse==TRUE) {
		code <- if (regexpr("x\\$|x\\[",xname) > 0) parse(text=xname) else parse(text=paste("x",xname,sep=""))
		r <- lapply(xlist, function(x) eval(code))				
	} else {
		r <- lapply(xlist, function(x) if (is.data.frame(x) || is.matrix(x)) x[,xname] else x[[xname]])
	}
	if (unlist) return (unlist(r))
	r
}

# list.bind, extract matrices or dataframes joining them in rows or columns
list.bind <- function(xlist, byrow=1) {
	x <- NULL
	for (i in 1:length(xlist)) {
		if (byrow==1) x <- rbind(x, xlist[[i]])
		else x <- cbind(x, xlist[[i]])
	}
	x
}





## http://www.jstatsoft.org/v11/c01/
## from venn package that could not been installed
venn <- function(id, category, cutoff=1, duplicates = FALSE, tab, func=function(x) sum(x) >= cutoff, main, labels=FALSE,
		labels.cex=.5, labels.width=40, radii=1, proportional=FALSE, labels.list=FALSE, list.numbers=FALSE,
		square.plot=TRUE) {
	
	if (missing(tab)) {
		# Create incidence table from id and category		
		tab <- incidence.table(id, category, cutoff=cutoff, 
					   duplicates = duplicates, func=func)
		if (missing(main)) main <- paste("Count of", deparse(substitute(id)), 
		 				 "by", deparse(substitute(category)))
	} else if (missing(main)) 
		main <- paste("Venn diagram of", deparse(substitute(tab)))
	
	# Convert rows to binary numbers and count them
	
	index <- tab %*% 2^(1:ncol(tab) - 1)
	itab <- table(index)
	
	save <- par(pty=ifelse(square.plot,"s","m"),	# Must be square to get labels right 
			mar=c(.5, 0, 1, 0)*par("mar")) # Don"t need side margins

	on.exit(par(save))
	
	if (ncol(tab) == 2) {
		# Set up coordinates.  xlim and ylim must be same length or coordinate
		# system will not be isometric
		plot(1, 1, xlim=c(-1.3, 2.3), ylim=c(-1.8, 1.8), bty="n", axes=FALSE, 
			 type="n", xlab="", ylab="", main=main)
		if (!is.na(zero <- itab[as.character(0)])) 
			title(sub=paste(zero, "not shown"))
		# Plot 2 circles
		cx <- c(0, 1.1)
		cy <- c(0, 0)
		mx <- mean(cx)
		my <- mean(cy)
		factor <- c(radii, radii)
		if (proportional) {
			s1 <- sum(tab[,1])
			s2 <- sum(tab[,2])
			factor <- sqrt((c(s1,s2)/max(s1,s2))/pi)
			factor <- factor/max(factor)
		}
		symbols(cx, cy, circles=factor, inches=FALSE, add=TRUE)
		# Put counts in the regions
		text(c(mx + 2*(cx[1] - mx)*sqrt(factor[1]), mx + 2*(cx[2] - mx)*sqrt(factor[2])),
			 c(my + 2*(cy[1] - my), my + 2*(cy[2] - my)),
			 itab[as.character(c(1, 2))])
		text(mx+ifelse(factor[1] > factor[2],1,-1)*0.75*(cx[2] - mx)*sqrt(min(factor)), my, itab["3"])
		# Label the circles
		# text(c(cx[1] -1, cx[2] + 1), c(0, 0),		 
		#	 pos=c(2, 4), colnames(tab))	
		text(c((mx + 3*(cx[1] - mx) + cx[1])/2,	  (mx + 3*(cx[2] - mx) + cx[2])/2),
			 c((my + 3*(cy[1] - my) + cy[1] -1.8)/2, (my + 3*(cy[2] - my) + cy[2] -1.8)/2),
			 pos = c(2, 4), colnames(tab))		 
	 } else if (ncol(tab) == 3) {
		# Set up coordinates.  xlim and ylim must be same length or coordinate
		# system will not be isometric
		plot(1, 1, xlim=c(-1.5, 2.6), ylim=c(-1.5, 2.6), bty="n", axes=FALSE, 
			 type="n", xlab="", ylab="", main=main)
		if (!is.na(zero <- itab[as.character(0)])) 
			mtext(paste(zero, "not shown"), side=1)
		# Plot 3 circles
		cx <- c(0, 1.1, 0.55)
		cy <- c(0, 0, 1.1*sqrt(3)/2)
		mx <- mean(cx)
		my <- mean(cy)
		symbols(cx, cy, circles=rep(radii, 3), inches=FALSE, add=TRUE)
		# Put counts in the regions
		text(c(mx + 2*(cx[3] - mx), mx + 2*(cx[1] - mx), mx + 2*(cx[2] - mx)),
			 c(my + 2*(cy[3] - my), my + 2*(cy[1] - my), my + 2*(cy[2] - my)),
			 itab[as.character(c(1, 2, 4))], col=c(1,2,4))
		text(c(mx + (cx[1] + cx[3] -2*mx), mx + (cx[2] + cx[3] -2*mx), 	
			   mx + (cx[2] + cx[1] -2*mx)),
			 c(my + (cy[1] + cy[3] -2*my), my + (cy[2] + cy[3] -2*my), 	
			   my + (cy[2] + cy[1] -3*my)),
			 itab[as.character(c(3, 5, 6))], col=c(3,5,6))
		text(mx, my, itab["7"], col=7)
		# Label the circles
		text(c(mx + 2.8*(cx[3] - mx), (mx + 3*(cx[1] - mx) + cx[1])/2,	  
			   (mx + 3.3*(cx[2] - mx) + cx[2])/2)*radii,
			 c(my + 2.8*(cy[3] - my), (my + 3*(cy[1] - my) + cy[1] - 1.2)/2, 
			   (my + 3.3*(cy[2] - my) + cy[2] - 1.2)/2)*radii,
			 pos = c(3, 2, 4), colnames(tab))
		## VT:
		tres <- length(rownames(index)[which(index[,1]==3)])
		cinco <- length(rownames(index)[which(index[,1]==5)])
		seis <- length(rownames(index)[which(index[,1]==6)])
		siete <- c(3,5,6)[which.min(c(tres,cinco))[1]] #tres,cinco,seis
		## :VT
		x <- c(mx , -1.6,-1.6,   2.7, 2.7,  mx, mx)*radii
		y <- c(2.5,-1.25, 1.75, -1.25, 1.75,-1.5,-1.5)*radii
		adjx <- c(.5, 0, 0, 1, 1, 0.5, 0.5)
		## VT:
		x[7] <- x[siete]
		y[7] <- 1.1#y[siete]
		adjx[7] <- adjx[siete]
		## :VT
		listnames <- list()
		for (i in 1:7) {
			#3,5,6 : 7
			listnames[[i]] <- rownames(index)[which(index[,1]==i)]
			rn <- paste(listnames[[i]],collapse=", ")
			if (labels  && !labels.list) {
				cr <- trunc(nchar(rn) / labels.width + 0.999)
				for (j in 1:cr) {
					rn <- paste(substr(rn, 1, labels.width*j+j-1),
								"\n",
								substr(rn,labels.width*j+j,nchar(rn)),sep="")
				}
				text(x[i], y[i], rn, col=i, cex=labels.cex, adj=c(adjx[i],0.5))
			}
		}
		if (labels && labels.list) {
			len <- unlist(lapply(listnames,length))
			if (list.numbers) {
				ul <- paste(unlist(sapply(len, function(i) if(i > 0) 1:i else NULL)), ":",unlist(listnames),sep="")				
			} else {
				ul <- unlist(listnames)
			}
			col <- rep(1:7, times=len)
			i <- 1
			j <- 2
			while (j <= length(ul)) {
				if (nchar(ul[j])+nchar(ul[i]) <= labels.width && col[j] == col[i]) {
					ul[i] <- paste(ul[i], ", ", ul[j], sep="")
					ul[j] <- NA
					col[j] <- NA
				} else {
					i <- j
				}
				j <- j + 1
			}
			ul <- ul[!is.na(ul)]
			col <- col[!is.na(col)]
			pri <- 1:trunc(length(ul)/2)
			sec <- (length(pri)+1):length(ul)
			h <- strheight(ul[1],cex=labels.cex*1.2)
			text(-1.6, seq(my+h*length(pri)/2, my-h*length(pri)/2, length.out=length(pri)), ul[pri], col=col[pri], adj=c(0,0.5), cex=labels.cex)
			text( 2.7,  seq(my+h*length(sec)/2, my-h*length(sec)/2, length.out=length(sec)), ul[sec], col=col[sec], adj=c(1,0.5), cex=labels.cex)
		}
		l <- colnames(tab)
		names(listnames) <- c(l[1], l[2],  
							paste(l[1],l[2],sep="+"),
							l[3],
							paste(l[1],l[3],sep="+"),
							paste(l[2],l[3],sep="+"),
							paste(l[1],l[2], l[3],sep="+"))
		invisible(listnames)
	} else stop("Can only Venn 2 or 3 categories")
}

incidence.table <- function(id, category, names = NULL, cutoff = 1, 
				duplicates = FALSE, func = function(x) sum(x) >= cutoff) {
   	if (is.matrix(id)  || is.data.frame(id)) {
   		if (is.data.frame(id)) {
   			cl <- c()
   			for (i in 1:ncol(id)) cl <- class(id[,i])
   			id <- data.matrix(id)
   			if (all(cl) == "logical") id <- id > 0
   		}
   		if ((is.matrix(id) && mode(id) != "logical")) {
   			id <- id >= cutoff   			
   		}
   		#id is a LOGICAL MATRIX, category is a vector of column classes
		e <- matrix(FALSE, nrow=nrow(id), ncol=nlevels(category))
		colnames(e) <- levels(category)
		rownames(e) <- rownames(id)
		for (i in 1:nrow(id)) {
			k = 0
			for (j in levels(category)) {
				k = k + 1
				w <- which(category==j)
				e[i,k] <- func(id[i,w])
			}
		}
		e
   	} else if (!duplicates) {
		# Count combinations and convert to TRUE/FALSE
		tab <- table(as.character(id), category)
		tab >= cutoff
	} else {
		# Count combinations
		tab <- table(as.character(id), category)
		# Set up matrix with one row per id
		result <- matrix(FALSE, length(id), ncol(tab))
		# Set appropriate entries TRUE
		for (i in 1:ncol(tab)) 
			result[, i] <- tab[as.character(id), i] >= cutoff
		# Return nice looking matrix
		rownames(result) <- as.character(names)
		colnames(result) <- colnames(tab)
		result
	}
}

plot.chars <- function(cex=1,from=1,to=130,col=from:to,pch=from:to,cols=10,rows=NULL,...) {
	d <- (to-from+1)
	if (missing(rows)) rows <- round(d/cols)
	if (missing(cols)) cols <- round(d/rows)
	df <- data.frame(x=rep(1:cols,times=rows)[1:d], y=rep(1:rows,each=cols)[1:d],pch=circularize(pch,to-from+1),col=circularize(col,to-from+1),cex=circularize(cex,to-from+1),stringsAsFactors=FALSE)
	with(df, plot(x,y,pch=pch,cex=cex,col=col,...))
}

# Draw sample stacked points similar to those in publications
# x.round should be a vector rounded to obtain equal values or breaked/cut has an histogram
# class is the class of each value
# x.levels is the levels perhaps sorted by user specified order
plot.vertical.densities <- function(x.round, cls, pch=c(20,17,18,15,1:25), col=1:nl, ylim=c(min(x.round),max(x.round)), xlim=c(0,nl), x.levels=levels(cls),xlab="Groups",nperline=10,cex=1,...) {
	nl <- length(x.levels)
	plot(0,0,xlim=xlim,ylim=ylim,type="n",xaxt="n",xlab=xlab,...)
	for (i in 1:nl) {
		x <- x.round[cls==x.levels[i]]
		xt <- table(x)
		xn <- as.numeric(names(xt))
		for (j in 1:length(xt)) {
			s <- paste(rep("@",xt[j]),collapse="")
			w <- strwidth(s,cex=cex)
			h <- strheight(s,cex=cex)
			z <- seq(-w/4,w/4,length.out=xt[j])
			points(rep(i-0.5,xt[j]) + if (xt[j] > 1) z else 0, rep(xn[j],xt[j]),pch=pch[i],col=col[i],cex=cex)
		}
		lines(c(i-0.5-cex/4,i-0.5+cex/4),c(mean(x),mean(x)), col=col[i], lwd=2*cex,lty=3)
		axis(1, i-0.5, paste(x.levels[i], "\n(n=", table(cls)[x.levels[i]],")",sep=""), col=col[i], col.axis=col[i],...)
	}
	
}


plot.contour <- function(x, y, z, xbreaks="Sturges", ybreaks="Sturges", ...) {
	h <- hist(x, plot=FALSE, breaks=xbreaks)
	v <- hist(y, plot=FALSE, breaks=ybreaks)
	m <- matrix(0, nrow=length(h$mids), ncol=length(v$mids))
	colnames(m) <- v$mids
	rownames(m) <- h$mids
	x1 <- h$breaks[1]
	for (i in 2:length(h$breaks)) {
		x2 <- h$breaks[i]
		y1 <- v$breaks[1]		
		for (j in 2:length(v$breaks)) {
			y2 <- v$breaks[j]
			m[i-1,j-1] <- sum(x > x1 & x <= x2 & y > y1 & y <= y2)
			y1 <- y2
		}
		x1 <- x2
	}
	contour(h$mids, v$mids, m, ...)
	m
}



rad2deg <- function(rad) {
	rad*180/pi
}

deg2rad <- function(deg) {
	deg*pi/180
}

polar2xy <- function(z, rad) {
	if (length(rad) < length(z)) rad <- rep(rad, len=length(z))
	if (length(z) < length(rad)) z <- rep(z, len=length(rad))
	matrix(z * c(cos(rad), sin(rad)), ncol=2)
}

# m - matrix, columns are variables, rows are samples
# radius = 1 (all variables are scaled)
#rx <- matrix(runif(30, max=100), ncol=5)
#rownames(rx) <- paste("Row ",1:nrow(rx))
#colnames(rx) <- paste("Col ",1:ncol(rx))
#plot.radar(rx)
#plot.radar(rx, density=c(2,1,1,1,0.5)*30, fill=1:5)
plot.radar <- function(m, alfa.off=90, minpoint=0.5, col=1:s, pch=1:s, lty=rep(1,s), lwd=rep(1,s), type=rep("l",s),
	xlab="", ylab="", main="", density=rep(NULL,s), fill=rep(NULL,s), angle=seq(15, 180-15, len=s), fill.lwd=rep(1,s),
	ticks.deg=4, digits=4, mar=c(2,1,3,1), rlab=rownames(m), clab=colnames(m), wx=0.5, wy=0.5, legend.cex=1, legend.cols=1, ...) {
	rn <- rlab
	cn <- clab
	n <- ncol(m)
	s <- nrow(m)
	alfa <- 360 / n # angle for each variable
	mx <- apply(m, 2, max)
	mn <- apply(m, 2, min)
	m <- apply(m, 2, function(x) ((x-min(x))/abs(max(x)))*(1-minpoint)+minpoint)
	angles <- seq(0, 360, length.out=n+1)[1:n]+alfa.off
	rad <- deg2rad(angles)
	cxy <- polar2xy(1, rad)
	up <- par("mar")
	on.exit(par(up))
	par(mar=mar)
	plot.new()
	f <- 1.075
	if (!is.null(cn)) f <- 1.25
	par(usr=c(min(cxy[,1]), max(cxy[,1]), min(cxy[,2]), max(cxy[,2]))*f)
	#plot(0,0,type="n",xlim=c(-0.95,0.95), ylim=c(-0.95,0.95), xaxt="n", xlab=xlab, ylab=ylab, main=main, yaxt="n")
	#abline(h=0)
	#abline(v=0)
	a <- pretty(m)
	#symbols(rep(0,length(a)), rep(0,length(a)), circles=a, fg=8, inches=FALSE, lty=3)
	rn.l <- length(rn)
	col <- circularize(col,rn.l)
	lty <- circularize(lty,rn.l)
	pch <- circularize(pch,rn.l)
	for (i in 1:length(a)) {
		xy <- polar2xy(a[i], rad)
		xy2 <- rbind(xy, xy[1,])
		lines(xy2, col=8, lty=3)
	}
	
	for (i in 1:n) {
		b <- polar2xy(max(a), rad[i])
		lines(c(0,b[1,1]), c(0,b[1,2]), col=1)
		if (!is.null(cn)) {
			b <- polar2xy(max(a)*1, rad[i])
			text(b[1,1], b[1,2], cn[i], col=1, adj=c(wx-cos(rad[i]), wy-sin(rad[i])))
		}
		b <- polar2xy(a, rad[i]-deg2rad(1))
		b2 <- polar2xy(a, rad[i]+deg2rad(1))
		b3 <- polar2xy(a, rad[i]-deg2rad(ticks.deg))
		for (j in 1:nrow(b)) {	
			lines(c(b[j,1], b2[j,1]), c(b[j,2],b2[j,2]), col=1)
			#text(b3[j,1], b3[j,2], format(a[j], digits=6), cex=0.66, adj=c(0.5,0.5))
			text(b3[j,1], b3[j,2], format((a[j]-minpoint)*(mx[i]-mn[i])/(1-minpoint)+mn[i], digits=digits), cex=0.66, adj=c(0.5,0.5))
		}
	}
	
	for (i in 1:s) {
		xy <- polar2xy(m[i,], rad)
		xy2 <- rbind(xy, xy[1,])
		for (j in 1:n) {
			lines(xy2, col=col[i], pch=pch[i], lty=lty[i], lwd=lwd[i], type=type[i])
		}
		polygon(xy2, border=col[i], lty=lty[i], lwd=fill.lwd[i], density=density[i], angle=angle[i], col=fill[i])
	}
	if (!is.null(rn)) {
		ph <- ifelse(type != "l", pch, "")
		plot.legend("topleft",rn,col=col,text.col=col,pch=ph,lty=lty,lwd=lwd,cex=legend.cex, ncol=legend.cols, ...)
	}
}


#FILTERING
# filter.by.variation <- function(data.n, data.sd=NULL, data.m=NULL, draw=TRUE, pch=20, cex=0.5, use=c("cv","sd"), 
# 	filter.th = 0.1, ...) 
#FILTERING
filter.by.variation <- function(data.n, data.sd=NULL, data.m=NULL, 
	draw=TRUE, pch=20, cex=0.5, ratios=FALSE, use=c("cv","sd"), 
	filter.th = 0.1,  percentage=FALSE, 
	minq=0.025, maxq=0.975, min.m=-Inf, max.m=+Inf, min.sd=0, max.sd=+Inf, min.cv=0, max.cv=+Inf,
	...) {
	if (ratios && missing(use)) use <- "sd"
	use = match.arg(use)
	if (is.null(data.sd)) data.sd <- apply(data.n, 1, sd, na.rm=TRUE)
	if (is.null(data.m)) data.m <- abs(apply(data.n, 1, mean, na.rm=TRUE))
	data.cv <- data.sd/data.m
	use.cv <- (use=="cv")
	var2use <- if(use.cv) data.cv else data.sd
	
	xp <- par()
	# FILTERING BY SD
	if (draw) {
		on.exit(par(xp))
		par(mfrow=c(2,3), mar=c(3,2,2,1))
		plot.densities(data.n, overall=TRUE, main="Original Data", xlab="", ylab="")
		plot.densities(qminmax(data.m,minq,maxq), overall=TRUE, xlab="Mean", ylab="", main="Mean")
		plot.densities(qminmax(var2use,minq,maxq), overall=TRUE, xlab=use, ylab="", main=use)		
		plot(qminmax(data.m,minq,maxq), qminmax(var2use,minq,maxq), pch=pch, cex=cex, main=paste("Mean vs",use), xlab="Mean", ylab=use, ...)
	}
	
	# FILTERING BY CV
	if (ratios) {
		xv <- sqrt(data.m^2+var2use^2)
		xa <- asin(var2use/xv) # rad
		data.loess <- loess(v~a, data.frame(a=xa, v=xv))
		s <- seq(min(xa), max(xa), length.out=100)
		p <- predict(data.loess, data.frame(a=s))
		lp <- predict(data.loess, data.frame(a=xa))
		if (percentage) {
			f <- var2use-lp
			filter.th <- quantile(f, 1-filter.th, na.rm=TRUE)
		}
		f <- xv - filter.th > lp & 
			(data.m >= min.m & data.m <= max.m) &
			(data.sd >= min.sd & data.sd <= max.sd) &
			(data.cv >= min.cv & data.cv <= max.cv) 
		s2 <- p * cos(s)
		p2 <- p * sin(s)
		if (draw) {
			lines(s2, p2, col=2, lwd=2)
			lines(-s2, p2, col=2, lwd=2)
			par(usr=c(0,max(xa)*2,-max(xv),max(xv)))
			symbols(max(xa)/2,max(xv)/2,rectangles=data.matrix(data.frame(x=max(xa),y=max(xv))), add=TRUE, inches=FALSE, col=4)
			points(xa, xv, pch=".", col=4)
			lines(s, p, col=2, lwd=1)
		}
		s <- c(s2,-s2)
		o <- order(s)
		s <- s[o]
		p <- c(p2, p2)[o]
		#filter.th : inscrease this to filter out more genes, decrease to filter in MORE
		data.f <- data.n[f,]
	} else {
		w <- if (length(data.m) > 10000) sort(sample(1:length(data.m),10000)) else 1:length(data.m)
		data.loess <- loess(v~m, data.frame(m=data.m, v=var2use)[w,])
		s <- seq(min(data.m, na.rm=TRUE), max(data.m, na.rm=TRUE),length.out=100)
		p <- predict(data.loess, data.frame(m=s))
		if (draw) {
			lines(s, p, col=2, lwd=2)
		}
		#filter.th : inscrease this to filter out more genes, decrease to filter in MORE
		lp <- predict(data.loess, data.frame(m=data.m))
		if (percentage) {
			f <- var2use-lp
			filter.th <- quantile(f, 1-filter.th, na.rm=TRUE)
		}
		f <- var2use - filter.th > lp & 
			(data.m >= min.m & data.m <= max.m) &
			(data.sd >= min.sd & data.sd <= max.sd) &
			(data.cv >= min.cv & data.cv <= max.cv)
		data.f <- data.n[f,]
		
	}
	if (draw) {
		plot(qminmax(data.m,minq,maxq), qminmax(var2use,minq,maxq), pch=pch, cex=cex, col=ifelse(f,1,8), main=paste(sum(f, na.rm=TRUE),"of",length(f),"Passing Filter=",round(filter.th,4)," (in black)"), xlab="Mean", ylab=use, ...)
		lines(s, p, col=2, lwd=2)
		plot.densities(qminmax(data.f,minq,maxq), overall=TRUE, main=paste("Passing Filter (",filter.th,")",sep=""), xlab="", ylab="")
	}
	
	invisible(list(filter=f, data=data.f))
}



plot.pca <- function(data, classes=factor(1:nrow(data)), col=1:nlevels(classes), npc=4, gap=0, pch=1:nlevels(classes), labels=rownames(data), legend.ncol=2, center=TRUE, scale.=TRUE, ...) {
	pca <- prcomp(data, center=center, scale.=scale.)
	ik <- unclass(classes)
	colo <- col
	pcho <- pch
	colfixed <- FALSE
	pchfixed <- FALSE
	if (length(col) != nrow(data)) {
		col <- col[ik]
		colfixed <- TRUE
	}
	if (length(pch) != nrow(data)) {
		pch <- pch[ik]
		pchfixed <- TRUE
	}
	upPanel <- function(x,y,...) {
		#points(x,y,pch=pch[ik],col=col[ik], ...)
		points(x,y,pch=pch,col=col, ...)
	}
	lowPanel <- function(x,y,...) {
		par(usr=c(range(y)+abs(range(y))*0.15*c(-1.15,1.15),range(x)+abs(range(x))*0.15*c(-1.15,1.15)))
		#text(y,x,labels,col=col[ik], adj=c(0,0),...)
		text(y,x,labels,col=col, adj=c(0,0),...)
	}
	.j <- 0
	txtPanel <- function(...) {
		j <- .j <<- .j + 1		
		usr <- par("usr"); on.exit(par(usr))
		par(usr = c(0,1,0,1) )
		text(0.5,0.1,paste("Comp ",j,"\n",round(pca$sdev[j]^2*100/sum(pca$sdev^2),3),"%\nCum=",round(sum(pca$sdev[1:j]^2)*100/sum(pca$sdev^2),3),"%",sep=""), adj=c(0.5,0))
		legend("top",legend=levels(classes),col=if(colfixed) colo else NA,pch=if(pchfixed) pcho else NA,ncol=legend.ncol)
	}
	if (npc < 1) {
		npc <- which(cumsum(pca$sdev^2)/sum(pca$sdev^2) > npc)[1]
	}
	pairs(pca$x[,1:min(npc,ncol(pca$x))], gap=gap, upper.panel=upPanel,
		lower.panel=if(!is.null(labels)) lowPanel else upPanel, text.panel=txtPanel, ...)
	invisible(pca)
}



# remove parts of the string that are repeated in a "fraction" (default=100%=all) of names
simplify.names <- function(x, pattern="[^A-Za-z0-9]", fraction=1, ...) {
	s <- strsplit(x, pattern, ...)
	u <- unique(unlist(s))
	k <- data.matrix(data.frame(sapply(u, function(i) lapply(s, function(j) sum(unique(j) == i)))))
	w <- which(apply(k, 2, sum) >= fraction*length(x))
	if (length(w) > 0) {
		n <- u[w]
		z <- unlist(lapply(s, function(i) paste(i[!(i %in% n)], collapse=".")))
	} else {
		z <- x
	}
	z
}

plot.string.rot <- function(x,y,xx,yy, stv, col=rep(1,length(stv)), srt=0, cex.min=.1, cex.max=3, word.sep=" ", plot=TRUE, adj=c(0, 0), f.adj=2, space.x=1, space.y=1.3, rect=FALSE, debug=FALSE, f.overlap=0.5) {
	r <- c(min(x,xx), max(x,xx), min(y, yy), max(y,yy))
	x <- r[1]
	xx <- r[2]
	y <- r[3]
	yy <- r[4]
	if (rect) {
		symbols(x+(xx-x)/2,y+(yy-y)/2,rectangles=matrix(c(xx-x, yy-y),ncol=2), add=TRUE, inches=FALSE)
	}
	col <- circularize(col, length(stv))
	pst <- paste(stv, collapse=word.sep)
	pstch <- character(nchar(pst))
	for (i in 1:nchar(pst)) pstch[i] <- substr(pst, i, i)
	colch <- rep(col, times=nchar(stv)+nchar(word.sep))
	pin <- par("pin")
	usr <- par("usr")
	xpin2usr <- function(xin) {
		usr[1] + (xin / pin[1]) * (usr[2]-usr[1])
	}
	ypin2usr <- function(yin) {
		usr[3] + (yin / pin[2]) * (usr[4]-usr[3])
	}
	xusr2pin <- function(x) {
		(x-usr[1]) / (usr[2]-usr[1]) * pin[1]
	}
	yusr2pin <- function(y) {
		(y-usr[3]) / (usr[4]-usr[3]) * pin[2]
	}
	xin <- xusr2pin(x)
	xxin <- xusr2pin(xx)
	yin <- yusr2pin(y)
	yyin <- yusr2pin(yy)
	inh <- pin[2]
	inw <- pin[1]
	inx <- (xx-x) * inw / (usr[2]-usr[1])
	iny <- (yy-y) * inh / (usr[4]-usr[3])
	cex.inf <- cex.min
	cex.sup <- cex.max
	a <- c(	abs(cos(srt*2*pi/360)),
			abs(cos((90-srt)*2*pi/360)), 
			abs(sin(srt*2*pi/360)), 
			abs(sin((90-srt)*2*pi/360)))
	area <- 0
	repeat {
		area0 <- area
		cex <- (cex.inf+cex.sup) / 2
		if (word.sep == "\n") {
			w <- max(inx, max(unlist(sapply(stv, strwidth, cex=cex, units="inches"))))
			h <- max(unlist(sapply(stv, strheight, cex=cex, units="inches")))
		} else {
			w <- strwidth(pst, cex=cex, units="inches") 
			h <- strheight(pst, cex=cex, units="inches")
			if (w > inx) {
				h <- h * w / inx
				w <- inx
			} else if (h > iny) {
				w <- w * h / iny
				h <- iny
			}
		}
		wr <- w #w * a[1] + h * a[2]
		hr <- h #w * a[3] + h * a[4]
		area <- sum(hr)*sum(wr) #*space.y*space.x
		if (debug) 
			cat("area=",area,",space=",inx*iny*f.adj,",sum(hr)=",sum(hr),",sum(wr)=",sum(wr),",cex=",cex,"\n")
		if (debug) cat("inx=",inx,", iny=",iny,",w=",w,",h=",h,"\n")
		# || strwidth("W",cex=cex,units="inches") > inx  || strheight("j",cex=cex,units="inches") > iny
		if (area*f.adj > (inx*iny)) {
			cex.sup <- cex
		} else {
			cex.inf <- cex
		}
		if (cex < cex.min || cex > cex.max || (abs(cex.sup - cex.inf) < 0.05)  ||  abs(area-area0) < 0.01) {
			break
		}
	}
	cex <- max(min((cex.inf+cex.sup) / 2, cex.max), cex.min)
	w <- sapply(pstch, strwidth, cex=cex, units="inches")
	h <- sapply(pstch, strheight, cex=cex, units="inches")
	wr <- (w * a[1] + h * a[2])
	hr <- (h * a[4] + w * a[3])
	x0 <- xin + h[1] * a[2] * 1
	if (srt > 75) {
		y0 <- yin
	} else {
		y0 <- yyin - max(h) * a[4] * space.y
	}
	xi <- x0
	yi <- y0
	jump <- TRUE
	topedy <- FALSE
	if (debug) cat("cex=",cex,",xin=",xin,",yin=",yin,",xxin=",xxin,",yyin=",yyin,"\n")
	if (plot) {
		for (i in 1:length(pstch)) {
			dx <- w[i] * a[1]
			dy <- w[i] * a[3]
			if (debug) cat("w[i]=",(w[i]),", h[i]=",(h[i]),", dy=",(dy),", dx=",(dx),"\n")
			#cat("wr[i]=",wr[i],", hr[i]=",hr[i],", dy=",dy,", dx=",dx,"\n")
			if ((yi+hr[i]*f.overlap > yyin) || (xi+w[i]*a[1]*f.overlap > xxin) || pstch[i] == "\n") {
				if (y0 >= yin  || x0 <= xxin) {
					jump <- TRUE
					yi <- y0 - (mean(h) * a[4] * space.y + mean(w) * a[3] * space.x) 
					xi <- x0
					if (yi < yin && xi < xxin) {
						yi <- yin
						xi <- x0 + (mean(w) * a[1] * space.x + mean(h) * a[2] * space.y) 
					}
					y0 <- yi
					x0 <- xi
				}
				#if (x0 > xin) {
				#	xin <- x0 + dx
				#}
				#if (y0 > yyin) {
				#	yyin <- y0 + dy
				#}
				#y <- min(y, yi-max(hr)*space.y/2)
			}
			if ((pstch[i] != " " && pstch[i] != "\n")  ||  !jump) {
				if (debug) cat("i=",i, ", xi=",(xi),", yi=",(yi),", ch=",pstch[i],", y=",y,", y0=",y0,",x0=",x0,"\n")
				text(xpin2usr(xi), ypin2usr(yi), pstch[i], col=colch[i], srt=srt, cex=cex, adj=adj)
				xi <- xi + dx
				yi <- yi + dy
				jump <- FALSE
				if (debug) cat("i=",i, ", xi=",(xi),", yi=",(yi),", ch=",pstch[i],", y=",y,", y0=",y0,",x0=",x0,"\n")
			}
		}
	}
	invisible(cex)
}


inschat <- function(st, ch="\n", at=80) {
	for (i in 1:length(st)) {
		if (nchar(st[i]) > at) {
			x <- nchar(st[i])/at
			for (j in 1:x) {
				st[i] <- paste(substr(st[i], 1, at*j+j-1),ch,substr(st[i],at*j+j,nchar(st[i])),sep="")
			}
		}
	}
	st
}


minmax <- function(x, min=0, max=1, quantiles=FALSE) {
	if (quantiles) {
		min = quantile(x, min, na.rm=TRUE)
		max = quantile(x, max, na.rm=TRUE)
	}
	x[x < min] <- min
	x[x > max] <- max
	x
}

qminmax <- function(x, min=0.01, max=0.99) {
	minmax(x, min, max, TRUE)
}

plot.text <- function(x, space=5, max.chars=40, autoscale=missing(cex), cex=1, precision.factor=1, col=NULL, box.col=NULL, add=FALSE, x0=NULL, y0=NULL, xx=NULL, yy=NULL, formatter=format, hlines=8, vlines=8, hadj=NULL, family=NULL, xlab="", ylab="", do.trim=TRUE, row.names=FALSE, mar=par("mar"), ...) {
	addedcols <- FALSE
	if (row.names && !is.null(rownames(x))) {
		x <- cbind("."=as.character(rownames(x)), x)
		addedcols <- TRUE
	}
	xmar <- par("mar")
	# x should be a matrix
	if (is.null(hadj)) {
		hadj <- rep(0,ncol(x))
		for (i in 1:ncol(x)) {
			if (any(class(x[,i]) %in% c("numeric","integer","double"))) {
				hadj[i] <- 1
			#} else if (any(class(x[,i]) %in% c("factor"))) {
			#	hadj[i] <- 0.5
			}
		}
	}
	for (i in 1:ncol(x)) {
		x[,i] <- formatter(x[,i])
	}
	xtrim <- if (do.trim) trim else { function(x) x }
	x <- rbind(colnames(x), x)
	txt <- apply(x, 2, function(y) {
		unlist(sapply(y,function(x) paste(substr(xtrim(formatter(x)),1,max.chars),ifelse(nchar(xtrim(formatter(x))) > max.chars, "..."," "),sep="")))		})
	nc <- apply(txt, 2, nchar)
	mnc <- pmin(max.chars+3, apply(nc, 2, max))
	cs <- c(1,cumsum(mnc)+space)[1:ncol(x)]
	x2 <- sum(mnc)-1
	x1 <- 1
	#y1 <- -0.25
	#y2 <- nrow(x)
	y1 <- -0.25 
	y2 <- nrow(txt)
	if (is.null(col)) {
		col <- matrix(c(4,rep(1,nrow(x)-1)), ncol=ncol(x), nrow=nrow(x))
		if (addedcols) { col[,1] <- 6; col[1,1] <- 0 }
	}
	if (!is.matrix(col) && !is.data.frame(col)) {
		if (length(col) == ncol(x)) 
			col <- matrix(col, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
		else if (length(col) == nrow(x)-1)
			col <- matrix(c(4,col), ncol=ncol(x), nrow=nrow(x))
		else {
			col <- matrix(col, ncol=ncol(x), nrow=nrow(x))
			col[1,] <- 4
			if (addedcols) { col[,1] <- 6; col[1,1] <- 0 }
		}
	}
	if ((is.matrix(col) || is.data.frame(col)) && (nrow(col) == nrow(x)-1)) {
		col <- rbind(rep(1,ncol(col)),col)
	}
	if (!is.null(box.col)  &&  !is.matrix(box.col) && !is.data.frame(box.col)) {
		if (length(box.col) == ncol(x)) 
			box.col <- matrix(box.col, ncol=ncol(x), nrow=nrow(x), byrow=TRUE)
		else if (length(box.col) == nrow(x)-1)
			box.col <- matrix(c(0,box.col), ncol=ncol(x), nrow=nrow(x))
		else 
			box.col <- matrix(box.col, ncol=ncol(x), nrow=nrow(x))
	}
	if (!is.null(box.col) && 
		(is.matrix(box.col) || is.data.frame(box.col)) && (nrow(box.col) == nrow(x)-1)) {
		box.col <- rbind(rep(1,ncol(box.col)),box.col)
	}
	prefam <- par("family")
	if (add) {
		xp <- par("usr")
		x0 <- if (is.null(x0)) xp[1] else x0
		xx <- if (is.null(xx)) xp[2] else xx
		y0 <- if (is.null(y0)) xp[4] else y0
		yy <- if (is.null(yy)) xp[3] else yy
		dx <- (x2-x1)*(xp[2]-xp[1])/(xx-x0)
		xb <- dx*(xp[2]-x0)/(xp[2]-xp[1])+x1
		xa <- xb - dx
		dy <- (y2-y1)*(xp[4]-xp[3])/(yy-y0)
		yb <- dy*(xp[4]-y0)/(xp[4]-xp[3])+y1
		ya <- yb - dy
		par(usr=c(xa, xb, ya, yb), family=if (is.null(family)) prefam else family, mar=mar)
	} else {
		par(mar=mar)
		plot(0,0,type="n",xlim=c(x1,x2*1.05), ylim=c(y2, y1), xlab=xlab, ylab=ylab, family=if (is.null(family)) prefam else family, ...)
		xp <- par("usr")
		x0 <- xp[1]
		xx <- xp[2]
		y0 <- xp[3]
		yy <- xp[4]
	}
		
	txt.w <- apply(txt, 2, function(x) abs(unlist(sapply(x, strwidth, cex=cex, units="user"))))
	txt.h <- apply(txt, 2, function(x) abs(unlist(sapply(x, strheight, cex=cex, units="user"))))
	txt.mw <- apply(txt.w, 2, which.max)
	txt.mh <- apply(txt.h, 2, which.max)
	#txt.area <- txt.w * txt.h
	#txt.ma <- apply(txt.area, 2, which.max)
	#wmx <- which.max(apply(txt, 1, function(x) sum(unlist(sapply(x, strwidth, cex=cex, units="user")))))[1]
	#ejemplo <- txt[wmx,]
	ejemplo <- unlist(txt[data.matrix(data.frame(r=txt.mw, c=1:ncol(txt)))])
	if (autoscale) {
		maxwrow <- which.max(apply(txt.w, 1, sum))[1]
		#cexw <- adjustCEXinches(paste(t(txt[maxwrow,]),collapse="  "))
		#cexw <- adjustCEXinches(paste(ejemplo,collapse=" "), width=par("pin")[1]/(ncol(txt)+row.names*1))
		cexw <- adjustCEXuser(paste(substr(ejemplo,1,nchar(ejemplo)-1),collapse=""), width=abs(x2-x1)) #/(ncol(txt)+row.names*1)
		maxhrow <- which.max(apply(txt.h, 1, max))[1]
		#cexh <- adjustCEXinches(paste(t(txt[maxhrow,]),collapse="  "), height=par("pin")[2]/nrow(txt))
		#cexh <- adjustCEXinches(paste(ejemplo,collapse=" "), height=par("pin")[2]/nrow(txt))
		cexh <- adjustCEXuser(paste(ejemplo,collapse=""), height=abs(y2-y1)/nrow(txt))
		cex <- min(cexh, cexw)
		if (FALSE) {
			#while (TRUE) {
			#	mx <- sum(unlist(sapply(ejemplo, strwidth, cex=cex, units="user")))
			#	if (mx*precision.factor >= x2) break
			#	cex <- cex * 1.5
			#}
			while (TRUE) {
				mx <- sum(unlist(sapply(ejemplo, strwidth, cex=cex, units="user")))
				if (mx*precision.factor < x2) break
				cex <- cex - 0.01
			} 
			ejemplo <- unlist(txt[data.matrix(data.frame(r=txt.mh, c=1:ncol(txt)))])
			while (TRUE) {
				mx <- max(abs(unlist(sapply(ejemplo, strheight, cex=cex, units="user"))))
				if (mx <= 1) break
				cex <- cex - 0.01
			}
		}
	}
	lw <- apply(txt, 2, function(x) max(sapply(x, strwidth, units="user", cex=cex)))
	lw <- lw+(x2-sum(lw))*lw/sum(lw)
	cs2 <- cumsum(c(x1,lw))
	cs <- cs2[1:ncol(x)]
	if (!add) {
		if (!is.null(family)) par(family=family)
		#plot(0,0,type="n",xlim=c(x1,max(x2,cs2[ncol(x)+1])), ylim=c(y2, y1), xlab=xlab, ylab=ylab, family=if (is.null(family)) prefam else family, ...)
		par(usr=c(x1*0.75,max(x2,cs2[ncol(x)+1])*1.005,y2*1.005,y1))
	}
	#if (!add) abline(h=0.01, col=8, lty=1, lwd=0.5)
	#if (hlines) lines(c(x1, x2),c(0.01, 0.01), col=8, lty=1, lwd=0.5)
	#if (!add) abline(v=cs, col=8, lty=1, lwd=0.5)
	#if (all(hadj==0))
	#	text(cs+0.25, 0, txt[1,], adj=c(0,0), cex=cex, col=col[1,], ...)
	
	#ys <- seq(y1, y2, length=nrow(x)+1)
	#dy <- 0*mean(diff(ys))/2
	mh <- max(abs(unlist(sapply(ejemplo, strheight, cex=cex, units="user"))))*1.5
	for (i in 0:nrow(txt)) {
		#lines(c(x1,cs2[ncol(x)+1]), c(i-mh,i-mh)+1, col=hlines, lty=1, lwd=0.5)
		lines(c(x1,cs2[ncol(x)+1]), c(i,i), col=hlines, lty=1, lwd=0.5)
		if (!is.null(box.col) && i < nrow(x)) {
			points(cs2[ncol(x)+1]+strwidth("."), i+0.5, pch=20, col=box.col[i+1,ncol(x)])
			points(cs2[1]-strwidth("."), i+0.5, pch=20, col=box.col[i+1,1])
		}
		#lines(c(x1,x2), c(ys[i+1]+dy,ys[i+1]+dy), col=hlines, lty=1, lwd=0.5)
	}

	for (i in 1:ncol(x)) {
		lines(c(cs2[i],cs2[i]), c(0, nrow(x)), col=vlines, lty=1, lwd=0.5)
		if (!is.null(box.col)) {
			symbols(rep((cs2[i+1]+cs2[i])/2,nrow(x)), 1:nrow(x)-0.5,
				rectangles=data.matrix(data.frame(width=rep(lw[i],nrow(x)),height=0.95)), 
				fg=box.col[,i], add=TRUE, inches=FALSE)
		}
		text(ifelse(hadj[i]==0,cs[i],ifelse(hadj[i]==1,cs2[i+1],(cs2[i+1]+cs2[i])/2))+0.25*0,
			 1:nrow(x), txt[,i], adj=c(hadj[i],0), cex=cex, col=col[,i], ...)
	}
	#lines(c(x2,x2), c(-mh, nrow(x)-mh), col=vlines, lty=1, lwd=0.5)
	lines(rep(cs2[ncol(x)+1],2), c(0, nrow(x)), col=vlines, lty=1, lwd=0.5)

	#for (i in 1:nrow(x)) {
	#	lines(c(x1,x2), c(i-0.9,i-0.9), col=hlines, lty=1, lwd=0.5)
	#	for (j in 1:ncol(x)) {
	#		text(ifelse(hadj[j]==0,cs[j],ifelse(hadj[j]==1,cs2[j+1],(cs2[j+1]+cs[j])/2))+0.25, i-1, txt[i,j], adj=c(hadj[j],0), cex=cex, col=col[i,j], ...)
	#	}
	#}
	
	if (add) par(usr=xp)
	if (!is.null(family)) par(family=prefam)
	par(mar=xmar)
	
	if (autoscale) {
		return (cex)
	}
}

plot.text.by.pages <- function(dframe, pageSize=50, numerate=TRUE, run.in.each.page=function(dframe, ...) { 0 }, ...) {
	for (page in 1:ceiling(nrow(dframe) / pageSize)) {
		from <- pageSize*(page-1)+1
		to   <- min(nrow(dframe),from+pageSize-1)
		dframeOK <- dframe[from:to,]
		if (numerate) { plot.text(data.frame(Num=from:to, dframeOK), ...) }
		else { plot.text(dframeOK, ...) }
		if (is.function(run.in.each.page)) run.in.each.page(dframeOK, ...)
	}
}



plot.string <- function(x0=NULL, y0=NULL, xx=NULL, yy=NULL, x, sep=ifelse(rotation,"  "," "), max.chars=40,  precision.factor=1.2, col=NULL, add=TRUE, wrap=TRUE, formatter=format, cex.min=.1, cex.max=10.1, autoscale=missing(cex), cex=cex.max, vfactor=1.25, srt=c(0,90,180,270)[1], adj=c(0,0), do.trim=TRUE, ...) {
	# x should be a character vector
	if (length(x) == 0) return (NULL)
	x2 <- sapply(x, formatter)
	if (!is.null(names(x))) x2 <- paste(names(x),"=",x2,sep="")
	xp <- par("usr")
	x0 <- if (is.null(x0)) xp[1] else x0[1]
	xx <- if (is.null(xx)) xp[2] else xx[1]
	y0 <- if (is.null(y0)) xp[3] else y0[1]
	yy <- if (is.null(yy)) xp[4] else yy[1]
	if (y0 > yy) {
		.y <- y0
		y0 <- yy
		yy <- .y
	}
	if (x0 > xx) {
		.x <- x0
		x0 <- xx
		xx <- .x
	}
	if (!add) {
		plot(0,0,type="n",xlim=c(x0,xx), ylim=c(nrow(x)-1, -0.5), xlab="", ylab="", ...)
	}
	#lines(c(x0,xx),c(y0,y0))
	#lines(c(x0,xx),c(yy,yy))
	#lines(c(x0,x0),c(y0,yy))
	#lines(c(xx,xx),c(y0,yy))

	rotation <- FALSE
	if (srt == 90  || srt == 270) {
		rotation <- TRUE
		a <- x0
		x0 <- y0
		y0 <- a
		a <- xx
		xx <- yy
		yy <- a
	}
	
	dx <- abs(xx-x0)*ifelse(rotation || srt == 180,.99995,1)
	dy <- abs(yy-y0)*ifelse(rotation || srt == 180,.99995,1)
	
	xtrim <- if (do.trim) trim else { function(x) x }
	if (xx > x0) {
		txt <- sapply(x2, function(y) {
		paste(substr(xtrim(y),1,max.chars),ifelse(nchar(xtrim(y)) > max.chars, "...",""),sep="")})
	} else {
		txt <- sapply(x2, function(y) {
		paste(substr(xtrim(y),1,max.chars),ifelse(nchar(xtrim(y)) > max.chars, "...",""),sep="")})
	}
	if (is.null(col)) {
		col <- rep(1,length(txt))
	}
	if (length(col) != length(txt)) col <- rep(col,length(txt))[1:length(txt)]
	#txt <- paste(sep,txt,sep="")
	txt2 <- rep(sep,length(txt)*2-1)
	txt2[1:length(txt)*2-1] <- txt
	txt <- txt2
	col2 <- rep(8,length(col)*2-1)
	col2[1:length(col)*2-1] <- col
	col <- col2
	rows <- 1

	#if (autoscale) {
		cex.sup <- cex.max
		cex.inf <- cex.min
		cexok <- NULL
		pos <- matrix(0,nrow=6,ncol=length(txt))
		while (TRUE) {
			#cat(cex,"\n")
			txt.w <- unlist(sapply(txt, strwidth, cex=cex, units="user"))
			txt.h <- unlist(sapply(txt, strheight, cex=cex, units="user"))
			if (rotation) {
				u <- par("usr")
				.y <- abs(u[4]-u[3])
				.x <- abs(u[2]-u[1])
				inch <- par("pin")		
				txt.w <- txt.w * (.y/.x) * (inch[1]/inch[2])
				txt.h <- txt.h * (.x/.y) * (inch[2]/inch[1])
			}
			pos[5,] <- txt.w
			pos[6,] <- txt.h
			y <- 0
			x <- 0
			rows <- 1
			overflow <- FALSE
			posok <- pos
			for (i in 1:length(txt)) {
				#cat(c(x,txt.w[i],xx),"\n")
				if (wrap) {
					ok <- 1
					if (x+txt.w[i] > dx && i > 1) {
						#pos[2,i-1] <- xx
						x <- 0
						y <- y + max(txt.h)*vfactor
						pos[4,i-1] <- y0+y
						rows <- rows + 1
						ok <- ifelse(txt[i] == sep,0,1)
					} else {
						if (i > 1) pos[2,i-1] <- x0+x
					}
					pos[1:4,i] <- c(x0+x,xx,y0+y,yy)
					x <- x + txt.w[i]*ok
					if (x > dx) overflow <- TRUE
				} else {
					for (j in 1:nchar(txt[i])) {
						wc <- strwidth(substr(txt[i],j,j),cex=cex,units="user")
						if (x+wc > dx) {
							x <- 0
							y <- y + max(txt.h)*vfactor
						}
						pos[,i] <- c(x,x+txt.w[i],y,y+txt.h[i])
						x <- x + wc
					}
				}
			}
			y <- y + max(txt.h)*min(1,vfactor)
			if (y > dy || x > dx || overflow) {
				cex.sup <- cex
			} else {
				posok <- pos
				cexok <- cex
				cex.inf <- cex
			}
			if (abs(cex.sup-cex.inf) < 0.01) {
				if (is.null(cexok)) cexok <- cex
				break
			}
			cex <- (cex.sup+cex.inf) / 2
		}
	#}
	cex <- cexok
	pos <- posok
	txt.w <- unlist(sapply(txt, strwidth, cex=cex, units="user"))
	txt.h <- unlist(sapply(txt, strheight, cex=cex, units="user"))
	yadj <- ifelse(yy > y0, ifelse(srt==270 || srt == 180,1,0), ifelse(srt==270 || srt == 180,0,1))
	xadj <- ifelse(xx > x0, ifelse(srt==90 || srt == 180,1,0), ifelse(srt==90 || srt == 180,0,1))
	ysign <- ifelse(yy > y0, 1, -1)
	xsign <- ifelse(xx > x0, 1, -1)
	y <- y0 #ifelse(yy > y0, y0, yy)
	x <- x0 #ifelse(xx > x0, x0, xx)
	sx <- 0
	if (rotation) {
		u <- par("usr")
		.y <- abs(u[4]-u[3])
		.x <- abs(u[2]-u[1])
		inch <- par("pin")		
		txt.w <- txt.w * (.y/.x) * (inch[1]/inch[2])
		txt.h <- txt.h * (.x/.y) * (inch[2]/inch[1])
	}
	ajuste <- if (rotation) c(yadj, xadj) else c(xadj,yadj) #ifelse(rotation,c(yadj,xadj),c(xadj,yadj))
	#cat("CEX=",cex,"------\n")
	xcoord <- numeric(length(txt))
	ycoord <- numeric(length(txt))
	if (adj[2] == 1) {
		y <- y+max(0,dy-max(txt.h)*rows*vfactor*ysign)
	} else if (abs(adj[2]-0.5) < 0.01) {
		y <- y+max(0,(dy-max(txt.h)*rows*vfactor*ysign)/2)
	} else {
		y <- y + 0
	}
	for (i in 1:length(txt)) {
		#cat(c(x,y,paste("\"",txt[i],"\"",sep=""),txt.w[i],txt.h[i],xx,yy),"\n")
		if (wrap) {
			ok <- 1
			if (sx+txt.w[i] > dx && i > 1) {
				sx <- 0
				x <- x0
				y <- y + max(txt.h)*vfactor*ysign
				ok <- ifelse(txt[i] == sep,0,1)
			}
			if (adj[1] == 1) {
				xcoord[i] <- ifelse(rotation,yy,xx)-txt.w[i]
			} else if (abs(adj[1]-0.5) < 0.01) {
				xcoord[i] <- (ifelse(rotation,pos[3,i]+pos[4,i],x+pos[2,i])-txt.w[i])/2
			} else {
				xcoord[i] <- ifelse(rotation,y,x)
			}
			ycoord[i] <- ifelse(rotation,x,y)
			text(xcoord[i], ycoord[i], txt[i], cex=cex, col=col[i], adj=ajuste,srt=srt)
			x <- x + txt.w[i]*xsign*ok
			sx <- sx + txt.w[i]*ok
		} else {
			jc <- 1:nchar(txt[i])
			if (xx < x0) jc <- rev(jc)
			for (j in 1:nchar(txt[i])) {
				wc <- strwidth(substr(txt[i],jc[j],jc[j]),cex=cex,units="user")
				ok <- 1
				if (sx+wc > dx) {
					sx <- 0
					x <- x0
					y <- y + max(txt.h)*vfactor*ysign
					ok <- ifelse(txt[i] == sep,0,1)
				}
				text(ifelse(rotation,y,x), ifelse(rotation,x,y), substr(txt[i],jc[j],jc[j]), cex=cex, col=col[i], adj=ajuste, srt=srt)
				x <- x + wc*xsign*ok
				sx <- sx + wc*ok
			}
		}
	}
	invisible(list(x=xcoord,y=ycoord,cex=cex))
}


p.cox <- function(data, time, status, method="breslow", value=c("logtest","sctest","waldtest","rsq","coefficient")[1], attr=c("pvalue", "rsq")[1]) {
	coef <- (value=="coefficient")
	pcox <- numeric(nrow(data))
	library(survival)
	for (i in 1:nrow(data)) {
		icox <- NULL
		try(icox <- coxph(Surv(time, status) ~ ., data.frame(t(data[i,,drop=FALSE])), method=method))
		if (!is.null(icox)) {
			if (coef)
				pcox[i] <- summary(icox)$coefficients[1,5]
			else
				pcox[i] <- summary(icox)[[value]][attr]
		} else pcox[i] <- 1
	}
	pcox
}

## just used to calibrate widths of chars
draw.width.chars <- function(pch=1,cex=1,from=2,to=from,n=11) {
	plot(0,0,type="n",xlim=c(0,n),ylim=c(0,n))
	w <- strwidth(pch,cex=cex)
	title(paste("reported width=",w))
	s <- seq(w/from,w*to,length=n)
	for (i in 1:n) {
		points(cumsum(rep(s[i],n)),rep(i,n),pch=pch,cex=cex, col=1:n)
	}
	text(par("usr")[2],1:n,paste("f=",round(s/w,4)),adj=c(1,0),cex=0.75,col=2)
	text(par("usr")[2],1:n,paste("w=",round(s,4)),adj=c(1,1),cex=0.75,col=4)
}


## just used to calibrate heights of chars
draw.height.chars <- function(pch=1,cex=1,from=2,to=from,n=11) {
	plot(0,0,type="n",ylim=c(0,n),xlim=c(0,n))
	h <- strheight(pch,cex=cex)
	title(paste("reported height=",h))
	s <- seq(h/from,h*to,length=n)
	for (i in 1:n) {
		points(rep(i,n),cumsum(rep(s[i],n)),pch=pch,cex=cex, col=1:n)
	}
	text(1:n,par("usr")[4],paste("f=",round(s/h,4)),adj=c(1,0),cex=0.75,col=2,srt=90)
	text(1:n,par("usr")[4],paste("h=",round(s,4)),adj=c(1,1),cex=0.75,col=4,srt=90)
}

## draw char symbols in circles to fill an area
draw.in.circle <- function(xcent, ycent, ipoints=100, cex=1, pch=1, col=1, xch=pch[1], anglef=8, radiusf=1, method=2, round.method=trunc) {
	w <- strwidth(xch, cex=cex)
	h <- strheight(xch, cex=cex)
	if (is.numeric(xch) && xch %in% 1:25) {
		widths <- c(0.95,  1.175, 1.325, 1.10,  1.25, 1.175, 0.95, 1.25,   1.15,  0.49,
					0.65,  0.5,   0.6,   0.55,  0.44, 0.44,  0.59, 0.4433, 0.55, 0.35,
					0.475, 0.50,  0.5,   0.6,   0.6)
		heights <- c(0.65, 0.8,   0.95,  0.9,   0.95, 0.8,   0.7,  0.95,   0.95, 0.8, 
					1.16,  0.9,   0.95,  0.85,  0.85, 0.65,  0.75, 0.65,   0.7,  0.55,
					0.7,   0.70,  0.85,  0.85,  0.85)
		w <- w * widths[xch]
		h <- h * heights[xch]
	}
	a <- w * h * 1.1 * ipoints
	r <- sqrt(a / pi)
	if (method==1) {
		rs <- seq(0, r*radiusf, length.out=ipoints)
		b <- seq(0, 2*pi*max(r/w,r/h)*anglef, length.out=ipoints)
		x <- rs * cos(b)
		y <- rs * sin(b)
	} else {
		if (ipoints > 4 || ipoints < 2) {
			x <- 0
			y <- 0
			n <- 1
		} else {
			x <- c()
			y <- c()
			n <- 0
		}
		rx <- 0
		ry <- 0
		r  <- 0
		l <- sqrt(w*w+h*h)/sqrt(2)
		lx <- w/sqrt(2)
		ly <- h/sqrt(2)
		while (n < ipoints) {
			r <- r + l
			rx <- rx + lx
			ry <- ry + ly
			if (method==3) {
				nx <- min(c(round.method(2*pi*rx / lx), round.method(2*pi*ry / ly), ipoints-n))
				b <- seq(0, 2*pi*(nx-1)/nx, length.out=nx)
				x <- c(x, rx * cos(b))
				y <- c(y, ry * sin(b))
			} else {
				nx <- min(round.method(2*pi*r / l), ipoints-n)
				b <- seq(0, 2*pi*(nx-1)/nx, length.out=nx)
				x <- c(x, r * cos(b))
				y <- c(y, r * sin(b))
			}
			n <- n + nx
		}
	}
	points(xcent+x,ycent+y,cex=cex,pch=pch,col=col)
}

#draw.in.circle(4,6,anglef=0.5,ipoints=200,radiusf=.5)


plot.xy.discrete <- function(x, y, col=1, pch=19, cex=1, jitters=TRUE, 
	use.swarm="none", dec=2, cex.labels=0.75, df=2, 
	lining=min(length(unique(x)),length(unique(y))) < 20, lining.col=8, xlim=NULL, ylim=NULL, 
	xaxt=ifelse(is.factor(x),"n","s"), yaxt=ifelse(is.factor(y),"n","s"), xper=TRUE, yper=TRUE, allper=TRUE, 
	alln=FALSE, top.cex.axis=1, right.cex.axis=1, in.circles=trunc, anglef=0.5, radiusf=0.5, ADJ=c(1,0), ...) {
	mar <- par("mar")
	par(mar=c(4,4,5,4))
	on.exit(par(mar))
	j <- if(is.function(jitters)) jitters else { if(jitters) { if (is.numeric(jitters)) function(x) jitter(x, factor=jitters) else jitter } else function(x) x }
	jx <- jy <- j
	xx <- if (is.factor(x)) as.numeric(x) else x
	yy <- if (is.factor(y)) as.numeric(y) else y
	if (is.null(xlim)) xlim <- range(xx)+(c(-1,1)*ifelse(is.function(jitters), 0.5, jitters/2))
	if (is.null(ylim)) ylim <- range(yy)+(c(-1,1)*ifelse(is.function(jitters), 0.5, jitters/2))
	if (use.swarm != "none") {
		library(beeswarm)
		plot(j(xx), j(yy), pch=pch, col=col, xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt, type="n", cex=cex, ...)
		for (i in unique(x)) {
			for (k in unique(y)) {
				w <- which(xx == i & yy == k)
				if (length(w) > 0) {
					.x <- swarmx(i, yy[w])
					.y <- swarmy(xx[w], k)
					if (use.swarm=="x") {
						points(j(xx)[w], .y[,2], pch=pch[pmin(length(pch),w)], col=col[pmin(length(col),w)], cex=cex)
					} else if (use.swarm=="y") {
						points(.x[,1], j(yy)[w], pch=pch[pmin(length(pch),w)], col=col[pmin(length(col),w)], cex=cex)
					} else {
						points(.x[,1], .y[,2], pch=pch[pmin(length(pch),w)], col=col[pmin(length(col),w)], cex=cex)
					}
				}
			}
		}
	} else if (!is.null(in.circles)) {
		plot(j(xx), j(yy), pch=pch, col=col, xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt, type="n", cex=cex, ...)
		for (i in unique(x)) {
			for (k in unique(y)) {
				w <- which(xx == i & yy == k)
				if (length(w) > 0) {
					draw.in.circle(i,k,anglef=anglef,radiusf=radiusf,cex=cex,ipoints=length(w),pch=pch[pmin(length(pch),w)],col=col[pmin(length(col),w)], round.method=if (is.function(in.circles)) in.circles else trunc)
				}
			}
		}
	} else {
		plot(j(xx), j(yy), pch=pch, col=col, cex=cex, xlim=xlim, ylim=ylim, xaxt=xaxt, yaxt=yaxt, ...)
	}
	tx <- table(xx)
	tx <- tx[order(as.numeric(names(tx)))]
	axis(3, at=as.numeric(names(tx)), labels=paste(tx,"\n(",round(tx*100/sum(tx,na.rm=TRUE),dec),"%)",sep=""), cex.axis=top.cex.axis)
	ty <- table(yy)
	ty <- ty[order(as.numeric(names(ty)))]
	axis(4, at=as.numeric(names(ty)), labels=paste(ty,"\n(",round(ty*100/sum(ty,na.rm=TRUE),dec),"%)",sep=""), cex.axis=right.cex.axis)
	z <- paste(xx,yy,sep=":")
	tz <- table(z)
	dx <- median(diff(as.numeric(names(tx))))
	dy <- median(diff(as.numeric(names(ty))))
	zx <- as.numeric(unlist(lapply(strsplit(names(tz),":"), function(x) x[[1]])))
	zy <- as.numeric(unlist(lapply(strsplit(names(tz),":"), function(x) x[[2]])))
	theper <- paste(round(tz*100/sum(tz),dec),"%",sep="")
	then   <- paste("",tz,sep="")#paste("n=",tz,sep="")
	res <- if (allper && !alln) theper else if (allper & alln) paste(theper,then,sep="\n") else then
	if (allper || alln) text(zx+dx/df-strwidth(res,cex=cex.labels)*0/2,zy-dy/df,res, cex=cex.labels, adj=ADJ, srt=45*0)
	if (lining) {
		abline(v=min(xx, na.rm=TRUE):(max(xx, na.rm=TRUE)+1)-ifelse(is.function(jitters), 0.5, jitters/2), col=lining.col)
		abline(h=min(yy, na.rm=TRUE):(max(yy, na.rm=TRUE)+1)-ifelse(is.function(jitters), 0.5, jitters/2), col=lining.col)
	}
	for (i in 1:length(zx)) {
		if (xper) text(zx[i]-dx/df,zy[i],paste(round(tz[i]*100/sum(yy == zy[i], na.rm=TRUE),dec),"%",sep=""), cex=cex.labels, srt=90, adj=c(0.5,0.9))
		if (yper) text(zx[i],zy[i]+dy/df,paste(round(tz[i]*100/sum(xx == zx[i],na.rm=TRUE),dec),"%",sep=""), cex=cex.labels, adj=c(0.5,0.9))
	}
	if (missing(xaxt) && is.factor(x)) axis(1,at=1:nlevels(x),levels(x),las=3)
	if (missing(yaxt) && is.factor(y)) axis(2,at=1:nlevels(y),levels(y),las=2)
}

quantize <- function(data, n=2, min.value=min(data, na.rm=TRUE), max.value=max(data, na.rm=TRUE)) {
	quantos <- n
	delta <- (max.value - min.value)
	brk <- c(-Inf, seq(min.value,max.value,length.out=quantos)+(delta/(2*quantos-2)))
	brk[length(brk)] <- Inf
	if (is.matrix(data) ||  is.data.frame(data)) {
		rown <- rownames(data)
		coln <- colnames(data)
		data <- t(apply(t(data), 2, function(x) cut(x, labels=FALSE, breaks=brk, include.lowest=TRUE)))
		func <- function(x, ...) { 
			y <- matrix(x, ...)
			rownames(y) <- rown
			colnames(y) <- coln
			y
		}
		
	} else {
		vn <- names(data)
		data <- cut(data, labels=FALSE, breaks=brk, include.lowest=TRUE)
		func <- function(x,...) {
			names(x) <- vn
			x
		}
	}
	if (any(is.na(data))) {
		data[is.na(data)] <- quantos
	}
	func(seq(min.value,max.value,length.out=quantos)[data],ncol=ncol(data))
}

print('getAnywhere("wilcox.test.default")')


soft.read.delim <- function(xFile, row.estimate=NULL, chunk.size=1000, wc="wc -l", skip=0, header=TRUE, sep="\t", nrows=-1, comment.char="", columns=NULL) {

	if (is.null(row.estimate)) {
		row.estimate <- if (nrows > 0) nrows else as.numeric(rev(strsplit(system(paste(wc,xFile,sep=" "),intern=TRUE), " +")[[1]])[2])
	}
	cnx <- file(xFile, "r", blocking = FALSE)
	eof <- FALSE
	if (skip > 0) {
		xL <- readLines(cnx, n=skip)
		eof <- (length(xL) < skip)
	}
	d.cols <- NULL
	t.cols <- NULL
	f.cols <- NULL
	h.cols <- NULL
	read.rows <- 0
	r <- 0
	catf("Reading",xFile,"by",chunk.size,"...")
	xH <- NULL
	while (!eof && (read.rows < nrows  ||  nrows <= 0)) {
		if (r %% 10 == 0) catf("\n", read.rows, "lines read ")
		else catf(read.rows,"")
		xL <- readLines(cnx, n=if (nrows > 0) min(chunk.size, nrows-read.rows) else chunk.size)
		eof <- (length(xL) < chunk.size)
		if (length(xL) == 0) break
		if (nchar(comment.char) > 0) {
			xL <- xL[substr(xL,1,1) != comment.char]
			if (length(xL) == 0) next
		}
		
		if (is.null(d.cols) && header) {
			xH <- xL[1]
			xL <- xL[-1]
		}
		
		lines.sep <- strsplit(xL, sep)
		if (! is.null(columns)) lines.sep <- lapply(lines.sep, function(l) l[columns])
		lines.len <- unlist(lapply(lines.sep, length))
		lines.cols <- max(lines.len)
		lines.n <- length(xL)
		
		if (is.null(d.cols)) {
			## first time reading
			wcon <- file("-tmpsft-.txt","w")
			writeLines(c(xH,xL), wcon)
			close(wcon)
			xrd <- read.delim("-tmpsft-.txt", header=header, sep=sep, stringsAsFactors=FALSE)
			#file.remove("-tmpsft-.txt")
			d.cols <- list()
			t.cols <- list()
			f.cols <- list()
			h.cols <- colnames(xrd)
			for (i in 1:lines.cols) {
				xf <- switch(class(xrd[,i]), 
							numeric=numeric, character=character,
							logical=logical, complex=complex,
							integer=integer, double=double, time=time,
							factor=factor)
				xc <- switch(class(xrd[,i]), 
							numeric=as.numeric, character=as.character, 
							logical=as.logical, complex=as.complex,
							integer=as.integer, double=as.double, time=as.time,
							factor=as.factor)
				
				if (is.null(xf)) {
					xf <- character
					xc <- as.character
				}
				xd <- xf(row.estimate)
				xd[] <- NA
				d.cols[[i]] <- xd
				f.cols[[i]] <- xf
				t.cols[[i]] <- xc
			}
		}
		
		## check for missing columns
		if (lines.cols > length(d.cols)) {
			cat("Adding ",lines.cols - d.cols, " columns as character.")
			for (i in (length(d.cols)+1):lines.col) {
				d.cols[[i]] <- character(length(d.cols[[1]]))
				t.cols[[i]] <- as.character
				f.cols[[i]] <- character
			}
		}
		
		if (length(d.cols[[1]]) < read.rows+lines.n) {			
			for (i in 1:length(d.cols)) {
				xd <- f.cols[[i]](read.rows + lines.n * 10)
				xd[] <- NA
				xd[1:length(d.cols[[i]])] <- d.cols[[i]]
				d.cols[[i]] <- xd
			}
		}
		for (i in 1:lines.cols) {
			y <- character(lines.n)
			for (j in which(lines.len >= i)) {
				y[j] <- lines.sep[[j]][i]
			}
			d.cols[[i]][read.rows + 1:lines.n] <- t.cols[[i]](y)
		}
		read.rows <- read.rows + lines.n
		r <- r + 1
	}
	catf("\n", read.rows, "lines read ")
	## cols processing
	close(cnx)
	read.rows <- max(read.rows, nrows, na.omit=TRUE)
	xdf <- NULL
	for (i in 1:length(d.cols)) {
		y <- d.cols[[i]][1:read.rows]
		d.cols[[i]] <- FALSE
		if (is.null(xdf)) {
			xdf <- data.frame(y)
		} else {
			xdf[,i] <- y
		}
	}
	colnames(xdf) <- rep(h.cols,1+ncol(xdf)/length(h.cols))[1:ncol(xdf)]
	xdf
}

#x <- soft.read.delim("GPL5188.txt", nrows=100000, comment.char="#")
#x <- soft.read.delim("GPL5188.txt", comment.char="#", chunk.size=10000)




align.score <- function(seq1, seq2, gap=-3, igap=-3, radio=min(10,min(nchar(c(seq1,seq2)))/2), match=1, mismatch=-0.5, ret.scores=FALSE) {
#radio=trunc(max(nchar(c(seq1,seq2)))-min(nchar(c(seq1,seq2))))
	#library(Biostrings)
	#x <- pairwiseAlignment(seq1, seq2)
	#score <- sum(unlist(strsplit(as.character(x@pattern),"")) == unlist(strsplit(as.character(x@subject),"")))
	#return (score)
	### below is depracated
	if (nchar(seq1) < nchar(seq2)) {
		x <- seq1
		seq1 <- seq2
		seq2 <- x
	}
	s1 <- unlist(strsplit(seq1,""), use.names=FALSE)
	s2 <- unlist(strsplit(seq2,""), use.names=FALSE)
	s <- matrix(NA, ncol=length(s1)+1, nrow=length(s2)+1)
	s[,1] <- igap
	s[1,] <- igap
	s[1,1] <- 0
	d <- s*0
	for (ii in 1:length(s1)) {
		i <- ii
		for (j in 1:length(s2)) {
			m <- ifelse(s1[i] == s2[j], match, mismatch)
			if (radio > 0) {
				minx <- i:max(1,i-radio)
				miny <- j:max(1,j-radio)
				mx <- c(s[j,i]+m, s[j+1,minx]+gap*1:length(minx), s[j:miny,i+1]+gap*1:length(miny))
				s[j+1,i+1] <- max(mx)
				d[j+1,i+1] <- which.max(mx)[1]
			} else {
				s[j+1,i+1] <- s[j,i]+m
				d[j+1,i+1] <- 1
			}	
			i <- i - 1
			if (i <= 0) break
		}
	}
	iii <- 2
	for (ii in rep(length(s1),length(s2)-1)) {
		i <- ii
		for (j in iii:length(s2)) {
			m <- ifelse(s1[i] == s2[j], match, mismatch)
			if (radio > 0) {
				minx <- i:max(1,i-radio)
				miny <- j:max(1,j-radio)
				mx <- c(s[j,i]+m, s[j+1,minx]+gap*1:length(minx), s[j:miny,i+1]+gap*1:length(miny))
				s[j+1,i+1] <- max(mx)
				d[j+1,i+1] <- which.max(mx)[1]
			} else {
				s[j+1,i+1] <- s[j,i]+m
				d[j+1,i+1] <- 1
			}	
			i <- i - 1
			if (i <= 0) break
		}
		iii <- iii + 1
	}
	if (ret.scores) 
		{
			colnames(s) <- c("gap",paste(s1,1:length(s1),sep=""))
			rownames(s) <- c("gap",paste(s2,1:length(s2),sep=""))
			names(d) <- names(s)
			list(scores=s, direction=d)
		} 
	else 
		max(s[,ncol(s)],s[nrow(s),])
}



best.align <- function(seq1, seq2, ret.scores=TRUE, ret.pos=FALSE, ...) { ## ... are parameters for align.score
	#library(Biostrings)
	#x <- pairwiseAlignment(seq1, seq2)
	#return (c(as.character(x@pattern), as.character(x@subject)))

	if (nchar(seq1) < nchar(seq2)) {
		x <- seq1
		seq1 <- seq2
		seq2 <- x
	}
	as <- align.score(seq1, seq2, ..., ret.scores=TRUE)
	s <- as$scores
	d <- as$direction
	wy <- which.max(s[,ncol(s)])
	wx <- which.max(s[nrow(s),])
	w <- which.max(c(s[wy,ncol(s)],s[nrow(s),wx]))[1]
	s1 <- ""
	s2 <- ""
	if (w == 2) {
		x <- wx
		y <- nrow(s)
		if (x <= nchar(seq1)) s1 <- substr(seq1,x,nchar(seq1)) #paste(rep("-",nchar(seq1)-x),collapse="")
	} else {
		x <- ncol(s)
		y <- wy
		if (y <= nchar(seq2)) s2 <- substr(seq2,y,nchar(seq2)) #paste(rep("-",nchar(seq2)-y),collapse="")
	}
	s1 <- paste(substr(seq1,x-1,x-1),s1,sep="")
	s2 <- paste(substr(seq2,y-1,y-1),s2,sep="")
	A <- list()
	B <- list()
	while (x > 2 && y > 2) {
		#w <- which.max(c(s[y-1,x-1],s[y-1,x],s[y,x-1]))[1]
		w <- d[y-1,x-1]
		#cat(x,y,s1,s2,w,"\n")
		if (w == 1) {
			x <- x - 1
			y <- y - 1
			n1 <- substr(seq1,x-1,x-1)
			n2 <- substr(seq2,y-1,y-1)			
		} else if (w == 2) {
			x <- x - 1
			n1 <- substr(seq1,x-1,x-1)
			n2 <- "-"
			B$gaps <- c(B$gaps, y)
		} else {
			y <- y - 1
			n2 <- substr(seq2,y-1,y-1)
			n1 <- "-"
			A$gaps <- c(A$gaps, x)
		}
		s1 <- paste(n1,s1,sep="")
		s2 <- paste(n2,s2,sep="")
	}
	#cat(x,y,s1,s2,n1,n2,"\n")
	#s1 <- paste(n1,s1,sep="")
	#s2 <- paste(n2,s2,sep="")
	A$start <- 1
	B$start <- 1
	if (x <= 2 && y > 1) {
		s1 <- paste(paste(rep("-",y-2),collapse=""),s1,sep="")
		s2 <- paste(substr(seq2,1,y-2),s2,sep="")
		A$start <- y
	}
	if (y <= 2 && x > 2) {
		s1 <- paste(substr(seq1,1,x-2),s1,sep="")
		s2 <- paste(paste(rep("-",x-2),collapse=""),s2,sep="")
		B$start <- y
	}
	if (nchar(s1) < nchar(s2)) s1 <- paste(s1,paste(rep("-",nchar(s2)-nchar(s1)),collapse=""),sep="")
	if (nchar(s2) < nchar(s1)) s2 <- paste(s2,paste(rep("-",nchar(s1)-nchar(s2)),collapse=""),sep="")
	if (ret.pos) list(A=A,B=B) else c(s1,s2)
}




consensus <- function(s, remove=c("-"," "), weight=rep(1,length(s)), align=FALSE, minrelfreq=NULL, removeends=FALSE) {
	# todas las seq del mismo numero de characteres
	s <- toupper(s)
	m <- max(nchar(s))
	cons <- character(m)
	cons[] <- "/"
	for (i in 1:m) {
		xs <- substr(s,i,i)
		for (j in 1:length(remove)) xs <- gsub(remove[j],"",xs)
		xs <- rep(xs,times=weight)
		xs <- xs[nchar(xs) > 0]
		if (length(xs) > 0) {
			txs <- table(xs)
			if (is.null(minrelfreq) || max(txs) >= length(x) * minrelfreq) {
				cons[i] <- names(which.max(txs)[1])
			} else {
				cons[i] <- "N"
			}
		}
	}
	if (align) {
		aligns <- ""
		for (i in 1:m) {
			xs <- substr(s,i,i)
			xs[xs == cons[i]] <- "."
			aligns <- paste(aligns, xs, sep="")
		}
		c(paste(cons,collapse=""), aligns)
	} else {
		paste(cons,collapse="")
	}
}



oligosignature <- function(s, len=6, normalize=FALSE, debug=FALSE) {
	library(gtools)
	idx <- apply(permutations(4,len,c("A","C","G","T"), repeats.allowed=TRUE), 1, paste,collapse="")
	m <- matrix(0,nrow=length(idx),ncol=length(s))
	rownames(m) <- idx
	colnames(m) <- names(s)
	for (i in 1:length(s)) {
		s[i] <- gsub("[^ACGT]","",s[i])
		if (debug) cat(i,"/",length(s),"...")
		n <- nchar(s[i])
		xs <- sapply(1:(n-len),function(j) substr(s[i],j,j+len-1))
		xx <- table(xs)
		m[names(xx),i] <- as.numeric(xx)/ifelse(normalize,n,1)
		if (debug) cat("\n")
	}
	m
}



# http://www.ncbi.nlm.nih.gov/blast/Doc/urlapi.html
# http://www.ncbi.nlm.nih.gov/BLAST/Doc/urlapi.pdf

ncbi.blastp <- 
function (seq, database = "pdb", hitlist=1000, wordsize=3, expect=10, matrix="BLOSUM62", params="") 
{
	if (!is.vector(seq)) {
		stop("Input 'seq' should be a single sequence as a single or multi element character vector")
	}
	seq <- paste(seq, collapse = "")
	if (!(database %in% c("pdb", "nr", "swissprot"))) 
		stop("Option database should be one of pdb, nr or swissprot")
		#		 "&MATRIX_NAME=PAM30&GAPCOSTS=9,1",

	urlput <- paste("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?CMD=Put&DATABASE=", 
		database, "&HITLIST_SIZE=",hitlist,
		"&PROGRAM=blastp",
		"&WORD_SIZE=",wordsize,
		"&CLIENT=web",
		"&EXPECT=",expect,
		"&MATRIX_NAME=",matrix,
		params,
		"&QUERY=",  
		paste(seq, collapse = ""), sep = "")
	txt <- scan(urlput, what = "raw", sep = "\n", quiet = TRUE)
	writeLines(txt,"blastout.txt")
	rid <- sub("^.*RID = ", "", txt[grep("RID =", txt)])
	cat(paste(" Searching ... please wait (updates every 5 seconds) RID =", 
		rid, "\n "))
	urlget <- paste("http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get", 
		"&FORMAT_OBJECT=Alignment", "&ALIGNMENT_VIEW=Tabular", 
		"&RESULTS_FILE=on", "&FORMAT_TYPE=CSV", "&ALIGNMENTS=20000", 
		"&RID=", rid, sep = "")
	raw <- read.csv(urlget, header = FALSE, sep = ",", quote = "\"", 
		dec = ".", fill = TRUE, comment.char = "")
	html <- 1
	while (length(html) == 1) {
		cat(".")
		Sys.sleep(5)
		raw <- read.csv(urlget, header = FALSE, sep = ",", quote = "\"", 
			dec = ".", fill = TRUE, comment.char = "")
		html <- grep("DOCTYPE", raw[1, ])
	}
	colnames(raw) <- c("queryid", "subjectids", "identity", "positives", 
		"alignmentlength", "mismatches", "gapopens", "q.start", 
		"q.end", "s.start", "s.end", "evalue", "bitscore")
	rawm <- as.matrix(raw)
	eachsubject <- strsplit(rawm[, "subjectids"], ";")
	subjectids <- unlist(eachsubject)
	n.subjects <- sapply(eachsubject, length)
	rawm <- apply(rawm, 2, rep, times = n.subjects)
	rawm[, "subjectids"] <- subjectids
	all.ids <- strsplit(subjectids, "\\|")
	gi.id <- sapply(all.ids, "[", 2)
	pdb.id <- paste(sapply(all.ids, "[", 4), "_", sapply(all.ids, 
		"[", 5), sep = "")
	mlog.evalue <- -log(as.numeric(rawm[, "evalue"]))
	mlog.evalue[is.infinite(mlog.evalue)] <- -log(1e-308)
	cat(paste("\n Reporting", length(pdb.id), "hits\n"))
	output <- list(bitscore = as.numeric(rawm[, "bitscore"]), 
		evalue = as.numeric(rawm[, "evalue"]), mlog.evalue = mlog.evalue, 
		gi.id = gi.id, pdb.id = pdb.id, hit.tbl = rawm, raw = raw, rid=rid)
	class(output) <- "blast"
	return(output)
}

# bl <- ncbi.blastp("DPFEEHED", database="nr", hitlist=100, wordsize=2, expect=50000, matrix="PAM30", params="&GAPCOSTS=9%201&ENTREZ_QUERY=Homo%20Sapiens%20%5BORGN%5D&COMPOSITION_BASED_STATISTICS=no")
#http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&RID=9W4543YZ012
# http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&RESULTS_FILE=on&FORMAT_TYPE=CSV&ALIGNMENTS=20000&RID=9W7M57BC013


#http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=Alignment&ALIGNMENT_VIEW=Tabular&FORMAT_TYPE=CSV&ALIGNMENTS=20000&RID=9W9NAS75013


p.cox <- function(data, time, status, method="breslow", value=c("logtest","sctest","waldtest","rsq","coefficient")[1], attr=c("pvalue", "rsq")[1]) {
	coef <- (value=="coefficient")
	pcox <- numeric(nrow(data))
	library(survival)
	for (i in 1:nrow(data)) {
		icox <- NULL
		try(icox <- coxph(Surv(time, status) ~ ., data.frame(t(data[i,,drop=FALSE])), method=method))
		if (!is.null(icox)) {
			if (coef)
				pcox[i] <- summary(icox)$coefficients[1,5]
			else
				pcox[i] <- summary(icox)[[value]][attr]
		} else 
			pcox[i] <- NA
	}
	pcox
}



plot.hist2d <- function(x, y, nbins=50, col = c("white",heat.colors(16)), ...) {
	require(gplots)		
	xh <- hist2d(x,y, nbins=nbins, col = col, ...)
	rug(x,side=1)
	rug(y,side=2)
	box()
	invisible(xh)
}


those <- function(pattern, x, ...) {
	g <- grep(pattern, x, ...)
	if (length(g) == 0) return (NULL)
	x[g]
}



oligosignature <- function(s, len=6, normalize=FALSE, debug=FALSE) {
	library(gtools)
	idx <- apply(permutations(4,len,c("A","C","G","T"), repeats.allowed=TRUE), 1, paste,collapse="")
	m <- matrix(0,nrow=length(idx),ncol=length(s))
	rownames(m) <- idx
	colnames(m) <- names(s)
	for (i in 1:length(s)) {
		s[i] <- gsub("[^ACGT]","",s[i])
		if (debug) cat(i,"/",length(s),"...")
		n <- nchar(s[i])
		xs <- sapply(1:(n-len),function(j) substr(s[i],j,j+len-1))
		xx <- table(xs)
		m[names(xx),i] <- as.numeric(xx)/ifelse(normalize,n,1)
		if (debug) cat("\n")
	}
	m
}

rgb2rgb <- function(xrgb) {
	apply(xrgb, 2, function(x) { 
		if (any(x > 1)) x <- x / 255
		rgb(x[1],x[2],x[3])
		})
}


ls.usage <- function(xls=ls(envir=.GlobalEnv)) {
	l <- length(xls)
	xd <- data.frame(name=xls, stringsAsFactors=FALSE, class=rep("",l), mode=rep("",l), length=rep(0,l), rows=rep(0,l), cols=rep(0,l))
	for (i in 1:length(xls)) {
		xv <- get(xls[i])
		xclass <- class(xv)
		xd$class[i] <- paste(xclass,collapse=".")
		xd$mode[i] <- mode(xv)
		xd$length[i] <- length(xv)
		if (any(xclass=="matrix")  || any(xclass=="data.frame")) {
			try(xd$rows[i] <- nrow(xv))
			try(xd$cols[i] <- ncol(xv))
		}	
	}
	xd
}





cluster.exemplars <- function(data, 
	threshold = 0.8,
	similarity = function(x) cor(x, method="spearman"),
	merge = function(m) apply(m, 1, mean),
	reorder=NA,
	verbose=FALSE) {
	#data <- t(poldata) 
	#threshold <- 0.8

	## columns of data will be clustered
	cluster <- rep(0, ncol(data)) ### 0 is no cluster
	exemplar <- rep(FALSE, ncol(data))
	names(cluster) <- names(exemplar) <- colnames(data)

	k <- 1

	set.exemplar <- function(tocluster, verbose=FALSE) {
		wclus <- which(cluster == tocluster)
		m <- merge(data[,wclus])
		cs <- similarity(cbind(data[,wclus],m))
		csi <- cs[,ncol(cs)][-ncol(cs)]
		csm <- max(csi)
		exemplar[wclus] <- FALSE
		if (sum(csi == csm) > 1) {
			### seleccionar al de mayor desviacion
			s <- apply(data[,wclus], 2, sd)
			w <- which.max(s)[1]
		} else {
			### seleccionar al de mayor similaridad
			w <- wclus[which.max(csi)[1]]
		}
		exemplar[w] <- TRUE
		if (verbose) cat("Exemplar",w,"\n")
		exemplar
	}


	cycle <- 0
	while (TRUE) {
		cycle <- cycle + 1
		if (verbose) cat("Starting cycle ",cycle,"\n")
		used <- which(cluster == 0 | exemplar)
		remdata <- data[, used, drop=FALSE]
		sm <- similarity(remdata)
		ii <- matrix(1:ncol(remdata),ncol=2,nrow=ncol(remdata))
		sm[ii] <- NA
		wmx <- apply(sm, 2, function(s) { which.max(s)[1] })
		smx <- sm[matrix(c(1:nrow(sm), wmx), ncol=2)]
		omx <- order(smx, decreasing=TRUE)
		w <- which(smx[omx] < threshold)[1]-1
		if (is.na(w)) w <- length(smx)
		if (verbose) {
			cat("Used=",length(used),"\n")
			cat("w=",w,"\n")
		}
		joined <- FALSE
		if (length(w) > 0 && w > 0) {
			if (verbose) {
				print(table(cluster))
				cat("Similarities ",smx[omx[1:w]],"\n")
				cat("Checking	 ",used[omx[1:w]],"\n")
				cat("With		 ",used[wmx[1:w]],"\n")
			}
			for (i in 1:w) {
				l <- used[omx[i]]
				r <- used[wmx[omx[i]]]
				if (verbose) {
					if (k < 100) print(table(cluster))
					cat("Comparing",l,r)
					cat(" (",omx[i],",",sep="")
					cat(wmx[omx[i]],") ",sep="")
					cat(sm[omx[i],wmx[omx[i]]],"\n")
				}
				if (cluster[l] == 0 && cluster[r] == 0) {
					if (verbose) cat("Joining ",l,r,"\n")
					# asignar cluster		
					cluster[c(l,r)] <- k
					k <- k + 1
					exemplar <- set.exemplar(cluster[r], verbose=verbose)
					joined <- TRUE
				} else if (cluster[l] == cluster[r]) {
					# ya forman un cluster, no hacer nada
				} else {
					if (cluster[l] == 0 || cluster[r] == 0) {
						# uno ya esta en un cluster y el otro no						
						# l debe ser el que tiene el 0 y r el otro, checar si hay que intercambiar...
						if (cluster[r] == 0) {
							.l <- l
							l <- r
							r <- .l
						}
						tojoin <- c(l, which(cluster == cluster[r]))						
					} else {
						tojoin <- c(which(cluster == cluster[r]), which(cluster == cluster[l]))
					}
					sm2 <- similarity(data[,tojoin])
					ii <- matrix(1:length(tojoin),ncol=2,nrow=length(tojoin))
					sm2[ii] <- NA
					if (all(sm2 >= threshold, na.rm=TRUE)) {
						if (verbose) cat("Joining ",tojoin,"\n")
						cluster[tojoin] <- max(cluster[tojoin])
						exemplar <- set.exemplar(cluster[tojoin[1]], verbose=verbose)
						joined <- TRUE
					} else {
						## no se pueden unir, ignorar
					}
				}
			}
		} else {
			break
		}
		if (!joined) break
	}

	w0 <- which(cluster == 0)
	if (length(w0) > 0) {
		cluster[w0] <- 1:length(w0) + k
		exemplar[w0] <- TRUE
	}
	if (!is.null(reorder) && !is.na(reorder)  && reorder != FALSE) {
		if (verbose) cat("Reordering...\n")
		w <- which(exemplar)
		sm <- similarity(data[,w,drop=FALSE])
		if (FALSE) {
			ii <- matrix(1:length(w),ncol=2,nrow=length(w))
			sm[ii] <- NA
			wmx <- which(sm == max(sm,na.rm=TRUE), arr.ind=TRUE)
			xw <- sort(unique(as.vector(wmx)))
			xmn <- apply(sm[,xw,drop=FALSE], 2, mean, na.rm=TRUE)
			o <- xw[which.min(xmn)[1]]
			while (length(o) < ncol(sm)) {
				xo <- o[length(o)]
				wmx <- which.max(sm[,xo])[1]
				sm[,xo] <- NA
				sm[xo,] <- NA
				o <- c(o, wmx)
			}
			newcluster <- cluster-cluster
			for (i in 1:length(w)) {
				newcluster[cluster == cluster[w[i]]] <- o[i]
			}
		} else {
			am <- as.dist(ceiling(max(sm))-sm)
			hc <- hclust(am, method=reorder)
			newcluster <- cluster-cluster
			for (i in 1:length(w)) {
				newcluster[cluster == cluster[w[i]]] <- hc$order[i]
			}
		}
		cluster <- newcluster
	}

	list(cluster=as.numeric(factor(cluster)), exemplar=exemplar)
}

panel.hist <- function(x, col="cyan", ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

panel.cor <- function(x, y, digits=2, prefix="", cor.min.cex=0.5, cor.max.cex=3, cor.cex, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- (cor(x, y, use="complete.obs")) #abs
    xs <- summary(lm(y ~ x))
    #ptxt <- format(c(xs$coefficients[2,4], 0.123456789), digits=digits)[1]
    p <- xs$coefficients[2,4]
    ptxt <- ifelse(p < 0.05, ifelse(p < 0.01, ifelse(p < 0.001, "***", "**"),"*"), "")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    if(missing(cor.cex) || !missing(cor.min.cex) || !missing(cor.max.cex)) cor.cex <- abs(r)*(cor.max.cex-cor.min.cex)+cor.min.cex #cex.cor <- cex.cor.scale/strwidth(txt)
    text(0.5, 0.5, txt, cex = cor.cex)
    text(0.5, par("usr")[3], paste0(ptxt), adj=c(0.5,0), cex=max(cor.min.cex, cor.cex/2))
}

panel.points <- function(x,y, col.diag="red", lty.diag=1, lwd.diag=1, cex.points=0.75, col.points=1, pch.points=1, ...) {
	points(x,y,cex=cex.points,col=col.points,pch=pch.points,...)
	abline(0,1,col=8,lwd=1,lty=3)
	abline(lm(y~x),col=col.diag,lwd=lwd.diag,lty=lty.diag)
}


pairs.3.panels <- function(data, upper.panel=panel.cor, diag.panel=panel.hist, lower.panel=panel.points, gap=0, ...) {
	pairs(data, upper.panel=upper.panel, diag.panel=diag.panel, lower.panel=lower.panel, gap=gap, ...)
}


panel.cor.robust.points <- function(...) {
	panel.cor.robust(..., draw.points=TRUE)
}

panel.cor.robust.text <- function(...) {
	panel.cor.robust(..., draw.points=FALSE, draw.cor=TRUE)
}

pairs.3.panels.2 <- function(data, upper.panel=panel.cor.robust.points, 
	diag.panel=panel.hist, 
	lower.panel=panel.cor.robust.text, pch=20, gap=0, ...) {
	pairs(data, upper.panel=upper.panel, diag.panel=diag.panel, lower.panel=lower.panel, pch=pch, gap=gap, ...)
}

order.by.ranks <- function(data, rank.func=mean, return.order=TRUE, cluster=1, ...) {
	# Order rows of data analyzing the average ranks of the columns, cluster is used to reverse order for odd clusters
	wr <- apply(data,2,rank)
	xcor <- cor(wr)
	#wm <- which.max(apply(abs(xcor),1,sum))[1]
	#wneg <- which(xcor[,wm] < 0)
	pc <- apply(xcor > 0,2,sum)
	wneg <- which(pc < nrow(xcor)/2)
	for (w in wneg) {
		wr[,w] <- nrow(data)+1-wr[,w]
	}
	xval <- apply(wr,1,rank.func)
	if (return.order) {
		return (order(xval, decreasing=(cluster %% 2 == 0)))
	} else {
		return (xval)
	}
}

order.by.correlations <- function(data, return.order=TRUE, cluster=1, cor.method="pearson", ...) {
	# Order rows of data depending on the correlation with column order
	xord <- 1:ncol(data)
	if (is.function(cor.method)) {
		wr <- apply(data,1,function(x) cor.method)
	} else {
		wr <- apply(data,1,function(x) cor(x,xord,method=cor.method))
	}
	if (return.order) {
		return (order(wr, decreasing=(cluster %% 2 == 0)))
	} else {
		return (wr)
	}
}

order.by.means <- function(data, return.order=TRUE, cluster=1, ...) {
	order(apply(data, 1, mean), decreasing=(cluster %% 2 == 0))
}

order.by.sums <- function(data, return.order=TRUE, cluster=1, ...) {
	order(apply(data, 1, sum), decreasing=(cluster %% 2 == 0))
}

#Order ROWS by tsne : rows for compatibility with order.by.ranks ...
order.by.tsne <- function(data, return.order=TRUE, cluster=1, ...) {
	library(tsne)
	z <- tsne(t(x), k=1, perplexity=max(min(5,nrow(x)),round(nrow(x)/5+1)))[,1]
	if (return.order) z <- order(z, decreasing=(cluster %% 2 == 0))
	z
}


as.HTML <- function(xdf, row.names=!is.null(rownames(xdf)), col.names=!is.null(colnames(xdf)), prefix="", posfix="") {
	#s[1] <- "<TABLE>"
	wordize <- function(s) gsub("([0-9])([A-Za-z])","\\1 \\2",gsub("([a-z])([A-Z0-9])","\\1 \\2",gsub("[_.]"," ",s)))
	m <- apply(xdf, 2, format)
	if (!is.matrix(m) && !is.data.frame(m)) {
		save(m, file="m.RData")
		m <- matrix(m, nrow=nrow(xdf), ncol=ncol(xdf))
	}
	colnames(m) <- if (!is.null(colnames(xdf))) colnames(xdf) else ""
	rownames(m) <- if (!is.null(rownames(xdf))) rownames(xdf) else as.character(1:nrow(m))
	if (row.names) {
		m <- cbind(rownames(m),m)
	}
	s <- character(nrow(m)+3)
	s[1] <- "<TABLE>"
	s[2] <- paste("<TR><TH>", paste(wordize(colnames(m)),collapse="</TH><TH>"),"</TD></TR>",sep="")
	s[-c(1,2,length(s))] <- apply(m,1,function(m) paste("<TR><TD>", prefix, paste(m,collapse=paste(posfix,"</TD><TD>",prefix,sep="")),posfix,"</TD></TR>",sep=""))
	#s[-c(1,2,length(s))] <- apply(m,1,function(m) paste("<TR><TD>", paste(m,collapse="</TD><TD>"),"</TD></TR>",sep=""))
	s[length(s)] <- "</TABLE>"
	return (if (col.names) s else s[-2])
}


plot.colors <- function(...) {
	thecol <- colors()[-grep("gray|grey",colors())]
	xcol <- length(thecol)
	xr <- trunc(sqrt(xcol)*1)
	xc <- ceiling(xcol/xr)
	plot(0,0,xlim=c(1,xc+1), ylim=c(1,xr+1),type="n")
	text(rep(1:xr,each=xc), rep(1:xc,times=xr), thecol, col=thecol, ...)
}


#brewer.pal(0,return.names=TRUE)
# [1] "Accent"   "Blues"    "BrBG"     "BuGn"     "BuPu"     "Dark2"    "GnBu"     
# 	  "Greens"   "Greys"    "Oranges"  "OrRd"     "Paired"  
#[13] "Pastel1"  "Pastel2"  "PiYG"     "PRGn"     "PuBu"     "PuBuGn"   "PuOr"     
#     "PuRd"     "Purples"  "RdBu"     "RdGy"     "RdPu"    
#[25] "Reds"     "RdYlBu"   "RdYlGn"   "Set1"     "Set2"     "Set3"     "Spectral" 
#     "YlGn"     "YlGnBu"   "YlOrBr"   "YlOrRd"  

# Accent : "palidos" de verde, morado, naranja, amarillo, azul, violeta, cafe, gris
# Blues  : "palidos" de azules
# BrBG   : palidos de cafe-->blanco-->aqua
# BuGn   : palidos de blanco --> verde obscuro
# BuGn   : palidos de blanco --> verde obscuro
# BuPu   : palidos de celeste --> purple obscuro
# Dark2  : Similar a Accent pero más obscuros
# GnBu   : Palidos de verde a azul
# Greens : Verdes Palidos
# Oranges: Palidos Narajas desde claro a osbcuro café
# OrRd   : Pálidos narajas como oranges terminando en rojos
# Greys  : Grises
# Paired : Claros y obscuros palidos alternados de Azul, Verde, Rojo, Naranja, Morado
# Pastel1: Claros palidos de rojo, azul, verde, morado, naraja, amarillo, café, violeta, gris
# Pastel2: Claros palidos de verde, naraja, azul, violeta, verde, amarillo, cafe, gris
# PRGn   : Violetas Obscuros Pálidos --> Blancos --> Verdes Pálidos obscuros
# PuBu   : purple to blue
# PuBuGn : blancos pálidos, celestes, verdes
# PuOr   : Oranges obscuros, blancos, púrpuras obscuros
# RdYlBu : "CHIDO" palidos de Rojo-->Amarillo-->Azul
# RdYlGn : "CHIDO" palidos de Rojo-->Amarillo-->Verde
# Set1   : Brillantes modificados de rojo, azul, verde, morado, naraja, amarillo, cafe, rosa, gris
# Set2   : Muy palidos claros de verde naraja, azul, rosa, verde, amarillo, arena, gris
# Set3   : Muy palidos claros de verde naraja, azul, rosa, verde, amarillo, arena, gris
# Set3   : Medio palidos claros de verde amarillo, morado, salmon, azul, naraja, verde, rosa
# Spectral: "CHIDO" Rojo-->Amarillo-->Azul parecido a RdYlBu
brewer.pal <- function (n, name, return.names=FALSE) 
{   

    brewer.colorlist <- list(
    	Accent = switch(n - 2, rgb(c(127, 190, 253), 
        c(201, 174, 192), c(127, 212, 134), maxColorValue = 255), 
        rgb(c(127, 190, 253, 255), c(201, 174, 192, 255), c(127, 
            212, 134, 153), maxColorValue = 255), rgb(c(127, 
            190, 253, 255, 56), c(201, 174, 192, 255, 108), c(127, 
            212, 134, 153, 176), maxColorValue = 255), rgb(c(127, 
            190, 253, 255, 56, 240), c(201, 174, 192, 255, 108, 
            2), c(127, 212, 134, 153, 176, 127), maxColorValue = 255), 
        rgb(c(127, 190, 253, 255, 56, 240, 191), c(201, 174, 
            192, 255, 108, 2, 91), c(127, 212, 134, 153, 176, 
            127, 23), maxColorValue = 255), rgb(c(127, 190, 253, 
            255, 56, 240, 191, 102), c(201, 174, 192, 255, 108, 
            2, 91, 102), c(127, 212, 134, 153, 176, 127, 23, 
            102), maxColorValue = 255)), Blues = switch(n - 2, 
        rgb(c(222, 158, 49), c(235, 202, 130), c(247, 225, 189), 
            maxColorValue = 255), rgb(c(239, 189, 107, 33), c(243, 
            215, 174, 113), c(255, 231, 214, 181), maxColorValue = 255), 
        rgb(c(239, 189, 107, 49, 8), c(243, 215, 174, 130, 81), 
            c(255, 231, 214, 189, 156), maxColorValue = 255), 
        rgb(c(239, 198, 158, 107, 49, 8), c(243, 219, 202, 174, 
            130, 81), c(255, 239, 225, 214, 189, 156), maxColorValue = 255), 
        rgb(c(239, 198, 158, 107, 66, 33, 8), c(243, 219, 202, 
            174, 146, 113, 69), c(255, 239, 225, 214, 198, 181, 
            148), maxColorValue = 255), rgb(c(247, 222, 198, 
            158, 107, 66, 33, 8), c(251, 235, 219, 202, 174, 
            146, 113, 69), c(255, 247, 239, 225, 214, 198, 181, 
            148), maxColorValue = 255), rgb(c(247, 222, 198, 
            158, 107, 66, 33, 8, 8), c(251, 235, 219, 202, 174, 
            146, 113, 81, 48), c(255, 247, 239, 225, 214, 198, 
            181, 156, 107), maxColorValue = 255)), BrBG = switch(n - 
        2, rgb(c(216, 245, 90), c(179, 245, 180), c(101, 245, 
        172), maxColorValue = 255), rgb(c(166, 223, 128, 1), 
        c(97, 194, 205, 133), c(26, 125, 193, 113), maxColorValue = 255), 
        rgb(c(166, 223, 245, 128, 1), c(97, 194, 245, 205, 133), 
            c(26, 125, 245, 193, 113), maxColorValue = 255), 
        rgb(c(140, 216, 246, 199, 90, 1), c(81, 179, 232, 234, 
            180, 102), c(10, 101, 195, 229, 172, 94), maxColorValue = 255), 
        rgb(c(140, 216, 246, 245, 199, 90, 1), c(81, 179, 232, 
            245, 234, 180, 102), c(10, 101, 195, 245, 229, 172, 
            94), maxColorValue = 255), rgb(c(140, 191, 223, 246, 
            199, 128, 53, 1), c(81, 129, 194, 232, 234, 205, 
            151, 102), c(10, 45, 125, 195, 229, 193, 143, 94), 
            maxColorValue = 255), rgb(c(140, 191, 223, 246, 245, 
            199, 128, 53, 1), c(81, 129, 194, 232, 245, 234, 
            205, 151, 102), c(10, 45, 125, 195, 245, 229, 193, 
            143, 94), maxColorValue = 255), rgb(c(84, 140, 191, 
            223, 246, 199, 128, 53, 1, 0), c(48, 81, 129, 194, 
            232, 234, 205, 151, 102, 60), c(5, 10, 45, 125, 195, 
            229, 193, 143, 94, 48), maxColorValue = 255), rgb(c(84, 
            140, 191, 223, 246, 245, 199, 128, 53, 1, 0), c(48, 
            81, 129, 194, 232, 245, 234, 205, 151, 102, 60), 
            c(5, 10, 45, 125, 195, 245, 229, 193, 143, 94, 48), 
            maxColorValue = 255)), BuGn = switch(n - 2, rgb(c(229, 
        153, 44), c(245, 216, 162), c(249, 201, 95), maxColorValue = 255), 
        rgb(c(237, 178, 102, 35), c(248, 226, 194, 139), c(251, 
            226, 164, 69), maxColorValue = 255), rgb(c(237, 178, 
            102, 44, 0), c(248, 226, 194, 162, 109), c(251, 226, 
            164, 95, 44), maxColorValue = 255), rgb(c(237, 204, 
            153, 102, 44, 0), c(248, 236, 216, 194, 162, 109), 
            c(251, 230, 201, 164, 95, 44), maxColorValue = 255), 
        rgb(c(237, 204, 153, 102, 65, 35, 0), c(248, 236, 216, 
            194, 174, 139, 88), c(251, 230, 201, 164, 118, 69, 
            36), maxColorValue = 255), rgb(c(247, 229, 204, 153, 
            102, 65, 35, 0), c(252, 245, 236, 216, 194, 174, 
            139, 88), c(253, 249, 230, 201, 164, 118, 69, 36), 
            maxColorValue = 255), rgb(c(247, 229, 204, 153, 102, 
            65, 35, 0, 0), c(252, 245, 236, 216, 194, 174, 139, 
            109, 68), c(253, 249, 230, 201, 164, 118, 69, 44, 
            27), maxColorValue = 255)), BuPu = switch(n - 2, 
        rgb(c(224, 158, 136), c(236, 188, 86), c(244, 218, 167), 
            maxColorValue = 255), rgb(c(237, 179, 140, 136), 
            c(248, 205, 150, 65), c(251, 227, 198, 157), maxColorValue = 255), 
        rgb(c(237, 179, 140, 136, 129), c(248, 205, 150, 86, 
            15), c(251, 227, 198, 167, 124), maxColorValue = 255), 
        rgb(c(237, 191, 158, 140, 136, 129), c(248, 211, 188, 
            150, 86, 15), c(251, 230, 218, 198, 167, 124), maxColorValue = 255), 
        rgb(c(237, 191, 158, 140, 140, 136, 110), c(248, 211, 
            188, 150, 107, 65, 1), c(251, 230, 218, 198, 177, 
            157, 107), maxColorValue = 255), rgb(c(247, 224, 
            191, 158, 140, 140, 136, 110), c(252, 236, 211, 188, 
            150, 107, 65, 1), c(253, 244, 230, 218, 198, 177, 
            157, 107), maxColorValue = 255), rgb(c(247, 224, 
            191, 158, 140, 140, 136, 129, 77), c(252, 236, 211, 
            188, 150, 107, 65, 15, 0), c(253, 244, 230, 218, 
            198, 177, 157, 124, 75), maxColorValue = 255)), Dark2 = switch(n - 
        2, rgb(c(27, 217, 117), c(158, 95, 112), c(119, 2, 179), 
        maxColorValue = 255), rgb(c(27, 217, 117, 231), c(158, 
        95, 112, 41), c(119, 2, 179, 138), maxColorValue = 255), 
        rgb(c(27, 217, 117, 231, 102), c(158, 95, 112, 41, 166), 
            c(119, 2, 179, 138, 30), maxColorValue = 255), rgb(c(27, 
            217, 117, 231, 102, 230), c(158, 95, 112, 41, 166, 
            171), c(119, 2, 179, 138, 30, 2), maxColorValue = 255), 
        rgb(c(27, 217, 117, 231, 102, 230, 166), c(158, 95, 112, 
            41, 166, 171, 118), c(119, 2, 179, 138, 30, 2, 29), 
            maxColorValue = 255), rgb(c(27, 217, 117, 231, 102, 
            230, 166, 102), c(158, 95, 112, 41, 166, 171, 118, 
            102), c(119, 2, 179, 138, 30, 2, 29, 102), maxColorValue = 255)), 
        GnBu = switch(n - 2, rgb(c(224, 168, 67), c(243, 221, 
            162), c(219, 181, 202), maxColorValue = 255), rgb(c(240, 
            186, 123, 43), c(249, 228, 204, 140), c(232, 188, 
            196, 190), maxColorValue = 255), rgb(c(240, 186, 
            123, 67, 8), c(249, 228, 204, 162, 104), c(232, 188, 
            196, 202, 172), maxColorValue = 255), rgb(c(240, 
            204, 168, 123, 67, 8), c(249, 235, 221, 204, 162, 
            104), c(232, 197, 181, 196, 202, 172), maxColorValue = 255), 
            rgb(c(240, 204, 168, 123, 78, 43, 8), c(249, 235, 
                221, 204, 179, 140, 88), c(232, 197, 181, 196, 
                211, 190, 158), maxColorValue = 255), rgb(c(247, 
                224, 204, 168, 123, 78, 43, 8), c(252, 243, 235, 
                221, 204, 179, 140, 88), c(240, 219, 197, 181, 
                196, 211, 190, 158), maxColorValue = 255), rgb(c(247, 
                224, 204, 168, 123, 78, 43, 8, 8), c(252, 243, 
                235, 221, 204, 179, 140, 104, 64), c(240, 219, 
                197, 181, 196, 211, 190, 172, 129), maxColorValue = 255)), 
        Greens = switch(n - 2, rgb(c(229, 161, 49), c(245, 217, 
            163), c(224, 155, 84), maxColorValue = 255), rgb(c(237, 
            186, 116, 35), c(248, 228, 196, 139), c(233, 179, 
            118, 69), maxColorValue = 255), rgb(c(237, 186, 116, 
            49, 0), c(248, 228, 196, 163, 109), c(233, 179, 118, 
            84, 44), maxColorValue = 255), rgb(c(237, 199, 161, 
            116, 49, 0), c(248, 233, 217, 196, 163, 109), c(233, 
            192, 155, 118, 84, 44), maxColorValue = 255), rgb(c(237, 
            199, 161, 116, 65, 35, 0), c(248, 233, 217, 196, 
            171, 139, 90), c(233, 192, 155, 118, 93, 69, 50), 
            maxColorValue = 255), rgb(c(247, 229, 199, 161, 116, 
            65, 35, 0), c(252, 245, 233, 217, 196, 171, 139, 
            90), c(245, 224, 192, 155, 118, 93, 69, 50), maxColorValue = 255), 
            rgb(c(247, 229, 199, 161, 116, 65, 35, 0, 0), c(252, 
                245, 233, 217, 196, 171, 139, 109, 68), c(245, 
                224, 192, 155, 118, 93, 69, 44, 27), maxColorValue = 255)), 
        Greys = switch(n - 2, rgb(c(240, 189, 99), c(240, 189, 
            99), c(240, 189, 99), maxColorValue = 255), rgb(c(247, 
            204, 150, 82), c(247, 204, 150, 82), c(247, 204, 
            150, 82), maxColorValue = 255), rgb(c(247, 204, 150, 
            99, 37), c(247, 204, 150, 99, 37), c(247, 204, 150, 
            99, 37), maxColorValue = 255), rgb(c(247, 217, 189, 
            150, 99, 37), c(247, 217, 189, 150, 99, 37), c(247, 
            217, 189, 150, 99, 37), maxColorValue = 255), rgb(c(247, 
            217, 189, 150, 115, 82, 37), c(247, 217, 189, 150, 
            115, 82, 37), c(247, 217, 189, 150, 115, 82, 37), 
            maxColorValue = 255), rgb(c(255, 240, 217, 189, 150, 
            115, 82, 37), c(255, 240, 217, 189, 150, 115, 82, 
            37), c(255, 240, 217, 189, 150, 115, 82, 37), maxColorValue = 255), 
            rgb(c(255, 240, 217, 189, 150, 115, 82, 37, 0), c(255, 
                240, 217, 189, 150, 115, 82, 37, 0), c(255, 240, 
                217, 189, 150, 115, 82, 37, 0), maxColorValue = 255)), 
        Oranges = switch(n - 2, rgb(c(254, 253, 230), c(230, 
            174, 85), c(206, 107, 13), maxColorValue = 255), 
            rgb(c(254, 253, 253, 217), c(237, 190, 141, 71), 
                c(222, 133, 60, 1), maxColorValue = 255), rgb(c(254, 
                253, 253, 230, 166), c(237, 190, 141, 85, 54), 
                c(222, 133, 60, 13, 3), maxColorValue = 255), 
            rgb(c(254, 253, 253, 253, 230, 166), c(237, 208, 
                174, 141, 85, 54), c(222, 162, 107, 60, 13, 3), 
                maxColorValue = 255), rgb(c(254, 253, 253, 253, 
                241, 217, 140), c(237, 208, 174, 141, 105, 72, 
                45), c(222, 162, 107, 60, 19, 1, 4), maxColorValue = 255), 
            rgb(c(255, 254, 253, 253, 253, 241, 217, 140), c(245, 
                230, 208, 174, 141, 105, 72, 45), c(235, 206, 
                162, 107, 60, 19, 1, 4), maxColorValue = 255), 
            rgb(c(255, 254, 253, 253, 253, 241, 217, 166, 127), 
                c(245, 230, 208, 174, 141, 105, 72, 54, 39), 
                c(235, 206, 162, 107, 60, 19, 1, 3, 4), maxColorValue = 255)), 
        OrRd = switch(n - 2, rgb(c(254, 253, 227), c(232, 187, 
            74), c(200, 132, 51), maxColorValue = 255), rgb(c(254, 
            253, 252, 215), c(240, 204, 141, 48), c(217, 138, 
            89, 31), maxColorValue = 255), rgb(c(254, 253, 252, 
            227, 179), c(240, 204, 141, 74, 0), c(217, 138, 89, 
            51, 0), maxColorValue = 255), rgb(c(254, 253, 253, 
            252, 227, 179), c(240, 212, 187, 141, 74, 0), c(217, 
            158, 132, 89, 51, 0), maxColorValue = 255), rgb(c(254, 
            253, 253, 252, 239, 215, 153), c(240, 212, 187, 141, 
            101, 48, 0), c(217, 158, 132, 89, 72, 31, 0), maxColorValue = 255), 
            rgb(c(255, 254, 253, 253, 252, 239, 215, 153), c(247, 
                232, 212, 187, 141, 101, 48, 0), c(236, 200, 
                158, 132, 89, 72, 31, 0), maxColorValue = 255), 
            rgb(c(255, 254, 253, 253, 252, 239, 215, 179, 127), 
                c(247, 232, 212, 187, 141, 101, 48, 0, 0), c(236, 
                  200, 158, 132, 89, 72, 31, 0, 0), maxColorValue = 255)), 
        Paired = switch(n - 2, rgb(c(166, 31, 178), c(206, 120, 
            223), c(227, 180, 138), maxColorValue = 255), rgb(c(166, 
            31, 178, 51), c(206, 120, 223, 160), c(227, 180, 
            138, 44), maxColorValue = 255), rgb(c(166, 31, 178, 
            51, 251), c(206, 120, 223, 160, 154), c(227, 180, 
            138, 44, 153), maxColorValue = 255), rgb(c(166, 31, 
            178, 51, 251, 227), c(206, 120, 223, 160, 154, 26), 
            c(227, 180, 138, 44, 153, 28), maxColorValue = 255), 
            rgb(c(166, 31, 178, 51, 251, 227, 253), c(206, 120, 
                223, 160, 154, 26, 191), c(227, 180, 138, 44, 
                153, 28, 111), maxColorValue = 255), rgb(c(166, 
                31, 178, 51, 251, 227, 253, 255), c(206, 120, 
                223, 160, 154, 26, 191, 127), c(227, 180, 138, 
                44, 153, 28, 111, 0), maxColorValue = 255), rgb(c(166, 
                31, 178, 51, 251, 227, 253, 255, 202), c(206, 
                120, 223, 160, 154, 26, 191, 127, 178), c(227, 
                180, 138, 44, 153, 28, 111, 0, 214), maxColorValue = 255), 
            rgb(c(166, 31, 178, 51, 251, 227, 253, 255, 202, 
                106), c(206, 120, 223, 160, 154, 26, 191, 127, 
                178, 61), c(227, 180, 138, 44, 153, 28, 111, 
                0, 214, 154), maxColorValue = 255), rgb(c(166, 
                31, 178, 51, 251, 227, 253, 255, 202, 106, 255), 
                c(206, 120, 223, 160, 154, 26, 191, 127, 178, 
                  61, 255), c(227, 180, 138, 44, 153, 28, 111, 
                  0, 214, 154, 153), maxColorValue = 255), rgb(c(166, 
                31, 178, 51, 251, 227, 253, 255, 202, 106, 255, 
                177), c(206, 120, 223, 160, 154, 26, 191, 127, 
                178, 61, 255, 89), c(227, 180, 138, 44, 153, 
                28, 111, 0, 214, 154, 153, 40), maxColorValue = 255)), 
        Pastel1 = switch(n - 2, rgb(c(251, 179, 204), c(180, 
            205, 235), c(174, 227, 197), maxColorValue = 255), 
            rgb(c(251, 179, 204, 222), c(180, 205, 235, 203), 
                c(174, 227, 197, 228), maxColorValue = 255), 
            rgb(c(251, 179, 204, 222, 254), c(180, 205, 235, 
                203, 217), c(174, 227, 197, 228, 166), maxColorValue = 255), 
            rgb(c(251, 179, 204, 222, 254, 255), c(180, 205, 
                235, 203, 217, 255), c(174, 227, 197, 228, 166, 
                204), maxColorValue = 255), rgb(c(251, 179, 204, 
                222, 254, 255, 229), c(180, 205, 235, 203, 217, 
                255, 216), c(174, 227, 197, 228, 166, 204, 189), 
                maxColorValue = 255), rgb(c(251, 179, 204, 222, 
                254, 255, 229, 253), c(180, 205, 235, 203, 217, 
                255, 216, 218), c(174, 227, 197, 228, 166, 204, 
                189, 236), maxColorValue = 255), rgb(c(251, 179, 
                204, 222, 254, 255, 229, 253, 242), c(180, 205, 
                235, 203, 217, 255, 216, 218, 242), c(174, 227, 
                197, 228, 166, 204, 189, 236, 242), maxColorValue = 255)), 
        Pastel2 = switch(n - 2, rgb(c(179, 253, 203), c(226, 
            205, 213), c(205, 172, 232), maxColorValue = 255), 
            rgb(c(179, 253, 203, 244), c(226, 205, 213, 202), 
                c(205, 172, 232, 228), maxColorValue = 255), 
            rgb(c(179, 253, 203, 244, 230), c(226, 205, 213, 
                202, 245), c(205, 172, 232, 228, 201), maxColorValue = 255), 
            rgb(c(179, 253, 203, 244, 230, 255), c(226, 205, 
                213, 202, 245, 242), c(205, 172, 232, 228, 201, 
                174), maxColorValue = 255), rgb(c(179, 253, 203, 
                244, 230, 255, 241), c(226, 205, 213, 202, 245, 
                242, 226), c(205, 172, 232, 228, 201, 174, 204), 
                maxColorValue = 255), rgb(c(179, 253, 203, 244, 
                230, 255, 241, 204), c(226, 205, 213, 202, 245, 
                242, 226, 204), c(205, 172, 232, 228, 201, 174, 
                204, 204), maxColorValue = 255)), PiYG = switch(n - 
            2, rgb(c(233, 247, 161), c(163, 247, 215), c(201, 
            247, 106), maxColorValue = 255), rgb(c(208, 241, 
            184, 77), c(28, 182, 225, 172), c(139, 218, 134, 
            38), maxColorValue = 255), rgb(c(208, 241, 247, 184, 
            77), c(28, 182, 247, 225, 172), c(139, 218, 247, 
            134, 38), maxColorValue = 255), rgb(c(197, 233, 253, 
            230, 161, 77), c(27, 163, 224, 245, 215, 146), c(125, 
            201, 239, 208, 106, 33), maxColorValue = 255), rgb(c(197, 
            233, 253, 247, 230, 161, 77), c(27, 163, 224, 247, 
            245, 215, 146), c(125, 201, 239, 247, 208, 106, 33), 
            maxColorValue = 255), rgb(c(197, 222, 241, 253, 230, 
            184, 127, 77), c(27, 119, 182, 224, 245, 225, 188, 
            146), c(125, 174, 218, 239, 208, 134, 65, 33), maxColorValue = 255), 
            rgb(c(197, 222, 241, 253, 247, 230, 184, 127, 77), 
                c(27, 119, 182, 224, 247, 245, 225, 188, 146), 
                c(125, 174, 218, 239, 247, 208, 134, 65, 33), 
                maxColorValue = 255), rgb(c(142, 197, 222, 241, 
                253, 230, 184, 127, 77, 39), c(1, 27, 119, 182, 
                224, 245, 225, 188, 146, 100), c(82, 125, 174, 
                218, 239, 208, 134, 65, 33, 25), maxColorValue = 255), 
            rgb(c(142, 197, 222, 241, 253, 247, 230, 184, 127, 
                77, 39), c(1, 27, 119, 182, 224, 247, 245, 225, 
                188, 146, 100), c(82, 125, 174, 218, 239, 247, 
                208, 134, 65, 33, 25), maxColorValue = 255)), 
        PRGn = switch(n - 2, rgb(c(175, 247, 127), c(141, 247, 
            191), c(195, 247, 123), maxColorValue = 255), rgb(c(123, 
            194, 166, 0), c(50, 165, 219, 136), c(148, 207, 160, 
            55), maxColorValue = 255), rgb(c(123, 194, 247, 166, 
            0), c(50, 165, 247, 219, 136), c(148, 207, 247, 160, 
            55), maxColorValue = 255), rgb(c(118, 175, 231, 217, 
            127, 27), c(42, 141, 212, 240, 191, 120), c(131, 
            195, 232, 211, 123, 55), maxColorValue = 255), rgb(c(118, 
            175, 231, 247, 217, 127, 27), c(42, 141, 212, 247, 
            240, 191, 120), c(131, 195, 232, 247, 211, 123, 55), 
            maxColorValue = 255), rgb(c(118, 153, 194, 231, 217, 
            166, 90, 27), c(42, 112, 165, 212, 240, 219, 174, 
            120), c(131, 171, 207, 232, 211, 160, 97, 55), maxColorValue = 255), 
            rgb(c(118, 153, 194, 231, 247, 217, 166, 90, 27), 
                c(42, 112, 165, 212, 247, 240, 219, 174, 120), 
                c(131, 171, 207, 232, 247, 211, 160, 97, 55), 
                maxColorValue = 255), rgb(c(64, 118, 153, 194, 
                231, 217, 166, 90, 27, 0), c(0, 42, 112, 165, 
                212, 240, 219, 174, 120, 68), c(75, 131, 171, 
                207, 232, 211, 160, 97, 55, 27), maxColorValue = 255), 
            rgb(c(64, 118, 153, 194, 231, 247, 217, 166, 90, 
                27, 0), c(0, 42, 112, 165, 212, 247, 240, 219, 
                174, 120, 68), c(75, 131, 171, 207, 232, 247, 
                211, 160, 97, 55, 27), maxColorValue = 255)), 
        PuBu = switch(n - 2, rgb(c(236, 166, 43), c(231, 189, 
            140), c(242, 219, 190), maxColorValue = 255), rgb(c(241, 
            189, 116, 5), c(238, 201, 169, 112), c(246, 225, 
            207, 176), maxColorValue = 255), rgb(c(241, 189, 
            116, 43, 4), c(238, 201, 169, 140, 90), c(246, 225, 
            207, 190, 141), maxColorValue = 255), rgb(c(241, 
            208, 166, 116, 43, 4), c(238, 209, 189, 169, 140, 
            90), c(246, 230, 219, 207, 190, 141), maxColorValue = 255), 
            rgb(c(241, 208, 166, 116, 54, 5, 3), c(238, 209, 
                189, 169, 144, 112, 78), c(246, 230, 219, 207, 
                192, 176, 123), maxColorValue = 255), rgb(c(255, 
                236, 208, 166, 116, 54, 5, 3), c(247, 231, 209, 
                189, 169, 144, 112, 78), c(251, 242, 230, 219, 
                207, 192, 176, 123), maxColorValue = 255), rgb(c(255, 
                236, 208, 166, 116, 54, 5, 4, 2), c(247, 231, 
                209, 189, 169, 144, 112, 90, 56), c(251, 242, 
                230, 219, 207, 192, 176, 141, 88), maxColorValue = 255)), 
        PuBuGn = switch(n - 2, rgb(c(236, 166, 28), c(226, 189, 
            144), c(240, 219, 153), maxColorValue = 255), rgb(c(246, 
            189, 103, 2), c(239, 201, 169, 129), c(247, 225, 
            207, 138), maxColorValue = 255), rgb(c(246, 189, 
            103, 28, 1), c(239, 201, 169, 144, 108), c(247, 225, 
            207, 153, 89), maxColorValue = 255), rgb(c(246, 208, 
            166, 103, 28, 1), c(239, 209, 189, 169, 144, 108), 
            c(247, 230, 219, 207, 153, 89), maxColorValue = 255), 
            rgb(c(246, 208, 166, 103, 54, 2, 1), c(239, 209, 
                189, 169, 144, 129, 100), c(247, 230, 219, 207, 
                192, 138, 80), maxColorValue = 255), rgb(c(255, 
                236, 208, 166, 103, 54, 2, 1), c(247, 226, 209, 
                189, 169, 144, 129, 100), c(251, 240, 230, 219, 
                207, 192, 138, 80), maxColorValue = 255), rgb(c(255, 
                236, 208, 166, 103, 54, 2, 1, 1), c(247, 226, 
                209, 189, 169, 144, 129, 108, 70), c(251, 240, 
                230, 219, 207, 192, 138, 89, 54), maxColorValue = 255)), 
        PuOr = switch(n - 2, rgb(c(241, 247, 153), c(163, 247, 
            142), c(64, 247, 195), maxColorValue = 255), rgb(c(230, 
            253, 178, 94), c(97, 184, 171, 60), c(1, 99, 210, 
            153), maxColorValue = 255), rgb(c(230, 253, 247, 
            178, 94), c(97, 184, 247, 171, 60), c(1, 99, 247, 
            210, 153), maxColorValue = 255), rgb(c(179, 241, 
            254, 216, 153, 84), c(88, 163, 224, 218, 142, 39), 
            c(6, 64, 182, 235, 195, 136), maxColorValue = 255), 
            rgb(c(179, 241, 254, 247, 216, 153, 84), c(88, 163, 
                224, 247, 218, 142, 39), c(6, 64, 182, 247, 235, 
                195, 136), maxColorValue = 255), rgb(c(179, 224, 
                253, 254, 216, 178, 128, 84), c(88, 130, 184, 
                224, 218, 171, 115, 39), c(6, 20, 99, 182, 235, 
                210, 172, 136), maxColorValue = 255), rgb(c(179, 
                224, 253, 254, 247, 216, 178, 128, 84), c(88, 
                130, 184, 224, 247, 218, 171, 115, 39), c(6, 
                20, 99, 182, 247, 235, 210, 172, 136), maxColorValue = 255), 
            rgb(c(127, 179, 224, 253, 254, 216, 178, 128, 84, 
                45), c(59, 88, 130, 184, 224, 218, 171, 115, 
                39, 0), c(8, 6, 20, 99, 182, 235, 210, 172, 136, 
                75), maxColorValue = 255), rgb(c(127, 179, 224, 
                253, 254, 247, 216, 178, 128, 84, 45), c(59, 
                88, 130, 184, 224, 247, 218, 171, 115, 39, 0), 
                c(8, 6, 20, 99, 182, 247, 235, 210, 172, 136, 
                  75), maxColorValue = 255)), PuRd = switch(n - 
            2, rgb(c(231, 201, 221), c(225, 148, 28), c(239, 
            199, 119), maxColorValue = 255), rgb(c(241, 215, 
            223, 206), c(238, 181, 101, 18), c(246, 216, 176, 
            86), maxColorValue = 255), rgb(c(241, 215, 223, 221, 
            152), c(238, 181, 101, 28, 0), c(246, 216, 176, 119, 
            67), maxColorValue = 255), rgb(c(241, 212, 201, 223, 
            221, 152), c(238, 185, 148, 101, 28, 0), c(246, 218, 
            199, 176, 119, 67), maxColorValue = 255), rgb(c(241, 
            212, 201, 223, 231, 206, 145), c(238, 185, 148, 101, 
            41, 18, 0), c(246, 218, 199, 176, 138, 86, 63), maxColorValue = 255), 
            rgb(c(247, 231, 212, 201, 223, 231, 206, 145), c(244, 
                225, 185, 148, 101, 41, 18, 0), c(249, 239, 218, 
                199, 176, 138, 86, 63), maxColorValue = 255), 
            rgb(c(247, 231, 212, 201, 223, 231, 206, 152, 103), 
                c(244, 225, 185, 148, 101, 41, 18, 0, 0), c(249, 
                  239, 218, 199, 176, 138, 86, 67, 31), maxColorValue = 255)), 
        Purples = switch(n - 2, rgb(c(239, 188, 117), c(237, 
            189, 107), c(245, 220, 177), maxColorValue = 255), 
            rgb(c(242, 203, 158, 106), c(240, 201, 154, 81), 
                c(247, 226, 200, 163), maxColorValue = 255), 
            rgb(c(242, 203, 158, 117, 84), c(240, 201, 154, 107, 
                39), c(247, 226, 200, 177, 143), maxColorValue = 255), 
            rgb(c(242, 218, 188, 158, 117, 84), c(240, 218, 189, 
                154, 107, 39), c(247, 235, 220, 200, 177, 143), 
                maxColorValue = 255), rgb(c(242, 218, 188, 158, 
                128, 106, 74), c(240, 218, 189, 154, 125, 81, 
                20), c(247, 235, 220, 200, 186, 163, 134), maxColorValue = 255), 
            rgb(c(252, 239, 218, 188, 158, 128, 106, 74), c(251, 
                237, 218, 189, 154, 125, 81, 20), c(253, 245, 
                235, 220, 200, 186, 163, 134), maxColorValue = 255), 
            rgb(c(252, 239, 218, 188, 158, 128, 106, 84, 63), 
                c(251, 237, 218, 189, 154, 125, 81, 39, 0), c(253, 
                  245, 235, 220, 200, 186, 163, 143, 125), maxColorValue = 255)), 
        RdBu = switch(n - 2, rgb(c(239, 247, 103), c(138, 247, 
            169), c(98, 247, 207), maxColorValue = 255), rgb(c(202, 
            244, 146, 5), c(0, 165, 197, 113), c(32, 130, 222, 
            176), maxColorValue = 255), rgb(c(202, 244, 247, 
            146, 5), c(0, 165, 247, 197, 113), c(32, 130, 247, 
            222, 176), maxColorValue = 255), rgb(c(178, 239, 
            253, 209, 103, 33), c(24, 138, 219, 229, 169, 102), 
            c(43, 98, 199, 240, 207, 172), maxColorValue = 255), 
            rgb(c(178, 239, 253, 247, 209, 103, 33), c(24, 138, 
                219, 247, 229, 169, 102), c(43, 98, 199, 247, 
                240, 207, 172), maxColorValue = 255), rgb(c(178, 
                214, 244, 253, 209, 146, 67, 33), c(24, 96, 165, 
                219, 229, 197, 147, 102), c(43, 77, 130, 199, 
                240, 222, 195, 172), maxColorValue = 255), rgb(c(178, 
                214, 244, 253, 247, 209, 146, 67, 33), c(24, 
                96, 165, 219, 247, 229, 197, 147, 102), c(43, 
                77, 130, 199, 247, 240, 222, 195, 172), maxColorValue = 255), 
            rgb(c(103, 178, 214, 244, 253, 209, 146, 67, 33, 
                5), c(0, 24, 96, 165, 219, 229, 197, 147, 102, 
                48), c(31, 43, 77, 130, 199, 240, 222, 195, 172, 
                97), maxColorValue = 255), rgb(c(103, 178, 214, 
                244, 253, 247, 209, 146, 67, 33, 5), c(0, 24, 
                96, 165, 219, 247, 229, 197, 147, 102, 48), c(31, 
                43, 77, 130, 199, 247, 240, 222, 195, 172, 97), 
                maxColorValue = 255)), RdGy = switch(n - 2, rgb(c(239, 
            255, 153), c(138, 255, 153), c(98, 255, 153), maxColorValue = 255), 
            rgb(c(202, 244, 186, 64), c(0, 165, 186, 64), c(32, 
                130, 186, 64), maxColorValue = 255), rgb(c(202, 
                244, 255, 186, 64), c(0, 165, 255, 186, 64), 
                c(32, 130, 255, 186, 64), maxColorValue = 255), 
            rgb(c(178, 239, 253, 224, 153, 77), c(24, 138, 219, 
                224, 153, 77), c(43, 98, 199, 224, 153, 77), 
                maxColorValue = 255), rgb(c(178, 239, 253, 255, 
                224, 153, 77), c(24, 138, 219, 255, 224, 153, 
                77), c(43, 98, 199, 255, 224, 153, 77), maxColorValue = 255), 
            rgb(c(178, 214, 244, 253, 224, 186, 135, 77), c(24, 
                96, 165, 219, 224, 186, 135, 77), c(43, 77, 130, 
                199, 224, 186, 135, 77), maxColorValue = 255), 
            rgb(c(178, 214, 244, 253, 255, 224, 186, 135, 77), 
                c(24, 96, 165, 219, 255, 224, 186, 135, 77), 
                c(43, 77, 130, 199, 255, 224, 186, 135, 77), 
                maxColorValue = 255), rgb(c(103, 178, 214, 244, 
                253, 224, 186, 135, 77, 26), c(0, 24, 96, 165, 
                219, 224, 186, 135, 77, 26), c(31, 43, 77, 130, 
                199, 224, 186, 135, 77, 26), maxColorValue = 255), 
            rgb(c(103, 178, 214, 244, 253, 255, 224, 186, 135, 
                77, 26), c(0, 24, 96, 165, 219, 255, 224, 186, 
                135, 77, 26), c(31, 43, 77, 130, 199, 255, 224, 
                186, 135, 77, 26), maxColorValue = 255)), RdPu = switch(n - 
            2, rgb(c(253, 250, 197), c(224, 159, 27), c(221, 
            181, 138), maxColorValue = 255), rgb(c(254, 251, 
            247, 174), c(235, 180, 104, 1), c(226, 185, 161, 
            126), maxColorValue = 255), rgb(c(254, 251, 247, 
            197, 122), c(235, 180, 104, 27, 1), c(226, 185, 161, 
            138, 119), maxColorValue = 255), rgb(c(254, 252, 
            250, 247, 197, 122), c(235, 197, 159, 104, 27, 1), 
            c(226, 192, 181, 161, 138, 119), maxColorValue = 255), 
            rgb(c(254, 252, 250, 247, 221, 174, 122), c(235, 
                197, 159, 104, 52, 1, 1), c(226, 192, 181, 161, 
                151, 126, 119), maxColorValue = 255), rgb(c(255, 
                253, 252, 250, 247, 221, 174, 122), c(247, 224, 
                197, 159, 104, 52, 1, 1), c(243, 221, 192, 181, 
                161, 151, 126, 119), maxColorValue = 255), rgb(c(255, 
                253, 252, 250, 247, 221, 174, 122, 73), c(247, 
                224, 197, 159, 104, 52, 1, 1, 0), c(243, 221, 
                192, 181, 161, 151, 126, 119, 106), maxColorValue = 255)), 
        Reds = switch(n - 2, rgb(c(254, 252, 222), c(224, 146, 
            45), c(210, 114, 38), maxColorValue = 255), rgb(c(254, 
            252, 251, 203), c(229, 174, 106, 24), c(217, 145, 
            74, 29), maxColorValue = 255), rgb(c(254, 252, 251, 
            222, 165), c(229, 174, 106, 45, 15), c(217, 145, 
            74, 38, 21), maxColorValue = 255), rgb(c(254, 252, 
            252, 251, 222, 165), c(229, 187, 146, 106, 45, 15), 
            c(217, 161, 114, 74, 38, 21), maxColorValue = 255), 
            rgb(c(254, 252, 252, 251, 239, 203, 153), c(229, 
                187, 146, 106, 59, 24, 0), c(217, 161, 114, 74, 
                44, 29, 13), maxColorValue = 255), rgb(c(255, 
                254, 252, 252, 251, 239, 203, 153), c(245, 224, 
                187, 146, 106, 59, 24, 0), c(240, 210, 161, 114, 
                74, 44, 29, 13), maxColorValue = 255), rgb(c(255, 
                254, 252, 252, 251, 239, 203, 165, 103), c(245, 
                224, 187, 146, 106, 59, 24, 15, 0), c(240, 210, 
                161, 114, 74, 44, 29, 21, 13), maxColorValue = 255)), 
        RdYlBu = switch(n - 2, rgb(c(252, 255, 145), c(141, 255, 
            191), c(89, 191, 219), maxColorValue = 255), rgb(c(215, 
            253, 171, 44), c(25, 174, 217, 123), c(28, 97, 233, 
            182), maxColorValue = 255), rgb(c(215, 253, 255, 
            171, 44), c(25, 174, 255, 217, 123), c(28, 97, 191, 
            233, 182), maxColorValue = 255), rgb(c(215, 252, 
            254, 224, 145, 69), c(48, 141, 224, 243, 191, 117), 
            c(39, 89, 144, 248, 219, 180), maxColorValue = 255), 
            rgb(c(215, 252, 254, 255, 224, 145, 69), c(48, 141, 
                224, 255, 243, 191, 117), c(39, 89, 144, 191, 
                248, 219, 180), maxColorValue = 255), rgb(c(215, 
                244, 253, 254, 224, 171, 116, 69), c(48, 109, 
                174, 224, 243, 217, 173, 117), c(39, 67, 97, 
                144, 248, 233, 209, 180), maxColorValue = 255), 
            rgb(c(215, 244, 253, 254, 255, 224, 171, 116, 69), 
                c(48, 109, 174, 224, 255, 243, 217, 173, 117), 
                c(39, 67, 97, 144, 191, 248, 233, 209, 180), 
                maxColorValue = 255), rgb(c(165, 215, 244, 253, 
                254, 224, 171, 116, 69, 49), c(0, 48, 109, 174, 
                224, 243, 217, 173, 117, 54), c(38, 39, 67, 97, 
                144, 248, 233, 209, 180, 149), maxColorValue = 255), 
            rgb(c(165, 215, 244, 253, 254, 255, 224, 171, 116, 
                69, 49), c(0, 48, 109, 174, 224, 255, 243, 217, 
                173, 117, 54), c(38, 39, 67, 97, 144, 191, 248, 
                233, 209, 180, 149), maxColorValue = 255)), RdYlGn = switch(n - 
            2, rgb(c(252, 255, 145), c(141, 255, 207), c(89, 
            191, 96), maxColorValue = 255), rgb(c(215, 253, 166, 
            26), c(25, 174, 217, 150), c(28, 97, 106, 65), maxColorValue = 255), 
            rgb(c(215, 253, 255, 166, 26), c(25, 174, 255, 217, 
                150), c(28, 97, 191, 106, 65), maxColorValue = 255), 
            rgb(c(215, 252, 254, 217, 145, 26), c(48, 141, 224, 
                239, 207, 152), c(39, 89, 139, 139, 96, 80), 
                maxColorValue = 255), rgb(c(215, 252, 254, 255, 
                217, 145, 26), c(48, 141, 224, 255, 239, 207, 
                152), c(39, 89, 139, 191, 139, 96, 80), maxColorValue = 255), 
            rgb(c(215, 244, 253, 254, 217, 166, 102, 26), c(48, 
                109, 174, 224, 239, 217, 189, 152), c(39, 67, 
                97, 139, 139, 106, 99, 80), maxColorValue = 255), 
            rgb(c(215, 244, 253, 254, 255, 217, 166, 102, 26), 
                c(48, 109, 174, 224, 255, 239, 217, 189, 152), 
                c(39, 67, 97, 139, 191, 139, 106, 99, 80), maxColorValue = 255), 
            rgb(c(165, 215, 244, 253, 254, 217, 166, 102, 26, 
                0), c(0, 48, 109, 174, 224, 239, 217, 189, 152, 
                104), c(38, 39, 67, 97, 139, 139, 106, 99, 80, 
                55), maxColorValue = 255), rgb(c(165, 215, 244, 
                253, 254, 255, 217, 166, 102, 26, 0), c(0, 48, 
                109, 174, 224, 255, 239, 217, 189, 152, 104), 
                c(38, 39, 67, 97, 139, 191, 139, 106, 99, 80, 
                  55), maxColorValue = 255)), Set1 = switch(n - 
            2, rgb(c(228, 55, 77), c(26, 126, 175), c(28, 184, 
            74), maxColorValue = 255), rgb(c(228, 55, 77, 152), 
            c(26, 126, 175, 78), c(28, 184, 74, 163), maxColorValue = 255), 
            rgb(c(228, 55, 77, 152, 255), c(26, 126, 175, 78, 
                127), c(28, 184, 74, 163, 0), maxColorValue = 255), 
            rgb(c(228, 55, 77, 152, 255, 255), c(26, 126, 175, 
                78, 127, 255), c(28, 184, 74, 163, 0, 51), maxColorValue = 255), 
            rgb(c(228, 55, 77, 152, 255, 255, 166), c(26, 126, 
                175, 78, 127, 255, 86), c(28, 184, 74, 163, 0, 
                51, 40), maxColorValue = 255), rgb(c(228, 55, 
                77, 152, 255, 255, 166, 247), c(26, 126, 175, 
                78, 127, 255, 86, 129), c(28, 184, 74, 163, 0, 
                51, 40, 191), maxColorValue = 255), rgb(c(228, 
                55, 77, 152, 255, 255, 166, 247, 153), c(26, 
                126, 175, 78, 127, 255, 86, 129, 153), c(28, 
                184, 74, 163, 0, 51, 40, 191, 153), maxColorValue = 255)), 
        Set2 = switch(n - 2, rgb(c(102, 252, 141), c(194, 141, 
            160), c(165, 98, 203), maxColorValue = 255), rgb(c(102, 
            252, 141, 231), c(194, 141, 160, 138), c(165, 98, 
            203, 195), maxColorValue = 255), rgb(c(102, 252, 
            141, 231, 166), c(194, 141, 160, 138, 216), c(165, 
            98, 203, 195, 84), maxColorValue = 255), rgb(c(102, 
            252, 141, 231, 166, 255), c(194, 141, 160, 138, 216, 
            217), c(165, 98, 203, 195, 84, 47), maxColorValue = 255), 
            rgb(c(102, 252, 141, 231, 166, 255, 229), c(194, 
                141, 160, 138, 216, 217, 196), c(165, 98, 203, 
                195, 84, 47, 148), maxColorValue = 255), rgb(c(102, 
                252, 141, 231, 166, 255, 229, 179), c(194, 141, 
                160, 138, 216, 217, 196, 179), c(165, 98, 203, 
                195, 84, 47, 148, 179), maxColorValue = 255)), 
        Set3 = switch(n - 2, rgb(c(141, 255, 190), c(211, 255, 
            186), c(199, 179, 218), maxColorValue = 255), rgb(c(141, 
            255, 190, 251), c(211, 255, 186, 128), c(199, 179, 
            218, 114), maxColorValue = 255), rgb(c(141, 255, 
            190, 251, 128), c(211, 255, 186, 128, 177), c(199, 
            179, 218, 114, 211), maxColorValue = 255), rgb(c(141, 
            255, 190, 251, 128, 253), c(211, 255, 186, 128, 177, 
            180), c(199, 179, 218, 114, 211, 98), maxColorValue = 255), 
            rgb(c(141, 255, 190, 251, 128, 253, 179), c(211, 
                255, 186, 128, 177, 180, 222), c(199, 179, 218, 
                114, 211, 98, 105), maxColorValue = 255), rgb(c(141, 
                255, 190, 251, 128, 253, 179, 252), c(211, 255, 
                186, 128, 177, 180, 222, 205), c(199, 179, 218, 
                114, 211, 98, 105, 229), maxColorValue = 255), 
            rgb(c(141, 255, 190, 251, 128, 253, 179, 252, 217), 
                c(211, 255, 186, 128, 177, 180, 222, 205, 217), 
                c(199, 179, 218, 114, 211, 98, 105, 229, 217), 
                maxColorValue = 255), rgb(c(141, 255, 190, 251, 
                128, 253, 179, 252, 217, 188), c(211, 255, 186, 
                128, 177, 180, 222, 205, 217, 128), c(199, 179, 
                218, 114, 211, 98, 105, 229, 217, 189), maxColorValue = 255), 
            rgb(c(141, 255, 190, 251, 128, 253, 179, 252, 217, 
                188, 204), c(211, 255, 186, 128, 177, 180, 222, 
                205, 217, 128, 235), c(199, 179, 218, 114, 211, 
                98, 105, 229, 217, 189, 197), maxColorValue = 255), 
            rgb(c(141, 255, 190, 251, 128, 253, 179, 252, 217, 
                188, 204, 255), c(211, 255, 186, 128, 177, 180, 
                222, 205, 217, 128, 235, 237), c(199, 179, 218, 
                114, 211, 98, 105, 229, 217, 189, 197, 111), 
                maxColorValue = 255)), Spectral = switch(n - 
            2, rgb(c(252, 255, 153), c(141, 255, 213), c(89, 
            191, 148), maxColorValue = 255), rgb(c(215, 253, 
            171, 43), c(25, 174, 221, 131), c(28, 97, 164, 186), 
            maxColorValue = 255), rgb(c(215, 253, 255, 171, 43), 
            c(25, 174, 255, 221, 131), c(28, 97, 191, 164, 186), 
            maxColorValue = 255), rgb(c(213, 252, 254, 230, 153, 
            50), c(62, 141, 224, 245, 213, 136), c(79, 89, 139, 
            152, 148, 189), maxColorValue = 255), rgb(c(213, 
            252, 254, 255, 230, 153, 50), c(62, 141, 224, 255, 
            245, 213, 136), c(79, 89, 139, 191, 152, 148, 189), 
            maxColorValue = 255), rgb(c(213, 244, 253, 254, 230, 
            171, 102, 50), c(62, 109, 174, 224, 245, 221, 194, 
            136), c(79, 67, 97, 139, 152, 164, 165, 189), maxColorValue = 255), 
            rgb(c(213, 244, 253, 254, 255, 230, 171, 102, 50), 
                c(62, 109, 174, 224, 255, 245, 221, 194, 136), 
                c(79, 67, 97, 139, 191, 152, 164, 165, 189), 
                maxColorValue = 255), rgb(c(158, 213, 244, 253, 
                254, 230, 171, 102, 50, 94), c(1, 62, 109, 174, 
                224, 245, 221, 194, 136, 79), c(66, 79, 67, 97, 
                139, 152, 164, 165, 189, 162), maxColorValue = 255), 
            rgb(c(158, 213, 244, 253, 254, 255, 230, 171, 102, 
                50, 94), c(1, 62, 109, 174, 224, 255, 245, 221, 
                194, 136, 79), c(66, 79, 67, 97, 139, 191, 152, 
                164, 165, 189, 162), maxColorValue = 255)), YlGn = switch(n - 
            2, rgb(c(247, 173, 49), c(252, 221, 163), c(185, 
            142, 84), maxColorValue = 255), rgb(c(255, 194, 120, 
            35), c(255, 230, 198, 132), c(204, 153, 121, 67), 
            maxColorValue = 255), rgb(c(255, 194, 120, 49, 0), 
            c(255, 230, 198, 163, 104), c(204, 153, 121, 84, 
                55), maxColorValue = 255), rgb(c(255, 217, 173, 
            120, 49, 0), c(255, 240, 221, 198, 163, 104), c(204, 
            163, 142, 121, 84, 55), maxColorValue = 255), rgb(c(255, 
            217, 173, 120, 65, 35, 0), c(255, 240, 221, 198, 
            171, 132, 90), c(204, 163, 142, 121, 93, 67, 50), 
            maxColorValue = 255), rgb(c(255, 247, 217, 173, 120, 
            65, 35, 0), c(255, 252, 240, 221, 198, 171, 132, 
            90), c(229, 185, 163, 142, 121, 93, 67, 50), maxColorValue = 255), 
            rgb(c(255, 247, 217, 173, 120, 65, 35, 0, 0), c(255, 
                252, 240, 221, 198, 171, 132, 104, 69), c(229, 
                185, 163, 142, 121, 93, 67, 55, 41), maxColorValue = 255)), 
        YlGnBu = switch(n - 2, rgb(c(237, 127, 44), c(248, 205, 
            127), c(177, 187, 184), maxColorValue = 255), rgb(c(255, 
            161, 65, 34), c(255, 218, 182, 94), c(204, 180, 196, 
            168), maxColorValue = 255), rgb(c(255, 161, 65, 44, 
            37), c(255, 218, 182, 127, 52), c(204, 180, 196, 
            184, 148), maxColorValue = 255), rgb(c(255, 199, 
            127, 65, 44, 37), c(255, 233, 205, 182, 127, 52), 
            c(204, 180, 187, 196, 184, 148), maxColorValue = 255), 
            rgb(c(255, 199, 127, 65, 29, 34, 12), c(255, 233, 
                205, 182, 145, 94, 44), c(204, 180, 187, 196, 
                192, 168, 132), maxColorValue = 255), rgb(c(255, 
                237, 199, 127, 65, 29, 34, 12), c(255, 248, 233, 
                205, 182, 145, 94, 44), c(217, 177, 180, 187, 
                196, 192, 168, 132), maxColorValue = 255), rgb(c(255, 
                237, 199, 127, 65, 29, 34, 37, 8), c(255, 248, 
                233, 205, 182, 145, 94, 52, 29), c(217, 177, 
                180, 187, 196, 192, 168, 148, 88), maxColorValue = 255)), 
        YlOrBr = switch(n - 2, rgb(c(255, 254, 217), c(247, 196, 
            95), c(188, 79, 14), maxColorValue = 255), rgb(c(255, 
            254, 254, 204), c(255, 217, 153, 76), c(212, 142, 
            41, 2), maxColorValue = 255), rgb(c(255, 254, 254, 
            217, 153), c(255, 217, 153, 95, 52), c(212, 142, 
            41, 14, 4), maxColorValue = 255), rgb(c(255, 254, 
            254, 254, 217, 153), c(255, 227, 196, 153, 95, 52), 
            c(212, 145, 79, 41, 14, 4), maxColorValue = 255), 
            rgb(c(255, 254, 254, 254, 236, 204, 140), c(255, 
                227, 196, 153, 112, 76, 45), c(212, 145, 79, 
                41, 20, 2, 4), maxColorValue = 255), rgb(c(255, 
                255, 254, 254, 254, 236, 204, 140), c(255, 247, 
                227, 196, 153, 112, 76, 45), c(229, 188, 145, 
                79, 41, 20, 2, 4), maxColorValue = 255), rgb(c(255, 
                255, 254, 254, 254, 236, 204, 153, 102), c(255, 
                247, 227, 196, 153, 112, 76, 52, 37), c(229, 
                188, 145, 79, 41, 20, 2, 4, 6), maxColorValue = 255)), 
        YlOrRd = switch(n - 2, rgb(c(255, 254, 240), c(237, 178, 
            59), c(160, 76, 32), maxColorValue = 255), rgb(c(255, 
            254, 253, 227), c(255, 204, 141, 26), c(178, 92, 
            60, 28), maxColorValue = 255), rgb(c(255, 254, 253, 
            240, 189), c(255, 204, 141, 59, 0), c(178, 92, 60, 
            32, 38), maxColorValue = 255), rgb(c(255, 254, 254, 
            253, 240, 189), c(255, 217, 178, 141, 59, 0), c(178, 
            118, 76, 60, 32, 38), maxColorValue = 255), rgb(c(255, 
            254, 254, 253, 252, 227, 177), c(255, 217, 178, 141, 
            78, 26, 0), c(178, 118, 76, 60, 42, 28, 38), maxColorValue = 255), 
            rgb(c(255, 255, 254, 254, 253, 252, 227, 177), c(255, 
                237, 217, 178, 141, 78, 26, 0), c(204, 160, 118, 
                76, 60, 42, 28, 38), maxColorValue = 255), rgb(c(255, 
                255, 254, 254, 253, 252, 227, 189, 128), c(255, 
                237, 217, 178, 141, 78, 26, 0, 0), c(204, 160, 
                118, 76, 60, 42, 28, 38, 38), maxColorValue = 255)))    
	namelist <- names(brewer.colorlist)
    if (return.names) return(namelist)
    if (!(name %in% namelist)) {
        stop(paste(name, "is not a valid palette name for brewer.pal\n"))
    }
    if (n < 3) {
        warning("minimal value for n is 3, returning requested palette with 3 different levels\n")
        return(brewer.pal(3, name))
    }
   # if (n > maxcolors[which(name == namelist)]) {
    #    warning(paste("n too large, allowed maximum for palette", 
     #       name, "is", maxcolors[which(name == namelist)]), 
      #      "\nReturning the palette you asked for with that many colors\n")
       # return(brewer.pal(maxcolors[which(name == namelist)], 
        #    name))
    #}
    #switch(name, brewer.colorlist)
    brewer.colorlist[[name]]
}
# [1] "Accent"   "Blues"    "BrBG"     "BuGn"     "BuPu"     "Dark2"    "GnBu"     "Greens"   "Greys"   
#[10] "Oranges"  "OrRd"     "Paired"   "Pastel1"  "Pastel2"  "PiYG"     "PRGn"     "PuBu"     "PuBuGn"  
#[19] "PuOr"     "PuRd"     "Purples"  "RdBu"     "RdGy"     "RdPu"     "Reds"     "RdYlBu"   "RdYlGn"  
#[28] "Set1"     "Set2"     "Set3"     "Spectral" "YlGn"     "YlGnBu"   "YlOrBr"   "YlOrRd"  

axis.col <- function(side, at, labels, col, ...) {
	for (u in unique(col)) {
		w <- which(col == u)
		axis(side=side, at=at[w], labels=labels[w], col=u, col.axis=u, ...)
	}
}


# f > 1 : brigther
# f < 1 : darker
# f=c(.5,.2,.3,1) for differents changes in red green blue and alpha inclusive
colorize <- function(xcol, f=1, alpha=FALSE, add=0, XOR=FALSE) {
	if (length(xcol) > 1) {
		xl <- length(xcol)
		fl <- length(f)
		al <- length(alpha)
		dl <- length(add)
		rl <- length(XOR)
		xc <- character(xl)
		for (i in 1:xl) xc[i] <- colorize(xcol[i], f[min(fl,i)], alpha[min(al,i)], add[min(dl,i)], XOR[min(rl,i)])
		xc
	} else {
		xf <- function(y,f,add) if (XOR) 255-y*f+add else y*f+add
		if (alpha) apply(col2rgb(xcol,alpha > 0) * c(1,1,1,alpha), 2, function(y) { x=pmin(255,pmax(0,xf(y,f,add))); rgb(x[1],x[2],x[3],x[4], maxColorValue=255) } )
		else apply(col2rgb(xcol), 2, function(y) { x=pmin(255,pmax(0,xf(y,f,add))); rgb(x[1],x[2],x[3], maxColorValue=255) })
	}
}


read.zip <- function(zipfile, internalfile=NA, read.function=read.delim, verbose=TRUE, ...) {
	# function based on http://stackoverflow.com/questions/8986818/automate-zip-file-reading-in-r

	# check the files within zip
	unzfiles <- unzip(zipfile, list=TRUE)
	if (is.na(internalfile) || is.numeric(internalfile)) {
		internalfile <- unzfiles$Name[ifelse(is.na(internalfile),1,internalfile[1])]
	}
    # Create a name for the dir where we'll unzip
    zipdir <- tempfile()
    # Create the dir using that name
    if (verbose) catf("Directory created:",zipdir,"\n")
    dir.create(zipdir)
    # Unzip the file into the dir
    if (verbose) catf("Unzipping file:",internalfile,"...")
    unzip(zipfile, file=internalfile, exdir=zipdir)
    if (verbose) catf("Done!\n")
    # Get the full name of the file
    file <- paste(zipdir, internalfile, sep="/")
    if (verbose) 
    	on.exit({ 
    		catf("Done!\nRemoving temporal files:",file,".\n") 
    		file.remove(file)
    		file.remove(zipdir)
    		}) 
    else
    	on.exit({file.remove(file); file.remove(zipdir);})
    # Read the file
    if (verbose) catf("Reading File...")
    read.function(file, ...)
}

right <- function(s,l) {
	return (substr(s,pmax(1,nchar(s)-l+1),nchar(s)))
}

read.z <- function(zfile, internalfile=NA, read.function=read.delim, quick.mode=TRUE, ...) {
	# THIS FUNCTION READ GZ/ZIP FILES IN MEMORY !!!!
	if (tolower(right(zfile,3)) %in% c(".gz",".xz") ||
		tolower(right(zfile,5)) %in% c(".lzma") ||
		tolower(right(zfile,6)) %in% c(".bzip2")) {
		# GZIP
		zf <- function(x) gzfile(x)
	} else {
		# ZIP
		if (is.na(internalfile) || is.numeric(internalfile)) {
			unzfiles <- unzip(zfile, list=TRUE)
			internalfile <- unzfiles$Name[ifelse(is.na(internalfile),1,internalfile[1])]
		}
		zf <- function(x) unz(x,internalfile)
	}
	cnx <- zf(zfile)
	if (quick.mode) {
		# nothing, read from cnx directly
	} else {
		raw <- readLines(cnx)
		close(cnx)
		cnx <- textConnection(raw)
	}
	dat <- read.function(cnx, ...)
	#if (isOpen(cnx)) try(close(cnx))
	return (dat)
}

# esta implementacion truena con valores altos de n
#perm <- function(n, x) {
#  return(factorial(n) / factorial(n-x))
#}

# esta implementacion truena con valores altos de n
#comb <- function(n, x) {
#  return(factorial(n) / (factorial(x) * factorial(n-x)))
#}

perm <- function(n, x) {
  return(exp(lfactorial(n) - lfactorial(n-x)))
}

comb <- function(n, x) { 
	return ( exp(lfactorial(n) - (lfactorial(x)+lfactorial(n-x))) )
}

#taken from http://druedin.com/2012/08/11/moving-averages-in-r/
# moving average
mav <- function(x,n=5,avoid.na=TRUE){
	y <- filter(x,rep(1/n,n), sides=2)
	if (avoid.na) {
		#VT:adapted to avoid NA
		wna <- which(is.na(y))
		if (length(wna)) {
			y[wna] <- x[wna]
		}	
	}
	y
}


## IN LOG SCALE THIS FUNCTION DOES NOT WORK
text.non.overlap <- function(x, y, labels, adj=c(0.5, 0.5), cex=1, col=1, 
	vfont=NULL, font=NULL, lines=FALSE, precision=3, marginx=1.1, marginy=1.1,
	background=NULL, erase=NULL, frame=NULL, method=c("circles","grid")[1], maxmag=1000,
	avoidpoints=FALSE, avoidx=NULL, avoidy=NULL, avoidjoin=0.30, show.avoid=FALSE, min.distance=0.5,
	verbose=FALSE, minmag=0, min.consistency.distance=0, line.aprox.factor=.2, initial.x=NULL, initial.y=NULL,
	max.time=60) {

	stopifnot(length(x) == length(y) && length(x) == length(labels))

	nona <- which(is.finite(x) & is.finite(y))
	x <- x[nona]
	y <- y[nona]
	labels <- labels[nona]
	if (length(cex) > 1) cex <- cex[nona]
	if (length(col) > 1) col <- col[nona]
	if (length(vfont) > 1) vfont <- vfont[nona]
	if (length(font) > 1) font <- font[nona]
	if (par("xlog")) { x <- log10(x); if (!is.null(avoidx)) avoidx <- log10(avoidx) }
	if (par("ylog")) { y <- log10(y); if (!is.null(avoidy)) avoidy <- log10(avoidy); }

	w <- strwidth(labels, font=font, vfont=vfont)*marginx*cex
	h <- strheight(labels, font=font, vfont=vfont)*marginy*cex

	# set positions
	xp <- x - w*adj[1]
	yp <- y - h*adj[2]
	if (min.distance > 0) {
		theta <- seq(0,2*pi,length=length(x)) #runif(length(x))*2*pi
		xp <- pmax(pmin(xp + cos(theta)*strwidth("W", font=font, vfont=vfont)*marginx*cex*min.distance, par("usr")[2]),par("usr")[1])
		yp <- pmax(pmin(yp + sin(theta)*strheight("W", font=font, vfont=vfont)*marginy*cex*min.distance, par("usr")[4]),par("usr")[3])
	}
	if (!is.null(initial.x)) {
		xp <- circularize(initial.x, length(x))
	}
	if (!is.null(initial.y)) {
		yp <- circularize(initial.y, length(y))
	}
	xxp <- xp + w
	yyp <- yp + h

	start <- as.numeric(format(Sys.time(), "%s"))

	elapsedTime <- function() {
		current <- as.numeric(format(Sys.time(), "%s"))
		return (current-start)
	}

	if (avoidpoints) {
		apxp <- (if (is.null(avoidx)) x else avoidx) - strwidth("x")*cex/2
		apxxp <- apxp + strwidth("x")*cex
		apyp <- (if (is.null(avoidy)) y else avoidy) - strheight("y")*cex/2
		apyyp <- apyp + strheight("y")*cex
		w2 <- which(is.na(apxp) | is.na(apyp))
		if (length(w2)) {
			apxp <- apxp[-w2]
			apxxp <- apxxp[-w2]
			apyp <- apyp[-w2]
			apyyp <- apyyp[-w2]
		}
		# optimize points
		while (is.finite(avoidjoin) && elapsedTime() < max.time) {
			if (verbose) { 
				cat("Avoid Points=",round(sum(!is.na(apxp))*100/length(apxp),1),"        \r")
				flush.console()
			}
			o <- order(apxp)
			changed <- FALSE
			for (i in 1:(length(o)-1)) {
				k <- o[i]
				if (is.na(apxp[k])) {
					break
				}
				for (j in o[(i+1):length(o)]) {
					if (!is.na(apxp[j])) {
						if (apxp[k] >= apxxp[j] || apxp[j] >= apxxp[k]) {
							break
						} else if (apyp[k] >= apyyp[j] || apyp[j] >= apyyp[k]) {
							break
						} else {
							# overlap
							kj <- c(k,j)
							l <- max(apxp[kj])
							r <- min(apxxp[kj])
							b <- max(apyp[kj])
							t <- min(apyyp[kj])
							a <- abs(t-b)*abs(r-l)
							aj <- abs(apxp[j]-apxxp[j])*abs(apyp[j]-apyyp[j])
							ak <- abs(apxp[k]-apxxp[k])*abs(apyp[k]-apyyp[k])
							if (a/max(aj,ak) > avoidjoin) {
								# join
								apxp[k] <- min(apxp[kj], na.rm=TRUE)
								apxxp[k] <- max(apxxp[kj], na.rm=TRUE)
								apyp[k] <- min(apyp[kj], na.rm=TRUE)
								apyyp[k] <- max(apyyp[kj], na.rm=TRUE)
								apxp[j] <- NA # <- apxxp[k] <- apyp[k] <- apyyp[k] <- NA
								changed <- TRUE
							}
						}
					}
				}
			}
			if (!changed) {
				break
			}
		}
		wh <- which(!is.na(apxp))
		apxp <- apxp[wh]
		apxxp <- apxxp[wh]
		apyp <- apyp[wh]
		apyyp <- apyyp[wh]
		if (show.avoid) symbols((apxp+apxxp)/2,(apyp+apyyp)/2,rectangles=data.matrix(data.frame(w=abs(apxxp-apxp),h=abs(apyyp-apyp))),inches=FALSE,add=TRUE)
	}


	# find overlaps
	n <- length(labels)
	ov <- numeric(n)
	ov[] <- 1

	overlap.region <- function(k, x,y,xx,yy) {
		wh <- which(!(x >= xxp | xp >= xx) & !(y >= yyp | yp >= yy))
		if (avoidpoints) {
			wap <- which(!(x >= apxxp | apxp >= xx) & !(y >= apyyp | apyp >= yy))
			return (c(wh[wh != k], -wap))
		}
		wh[wh != k]
	}

	overlaps <- function(k) {
		overlap.region(k, xp[k], yp[k], xxp[k], yyp[k])
	}

	xposs <- seq(par("usr")[1], par("usr")[2], by=strwidth("i",cex=cex)/precision)
	yposs <- seq(par("usr")[3], par("usr")[4], by=strheight("i",cex=cex)/precision)

	method <- match.arg(method, c("circles","grid"))
	#omega <- seq(0, 2*pi, len=72)
	xscale <- abs(par("usr")[2]-par("usr")[1])
	yscale <- abs(par("usr")[4]-par("usr")[3])
	mag <- seq(minmag,1,len=maxmag)[-1]

	pusr <- par("usr") + c(strwidth("x"),-strwidth("x"),strheight("y"),-strheight("y"))*cex
	prevt <- 0
	maxt <- 5
	ronda <- 0
	while (elapsedTime() < max.time) {
		ronda <- ronda + 1
		maxt <- maxt * 5
		wov <- which(ov > 0)
		ov[wov] <- sapply(wov, function(i) length(overlaps(i)))
		if (verbose) { 
			cat("**** Round ",ronda,", Left ",sum(ov > 0),"out of",length(ov),"\n")
			flush.console()
		}
		o <- order(ov)
		overlap <- FALSE
		for (i in o) {
			if (ov[i] > 0) {
				if (verbose) {
					cat("Checking ",labels[i],i,"out of",length(o),"                        \r")
					flush.console()
				}
				xc <- (xp[i]+xxp[i])/2
				yc <- (yp[i]+yyp[i])/2
				w.5 <- w[i]/2
				h.5 <- h[i]/2
				if (method == "circles") {
					#metodo de buscar en ángulos
					assigned <- FALSE
					for (im in 1:length(mag)) {
						omega <- seq(0, 2*pi, len=8+176*im/maxmag)
						ncx <- xc + cos(omega) * mag[im] * xscale
						ncy <- yc + sin(omega) * mag[im] * yscale
						wh <- which(ncx-w.5 > pusr[1] & ncx+w.5 < pusr[2] & ncy-h.5 > pusr[3] & ncy+h.5 < pusr[4])
						#cat(labels[i],length(wh),":")
						#flush.console()
						for (ic in sample(wh)) {
							#cat(ic,"")
							if (length(overlap.region(i, ncx[ic]-w.5, ncy[ic]-h.5, ncx[ic]+w.5, ncy[ic]+h.5)) == 0) {
								assigned <- TRUE
								ov[i] <- 0
								xp[i] <- ncx[ic]-w.5
								yp[i] <- ncy[ic]-h.5
								xxp[i] <- ncx[ic]+w.5
								yyp[i] <- ncy[ic]+h.5
								#text((xp[i]+xxp[i])/2,(yp[i]+yyp[i])/2,labels[i],col=col[min(i,length(col))],cex=cex,vfont=vfont,font=font,adj=c(0.5,0.5))
								break
							}
						}
						#cat("\n")
						if (assigned) {
							break
						}
					}
					if (!assigned) {
						overlap <- TRUE
					}
				} else {
					#metodo de buscar en malla
					wx <- which(xposs <= xp[i] | xposs >= xxp[i])
					wy <- which(yposs <= yp[i] | yposs >= yyp[i])
					if (length(wx) > 0 && length(wy) > 0) {
						overlap <- TRUE
						prev <- c(xp[i], xxp[i], yp[i], yyp[i])
						ox <- order((xc-xposs[wx])^2)
						oy <- order((yc-yposs[wy])^2)
						oxx <- xposs[wx[ox]]
						oyy <- yposs[wy[oy]]
						if (length(oxx) > maxt) {
							oxx <- oxx[1:maxt]
						}
						if (length(oyy) > maxt) {
							oyy <- oyy[1:maxt]
						}
						dx <- (xc-oxx)^2
						dy <- (yc-oyy)^2
						m <- matrix(0, nrow=length(dy), ncol=length(dx))
						for (j in 1:length(dx)) {
							m[,j] <- (dy+dx[j])
						}
						om <- order(m)
						nr <- nrow(m)
						assigned <- FALSE
						for (j in 1:length(om)) {
							rw <- (om[j]-1) %% nr + 1
							cl <- ceiling(om[j] / nr)
							oxi <- oxx[cl]
							oyi <- oyy[rw]
							if (length(overlap.region(i, oxi-w[i]/2, oyi-h[i]/2, oxi+w[i]/2, oyi+h[i]/2)) == 0) {
								assigned <- TRUE
								#cat(oxi,oyi,"")
								break
							}
						}
						#cat(assigned,"\n")
						if (assigned) {
							xp[i] <- oxi-w[i]/2
							xxp[i] <- oxi+w[i]/2
							yp[i] <- oyi-h[i]/2
							yyp[i] <- oyi+h[i]/2
						} else {
						#	xp[i] <- prev[1]
						#	xxp[i] <- prev[2]
						#	yp[i] <- prev[3]
						#	yyp[i] <- prev[4]
						}
					} else {
						cat("No es posible asignar alguna posicion a ",i,"\n")
					}
				}
			}
		}
		# Check consistency
		for (i in 1:length(ov)) {
			if (ov[i] == 0) {
				wh <- which(ov == 0 & (x > x[i] & xp < xp[i] | x < x[i] & xp > xp[i] | y > y[i] & yp < yp[i] | y < y[i] & yp > yp[i]))
				if (length(wh) > 0) {
					ed <- (x[wh]-x[i])^2+(y[wh]-y[i])^2
					if (min(ed) < min.consistency.distance) {
						wh <- wh[which.min(ed)[1]]
						if (verbose) cat("Exchaning ",i,labels[i],"vs",wh,labels[wh],", Distance=",min(ed),"<",min.consistency.distance,"\n")
						pp <- c(xp[i],yp[i],xxp[i],yyp[i])
						xp[i]  <- (xp[wh]+xxp[wh]-w[i])/2
						xxp[i] <- (xp[wh]+xxp[wh]+w[i])/2
						yp[i]  <- (yp[wh]+yyp[wh]-h[i])/2
						yyp[i] <- (yp[wh]+yyp[wh]+h[i])/2
						xp[wh]  <- (pp[1]+pp[3]-w[wh])/2
						xxp[wh] <- (pp[1]+pp[3]+w[wh])/2
						yp[wh]  <- (pp[2]+pp[4]-h[wh])/2
						yyp[wh] <- (pp[2]+pp[4]+h[wh])/2
						ov[i] <- 1
						ov[wh] <- 1
						overlap <- TRUE
					}
				} 
			}
		}
		if (!overlap || ronda > 100) {
			break
		}
		prevt <- maxt
	}

	if (any(lines != 0)) {
		if (length(lines) < length(xp)) lines <- circularize(lines, length(xp))
		if (is.logical(lines)) lines <- (lines * 1)
		laf <- line.aprox.factor
		wsp <- strwidth(" ", cex=cex,vfont=vfont,font=font)*0.333
		hsp <- strheight(" ", cex=cex,vfont=vfont,font=font)*0.1
		.y <- ifelse(yp > y, yp+hsp, ifelse(yyp < y, yyp-hsp, y))
		.x <- ifelse(xp > x, xp+wsp, ifelse(xxp < x, xxp-wsp, x))
		for (i in 1:length(xp)) {
			#.y <- ifelse(yp[i] > y[i], yp[i]+h[i]*laf, ifelse(yyp[i] < y[i], yyp[i]-h[i]*laf, y[i]))
			#.x <- ifelse(xp[i] > x[i], xp[i]+w[i]*laf, ifelse(xxp[i] < x[i], xxp[i]-w[i]*laf, x[i]))
			#if (lines[i]) lines(c(x[i],.x),c(y[i],.y),col=col[min(i,length(col))], lwd=lines[i])
			if (lines[i]) lines(c(x[i],.x[i]),c(y[i],.y[i]),col=col[min(i,length(col))], lwd=lines[i])
		}
	}
	if (!is.null(background)) {
		symbols((xp+xxp)/2,(yp+yyp)/2,rectangles=data.matrix(data.frame(w=w,h=h)),fg=background, bg=background, add=TRUE, inches=FALSE)
	}
	if (!is.null(frame)) {
		symbols((xp+xxp)/2,(yp+yyp)/2,rectangles=data.matrix(data.frame(w=w,h=h)),fg=frame, bg="transparent", add=TRUE, inches=FALSE)
	}
	if (!is.null(erase)) {
		dx <- abs(par("usr")[2]-par("usr")[1])/1000
		dy <- abs(par("usr")[4]-par("usr")[3])/1000
		for (i in 1:length(erase)) {
			text((xp+xxp)/2+rnorm(1,m=0,sd=dx),(yp+yyp)/2+rnorm(1,m=0,sd=dy),labels,col=erase[i],cex=cex,vfont=vfont,font=font,adj=c(0.5,0.5))
		} 
	}
	text((xp+xxp)/2,(yp+yyp)/2,labels,col=col,cex=cex,vfont=vfont,font=font,adj=c(0.5,0.5))
}

# taken from http://stackoverflow.com/questions/12995683/any-way-to-make-plot-points-in-scatterplot-more-transparent-in-r
try(library(scales))
# it contains alpha, a function to control transparency
cat("alpha::scales - Control of transparency\n")

mfrowrect <- function(n) c(trunc(sqrt(n)), ceiling(n/trunc(sqrt(n))))

finites <- function(x) x[is.finite(x)]

change <- function(x, w, r) { x[w] <- r; x }




pubmed.to.data <- function(xfile, records="^[0-9]+\\. ", field="^[^:.]+:", use.period=FALSE) {
	cat("Reading PubMed export file:",xfile,"\n")
	x <- readLines(xfile)
	cat("Looking for records... ")
	rec <- grep(records, x)
	cat(rec,"records found.\n")
	recpos <- c(rec, length(rec)+1)
	fields <- c()
	for (j in 1:2) {
		saving <- (j == 2)
		if (saving) {
			res <- matrix("", nrow=length(rec), ncol=length(fields))
			colnames(res) <- fields
		}
		cat(ifelse(saving,"Storing...","Analyzing..."))
		for (i in 1:length(rec)) {
			if (i %% 100 == 0) cat(i,"")
			irec <- x[recpos[i]:(recpos[i+1]-1)]
			irec[1] <- sub("^[0-9]+\\. ","",irec[1])
			if (use.period) {
				irec <- paste(irec, collapse=" ")
				irec <- unlist(sapply(irec, strsplit, "\\. "))
			}
			w <- 1
			while (length(w)) {
				rex <- regexpr(field, irec)
				rem <- regmatches(irec, rex)
				w <- which(rex <= 0)
				if (length(w)) {
					irec[w-1] <- paste(irec[w-1],irec[w],sep="\n")
					irec <- gsub("\n$","",irec[-w])
				}
			}

			if (saving) {
				remval <- regmatches(irec, rex, invert=TRUE)
				res[i,rem] <- unlist(lapply(remval, function(x) trim(x[2])))
			} else {
				if (!(all(rem) %in% fields)) {
					fields <- unique(c(fields,rem))
				}
			}
		}
		cat("\n")
	}
	colnames(res) <- gsub(":","",colnames(res))
	res
}



# function to perform a grep to filter input lines but with an expression that finalize the concatenation of lines
readLines.and.grep <- function(xfile, grepbegin, grepend=rep("",length(grepbegin)), 
	postreadfunc=trim, collapse="", block=10000, maxLines=-1, filter.func=NULL, header=1, ...) {
	if (length(grep("gz$", xfile))) {
		cnxf <- gzfile(xfile, open="r")
	} else {
		cnxf <- file(xfile, open="r")
	}
	linesHeader <- NULL
	if (header) {
		linesHeader <- readLines(cnxf, n=header)
	}
	L <- block # lines
	db <- list()
	idb <- 0
	xL <- 1:L
	linesRead <- 0
	if (is.null(postreadfunc) || !is.function(postreadfunc)) {
		postreadfunc <- function(x) { x }
	}
	add <- NULL
	wend <- which(nchar(grepend) > 0)
	while (length(xL) == L && (maxLines == -1 || linesRead < maxLines)) {
		cat("[[[***]]] Reading", L, "lines, (current",linesRead,")...\n")
		flush.console()
		xL <- readLines(cnxf, n=L)
		linesRead <- linesRead + length(xL)
		xLS <- postreadfunc(xL)
		if (!is.null(add)) {
			xLS <- c(add, xLS)
			add <- NULL
		}
		maxProc <- length(xLS)
		if (!is.null(grepbegin)) {			
			gbegin <- list()
			gend <- list()
			for (i in 1:length(grepbegin)) {
				gbegin[[i]] <- grep(grepbegin[i], xLS, ...)
				gend[[i]] <- gbegin[[i]]
				if (length(gbegin[[i]])) {
					if  (i %in% wend) {
						gxe <- grep(grepend[i], xLS, ...)
						if (length(gxe)) {
							gbi <- gbegin[[i]]
							geok <- c()
							for (j in 1:length(gbi)) {
								wg <- which(gxe >= gbi[j])
								if (length(wg)) {
									geok[j] <- gxe[wg[1]]
								} else {
									maxProc <- min(maxProc, gbi[j]-1)
									gbegin[[i]] <- gbi[1:(j-1)]
									break
								}
							}
							gend[[i]] <- geok
						} else {
							cat("WARNING: Terminator for expresion ["+grepend[i]+"] was not found within a block. The block needs to be larger. Possible memory issues.\n")
							flush.console()
							maxProc <- 0
							gbegin[[i]] <- c()
							gend[[i]] <- c()
						}
					} else {
						gend[[i]] <- gbegin[[i]]
					}
				}
			}
			if (maxProc) {
				toCopy <- numeric(length(xLS))
				for (i in 1:length(grepbegin)) {
					gb <- gbegin[[i]]
					ge <- gend[[i]]
					if (length(gb) && length(ge)) {
						for (j in 1:length(gb)) {
							if (ge[j] <= maxProc &&  gb[j] <= maxProc) {
								if (gb[j] == ge[j]) {
									toCopy[gb[j]:ge[j]] <- 1
								} else {
									toCopy[gb[j]:ge[j]] <- c(1,rep(-1,ge[j]-gb[j]))
								}
							}
						}
					
					}
				}
				wok <- which(toCopy != 0)
				if (length(wok)) {
					xl <- character(length(wok))
					j <- 0
					for (i in 1:length(wok)) {
						if (toCopy[wok[i]] == 1) {
							j <- j + 1
							xl[j] <- xLS[wok[i]]
						} else {
							xl[j] <- paste(xl[j], xLS[wok[i]], sep=collapse)
						}
					}
					if (j < length(xl)) xl <- xl[1:j]
					db[[length(db)+1]] <- xl
				}
			}
			if (maxProc < length(xLS)) {
				cat("Leaving ",(maxProc+1),":",length(xLS),"for next block.\n")
				flush.console()
				add <- xLS[(maxProc+1):length(xLS)]
			}
		} else if (is.function(filter.function)) {
			isok <- sapply(xLS, filter.function) # debe regresar TRUE/FALSE 1/0
			if (any(!is.na(isok) & isok)) {
				db[[length(db)+1]] <- xLS[!is.na(isok) & isok]
			}
		}
	}
	close(cnxf)
	if (header) {
		list(header=linesHeader, data=unlist(db))
	} else unlist(db)
}

#### Tomado de demo(error.catching)
tryCatch.W.E <- function(expr) {
     W <- NULL
     w.handler <- function(w){ # warning handler
       W <<- w
       invokeRestart("muffleWarning")
     }
     list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                    warning = w.handler),
        warning = W)
 }


arcoiris <- function(n, first="orange", last="red") {
	redgreenblue(c(first,"yellow","cyan","deepskyblue","blue","violet","pink",last), len=n)
}


gz.write.table <- function(dm, file="xf.gz", ...) {
	gz <- gzfile(file, "w")
	write.table(dm, file=gz, ...)
	close(gz)
}

'%nin%' <- Negate('%in%')


## This function 'put' set values in a data frame / matrix / list operating from a named list 
## This is similar to transform and within BUT new variables can be created THAT DEPEND on other
## also created in the same expression if the expression value does not start with "."
## For example if D is a data frame that contains X, Y, and Z columns:
## put(D, A=Y, Y=Y*2, B=Y+A+1, .W=A+1, C=.W-2)
## it will:
## (1) Create A in D with the current values of Y
## (2) Update Y in D with Y*2
## (3) Create B in D
## (4) Create C in D
## NOTE that .W is not created in D because it start with dot "."
put <- function(data, ...) {
	xarg <- match.call(expand.dots=FALSE)
	env <- globalenv() 
	if (is.data.frame(data) || is.matrix(data) || is.array(data)) {
		vars <- colnames(data)
		if (is.matrix(data)) {
			data <- as.data.frame(data)
		}
	} else {
		# List, vector
		vars <- names(data)
	}
	xls <- ls(envir=env)
	# backup previous vars from environment
	restore <- list()
	for (i in vars) {
		if (i %in% xls) {
			restore[[i]] <- get(i, envir=env)
		}
	}
	# Assign vars from data
	for (i in vars) {
		assign(i, data[[i]], envir=env)
	}
	# Assign values to new data
	x3 <- xarg[[3]]
	nx3 <- names(x3) 
	for (j in 1:length(x3)) {
		r <- eval(x3[[j]])
		assign(nx3[j], r, envir=env)
		if (substr(nx3[j],1,1) != ".") {
			data[,nx3[j]] <- r
		}
	}
	# Restore original values in the environment
	for (i in names(restore)) {
		assign(i, restore[[i]], envir=env)
	}
	# return updated data
	data
}

macCopyToClipboard <- function(df, quote=FALSE, sep="\t", ...) {
	cb <- pipe("pbcopy","w")
	write.table(df,cb,quote=quote,sep=sep, ...)
	close(cb)
}

oprintf <- function(s, ..., collapse=",") {
	xl <- list(...)
	for (i in length(xl):1) {
		s <- gsub(paste0("%",i),paste(xl[[i]],collapse=collapse),s)
	}
	s
}

dfapply <- function(df, margin, func) {
	r <- list()
	if (margin==2) {
		for (i in 1:ncol(df)) {
			r[[i]] <- func(df[,i])
		}
		names(r) <- colnames(df)
	} else {
		for (i in 1:nrow(df)) {
			r[[i]] <- func(df[i,])
		}
		names(r) <- rownames(df)
	}
	ml <- unlist(lapply(r, length))
	if (all(ml==ml[1])) {
		if (ml[1] == 1) {
			nr <- names(r)
			r <- unlist(r)
			names(r) <- nr
		} else {
			m <- matrix(NA, nrow=ml[1], ncol=length(r))
			for (i in 1:length(r)) {
				m[,i] <- r[[i]]
			}
			colnames(m) <- names(r)
			r <- m
		}
	}
	r
}


cexHeightAdjustFactor <- 0.91 # Depende del monitor puede ir entre 0.9 y 0.95
cexWidthAdjustFactor <- 0.91  # Depende del monitor

adjustCEX <- function(st, maxCex=3, minCex=0.1, width=par("pin")[1], height=par("pin")[2], hf=1, wf=1, units="inches", ...) {
	ih <- which.max(strheight(st, units=units))[1]
	iw <- which.max(strwidth(st, units=units))[1]
	len <- length(st)
	while (maxCex-minCex > 0.01) {
		midCex <- (maxCex+minCex)/2
		h <- abs(strheight(st[ih], units=units, cex=midCex, ...))
		w <- abs(strwidth(st[iw], units=units, cex=midCex, ...))
		if (h*len > height*hf*cexHeightAdjustFactor || w > width*wf*cexWidthAdjustFactor) {
			maxCex <- midCex
		} else {
			minCex <- midCex
		}
	}
	(maxCex+minCex)/2
}


adjustCEXinches <- function(st, maxCex=3, minCex=0.1, width=par("pin")[1], height=par("pin")[2], hf=1, wf=1, ...) {
	adjustCEX(st, maxCex, minCex, width, height, hf, wf, units="inches")
}


adjustCEXuser <- function(st, maxCex=3, minCex=0.1, width=abs(diff(par("usr")[1:2])), height=abs(diff(par("usr")[3:4])), hf=1, wf=1, ...) {
	adjustCEX(st, maxCex, minCex, width, height, hf, wf, units="user")
}


strinsert <- function(str, what="\n", every, at) {
	if (!missing(every)) {
		gsub(paste0("(.{",every[1],",",every[1],"})"),paste0("\\1",what),str)
	} else if (!missing(at)) {
		what <- circularize(what, length(at))
		nwhat <- c(0,cumsum(nchar(what)))
		xstr <- str
		for (i in 1:length(at)) {
			xstr <- paste0(substr(xstr,1,at[i]+nwhat[i]),what[i],substr(xstr,at[i]+nwhat[i]+1,nchar(xstr)))
		}
		xstr
	} else {
		stop("every or at must be specified.")
	}
}

#ltable(okterms$Term,strsplit(okterms$Genes,", ?"))
# f, list or vector
# l, list (if vector table can be used directly)
ltable <- function(f, l=NULL) {
	if (is.null(l)) {
		table(unlist(f))
	} else if (length(f) != length(l)) {
		stop("Lengths of f and l are different (",length(f)," vs ",length(l),").")
	}
	ul <- unique(unlist(l))
	uf <- unique(f)
	m <- matrix(0, nrow=length(uf), ncol=length(ul))
	colnames(m) <- ul
	rownames(m) <- uf
	for (i in uf) {
		w <- which(f==i)
		xt <- table(unlist(l[w]))
		m[i, names(xt)] <- xt[]
	}
	m
}

# Get Sequence named cont from FASTA file fa reading by blocks of size block lines
getSeqFromFasta <- function(fa, cont, block=100000, suffix=" ",verbose=FALSE) {
	xfile <- file(fa, open="r")
	f <- 0
	xexp <- paste0("^>",cont,suffix)
	n <- 0
	seq <- list()
	if (verbose) {
		catf("Searching `", xexp, "` in file `",fa,"`...")
	}
	read <- TRUE
	nlines <- block
	while (length(seq) < length(cont) && nlines == block) {
		if (read) {
			lines <- readLines(xfile, n=block)
			nlines <- length(lines)
			n <- n + nlines
			if (verbose) {
				catf(".")
			}
		}
		read <- TRUE
		lg <- sapply(xexp, grep, lines) #grep(xexp, lines)
		ll <- unlist(lapply(lg, length))
		if (any(ll > 0)) {
			wl <- which(ll > 0)
			xg <- lg[[wl[1]]]
			f <- n + xg[1]
			if (xg[1] < length(lines)) {
				lines <- lines[(xg[1]+1) : length(lines)]
			} else {
				lines <- c()
			}
			xseq <- c()
			while (nlines > 0) {
				xg <- grep("^>",lines)
				if (length(xg)) {
					xseq <- c(xseq, lines[1:(xg[1]-1)])
					lines <- lines[xg[1]:length(lines)]
					read <- FALSE
					break
				} else {
					xseq <- c(xseq, lines)
					lines <- readLines(xfile, n=block)
					nlines <- length(lines)
				}
			}
			s <- paste(xseq, collapse="")
			attr(s, "line") <- f
			seq[[cont[wl]]] <- s
			if (verbose) {
				catf(length(seq),":",cont[wl],"")
			}
		}
	}
	if (verbose) {
		catf("Done.\n")
	}
	close(xfile)
	seq
}


getIndexFromFasta <- function(fa, block=100e6, del=" .+", verbose=FALSE) {
	read <- TRUE
	bytes <- 0
	n <- 0
	nreads <- 0
	cnx <- file(fa, open="rb")
	if (verbose) {
		catf("Reading [",fa,"]... ")
	}
	fnd <- data.frame(n=0, id="", start=0, end=0, stringsAsFactors=FALSE)
	while (TRUE) {
		if (read) {
			n <- n + bytes
			buf <- readChar(cnx, n=block, useBytes=TRUE)
			bytes <- nchar(buf)
			nreads <- nreads + 1
			if (verbose) {
				catf(nreads,paste0("{",nrow(fnd)-1,"}"),"")
			}
		}
		pos <- grepRaw(">", buf, all=TRUE)
		if (length(pos)) {
			if (nrow(fnd) > 1) {
				if (fnd$end[nrow(fnd)] == 0) {
					fnd$end[nrow(fnd)] <- n+pos[1]
				}
			}
			ids <- gsub(del,"",sapply(pos, function(k) substr(buf, k+1, k+100)))
			if (pos[length(pos)]+50 > bytes) {
				ids[length(pos)] <- ""
			}
			xfnd <- data.frame(n=1:length(pos)+nrow(fnd)-1, id=ids, start=n+pos, end=c((n+pos-1)[-1],0), stringsAsFactors=FALSE)
			fnd <- rbind(fnd, xfnd)
		}
		if (bytes < block) {
			break
		}
	}
	close(cnx)
	n <- n + bytes
	if (nrow(fnd) > 1) {
		if (fnd$end[nrow(fnd)] == 0) {
			fnd$end[nrow(fnd)] <- n
		}
		fnd <- fnd[-1,]
		# process empty id in fnd
		w <- which(nchar(fnd$id)==0)
		if (length(w)) {
			if (verbose) catf("\nFilling ids : ", w, "... ")
			cnx <- file(fa, open="rb")
			for (i in w) {
				seek(cnx, fnd[i, start]+1)
				buf <- readChar(cnx, n=100, useBytes=TRUE)
				fnd$id[i] <- gsub(del,"",buf)
				if (verbose) catf(i,"")
			}
			close(cnx)
		}
		rownames(fnd) <- fnd$id
		return (fnd)
	} else {
		return (NULL)
	}
}


getSequencesFromIndexedFasta <- function(fa, idx, seqid, removeName=TRUE) {
	seq <- c()
	cnx <- file(fa, open="rb")
	for (sid in seqid) {
		s <- ""
		rec <- idx[sid, ]
		if (nrow(rec)) {
			seek(cnx, rec$start-1)
			buf <- readChar(cnx, n=rec$end-rec$start+1, useBytes=TRUE)
			s <- strsplit(buf,"\n")[[1]]
			if (removeName) s <- s[-1]
		}
		seq[sid] <- paste(s, collapse="")
	}
	close(cnx)
	seq
}



set.names <- function(x, idxs, xnames, funcnames=ifelse(is.matrix(x) || is.data.frame(x), colnames, names) ) {
	stopifnot(length(idxs) == length(xnames))
	cn <- funcnames(x)
	cn[idxs] <- xnames
	funcnames(x) <- cn
	x
}

set.colnames <- function(x, idxs, xnames) {
	stopifnot(length(idxs) == length(xnames))
	cn <- colnames(x)
	if (is.null(cn)) cn <- paste0("col",1:ncol(x))
	cn[idxs] <- xnames
	colnames(x) <- cn
	x
}


set.rownames <- function(x, idxs, xnames) {
	stopifnot(length(idxs) == length(xnames))
	cn <- rownames(x)
	if (is.null(cn)) cn <- paste0("row",1:nrow(x))
	cn[idxs] <- xnames
	rownames(x) <- cn
	x
}

set.names <- function(x, idxs, xnames) {
	stopifnot(length(idxs) == length(xnames))
	cn <- names(x)
	if (is.null(cn)) cn <- paste0("name",1:length(x))
	cn[idxs] <- xnames
	names(x) <- cn
	x
}

cumsumreset <- function(x, v=0, z=0) {
	cs <- numeric(length(x))
	s <- z
	for (i in 1:length(x)) {
		if (x[i]==v) {
			s <- z
		}
		s <- s + x[i]
		cs[i] <- s
	}
	cs
}

cumonequal <- function(x, by=1) {
	cumsumreset(c(0,(1-(diff(smut) > 0)))) * by
}


# Order cols of mat by applying a function f1 (f1 must return 'which' is the best to be ordered)
order.cols.by.func <- function(mat, f1=which.min, f2.decreasing=FALSE, level=0, max.level=1000, verbose=FALSE) {
	#cat(level)
	if (nrow(mat) < 2) order(mat, decreasing=f2.decreasing)
	wm <- apply(mat, 2, f1)
	wmt <- sort(table(wm),decreasing=TRUE)
	o <- numeric(length(wm))
	prev <- 0
	for (j in names(wmt)) {
		i <- as.numeric(j)
		w <- which(wm == i)
		lw <- length(w)
		xt <- table(mat[i,w])
		xo <- order(mat[i,w], decreasing=f2.decreasing)
		xi <- 1:lw+prev
		o[xi] <- w[xo] ## default order
		if (any(xt > 1) && max.level > level && nrow(mat) > 1) {
			w <- w[xo]
			xu <- unique(mat[i,w])
			for (u in xu) {
				wu <- which(mat[i,w] == u)
				if (length(wu) > 1) {
					#cat(".")
					if (verbose) cat("Level:",level,", sorting:",length(wu),", value:",u,", from:",length(w),", given:",ncol(mat),", row:",i,"\n")
					yo <- order.cols.by.func(mat[-i,w[wu],drop=FALSE], f1, f2.decreasing, level=level+1, max.level)
					#cat(level)
					#print(dim(mat[-i,w[wu],drop=FALSE]))
					w[wu] <- w[wu[yo]]
				#} else {
				#	o[prev+w[wu]] 
				}
			}
			o[xi] <- w
		}
		prev <- prev + lw
	}
	#cat("\n")
	o
}


#applying switch for a vector
swapply <- function(x, ...) {
	unlist(sapply(x, function(y) switch(y, ...)))
}


l2m <- function(xlist) {
	rows <- length(xlist)
	cols <- length(xlist[[1]])
	m <- matrix(NA, ncol=cols, nrow=rows)
	for (i in 1:length(xlist)) {
		m[i,] <- xlist[[i]]
	}
	rownames(m) <- names(xlist)
	colnames(m) <- names(xlist[[1]])
	m
}



#left & right are dendrograms
reorder.leafs <- function(left, right, orders, distmat, func=mean) {
	lo <- order.dendrogram(left)
	ro <- order.dendrogram(right)
	R <- length(ro)
	L <- length(lo)
	O <- length(orders)
	lFirst <- orders[1] == lo[1]
	rLast <- orders[O] == ro[R]
	#prev neighbours left
	pnlo <- list(c(lo[L], ro[1]))
	if (!lFirst) { # not first
		wlo <- which(orders == lo[1])
		pnlo[[2]] <- c(lo[1], orders[wlo-1])
	}
	#prev neighbours right
	wro <- which(ro[1] == orders)
	pnro <- list(c(ro[1], orders[wro-1]))
	if (!rLast) { # not last
		wro <- which(ro[R] == orders)
		pnro[[2]] <- c(ro[R], orders[wro+1])
	}
	#new neighbours left
	nnlo <- list(c(lo[1], ro[R]))
	if (!rLast) { # not last
		nnlo[[2]] <- c(lo[L], pnro[[2]][2])
	}
	nnro <- list()
	if (!lFirst) {
		nnro <- list(c(ro[1], pnlo[[2]][2]))
	}
	pno <- c(pnlo, pnro)
	nno <- c(nnlo, nnro)
	#cat("lo:",lo,"\n")
	#cat("ro:",ro,"\n")
	#print(pno)
	#print(nno)
	preval <- func(unlist(lapply(pno, function(x) distmat[x[1],x[2]])))
	newval <- func(unlist(lapply(nno, function(x) distmat[x[1],x[2]])))
	if (is.na(preval)) preval <- 0
	if (is.na(newval)) newval <- 0
	#print(c(preval, newval))
	return (newval < preval)
}

reorder.subtree <- function(ddg, ord, distmat, func=mean, method="top-down") {
	l <- ddg[[1]]
	r <- ddg[[2]]
	if (method == "top-down") {
		if (!is.leaf(l)) l <- reorder.subtree(l, ord, distmat, mean, method)
		if (!is.leaf(r)) r <- reorder.subtree(r, ord, distmat, mean, method)
		ddg <- merge(l, r, height=attr(ddg,"height"))
	}
	reo <- reorder.leafs(l, r, ord, distmat, func)
	if (reo) {
		ddg <- merge(r, l)
		x <- r
		r <- l
		l <- x
	}
	if (method == "bottom-up") {
		if (!is.leaf(l)) l <- reorder.subtree(l, ord, distmat, mean, method)
		if (!is.leaf(r)) r <- reorder.subtree(r, ord, distmat, mean, method)
		ddg <- merge(l, r, height=attr(ddg,"height"))
	}
	ddg
}

reorder.tree <- function(ddg, distmat, func=mean, method=c("top-down","bottom-up")[1], iterations=20) {
	repeat {
		ord <- order.dendrogram(ddg)
		ddg <- reorder.subtree(ddg, ord, distmat, func, method)
		iterations <- iterations - 1
		if (any(ord != order.dendrogram(ddg)) && iterations > 0) break
	} 
	ddg
}

#data <- t(scale(t(data.lq[markers,])))
#dm <- dist(data)
#xdm <- as.matrix(dm)
#hc <- hclust(dm, method="ward")
#hcr <- as.dendrogram(hc)
#reorder.leafs(hcr[[1]], hcr[[2]], order.dendrogram(hcr), xdm, mean)
#reorder.leafs(hcr[[1]][[1]], hcr[[1]][[2]], order.dendrogram(hcr), xdm, mean)
#hcr2 <- reorder.tree(hcr, xdm)


# Based on pie but also considering lwd, add, labels.col
plot.pie <- function (x, labels = names(x), edges = 200, radius = 0.8, clockwise = FALSE, 
    init.angle = if (clockwise) 90 else 0, density = NULL, angle = 45, 
    col = NULL, border = NULL, lty = NULL, main = NULL, lwd=NULL, add=FALSE, labels.col=NULL, 
    offset.x=0, offset.y=0, radiusX=radius, radiusY=radius, ...) 
{
    if (!is.numeric(x) || any(is.na(x) | x < 0)) 
        stop("'x' values must be positive.")
    if (is.null(labels)) 
        labels <- as.character(seq_along(x))
    else labels <- as.graphicsAnnot(labels)
    x <- c(0, cumsum(x)/sum(x))
    dx <- diff(x)
    nx <- length(dx)
    marg <- par("mar")
    if (!add) plot.new()
    pin <- par("pin")
    xlim <- ylim <- c(-1, 1)
    if (pin[1L] > pin[2L]) 
        xlim <- (pin[1L]/pin[2L]) * xlim
    else ylim <- (pin[2L]/pin[1L]) * ylim
    dev.hold()
    on.exit({ par(mar=marg); dev.flush() })
    if (!add) plot.window(xlim, ylim, "", asp = 1)
    if (!add) par(mar=c(0,0,0,0))
    if (is.null(col)) 
        col <- if (is.null(density)) 
            c("white", "lightblue", "mistyrose", "lightcyan", 
                "lavender", "cornsilk")
        else par("fg")
    if (!is.null(col)) 
        col <- rep_len(col, nx)
    if (!is.null(border)) 
        border <- rep_len(border, nx)
    if (!is.null(lty)) 
        lty <- rep_len(lty, nx)
    if (!is.null(lwd)) 
        lwd <- rep_len(lwd, nx)
    labels.col <- rep_len(if (!is.null(labels.col)) labels.col else col, nx)
    angle <- rep(angle, nx)
    if (!is.null(density)) 
        density <- rep_len(density, nx)
    twopi <- if (clockwise) 
        -2 * pi
    else 2 * pi
    t2xy <- function(t) {
        t2p <- twopi * t + init.angle * pi/180
        list(x = radiusX * cos(t2p), y = radiusY * sin(t2p))
    }
    x_ <- 1:nx
    y_ <- 1:nx
    for (i in 1L:nx) {
        n <- max(2, floor(edges * dx[i]))
        P <- t2xy(seq.int(x[i], x[i + 1], length.out = n))
        polygon(c(P$x, 0)+offset.x, c(P$y, 0)+offset.y, density = density[i], angle = angle[i], 
            border = border[i], col = col[i], lty = lty[i], lwd=lwd[i])
        P <- t2xy(mean(x[i + 0:1]))
        lab <- as.character(labels[i])
        if (!is.na(lab) && nzchar(lab)) {
            lines(c(1, 1.05) * P$x+offset.x, c(1, 1.05) * P$y+offset.y)
            text(1.1 * P$x + offset.x, 1.1 * P$y + offset.y, labels[i], xpd = TRUE, 
                adj = ifelse(P$x < 0, 1, 0), col=labels.col[i], ...)
        }
        x_[i] <- P$x+offset.x
        y_[i] <- P$y+offset.y
    }
    title(main = main, ...)
    invisible(list(x=x_, y=y_, angle=(twopi*(x[-1]+x[-length(x)])/2+init.angle * pi/180)*180/pi))
}


# Dibuja donas jerarquicas... las mas externas son más detalladas y las mas internas las mas generales
# d2draw - columnas de todas las donas a dibujar
# c2draw - columna con el valor de las donas
# c2chk  - columnas a considerar para dibujar, esto sirve para que si se van a dibujar varias donas, sean comparables en rubros y colores
# fademode - level : color fading is used per level of the dona, root, color fading is per first dona at the center
draw.doughnut <- function(df, c2chk, c2draw=1, d2draw, hierachy=rev(d2draw), hole=1, radius=rep(1,length(d2draw)), col=1:length(hList[[1]]), cex=1, lwd=NULL, lty=NULL, border=1, labels="$L\n$%", minimum=1, fadeto="white", fademode=c("level","root")[1], total=100, maxcolperc=70, expfactor=0.03, labels.col="white", fadelen=100, print=TRUE) {

	nx <- length(d2draw)
	if (!is.null(border)) border <- rep_len(border, nx)
	if (!is.null(lty)) lty <- rep_len(lty, nx)
	if (!is.null(lwd)) lwd <- rep_len(lwd, nx)
	if (!is.null(cex)) cex <- rep_len(cex, nx)
	if (!is.null(labels.col)) labels.col <- rep_len(labels.col, nx)

	if (is.numeric(c2draw)) c2draw <- c2chk[c2draw]

	#df$KEY <- paste(df[,hierachy,drop=FALSE],sep=".")
	#tapply(df[c2draw,], df$KEY)

	p2c <-function(p) fadelen*(1-exp(-expfactor*(p)))

	hList <- list()
	tList <- list()
	cList <- list()
	others <- c()
	for (i in 1:length(hierachy)) {
		h <- hierachy[i]
		k <- paste0("KEY__",i)
		df[,k] <- apply(df[, hierachy[1:i],drop=FALSE],1,paste,collapse="/")
		ts <- sapply(c2chk, function(column) tapply(df[,column], df[,k], sum, na.rm=TRUE))
		ts <- cbind(ts,TO_DRAW=apply(ts >= minimum,1,any))
		#print(ts)
		tList[[i]] <- ts
		l <- ts[ts[,"TO_DRAW"]==1,c2draw]
		if (i==1 && sum(l) != total) others <- total-sum(l)
		if (others) l <- c(l, Others=others)
		if (i > 1) {
			hlant <- hList[[i-1]]
			kant <- paste0("KEY__",i-1)
			kroot <- paste0("KEY__",1)
			l <- c()
			ch <- c()
			for (j in names(hlant)) {
				dfok <- df[df[,kant] == j,]
				r <- tapply(dfok[,c2draw], dfok[,k], sum, na.rm=TRUE)
				r <- r[r > minimum]
				if (abs(sum(r)-hlant[j]) > 1e-6) { 
					r <- c(r, round(hlant[j]-sum(r),6))
					names(r) <- c(names(r)[-length(r)], paste0(j, ifelse(extract("Others", j)=="Others",""," (Others)")))
				}
				r <- r[order(r,decreasing=TRUE)]
				l <- c(l, r)
				if (i > 1) {
					if (fademode == "level") xcol <- redgreenblue(c(cList[[i-1]][j], fadeto[i]), len=fadelen)
					if (fademode == "root") {
						if (nrow(dfok) == 0) ft <- length(fadeto)
						else ft <- which(names(cList[[1]])==dfok[1,kroot])
						xcol <- redgreenblue(c(cList[[i-1]][j], fadeto[ft]), len=fadelen)
					}
					#print(10*r/(hList[[i-1]][j]))
					#print(xcol)
					ch <- c(ch, xcol[pmin(maxcolperc,pmax(1,fadelen-p2c(fadelen*r/hlant[j])))])
				}
			}
			cList[[i]] <- ch
			#print(ch)
		}
		hList[[i]] <- l
		#print(l)
		if (i == 1) {
			if (!is.null(col)) col <- rep_len(col, length(hList[[i]]))
			if (!is.null(fadeto)) fadeto <- rep_len(fadeto, length(col))
			#print(col)
			#print(fadeto)
			cList[[i]] <- apply(col2rgb(col)/255,2,function(x) rgb(x[1],x[2],x[3]))
		}
		#print(cList[[i]])
		#print(hList[[i]])
		names(cList[[i]]) <- names(hList[[i]])
	}

	if(sum(radius) > 1) {
		sr <- sum(c(radius,hole))
		radius <- rev(cumsum(rev(c(radius,hole))))/sr
		hole <- hole/sr
	}
	for (i in 1:length(d2draw)) {
		c2d <- d2draw[i]
		#r <- sum(radius[i:nx])
		w <- which(hierachy==c2d)
		tod <- hList[[w]]
		#print(tod)
		#pp <- plot.pie(tod, labels=if (i==1) paste0(gsub("[^/]+/","",names(tod)),"\n",tod,"%") else NA, radius=radius[i], add=i > 1, col=cList[[w]], border=border[i], lwd=lwd[i], lty=lty[i], cex=cex[i])
		if (print) {
			print(tod)
		}
		pp <- plot.pie(tod, labels=NA, radius=radius[i], add=i > 1, col=cList[[w]], border=border[i], lwd=lwd[i], lty=lty[i], cex=cex[i])
		for (k in 1:length(pp$x)) {
			a <- pp$angle[k]
			#|| tod[k] > minimum/2 && i == 1
			if (tod[k] > minimum) text(pp$x[k]*0.99, pp$y[k]*.99, paste0(round(tod[k],2),"%"), col=labels.col[i], cex=cex[i], srt=ifelse(a > 90 && a < 270,a+180,a), adj=c(ifelse(a > 90 && a < 270,0,1),0.5))
		}
		if (i == 1) {
			for (k in 1:length(pp$x)) {
				a <- pp$angle[k]
				#if (tod[k] > minimum/2) text(pp$x[k]*1.01, pp$y[k]*1.01, paste0(gsub("[^/]+/","",names(tod[k]))," ",round(tod[k],2),"%"), col=if (i==1) cList[[w]][k] else "white", cex=cex[i], srt=ifelse(a > 90 && a < 270,a+180,a), adj=c(ifelse(a > 90 && a < 270,1,0),0.5))
				if (tod[k] > minimum) text(pp$x[k]*1.01, pp$y[k]*1.01, paste0(gsub("[^/]+/","",names(tod[k]))), col=if (i==1) cList[[w]][k] else "white", cex=cex[i], srt=ifelse(a > 90 && a < 270,a+180,a), adj=c(ifelse(a > 90 && a < 270,1,0),0.5)) 
				else text(pp$x[k]*1.01, pp$y[k]*1.01, paste0(ifelse(a > 90 && a < 270,paste0(round(tod[k],2),"% "),""),gsub("[^/]+/","",names(tod[k])),ifelse(a > 90 && a < 270,"",paste0(" ",round(tod[k],2),"%"))), col=if (i==1) cList[[w]][k] else labels.col[i], cex=cex[i], srt=ifelse(a > 90 && a < 270,a+180,a), adj=c(ifelse(a > 90 && a < 270,1,0),0.5)) 
				#print(tod[k])
			}
		}
	}
	if (hole > 0) symbols(0,0,circles=hole,add=TRUE,fg="white",bg="white", inches=FALSE)
	i = length(d2draw)

}

invertColor <- function(col) {
	kol <- 255-col2rgb(col, alpha=TRUE)
	apply(kol, 2, function(x) rgb(x[1],x[2],x[3],alpha=255-x[4],maxColorValue=255))
}

cat("https://www.nceas.ucsb.edu/~frazier/RSpatialGuides/colorPaletteCheatsheet.pdf\n")
cat('library("colorspace")\n')
cat('pal <- choose_palette()\n')
cat('library(RColorBrewer)\n')
cat('display.brewer.all()\n')


