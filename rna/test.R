
library(scde)
library(parallel)
# load example dataset
data(es.mef.small)
# factor determining cell types
sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(es.mef.small)), levels = c("ESC", "MEF"))
# the group factor should be named accordingly
names(sg) <- colnames(es.mef.small)  
table(sg)

# clean up the dataset
cd <- es.mef.small
# omit genes that are never detected
cd <- cd[rowSums(cd)>0, ]
# omit cells with very poor coverage
cd <- cd[, colSums(cd)>1e4]

# EVALUATION NOT NEEDED
# calculate models
o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE, save.model.plots = FALSE, verbose = 1)
devtools::use_data(o.ifm)  # save for later since this step takes a long time

data(o.ifm)

# filter out cells that don't show positive correlation with
# the expected expression magnitudes (very poor fits)
valid.cells <- o.ifm$corr.a > 0
table(valid.cells)
o.ifm <- o.ifm[valid.cells, ]

# estimate gene expression prior
o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)

# define two groups of cells
groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels  =  c("ESC", "MEF"))
names(groups) <- row.names(o.ifm)
# run differential expression tests on all genes.
ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups  =  groups, n.randomizations  =  100, n.cores  =  1, verbose  =  1)

# top upregulated genes (tail would show top downregulated ones)
head(ediff[order(ediff$Z, decreasing  =  TRUE), ])

# write out a table with all the results, showing most significantly different genes (in both directions) on top
write.table(ediff[order(abs(ediff$Z), decreasing = TRUE), ], file = "results.txt", row.names = TRUE, col.names = TRUE, sep = "\t", quote = FALSE)

scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior)


batch <- as.factor(ifelse(rbinom(nrow(o.ifm), 1, 0.5) == 1, "batch1", "batch2"))
# check the interaction between batches and cell types (shouldn't be any)
table(groups, batch)

# test the Tdh gene again
scde.test.gene.expression.difference("Tdh", models = o.ifm, counts = cd, prior = o.prior, batch = batch)

# test for all of the genes
ediff.batch <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, batch = batch, n.randomizations = 100, n.cores = 1, return.posteriors = TRUE, verbose = 1)


# calculate joint posterior for ESCs (set return.individual.posterior.modes=T if you need p.modes)
jp <- scde.posteriors(models = o.ifm[grep("ESC",rownames(o.ifm)), ], cd, o.prior, n.cores = 1)

# get expression magntiude estimates
o.fpm <- scde.expression.magnitude(o.ifm, counts = cd)
# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
par(mfrow = c(1,1), mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1)
plot(c(), c(), xlim=range(o.prior$x), ylim=c(0,1), xlab="expression magnitude (log10)", ylab="drop-out probability")
invisible(apply(o.fail.curves[, grep("ES",colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y,col = "orange")))
invisible(apply(o.fail.curves[, grep("MEF", colnames(o.fail.curves))], 2, function(y) lines(x = o.prior$x, y = y, col = "dodgerblue")))

# get failure probabilities on the expresison range
o.fail.curves <- scde.failure.probability(o.ifm, magnitudes = log((10^o.prior$x)-1))
# get self-fail probabilities (at a given observed count)
p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)


p.self.fail <- scde.failure.probability(models = o.ifm, counts = cd)
# simulate drop-outs
# note: using 10 sampling rounds for illustration here. ~500 or more should be used.
n.simulations <- 10; k <- 0.9;
cell.names <- colnames(cd); names(cell.names) <- cell.names;
dl <- mclapply(1:n.simulations,function(i) {
  scd1 <- do.call(cbind,lapply(cell.names,function(nam) {
    x <- cd[,nam];
    # replace predicted drop outs with NAs
    x[!as.logical(rbinom(length(x),1,1-p.self.fail[,nam]*k))] <- NA;
    x;
  }))
  rownames(scd1) <- rownames(cd); 
  # calculate correlation on the complete observation pairs
  cor(log10(scd1+1),use="pairwise.complete.obs");
}, mc.cores = 1)
# calculate average distance across sampling rounds
direct.dist <- as.dist(1-Reduce("+",dl)/length(dl))



# load boot package for the weighted correlation implementation
require(boot)
k <- 0.95;
reciprocal.dist <- as.dist(1 - do.call(rbind, mclapply(cell.names, function(nam1) {
  unlist(lapply(cell.names, function(nam2) {
    # reciprocal probabilities
    f1 <- scde.failure.probability(models = o.ifm[nam1,,drop = FALSE], magnitudes = o.fpm[, nam2])
    f2 <- scde.failure.probability(models = o.ifm[nam2,,drop = FALSE], magnitudes = o.fpm[, nam1])
    # weight factor
    pnf <- sqrt((1-f1)*(1-f2))*k +(1-k); 
    boot::corr(log10(cbind(cd[, nam1], cd[, nam2])+1), w = pnf)
  }))
},mc.cores = 1)), upper = FALSE)


# reclculate posteriors with the individual posterior modes 
jp <- scde.posteriors(models = o.ifm, cd, o.prior, return.individual.posterior.modes = TRUE, n.cores = 1)
# find joint posterior modes for each gene - a measure of MLE of group-average expression
jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
p.mode.fail <- scde.failure.probability(models = o.ifm, magnitudes = jp$jp.modes)
# weight matrix
matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
# magnitude matrix (using individual posterior modes here)
mat <- log10(exp(jp$modes)+1);
# weighted distance
mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
  unlist(lapply(cell.names,function(nam2) {
    corr(cbind(mat[, nam1], mat[, nam2]), w = sqrt(sqrt(matw[, nam1]*matw[, nam2])))
  }))
}, mc.cores = 1)), upper = FALSE)



############### Pathway and Gene Set Overdispersion Analysis

# read in the expression matrix
data(pollen)
cd <- pollen
# filter data
# filter out low-gene cells
cd <- cd[, colSums(cd>0)>1.8e3]
# remove genes that don't have many reads
cd <- cd[rowSums(cd)>10, ]
# remove genes that are not seen in a sufficient number of cells
cd <- cd[rowSums(cd>0)>5, ]
# check the final dimensions of the read count matrix
dim(cd)

x <- gsub("^Hi_(.*)_.*", "\\1", colnames(cd))
l2cols <- c("coral4", "olivedrab3", "skyblue2", "slateblue3")[as.integer(factor(x, levels = c("NPC", "GW16", "GW21", "GW21+3")))]


# EVALUATION NOT NEEDED
knn <- knn.error.models(cd, k = ncol(cd)/4, n.cores = 1, min.count.threshold = 2, min.nonfailed = 5, max.model.plots = 10)
devtools::use_data(knn)  # save for later since this step takes a long time

data(knn)

varinfo <- pagoda.varnorm(knn, counts = cd, trim = 3/ncol(cd), max.adj.var = 5, n.cores = 1, plot = TRUE)
# list top overdispersed genes
sort(varinfo$arv, decreasing = TRUE)[1:10]
varinfo <- pagoda.subtract.aspect(varinfo, colSums(cd[, rownames(knn)]>0))


library(org.Hs.eg.db)
# translate gene names to ids
ids <- unlist(lapply(mget(rownames(cd), org.Hs.egALIAS2EG, ifnotfound = NA), function(x) x[1]))
# reverse map
rids <- names(ids)
names(rids) <- ids
# list all the ids per GO category
go.env <- eapply(org.Hs.egGO2ALLEGS, function(x) as.character(na.omit(rids[x])))
# omit categories with too few genes
go.env <- go.env[unlist(lapply(go.env, length))>5]

# append descriptions to the GO names
library(GO.db)
desc <- unlist(lapply(mget(names(go.env), GOTERM, ifnotfound = NA), function(x) if(is.logical(x)) { return("") } else { slot(x, "Term")}))
names(go.env) <- paste(names(go.env), desc)  # append description to the names

go.env <- list2env(go.env)  # convert to an environment

pwpca <- pagoda.pathway.wPCA(varinfo, go.env, n.components = 1, n.cores = 1, n.internal.shuffles = 0)
df <- pagoda.top.aspects(pwpca, return.table = TRUE, plot = TRUE, z.score = 1.96)

head(df)

clpca <- pagoda.gene.clusters(varinfo, trim = 7.1/ncol(varinfo$mat), n.clusters = 150, n.cores = 1, plot = TRUE)
df <- pagoda.top.aspects(pwpca, clpca, return.table = TRUE, plot = TRUE, z.score = 1.96)
head(df)

# get full info on the top aspects
tam <- pagoda.top.aspects(pwpca, clpca, n.cells = NULL, z.score = qnorm(0.01/2, lower.tail = FALSE))
# determine overall cell clustering
hc <- pagoda.cluster.cells(tam, varinfo)

tamr <- pagoda.reduce.loading.redundancy(tam, pwpca, clpca)

tamr2 <- pagoda.reduce.redundancy(tamr, distance.threshold = 0.9, plot = TRUE, cell.clustering = hc, labRow = NA, labCol = NA, box = TRUE, margins = c(0.5, 0.5), trim = 0)

col.cols <- rbind(groups = cutree(hc, 3))
pagoda.view.aspects(tamr2, cell.clustering = hc, box = TRUE, labCol = NA, margins = c(0.5, 20), col.cols = rbind(l2cols))


# compile a browsable app, showing top three clusters with the top color bar
app <- make.pagoda.app(tamr2, tam, varinfo, go.env, pwpca, clpca, col.cols = col.cols, cell.clustering = hc, title = "NPCs")
# show app in the browser (port 1468)
show.app(app, "pollen", browse = TRUE, port = 1468)  
saveRDS(app, file = "pollen.app.rds")

pagoda.show.pathways(c("GO:0022008 neurogenesis","GO:0048699 generation of neurons"), varinfo, go.env, cell.clustering = hc, margins = c(1,5), show.cell.dendrogram = TRUE, showRowLabels = TRUE, showPC = TRUE)

# get cell cycle signature and view the top genes
cc.pattern <- pagoda.show.pathways(c("GO:0000280 nuclear division", "GO:0007067 mitotic nuclear division"), varinfo, go.env, show.cell.dendrogram = TRUE, cell.clustering = hc, showRowLabels = TRUE)
# subtract the pattern
varinfo.cc <- pagoda.subtract.aspect(varinfo, cc.pattern)














