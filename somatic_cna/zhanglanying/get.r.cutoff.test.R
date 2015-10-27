get.specificity <-
function(chi, r, W, l, rho, overdisperse="no", phi=1, od.alpha=0) {
	1 - get.alpha(chi, r, W, l, rho, overdisperse, phi, od.alpha)
}

get.alpha <-
function(chi, r, W, l, rho, overdisperse="no", phi=1, od.alpha=0) {
	C = W*chi/l
	if (overdisperse == "no") {
		T = (r-1)*sqrt(C /(1+r^2))
	} else if (overdisperse %in% c("ql", "quasi-likelihood")) {
		T = (r-1)*sqrt(C /((1+r^2)*phi))
	} else if (overdisperse %in% c("nb", "negative binomial")) {
		T = (r-1)*sqrt(C /((1+r^2)*(1+od.alpha*C)))
	}

	alpha = if (rho > 1) { (1-pnorm(T)) } else { pnorm(T) }
	if (any(alpha > 1)) alpha[alpha>1] = 1
	if (any(alpha < 0)) alpha[alpha<0] = 0
	return(alpha)
}

get.sensitivity <-
function(chi, r, W, l, rho, overdisperse="no", phi=1, od.alpha=0) {
	get.power(chi, r, W, l, rho, overdisperse, phi, od.alpha)
}

get.power <-
function(chi, r, W, l, rho, overdisperse="no", phi=1, od.alpha=0) {
	C = W*chi/l
	if (overdisperse == "no") {
		T = (r-rho)*sqrt(C /(rho+r^2))
	} else if (overdisperse %in% c("ql", "quasi-likelihood")) {
		T = (r-rho)*sqrt(C /((rho+r^2)*phi))
	} else if (overdisperse %in% c("nb", "negative binomial")) {
		T = (r-rho)*sqrt(C /((rho+r^2)+(rho^2+r^2)*od.alpha*C))
	}

	power = if (rho > 1) { 1-pnorm(T) } else { pnorm(T) }	
	if (any(is.nan(power))) power[is.nan(power)] = 0
	return(power)
}

read.coverage.gatk <- 
function(file) {
	gatk = read.table(file, header=TRUE)
	chrpos = matrix(unlist(strsplit(as.character(gatk$Target),":")), ncol=2, byrow=TRUE)
	chr = factor(paste("chr",chrpos[,1],sep=""))
	pos = matrix(as.integer(unlist(strsplit(chrpos[,2],"-"))), ncol=2, byrow=TRUE)
	start = pos[,1]
	end = pos[,2]
	return(data.frame( probe=gatk$Target, 
			   chr=chr, 
			   probe_start=start, 
			   probe_end=end, 
			   targeted.base=end-start+1, 
			   sequenced.base=NA, 
			   coverage=as.numeric(gatk$total_coverage), 
			   average.coverage=as.numeric(gatk$average_coverage), 
			   base.with..10.coverage=NA))
}

get.alpha.inv <-
function(chi, W, l, rho, alpha, overdisperse="no", phi=1, od.alpha=0) {
	t = if (rho < 1) { qnorm(alpha) } else { qnorm(1-alpha) }
	if (overdisperse %in% c("ql", "quasi-likelihood")) { t = t*sqrt(phi) }
	c = W*chi/l
	# k = c/(c-t^2)
	# r = k + sqrt(k^2-1)
	oldWarn = getOption("warn")
	options(warn=-1)
	if (overdisperse %in% c("nb", "negative binomial")) {
		A = 1 + c*od.alpha
		r = (c + t*sqrt(2*c*A - A^2*t^2))/(c-A*t^2)
	} else {
		r = (c + t*sqrt(2*c-t^2))/(c-t^2) 
	}
	options(warn=oldWarn)
	if (any(r < 0 | is.nan(r))) { r[r<0] = NaN }
	return(r)
}

get.power.inv <-
function(chi, W, l, rho, power, overdisperse="no", phi=1, od.alpha=0) {
	t = if (rho < 1) { qnorm(power) } else { qnorm(1-power) }
	if (overdisperse %in% c("ql", "quasi-likelihood")) { t = t*sqrt(phi) }
	c = W*chi/l
	# k = rho*c/(c-t^2)
	# r = k + t*sqrt((k^2+k)/c)
	oldWarn = getOption("warn")
	options(warn=-1)
	if (overdisperse %in% c("nb", "negative binomial")) {
		A = rho + c*od.alpha*rho^2
		B = 1 + c*od.alpha
		r = (rho*c + t*sqrt(c*B*rho^2 - A*B*t^2 + c*A))/(c-B*t^2)
	} else {
		r = (rho*c + t*sqrt(rho*(c - t^2 + rho*c)))/(c-t^2) 
	}
	options(warn=oldWarn)
	if (any(r < 0 | is.nan(r))) { r[r<0] = NaN }
	return(r)
}

get.r.cutoff <-
function(chi, W, l, rho, min.spec, min.sens, option) {
	stopifnot(option %in% c("auc", "spec", "sens"))
	r.min.spec = get.alpha.inv(chi=chi, W=W, l=l, rho=rho, alpha=1-min.spec)
	# r.min.sens = get.power.inv(chi=chi, W=W, l=l, rho=rho, power=min.sens)
	power = min.sens
	t = if (rho < 1) { qnorm(power) } else { qnorm(1-power) }
	c = W*chi/l
	r.min.sens = (rho*c + t*sqrt(rho*(c - t^2 + rho*c)))/(c-t^2)
	opt.lower = if (rho < 1) { r.min.sens } else { r.min.spec }
	opt.upper = if (rho < 1) { r.min.spec } else { r.min.sens }
	#if (is.nan(opt.lower) || is.nan(opt.upper) || opt.lower > opt.upper) {
		# print("not enough power!")
		#return(c(cutoff=NaN, spec=NaN, sens=NaN))
	#}
	if (option == "auc") {
		old.warn = options()$warn
		options(warn=-1) # suppress warning
		# note optimize -auc because optim do minimization
		r.opt = optim(rho, function(this.r){-get.AUC(chi=chi, r=this.r, W=W, l=l, rho=rho)}, method="L-BFGS-B", lower=opt.lower, upper=opt.upper, control=list(abstol=1e-10))
		r.cutoff = r.opt$par
		options(warn=old.warn) # setback old value
	} else if (option == "spec") {
		r.cutoff = r.min.sens
	} else if (option == "sens") {
		r.cutoff = r.min.spec
	}
	return(c(cutoff=r.cutoff))
		 # spec=get.specificity(chi, r.cutoff, W, l, rho),
		 # sens=get.sensitivity(chi, r.cutoff, W, l, rho)))
}



normal_depth = "lu05.bed.HER2.normal.sample_interval_summary"
normal = read.coverage.gatk(normal_depth)
tumor_depth = "lu05.bed.HER2.tumor.sample_interval_summary"
tumor = read.coverage.gatk(tumor_depth)
n_sum = sum(as.numeric(normal$coverage), na.rm=TRUE)
t_sum = sum(as.numeric(tumor$coverage), na.rm=TRUE)
mean_sum = ( n_sum + t_sum ) / 2 
normal$coverage = ( normal$coverage / n_sum ) * mean_sum
normal$average.coverage = normal$coverage / normal$targeted.base

covered.exon = (normal$average.coverage > 0 & tumor$average.coverage > 0)
covered.exon[is.na(covered.exon)] = FALSE
eCNV = normal[,c("probe", "chr", "probe_start", "probe_end", "coverage", "average.coverage", "targeted.base")]

c = 0.3
test.num.copy=c(1,3)
rho = c + (1-c)*test.num.copy/2
min.spec=0.9999
min.sens=0.9999
option="spec"
l = read.len = 100
for (i in 1:nrow(normal)) {
		if (covered.exon[i]) {
			chi = normal$average.coverage[i]
			W = normal$targeted.base[i]
			r = matrix(NA, nrow=length(rho), ncol=3, dimnames=list(rho,c("cutoff","spec","sens")))
			for (j in 1:length(rho)) {
				r[j,] = get.r.cutoff(chi, W, l, rho[j], min.spec, min.sens, option)
			}
			r = data.frame(r)
			# cn = classify.logR(norm.log.ratio[i], log2(r$cutoff))
			# eCNV$copy.number[i] = cn
			eCNV$lower.cutoff[i] = r[1,1]
			eCNV$upper.cutoff[i] = r[2,1]
		} else {
			# eCNV$copy.number[i] = 0 # zero to signify "no coverage"
			eCNV$lower.cutoff[i] = NA
			eCNV$upper.cutoff[i] = NA
		}
}
write.table(eCNV, "lu05.bed.get.r.cutoff.txt",sep="\t", eol="\n", quote=F, col.names=F, row.names=F)
