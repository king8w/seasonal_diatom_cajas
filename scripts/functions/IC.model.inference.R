`IC.selection` <- function(IC.model, IC.sel = "AIC", C.hat = 1000, IC.max.var = NA, IC.fix.var = NULL, beta = FALSE, IC.save = "IC.model.results.txt", IC.coef = FALSE, ...) {

	#Delete results generated with script
	  unlink(IC.save, recursive = FALSE)
	  unlink("IC.coefficients.txt", recursive = FALSE)

	#Evaluate IC.sel, and if is different of AIC and is in list (AIC, AICc, QAIC, QAICc, BIC & DIC) 
	  if(IC.sel != "AIC" & IC.sel != "AICc" & IC.sel != "QAIC" & IC.sel != "QAICc" & IC.sel != "BIC") {
		stop(paste(IC.sel, " not recognised, possible types are ","\"AIC\", \"AICc\", \"QAIC\", \"QAICc\" and \"BIC\""))
	  } else {if(IC.sel == "AIC") IC.rank = 1 
	          if(IC.sel == "AICc") IC.rank = 2
	          if(IC.sel == "QAIC") IC.rank = 3
	          if(IC.sel == "QAICc") IC.rank = 4
	          if(IC.sel == "BIC") IC.rank = 5}


	#Define model intercept and extract all model terms
	  ICm.terms <- IC.get.model.terms(IC.model)
	  ICm.intercept <- attr(ICm.terms, "intercept")

	  ICm.n.variable <- length(ICm.terms)
	  ICm.table <- numeric(0)
	  ICm.formulas <- numeric(0)

	  if(IC.coef != FALSE)ICc.table <- numeric(0)


	  if (is.na(IC.max.var)) {IC.max.var <- ICm.n.variable} 
	 	 else {IC.max.var <- min(ICm.n.variable, IC.max.var)}
 
	  if (!is.null(IC.fix.var)) {
	 	if (inherits(IC.fix.var, "formula")) {
	     		if (IC.fix.var[[1]] != "~" || length(IC.fix.var) != 2) warning(sQuote("IC.fix.var"), " formula should be of form ", dQuote("~ a + b + c"))
			IC.fix.var <- c(IC.get.model.terms(IC.fix.var))}
		else if (!is.character(IC.fix.var)) {
			stop (sQuote("IC.fix.var"), " should be either a character vector with names of variables or a one-sided formula")}
		
	 	if (!all(IC.fix.var %in% ICm.terms)) {
			warning("Not all terms of ", sQuote("IC.fix.var"), " exist in ", sQuote("IC.model"))
			IC.fix.var <- IC.fix.var[IC.fix.var %in% ICm.terms]}
			IC.max.var <- IC.max.var - length(IC.fix.var)
	  }

	  if (IC.max.var > 0) {
		ICm.n.opt.vars <- 1:ICm.n.variable
			if (!is.null(IC.fix.var)) ICm.n.opt.vars <- ICm.n.opt.vars[!(ICm.terms %in% IC.fix.var)]

		ICm.comb <- lapply(seq(IC.max.var), combn, x = ICm.n.opt.vars)
		ICm.comb <- unlist(lapply(ICm.comb, function(.x) split(.x, col(.x))), recursive = FALSE)
		ICm.comb <- c(`0` = list(0), ICm.comb)} 
	  else {ICm.comb <- list (0)}

	  if (!is.null(IC.fix.var)) ICm.comb <- lapply(ICm.comb, append, (1:ICm.n.variable)[ICm.terms %in% IC.fix.var])

	  ICm.formulas <- lapply(ICm.comb, function(.x) reformulate(c("1", ICm.terms[.x]), response = "." ))
	  IC.supported <- sapply(ICm.formulas, Allowed.formula)
	  ICm.comb <- ICm.comb[IC.supported]
	  ICm.formulas <- ICm.formulas[IC.supported]
	  names(ICm.formulas) <- seq(ICm.formulas)

	  for(ICm in seq(length(ICm.comb))) {
		m.number <- ICm
		m.terms <- ICm.terms[ICm.comb[[ICm]]]
		m.formula <- ICm.formulas[[ICm]]
		m.row <- rep(NA, ICm.n.variable)
		m.row[match(m.terms, ICm.terms)] <- rep(1, length(m.terms))
		m.update.call <- call("update", substitute(IC.model), m.formula)

		m.run <- try(eval(m.update.call, parent.frame()))
			if (inherits(m.run, "try-error")) {ICm.formulas[[as.character(ICm)]] <- NA
			next;
			}

		m.coefficients <- c(na.omit(match(ICm.terms, names(Extract.m.coefficients(m.run)))))
		m.intercept <- if (attr(ICm.terms, "intercept")) Extract.m.coefficients(m.run)["(Intercept)"] else NA
		m.beta.coefficients <- if (beta) Beta.coefficients(m.run)[,3] else Extract.m.coefficients(m.run)
		m.coefficients <- m.beta.coefficients[m.coefficients]
		m.row[match(names(m.coefficients), ICm.terms)] <- m.coefficients

		if(beta == FALSE){m.beta <- as.numeric(0)} else {m.beta <- as.numeric(1)}


		m.N = as.numeric(length(resid(m.run)))
		m.K = as.numeric(attr(logLik(m.run), "df"))
		m.independent.var <- as.numeric(m.K - 2)
		
		m.residual.df <- as.numeric(m.run$df.residual)
		m.deviance <- as.numeric(deviance(m.run))
		m.deviance.p <- as.numeric((1 - pchisq(m.deviance, m.N - (m.K - 1))))
		m.RSS <- as.numeric((summary(m.run)$sigma ^ 2) * (m.N - (m.K - 2 + 1)))
		m.sigma <- as.numeric(summary(m.run)$sigma)
		m.variance <- as.numeric(summary(m.run)$sigma ^ 2)
		m.r <- as.numeric(sqrt(summary(m.run)$r.squared))
		m.r2 <- as.numeric(summary(m.run)$r.squared)
		m.r2.adj <- as.numeric(summary(m.run)$adj.r.squared)

			if(is.null(summary(m.run)$fstatistic[1])) {m.F.value <- NA} else {m.F.value <- as.numeric(summary(m.run)$fstatistic[1])}
			if(is.null(summary(m.run)$fstatistic[2])) {m.F.df1 <- NA} else {m.F.df1 <- as.numeric(summary(m.run)$fstatistic[2])}
			if(is.null(summary(m.run)$fstatistic[3])) {m.F.df2 <- NA} else {m.F.df2 <- as.numeric(summary(m.run)$fstatistic[3])}
			if(is.null(summary(m.run)$fstatistic[1])) {m.F.p <- NA} else {m.F.p <- 1 - as.numeric(pf(m.F.value, m.F.df1, m.F.df2))}

		m.AIC <- IC.AIC(m.run)
		m.AICc <- IC.AICc(m.run)
		m.QAIC <- IC.QAIC(m.run, C.hat)
		m.QAICc <- IC.QAICc(m.run, C.hat)
		m.BIC	<- IC.BIC(m.run)

		m.row <- c(m.number, m.intercept, m.row, m.beta, m.independent.var, m.N, m.K, (m.N / m.K), m.residual.df, m.deviance, m.deviance.p, m.RSS, m.sigma, 
			   m.variance, m.r, m.r2, m.r2.adj, m.F.value, m.F.df1, m.F.df2, m.F.p, m.AIC, m.AICc, m.QAIC, m.QAICc, m.BIC, IC.rank) 

		ICm.table <- rbind(ICm.table, m.row)

			if(IC.coef != FALSE){
				m.coeff <- cbind(m.number, Term = row.names(summary(m.run)$coefficients), summary(m.run)$coefficients, Beta.coefficients(m.run))
				ICc.table <- rbind(ICc.table, m.coeff)}	
	  }
	 	  
	  ICm.formulas[is.na(ICm.formulas)] <- NULL
	  ICm.table <- data.frame(ICm.table, row.names=1:NROW(ICm.table))
	  colnames (ICm.table) <- c("Model.number", "(Intercept)", ICm.terms, "Beta.coeff", "N.indep.var", "N.data", "K.par", "N/K", "Residual.df", "Deviance", "Deviance.pvalue", 
				    "RSS", "Sigma", "Variance", "r", "r2", "Adj.r2", "F.value", "F.df1", "F.df2", "F.pvalue", "AIC", "AICc", "QAIC", "QAICc", "BIC", "Sel.by")

	  o <- order(ICm.table[, IC.sel], decreasing = FALSE)
	  ICm.table <- ICm.table[o, ]
	  ICm.table$Delta <- as.numeric(ICm.table[, IC.sel] - min(ICm.table[, IC.sel]))
	  ICm.table$Weight <- as.numeric(exp(-ICm.table$Delta / 2) / sum(exp(-ICm.table$Delta / 2)))

	  class(ICm.table) = c("m.selection", "data.frame")
	  attr(ICm.table, "formulas") <- ICm.formulas[o]
	  attr(ICm.table, "global") <- IC.model
	  attr(ICm.table, "terms") <- c("(Intercept)", ICm.terms)

	  if(!is.null(IC.save)){write.table(ICm.table, file = IC.save, sep = "\t", col.names = TRUE, row.names=FALSE, na="")}


		if(IC.coef != FALSE){
			colnames (ICc.table) <- c("Model.number", "Term", "B", "SE", "t.value", "t.pvalue", "Stand.B", "Stand.B.SE", "Stand.B.t.value", "Stand.B.t.pvalue")
			write.table(ICc.table, file = "IC.coefficients.txt", sep = "\t", col.names = TRUE, row.names=FALSE, na="")}

	  return(ICm.table)
	 }

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.get.model.terms.default` <- function(x, ...) {return(IC.get.model.terms(as.formula(formula(x))))}

	`IC.get.model.terms.formula` <- function(x, ...) {
		My.model.terms <- terms(x)
		Return.terms <- attr(terms(x),"term.labels")
		if (length(Return.terms) > 0) {
			Return.terms <- Return.terms[order(Return.terms)]
			i <- grep(" ", Return.terms)
			Return.terms[i] <- paste("(", Return.terms[i] , ")")

			My.model.terms <- terms(as.formula(paste(". ~", paste(Return.terms, sep=" ", collapse=" + "))))
			Return.terms <- attr(My.model.terms, "term.labels")
		}

		attr(Return.terms, "intercept") <- attr(My.model.terms, "intercept")
		Return.terms
	}

`IC.get.model.terms` <- function (x, ...) UseMethod("IC.get.model.terms")

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Allowed.formula` <- function(m.formula) {
	Factors <- attr(terms(m.formula), "factors")
	return(all(Factors < 2))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Beta.coefficients` <- function(model) {

	summ <- summary(model)
	m.coef <- summ$coefficients[,1]
	m.se <- summ$coefficients[,2]

	response.sd <- sd(eval(attr(model$terms ,"variables"), envir=model$model)[[attr(model$terms ,"response")]])
	m.terms.sd <- sd(model.matrix(model))
	bx <- m.terms.sd / response.sd

	m.b.coef <- m.coef * bx
	m.b.se <- m.se * bx

	ret <- data.frame(m.coef, m.se, m.b.coef, m.b.se)

	colnames(ret) <- c("Estimate", "Std. Err.", "Beta", "Std. Err. Beta")
	rownames(ret) <- names(model$coefficients)

	ret <- as.matrix(ret)
	return (ret)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Extract.m.coefficients.default` <- function(model) {model$coefficients}
`Extract.m.coefficients` <- function (model) UseMethod("Extract.m.coefficients")
`Extract.m.coefficients.spautolm` <- function(model) { model$fit$coefficients}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`print.m.selection` <- function(x, ...) {
	x$Weight <- round(x$Weight, 3)
	nn <- attr(x, "terms")
	names(x)[seq(along=nn)] <- sapply( strsplit(nn, ":"), function(xx) paste(sapply(xx, abbreviate, 15 )  , collapse=":") )
	cat ("Model selection table", "\n")
	print(signif(as.matrix(x), digits=4), na.print="")
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Select.models` <- function(m.data, subset = cumsum(Weight) <= 0.95, ...) {

	 #subset = evalution formula, e.g. = cumsum(Weight) <= 0.95 or delta <= 4
	 o <- order(m.data[, "Weight"], decreasing = TRUE)
	 m.data <- m.data[o, ]

	 subset <- eval(substitute(subset), envir = m.data, enclos = parent.frame())
	 gmod <- attr(m.data, "global")
	 frm <- attr(m.data, "formulas")[subset]
	 sgmod <- substitute(gmod)
	 models <- lapply(frm, function(.x) eval(call("update", sgmod, .x), sys.parent(3)))
	 if (!is.null(attr(m.data, "rank.call"))) {attr(models, "rank.call") <- attr(m.data, "rank.call")}

	 return(models)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Average.models` <- function(m.selected, IC.sel = "AIC", C.hat = 1000, beta = FALSE, na.values = c("NA", "0"), alpha = 0.05, Avm.save = "Average.model.results.txt") {

	 unlink(Avm.save, recursive = FALSE)

	 na.values <- match.arg(na.values)
 	 Models <- m.selected
	 C.hat <- C.hat
	 if (length(Models) == 1) {stop("Only one model supplied. Thus nothing to do")}
	 m.selected <- Models[[1]] #Select the best model


	 if(IC.sel %in% c("AIC", "AICc", "QAIC", "QAICc", "BIC")) {
	 	cl <- as.call(c(as.name("sapply"), quote(Models), quote(paste("IC.", IC.sel, sep = "")), quote(C.hat)))
	 	M.rank <- eval(cl)
	 }else {stop(paste(IC.sel, " not recognised, possible types are ","\"AIC\", \"AICc\", \"QAIC\", \"QAICc\" and \"BIC\""))} 


	 Deviance <- sapply (Models, deviance)
	 Delta <- M.rank - min(M.rank)
	 Weight <- exp(-Delta / 2) / sum(exp(-Delta / 2))
	 M.order <- order(Weight, decreasing=TRUE)
	 M.rank <- M.rank[M.order]
	 Delta <- Delta[M.order]
	 Weight <- Weight[M.order]
	 Models <- Models[M.order]
	 Deviance <- Deviance[M.order]
	 M.table <- data.frame (Deviance, M.rank, Delta, Weight)
	 colnames(M.table) <- c("Deviance", IC.sel, "Delta", "Weight")

	 M.parameters <- unique(unlist(lapply(Models, function(m) names(Extract.m.coefficients(m)))))
	 M.terms <- unique(unlist(lapply(Models, IC.get.model.terms)))
	 M.terms <- M.terms[order(sapply(gregexpr(":", M.terms), function(x) if(x[1] == -1) 0 else length(x)), M.terms)]
	 M.parameters <- M.parameters[order(sapply(gregexpr(":", M.parameters), function(x) if(x[1] == -1) 0 else length(x)), M.parameters)]
	 M.coef <- M.var <- M.df <- numeric(0)
	 ac <- rep(0, length = length(M.parameters))

 	 for (m in Models) {
		m.tTable <- summary(m)$coefficients
		N <- length(resid(m))
		m.coef <- m.tTable[,1]
		m.var <- m.tTable[,2]
 		m.df <- N - length(m.coef)
		m.names <- paste(match(Extract.m.coefficients(m), M.terms), collapse="+")

		if (beta) {
			response.sd <- sd(model.frame(m.selected)[, attr(terms(m.selected) ,"response")])
			m.vars.sd <- sd(model.matrix(m))
			bx <- m.vars.sd / response.sd
			m.coef <- m.coef * bx
			m.var <- m.var * bx
		}
		m.vars <- match(M.parameters, rownames(m.tTable))
		M.coef <- rbind(M.coef, model = c(m.coef[m.vars]))
		M.var <- rbind(M.var, model = c(m.var[m.vars]))
		M.df <- append(M.df, m.df)
	 }

	 M.names <- sapply(Models, function(x) paste(match(IC.get.model.terms(x), M.terms), collapse="+"))
	 M.importance <- apply(Weight * t(sapply(Models, function(x) M.terms %in% IC.get.model.terms(x))), 2, sum)
	 names(M.importance) <- M.terms

	#check if models are unique:
	 dup <- duplicated(M.names)
	 	if (any(dup)) {dup <- table(M.names)
				   dup <- seq(M.names)[M.names %in% names(dup[dup > 1])]
				   stop("Models are not unique. Duplicates: ", paste(dup, collapse=", "))}

	 rownames(M.var) <- rownames(M.coef) <- rownames(M.table) <- M.names
	 if (na.values == "0") {all.coef[is.na(M.coef)] <- 0
					all.var[is.na(M.var)] <- 0}

	 M.average <- t(sapply(seq_along(M.parameters), function(i) Average.parameters(M.coef[,i], M.var[,i], M.df, Weight, alpha)))
	 M.coef[M.coef == 0] <- NA
	 M.var[M.var == 0] <- NA
 
	 M.importance <- sort(M.importance, decreasing=T)
	 colnames(M.coef) <- colnames(M.var) <- rownames(M.average) <-  M.parameters
       	 names(M.terms) <- seq_along(M.terms)
	
	 if (!is.null(IC.sel)) {colnames(M.table)[2] <- as.character(IC.sel)}

	 M.result <- list(summary = M.table,
			 Coefficients = M.coef,
			 Variable.codes = M.terms,
			 Variance = M.var,
			 Averaged.model = M.average,
			 Relative.importance = M.importance,
			 Weights = Weight,
			 Beta = beta,
			 Terms = M.parameters)

	 class(M.result) <- "Average"

		if(!is.null(Avm.save)){
			M.res <- c("Selected models summary:")
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = FALSE, row.names=FALSE, append = TRUE)
			M.res <- data.frame(M.result[1])
			M.res <- cbind(Selected.model = row.names(M.res), M.res)
			colnames(M.res) <- c("Regression.model", "Deviance", IC.sel, "Delta", "Weight")
			empty.row <- rep(NA, length(M.res))
			M.res <- rbind(M.res, empty.row)
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = TRUE, row.names=FALSE, append = TRUE)

			M.res <- c("Selected models - Variables codes:")
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = FALSE, row.names=FALSE, append = TRUE)
			M.res <- data.frame(M.result[3])
			M.res <- cbind(Number = row.names(M.res), M.res)
			colnames(M.res) <- c("Code", "Model.term")
			empty.row <- rep(NA, length(M.res))
			M.res <- rbind(M.res, empty.row)
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = TRUE, row.names=FALSE, append = TRUE)

			if(beta == FALSE) {M.res <- c("Averaged model parameters - Beta unstandarized:")}
			if(beta == TRUE) {M.res <- c("Averaged model parameters - Beta standarized:")}
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = FALSE, row.names=FALSE, append = TRUE)
			M.res <- data.frame(M.result[5])
			M.res <- cbind(Model.term = row.names(M.res), M.res)
			colnames(M.res) <- c("Model.term", "Coefficient", "Variance", "SE", "Unconditional SE", "Lower CI", "Upper CI")
			empty.row <- rep(NA, length(M.res))
			M.res <- rbind(M.res, empty.row)
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = TRUE, row.names=FALSE, append = TRUE)

			M.res <- c("Selected models - Relative variable importance:")
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = FALSE, row.names=FALSE, append = TRUE)
			M.res <- data.frame(M.result[6])
			M.res <- cbind(Model.term = row.names(M.res), M.res)
			colnames(M.res) <- c("Model.term", "Relative.importance")
			empty.row <- rep(NA, length(M.res))
			M.res <- rbind(M.res, empty.row)
			write.table(M.res, file = Avm.save, sep = "\t", na = "", col.names = TRUE, row.names=FALSE, append = TRUE)
		}
	 return(M.result)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.AIC` <- function(object, ...) {
		Link <- logLik(object)
		k <- attr(Link,"df")
		IC.aic <- as.numeric(-2*logLik(object)+2*k)

		return(IC.aic)	
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.AICc` <- function(object, ...) {
		Link <- logLik(object)
		k <- attr(Link,"df")
  		n <- length(resid(object))
		IC.aicc <- as.numeric(-2*logLik(object) + 2*k*n/max(n-k-1,0))

		return(IC.aicc)	
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.QAIC` <- function(object, C.hat, ...) {
		Link <- logLik(object)
		k <- attr(Link,"df")
		IC.qaic <- as.numeric(return(Link / C.hat + 2 * k))

		return(IC.qaic)	
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.QAICc` <- function(object, C.hat, ...) {
		Link <- logLik(object)
		k <- attr(Link,"df")
		n= length(resid(object))
		IC.qaicc <- as.numeric(Link / C.hat + 2*k*n/max(n-k-1,0))

		return(IC.qaicc)	
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`IC.BIC` <- function(object, ...) {
		Link <- logLik(object)
		k <- attr(Link,"df")
		n= length(resid(object))
		IC.bic <- as.numeric(-2*logLik(object) + k*log(n))

		return(IC.bic)	
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Average.parameters` <- function(x, se, npar, weight, alpha = 0.05) {
	
	 weight[is.na(weight)] <- 0 # not really necessary
	 weight.x <- weighted.mean(x, weight, na.rm = TRUE)
	 x.sqdiff <- (x - weight.x) ^ 2
	 variance.x <- se ^ 2
	 variance.coef <- weighted.mean(variance.x + x.sqdiff, weight, na.rm = TRUE) ^ 2
	 se.coef <- weighted.mean(sqrt(variance.x + x.sqdiff), weight, na.rm = TRUE)
	 z <- c((qt(1 - (alpha / 2), npar) / qnorm(1 - (alpha / 2))) ^ 2)
	 unconditional.se <- weighted.mean(sqrt((variance.x * z) + x.sqdiff), weight, na.rm = TRUE)
	 confidence.interval <- qnorm(1 - (alpha / 2)) * unconditional.se
	 return(c(`Coefficient` = weight.x, `Variance` = variance.coef,  `SE` = se.coef, `Unconditional SE` = unconditional.se, 
		 `Lower CI` = weight.x - confidence.interval, `Upper CI` = weight.x + confidence.interval))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`print.Average` <- function(x, ...) {
	 cat("\nModel summary:\n")
	 print(signif(x$summary,3))

	 cat("\nVariables:\n")
	 print(x$Variable.codes, quote= F)

	 if(x$Beta == FALSE){cat("\nAveraged model parameters - Beta unstandarized:\n")} else {cat("\nAveraged model parameters - Beta standarized:\n")}
	 print(signif(x$Averaged.model, 3))

	 cat("\nRelative variable importance:\n")
	 print(round(x$Relative.importance, 2))
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`VIF.evaluation` <- function(m.data, g.model, beta = FALSE, VIF.deletion = TRUE, VIF.value = 5, VIF.save = "VIF.selected.models.txt", VIF.analysis = "VIF.results.txt", ...) {

	unlink(VIF.save, recursive = FALSE)
	unlink(VIF.analysis, recursive = FALSE)

	M.class <- class(m.data)
	M.global <- attr(m.data, "global")
	M.terms <- attr(m.data, "terms")
	M.formulas <- attr(m.data, "formulas")

	m.terms <- IC.get.model.terms(g.model)
	m.intercept <- attr(m.terms, "intercept")
	n.variable <- length(m.terms)
	vif.terms <- paste("vif.", m.terms, sep = "")
	vif.sel <- paste("VIF>",VIF.value, sep="")

	m.table <- numeric(0)

	for(mod in 1:NROW(m.data)) {
		m.number <- m.data$Model.number[[mod]]
		m.formula <- M.formulas[[mod]]
		m.update.call <- call("update", substitute(g.model), m.formula)
		m.run <- try(eval(m.update.call, parent.frame()))

		m.coefficients <- c(na.omit(match(m.terms, names(Extract.m.coefficients(m.run)))))
		m.intercept <- if (attr(m.terms, "intercept")) Extract.m.coefficients(m.run)["(Intercept)"] else NA
		m.beta.coefficients <- if (beta) Beta.coefficients(m.run)[,3] else Extract.m.coefficients(m.run)

		m.coefficients <- m.beta.coefficients[m.coefficients]
		m.row <- rep(NA, n.variable)
		m.row[match(names(m.coefficients), m.terms)] <- m.coefficients

		if(beta == FALSE){m.beta <- as.numeric(0)} else {m.beta <- as.numeric(1)}


		if(length(m.coefficients) > 1){vif.coefficients <- vif(m.run)} else {vif.coefficients <- is.na(m.coefficients)} 
		vif.row <- rep(NA, n.variable)
		vif.row[match(names(vif.coefficients), m.terms)] <- vif.coefficients
		if(na.omit(any(vif.coefficients > VIF.value))) {vif.del <- 1} else {vif.del <- 0}

		m.row <- c(m.number, m.intercept, m.row, m.beta, vif.row, vif.del)
		m.table <- rbind(m.table, m.row)
	}

	colnames(m.table) <- c("Model.number", "(Intercept)", m.terms, "Beta.coeff", vif.terms, vif.sel)
	if(!is.null(VIF.analysis)){write.table(m.table, file = VIF.analysis, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}

	NAMES <- names(m.data)
	m.data <- cbind(m.data, m.table[ ,vif.sel])
	colnames(m.data) <- c(NAMES, vif.sel)

	if(VIF.deletion == TRUE){
		colnames(m.data) <- c(NAMES, "VIF")
		selection <- eval(substitute(VIF != 1), envir = m.data, enclos = parent.frame())
		m.data <- m.data[selection == TRUE, ]
		M.formulas <- M.formulas[selection]

		if(NROW(m.data) == 0){stop(paste("Upsss!!! Only one model remains, probably something happened"))}
		m.data <- m.data[,!names(m.data)%in% c("VIF")]

		if(m.data$Sel.by[1] == 1){IC.sel <- c("AIC")}
		if(m.data$Sel.by[1] == 2){IC.sel <- c("AICc")}
		if(m.data$Sel.by[1] == 3){IC.sel <- c("QAIC")}
		if(m.data$Sel.by[1] == 4){IC.sel <- c("QAICc")}
		if(m.data$Sel.by[1] == 5){IC.sel <- c("BIC")}

		o <- order(m.data[, IC.sel], decreasing = FALSE)
		m.data <- m.data[o, ]
		M.formulas <- M.formulas[o]
		m.data$Delta <- as.numeric(m.data[, IC.sel] - min(m.data[, IC.sel]))
		m.data$Weight <- as.numeric(exp(-m.data$Delta / 2) / sum(exp(-m.data$Delta / 2)))
	}

	class(m.data) = M.class
	attr(m.data, "global") <- M.global
	attr(m.data, "formulas") <- M.formulas
	attr(m.data, "terms") <- M.terms

	if(!is.null(VIF.save)){write.table(m.data, file = VIF.save, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}

	return(m.data)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Likelihood.test` <- function(m.data, g.model, LH.deletion = TRUE, LH.alpha = 0.05, LH.save = "Likelihood.selected.models.txt", LH.analysis = "Likelihood.results.txt", ...) {

	unlink(LH.save, recursive = FALSE)
	unlink(LH.analysis, recursive = FALSE)

	M.class <- class(m.data)
	M.global <- attr(m.data, "global")
	M.terms <- attr(m.data, "terms")
	M.formulas <- attr(m.data, "formulas")

	Null.formula <- c(". ~ 1")
	Null.update.call <- call("update", substitute(g.model), Null.formula)
	Null.run <- try(eval(Null.update.call, parent.frame()))
	Null.link <- logLik(Null.run)
	Null.df <- attr(Null.link,"df")

	m.table <- numeric(0)

	for(mod in 1:NROW(m.data)) {
		m.number <- m.data$Model.number[[mod]]
		m.formula <- M.formulas[[mod]]
		m.update.call <- call("update", substitute(g.model), m.formula)
		m.run <- try(eval(m.update.call, parent.frame()))

		m.link <- logLik(m.run)
		m.df <- attr(m.link,"df")
		m.chi <- as.numeric(-2 * (Null.link - m.link))
		m.chi.df <- as.numeric(m.df  - Null.df)
		if(m.chi != 0 & m.chi.df != 0) {m.chi.p <- as.numeric(pchisq(m.chi, m.chi.df, lower.tail = FALSE))} else m.chi.p <- NA
		if(!is.na(m.chi.p) & m.chi.p > LH.alpha) {LH.del <- 1} else {LH.del <- 0}

		m.row <- c(m.number, Null.link, m.link, m.chi, m.chi.df, m.chi.p, LH.del)
		m.table <- rbind(m.table, m.row)
	}

	
	colnames(m.table) <- c("Model.number", "Null.Likelihood", "Model.Likelihood", "Chi.square", "df", "p.value", paste("LH.test>",LH.alpha, sep=""))
	if(!is.null(LH.analysis)){write.table(m.table, file = LH.analysis, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}

	if(LH.deletion == TRUE) {colnames(m.table) <- c("Model.number", "Null.Likelihood", "Model.Likelihood", "Chi.square", "df", "p.value", "LH.test")}
	m.data <- cbind(m.data, m.table[ ,!colnames(m.table) %in% c("Model.number")])

	if(LH.deletion == TRUE){
		selection <- eval(substitute(LH.test != 1), envir = m.data, enclos = parent.frame())
		m.data <- m.data[selection == TRUE, ]
		M.formulas <- M.formulas[selection]

		if(NROW(m.data) == 0){stop(paste("Upsss!!! Only one model remains, probably something happened"))}
		m.data <- m.data[,!names(m.data)%in% c("LH.test")]

		if(m.data$Sel.by[1] == 1){IC.sel <- c("AIC")}
		if(m.data$Sel.by[1] == 2){IC.sel <- c("AICc")}
		if(m.data$Sel.by[1] == 3){IC.sel <- c("QAIC")}
		if(m.data$Sel.by[1] == 4){IC.sel <- c("QAICc")}
		if(m.data$Sel.by[1] == 5){IC.sel <- c("BIC")}

		o <- order(m.data[, IC.sel], decreasing = FALSE)
		m.data <- m.data[o, ]
		M.formulas <- M.formulas[o]
		m.data$Delta <- as.numeric(m.data[, IC.sel] - min(m.data[, IC.sel]))
		m.data$Weight <- as.numeric(exp(-m.data$Delta / 2) / sum(exp(-m.data$Delta / 2)))
	}

	class(m.data) = M.class
	attr(m.data, "global") <- M.global
	attr(m.data, "formulas") <- M.formulas
	attr(m.data, "terms") <- M.terms

	if(!is.null(LH.save)){write.table(m.data, file = LH.save, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}

	return(m.data)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Subset.models` <- function(m.data, subset = cumsum(Weight) <= 0.95, sub.save = "IC.models.subset.txt", ...) {

	M.class <- class(m.data)
	M.global <- attr(m.data, "global")
	M.terms <- attr(m.data, "terms")

	o <- order(m.data[, "Weight"], decreasing = TRUE)
	m.data <- m.data[o, ]
	
	selection <- eval(substitute(subset), envir = m.data, enclos = parent.frame())
	m.data <- m.data[selection == TRUE, ]
	M.formulas <- attr(m.data, "formulas")[selection]

	class(m.data) <- M.class
	attr(m.data, "global") <- M.global
	attr(m.data, "formulas") <- M.formulas
	attr(m.data, "terms") <- M.terms	 

	if(!is.null(sub.save)){write.table(m.data, file = sub.save, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}
	return(m.data)
}
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Extract.models` <- function(m.data, ...) {

	o <- order(m.data[, "Weight"], decreasing = TRUE)
	m.data <- m.data[o, ]

	gmod <- attr(m.data, "global")
	frm <- attr(m.data, "formulas")
	sgmod <- substitute(gmod)
	models <- lapply(frm, function(.x) eval(call("update", sgmod, .x), sys.parent(3)))
	if (!is.null(attr(m.data, "rank.call"))) {attr(models, "rank.call") <- attr(m.data, "rank.call")}
	return(models)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Bias.models` <- function(m.data, g.model, g.Beta = FALSE, Bias.save = "Models.parameters.bias.txt",...) {

	unlink(Bias.save, recursive = FALSE)

	IC.sel <- colnames(m.data$summary[2])

	if(g.Beta == FALSE) {
		Beta = 0
		IC.beta <- as.numeric(0)
		M.parameters <- Extract.m.coefficients(g.model)}
	if(g.Beta == TRUE) {
		Beta = 1
		IC.beta <- as.numeric(1)
		M.parameters <- Beta.coefficients(g.model)[,1]}

	B.table <- matrix(data = NA, nrow = length(M.parameters), ncol = 6, byrow = FALSE, dimnames = NULL)
	colnames(B.table) <- c("Parameter", "Beta.coeff", "Sel.by", "Model.coefficient", "IC.Estimate", "Bias")
	
	B.table[, "Parameter"] <- names(M.parameters)
	B.table[, "Beta.coeff"] <- rep(Beta, NROW(B.table))
	B.table[, "Sel.by"] <- rep(IC.sel, NROW(B.table))
	B.table[match(names(M.parameters), B.table[ , "Parameter"]), "Model.coefficient"] <- M.parameters
	B.table[match(row.names(m.data$Averaged.model), B.table[ , "Parameter"]), "IC.Estimate"] <- m.data$Averaged.model[, "Coefficient"]
	B.table[, "Bias"] <- (as.numeric(B.table[, "IC.Estimate"]) - as.numeric(B.table[, "Model.coefficient"])) / as.numeric(B.table[, "IC.Estimate"])

	if(!is.null(Bias.save)){write.table(B.table, file = Bias.save, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}

	print(B.table, na.print="", quote = FALSE)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Models.evaluation` <- function(Full.model, m.selected, m.averaged, IC.sel = "AIC", C.hat = 1000, na.values = c("NA", "0"), Eval.save = "Models.evaluation.results.txt") {

	#With estimating purpouses beta coefficients have no sense
	  unlink(Eval.save, recursive = FALSE)

	  Est.table <- numeric(0)
	  Est.table <- row.names(Full.model$model)
	  Est.table <- cbind(Est.table, Full.model$model, Full.model$residuals, Full.model$fitted.values)

	 na.values <- match.arg(na.values)
 	 Models <- m.selected
	 C.hat <- C.hat
	 if (length(Models) == 1) {stop("Only one model supplied. Thus nothing to do")}

	 Best.model <- Models[[1]] #Select the best model
	 Est.table <- cbind(Est.table, Best.model$residuals, Best.model$fitted.values)


	 if(IC.sel %in% c("AIC", "AICc", "QAIC", "QAICc", "BIC")) {
	 	cl <- as.call(c(as.name("sapply"), quote(Models), quote(paste("IC.", IC.sel, sep = "")), quote(C.hat)))
	 	M.rank <- eval(cl)
	 }else {stop(paste(IC.sel, " not recognised, possible types are ","\"AIC\", \"AICc\", \"QAIC\", \"QAICc\" and \"BIC\""))} 
	
	Delta <- M.rank - min(M.rank)
	Sum.Delta <- sum(exp(-Delta / 2))
	
		Avg.values <- numeric(0)
		Avg.residuals <- numeric(0)

 		for(m in 1:NROW(Models)){
			m.Delta <- M.rank[m] - min(M.rank)	
			m.Weight <- exp(-m.Delta / 2) / Sum.Delta
			m.fitted <- Models[[m]]$fitted.values * m.Weight
			m.residuals <- Models[[m]]$residuals * m.Weight
			Avg.values <- cbind(Avg.values, m.fitted)
			Avg.residuals <- cbind(Avg.residuals, m.residuals)
		}

		Avg.model <- matrix(data = NA, nrow = NROW(Avg.values), ncol = 2, byrow = FALSE, dimnames = NULL)

	 	for(r in 1:NROW(Avg.values)){
			Avg.model[r, 1] <- sum(Avg.residuals[r, ])
			Avg.model[r, 2] <- sum(Avg.values[r, ])
		}

	Est.table <- cbind(Est.table, Avg.model)

	#Estimate model average from averaged coefficients
	  Avg.mod <- m.averaged$Averaged.model[,1]
	  Avg.data <- Full.model$model[, (names(Full.model$model) %in% names(Avg.mod))]
	  Avg.est <- matrix(data = NA, nrow = NROW(Avg.data), ncol = length(Avg.mod), byrow = FALSE, dimnames = NULL)

	  for(r in 1:NROW(Avg.data)){
		for(m in 2:length(Avg.mod)){
			Avg.est[r, (m - 1)] <- Avg.mod[[m]] * Avg.data[r, as.numeric(match(names(Avg.mod[m]), names(Avg.data)))]
		}
		
		Avg.est[r, m] <- sum(Avg.est[r, 1:(m - 1)]) + Avg.mod[[1]]
		}	

	Est.table <- cbind(Est.table, Avg.est[, m])
	colnames(Est.table) <- c("Index", colnames(Full.model$model),"Global.mod.Residuals", "Global.mod.Fitted", "Best.mod.Residuals","Best.model.Fitted", 
			         "Avg.mod.Residuals", "Avg.mod.Fitted", "Avg.coeff.mod.fitted")
	
	if(!is.null(Eval.save)){write.table(Est.table, file = Eval.save, sep = "\t", col.names = TRUE, row.names = FALSE, na = "")}
	return(Est.table)
}
	
#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

`Evaluation.plot` <- function(X.var = "Xvar", Y.var = "Yvar", Database, m.data, Graphs = c(1, 2, 3, 4, 5), Output = "Plot.png", Title = "Title", X.Title = "Title", Y.Title = "Title", Legend = c(1, 2, 3, 4), Xv = 2, Yv = 0.15,   ... ) {

	#Graphs types: 1 = "Observed.data"; 2 = "Global.model"; 3 = "Best.model"; 4 = "Averaged model"; 5 = "Averaged coefficients model"
	 
	unlink(Output, recursive = FALSE)

	if(any(is.numeric(Graphs))) {Graphs <- Graphs} else {Graphs <- c(1, 2, 3, 4, 5)}
	if(any(!is.na(Graphs))) {Graphs <- Graphs} else {Graphs <- c(1, 2, 3, 4, 5)}
	if(is.null(Graphs)) {Graphs <- c(1, 2, 3, 4, 5)}
	if(!any(Graphs < 6)) {Graphs <- c(1, 2, 3, 4, 5)}

	Plot.table <- matrix(data = NA, nrow = NROW(Database), ncol = 6, byrow = FALSE, dimnames = NULL)
	colnames(Plot.table) <- c("X.data", "Y.data", "Global", "Best", "Averaged", "Avg.coeff")

		if(!is.na(as.numeric(match(X.var, names(Database))))){
			Plot.table[, "X.data"] <- Database[, X.var]
		} else {stop(paste(X.var, " not recognised or not found in the Database provided"))} 

		if(!is.na(as.numeric(match(Y.var, names(Database))))){
			Plot.table[, "Y.data"] <- Database[, Y.var]
		} else {stop(paste(Y.var, " not recognised or not found in the Database provided"))} 

	 	if(any(Graphs == 2)) {Plot.table[,  "Global"] <- m.data[, "Global.mod.Fitted"]}
	 	if(any(Graphs == 3)) {Plot.table[,  "Best"] <- m.data[, "Best.model.Fitted"]}
	 	if(any(Graphs == 4)) {Plot.table[,  "Averaged"] <- m.data[, "Avg.mod.Fitted"]}
	 	if(any(Graphs == 5)) {Plot.table[,  "Avg.coeff"] <- m.data[, "Avg.coeff.mod.fitted"]}

	X.data <- Plot.table[, "X.data"]
	Y.data<- Plot.table[, 2:6]

	Legend.text <- numeric(0)
	Legend.color <- numeric(0)
	Legend.circle <- numeric(0)

	#Remove extra top margin:
	  par(mar = c(4, 4, 2, 1))  # Trim margin around plot [bottom,left,top,right]
	  par(tcl = -0.5)  # Switch tick marks
	  par(mgp = c(2.5, 0.5, 0))  # Set other margins; default c(3,1,0) [ticks with title, ticks with ticks labels, ticks with graph margin line]
	  par(xaxs = "i", yaxs = "i")  # Extend axis limits by 4% ("i" does no extension)
	  par(family = "sans")

	#Compute graph axes margins:
	  X.inf <- round(min(X.data, na.rm = TRUE) / 0.1) * 0.1 
	  X.sup <- round((max(X.data, na.rm = TRUE) + 0.1) / 0.1) * 0.1
	  X.range <- (X.sup - X.inf) / Xv

	  Y.inf <- round(min(Y.data, na.rm = TRUE) / 0.1) * 0.1  
	  Y.sup <- round((max(Y.data, na.rm = TRUE) + 0.1) / 0.1) * 0.1
	  Y.range <- (Y.sup - Y.inf) / Yv

	#Define graph:
 	  matplot(X.data, rep(NA, NROW(X.data)), axes = T, type = "p", lty = 1, las = 1, lwd = 2, col = "grey",
		ylim = c(Y.inf, Y.sup), yaxp = c(Y.inf, Y.sup, Y.range),
		xlim = c(X.inf, X.sup), xaxp = c(X.inf, X.sup, X.range),
		xlab = X.Title, ylab = Y.Title, font.lab = 2, cex.lab = 1.5, main = Title, font.main = 2, cex.main = 2)

	#Define graph-series:
	  if(any(Graphs == 1)) {
		 points(X.data, Y.data[, 1], col = "black", bg = "white", cex =1.25, pch = 21, lwd=2)
		 lmfit1 <- lm(Y.data[, 1] ~ X.data)
		 abline(lmfit1, col = "black", lty = 4, lwd=2)
	 	 Legend.text <- c(Legend.text, paste("Observed: r = ", as.numeric(signif(sqrt(summary(lmfit1)$r.squared), 2)), sep = ""))
		 Legend.color <- c(Legend.color, "black")
		 Legend.circle <- c(Legend.circle, "white") }

	  if(any(Graphs == 2)) {
		 points(X.data, Y.data[, 2], col = "black", bg = "black", cex =1.25, pch = 21, lwd=2)
		 lmfit2 <- lm(Y.data[, 2] ~ X.data)
		 abline(lmfit2, col = "black", lty = 1, lwd=2)
 	 	 Legend.text <- c(Legend.text, paste("Global model: r = ", as.numeric(signif(sqrt(summary(lmfit2)$r.squared), 2)), sep = ""))
		 Legend.color <- c(Legend.color, "black")
		 Legend.circle <- c(Legend.circle, "black") }

	  if(any(Graphs == 3)) {
		 points(X.data, Y.data[, 3], col = "black", bg = gray(.85), cex =1.25, pch = 21, lwd=2)
		 lmfit3 <- lm(Y.data[, 3] ~ X.data)
		 abline(lmfit3, col = gray(.85), lty = 1, lwd=2) 
	 	 Legend.text <- c(Legend.text, paste("Best model: r = ", as.numeric(signif(sqrt(summary(lmfit3)$r.squared), 2)), sep = ""))
		 Legend.color <- c(Legend.color, "black")
		 Legend.circle <- c(Legend.circle, gray(.85)) }

	  if(any(Graphs == 4)) {
		 points(X.data, Y.data[, 4], col = "black", bg = gray(.55), cex =1.25, pch = 21, lwd=2)
		 lmfit4 <- lm(Y.data[, 4] ~ X.data)
		 abline(lmfit4, col = gray(.55), lty = 1, lwd=2) 
	 	 Legend.text <- c(Legend.text, paste("Averaged model: r = ", as.numeric(signif(sqrt(summary(lmfit4)$r.squared), 2)), sep = ""))
		 Legend.color <- c(Legend.color, "black")
		 Legend.circle <- c(Legend.circle, gray(.55)) }

	  if(any(Graphs == 5)) {
		 points(X.data, Y.data[, 5], col = "black", bg = gray(.40), cex =1.25, pch = 21, lwd=2)
		 lmfit5 <- lm(Y.data[, 5] ~ X.data)
		 abline(lmfit5, col = gray(.40), lty = 2, lwd=2) 
	 	 Legend.text <- c(Legend.text, paste("Avg. coef. model: r = ", as.numeric(signif(sqrt(summary(lmfit5)$r.squared), 2)), sep = ""))
		 Legend.color <- c(Legend.color, "black")
		 Legend.circle <- c(Legend.circle, gray(.40)) }

	#Define legend.position:
	 if(length(Legend) != 1) {Legend <- c(1)}
	 if(!is.numeric(Legend)) {Legend <- c(1)}
	 if(is.na(Legend)) {Legend <- c(1)}
	 if(is.null(Legend)) {Legend <- c(1)}
	 if(!any(Legend < 5)) {Legend <- c(1)}

	  if(Legend == 1){Position = c("topleft")}
	  if(Legend == 2){Position = c("topright")}
	  if(Legend == 3){Position = c("bottomright")}
	  if(Legend == 4){Position = c("bottomleft")}
	  legend(Position, legend = Legend.text, cex=1.25, col = Legend.color, pt.bg = Legend.circle, pch=21)

	#Draw and export graph:
	  dev.print(device=png, Output, height = 480*5, width = 480*5, res = 72*5)
}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#`New.function` <- function(, ...) {}

#----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
