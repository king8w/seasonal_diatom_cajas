#GAMs___________________________________________________________
## Perform GAMs across species (columns), predict, and plot 

models.gam<-function(taxa,env, title, gam.save = "gam.results.txt"){
  
  require(mgcv)
  
  #dev.new(width=4,height=6)
  par(mfrow = c(6, 6))
  par(mar = c(2.5, 3.5, 1, 0.5))
  par(mgp = c(1.5, 0.5, 0))
  par(oma = c(0, 0, 3, 0))
  
  #par(mfrow=c(length(taxa),1),mar=c(0,2,0,1),oma=c(4,4,2,1))
  
  for(i in 1:length(taxa)){
    spp <- taxa[,i]
    env_var<-env
    spp.data<-as.data.frame(cbind(spp,env_var))
    names(spp.data)<-c("spp","env_var")
    
    mod<-gam(spp~s(env_var,k=5), data=spp.data)

    pdat<-with(spp.data,data.frame(env_var=seq(min(env_var),max(env_var),length=200)))
    pred<-predict(mod,newdata=pdat)
    
    m.d<-Deriv(mod,n=200)
    
    pred1 <- with(spp.data, data.frame(env_var = seq(min(env_var), max(env_var), length = 200)))
    pred1 <- cbind(pred1, as.data.frame(predict(mod, pred1, se.fit = TRUE, unconditional = TRUE)))
    pred1 <- transform(pred1,Fitted = fit,Upper = fit + (2 * se.fit),Lower = fit - (2 * se.fit),env_var = env_var)
    
    
    CI<-confint(m.d,alpha=0.01)
    S<-signifD(pred,m.d$env_var$deriv,CI$env_var$upper,CI$env_var$lower,eval=0)
    
    plot(spp~env_var,data=spp.data,type="p",pch=16,ylab="", xlab="", las=1)

    points(env_var,spp,pch=16,col="gray50")
    lines(pred1$env_var,pred1$Upper,lty=3,lwd=2)
    lines(pred1$env_var,pred1$Lower,lty=3,lwd=2)
    lines(pred1$env_var,pred1$fit,lty=1,lwd=2)
   
    lines(S$incr~env_var,data=pdat,lwd=3,col="blue")
    lines(S$decr~env_var,data=pdat,lwd=3,col="red")	
    
    mtext(colnames(diat)[i], side = 3, line = 0.2, cex = 0.6)
    
    for(r in 1:ncol(spp.data)){
      gam.results <- matrix()
      gam.results <- summary(mod)["s.pv"][r]
    }
    

    # if(!is.null(gam.save)){write.table(gam.results, file = gam.save, sep = "\t", col.names = TRUE, row.names=FALSE, na="")}
    
    return(gam.results)
    
  }
  mtext(title, outer = TRUE,
        side = 3, cex = 1.2, line = 1)

  
}
##___________________________________________________________

Deriv <- function(mod, n = 200, eps = 1e-7, newdata, term) {
  if(inherits(mod, "gamm"))
    mod <- mod$gam
  m.terms <- attr(terms(mod), "term.labels")
  if(missing(newdata)) {
    newD <- sapply(model.frame(mod)[, m.terms, drop = FALSE],
                   function(x) seq(min(x), max(x), length = n))
    names(newD) <- m.terms
  } else {
    newD <- newdata
  }
  X0 <- predict(mod, data.frame(newD), type = "lpmatrix")
  newD <- newD + eps
  X1 <- predict(mod, data.frame(newD), type = "lpmatrix")
  Xp <- (X1 - X0) / eps
  Xp.r <- NROW(Xp)
  Xp.c <- NCOL(Xp)
  ## dims of bs
  bs.dims <- sapply(mod$smooth, "[[", "bs.dim") - 1
  ## number of smooth terms
  t.labs <- attr(mod$terms, "term.labels")
  ## match the term with the the terms in the model
  if(!missing(term)) {
    want <- grep(term, t.labs)
    if(!identical(length(want), length(term)))
      stop("One or more 'term's not found in model!")
    t.labs <- t.labs[want]
  }
  nt <- length(t.labs)
  ## list to hold the derivatives
  lD <- vector(mode = "list", length = nt)
  names(lD) <- t.labs
  for(i in seq_len(nt)) {
    Xi <- Xp * 0
    want <- grep(t.labs[i], colnames(X1))
    Xi[, want] <- Xp[, want]
    df <- Xi %*% coef(mod)
    df.sd <- rowSums(Xi %*% mod$Vp * Xi)^.5
    lD[[i]] <- list(deriv = df, se.deriv = df.sd)
  }
  class(lD) <- "Deriv"
  lD$gamModel <- mod
  lD$eps <- eps
  lD$eval <- newD - eps
  lD ##return
}

confint.Deriv <- function(object, term, alpha = 0.05, ...) {
  l <- length(object) - 3
  term.labs <- names(object[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else { ## how many attempts to get this right!?!?
    ##term <- match(term, term.labs)
    ##term <- term[match(term, term.labs)]
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  res <- vector(mode = "list", length = length(term))
  names(res) <- term
  residual.df <- df.residual(object$gamModel)
  tVal <- qt(1 - (alpha/2), residual.df)
  ##for(i in term.labs[term]) {
  for(i in term) {
    upr <- object[[i]]$deriv + tVal * object[[i]]$se.deriv
    lwr <- object[[i]]$deriv - tVal * object[[i]]$se.deriv
    res[[i]] <- list(upper = drop(upr), lower = drop(lwr))
  }
  res$alpha = alpha
  res
}

signifD <- function(x, d, upper, lower, eval = 0) {
  miss <- upper > eval & lower < eval
  incr <- decr <- x
  want <- d > eval
  incr[!want | miss] <- NA
  want <- d < eval
  decr[!want | miss] <- NA
  list(incr = incr, decr = decr)
}

plot.Deriv <- function(x, alpha = 0.05, polygon = TRUE,
                       sizer = FALSE, term,
                       eval = 0, lwd = 3,
                       col = "lightgrey", border = col,
                       ylab, xlab, main, ...) {
  l <- length(x) - 3
  ## get terms and check specified (if any) are in model
  term.labs <- names(x[seq_len(l)])
  if(missing(term)) {
    term <- term.labs
  } else {
    term <- term.labs[match(term, term.labs)]
  }
  if(any(miss <- is.na(term)))
    stop(paste("'term'", term[miss], "not a valid model term."))
  if(all(miss))
    stop("All terms in 'term' not found in model.")
  l <- sum(!miss)
  nplt <- n2mfrow(l)
  tVal <- qt(1 - (alpha/2), df.residual(x$gamModel))
  if(missing(ylab))
    ylab <- expression(italic(hat(f)*"'"*(x)))
  if(missing(xlab)) {
    xlab <- attr(terms(x$gamModel), "term.labels")
    names(xlab) <- xlab
  }
  if (missing(main)) {
    main <- term
    names(main) <- term
  }
  ## compute confidence interval
  CI <- confint(x, term = term)
  ## plots
  layout(matrix(seq_len(l), nrow = nplt[1], ncol = nplt[2]))
  for(i in term) {
    upr <- CI[[i]]$upper
    lwr <- CI[[i]]$lower
    ylim <- range(upr, lwr)
    plot(x$eval[,i], x[[i]]$deriv, type = "n",
         ylim = ylim, ylab = ylab, xlab = xlab[i], main = main[i], ...)
    if(isTRUE(polygon)) {
      polygon(c(x$eval[,i], rev(x$eval[,i])),
              c(upr, rev(lwr)), col = col, border = border)
    } else {
      lines(x$eval[,i], upr, lty = "dashed")
      lines(x$eval[,i], lwr, lty = "dashed")
    }
    abline(h = 0, ...)
    if(isTRUE(sizer)) {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 1)
      S <- signifD(x[[i]]$deriv, x[[i]]$deriv, upr, lwr,
                   eval = eval)
      lines(x$eval[,i], S$incr, lwd = lwd, col = "blue")
      lines(x$eval[,i], S$decr, lwd = lwd, col = "red")
    } else {
      lines(x$eval[,i], x[[i]]$deriv, lwd = 2)
    }
  }
  layout(1)
  invisible(x)
}
