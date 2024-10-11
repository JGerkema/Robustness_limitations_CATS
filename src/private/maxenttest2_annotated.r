maxent.test2 <- function (model, obs, sub.c, nperm = 999, quick = TRUE, alpha = 0.05, 
    plot = TRUE, minperms = 20) {
  
   set.seed(19970606)
  
    if (is.vector(obs)) {
        s.names <- names(obs)
        obs <- matrix(obs, 1, length(obs))
        dimnames(obs) <- list("set1", s.names)
    }
    if (is.data.frame(obs)) 
        obs <- as.matrix(obs)
    obs.names <- dimnames(obs)
    if (!is.numeric(obs)) 
        stop("obs must only contain numeric values\n")
    if (dim(obs)[2] == 1 && dim(obs)[1] > 1)  
        obs <- t(obs) # If the matrix has 1 column and multiple rows, it is flipped
    if (any(is.na(obs))) 
        stop("no NA's allowed\n")
    if (any(obs < 0 | obs > 1)) 
        stop("obs must contain probabilities between 0 and 1\n")
    obs <- t(apply(obs, 1, function(x) x/sum(x))) 
    # This function takes the abundance of each species and divides it by the total 
    # abundance. When the relative abundances are used as input, the sum is 1, so nothing
    # changes. 
    
    dimnames(obs) <- obs.names
    states <- model$states
    if (sum(states <= 0) > 0) 
        stop("trait values must not be negative or zero")
    
 
    
    s.names <- dimnames(states)[[2]] # Extracts the species names
    constr <- model$constr # Extracts the cwm's
    prob <- model$prob # Extracts the predicted probabilities
    prior <- model$prior
    n.states <- dim(states)[2] # Extracts the number of species
    KLR2.fit <- function(oi, pi, n.states) {
        sel <- oi > 0
        qi.uniform <- rep(1/n.states, length(oi))
        1 - sum(oi[sel] * log(oi[sel]/pi[sel]))/sum(oi[sel] * 
            log(oi[sel]/qi.uniform[sel]))
    }
    sel <- obs > 0
    oi <- as.double(as.vector(as.matrix(obs[sel]))) #Takes only the abundances higher than 0
    pi <- as.double(as.vector(as.matrix(prob[sel])))
    KLR2.prior.plus.traits <- KLR2.fit(oi, pi, n.states)
   
    if (is.vector(constr)) {
        means.names <- names(constr)
        constr <- matrix(constr, 1, length(constr))
        dimnames(constr) <- list("set1", means.names)
        prior <- matrix(prior, 1, length(prior))
        dimnames(prior) <- list("set1", s.names)
        prob <- matrix(prob, 1, length(prob))
        dimnames(prob) <- list("set1", s.names)
    }
    n.sets <- dim(constr)[1]
    if (n.sets != dim(obs)[1]) 
        stop("number of rows in obs and constr should be equal\n")
    stat <- function(o, p, q) sum(o * log(p/q))
    values <- rep(NA, nperm)
    KLR2.null <- rep(NA, nperm)
    oi <- as.double(as.vector(as.matrix(obs)))
    qi <- as.double(as.vector(as.matrix(prior)))
    sel.oi <- (oi > 0)
    if (missing(sub.c)) {
        obs.stat <- stat(obs, prob, prior)
        count <- 0
        outside.ci <- FALSE
        while ((!outside.ci && count < nperm) | count < minperms) {
# The loop will continue until one of the following conditions is met:
# The computed statistic falls outside the confidence interval and the number of iterations exceeds 
# or equals nperm. 
# The number of iterations exceeds or equals minperms.
         
            count <- count + 1
            prob.temp <- matrix(NA, n.sets, n.states)
            for (j in 1:n.sets) {
                shuffled <- sample(1:n.states, n.states)
                states.perm <- states[, shuffled, drop = F]
                colnames(states.perm) <- s.names
                constr.perm <- functcomp(t(states.perm), obs[j, 
                  , drop = F])
                prob.temp[j, ] <- maxent2(constr.perm, states.perm, 
                  prior[j, ])$prob
            }
            pi <- as.double(as.vector(as.matrix(prob.temp)))
            KLR2.null[count] <- KLR2.fit(oi, pi, n.states)
            values[count] <- stat(obs, prob.temp, prior)
            val.temp <- values[1:count]
            p.temp <- (length(val.temp[val.temp >= obs.stat]) + 
                1)/(length(val.temp) + 1)
            if (quick && p.temp < 1) {
                ci.hi <- p.temp + 1.96 * sqrt((p.temp * (1 - 
                  p.temp))/count)
                ci.lo <- p.temp - 1.96 * sqrt((p.temp * (1 - 
                  p.temp))/count)
                outside.ci <- ci.hi <= alpha || ci.lo >= alpha
                if (outside.ci && count > minperms) 
                  nperm = count
            }
        }
    }
    else {
        if (length(sub.c) >= n.states) 
            stop("sub.c contains as many or more elements than there are states\n")
        if (is.character(sub.c) && !all(sub.c %in% dimnames(states)[[1]])) 
            stop("sub.c does not match the names of the state attributes\n")
        if (is.character(sub.c) && !all(sub.c %in% dimnames(constr)[[2]])) 
            stop("sub.c does not match the constraint names\n")
        if (is.character(sub.c)) 
            sub.c <- which(dimnames(states)[[1]] %in% sub.c)
        if (!is.vector(sub.c)) 
            stop("sub.c must be a vector\n")
        prob.a <- maxent(constr[, -sub.c, drop = F], states[-sub.c, 
            , drop = F], prior)$prob
        obs.stat <- stat(obs, prob, prob.a)
        count <- 0
        outside.ci <- FALSE
        while ((!outside.ci && count < nperm) | count < minperms) {
            count <- count + 1
            prob.temp <- matrix(NA, n.sets, n.states)
            for (j in 1:n.sets) {
                shuffled <- sample(1:n.states, n.states)
                states.perm <- states
                states.perm[sub.c, ] <- states[sub.c, shuffled, 
                  drop = F]
                colnames(states.perm) <- s.names
                constr.perm <- functcomp(t(states.perm), obs[j, 
                  , drop = F])
                prob.temp[j, ] <- maxent2(constr.perm, states.perm, 
                  prior[j, ])$prob
            }
            pi <- as.double(as.vector(as.matrix(prob.temp)))
            KLR2.null[count] <- KLR2.fit(oi, pi, n.states)
            values[count] <- stat(obs, prob.temp, prob.a)
            val.temp <- values[1:count]
            p.temp <- (length(val.temp[val.temp >= obs.stat]) + 
                1)/(length(val.temp) + 1)
            if (quick && p.temp < 1) {
                ci.hi <- p.temp + 1.96 * sqrt((p.temp * (1 - 
                  p.temp))/count)
                ci.lo <- p.temp - 1.96 * sqrt((p.temp * (1 - 
                  p.temp))/count)
                outside.ci <- ci.hi <= alpha || ci.lo >= alpha
                if (outside.ci && count > minperms) 
                  nperm = count
            }
        }
    }
    values <- values[!is.na(values)]
    p.perm <- (length(values[values >= obs.stat]) + 1)/(length(values) + 
        1)
    p.perm.hi <- p.perm + 1.96 * sqrt((p.perm * (1 - p.perm))/nperm)
    p.perm.lo <- p.perm - 1.96 * sqrt((p.perm * (1 - p.perm))/nperm)
    opqfit <- function(o, p, q) 1 - (sum((o - p)^2)/sum((o - 
        q)^2))
    r2.op <- cor(as.double(prob), as.double(obs))^2
    if (length(unique(as.double(prior))) != 1) 
        r2.oq <- cor(as.double(prior), as.double(obs))^2
    else r2.oq = 0
    fit <- opqfit(obs, prob, prior)
    if (!missing(sub.c)) {
        r2.opa <- cor(as.double(prob.a), as.double(obs))^2
        fit.a <- opqfit(obs, prob.a, prior)
    }
    mean.KLR2.null <- mean(KLR2.null[!is.na(KLR2.null)])
    if (plot) {
        if (missing(sub.c)) {
            par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
            fit.text <- bquote(fit[bold(o * "," * p) * "|" * 
                bold(q)] == .(round(fit, 3)))
            r2.op.text <- bquote(italic(KLR)^2 == .(round(KLR2.prior.plus.traits, 
                3)))
            r2.oq.text <- bquote(italic(KLR)^2 == .(round(mean.KLR2.null, 
                3)))
            plot(as.double(obs), as.double(prob), xlim = c(0, 
                1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities", 
                main = expression(bold(p)(bold(C) * "," * ~bold(q))))
            abline(0, 1, col = "grey25")
            text(0.1, 0.9, fit.text, cex = 1.2, pos = 4)
            text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4)
            lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", 
                col = "grey50")
            lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", 
                col = "grey50")
            mtext("arithmetic scale", line = 0)
            plot(as.double(obs) + 1e-04, as.double(prob) + 1e-04, 
                xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed probabilities", 
                ylab = "predicted probabilities", log = "xy", 
                main = expression(bold(p)(bold(C) * "," * ~bold(q))))
            abline(0, 1, col = "grey25")
            lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04, 
                1 + 1e-04), lty = "dashed", col = "grey50")
            lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 + 
                1e-04), lty = "dashed", col = "grey50")
            mtext("log10 scale, + 1e-4", line = 0)
            plot(as.double(obs), as.double(prior), xlim = c(0, 
                1), ylim = c(0, 1), xlab = "", ylab = "prior", 
                main = expression(bold(q)))
            abline(0, 1, col = "grey25")
            text(0.1, 0.9, r2.oq.text, cex = 1.2, pos = 4)
            lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", 
                col = "grey50")
            lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", 
                col = "grey50")
            mtext("arithmetic scale", line = 0)
            plot(as.double(obs) + 1e-04, as.double(prior) + 1e-04, 
                xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed probabilities", 
                ylab = "prior", log = "xy", main = expression(bold(q)))
            abline(0, 1, col = "grey25")
            lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04, 
                1 + 1e-04), lty = "dashed", col = "grey50")
            lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 + 
                1e-04), lty = "dashed", col = "grey50")
            mtext("log10 scale, + 1e-4", line = 0)
            mtext(expression(H[0] %->% bold(p)(bold(C) * "," * 
                ~bold(q)) == bold(q)), cex = 1.2, outer = T, 
                las = 1, line = 1)
            iter.text <- bquote(italic(Nperms) == .(round(count, 
                0)))
            p.text <- bquote(italic(P) == .(round(p.perm, 3)))
            mtext(p.text, cex = 1.2, outer = T, las = 1, line = -0.5)
            mtext(iter.text, cex = 0.75, outer = T, las = 1, 
                line = -1.5)
        }
        else {
            par(las = 1, mfcol = c(2, 2), oma = c(0, 0, 3, 0))
            fit.text <- bquote(fit[bold(o * "," * p) * "|" * 
                bold(q)] == .(round(fit, 3)))
            fit.a.text <- bquote(fit[bold(o * "," * p) * "|" * 
                bold(q)] == .(round(fit.a, 3)))
            r2.op.text <- bquote(italic(KLR)^2 == .(round(KLR2.prior.plus.traits, 
                3)))
            r2.opa.text <- bquote(italic(KLR)^2 == .(round(KLR2.null, 
                3)))
            plot(as.double(obs), as.double(prob), xlim = c(0, 
                1), ylim = c(0, 1), xlab = "", ylab = "predicted probabilities", 
                main = expression(bold(p)(bold(A) ~ union(bold(B) * 
                  "," * ~bold(q)))))
            abline(0, 1, col = "grey25")
            text(0.1, 0.9, fit.text, cex = 1.2, pos = 4)
            text(0.1, 0.75, r2.op.text, cex = 1.2, pos = 4)
            lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", 
                col = "grey50")
            lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", 
                col = "grey50")
            mtext("arithmetic scale", line = 0)
            plot(as.double(obs) + 1e-04, as.double(prob) + 1e-04, 
                xlim = c(1e-04, 1), ylim = c(1e-04, 1), xlab = "observed probabilities", 
                ylab = "predicted probabilities", log = "xy", 
                main = expression(bold(p)(bold(A) ~ union(bold(B) * 
                  "," * ~bold(q)))))
            abline(0, 1, col = "grey25")
            lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04, 
                1 + 1e-04), lty = "dashed", col = "grey50")
            lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 + 
                1e-04), lty = "dashed", col = "grey50")
            mtext("log10 scale, + 1e-4", line = 0)
            plot(as.double(obs), as.double(prob.a), xlim = c(0, 
                1), ylim = c(0, 1), xlab = "", ylab = "", main = expression(bold(p)(bold(A) * 
                "," * ~bold(q))))
            abline(0, 1, col = "grey25")
            text(0.1, 0.9, fit.a.text, cex = 1.2, pos = 4)
            text(0.1, 0.75, r2.opa.text, cex = 1.2, pos = 4)
            lines(x = c(0.05, 0.05), c(0, 1), lty = "dashed", 
                col = "grey50")
            lines(x = c(0, 1), c(0.05, 0.05), lty = "dashed", 
                col = "grey50")
            mtext("arithmetic scale", line = 0)
            plot(as.double(obs) + 1e-04, as.double(prob.a) + 
                1e-04, xlim = c(1e-04, 1), ylim = c(1e-04, 1), 
                xlab = "observed probabilities", ylab = "", log = "xy", 
                main = expression(bold(p)(bold(A) * "," * ~bold(q))))
            abline(0, 1, col = "grey25")
            lines(x = c(0.05 + 1e-04, 0.05 + 1e-04), c(0 + 1e-04, 
                1 + 1e-04), lty = "dashed", col = "grey50")
            lines(x = c(1e-04, 1 + 1e-04), c(0.05 + 1e-04, 0.05 + 
                1e-04), lty = "dashed", col = "grey50")
            mtext("log10 scale, + 1e-4", line = 0)
            mtext(expression(H[0] %->% bold(p)(bold(A) ~ union(bold(B)) * 
                "," * ~bold(q)) == bold(p)(bold(A) * "," * ~bold(q))), 
                cex = 1.2, outer = T, las = 1, line = 1)
            iter.text <- bquote(italic(Nperms) == .(round(count, 
                0)))
            p.text <- bquote(italic(P) == .(round(p.perm, 3)))
            mtext(p.text, cex = 1.2, outer = T, las = 1, line = -0.5)
            mtext(iter.text, cex = 0.75, outer = T, las = 1, 
                line = -1.5)
        }
    }
    res <- list()
    res$fit <- fit
    if (!missing(sub.c)) {
        res$fit.a <- fit.a
        res$r2.a <- r2.opa
    }
    res$obs.stat <- obs.stat
    res$nperm <- nperm
    res$pval <- p.perm
    res$ci.pval <- c(p.perm.lo, p.perm.hi)
    res$KLR2.null <- KLR2.null[!is.na(KLR2.null)]
    res$mean.KLR2.null <- mean(KLR2.null[!is.na(KLR2.null)])
    res$KLR2.prior.plus.traits <- KLR2.prior.plus.traits
    res$values <- values
    return(res)
}
