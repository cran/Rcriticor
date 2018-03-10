critic<-function (t, Y, fac = NULL, dinf = 10, durinf = 2, dsup = 90, 
    dursup = 90, nperm = 0, nboot = 0, period = 365, dt = 1, 
    seriesName = "year", grType = "image", roll = FALSE, alpha = 0.05, 
    ps.print = FALSE) 
{
    Rep = if (is.null(fac)) 
        F
    else T
    m = dsup - dinf + 1
    p = dursup - durinf + 1
    z = matrix(0, m, p)
    n = length(Y)
    fac = if (Rep) 
        fac
    else as.factor(1:n)
    cat("Critical correlogram\n")
    cat("====================\n\n")
    cat("J.S.Pierre - 2015\n\n")
    cat("Data type: ", if (Rep) 
        "With replicates\n\n"
    else "Without replicates\n\n")
    cat("Name of the series : ", seriesName, "\n")
    cat("length of the basic series : ", period, "\n\n")
    cat("Lower beginning\t", dinf, "\n")
    cat("Upper beginning\t", dsup, "\n")
    cat("Minimal span\t", durinf, "\n")
    cat("Maximal span\t", dursup, "\n")
    cat("Matrix dimension : ", dim(z), "\n")
    S = array(0, dim = c(m, p, n))
    if (!Rep) 
        for (i in dinf:dsup) for (j in durinf:dursup) {
            begin = i - dinf + 1
            dur = j - durinf + 1
            for (k in 1:n) {
                start = i + period * (k - 1)
                stop = (i + j + period * (k - 1))
                S[begin, dur, k] = sum(t[start:stop])
            }
            z[begin, dur] = cor(Y, S[begin, dur, ])
        }
    else {
        num = as.numeric(fac)
        if (!all(sort(num) == num)) 
            stop("******* ERROR: years not in increasing order ********\n")
        num = num - min(num) + 1
        for (i in dinf:dsup) for (j in durinf:dursup) {
            for (k in 1:n) S[i - dinf + 1, j - durinf + 1, k] = sum(t[(i + 
                period * (num[k] - 1)):(i + j + period * (num[k] - 
                1))])
            z[i - dinf + 1, j - durinf + 1] = cor(Y, S[i - dinf + 
                1, j - durinf + 1, ])
        }
    }
    cat("\nMaximum correlation :", max(z), "\n")
    v = which(z == max(z), arr.ind = TRUE)
    debstar = v[1]
    durstar = v[2]
    intercept = debstar + dinf + durstar + durinf
    cat("Beginning of the period: ", debstar + dinf, "\n")
    cat("Span of the period: ", durstar + durinf, "\n")
    cat("Ordinary confidenceInterval:", cor.test(Y, S[debstar, 
        durstar, ])$conf.int, "\n")
    cat("p-value :", cor.test(Y, S[debstar, durstar, ])$p.value, 
        "\n")
    zB = z * 0
    for (k in 1:nperm) 
        {
        for (i in dinf:dsup) for (j in durinf:dursup) {
            YB = sample(Y)
            rho1 = cor(YB, S[i - dinf + 1, j - durinf + 1, ])
            rho2 = z[i - dinf + 1, j - durinf + 1]
            if (abs(rho1) < abs(rho2)) 
                zB[i - dinf + 1, j - durinf + 1] = zB[i - dinf + 
                  1, j - durinf + 1] + 1
        }
    #}
    if (nperm != 0) 
        zB = 1 - zB/nperm
    else zB = z * 0
    cat("\nNumber of random permutations: ", nperm, "\n")
    cat("Estimated p-value (random permutations):", if (nperm != 
        0) 
        zB[debstar, durstar]
    else "NA", "\n")
    if (nboot != 0) {
        zBoot = z * 0
        cat("\n\nBootstrap estimation\n")
        cat("====================\n\n")
        Sboot = S
        pseudo = NULL
        navalue = FALSE
        for (k in 1:nboot) {
            zBoot = z * 0
            fc = sample(fac, replace = T)
            YB = Y[as.numeric(fc)]
            Sboot = S[, , as.numeric(fc)]
            for (i in dinf:(dsup - 1)) for (j in durinf:(dursup - 
                1)) zBoot[i - dinf + 1, j - durinf + 1] = cor(YB, 
                Sboot[i - dinf + 1, j - durinf + 1, ])
            vb = which(zBoot == max(zBoot), arr.ind = TRUE)
            debstarb = vb[1] + dinf
            durstarb = vb[2] + durinf
            pseudo = rbind(pseudo, c(max(zBoot), debstarb, durstarb))
        }
        cat("\npseudovalues \n============\n\n")
        pseudo = data.frame(pseudo)
        names(pseudo) = c("rho", "begin", "duration")
        if (ps.print) 
            print(pseudo)
        rhoBoot = mean(pseudo[, 1], na.rm = TRUE)
        debBoot = mean(pseudo[, 2], na.rm = TRUE)
        durBoot = mean(pseudo[, 3], na.rm = TRUE)
        cat("\nBootstrap estimator for best correlation :", rhoBoot, 
            "\n")
        cat("Bootstrap estimator for  beginning of the period:", 
            debBoot, "\n")
        cat("Bootstrap estimator for  length of the period:", 
            durBoot, "\n")
        serho = sqrt(var(pseudo[, 1], na.rm = TRUE))
        cat("Standard error on rhomax :", serho, "\n")
        cat("\n Standard error for best correlation :", serho, 
            "\n")
        cat("Confidence interval for best correlation (alpha =", 
            alpha, "):\nLower : ", rhoBoot - qnorm(1 - alpha/2) * 
                serho, "\nUpper : ", rhoBoot + qnorm(1 - alpha/2) * 
                serho, "\n")
        cat("\nCovariance and correlation matrix of the two dimensions of the period :\n")
        cat("======================================================================\n\n")
        cat("Number of calculable pseudovalues :", dim(!is.na(pseudo))[1], 
            "\n")
        MVCV = var(cbind(pseudo$begin, pseudo$duration), na.rm = TRUE)
        MVCV = data.frame(MVCV)
        names(MVCV) = c("begin", "dur")
        row.names(MVCV) = c("begin", "dur")
        print(MVCV)
        cat("Correlation beginning-duration: \n")
        print(cov2cor(as.matrix(MVCV)))
        mu = c(debBoot, durBoot)
        A = if (nboot < 50) 
            qf(1 - alpha, 2, nboot - 2) * 2
        else qchisq(1 - alpha, 2)
        xmin = debBoot - sqrt(A * MVCV[1, 1])
        xmax = debBoot + sqrt(A * MVCV[1, 1])
        dx = (xmax - xmin)/500
        D = MVCV[1, 1] * MVCV[2, 2] - MVCV[1, 2]^2
        ellips = NULL
        for (x in c(seq(xmin, xmax, dx), xmax)) {
            X = x - mu[1]
            disc = X^2 * MVCV[1, 2]^2 - MVCV[1, 1] * (X^2 * MVCV[2, 
                2] - A * D)
            if (disc < 0) 
                disc = 0
            y1 = mu[2] + (X * MVCV[1, 2] - sqrt(disc))/MVCV[1, 
                1]
            y2 = mu[2] + (X * MVCV[1, 2] + sqrt(disc))/MVCV[1, 
                1]
            ellips = rbind(ellips, c(x, y1, y2))
        }
        ellips = data.frame(ellips)
        names(ellips) = c("x", "y1", "y2")
    }
    if (grType == "image") {
        dev.new()
        image(durinf:dursup, dinf:dsup, t(z), xlab = "span", 
            ylab = "start", main = "Correlogram", col = rainbow(21))
        abline(v = durstar + durinf)
        abline(debstar + dinf, 0)
        abline(intercept, -1, -1)
        points(durstar + durinf, debstar + dinf, pch = 19, bg = "red", 
            cex = 2)
        if (nboot != 0) {
            lines(ellips$x, ellips$y1)
            lines(ellips$x, ellips$y2)
            points(mu[1], mu[2], pch = 20)
        }
    }
    if (grType == "contour") {
        dev.new()
        contour(durinf:dursup, dinf:dsup, t(z), xlab = "span", 
            ylab = "start", main = "Correlogram", xaxs = "i", 
            yaxs = "i")
        abline(v = durstar + durinf)
        abline(debstar + dinf, 0)
        abline(intercept, -1, -1)
        points(durstar + durinf, debstar + dinf, pch = 19, bg = "red", 
            cex = 2)
        if (nboot != 0) {
            lines(ellips$x, ellips$y1)
            lines(ellips$x, ellips$y2)
            points(mu[1], mu[2], pch = 20, cex = 2)
            points(pseudo$duration, pseudo$begin, pch = 20, col = "red")
        }
    }
    if (grType == "filledcontour") {
        dev.new()
        filled.contour(durinf:dursup, dinf:dsup, t(z), xlab = "span", 
            ylab = "start", plot.axes = if (nboot != 0) {
                axis(1)
                axis(2)
                abline(v = durstar + durinf)
                abline(debstar + dinf, 0)
                abline(intercept, -1, -1)
                points(durstar + durinf, debstar + dinf, pch = 19, 
                  bg = "red", cex = 2)
            }
            else {
                axis(1)
                axis(2)
                abline(v = durstar + durinf)
                abline(debstar + dinf, 0)
                abline(intercept, -1, -1)
                points(durstar + durinf, debstar + dinf, pch = 19, 
                  bg = "red", cex = 2)
            }, main = "Correlogram",col=rainbow(21))
          }
    }

    if (grType == "persp") {
        dev.new()
        if (!roll) 
            imax = 15
        else imax = 360
        for (i in 15:imax) {
            u = persp(durinf:dursup, dinf:dsup, t(z), theta = i, 
                phi = 20, xlab = "span", ylab = "start", zlab = "rho", 
                shade = 1, col = "lightblue", ticktype = "detailed", 
                main = "Correlogram")
            points(trans3d(durstar + durinf, debstar + dinf, 
                min(z), u), pch = 19, bg = "red", cex = 2)
            v1 = trans3d(durstar + durinf, dinf, min(z), u)
            x1 = v1[1]
            y1 = v1[2]
            v2 = trans3d(durstar + durinf, dsup, min(z), u)
            x2 = v2[1]
            y2 = v2[2]
            lines(c(x1, x2), c(y1, y2))
            v1 = trans3d(durinf, debstar + dinf, min(z), u)
            x1 = v1[1]
            y1 = v1[2]
            v2 = trans3d(dursup + durinf, debstar + dinf, min(z), 
                u)
            x2 = v2[1]
            y2 = v2[2]
            lines(c(x1, x2), c(y1, y2))
            v1 = trans3d(durstar + durinf, debstar + dinf, min(z), 
                u)
            x1 = v1[1]
            y1 = v1[2]
            v2 = trans3d(durstar + durinf, debstar + dinf, max(z), 
                u)
            x2 = v2[1]
            y2 = v2[2]
            lines(c(x1, x2), c(y1, y2))
        }
    }
    if (nperm != 0) {
        par(mfrow = c(1, 2))
        dev.new()
        filled.contour(durinf:dursup, dinf:dsup, t(zB), nlevels = 4, 
            levels = c(1e-04, 0.001, 0.01, 0.05), xlab = "span", 
            ylab = "start", plot.axes = {
                axis(1)
                axis(2)
                abline(v = durstar + durinf)
                abline(debstar + dinf, 0)
                abline(intercept, -1, -1)
                points(durstar + durinf, debstar + dinf, pch = 19, 
                  bg = "red", cex = 2)
            }, main = "Significance map",col=rainbow(21))
        dev.new()
        filled.contour(durinf:dursup, dinf:dsup, 1 - t(zB), xlab = "span", 
            ylab = "start", plot.axes = {
                axis(1)
                axis(2)
                abline(v = durstar + durinf)
                abline(debstar + dinf, 0)
                abline(intercept, -1, -1)
                points(durstar + durinf, debstar + dinf, pch = 19, 
                  bg = "red", cex = 2)
            }, main = "Significance map",col=rainbow(21))
        par(mfrow = c(1, 1))
    }
    u = list()
    u$dinf = dinf
    u$dsup = dsup
    u$durinf = durinf
    u$dursup = dursup
    u$S = S
    u$z = z
    u$grType = grType
    class(u) = "criticor"
    invisible(u)
}
