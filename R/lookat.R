lookat<-function (u, Y) 
{
    if (class(u) != "criticor") 
        stop("ERROR: object not of class criticor\n")
    if (u$grType == "persp") 
        stop("ERROR: graphic type persp not convenient. Use image or contour in critic()\n")
    S<-u$S
    dinf<-u$dinf
    dsup<-u$dsup
    durinf<-u$durinf
    dursup<-u$dursup 
    grType<-u$grType
    device = dev.cur()
    existsModel = FALSE
    repeat {
        cat("\nlookat: function to inspect point correlations.\n\n")
        cat("Left-click on the point you want to illustrate\nclick outside the frame to close the function\n")
        pt = locator(1)
        if (grType == "filledcontour") {
            begin = as.integer(pt$y)
            dur = round(as.integer(pt$x)/0.7) + durinf
        }
        else {
            begin = as.integer(pt$y)
            dur = as.integer(pt$x)
        }
        outFrame = (dur < 0) | (begin < dinf) | (begin > dsup) | 
            (dur > dursup)
        if (outFrame) {
            cat("Exploration finished\n")
            (break)()
        }
        cat("\nbegin = ", begin, "dur = ", dur, "\n")
        i = begin - dinf + 1
        j = dur - durinf + 1
        dev.new(xpos=-25-7*100)
        plot(S[i, j, ], Y, pch = 20, main = paste("regression:\n begin = ", 
            begin, "\nspan = ", dur), xlab = "Sum")
        m = lm(Y ~ S[i, j, ])
        existsModel = TRUE
        abline(m)
        rho = cor(S[i, j, ], Y)
        rhos = cor(S[i, j, ], Y, method = "spearman")
        cat("correlation coefficient (Pearson): ", rho, "\n")
        cat("correlation coefficient (Spearman): ", rhos, "\n")
        cat("R-square: ", rho^2, "\n")
        print(anova(m))
        dev.set(device)
    }
    if (existsModel) 
        return(invisible(m))
    else return(invisible())
}
