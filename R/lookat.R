function(S, dinf = dinf, dsup = dsup, durinf = durinf, 
        dursup = dursup) 
{
        repeat {
            pt = locator(1)
            begin = as.integer(pt$y)
            dur = as.integer(pt$x)
            cat("\nlookat: function to inspect point correlations.\n\n")
            print(begin)
            print(dur)
            i = begin - dinf + 1
            j = dur - durinf + 1
            plot(S[i, j, ], Y, pch = 20)
            abline(lm(Y ~ S[i, j, ]))
            cat("correlation coefficient: ", cor(S[i, j, ], Y), 
                "\n")
            if (dur < 0) 
                stop("finished")
        }
}
