library(shiny)
library(compute.es)


shinyServer(function(input, output) {
    
    options(warn=-1)


#==========================================================
# Standardized Mean Difference
#==========================================================

    sliderValues <- reactive ({
        n1 <- as.integer(input$nx)
        n2 <- as.integer(input$ny)
        
        data.frame(
            n = c(n1, n2),
            Mean = c(input$mx, input$my),
            SD = c(input$sdx, input$sdy),
            stringsAsFactors=FALSE)
    })
    # Show the values using an HTML table
    output$values <- renderTable({
        sliderValues()
    })
    
    
    difference <- reactive({
            nx <- input$nx
            mx <- input$mx
            sdx <- input$sdx
            ny <- input$ny
            my <- input$my
            sdy <- input$sdy
            
            if (input$varequal) {
                df <- nx+ny-2
                v <- ((nx-1)*sdx^2+(ny-1)*sdy^2)/df
                diff <- round((mx - my), 3)
                diff.std <- sqrt(v * (1/nx + 1/ny))
                diff.lower <- round(diff + diff.std * qt(0.05/2, df),3)
                diff.upper <- round(diff + diff.std * qt(0.05/2, df, lower.tail = FALSE),3)
            } else {
                stderrx <- sqrt(sdx^2/nx)
                stderry <- sqrt(sdy^2/ny)
                stderr <- sqrt(stderrx^2 + stderry^2)
                df <- round(stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny - 1)),3)
                tstat <- round(abs(mx - my)/stderr,3)
                diff <- round((mx - my), 3)
                cint <- qt(1 - 0.05/2, df)
                diff.lower <- round(((tstat - cint) * stderr),3)
                diff.upper <- round(((tstat + cint) * stderr),3)
            }
            
            cat("Mean of the differences [95% CI] =", diff, "[", diff.lower,",", diff.upper,"]", "\n")
    })
    output$difference.out <- renderPrint({
        difference()
    })


    es <- reactive({
        nx <- input$nx
        mx <- input$mx
        sdx <- input$sdx
        ny <- input$ny
        my <- input$my
        sdy <- input$sdy
    
        mes(mx, my, sdx, sdy, nx, ny)
    })
    output$es.out <- renderPrint({
        es()
    })
    
    
    ttest <- reactive({
        nx <- input$nx
        mx <- input$mx
        sdx <- input$sdx
        ny <- input$ny
        my <- input$my
        sdy <- input$sdy
        
        if (input$varequal) {
            df1 <- input$nx+input$ny-2
            v1 <- ((input$nx-1)*input$sdx^2+(input$ny-1)*input$sdy^2)/df1
            tstat1 <- round(abs(input$mx-input$my)/sqrt(v1*(1/input$nx+1/input$ny)),3)
            diff <- round((input$mx - input$my), 3)
            P1 <- 2 * pt(-abs(tstat1), df1)
        
            cat("Independent t-test (equal variances assumed)", "\n",
            " t =", tstat1, ",", "df =", df1, ",", "p-value =", P1, "\n")
        
        } else {

            stderrx <- sqrt(input$sdx^2/input$nx)
            stderry <- sqrt(input$sdy^2/input$ny)
            stderr <- sqrt(stderrx^2 + stderry^2)
            df2 <- round(stderr^4/(stderrx^4/(input$nx - 1) + stderry^4/(input$ny - 1)),3)
            tstat2 <- round(abs(input$mx - input$my)/stderr,3)
            P2 <- 2 * pt(-abs(tstat2), df2)
        
            cat("Welch's t-test (equal variances not assumed)", "\n",
            " t =", tstat2, ",", "df =", df2, ",", "p-value =", P2, "\n")
        }
    })
    output$ttest.out <- renderPrint({
        ttest()
    })
    
    
    vartest <- reactive({
        if (input$vartest) {
            nx <- input$nx
            sdx <- input$sdx
            vx <- sdx^2
            ny <- input$ny
            sdy <- input$sdy
            vy <- sdy^2
            
            if (vx > vy) {
                f <- vx/vy
                df1 <- nx-1
                df2 <- ny-1
            } else {
                f <- vy/vx
                df1 <- ny-1
                df2 <- nx-1
            }
            
            p <- 2*pf(f, df1, df2, lower.tail=FALSE)
            dfs <- c("num df"=df1, "denom df"=df2)
            
            cat(" Test for equality of variances", "\n",
                "  F =", f, ",", "num df =", df1, ",", "denom df =", df2, "\n",
                "  p-value = ", p, "\n"
                )
        
        } else {
            cat("Test for equality of variances will be displayed if the option is selected.")
        }
    })
    output$vartest.out <- renderPrint({
        vartest()
    })





#==========================================================
# Chi-Squared Statistic
#==========================================================

    chies.value <- reactive({
                    Chisq <- input$chi.v
                    n <- input$chi.n
                    rowcoln <- input$chi.rcn
                    
                    r <- sqrt(Chisq/(n*(rowcoln-1)))
                    var.r <- (1 - r^2)^2/(n - 1)
                    z <- 0.5 * log((1 + r)/(1 - r))
                    var.z <- 1/(n - 3)
                    df <- (n) - 2
                    j <- 1 - (3/(4 * df - 1))
                    var.z <- 1/(n - 3)
                    
                    level <- 95
                    alpha <- (100 - level)/100
                    crit <- qt(alpha/2, df, lower.tail = FALSE)
                    zval.z <- z/sqrt(var.z)
                    pval.z <- 2 * pt(abs(zval.z), df, lower.tail = FALSE)
                    lower.z <- z - crit * sqrt(var.z)
                    upper.z <- z + crit * sqrt(var.z)
                    pval.r <- pval.z
                    lower.r <- ztor(lower.z)
                    upper.r <- ztor(upper.z)
                    
                    cat("Correlation ES:", "\n", "\n", "r [", level, 
                    "%CI] =", round(r, 3), "[", round(lower.r, 3), 
                    ",", round(upper.r, 3), "]", "\n", " var(r) =", 
                    round(var.r, 3), "\n", " p-value(r) =", round(pval.r, 
                    3), "\n", "\n")
    })
    output$chies.out <- renderPrint({
        chies.value()
    })




#==========================================================
# R Session Info (tabPanel "Standardized Mean Difference")
#==========================================================
    info1 <- reactive({
        info1 <- paste("This analysis was conducted with ", strsplit(R.version$version.string, " \\(")[[1]][1], ".", sep = "")# バージョン情報
        info2 <- paste("It was executed on ", date(), ".", sep = "")# 実行日時
        cat(sprintf(info1), "\n")
        cat(sprintf(info2), "\n")
    })

    output$info1.out <- renderPrint({
        info1()
    })

#==========================================================
# R Session Info (tabPanel "Chi-Squared Statistic")
#==========================================================
    info2 <- reactive({
        info1 <- paste("This analysis was conducted with ", strsplit(R.version$version.string, " \\(")[[1]][1], ".", sep = "")# バージョン情報
        info2 <- paste("It was executed on ", date(), ".", sep = "")# 実行日時
        cat(sprintf(info1), "\n")
        cat(sprintf(info2), "\n")
    })

    output$info2.out <- renderPrint({
        info2()
    })



}) # shinyServer(function(input, output) {