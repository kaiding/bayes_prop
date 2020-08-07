library(shiny)
library(shinybusy)
library(markdown)
library(mosaic)
library(flextable)
library(officer)

# Define UI for dataset viewer app ----
ui <- navbarPage(title = "Simple R ShinyApp for binary endpoints using Bayesian Posterior Probability",
                 
                 tabsetPanel(
                              tabPanel("One-Arm",
                                       sidebarLayout(
                                           sidebarPanel(
                                               p("Beta prior (in blue)", style="color: blue"),
                                               p("Binomial likelihood (in green)", style="color: green"),
                                               p("Beta posterior (in red)", style="color: red"),
                                               numericInput("alphaG", "alpha, shape parameter for Beta prior:", 
                                                            value = 1,min = 0.01, max = 100),
                                               numericInput("betaG", "beta, shape parameter for Beta prior:", 
                                                            value = 1,min = 0.01, max = 100),
                                               numericInput("rG", "Number of observed responses:", 15,
                                                            min = 1, max = 100),
                                               numericInput("nG", "Number of participants enrolled:", 20,
                                                            min = 1, max = 100),
                                               numericInput("pG", "Target response rate:", 0.5,
                                                            min = 0, max = 1),
                                               actionButton(inputId = "plotG",
                                                            label = "Draw")
                                           ),
                                           
                                           mainPanel(plotOutput("plotG"),
                                                     
                                                     textOutput("description_G"))
                                           )
                                       ),
                              
                              tabPanel("Two-Arm",
                                       sidebarLayout(
                                         sidebarPanel(
                                           p("Beta prior (in blue)", style="color: blue"),
                                           p("Binomial likelihood (in green)", style="color: green"),
                                           p("Beta posterior (in red)", style="color: red"),
                                           numericInput("alphaG1", "alpha, shape parameter for Beta prior in group 1:", 
                                                        value = 1,min = 0.01, max = 200),
                                           numericInput("betaG1", "beta, shape parameter for Beta prior in group 1:", 
                                                        value = 1,min = 0.01, max = 200),
                                           numericInput("alphaG2", "alpha, shape parameter for Beta prior in group 2:", 
                                                        value = 1,min = 0.01, max = 200),
                                           numericInput("betaG2", "beta, shape parameter for Beta prior in group 2:", 
                                                        value = 1,min = 0.01, max = 200),
                                           numericInput("rG1", "Number of observed responses in group 1:", 25,
                                                        min = 1, max = 100),
                                           numericInput("nG1", "Number of participants enrolled in group 1:", 35,
                                                        min = 1, max = 100),
                                           numericInput("rG2", "Number of observed responses in group 2:", 15,
                                                        min = 1, max = 100),
                                           numericInput("nG2", "Number of participants enrolled in group 2:", 35,
                                                        min = 1, max = 100),
                                           numericInput("pG2", "Target treatment difference between groups 
                                                        (eg., p1 - p2 >= delta):", 0.15,
                                                        min = 0, max = 1),
                                           numericInput("nsim", "Number of replicas for Monte Carlo sampling 
                                                        from the two posteriors", 10000),
                                           actionButton(inputId = "plotG2",
                                                        label = "Draw")
                                         ),
                                         
                                         mainPanel(#plotOutput("plotG2"),
                                                   
                                                   textOutput("description_G2"))
                                         )
                                       )
                              )
                 )


# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {
    
    ###Bayesian Graphica Illustration
    
    
    drawG <- eventReactive(input$plotG, {
        pG = input$pG
        
        priorG = list(input$alphaG, input$betaG)
        posteriorG = list((input$alphaG + input$rG), (input$betaG + input$nG - input$rG))
        
        xG = seq(0, 1, 0.01)
        max1 = max(dbeta(xG, priorG[[1]], priorG[[2]]))
        max2 = max(dbeta(xG, posteriorG[[1]], posteriorG[[2]]))
        maxy = min(8, max(max1, max2) + 0.1)
        
        # pesky detail to rescale the likelihood
        maxbinom = dbinom(input$rG, input$nG, input$rG/input$nG)
        mle = makeFun(0.99*maxy/maxbinom*choose(nG, rG)*p^rG*(1-p)^(nG-rG) ~ p, nG=input$nG, rG=input$rG)
        
        # display the results
        
        l1 = stat_function(fun = function(x) dbeta(x, posteriorG[[1]], posteriorG[[2]]), 
                          col = "red", lwd = 1.5)
        
        l2 = stat_function(fun = function(x) dbeta(x, posteriorG[[1]], posteriorG[[2]]), 
                          geom = "area", xlim = c(pG, 1), col = "red", fill = "red", alpha = 0.5)
            
        l3 = stat_function(fun = function(x) dbeta(x,priorG[[1]], priorG[[2]]), 
                          col = "blue", lwd = 1, alpha = 0.6)
        
        l4 =  stat_function(fun = function(x) mle(x), col = "green", lwd = 1, alpha = 0.6) 
        
        l5 = geom_vline(xintercept = pG, col = "black", lty = 2, lwd = 1.2)
        
        l6 = geom_text(aes(x = pG, y = maxy/2, label = "target response rate"),
                       col = "red", angle = 90, vjust = -0.4)
        
        ggplot(data.frame(x = c(0, 1)), ylim = c(0, maxy), aes(x)) + 
            theme_bw() + l1 + l2 + l3 + l4 + l5 + l6
    })
    
    output$plotG = renderPlot({drawG()})
    
    description_G = eventReactive(input$plotG, {
        pG = input$pG
        
        priorG = list(input$alphaG, input$betaG)
        posteriorG = list((input$alphaG + input$rG), (input$betaG + input$nG - input$rG))
        postPR = 1 - pbeta(pG, posteriorG[[1]], posteriorG[[2]])
        
        d1 = paste0("With prior Beta(", 
                    priorG[[1]], ",", priorG[[2]], "), and  ", 
                    input$rG, " responder(s) out of ", input$nG, " subjects, the posterior is Beta(",
                    posteriorG[[1]], ",", posteriorG[[2]], "). ",
                    "The posterior probability that the true response rate is greater than target rate (",
                    pG, ") is ", postPR, ".")
        
        d1
    })
    
    output$description_G = renderText({description_G()})
    
    description_G2 = eventReactive(input$plotG2, {
      pG2 = input$pG2
      
      priorG1 = list(input$alphaG1, input$betaG1)
      priorG2 = list(input$alphaG2, input$betaG2)
      
      posteriorG1 = list((input$alphaG1 + input$rG1), (input$betaG1 + input$nG1 - input$rG1))
      posteriorG2 = list((input$alphaG2 + input$rG2), (input$betaG2 + input$nG2 - input$rG2))
      
      theta1 = rbeta(input$nsim, posteriorG1[[1]], posteriorG1[[2]])
      theta2 = rbeta(input$nsim, posteriorG2[[1]], posteriorG2[[2]])
      
      diff = theta1 - theta2
      
      postPR2 = mean(diff >= pG2)
      
      
      d2 = paste0("With priors Beta(", 
                  priorG1[[1]], ",", priorG1[[2]], ") and  Beta(", 
                  priorG2[[1]], ",", priorG2[[2]], ") for group 1 and 2, respectively. ",
                  "After ", input$rG1, " responder(s) out of ", input$nG1, " subjects in group 1, and ",
                  input$rG2, " responder(s) out of ", input$nG2, " subjects in group 2, ",
                  "the posteriors are Beta(",
                  posteriorG1[[1]], ",", posteriorG1[[2]], ") and Beta(",
                  posteriorG2[[1]], ",", posteriorG2[[2]], ") for group 1 and 2, respectively. ",
                  "The posterior probability that the true treatment different between groups is at least ",
                  pG2, " is ", postPR2, ", under ", input$nsim, " replicas from Monte Carlo sampling from the two posteriors.")
      
      d2
    })
    
    output$description_G2 = renderText({description_G2()})
}

# Create Shiny app ----
shinyApp(ui, server)