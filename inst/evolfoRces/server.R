library(tidyverse)
library(cowplot)
library(shinyWidgets)
library(colourpicker)
library(DT)
library(promises)
library(future)
plan(multiprocess)

options(htmlwidgets.TOJSON_ARGS = list(na = 'string')) ## Make NA in DT show as NA instead of blank cell

shinyServer(function(input,output,session){
    source("./simul.R")
    

    current_simulation <- reactiveValues()
    stored_simulation <- reactiveValues()
           
                                          
    ####################### SINGLE POPULATION TAB ########################################
    ######################################################################################


    ### Simulate, plot, summarize reactive to "Run Simulation!"
    observeEvent(input$go_s,  {

        # isolate simulation parameters
        gen <- as.numeric(isolate(input$gen_s))
        p   <- as.numeric(isolate(input$p))
        Waa <- as.numeric(isolate(input$Waa_s))
        Wab <- as.numeric(isolate(input$Wab_s))
        Wbb <- as.numeric(isolate(input$Wbb_s))
        Uab <- as.numeric(isolate(input$Uab))
        Uba <- as.numeric(isolate(input$Uba))


        usedrift <- ifelse( isolate(input$usedrift) == "drift", TRUE, FALSE )  
        infinitePop <- ifelse(usedrift==FALSE,TRUE,FALSE)
        if (infinitePop)
        {   
            Neff <- 1  ## doesnt matter
            nRep <- 1
        } else
        {
            Neff <- as.numeric(isolate(input$N))
            nRep <- as.numeric(isolate(input$nRep))
        }
        
        ## Reset some limits ##
        if (gen > max.gen) gen = max.gen
        if (nRep > max.rep) nrep = max.rep
        if (Neff > max.n) Neff = max.n
    
            

       simulation_data <- future( simulatePopulations.single(gen=gen,
                                                             p=p,
                                                             Waa=Waa,
                                                             Wab=Wab,
                                                             Wbb=Wbb,
                                                             Uab=Uab,
                                                             Uba=Uba,
                                                             Neff=Neff,
                                                             infinitePop=infinitePop,
                                                             nRep=nRep))
                                

        current_simulation$data <- simulation_data
        current_simulation$gen <- gen
        current_simulation$infinitePop <- infinitePop


        output$single_table <- renderDT(rownames= FALSE, server=FALSE, 
                                options = list(dom = 'tp', columnDefs = list(list(className = 'dt-left', targets = "_all"))),
            simulation_data %...>% process_simulation(gen)
        )
        
          
        output$single_plot <- renderPlot( { 
            simulation_data %...>% {
                simulation_data_df <- .
                plot_simulation(simulation_data_df, gen, input$line_color_s, infinitePop)
            }
        })

    })
    
    
    observeEvent(input$store_single_btn,  {
        inputSweetAlert(
            session = session, inputId = "store_single_name",
            title = "Name for your stored results?",
            type = "warning"
        )    
       stored_simulation$data        <- current_simulation$data
       stored_simulation$gen         <- current_simulation$gen
       stored_simulation$infinitePop <- current_simulation$infinitePop     
       stored_simulation$line_color  <- input$line_color_s
    })
        
    observeEvent(input$store_single_name, {
    
        output$single_plot_stored <- renderPlot({
            stored_simulation$data %...>% {
                stored_simulation_df <- .
                plot_simulation(stored_simulation_df, stored_simulation$gen, stored_simulation$line_color, stored_simulation$infinitePop)
            }
        })        

        output$single_table_beforeafter_stored <- renderDT(rownames= FALSE, server=FALSE, 
                                                    options = list(dom = 'tp', columnDefs = list(list(className = 'dt-left', targets = "_all"))),
            stored_simulation$data %...>% {
                stored_simulation_df <- .
                process_simulation(stored_simulation_df, stored_simulation$gen)
            }
        )
                
        output$single_table_fixation_stored <- renderDT(rownames= FALSE, server=FALSE, 
                                                    options = list(dom = 'tp', columnDefs = list(list(className = 'dt-left', targets = "_all"))),
            stored_simulation$data %...>% {
                stored_simulation_df <- .
                process_simulation(stored_simulation_df, stored_simulation$gen)
            }
        )
        
        output$single_name_stored <- renderText({input$store_single_name})

    })    
    
    observeEvent(input$clear_single_btn, {
        stored_simulation <- reactiveValues()
        output$single_name_stored  <- renderText({})
        output$single_table_stored <- renderDT({})
        output$single_plot_stored <- renderPlot({})
    })
    
    
    ######################################################################################
    ############################### MIGRATION TAB ########################################
    ######################################################################################
    observeEvent(input$go_m,  {

                                  
                                  
        # isolate simulation parameters..and color!!
        gen <- as.numeric(isolate(input$gen_m))
        p.main <- as.numeric(isolate(input$p.main))
        p.island <- as.numeric(isolate(input$p.island))
        m <- as.numeric(isolate(input$m))
        Waa <- as.numeric(isolate(input$Waa_m))
        Wab <- as.numeric(isolate(input$Wab_m))
        Wbb <- as.numeric(isolate(input$Wbb_m))
    

        ## Reset some limits ##
        if (gen > max.gen) gen = max.gen
        infinitePop <- TRUE

        # simulate
        simulation_data <- future( simulatePopulations.migration(gen=gen,
                                                  p.main = p.main,
                                                  p.island = p.island,
                                                  m = m,
                                                  Waa=Waa,
                                                  Wab=Wab,
                                                  Wbb=Wbb) %>% 
                                     mutate("population" = 1))


        current_simulation$data <- simulation_data
        current_simulation$gen <- gen
        current_simulation$infinitePop <- infinitePop


        output$migration_table <- renderDT(rownames= FALSE, server=FALSE, 
                                options = list(dom = 'tp', columnDefs = list(list(className = 'dt-center', targets = "_all"))),
            simulation_data %...>% process_simulation(gen)
        )
        
          
        output$migration_plot <- renderPlot( { 
            simulation_data %...>% {
                simulation_data_df <- .
                plot_simulation(simulation_data_df, gen, input$line_color_m, infinitePop)
            }
        })

    })
    
    
    observeEvent(input$store_migration_btn,  {
        inputSweetAlert(
            session = session, inputId = "store_migration_name",
            title = "Name for your stored results?",
            type = "warning"
        )    
       stored_simulation$data        <- current_simulation$data
       stored_simulation$gen         <- current_simulation$gen
       stored_simulation$infinitePop <- current_simulation$infinitePop     
       stored_simulation$line_color  <- input$line_color_m
    })
        
    observeEvent(input$store_migration_name, {
    
        output$migration_plot_stored <- renderPlot({
            stored_simulation$data %...>% {
                stored_simulation_df <- .
                plot_simulation(stored_simulation_df, stored_simulation$gen, stored_simulation$line_color, stored_simulation$infinitePop)
            }
        })        
        
        output$migration_table_stored <- renderDT(rownames= FALSE, server=FALSE, 
                                                    options = list(dom = 'tp', columnDefs = list(list(className = 'dt-center', targets = "_all"))),
            stored_simulation$data %...>% {
                stored_simulation_df <- .
                process_simulation(stored_simulation_df, stored_simulation$gen)
            }
        )
        
        output$migration_name_stored <- renderText({input$store_migration_name})

    })    
    
    observeEvent(input$clear_migration_btn, {
        stored_simulation <- reactiveValues()
        output$migration_name_stored  <- renderText({})
        output$migration_table_stored <- renderDT({})
        output$migration_plot_stored <- renderPlot({})
    })
    


})

