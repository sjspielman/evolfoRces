library(tidyverse)
library(cowplot)
library(viridis) 

shinyServer(function(input,output,session){
  
    source("./simul.R")
  
    ####################### SINGLE POPULATION TAB ########################################
    ######################################################################################


    plot.simulation.single <- function(sim.data, gen, line.size, t1, t2)
    {
        theme_set(theme_cowplot() + theme(legend.position = "None", 
                                            axis.text = element_text(size=t1),
                                            axis.title = element_text(size=t2)))

        freq.plot <- ggplot(sim.data, aes(x = generation, y = p, group = population, color = population)) + 
            geom_path(size=line.size) + 
            scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen), expand=c(0,0)) + 
            scale_color_viridis() + 
            background_grid() + 
            xlab("Generation") + ylab("Frequency of allele A")
        
        fit.plot <- ggplot(sim.data, aes(x = generation, y = w, group = population, color = population)) + 
            geom_path(size=line.size) + 
            scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen), expand=c(0,0)) + 
            scale_color_viridis() + 
            background_grid() + 
            xlab("Generation") + ylab("Mean population fitness")
    
        p <- plot_grid(freq.plot, fit.plot, nrow=2, scale=0.95)
        print(p)
    }

    plot.simulation.migration <- function(sim.data, gen, line.color, line.size, t1, t2)
    {
        theme_set(theme_cowplot() + theme(axis.text = element_text(size=t1),
                                            axis.title = element_text(size=t2),
                                            legend.position = "none"))


        freq.plot <- ggplot(sim.data, aes(x = generation, y = p)) + 
            geom_path(size = line.size, color = line.color) + 
            scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen), expand=c(0,0)) +  
            background_grid() + 
            xlab("Generation") + ylab("Frequency of allele A (island)")
        
        fit.plot <- ggplot(sim.data, aes(x = generation, y = w)) + 
            geom_path(size = line.size, color = line.color) +   
            scale_y_continuous(limits=c(0,1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen), expand=c(0,0)) + 
            background_grid() + 
            xlab("Generation") + ylab("Mean population fitness (island)")
    
        p <- plot_grid(freq.plot, fit.plot, nrow=2, scale=0.95)
        print(p)
    }        



    ### Simulate, plot, summarize reactive to "Run Simulation!"
    observeEvent(input$go_s,  {

        # isolate simulation parameters
        gen2 <- as.numeric(isolate(input$gen_s))
        p2 <- as.numeric(isolate(input$p))
        Waa2 <- as.numeric(isolate(input$Waa_s))
        Wab2 <- as.numeric(isolate(input$Wab_s))
        Wbb2 <- as.numeric(isolate(input$Wbb_s))
        Uab2 <- as.numeric(isolate(input$Uab))
        Uba2 <- as.numeric(isolate(input$Uba))
        Neff2 <- as.numeric(isolate(input$Neff))
        nRep2 <- as.numeric(isolate(input$nRep))
        infinitePop2 <- as.numeric(isolate(input$infinitePop))   
    
        # Override for infinite populations to one replicate, reduces comp load for deterministic.
        if(infinitePop2) nRep2 = 1
        ## Reset some limits ##
        if (gen2 > 1000) gen2 = 1000
        if (nRep2 > 100) nrep2 = 100
        if (Neff2 > 1000000) Neff2 = 1000000
    
            
        # simulate
        sim.data <- simulatePopulations.single(gen=gen2,
                                          p=p2,
                                          Waa=Waa2,
                                          Wab=Wab2,
                                          Wbb=Wbb2,
                                          Uab=Uab2,
                                          Uba=Uba2,
                                          Neff=Neff2,
                                          infinitePop=infinitePop2,
                                          nRep=nRep2)

        ################ If NO DRIFT, what are the final frequencies? #####################
        if (infinitePop2)
        {
            sim.data %>% filter(generation == gen2) -> final.generation
        
            fix.loss.table <- tibble("Allele" = "A", "Final population frequency" = final.generation$p, "Final population fitness" = final.generation$w)

            if ( final.generation$p > ZERO && (1-final.generation$p) > ZERO) 
            {
                output$fixation_s <- renderText({
                    "Neither allele A nor B fixed in your simulation."
                })
                output$fixlosstable <- renderTable({as.data.frame( fix.loss.table )}, striped=TRUE, hover=TRUE, bordered=TRUE, align="l", digits=8)
            }
            else 
            {
            
                if ( final.generation$p <= ZERO){
                    allele = "B"
                }
                else {
                    allele = "A"
                }
                sim.data %>% filter(p == final.generation$p) -> sim.data.fixed
                fixedgen <- min(sim.data.fixed$generation)
                output$fixation_s <- renderText({
                    paste0("Allele ", allele, " fixed in generation ", fixedgen, " of your simulation, for a final population mean fitness of ", final.generation$w, ".")
                })
                output$fixlosstable <- renderTable({})
            }
        }   
        else
        { 
        
            ################ If DRIFT, fixation events at END of simulation?
            sim.data %>% filter(generation == gen2, p == 1) -> fixed.data
            sim.data %>% filter(generation == gen2, p == 0) -> lost.data
            nfixed <- nrow(fixed.data)
            nlost  <- nrow(lost.data)

        
            wefixed <- sum(nfixed + nlost) > 0
            # If we have fixation, make a table
            if (wefixed)
            {
        
                #### There is about a 100% chance this can be vectorized. Alas, it is dinner time. ####
                sim.data %>%
                    filter(generation > 0, population %in% fixed.data$population) -> sim.fixed
                sim.data %>%
                    filter(generation > 0, population %in% lost.data$population) -> sim.lost
        
                fix.loss.table <- tibble("allele" = as.character(), "gen" = as.integer())
        
        
        
                for (i in unique(sim.fixed$population))
                {
                    this <- sim.fixed %>% filter(population == i) %>% select(p)
                    when.fixed <- max(which(this$p != 1))
                    fix.loss.table <- bind_rows(fix.loss.table, tibble("allele" = "A", "gen" = when.fixed))
                }
                for (i in unique(sim.lost$population))
                {
                    this <- sim.lost %>% filter(population == i) %>% select(p) %>% mutate(p = 1-p)
                    when.fixed <- max(which(this$p != 1))
                    fix.loss.table <- bind_rows(fix.loss.table, tibble("allele" = "B", "gen" = when.fixed))
                }       
                fix.loss.table %>%
                    arrange(allele, gen) %>% 
                    group_by(allele) %>%
                    mutate("Replicate" = 1:n()) -> fix.loss.table
                fix.loss.table <- fix.loss.table[,c(3,1,2)]
                names(fix.loss.table) <- c("Replicate", "Allele Fixed", "Generation Fixed")

                
                output$fixation_s <- renderText({
                paste0("At the end of your simulation, \n Allele A fixed ", nfixed, " times.\n Allele B fixed ", nlost, " times.")})
            }
            else {
                output$fixation_s <- renderText({
                    "Neither allele A nor B fixed in any simulation replicate."
                })

                sim.data %>% filter(generation == gen2) -> final.generation
                fix.loss.table <- tibble("Allele" = "A", "Final population frequency" = final.generation$p, "Final population fitness" = final.generation$w) %>% 
                    arrange(`Final population frequency`) %>%
                    mutate("Replicate" = 1:n())
                fix.loss.table <- fix.loss.table[,c(4,1,2,3)]
            }
            output$fixlosstable_s <- renderTable(
                    {as.data.frame( fix.loss.table )}, striped=TRUE, hover=TRUE, bordered=TRUE, align="l", digits=8)
        }

        output$singleplot <- renderPlot( { 
            plot.simulation.single(sim.data, gen2, 1.15, 16, 19)
        })
        
        output$downloadp_s <- renderUI({
            downloadButton('download_plot_s', 'Download plots')
        })
        output$download_plot_s <- downloadHandler(
            filename = function() {
                "download.png"
            },
            content = function(file) {
                save_plot(file, plot.simulation.single(sim.data, gen2, 1, 13, 14), base_width = 8, base_height=6)
            }
        )
        
        output$downloaddata_s <- renderUI({
            downloadButton('download_data_s', 'Download simulation data')
        })
        output$download_data_s <- downloadHandler(
            filename = function() {
                "simulation_data.csv"
            },
            content = function(file) {
                write_csv(sim.data, file) 
            }
        )

    })
    
    
    ############################### MIGRATION TAB ########################################
    ######################################################################################
    observeEvent(input$go_m,  {

                                  
                                  
        # isolate simulation parameters..and color!!
        gen2 <- as.numeric(isolate(input$gen_m))
        p.main2 <- as.numeric(isolate(input$p.main))
        p.island2 <- as.numeric(isolate(input$p.island))
        m2 <- as.numeric(isolate(input$m))
        Waa2 <- as.numeric(isolate(input$Waa_m))
        Wab2 <- as.numeric(isolate(input$Wab_m))
        Wbb2 <- as.numeric(isolate(input$Wbb_m))
    
        line.color <- isolate(input$yaycolor)
        # simulate
        sim.data <- simulatePopulations.migration(gen=gen2,
                                          p.main = p.main2,
                                          p.island = p.island2,
                                          m = m2,
                                          Waa=Waa2,
                                          Wab=Wab2,
                                          Wbb=Wbb2)

        ## how many fixation events at END of simulation? We are only looking into dynamics **on the island**
        sim.data %>% filter(generation > 0) -> sim.data.nozero
        
        
        ## Deterministic simulation w/out drift so much easier to find generations.
        sim.data.nozero %>% filter(generation == gen2) -> final.generation
        
        fixed <- FALSE
        lost <- FALSE
        nofixation <- FALSE
        if (final.generation$p == 1)
        {
            fixed <- TRUE
            finalgen <- max(which(sim.data.nozero$p != 1)) + 1
        }
        else if (final.generation$p == 0)
        {
            lost <- TRUE
            this <- sim.data.nozero %>% mutate(p = 1 - p)
            finalgen <- max(which(this$p != 1)) + 1
        }
        else {
            nofixation <- TRUE
        }
        
        


        output$migrationplot <- renderPlot( { 
            plot.simulation.migration(sim.data, gen2, line.color, 1.25, 16, 19)
        })

        if (fixed) output$fixation <- renderText({ paste0("Allele A fixed in generation ", finalgen, " on the island. Allele frequencies up until allele A fixation are in the table below.")})
        if (lost) output$fixation <- renderText({ paste0("Allele B fixed in generation ", finalgen, " on the island. Allele frequencies up until allele B fixation are in the table below.")})
        if (nofixation) output$fixation <- renderText({ "No allele went to fixation on the island. Allele A frequencies in the final generation are in table below."})
        
        if (!(nofixation)){
            fix.loss.table <- sim.data.nozero %>% filter(generation <= finalgen)
            fix.loss.table$generation <- as.integer(fix.loss.table$generation)
            names(fix.loss.table) <- c("Generation", "Frequency of allele A", "Mean population fitness")
            output$fixlosstable <- renderTable(
                {as.data.frame( fix.loss.table )}, striped=TRUE, hover=TRUE, bordered=TRUE, align="l", digits=8)
        }
        else {           
            final.generation$generation <- as.integer(final.generation$generation)
            names(final.generation) <- c("Generation", "Frequency of allele A", "Mean population fitness")
            output$fixlosstable_m <- renderTable({final.generation}, hover=TRUE, bordered=TRUE, align="l", digits=8)
        }

        output$downloadp_m <- renderUI({
            downloadButton('download_plot_m', 'Download plots')
        })
        output$download_plot_m <- downloadHandler(
            filename = function() {
                "download.png"
            },
            content = function(file) {
                save_plot(file, plot.simulation.migration(sim.data, gen2, line.color, 1.15, 13, 14), base_width = 8, base_height=6)
            }
        )

        output$downloaddata_m  <- renderUI({
            downloadButton('download_data_m', 'Download simulation data')
        })
        output$download_data_m <- downloadHandler(
            filename = function() {
                "simulation_data."
            },
            content = function(file) {
                write_csv(sim.data, file) 
            }
        )



    })


})


