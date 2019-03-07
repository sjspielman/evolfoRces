library(tidyverse)
library(cowplot)
library(viridis) 
library(plotly)
library(colourpicker)
library(shinythemes)
options(digits = 8)

shinyServer(function(input,output,session){
    source("./simul.R")
    
    max.gen <- 1000
    max.rep <- 50
    max.n   <- 100000                    
                        
    process.simulation <- function(sim.data, gen)
    {
        ## Table of final generation
        sim.data %>% filter(generation == gen) -> final.generation

        final.generation %>% mutate("het" = 2*p*(1-p)) -> final.generation


        result.header <- paste("Simulation results after", gen, "generations of evolution:")
        result.table <- tibble("Simulation Replicate" = final.generation$population,
                               "Allele A frequency" = final.generation$p, 
                               "Population fitness" = final.generation$w,
                               "Population heterozygosity" = final.generation$het)


        #write_csv(sim.data, "sim.csv")
        
        sim.data %>% 
            group_by(population) %>%
            filter(generation == gen) %>%
            mutate(fixed = ifelse(p == 1, TRUE, FALSE),
                   lost = ifelse(p == 0, TRUE, FALSE)) -> sim.data.fixation
        
                            
        #### Check for fixation ####
        wefixed <- sum(sim.data.fixation$fixed + sim.data.fixation$lost) > 0


        if (wefixed)
        {
        
            fix.loss.table <- tibble("rep" = as.integer(), "allele" = as.character(), "gen" = as.integer())
        
            if (sum(sim.data.fixation$fixed) > 0){
    
                sim.data.fixation %>%
                    ungroup() %>%
                    filter(fixed == TRUE) %>%
                    select(population) -> fixed.pops
                
                for (i in unique(fixed.pops$population))
                {
                    sim.data %>% 
                        filter(population == i, p == 1) %>%
                        select(generation) -> s2
                    when.fixed <- min(s2$generation)
                    fix.loss.table <- bind_rows(fix.loss.table, tibble("rep" = i, "allele" = "A", "gen" = when.fixed))
                }
            }
            if (sum(sim.data.fixation$lost) > 0){
    
                sim.data.fixation %>%
                    ungroup() %>%
                    filter(lost == TRUE) %>%
                    select(population) -> lost.pops
                
                for (i in unique(lost.pops$population))
                {
                    sim.data %>% 
                        filter(population == i, p == 0) %>%
                        select(generation) -> s2
                    when.lost <- min(s2$generation)
                    fix.loss.table <- bind_rows(fix.loss.table, tibble("rep" = i, "allele" = "a", "gen" = when.lost))
                }
            }
            
            
     
            fix.loss.table$allele <- factor(fix.loss.table$allele, levels=c("A", "a"))
            fix.loss.table %>%
                arrange(allele, gen) %>% 
                rename("Simulation Replicate" = rep, 
                       "Allele Fixed" = allele, 
                       "Generation Fixed" = gen) -> fix.loss.table

            result.table <- left_join(result.table, fix.loss.table) %>% arrange(`Allele Fixed`)
        } else
        {
            result.table <- result.table %>% mutate("Allele Fixed" = "NA", "Generation Fixed" = "NA")
        }
    
        processed <- list(result.header, result.table)
        
    }                                
                
                
                                          
    ####################### SINGLE POPULATION TAB ########################################
    ######################################################################################



    line.size <- 1
    t1 <- 12
    t2 <- 14
    plot.simulation.single <- function(sim.data, gen, plottype)
    {

        theme_set(theme_cowplot() + theme(legend.position = "none", 
                                            axis.text = element_text(size=t1),
                                            axis.title = element_text(size=t2)))
                                            
        sim.data %>% mutate(display_text_freq = paste("Frequency of A:", round(sim.data$p, 8)),
                            display_text_fit = paste("Mean population fitness:", round(sim.data$w, 8))) -> sim.data
        
        

        sim.data %>% rename("Simulation Replicate" = population) -> sim.data
        if (plottype == "frequency")
        {
        p1 <- ggplot(sim.data, aes(x = generation, y = p, group = `Simulation Replicate`, color = `Simulation Replicate`)) + 
            geom_path(size=line.size) + 
            scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) + 
            scale_color_viridis() + 
            background_grid() + 
            xlab("Generation") + ylab("Frequency of allele A")
        }
        if (plottype == "fitness")
        {
        p1 <- ggplot(sim.data, aes(x = generation, y = w, group = `Simulation Replicate`, color = `Simulation Replicate`)) + 
            geom_path(size=line.size) + 
            scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) + 
            scale_color_viridis() + 
            background_grid() + 
            xlab("Generation") + ylab("Mean population fitness")
        }
        ggplotly(p1, tooltip = c("x", "y", "group"))       
    }



    plot.simulation.migration <- function(sim.data, gen, line.color, plottype)
    {
        theme_set(theme_cowplot() + theme(axis.text = element_text(size=t1),
                                            axis.title = element_text(size=t2),
                                            legend.position = "none"))


        if (plottype == "frequency")
        {
         p2 <- ggplot(sim.data, aes(x = generation, y = p)) + 
            geom_path(size = line.size, color = line.color) + 
            scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) +  
            background_grid() + 
            xlab("Generation") + ylab("Frequency of allele A (island)")
        }
        if (plottype == "fitness")
        {
        p2 <- ggplot(sim.data, aes(x = generation, y = w)) + 
            geom_path(size = line.size, color = line.color) +   
            scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) + 
            background_grid() + 
            xlab("Generation") + ylab("Mean population fitness (island)")
        }
        ggplotly(p2)
    }        



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
        if (gen > 1000) gen = max.gen
        if (nRep > 100) nrep = max.rep
        if (Neff > 100000) Neff = max.n
    
            

        sim.data <- simulatePopulations.single(gen=gen,
                                          p=p,
                                          Waa=Waa,
                                          Wab=Wab,
                                          Wbb=Wbb,
                                          Uab=Uab,
                                          Uba=Uba,
                                          Neff=Neff,
                                          infinitePop=infinitePop,
                                          nRep=nRep)
        
        
        
        processed <- process.simulation(sim.data, gen)
        
        
         
        result.header <- processed[[1]]
        result.table  <- processed[[2]]

        result.table$`Simulation Replicate` <- as.integer(result.table$`Simulation Replicate`)
        result.table$`Generation Fixed` <- as.integer(result.table$`Generation Fixed`)
        
        result.table %>%
            select(-`Simulation Replicate`, everything()) -> result.table
        output$result_header_s <- renderText({result.header})

        output$result_table_s <- renderTable(
            {as.data.frame( result.table )}, 
            striped=TRUE, hover=TRUE, bordered=TRUE, align="l", digits=8, rownames=TRUE)
            
        
        output$singleplot.frequency_s <- renderPlotly( { 
            plot.simulation.single(sim.data, gen, "frequency")
        })

        output$singleplot.fitness_s <- renderPlotly( { 
            plot.simulation.single(sim.data, gen, "fitness")
        })

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
    
        line.color <- isolate(input$yaycolor)
        # simulate
        sim.data <- simulatePopulations.migration(gen=gen,
                                                  p.main = p.main,
                                                  p.island = p.island,
                                                  m = m,
                                                  Waa=Waa,
                                                  Wab=Wab,
                                                  Wbb=Wbb) %>% mutate("population" = 1)

        processed <- process.simulation(sim.data, gen)
             
        result.header <- processed[[1]]
        result.table  <- processed[[2]]
        result.table$`Simulation Replicate` <- as.integer(result.table$`Simulation Replicate`)
        result.table$`Generation Fixed`     <- as.integer(result.table$`Generation Fixed`)
        result.table %>%
            select(-`Simulation Replicate`, everything()) -> result.table
        
        output$result_header_m <- renderText({result.header})

        output$result_table_m <- renderTable(
            {as.data.frame( result.table )}, 
            striped=TRUE, hover=TRUE, bordered=TRUE, align="l", digits=8, rownames=TRUE)
            
        
        output$singleplot.frequency_m <- renderPlotly( { 
            plot.simulation.migration(sim.data, gen, line.color, "frequency")
        })

        output$singleplot.fitness_m <- renderPlotly( { 
            plot.simulation.migration(sim.data, gen, line.color, "fitness")
        })




    })


})


