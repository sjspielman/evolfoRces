library(tidyverse)
library(ggrepel)
library(cowplot)
library(viridis) 
library(plotly)
library(colourpicker)
library(shinythemes)
library(DT)
library(promises)
library(future)
plan(multiprocess)

max.gen <- 1000
max.rep <- 50
max.n   <- 100000   
round_table <- 4
source("./simul.R")    
shinyServer(function(input,output,session){

    
                 
                                                         
    ####################### SINGLE POPULATION TAB ########################################
    ######################################################################################



    line.size <- 1
    t1 <- 12
    t2 <- 14
    plot.simulation.single <- function(sim.data, gen, plottype)
    {

        theme_set(theme_cowplot() + theme(legend.position = "none",
                                          plot.margin = margin(.1, 2, .1, .1, "cm"),))#, 
                                           # axis.text = element_text(size=t1),
                                           # axis.title = element_text(size=t2)))
                                            
       # sim.data %>% mutate(display_text_freq = paste("Frequency of A:", round(sim.data$p, 8)),
       #                     display_text_fit = paste("Mean population fitness:", round(sim.data$w, 8))) -> sim.data

        sim.data %>% 
            rename("Simulation Replicate" = population) %>%
            mutate(repel_label = paste0("p = ", round(p, 4), "; w = ", round(w, 4))) -> sim.data
        
        if (plottype == "frequency")
        {
        p1 <- ggplot(sim.data, aes(x = generation, y = p, group = `Simulation Replicate`, color = `Simulation Replicate`)) + 
            geom_path(size=line.size) + 
            scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
            scale_x_continuous(limits=c(0,gen+10), expand=c(0,0)) + 
            scale_color_viridis() + 
            background_grid() + 
            xlab("Generation") + ylab("Frequency of allele A") +
            geom_label_repel(data = subset(sim.data, generation == max(generation)), 
                  aes(label = repel_label),
                  nudge_x = 1,
                  na.rm = TRUE)
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
        p1
        #ggplotly(p1, tooltip = c("x", "y", "group"))       
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


    ###################### Generate simulation data upon go button ########################
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
        
        sim_data <- future(
          simulatePopulations.single(gen=gen,
                                     p=p,
                                     Waa=Waa,
                                     Wab=Wab,
                                     Wbb=Wbb,
                                     Uab=Uab,
                                     Uba=Uba,
                                     Neff=Neff,
                                     infinitePop=infinitePop,
                                     nRep=nRep))
      
        
        output$result_header_s <- renderText({
            paste("Simulation results after", gen, "generations of evolution:")
        })

        output$result_table_s <- renderDT(rownames= FALSE, server=FALSE, options = list(dom = 't'),
            sim_data %...>% process.simulation(round_table)
        )
            
        output$singleplot.frequency_s <- renderPlot( { 
            sim_data %...>% plot.simulation.single(gen, "frequency")
        })

        output$singleplot.fitness_s <- renderPlot( { 
            sim_data %...>% plot.simulation.single(gen, "fitness")
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


