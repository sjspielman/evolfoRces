# SJS modification of original `dev.R` from driftR package
library(tidyverse)
library(viridis)

ZERO <- 1e-7
max.gen <- 1000
max.rep <- 50
max.n   <- 100000                    


line.size <- 1
t1 <- 12
t2 <- 14
theme_set(theme_classic() + theme(legend.position = "none", 
                                  axis.text = element_text(size=t1),
                                  axis.title = element_text(size=t2)))

plot_simulation <- function(sim.data, gen, line_color, is_infinite)
{

                                        
    sim.data %>% mutate(display_text_freq = paste("Frequency of A:", round(sim.data$p, 8)),
                        display_text_fit = paste("Mean population fitness:", round(sim.data$w, 8))) %>% 
                        rename("Simulation Replicate" = population) -> sim.data
    
    if (is_infinite) sim.data$`Simulation Replicate` <- factor(sim.data$`Simulation Replicate`)


    p1 <- ggplot(sim.data, aes(x = generation, y = p, group = `Simulation Replicate`, color = `Simulation Replicate`)) + 
        geom_path(size=line.size) + 
        scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
        scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) + 
        background_grid() + 
        xlab("Generation") + ylab("Frequency of allele A")

    p2 <- ggplot(sim.data, aes(x = generation, y = w, group = `Simulation Replicate`, color = `Simulation Replicate`)) + 
        geom_path(size=line.size) + 
        scale_y_continuous(limits=c(0,1.1), expand=c(0,0)) + 
        scale_x_continuous(limits=c(0,gen+5), expand=c(0,0)) + 
        background_grid() + 
        xlab("Generation") + ylab("Mean population fitness")


    if (is_infinite) {
        p1 <- p1 + scale_color_manual(values = line_color)
        p2 <- p2 + scale_color_manual(values = line_color)
    } else {
        p1 <- p1 + scale_color_viridis()
        p2 <- p2 + scale_color_viridis()
    }    
    plot_grid(p1, p2, nrow=1, scale = 0.9)
}



process_simulation <- function(sim.data, gen, infinitePop)
{
    ## before/after basics
    sim.data %>% 
        filter(generation %in% c(0,gen)) %>% 
        mutate(generation = ifelse(generation == 0, "Before simulation", "After simulation"),
               het = 2*p*(1-p)) %>%
        rename("Simulation Replicate" = population,
               "Time Point" = generation,
               "Allele A frequency" = p,
               "Population fitness" = w,
               "Population heterozygosity" = het) -> result.table
    
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
                   "Generation Fixed" = gen) %>%
            mutate("Time Point" = "After simulation") -> fix.loss.table

        result.table <- left_join(result.table, fix.loss.table) 
    } else
    {
        result.table <- result.table %>% mutate("Allele Fixed" = "NA", "Generation Fixed" = "NA")
    }

    result.table %>%
        mutate(`Time Point`         = factor(result.table$`Time Point`, levels=c("Before simulation", "After simulation")),
               `Allele A frequency` = round(`Allele A frequency`, 5), 
               `Population fitness` = round(`Population fitness`, 5), 
               `Population heterozygosity` = round(`Population heterozygosity`, 5),
               `Allele Fixed` = case_when(`Time Point` == "Before simulation" ~ " ", 
                                          `Time Point` == "After simulation" & `Allele Fixed` == "A" ~ "A",
                                          `Time Point` == "After simulation" & `Allele Fixed` == "a" ~ "a"),
               `Generation Fixed` = ifelse(`Time Point` == "Before simulation", " ", `Generation Fixed`)) %>%
        arrange(`Simulation Replicate`, `Time Point`) -> result.table
    
    if (!(infinitePop)){
        result.table %>%
            filter(`Time Point` == "After simulation") -> result.table
    }
    result.table
    
}                                
            
     

### WE ONLY ARE KEEPING TRACK OF THE ISLAND!!!!!!!!!!!!! ###
simulatePopulations.migration <- function(gen,
                                  Waa,
                                  Wab,
                                  Wbb,
                                  p.main,
                                  p.island,
                                  m)
{  

    
    ## Calculate starting fitness
    w.start <- (p.island**2)*Waa + (2*p.island*(1-p.island)*Wab) + ((1-p.island)**2)*Wbb
    
    full.sim <- tibble("generation" = 0,
                       "p" = p.island,
                       "w" = w.start)    
    
    
    withProgress(message="simulating populations...",value=0,{

        for(i in 1:gen)
        { 

            full.sim %>% filter(generation == i-1) -> last.generation  
    
            # migrate mainland to island. 
            p <- last.generation$p*(1-m) + p.main*m 
            
            # fitness after migration
            mean.w <- (p**2)*Waa + (2*p*(1-p)*Wab) + ((1-p)**2)*Wbb
            
            ## Ok seriously, rounding. Get it together R flops.
            if (1 - p <= ZERO) p <- 1.0
            if (p <= ZERO)    p <- 0.0
            
            
            # selection (reproduction) for next generation
            freq.aa <- (p*p*Waa)/mean.w
            freq.ab <- (2*p*(1-p)*Wab)/mean.w
            p <- freq.aa+(freq.ab/2)
            if (1 - p <= ZERO) p <- 1.0
            if ( p <= ZERO)    p <- 0.0
            mean.w <- (p**2)*Waa + (2*p*(1-p)*Wab) + ((1-p)**2)*Wbb
        

            ## Add new row
            full.sim <- bind_rows(full.sim,  tibble("generation" = i,
                                                    "p" = p,
                                                    "w" = mean.w))                                           
            incProgress(1/gen)
        
        } #end generations loop
        
    }) #end progress bar

    return(full.sim)
}



#main simulation function
simulatePopulations.single <- function(gen,
                                  p,
                                  Waa,
                                  Wab,
                                  Wbb,
                                  Uab,
                                  Uba,                      
                                  Neff,
                                  infinitePop,
                                  nRep)
{  
    
    ## Calculate starting fitness
    q <- 1 - p
    w.start <- (p**2)*Waa + (2*p*q*Wab) + (q**2)*Wbb
   
    full.sim <- tibble( "population" = 1:nRep,
                          "generation" = rep(0,nRep),
                          "p" = rep(p, nRep),
                          "w" = rep(w.start, nRep))
                          
    withProgress(message="simulating populations...",value=0,{

        for(i in 1:gen)
        { 

            full.sim %>% 
                filter(generation == i-1) -> last.generation
        
            p <- mean(last.generation$p)
    
            for (j in 1:nRep)
            {
                # starting frequency
                p <- last.generation$p[last.generation$population == j]
            
                # mutate
                p <- (1 - Uab) * p + Uba * (1 - p)
                q <- 1 - p

                # NOT fixed
                if( p > ZERO && (1 - p) > ZERO){
           
                    # last mean population fitness
                    mean.w <- p * p * Waa + (2*p*q*Wab) + q * q * Wbb

                    # post-selection genotype frequencies, weighted by relative fitness
                    freq.aa <- (p*p*Waa)/mean.w
                    freq.ab <- (2*p*q*Wab)/mean.w
                    ## Reproduction, with drift
                    if(infinitePop == FALSE)
                    { 
                        Naa <- rbinom(1, Neff, freq.aa)
                        if(freq.aa<1)
                        { 
                            Nab <- rbinom(1, (Neff - Naa), (freq.ab/(1 - freq.aa)) )
                        }
                        else {
                            Nab <- 0
                        }
                        p <- ((2*Naa)+Nab)/(2*Neff)
                        q <- 1 - p
              
                    } 
            
                    else { ## Reproduction, without drift (infinite population) 
                        p <- freq.aa+(freq.ab/2)
                        q <- 1 - p
                    }
                    
                    # UPDATE mean population fitness
                    mean.w <- p * p * Waa + (2*p*q*Wab) + q * q * Wbb
                }
                else { # Fixed at this round
                    if (p <= ZERO){
                        p <- 0.0
                        q <- 1.0
                        mean.w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
                    } 
                    else {
                        p <- 1.0
                        q <- 0.0
                        mean.w <- p*p*Waa+2*p*q*Wab+q*q*Wbb
                    }
                }
    
                ## Add new row
                full.sim <- bind_rows(full.sim,  tibble("population" = j,
                                                        "generation" = i,
                                                        "p" = p,
                                                        "w" = mean.w))                                        
            } #end populations loop
        
            incProgress(1/gen)
        
        } #end generations loop
        
    }) #end progress bar
  
  
    return(full.sim)
}