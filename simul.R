# SJS modification of original `dev.R` from driftR package
library(tidyverse)
library(viridis) 
library(plotly)

library(shinythemes)
ZERO <- 1e-7


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
