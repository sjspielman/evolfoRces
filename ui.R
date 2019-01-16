library(shiny)
library(tidyverse)
library(colourpicker)
library(shinythemes)
library(plotly)
library(cowplot)


ui <- shinyUI(navbarPage(theme = shinytheme("sandstone"), "Evolutionary Forces Simulation",
    tabPanel("Single population",
        sidebarPanel(
            sliderInput("p","Starting allele frequency A",value=0.5,min=0,max=1),
            sliderInput("Waa_s","Fitness (W) of genotype AA",value=1,min=0,max=1),
            sliderInput("Wab_s","Fitness (W) of genotype Aa",value=0.95,min=0,max=1),
            sliderInput("Wbb_s","Fitness (W) of genotype aa",value=0.90,min=0,max=1),
            sliderInput("Uab","Mutation Rate from A to a",value=0,min=0,max=0.5),
            sliderInput("Uba","Mutation Rate from a to A",value=0,min=0,max=0.5),
            numericInput("gen_s","Number of Generations (maximum 1000)",100,min=1,max=1000),
            radioButtons("usedrift", "Genetic drift", 
                          choices = c("Turn off genetic drift (INFINITE population size)" = "infinite",
                                      "Turn on genetic drift (FINITE population size)" = "drift"), 
                          selected = "infinite"
                        ),
            conditionalPanel(condition = "input.usedrift == 'drift'",
                {textInput("N","Population Size (maximum 100000)",value=100)}
            ),
            conditionalPanel(condition = "input.usedrift == 'drift'",
                {numericInput("nRep","Number of Replicate Populations (maximum 50)",value = 10)}
            ),
            br(),
            actionButton("go_s","Run Simulation!",width="100%")
        ),
        mainPanel(

            br(), 
            div(style = "float:center; font-size:18px",
                tags$b(textOutput("result_header_s"))
            ),
            br(),
            div(style = "float:left; font-size:16px",
                tags$b(tableOutput("result_table_s"))
            ),
            br(),br(),br(),br(),br(),
            div(style = "display:inline-block; height: 300px; width=300px; ",
                plotlyOutput("singleplot.frequency_s", height = "100%", width = "100%"),
                br(),br(),
                plotlyOutput("singleplot.fitness_s", height = "100%", width = "100%")
            ),
            br(),br(),
            div(style = "float:left",
                uiOutput("downloaddata_s")
            )
        )
    ),
    tabPanel("Migration",
        sidebarPanel(sliderInput("p.main","CONTINENT fixed allele A frequency",value=0.5,min=0,max=1),
            sliderInput("p.island","ISLAND starting allele A frequency",value=0.5,min=0,max=1),
            sliderInput("m","Migration rate, from continent to island. (Fraction of island individuals who are migrants, per generation).",value=0.1,min=0,max=1),
            sliderInput("Waa_m","Fitness (W) of genotype AA, on ISLAND",value=1,min=0,max=1),
            sliderInput("Wab_m","Fitness (W) of genotype Aa, on ISLAND",value=0.95,min=0,max=1),
            sliderInput("Wbb_m","Fitness (W) of genotype aa, on ISLAND",value=0.90,min=0,max=1),
            numericInput("gen_m","Number of Generations",100,min=1,max=1000),
            colourInput("yaycolor", "Select color:", value = "seagreen"),
            br(),br(),
            actionButton("go_m","Run Simulation!",width="100%")
        ),
       mainPanel(

            br(), 
            div(style = "float:center; font-size:18px",
                tags$b(textOutput("result_header_m"))
            ),
            br(),
            div(style = "float:left; font-size:16px",
                tags$b(tableOutput("result_table_m"))
            ),
            br(),br(),br(),br(),br(),
            div(style = "display:inline-block; height: 325px; width=325px; ",
                plotlyOutput("singleplot.frequency_m", height = "100%", width = "100%"),
                br(),br(),
                plotlyOutput("singleplot.fitness_m", height = "100%", width = "100%")
            ),
            br(),br(),
            div(style = "float:left",
                uiOutput("downloaddata_m")
            )
        )
    ),
    tabPanel("About and Help",            
        div(style = "width:80%;  font-size:16px; float:center; text-align:left;",
            p("This Shiny application allows users to visualize how evolutionary forces affect i) allele frequencies and b) population fitness over time. Source code and licensing for this app is available from", tags$a(href = "https://github.com/spielmanlab/evolforRces", "https://github.com/spielmanlab/evolforRces"),".This application is primarily intended for classroom use."),
            br(),
            p("Written and maintained by ", tags$a(href = "http://spielmanlab.io", "Stephanie J. Spielman, PhD"), ", with many thanks to CJ Battey's excellent", tags$a(href = "https://github.com/cjbattey/driftR", "driftr Shiny app"), " for inspiration."),
            h3(tags$b("Tab One: Single Population")),
            hr(),
            p("The ", tags$b("Single Population"), "tab allows you to simulate allele frequencies over time for a single population with two alleles (A/a,) with NO migration. Evolutionary forces you can specify for this tab include the following:"),
        tags$ul(
            tags$li(tags$b("Natural selection"), "by setting the fitness for each genotype"), 
            tags$li(tags$b("Mutation rate"), "both forward (to allele a) and backwards (to allele A)"), 
            tags$li(tags$b("Genetic drift"))
        ),
        br(),
        p("You can also specify the following simulation settings:"),
        tags$ul(
            tags$li(tags$b("Starting allele 'A' frequency."), "This represents the proportion of alleles in the population that are 'A' at time 0, i.e. before the simulation begins."), 
            tags$li(tags$b("Number of generations"), "telling the program how many generations to run the simulation for"), 
            tags$li(tags$b("Population size"), "if drift is turned OFF"),
            tags$li(tags$b("Number of simulation replicates."), "For any simulation without drift, only one replicate will be run.")
        ),
     
        
        h3(tags$b("Tab Two: Migration Population")),
            hr(),
            p("The", tags$b("Migration"), "tab allows you to simulate allele frequencies over time for an Island-Continent Model of two populations, with two alleles (A/a).  Evolutionary forces you can specify for this tab include the following:"),
        tags$ul(
            tags$li(tags$b("Natural selection"), "by setting the fitness for each genotype, for each population (island and continent)"), 
            tags$li(tags$b("Migration rate"), "the proportion of island individuals arriving from continent each generation"),
            tags$li(tags$b("There is NO GENETIC DRIFT and NO MUTATION in migration simulations."))
        ),
            p("You can also specify the following simulation settings:"),
            tags$ul(
                tags$li(tags$b("Starting allele 'A' frequency for each population (island and continent)."), "This represents the proportion of alleles in the population that are 'A' at time 0, i.e. before the simulation begins. The allele frequency for the island will change with simulations, but the continent will not change."), 
                tags$li(tags$b("Number of generations"), "for the simulation to run for"),
                tags$li(tags$b("Choose your own plot color!"))
            ) ,
        br(),br(),br()       
        ) 
    
    )
))



  


