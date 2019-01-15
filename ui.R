library(shiny)
library(tidyverse)
library(colourpicker)
library(shinythemes)

# shinyUI(navbarPage("App Title",
#   tabPanel("Tab Name",
#     sidebarPanel([inputs for the first tab]),
#     mainPanel([outputs for the first tab])
#   )


ui <- shinyUI(navbarPage(theme = shinytheme("sandstone"), "Evolutionary Forces Simulation",
    tabPanel("About and Help",
    
        div(style = "width:80%;  font-size:16px; float:center; text-align:left;",
            p("This Shiny application, inspired by CJ Battey's",
            tags$a(href = "https://github.com/cjbattey/driftR", "driftr Shiny app"), "allows users to visualize how evolutionary forces affect i) allele frequencies and b) population fitness over time. Source code and licensing for this app is available from", tags$a(href = "https://github.com/spielmanlab/driftr-lite", "https://github.com/spielmanlab/driftr-lite"),"."),
            br(),
            h3(tags$b("Tab One: Single Population")),
            br(),
            p("The ", tags$b("Single Population"), "tab allows you to simulate allele frequencies over time for a single population with two alleles (A/B,) with NO migration. Evolutionary forces you can specify for this tab include the following:"),
        tags$ul(
            tags$li(tags$b("Natural selection"), "by setting the fitness for each genotype"), 
            tags$li(tags$b("Mutation rate"), "both forward (to allele B) and backwards (to allele A)"), 
            tags$li(tags$b("Genetic drift"), "by setting the population size (or clicking the Infinite Population setting to turn off drift")
        ),
        br(),
        p("You can also specify the following simulation settings:"),
        tags$ul(
            tags$li(tags$b("Starting allele 'A' frequency."), "This represents the proportion of alleles in the population that are 'A' at time 0, i.e. before the simulation begins."), 
            tags$li(tags$b("Number of generations"), "for the simulation to run for"), 
            tags$li(tags$b("Number of simulation replicates."), "For any simulation without drift, only one replicate will be run.")
        ),
         
        
        h3(tags$b("Tab Two: Migration Population")),
            br(),
            p("The", tags$b("Migration"), "tab allows you to simulate allele frequencies over time for an Island-Continent Model of two populations, with two alleles (A/B).  Evolutionary forces you can specify for this tab include the following:"),
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
            )        
       
              
        )
    ),
    tabPanel("Single population",
        sidebarPanel(
            sliderInput("p","Starting allele frequency A",value=0.5,min=0,max=1),
            sliderInput("Waa_s","Relative fitness of genotype AA",value=1,min=0,max=1),
            sliderInput("Wab_s","Relative fitness of genotype AB",value=0.95,min=0,max=1),
            sliderInput("Wbb_s","Relative fitness of genotype BB",value=0.90,min=0,max=1),
            sliderInput("Uab","Mutation Rate from A to B",value=0,min=0,max=0.5),
            sliderInput("Uba","Mutation Rate from B to A",value=0,min=0,max=0.5),
            numericInput("gen_s","Number of Generations",100,min=1,max=5000),
            textInput("Neff","Population Size",value=100),
            checkboxInput("infinitePop","Infinite Population (no drift). Checking this box will ignore population size above.",value = F),
            numericInput("nRep","Number of Replicate Populations",10,min=1,max=100),
            br(),
            actionButton("go_s","Run Simulation!",width="100%")
        ),
        mainPanel(
            div(style = "height: 700px; width=600px; font-size:18px; float:center; text-align:center;",
                plotOutput("singleplot", height = "90%", width = "75%"),
                tags$b(textOutput("fixation_s"))
            ),
            div(style = "float:left; font-size:16px", 
                tableOutput("fixlosstable_s")
            ),
            div(style = "float:right",
                uiOutput("downloadp_s"),
                br(), br(),
                uiOutput("downloaddata_s")
            )
        )
    ),
    tabPanel("Migration",
        sidebarPanel(sliderInput("p.main","CONTINENT fixed allele A frequency",value=0.5,min=0,max=1),
            sliderInput("p.island","ISLAND starting allele A frequency",value=0.5,min=0,max=1),
            sliderInput("m","Migration rate, from continent to island. (Fraction of island individuals who are migrants, per generation).",value=0.1,min=0,max=1),
            sliderInput("Waa_m","Relative fitness of genotype AA, on ISLAND",value=1,min=0,max=1),
            sliderInput("Wab_m","Relative fitness of genotype AB, on ISLAND",value=0.95,min=0,max=1),
            sliderInput("Wbb_m","Relative fitness of genotype BB, on ISLAND",value=0.90,min=0,max=1),
            numericInput("gen_m","Number of Generations",100,min=1,max=5000),
            colourInput("yaycolor", "Select color:", value = "seagreen"),
            br(),br(),
            actionButton("go_m","Run Simulation!",width="100%")
        ),
        mainPanel(
            div(style = "height: 700px; width=600px; font-size:18px; float:center; text-align:left;",
                plotOutput("migrationplot", height = "90%", width = "75%"),
                tags$b(textOutput("fixation_m"))
            ),
            div(style = "float:left",
                tableOutput("fixlosstable_m")
            ),
            div(style = "float:right",
                uiOutput("downloadp_m"),
                br(),
                uiOutput("downloaddata_m")
            )
        )
    )
))



  


