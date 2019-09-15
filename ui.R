library(shiny)
library(shinythemes)
library(shinyWidgets)
library(tidyverse)
library(colourpicker)
library(DT)


#################################################################################################



ui <- shinyUI(navbarPage(theme = shinytheme("flatly"), "evolfoRces: Two-allele population genetics simulations",
    tabPanel("Single population",
        sidebarPanel(width=3,
            sliderInput("p","Starting allele frequency A",value=0.5,min=0,max=1),
            sliderInput("Waa_s","Fitness (W) of genotype AA",value=1,min=0,max=1),
            sliderInput("Wab_s","Fitness (W) of genotype Aa",value=0.95,min=0,max=1),
            sliderInput("Wbb_s","Fitness (W) of genotype aa",value=0.90,min=0,max=1),
            sliderInput("Uab","Mutation Rate from A to a",value=0,min=0,max=0.1),
            sliderInput("Uba","Mutation Rate from a to A",value=0,min=0,max=0.1),
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
            
            
            conditionalPanel(condition = "input.usedrift == 'infinite'",   
                {colourpicker::colourInput("line_color_s", "Color:", value = "blue")}
            ),               
            br(), 
            actionButton("go_s","Run Simulation!",width="100%")
        ),
        mainPanel(width=9,
            h3("Single population simulation results"),
            br(),
            div(style = "float:center; font-size:12px",
                DTOutput("single_table")
            ),
            br(),
            div(style = "display:inline-block;",
                plotOutput("single_plot", height = "300px", width = "1000px"),
                br(),
                div(style = "float:right;",
                    actionButton(inputId = "store_single_btn", label = "Click to store results.")
                )
            ),
            br(), 
            tags$hr(),

            h3("Stored simulation"),
            tags$b(uiOutput("single_name_stored")),
            br(),
            div(style = "float:center; font-size:12px",
                DTOutput("single_table_stored")
            ),
            br(),br(),br(),
            div(style = "display:inline-block;",
                plotOutput("single_plot_stored", height = "300px", width = "1000px"),
                br(),
                div(style = "float:right;",
                    actionButton(inputId = "clear_single_btn", label = "Click to clear stored results.")
                )
            )
        )

    ),
    tabPanel("Migration",
        sidebarPanel(width = 3, 
            sliderInput("p.main","CONTINENT fixed allele A frequency",value=0.5,min=0,max=1),
            sliderInput("p.island","ISLAND starting allele A frequency",value=0.5,min=0,max=1),
            sliderInput("m","Migration rate to island",value=0.05,min=0,max=0.5),
            sliderInput("Waa_m","Fitness (W) of genotype AA, on ISLAND",value=1,min=0,max=1),
            sliderInput("Wab_m","Fitness (W) of genotype Aa, on ISLAND",value=0.95,min=0,max=1),
            sliderInput("Wbb_m","Fitness (W) of genotype aa, on ISLAND",value=0.90,min=0,max=1),
            numericInput("gen_m","Number of Generations (maximum 1000)",100,min=1,max=1000),
            colourInput("line_color_m", "Color:", value = "blue"),
            br(),br(),
            actionButton("go_m","Run Simulation!",width="100%")
        ),
       mainPanel(width = 9, 
            h3("Island-continent simulation results"),
            br(),
            div(style = "float:center; font-size:12px",
                DTOutput("migration_table")
            ),
            br(),
            div(style = "display:inline-block;",
                plotOutput("migration_plot", height = "300px", width = "1000px"),
                br(),
                div(style = "float:right;",
                    actionButton(inputId = "store_migration_btn", label = "Click to store results.")
                )
            ),
            br(), 
            tags$hr(),

            h3("Stored simulation"),
            tags$b(uiOutput("migration_name_stored")),
            br(),
            div(style = "float:center; font-size:12px",
                DTOutput("migration_table_stored")
            ),
            br(),br(),br(),
            div(style = "display:inline-block;",
                plotOutput("migration_plot_stored", height = "300px", width = "1000px"),
                br(),
                div(style = "float:right;",
                    actionButton(inputId = "clear_migration_btn", label = "Click to clear stored results.")
                )
            )
        )

    ),
    tabPanel("About and Help",            
        div(style = "width:80%;  font-size:16px; float:center; text-align:left;",
            p("This Shiny application allows users to visualize how evolutionary forces affect i) allele frequencies and b) population fitness over time. Source code and licensing for this app is available from", tags$a(href = "https://github.com/spielmanlab/evolforRces", "https://github.com/spielmanlab/evolforRces"),".This application is primarily intended for classroom use."),
            br(),
            p("Written and maintained by ", tags$a(href = "http://spielmanlab.io", "Stephanie J. Spielman, PhD"), ", with many thanks to CJ Battey's excellent", tags$a(href = "https://github.com/cjbattey/driftR", "driftr Shiny app"), " for inspiration."),
            br(),br(),
            
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
                tags$li(tags$b("Number of generations"), "for which the simulation will run"),
                tags$li(tags$b("Population size"), "if drift is turned OFF"),
                tags$li(tags$b("Number of simulation replicates."), "For any simulation without drift, only one replicate will be run."),
                tags$li(tags$b("Choose your own plot color!"))
            ),
     
            br(), br(),
            h3(tags$b("Tab Two: Migration Population")),
            hr(),
            p("The", tags$b("Migration"), "tab allows you to simulate allele frequencies over time for an Island-Continent Model of two populations, with two alleles (A/a).  Evolutionary forces you can specify for this tab include the following:"),
            tags$ul(
                tags$li(tags$b("Natural selection"), "by setting the fitness for each genotype, for each population (island and continent)"), 
                tags$li(tags$b("Migration rate"), "the proportion of island individuals arriving from continent each generation"),
                tags$li(tags$b("There is NO GENETIC DRIFT and NO MUTATION in migration simulations."))
            ),
            br(),
            p("You can also specify the following simulation settings:"),
            tags$ul(
                tags$li(tags$b("Starting allele 'A' frequency for each population (island and continent)."), "This represents the proportion of alleles in the population that are 'A' at time 0, i.e. before the simulation begins. The allele frequency for the island will change with simulations, but the continent will not change."), 
                tags$li(tags$b("Number of generations"), "for which the simulation will run"),
                tags$li(tags$b("Choose your own plot color!"))
            ),
            br(),br(),br()       
        ) 
    
    )
))



  


