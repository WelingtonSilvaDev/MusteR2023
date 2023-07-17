# ----------- Interface MusteR 2022 ----------------------
#    __  __           _       _____  
#   |  \/  |         | |     |  __ \ 
#   | \  / |_   _ ___| |_ ___| |__) |
#   | |\/| | | | / __| __/ _ \  _  / 
#   | |  | | |_| \__ \ ||  __/ | \ \ 
#   |_|  |_|\__,_|___/\__\___|_|  \_\
# Created by Welington Goncalves Silva
# Last update on April 8, 2023
# My github: https://github.com/WelingtonSilvaDev

#-----------------------------------------Packages-----------------------------------------------------------------------------------------------------------------------
# Lista de pacotes necessários
 pacotes <- c("bio3d","flexclust", "Matrix","cluster", "stringr","prodlim","stringi","RColorBrewer","rgl", "doMC", "bio3d", "Rpdb","rmarkdown", "pracma","geometry",
                  "deldir", "caret", "graphkernels", "shape", "rARPACK", "mongolite", "tnet",
                  "markdown", "DT", "shiny", "shinythemes", "shinycssloaders", "dplyr","r3dmol",
                  "plotly","bslib", "ragg", "systemfonts", "textshaping", "tidyverse") 
# Verificar se cada pacote está instalado e instalá-lo se necessário
for (pacote in pacotes) {
  if (!require(pacote, character.only = TRUE)) {
    install.packages(pacote, repos = "https://cran.rstudio.com/", dependencies=TRUE)
  }
}
if (!require("Rpdb", character.only = TRUE)) {
   install.packages(pacote, repos = "https://cran.r-project.org/src/contrib/Archive/Rpdb/Rpdb_2.3.tar.gz", dependencies=TRUE)
}
 #-----------------------------------------Sources and Libraries-----------------------------------------------------------------------------------------------------------------------
# source("bib/common-v2.R")
# source("bib/fun-base-v46.R")
# source("bib/fun-mass-v46.R")
# source("main/fun-ligs-v1.R")
source("main/main-ligs.R")
source("show3d/3dmol.R")

library(shiny)
library(bslib)
library(rmarkdown)
library(shinythemes)
library(DT)
library(shinycssloaders)
library(dplyr)
library(r3dmol)
library(plotly)
library(bio3d)
library(stringr)
library(shinyjs)
#-----------------------------------------Spinner options-----------------------------------------------------------------------------------------------------------------------
options(spinner.color = "#0275D8", spinner.color.background = "#ffffff", spinner.size = 2)
 
#----------------------------------------- user interface -----------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(

  theme = bslib::bs_theme(bootswatch = "yeti", base_font = font_google("Montserrat")), # united, yeti,zephyr ou cerulean
  navbarPage(
    title = "MusteR",
    id = "navbarID",
    # imageOutput("home_img",   #isn't being used
    #             width = "100%",
    #             height = "60%",
    # ),
    # tabsetPanel(id="tabset",
    tabPanel(
      title = "Home",
      value = "sweethome",
      splitLayout(
        sidebarPanel(
          width = "100%",
          tags$style(type="text/css", "
           #loadmessage {
             position: fixed;
             top: 590px;
             left: 0px;
             width: 100%;
             padding: 5px 0px 5px 0px;
             text-align: center;
             font-weight: bold;
             font-size: 100%;
             color: #000000;
             background-color: #00FA9A;
             z-index: 105;
           }
          "),
          conditionalPanel(condition="$('html').hasClass('shiny-busy')",
             tags$div("Loading... (it can take some seconds/minutes)",id="loadmessage"),
          ),
          h2(strong("Welcome to MusteR!")),
          h6("Mining System of Contacts in Protein-Ligand Complexes"),   hr(),
          #h4(strong("Fill in the table below separating by comma")), br(),
          radioButtons(inputId = "check_pdb", label = "Download or upload a PDB file",
                       choices = c("Download from PDB site", "Upload local PDB file")),
          
          conditionalPanel(
            condition = "input.check_pdb == 'Download from PDB site'",
            tags$blockquote("You can type in the field below if you want to download directly from RCSB - Protein Data Bank source."),
            textInput("file1", "PDB Code,Ligand Code,Ligand ID", placeholder = "Example:\n6XQU,U5G,401", value = "6XQU,U5G,401"), 
            actionButton("action", "Submit", class = "btn-success"), 
            br(), 
          ),
          
          conditionalPanel(
            condition = "input.check_pdb == 'Upload local PDB file'",
            tags$blockquote("You can upload and type in the fields below the informations about your local upload."), br(),
            fileInput(inputId = "pdb_file", label = "Choose a local PDB file",  accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            ), placeholder = "No file selected"),
            textInput("file2", "PDB Code,Ligand Code,Ligand ID", placeholder = "Info about your upload, example: \n1A2B,LIG,401"),
            actionButton("pdb_archive", "Submit upload", class = "btn-success")
          ),
        #),
         
        ),
        
      ), 
    ),
   
        tabPanel(
          title = "Interactive Map",
          #id = "tabmaps",
          value = "Mapid",
          sidebarLayout(position = "left", 
              sidebarPanel(#h5(strong("Molecular Visualization")),                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            verbatimTextOutput("titlepdb"),
                           textOutput("pdbtitle"),
                           # verbatimTextOutput("pdbtitle"),
                           #textOutput("pdbinfo"),
                           width = 12
                           ), 
          mainPanel(width=12,
                fluidRow(
                  column(6,
                         withSpinner(r3dmolOutput("graph1", width = "auto", height = "500px"), type = 3)
                  ),  
                  column(6,
                         withSpinner(plotlyOutput("graph2", width = "auto", height = "500px"), type = 3)
                  ),
                  column(5,
                         sliderInput("slider", "Choose the number of the largest residues areas", #"Choose the number of the largest residues areas",
                                     min = 2, max = 7, value = 7, step = 1)
                  ),
                  column(3,
                         radioButtons("molecule_style", label = "Molecular Visualization:",
                                      choices = list("Sticks" = "sticks", "Spheres" = "spheres"),
                                      selected = "sticks"
                         )
                  ),
                  column(3,
                         downloadButton("downloadBtn", "Download Tags Map Csv File")
                  ),

            ),
   
           )
      ),

    ),
   
      tabPanel(
        title = "About",
        # Markdown sobre o MusteR
        mainPanel(
          includeMarkdown("md/about.md"),
        ),
       
      ),
    inverse = T

  ),

)


#----------------------------------------------------------------------------------Funcao Server-----------------------------------------------------------------------------------

server <- function(input, output, session) {
  
#----------------------------------------- atualizar visualizacao R3DMOL a cada vez que o slider atualizar  ----------------------------------------------------------------  

    observeEvent(input$molecule_style, {
      output$graph1 <- renderR3dmol({
        mol3d(paste0(lig.table$pdb_title, ".pdb"), lig.table$lig1id, check_style = input$molecule_style)
      })
    })

observeEvent(input$action, {
  observeEvent(input$slider, {
      nn <- input$slider
      mol_name = list.files(path = "pdb-mol/", pattern = ".pdb", all.files = T)
      # browser()
      area <- m.ligs1 %>% filter(resn != lig.table$lig1n) %>% arrange(desc(area)) %>% slice_head(n=nn) %>% 
        select(area) %>% slice_tail() %>% as.numeric()
      area = area - 0.1
      x = write_molecule_muster_version2(pdb1[[1]], dba2[[1]][[2]], area.lim = area)
      output$graph1 <- renderR3dmol({
         mol3d(paste0(lig.table$pdb_title, ".pdb"), lig.table$lig1id, check_style = input$molecule_style)
      })
      output$graph2 = renderPlotly({
        p.ligs1 = plot_lig_res_interactions(m.ligs1,scene="scene1",area.lim=area,alpha=0.1,n=1200)
        p.ligs1
      })
    })
})
observeEvent(input$pdb_archive, {
  observeEvent(input$slider, {
    nn <- input$slider
    mol_name = list.files(path = "pdb-mol/", pattern = ".pdb", all.files = T)
    
    area <- m.ligs1 %>% filter(resn != lig.table$lig1n) %>% arrange(desc(area)) %>% slice_head(n=nn) %>% 
      select(area) %>% slice_tail() %>% as.numeric()
    area = area - 0.1
    x = write_molecule_muster_version2(pdb1[[1]], dba2[[1]][[2]], area.lim = area)
    output$graph1 <- renderR3dmol({
      mol3d(paste0(lig.table$pdb_title, ".pdb"), lig.table$lig1id, check_style = input$molecule_style)
    })
    output$graph2 = renderPlotly({
      p.ligs1 = plot_lig_res_interactions(m.ligs1,scene="scene1",area.lim=area,alpha=0.1,n=1200)
      p.ligs1
    })
  })
})
# }})
write_molecule_muster_version2 = function(pdb, dba, area.lim = 500,
                                           ids = c(1,1), force = T,
                                           title = T, pdb_out = T, 
                                           pdbpath = "pdb-mol", 
                                           pdbpathout = pdbpath){

    if (!dir.exists(pdbpathout)){
      dir.create(pdbpathout)
    }
    
    res.tab = dba$exp$residues[[ids[1]]][[ids[2]]]$res.table
    res.tab = res.tab %>% filter((area > area.lim & chain == "A") |
                                   (chain == "B"))
    
    
    aux1 = atom.select(pdb, resid = res.tab$resn)
    aux2 = atom.select(pdb, resno = as.numeric(res.tab$resi))
    aux12 = combine.select(aux1, aux2, operator="AND")
    
    pdb.out = pdb
    pdb.out$atom = pdb$atom[aux12$atom, ]
    pdb.out$xyz = pdb.out$xyz[aux12$xyz]
    pdb.out$calpha = pdb.out$atom$elety == "CA"
    #browser()
    res.out = write_pdb(pdb.out,
                        pdbpath = pdbpath,
                        force = force,
                        title = title,
                        pdb_out = pdb_out,
                        pdbpathout = pdbpathout)
    
    return(res.out)
  }
#----------------------------------------- Verifica se eh pra usar LOCAL = 0 ou 1  ----------------------------------------------------------------
  observeEvent(input$check_pdb, {
    if(input$check_pdb == "Upload local PDB file") {
      check_pdb <- 1
    } else {
      check_pdb <- 0
    }
  })
#----------------------------------------- Realizar o download do tagmaps  -------------------------------------------------------------------------------------------------------- 
  output$downloadBtn <- downloadHandler(
    filename = function() {
      # "PDB Name,Ligand Name,Ligand ID"
      filename_ = paste0(lig.table$pdb_title,"_", lig.table$lig1n,"_",lig.table$lig1id, ".csv")
    },
    content = function(file) {
      # Aqui é o caminho para o arquivo CSV que deseja baixar
      file_path <- "tags_out/tags.csv"
      # Lê o arquivo CSV e escreve no arquivo que será baixado
      data <- read.csv(file_path)
      write.csv(data, file, row.names = F)
      headers = function() {
        addHeader("Content-Disposition", "attachment; filename=filename_")
      }
    },

  )
#----------------------------------------- metodo para pegar PDB Name, Ligand Name e Ligand ID ------------------------------------------------------------------------------------

observeEvent(input$action,
                          {
                            pdb_data <- reactive({
                              read.pdb(paste0("pdb/", lig.table$pdb_title, ".pdb"))
                            })

                            pdb_title <- reactive({
                              pdb_data()$title # nao estava aparecendo o titulo do pdb por esse metodo
                            })
                            output$pdbtitle <- renderPrint({
                              paste0("PDB Name: ", lig.table$pdb_title, "; Ligand Name: ", lig.table$lig1n, "; Ligand ID:  ", lig.table$lig1id)
                              #browser()
                            })
})

observeEvent(input$pdb_archive,
               {
               pdb_data <- reactive({
                 read.pdb(paste0("pdb/", lig.table$pdb_title, ".pdb"))
               })

               pdb_title <- reactive({
                 pdb_data()$title # nao estava aparecendo o titulo do pdb por esse metodo
               })
               output$pdbtitle <- renderPrint({
                 paste0("PDB Name: ", lig.table$pdb_title, "; Ligand Name: ", lig.table$lig1n, "; Ligand ID:  ", lig.table$lig1id)
                 #browser()
               })
})
  
#----------------------------------------- atualiza a visualizacao para pagina de Interactive Map------------------------------------------------------------------------------------
  observeEvent(input$action, {
    updateTabsetPanel(session, "navbarID",
                      selected = "Mapid")
  })
  observeEvent(input$pdb_archive, {
    updateTabsetPanel(session, "navbarID",
                      selected ="Mapid")
  })
  
#----------------------------------------- pega arquivo de download e info------------------------------------------------------------------------------------
  observeEvent(input$action,{
    
    file_df <- read.table(text = gsub(" ","",input$file1), sep = ",", header = FALSE, fileEncoding = "UTF-8")
    colnames(file_df) <- c("pdb_title", "lig1n", "lig1id")
    file_df <- data.frame(lapply(file_df, function(v) {
      if (is.character(v)) {
        return(toupper(v))
      } else {
        return(v)
      }
    }))
    write.csv(file_df, "file_df.csv", row.names = FALSE)
    check_pdb = 0
    backend_MusteR("file_df.csv", check_pdb)
    write_csv(m.ligs1, file = "tags_out/tags.csv")
  })
#----------------------------------------- pega arquivo de upload e info----------------------------------------------------------------------------------------------------
  observeEvent(input$pdb_archive, {
    file_df <- read.table(text = gsub(" ","",input$file2), sep = ",", header = FALSE, fileEncoding = "UTF-8")
    colnames(file_df) <- c("pdb_title", "lig1n", "lig1id")
    file_df <- data.frame(lapply(file_df, function(v) {
      if (is.character(v)) {
        return(toupper(v))
      } else {
        return(v)
      }
    }))
    check_pdb <- 1
    write.csv(file_df, "file_df.csv", row.names = FALSE)
    backend_MusteR("file_df.csv", check_pdb)
    write_csv(m.ligs1, file = "tags_out/tags.csv")
  })
}

shinyApp(ui = ui, server = server)
