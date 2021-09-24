# SpliceR UI
options(shiny.sanitize.errors = FALSE)
version = "1.2.0"
update = "May 24th, 2021"

shinyUI(
  pageWithSidebar(
    
    # App title ----
    headerPanel(paste0("SpliceR v", version)),
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      p(""),
      p(""),
      em(paste0("Updated ", update, ".")),
      em(paste0("Please email klues009@umn.edu with any questions, comments, or concerns.")),
      p(""),
      tags$a(href="https://www.biorxiv.org/content/biorxiv/early/2019/05/09/633685.full.pdf", em("Check out our Nature Communications article")),
      # actionButton(inputId = "actionButton", label = "Update"),
      p(""),
      checkboxInput(inputId = "use_example", label = "See example.", value = FALSE, width = NULL),
      p(""),
      p("Please enter an Ensembl ID from your transcript of interest from www.ensembl.org"),
      textInput(inputId = 'id',
                label = 'Enter Ensembl id'),
      selectInput(inputId = 'pam',
                  label = "Choose PAM",
                  choices = c("Choose a PAM" = "",
                                "NGG - SpCas9" = "NGG",
                                "NGA - SpCas9-VQR" = "NGA",
                                "NG - xCas9 / Cas9-NG / SpG Cas9" = "NGN",
                                "N[R|Y] - SpRY Cas9" = "NNN", 
                                "NGCG - SpCas9-VRER" = "NGCG",
                                "NNNRRT - SaCas9-KKH" = "NNNRRT",
                                "NRCH — SpCas9-NRCH" = "NRCH",
                                "NRRH — SpCas9-NRRH" = "NRRH",
                                "NRTH — SpCas9-NRTH" = "NRTH"
                              ),
                  selected = "NGG - Sp Cas9",
                  multiple = FALSE),
      selectInput(inputId = 'species',
                  label = "Choose species",
                  choices = c("Choose a species" = "",
                              "Human" = "Homo_sapiens",
                              "Mouse" = "Mus_musculus",
                              "Zebrafish" = "Danio_rerio",
                              "Pig" = "Sus_scrofa",
                              "Cow" = "Bos_taurus",
                              "Mouse lemur" = "Microcebus_murinus",
                              "Rat" = "Rattus_norvegicus",
                              "Rabbit" = "Oryctolagus_cuniculus",
                              "Cat" = "Felis_catus",
                              "Dog" = "Canis_familiaris",
                              #"Chimpanzee",
                              "Rhesus Macaque" = "Macaca_mulatta",
                              "C. elegans" = "Caenorhabditis_elegans",
                              "D. melanogaster" = "Drosophila_melanogaster",
                              "Yeast" = "Saccharomyces_cerevisiae"
                              #"Hamster",
                              # Also include plants via plants.ensembl.org
                              ),
                  selected = "Human",
                  multiple = FALSE),
      selectInput(inputId = 'enzyme',
                  label = "Choose class of base editor",
                  choices = c("Choose base editors class" = "",
                              "CBE",
                              "ABE",
                              "CBE and ABE",
                              "CBE and-or ABE"
                              ),
                  selected = "CBE"
                  ),
      selectInput(inputId = 'splice_site',
                  label = "Choose target splice-sites",
                  choices = c("Choose splice-sites to target" = "",
                              "splice-donors",
                              "splice-acceptors",
                              "splice-donors and splice-acceptors"),
                  selected = "splice-donors and splice-acceptors"
      ),
      # drop-down for different deaminase options
      # 
      
      sliderInput(inputId = "editing_window",
                  label = "Choose base editing window",
                  min = 1,
                  max = 20,
                  value = c(4,9),
                  step = 1),
      checkboxInput(inputId = "strictFilter", 
                    label = p("Only include splice-donor guides within editing window for CBE ", strong(em("and")), "ABE.")
                    ),
      downloadButton("downloadAllData", "Download all guides"),
      p(),
      downloadButton("downloadSelectedData", "Download selected guides")
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Instructions",
                           includeMarkdown("instructions.md")
                  ),
                  tabPanel("Predicted guides",
                           
      fluidRow(
        column(12,
               h2("Predicted guides"),
               br(),
               DT::DTOutput("df"),
               br(),
               h2("Parameters used"),
               verbatimTextOutput("print"),
        )
      )
    ),
    
    tabPanel("Ensembl gene map",
             fixedPage(
               br(),
               htmlOutput("ensembl")
             )
    )# ,
    # This function is under development
    # tabPanel("BE-Hive predictions",
    #          fixedPage(
    #            br(),
    #            htmlOutput("behive")
    #          )
    # )
  )
)
))
