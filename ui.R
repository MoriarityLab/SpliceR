# EditR UI
version = "1.0.7"
update = "April 2nd, 2020"

shinyUI(
  pageWithSidebar(
    
    # App title ----
    headerPanel(paste0("SpliceR v", version)),
    
    # Sidebar panel for inputs ----
    sidebarPanel(
      # tags$img(href = "https://www.cancer.umn.edu/bio/osteosarcoma-staff/branden-moriarity", src="https://upload.wikimedia.org/wikipedia/commons/thumb/f/f9/Minnesota_Golden_Gophers_logo.svg/1200px-Minnesota_Golden_Gophers_logo.svg.png", alt="SpliceR", width="100%", class="unframed", align="center"),
      p(""),
      p(""),
      em(paste0("Please note this application is under development.")),
      em(paste0("Updated ", update, ".")),
      em(paste0("Please email klues009@umn.edu with any questions, comments, or concerns.")),
      p(""),
      checkboxInput(inputId = "use_example", label = "See example.", value = FALSE, width = NULL),
      p(""),
      p("Please enter an Ensembl ID from your transcript of interest from www.ensembl.org"),
      textInput(inputId = 'id',
                label = 'Enter Ensembl id'),
      selectInput(inputId = 'pam',
                  label = "Choose PAM",
                  choices = c("Choose a PAM" = "",
                                "NGN - xCas9 / Cas9-NG" = "NGN",
                                "NGG - SpCas9" = "NGG",
                                "NGA - SpCas9-VQR" = "NGA"#,
                                # "NGCG - SpCas9-VRER" = "NGCG",
                                #"NNNRRT - SaCas9-KKH" = "NNNRRT"
                              ),
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
                              #"Rhesus",
                              "C. elegans" = "Caenorhabditis_elegans",
                              "D. melanogaster" = "Drosophila_melanogaster",
                              "Yeast" = "Saccharomyces_cerevisiae"
                              #"Hamster",
                              # Also include plants via plants.ensembl.org
                              ),
                  multiple = FALSE),
      selectInput(inputId = 'enzyme',
                  label = "Choose class of base editor",
                  choices = c("Choose a class of base editor" = "",
                              "CBE",
                              "ABE",
                              "CBE or ABE"),
                  selected = "CBE"
                  ),
      selectInput(inputId = 'splice_site',
                  label = "Choose target splice-sites",
                  choices = c("Choose splice-sites to target" = "",
                              "splice-donor",
                              "splice-acceptor",
                              "splice-donors and splice-acceptors"),
                  selected = "splice-donors and splice-acceptors"
      ),
      p(),
      downloadButton("downloadData", "Download")
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
               br(),
               h2("Parameters used"),
               verbatimTextOutput("print")
        )
      )
    ),
    
    tabPanel("Ensembl gene map",
             fixedPage(
               br(),
               htmlOutput("frame")
             )
    )
  )
)
))
