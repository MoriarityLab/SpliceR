options(shiny.sanitize.errors = FALSE)
# Example data
example_id <- "ENST00000380152.8"
example_pam <- "NGG"
example_species <- "Homo_sapiens"
example_enzyme = "CBE or ABE"
example_splice_site = "splice-donors and splice-acceptors"




### To do:
# Update with note to contact about other species.
# Add flank sequences to be used with BE-Hive
# Add notes about how the scoring works
# Add BE-Hive page
# Add figures to explain the algorithm
# In general more verbose explanations of the application and theory
# Have a figure for showing where guides are mapping to
# Have the ability to choose a guide set?
# Have suggested for how we work up experiments
# No guides available? Try these other resources
# Have an adjustable base editing window
# Change Aneesha's name on the application
# Error handling...
#
## Reaches
# Error handling

shinyServer(
  function(input, output) {
    
    id.Reactive = reactive({
      if(input$id != "") {
        input$id
      } else {
        if(input$use_example) {
          example_id
        } else {
          validate(
            need(input$id != "", "Please enter an ensembl id")
          )
          return(input$id)
        }
      }
    })
    
    pam.Reactive = reactive({
      if(input$pam != "") {
        input$pam
      } else {
        if(input$use_example) {
          example_pam
        } else {
          validate(
            need(input$pam != "", "Please enter a PAM identity")
          )
          return(input$pam)
        }
      }
    })
    
    species.Reactive = reactive({
      if(input$species != "") {
        input$species
      } else {
        if(input$use_example) {
          example_species
        } else {
          validate(
            need(input$species != "", "Please choose a species")
          )
          return(input$species)
        }
      }
    })
    
    enzymeClass.Reactive = reactive({
      if(input$enzyme != "") {
        input$enzyme
      } else {
        if(input$use_example) {
          example_enzyme
        } else {
          validate(
            need(input$enzyme != "", "Please choose a base editor(s)")
          )
          return(input$enzyme)
        }
      }
    })
    
    splice_site.Reactive = reactive({
      if(input$splice_site != "") {
          input$splice_site
      } else {
        if(input$use_example) {
          example_splice_site
        } else {
          validate(
            need(input$splice_site != "", "Please enter a target splice-site(s)")
          )
          return(input$splice_site)
        }
      }
    })
    
    reactive({ validate(
      need(input$splice_site == "splice-acceptor" & input$enzyme == "CBE or ABE", "Please enter a target splice-site(s)")
    )})
    
    # console for debugging
    output$console = renderPrint({
      id.Reactive()
    })
    
    
    # Render table with gRNAs
    output.df = reactive({
      
      runSpliceR(ensembl_transcript_id = id.Reactive(), pam = pam.Reactive(), species = species.Reactive() #,
                 # min_editing_window = input$editing_window[1],
                 # max_editing_window = input$editing_window[2]
                 )[[1]] %>%
        
        filterGuides(runSpliceR.React = .,
                     enzymeClass.React = enzymeClass.Reactive(),
                     min_editing_window = input$editing_window[1],
                     max_editing_window = input$editing_window[2],
                     splice_site.React = splice_site.Reactive(),
                     strictFilter = input$strictFilter
                     )
        
      #   # IF CBE or ABE is selected, THEN do not filter, else
      #   {if(enzymeClass.Reactive() == "CBE and-or ABE") {
      #     {.} %>%
      #       filter(.,
      #       (abe_position >= input$editing_window[1] | cbe_position >= input$editing_window[1]) &
      #         (abe_position <= input$editing_window[2] | cbe_position <= input$editing_window[2])
      #     )
      #   } else {
      #     
      #     # IF CBE and ABE is selected, THEN include only CBE and ABE
      #     if(enzymeClass.Reactive() == "CBE and ABE") {
      #       filter(., enzyme == "CBE and ABE") %>%
      #         filter(
      #           (abe_position >= input$editing_window[1]| cbe_position >= input$editing_window[1]) &
      #             (abe_position <= input$editing_window[2] | cbe_position <= input$editing_window[2])
      #         )
      #     } else {
      #       
      #       # IF CBE is only base editor selected, THEN include only CBE AND CBE and ABE
      #       if(enzymeClass.Reactive() == "CBE") {
      #         filter(., enzyme == "CBE" | enzyme == "CBE and ABE") %>%
      #           filter(
      #             (cbe_position >= input$editing_window[1]) & (cbe_position <= input$editing_window[2])
      #           )
      #       } else {
      #         
      #         # IF ABE is only base editor selected, THEN include only ABE AND CBE and ABE
      #         filter(., enzyme == "ABE" | enzyme == "CBE and ABE") %>%
      #           filter(
      #             (abe_position >= input$editing_window[1]) & (abe_position <= input$editing_window[2])
      #           )
      #       }
      #     }
      #   }
      #     } %>%
      #   {
      #     if(splice_site.Reactive() == "splice-donors") {
      #       filter(., splice_site == "donor")
      #     } else {
      #       if(splice_site.Reactive() == "splice-acceptors") {
      #         filter(., splice_site == "acceptor")
      #       } else {
      #       .
      #     }
      #     }
      #   } %>%
      # dplyr::rename(Exon = 1, `Splice Site` = 2, `Protospacer` = 3, `PAM` = 4, `Enzyme` = 5, `cDNA Disruption` = 6,
      #               `CBE Position` = 7, `CBE Score` = 8, `ABE Position` = 9, `ABE Score` = 10, `Transcript ID` = 11, `Gene ID` = 12) #%>%

      # # Remove guides that have targets outside of 1-20 bases
      # filter(
      #   (`ABE Position` >= input$editing_window[1] | `CBE Position` >= input$editing_window[1]) &
      #     (`ABE Position` <= input$editing_window[2] | `CBE Position` <= input$editing_window[2])
      # )




      # runSpliceR(ensembl_transcript_id = id.Reactive(), pam = pam.Reactive(), species = species.Reactive())[[1]] %>%
        # mutate(ABE = grepl("ABE", enzyme)) %>%
        # mutate(CBE = grepl("CBE", enzyme)) %>%
        # mutate(splice_site = paste0("splice-",splice_site)) %>%
        # 
        # {if(enzyme.Reactive() == "CBE") {filter(., CBE)} else {
        #   if(enzyme.Reactive() == "ABE") {filter(., ABE)} else {
        #     if(enzyme.Reactive() == "CBE and ABE") {filter(., CBE & ABE)} else {.}
        #   }
        # }
        #   } %>%
        # 
        # {if(splice_site.Reactive() == "splice-donor") {filter(., splice_site == "splice-donor")} else {
        #   if(splice_site.Reactive() == "splice-acceptor") {filter(., splice_site == "splice-acceptor")} else {.}
        # }
        #   }
        # 
      }
        )
    
    output$print = renderPrint({
      cat(
          paste(
            c("Ensembl transcript ID:", "PAM:\t\t", "Species:\t", "Enzyme:\t\t", "Splice-site:\t", "Editing window:\t", "Filtering:\t", "Scoring enzyme:\t"),
            c(id.Reactive(), pam.Reactive(), species.Reactive(), enzymeClass.Reactive(), splice_site.Reactive(),
              paste0(input$editing_window, collapse = " - "), 
              if(input$strictFilter) {"CBE AND ABE within window"} else {"CBE OR ABE within window"},
              "rAPOBEC1 and evoTadA 7.10"),
        sep = "\t"),
        sep = "\n"
      )
    })

    output$ensembl <- renderUI({
      withProgress(tags$iframe(src = runSpliceR(ensembl_transcript_id = id.Reactive(), pam = pam.Reactive(), species = species.Reactive()
                                                )[[3]],
                             height = 1000, width = 1500),
      message = "Rendering webpage"
      )
    })
    
    output$behive <- renderUI({
      withProgress(tags$iframe(src = "https://www.crisprbehive.design/single",
                               height = 1000, width = 1500),
                   message = "Rendering BE-hive webpage"
      )
    })
    
    

    output$df <- DT::renderDT(
      withProgress(output.df(),
                   message = 'Generating guide RNAs'),
      options = list(lengthChange = FALSE)
      )

    # Render .csv file from dataTable
    output$downloadAllData <- downloadHandler(
      filename = function() {
        paste(id.Reactive(), pam.Reactive(), "SpliceR_gRNAs.csv", sep = "_")
      },
      content = function(file) {
        write.csv(output.df(), file, row.names = FALSE)
      }
    )
    
   # Render .csv file for selected rows from dataTable
    output$downloadSelectedData <- downloadHandler(
      filename = function() {
        paste(id.Reactive(), pam.Reactive(), "SpliceR_gRNAs.csv", sep = "_")
      },
      content = function(file) {
        write.csv(output.df()[input$df_rows_selected,], file, row.names = FALSE)
      }
    )
    
  }
)
