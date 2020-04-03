# Example data
example_id <- "ENST00000380152.7"
example_pam <- "NGG"
example_species <- "Homo_sapiens"
example_enzyme = "CBE or ABE"
example_splice_site = "splice-donors and splice-acceptors"

## To do:
# Update with note to contact about other species.

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
    
    enzyme.Reactive = reactive({
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
    
    
    # Render table with gRNAs
    output_df = reactive({
      generate_guides(ensembl_id = id.Reactive(), pam = pam.Reactive(), species = species.Reactive())[[1]] %>%
        filter(enzyme == enzyme.Reactive()) %>%
        mutate(splice_site = paste0("splice-", splice_site)) %>%
        filter(mapply(FUN = grepl, pattern = splice_site, x = splice_site.Reactive())) %>%
        dplyr::select(exon, splice_site, guide, pam, be_efficiency, be_dinucleotide, be_position, abe_efficiency, abe_dinucleotide, abe_position, full_motif, gene_id, id) %>%
        dplyr::rename(Exon = 1, `Splice-site` = 2, Protospacer = 3, PAM = 4,
                      `CBE efficiency` = 5, `CBE dinucleotide` = 6, `CBE position` = 7,
                      `ABE efficiency` = 8, `ABE dinucleotide` = 9, `ABE position` = 10,
                      `Motif` = 11, `Gene ID` = 12, `Transcript ID` = 13)
                      
        # Columns to include
        # Exon, Splice Site, Guide, PAM, CBE Efficiency, ABE Efficiency, Dinucleotide, Position, Motif, Gene ID, transcript ID
        
      }
        )
    
    output$print = renderPrint({
      cat(
          paste(
            c("Ensembl transcript ID:", "PAM:", "Species:", "Enzyme:", "Splice-site:"),
            c(id.Reactive(), pam.Reactive(), species.Reactive(), enzyme.Reactive(), splice_site.Reactive()),
        sep = "\t"),
        sep = "\n"
      )
    })

    output$frame <- renderUI({
      withProgress(tags$iframe(src = generate_guides(ensembl_id = id.Reactive(), pam = pam.Reactive(), species = species.Reactive())[[2]],
                             height = 1000, width = 1500),
      message = "Rendering webpage"
      )
    })

    output$df <- DT::renderDT(
      withProgress(output_df(),
                   message = 'Generating guide RNAs')
      )

    # Render .csv file from dataTable
    output$downloadData <- downloadHandler(
      filename = function() {
        paste(input$id, input$pam, "gRNAs.csv", sep = "_")
      },
      content = function(file) {
        write.csv(output_df(), file, row.names = FALSE)
      }
    )
    
  }
)
