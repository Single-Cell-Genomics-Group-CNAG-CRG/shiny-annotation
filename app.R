library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(sctransform)
library(Seurat)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(profvis)
library(scattermore)
library(DT)
library(glue)
library(viridis)

# Metadata dataframe set as NULL at the beginning to avoid showing error
if (! exists("metadata_df")) metadata_df <- NULL
if (! exists("expr_mtrx")) expr_mtrx <- NULL

## Only run examples in interactive R sessions
# if (interactive()) {
##################################################################################################################
##################################################################################################################
################################################## INTERFACE #####################################################
##################################################################################################################
##################################################################################################################
ui <- fluidPage(
  sidebarLayout(fluid = TRUE,
    
    #########################
    #### Parameter entry ####
    #########################
    sidebarPanel(width = 3,
     # CSS formatting for description text
     tags$head(tags$style("#text1{color: black;
                 font-size: 40px;
             font-style: bold;
             }"
     )
     ),
     tags$head(tags$style("#text2{color: #5e5e5e;
                 font-size: 20px;
             font-style: italic;
             }"
     )
     ),
      # Load Metadata
      fileInput(inputId = "metadata",
                label = "Metadata + 2D embedings (coord_x, coord_y)",
                multiple = FALSE),
      
      # Load Expression data
      fileInput(inputId = "data",
                label = "Expression matrix",
                multiple = FALSE),
      
      # Select genes to use as markers
      selectizeInput("gene_ls", "Select marker genes:",
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = TRUE),
      actionButton(inputId = "apply_markers", label = "Update markers"),
      
      # Enter value to group by
      selectizeInput("groupby", "Coloring feature:",
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     # multiple = TRUE,
                     multiple = FALSE),
      actionButton(inputId = "apply_groupby", label = "Update Grouping"),
      
      # Slider to determine point size of UMAP plots
      sliderInput("size", "Dot size:",
                  min = 0,
                  max = 20,
                  value = 3,
                  step = 0.1),
     
     # Select gene to to filter expression
     selectizeInput("gene_filt", "Gene to filter:",
                    selected = NULL,
                    choices = NULL,
                    options = list(create = TRUE),
                    multiple = FALSE),
     actionButton(inputId = "apply_expr_filt",
                  label = "Update gene slider"),
     
     # Slider to filter by gene expression
     sliderInput("expression_slider",
                 label = "Gene expression filter",
                 min = 0, 
                 max = 10,
                 step = 0.5,
                 value = c(-5, 10)),
      
      # Which labels to add to interactive output
      selectizeInput("filter_var", "Filtering group:",
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = FALSE),
      
      actionButton(inputId = "apply_filter", label = "Update variable"),
      
      # Which groups to keep
      checkboxGroupInput(inputId = "filter_grp",
                         label = "Groups to include:",
                         choices = "",
                         selected = NULL),
      
      actionButton(inputId = "apply_grp", label = "Update filter"),
      
    ),
    
    ############################
    #### Plot visualization ####
    ############################
    mainPanel(
      tabsetPanel(
        tabPanel(title = "Description",
                 textOutput("text1"),
                 textOutput("text2"),
                 htmlOutput("text3"),
                 htmlOutput("text4")),
        
        tabPanel(title = "UMAP Plot",
                 plotOutput("dimPlot", height = "800")),
        
        tabPanel(title = "Feature Plots",
                 radioButtons(inputId = "feat_col",
                                    label = "Color Palette:",
                                    choices = c("GreyBlue",
                                                "viridis",
                                                "heat",
                                                "magma"),
                                    selected = "GreyBlue",
                                    inline = TRUE
                                    # width = "100%"
                              ),
                 plotOutput("FeaturePlot", height = "800")
                 ),
        
        tabPanel(title = "Violin Plots",
                 plotOutput("ViolinPlot", height = "800")),
        
        # Adding Differential Expression Tab
        tabPanel(title = "DE",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                # Variables specifying groups
                                selectizeInput("de_var", "DE group:",
                                               selected = NULL,
                                               choices = NULL,
                                               options = list(create = TRUE),
                                               multiple = FALSE),
                                actionButton(inputId = "update_indents",
                                             label = "Update groups"),
                                
                                # Variables specifying groups
                                selectizeInput("de_g1", "Group 1:",
                                               selected = NULL,
                                               choices = NULL,
                                               options = list(create = TRUE),
                                               multiple = FALSE),
                                selectizeInput("de_g2", "Group 2:",
                                               selected = NULL,
                                               choices = NULL,
                                               options = list(create = TRUE),
                                               multiple = FALSE),
                                actionButton(inputId = "do_de",
                                             label = "Run DE"),
                                # Button
                                downloadButton("downloadData", "Download")
                                ),
                   mainPanel(
                     DTOutput("de_table"))
                   )
                 ),
        # Adding Differential Expression Tab
        tabPanel(title = "Cell Selection",
                 sidebarLayout(
                   sidebarPanel(width = 3,
                                # Variables specifying groups
                                selectizeInput("sel_gene", "Gene expression:",
                                               selected = NULL,
                                               choices = NULL,
                                               options = list(create = TRUE),
                                               multiple = FALSE),
                                actionButton(inputId = "update_gene_sel",
                                             label = "Update gene"),
                   ),
                   mainPanel(
                     plotlyOutput("sel_plot"),
                     DTOutput("barcode_table")
                     )
                 )
        )
        )
      )
    )
  )

##################################################################################################################
##################################################################################################################
################################################ SERVER ##########################################################
##################################################################################################################
##################################################################################################################

server <- function(input, output, session) {
  # Setting maximum file size to 8GB
  options(shiny.maxRequestSize = 8000 * 1024 ^ 2)
  
  observeEvent(input$apply_checkbox, {
    shinyjs::toggle("interactive_filter")
  })
  
  ##########################################
  ##### Defining environment variables #####
  ##########################################
  # None
  
  ##############################
  ##### Update input fields ####
  ##############################
  observe({
    
    # Load data
    # Read marker list
    file1 <- input$metadata
    if (is.null(file1) ) { return() }
    tmp1_ds <- readRDS(file1$datapath)
    metadata_df <<- tmp1_ds

    file2 <- input$data
    if( is.null( file2 ) ) { return() }
    tmp2_ds <- readRDS(file2$datapath)
    expr_mtrx <<- as.matrix(tmp2_ds)
    
    # Update marker selection
    updateSelectizeInput(session,
                         inputId = "gene_ls",
                         choices = c("", rownames(expr_mtrx)),
                         selected = rownames(expr_mtrx)[1])
    
    # Update filter gene
    updateSelectizeInput(session,
                         inputId = "gene_filt",
                         choices = rownames(expr_mtrx),
                         selected = rownames(expr_mtrx)[1])
    
    
    # Subset character/factor columns with <= 50 unique values
    feat_ch1 <- sapply(metadata_df, function(x) length(unique(x)) <= 50)
    feat_ch2 <- sapply(colnames(metadata_df), function(x) !is.numeric(metadata_df[, x]))
    feat_ch <- colnames(metadata_df)[feat_ch1 & feat_ch2]
    
    # Subset numeric metadata columns
    feat_num1 <- sapply(colnames(metadata_df), function(x) is.numeric(metadata_df[, x]))
    feat_num <- colnames(metadata_df)[feat_num1]
    
    # Join metadata variables of interest subset
    feat_sub <- c(feat_ch, feat_num)

    # Update groupby selection
    updateSelectizeInput(session,
                         inputId = "groupby",
                         choices = c("", 
                                     feat_sub),
                         selected = feat_sub[1])
    
    # Update Filtering variable selection
    updateSelectizeInput(session,
                         inputId = "filter_var",
                         choices = c("", 
                                     feat_ch),
                         selected = feat_ch[1])
    
    # Update gene_sel selection
    updateSelectizeInput(session,
                         inputId = "sel_gene",
                         choices = c("", rownames(expr_mtrx)),
                         selected = rownames(expr_mtrx)[1])
  })
  
  # Differential expression analysis update
  observe({
    
    # Load data
    file1 <- input$metadata
    if (is.null(file1) ) { return() }
    tmp1_ds <- readRDS(file1$datapath)
    metadata_df <<- tmp1_ds
    
    # Exrtact metadata features + groups
    # Subset colnames
    mask <- sapply(metadata_df, function(x) length(unique(x)) <= 50)
    de_vr <- names(mask)[mask]
    
    # Update Filtering variable selection
    updateSelectizeInput(session,
                         inputId = "de_var",
                         choices = c("", 
                                     de_vr),
                         selected = de_vr[1])
  })
  
  observe({
    vr_groups <- unique(as.character(metadata_df[, de_var()]))
    # Update Filtering variable selection
    updateSelectizeInput(session,
                         inputId = "de_g1",
                         choices = c("",
                                     vr_groups),
                         selected = vr_groups[1])
    
    # Update Filtering variable selection
    updateSelectizeInput(session,
                         inputId = "de_g2",
                         choices = c("", 
                                     vr_groups),
                         selected = vr_groups[2])
  })
  
  # Independent observe event for Gene expression slider event input so it doesn't update all the rest  
  observe({
    # Here we want to update the slider to filter by gene expression
    updateSliderInput(session, 
                      inputId = "expression_slider",
                      value = min(expr_mtrx[gene_filt(), ]),
                      min = min(expr_mtrx[gene_filt(), ]),
                      max = max(expr_mtrx[gene_filt(), ]),
                      step = 0.1)
  })
  
  # Independent observe event for checkbox input so it doesn't update all the rest
  observe({
    # Update Filtering groups
    updateCheckboxGroupInput(session,
                             inputId = "filter_grp",
                             choices = unique(metadata_df[, filter_var()]),
                             selected = unique(metadata_df[, filter_var()]))
  })
  
  
  ########################################
  ######## Setting reactive events #######
  ########################################
  gene_list <- eventReactive(input$apply_markers, {
    input$gene_ls
  })
  
  groupby_var <- eventReactive(input$apply_groupby, {
    input$groupby
  })
  
  filter_var <- eventReactive(input$apply_filter, {
    input$filter_var
  })
  
  gene_filt <- eventReactive(input$apply_expr_filt, {
    input$gene_filt
  })
  
  apply_grp <- eventReactive(input$apply_grp, {
    input$filter_grp
  })
  
  dfInput <- reactive({
    ## subsetting is a bit tricky here to id the column on which to subset        
    metadata_df[metadata_df[, filter_var()] %in% apply_grp(), ]
  })

  exprInput <- reactive({
    ## subsetting is a bit tricky here to id the column on which to subset
    keep_id <- metadata_df[metadata_df[, filter_var()] %in% apply_grp(), "barcode"]
    expr_mtrx[, keep_id]
  })
  
  # Differential expression
  de_var <- eventReactive(input$update_indents, {
    input$de_var
  })
  
  de_g1 <- eventReactive(input$do_de, {
    input$de_g1
  })

  de_g2 <- eventReactive(input$do_de, {
    input$de_g2
  })
  
  de_table <- reactive({
    # Some old verision metadatas don't have column names but they are save in the "barcode" column
    # I've updated seurat2shiny so now the metadata has rownames. In time this If can disappear.
    
    # If all the rownames in metadata are in the columns of expression matrix go ahead,
    # if not and barcode is a column in metadata assign it to rownames
    if(sum(rownames(metadata_df) %in% colnames(expr_mtrx)) != nrow(metadata_df) &
       "barcode" %in% colnames(metadata_df)) {
      rownames(metadata_df) <- metadata_df$barcode
    }
    
    se_obj <- Seurat::CreateSeuratObject(counts = expr_mtrx,
                                         meta.data = metadata_df)
    
    Seurat::Idents(se_obj) <- metadata_df[, de_var()]
    
    markers <- Seurat::FindMarkers(object = se_obj,
                                   ident.1 = de_g1(),
                                   ident.2 = de_g2())
    markers <- markers %>% tibble::rownames_to_column("gene")
  })
  
  sel_gene <- eventReactive(input$update_gene_sel, {
    input$sel_gene
  })
  
  ###########################################
  ######### 1st tab App description #########
  ###########################################
  output$text1 <- renderText({
    "Introduction"
  })
  
  output$text2 <- renderText({
    "Please read this carefully before using the app. Here we explain the purpose of the app as well as the file format requirements."
  })
  
  output$text3 <- renderUI({
    HTML("<p>This App is designed to take in 2 RDS files, one containing the metadata of the cells and the second containing the gene expression matrix of choice.<br/>
    These RDS objects can be obtained using the function found <a href='https://github.com/Single-Cell-Genomics-Group-CNAG-CRG/shiny-annotation/blob/master/seurat2shiny.R'> <B>here</B></a>!<br/>
    <B>Before visualizing the plots</B> for the 1st time one must click the update buttons selecting <i>genes of interest</i>, <i>grouping variable</i>, <i>filtering variable</i> and <i>filtering selection</i></p>")
  })
  
  output$text4 <- renderUI({
    HTML("&#8226;<B>Metadata file</B>: this file contains as many rows as cells are in the dataset with information regarding each one.
    Please check that it contains the variables:<br/>
    &nbsp;&#8212;<B>coord_x, coord_y</B>: containing the 2D embedding of the cells<br/>
    &nbsp;&#8212;<B>barcode</B>: containing the cell barcode matching the colnames of the expression matrix<br/>
    &#8226;<B>Expression matrix</B>: this file is a GENExCELL expression matrix with gene names as rownames and cell barcodes as colnames.")
  })
  
  ########################################
  ######### Plot visualization ###########
  ########################################
  
  ### UMAP coloring by metadata features ### 
  output$dimPlot <- renderPlot({

    metadata_df <- dfInput()

    dim_plot <- ggplot2::ggplot(data.frame(x = metadata_df[, "coord_x"],
                                           y = metadata_df[, "coord_y"])) +
      scattermore::geom_scattermore(aes(x,
                                        y,
                                        color = metadata_df[, groupby_var()]),
                                    pointsize = as.numeric(input$size),
                                    alpha = 1,
                                    pixels = c(2000, 2000),
                                    interpolate = TRUE) +
      ggplot2::theme_classic() +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 5))) +
      ggplot2::labs(
        title = "",
        x = "DIM-1",
        y = "DIM-2",
        color = groupby_var())
    
    # Define plot color differences when character coloring variable passed
    if (! is.numeric(metadata_df[, groupby_var()])) {
      # Define the number of colors you want
      nb.cols <- length(unique(metadata_df[, groupby_var()]))
      set2_expand <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nb.cols)
      
      dim_plot <- dim_plot + ggplot2::scale_color_manual(values = set2_expand)
    } else {
      # Define coloring for numerical variable
      dim_plot <- dim_plot + ggplot2::scale_color_gradient(low = "lightgrey",
                                    high = "blue")
        
    }

    return(dim_plot)
  })

  
  ### Feature plot and Violin plot are dependent on apply_markers button to be pressed ###
  
  # Feature plot
  output$FeaturePlot <- renderPlot({
    
    # Read data from reactive observed slots
    metadata_df <- dfInput()
    expr_mtrx <- exprInput()
    
    # Subset by gene expression
    ## Added -0.1 and +0.1 since I think the slider input rounds the number and the extreme numbers might be ignored
    mask <- expr_mtrx[gene_filt(), ] >= input$expression_slider[1] - 0.1 &
            expr_mtrx[gene_filt(), ] <= input$expression_slider[2] + 0.1
    
    expr_mtrx <- expr_mtrx[, mask]
    metadata_df <- metadata_df[mask, ]

    ## Plot all genes
    plt_ls <- lapply(gene_list(), function(gene) {
      if (gene %in% rownames(expr_mtrx)) {
        ## Plot
        feat_plt <- ggplot2::ggplot(data.frame(x = metadata_df[, "coord_x"],
                                               y = metadata_df[, "coord_y"])) +
          scattermore::geom_scattermore(aes(x,
                                            y,
                                            color = expr_mtrx[gene, ]),
                                        pointsize = as.numeric(input$size),
                                        alpha = 0.7,
                                        pixels = c(1000,1000),
                                        interpolate = TRUE) +
          # ggplot2::scale_color_gradient(low = "lightgrey",
          #                               high = "blue",
          #                               limits = c(min(expr_mtrx[gene, ]),
          #                                          max(expr_mtrx[gene, ]))) +
          ggplot2::theme_classic() +
          ggplot2::labs(
            title = gene,
            x = "DIM-1",
            y = "DIM-2",
            color = "Expression") +
          ggplot2::theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        
        # Set color palette
        if (input$feat_col == "GreyBlue") {
          feat_plt <- feat_plt +
            ggplot2::scale_color_gradient(low = "lightgrey",
                                          high = "blue",
                                          limits = c(min(expr_mtrx[gene, ]),
                                                     max(expr_mtrx[gene, ])))
        } else if (input$feat_col %in% c("viridis", "magma")) {
          feat_plt <- feat_plt +
            ggplot2::scale_colour_viridis_c(option = input$feat_col,
                                            limits = c(min(expr_mtrx[gene, ]),
                                                       max(expr_mtrx[gene, ])))
        } else if (input$feat_col == "heat") {
          feat_plt <- feat_plt +
            ggplot2::scale_color_gradient(low = "#FFFF80FF",
                                          high = "#FF0000FF",
                                          limits = c(min(expr_mtrx[gene, ]),
                                                     max(expr_mtrx[gene, ])))
        }
        
      } else {
        warning(sprintf("Gene %s not found in the expression matrix.", gene))
        return(NULL)
      }
    })

    # Arrange all the plots
    feat_arr <- cowplot::plot_grid(plotlist = plt_ls,
                                   align = "hv",
                                   axis = "tbrl")
    
    return(feat_arr)
  })
  
  # Violin plots
  output$ViolinPlot <- renderPlot({
    
    # Read data from reactive observed slots
    metadata_df <- dfInput()
    expr_mtrx <- exprInput()
    
    # Subset by gene expression
    ## Added -0.1 and +0.1 since I think the slider input rounds the number and the extreme numbers might be ignored
    mask1 <- expr_mtrx[gene_filt(), ] >= input$expression_slider[1] - 0.1
    mask2 <- expr_mtrx[gene_filt(), ] <= input$expression_slider[2] + 0.1
    mask <- mask1 & mask2
    
    expr_mtrx <- expr_mtrx[, mask]
    metadata_df <- metadata_df[mask, ]
    
    
    nb.cols <- length(unique(metadata_df[, groupby_var()]))
    set2_expand <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nb.cols)
    
    ## Plot all genes
    vln_ls <- lapply(gene_list(), function(gene) {
      vln_plt <- ggplot2::ggplot(data.frame(x = metadata_df[, groupby_var()],
                                            y = expr_mtrx[gene, ])) +
        ggplot2::geom_violin(aes(x = metadata_df[, groupby_var()],
                                 y = expr_mtrx[gene, ],
                                 color = metadata_df[, groupby_var()],
                                 fill = metadata_df[, groupby_var()]),
                             alpha = 0.6) +
        ggplot2::ylim(c(min(expr_mtrx[gene, ]),
                        max(expr_mtrx[gene, ]))) +
        ggplot2::scale_color_manual(values = set2_expand) +
        ggplot2::scale_fill_manual(values = set2_expand) +
        ggplot2::theme_classic() +
        ggplot2::labs(
          title = gene,
          x = groupby_var(),
          y = sprintf("%s expression", gene),
          color = groupby_var(),
          fill  = groupby_var()) +
        ggplot2::theme(
          plot.title = ggplot2::element_text(
            hjust = 0.5,
            face = "bold")
        )
      
      return(vln_plt)
      })
    
    # Arrange all the plots
    vln_arr <- cowplot::plot_grid(plotlist = vln_ls,
                                   align = "hv",
                                   axis = "tbrl")
    
    return(vln_arr)
  })
  
  # Differential Expression Table, table is generated in a reactive event so it can be picked up by de_table and the download button
  output$de_table <- renderDT({
    # Some old verision metadatas don't have column names but they are save in the "barcode" column
    # I've updated seurat2shiny so now the metadata has rownames. In time this If can disappear.
    
    # If all the rownames in metadata are in the columns of expression matrix go ahead,
    # if not and barcode is a column in metadata assign it to rownames
    # if(sum(rownames(metadata_df) %in% colnames(expr_mtrx)) != nrow(metadata_df) &
    #    "barcode" %in% colnames(metadata_df)) {
    #   rownames(metadata_df) <- metadata_df$barcode
    # }
    # 
    # se_obj <- Seurat::CreateSeuratObject(counts = expr_mtrx,
    #                                      meta.data = metadata_df)
    # 
    # Seurat::Idents(se_obj) <- metadata_df[, de_var()]
    # 
    # markers <- Seurat::FindMarkers(object = se_obj,
    #                                ident.1 = de_g1(),
    #                                ident.2 = de_g2())
    DT::datatable(de_table(),
                  filter = "top",
                  options = list(
                    lengthMenu = c(10, 25, 50),
                    pageLength = 5)
                  )
  })
  
  # Download csv of DE table generated
  output$downloadData <- downloadHandler(
    filename = "DE_shiny.csv",
    content = function(file) {
      write.csv(de_table(),
                file = file,
                row.names = FALSE)
    })
  
  # Lasso selection plot
  output$sel_plot <- renderPlotly({
    
    # Read data from reactive observed slots
    metadata_df <- dfInput()
    expr_mtrx <- exprInput()
    
    p <- ggplot2::ggplot(data = metadata_df) +
      ggplot2::geom_point(ggplot2::aes(x = coord_x,
                                       y = coord_y,
                                       color = expr_mtrx[sel_gene(), ],
                                       text = barcode),
                          size = as.numeric(input$size)) +
      ggplot2::theme_classic() +
      ggplot2::labs(
        title = glue::glue("Gene: {sel_gene()}"),
        x = "DIM-1",
        y = "DIM-2",
        color = "Expression") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold")) +
      ggplot2::scale_color_gradient(low = "lightgrey",
                                    high = "blue",
                                    limits = c(min(expr_mtrx[sel_gene(), ]),
                                               max(expr_mtrx[sel_gene(), ])))
    
    ggplotly(p, tooltip = "text") %>%
      layout(dragmode = "lasso")
  })
  
  # Lasso selection table
  output$barcode_table <- DT::renderDT(
    # It is important here to set server = FALSE so that when we save the table in CSV format i saves ALL the entries. By default it is TRUE and it only saves those entries currently shown in the table!
    # Note that your users might run into performance and memory issues using server=FALSE if your data table is very large.
    # https://stackoverflow.com/questions/50508854/button-extension-to-download-all-data-or-only-visible-data
    server = FALSE,
    {
    # Read data from reactive observed slots
    metadata_df <- dfInput()
    metadata_df <- metadata_df %>%
      # Round and make as character so we can then join with the d dataframe representing the coordinates in the shiny app.
      dplyr::mutate(
        coord_x = as.character(round(coord_x, 6)),
        coord_y = as.character(round(coord_y, 6))
      )
      

    d <- event_data(event = "plotly_selected")
    if (!is.null(d)) {
      d <- d %>%
        # we need to round the numbers to 6 since those are the coord that the metadata gives, shiny app returns up tp 14 decimals so the left join doesn't work
        dplyr::mutate(
          x = as.character(round(x, 6)),
          y = as.character(round(y, 6))
        ) %>%
        # dplyr::mutate(y = -y) %>%
        dplyr::left_join(metadata_df,
                         by = c("y" = "coord_y",
                                "x" = "coord_x")) %>%
        dplyr::select("pointNumber", "x", "y", "barcode")
      
      filename <- glue::glue("{sel_gene()}-shinyapp")
      # Return datatable with the csv option to save the table directly
      # Download options following this -> https://rstudio.github.io/DT/003-tabletools-buttons.html
      # Also this -> https://github.com/rstudio/DT/issues/409
      DT::datatable(
        data = d,
        extensions = c("Buttons"),
        options = list(
          dom = 'Bfrtip',
          buttons =  list(list(extend = "csv", filename = filename))
          # buttons = list(list(
          #   extend = 'collection',
          #   buttons = c('csv', 'excel', 'pdf'),
          #   text = 'Download'
          # ))
            # list("copy", "print", list(
            #   extend = "collection",
            #   buttons = c("csv", "excel", "pdf"),
            #   filename = filename,
            #   text = "Download"
            # ))
        )
      )
    }
  })
}

shinyApp(ui, server)

# Code profiling
# profvis(runApp())
