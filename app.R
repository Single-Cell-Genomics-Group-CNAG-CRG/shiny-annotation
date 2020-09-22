library(BiocManager)
options(repos = BiocManager::repositories())
library(shiny)
library(shinyjs)
library(Seurat)
library(ggplot2)
library(plotly)
library(RColorBrewer)
library(profvis)
library(scattermore)

# Metadata dataframe set as NULL at the beginning to avoid showing error
if (! exists("metadata_df")) metadata_df <- NULL

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
                  max = 10,
                  value = 3,
                  step = 0.1),
      
      # Which labels to add to interactive output
      selectizeInput("filter_var", "Filtering group:",
                     selected = NULL,
                     choices = NULL,
                     options = list(create = TRUE),
                     multiple = FALSE),
      
      actionButton(inputId = "apply_filter", label = "Update variable"),
      
      # Which groups to keep
      checkboxGroupInput("filter_grp", "Groups to include:",
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
                 plotOutput("FeaturePlot", height = "800")),
        tabPanel(title = "Violin Plots",
                 plotOutput("ViolinPlot", height = "800"))
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
    These RDS objects can be obtained using the function found <a href='https://github.com/MarcElosua/Shiny_annotation/blob/master/seurat_preprocess_fun.R'> <B>here</B></a>!<br/>
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
                                    alpha = 0.7,
                                    pixels = c(2000, 2000),
                                    interpolate = TRUE) +
      ggplot2::theme_classic() +
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
          ggplot2::scale_color_gradient(low = "lightgrey",
                                        high = "blue") +
          ggplot2::theme_classic() +
          ggplot2::labs(
            title = gene,
            x = "DIM-1",
            y = "DIM-2",
            color = "Expression") +
          ggplot2::theme(
            plot.title = element_text(hjust = 0.5, face = "bold")
          )
        
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
          plot.title = ggplot2::element_text(hjust = 0.5,
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
  
}

shinyApp(ui, server)

# Code profiling
# profvis(runApp())
