ui <- bootstrapPage(
  mainPanel(
    titlePanel("Classification app"),
    
    tabsetPanel(
      
      tabPanel("Data input", 
               p("Before uploading your data, check that it is clean, especially ensure that the the numeric variables contain only the digits 0-9 or NA (to indicate missing data)."),
               p("Rows that contain one or more NAs will be excluded from the PCA."),
               p("Columns that contain a mixture of numbers and text will not be included in the computation of the PCA results."),
               tags$hr(),
               p("Select the options that match your CSV file, then upload your file:"),
               
               
               radioButtons(inputId = 'header',  
                            label = 'Header',
                            choices = c('Columns have headers'='Yes',
                                        'Columns do not have headers'='No'), 
                            selected = 'Yes'),
               
               radioButtons('sep', 'Separator',
                            c(Comma=',',
                              Semicolon=';',
                              Tab='\t'),
                            ','),
               
               radioButtons('quote', 'Quote',
                            c(None='',
                              'Double Quote'='"',
                              'Single Quote'="'"),
                            '"'),
               
               tags$hr(),
               
               fileInput('file1', 'Choose a CSV file to upload:',
                         accept = c(
                           'text/csv',
                           'text/comma-separated-values',
                           'text/tab-separated-values',
                           'text/plain',
                           '.csv',
                           '.tsv'
                         )),
               p("After uploading your CSV file, click on the 'Data Exploration' tab")
               
      ), # end file  tab
      
      tabPanel("Data Exploration",
               
               p("Here is a summary of the data"),
               tableOutput('summary'),
               tags$hr(),
               p("Here is the raw data from the CSV file"),
               DT::dataTableOutput('contents')
      ), # end  tab
      
      
      tabPanel("Correlation Plots",
               uiOutput("choose_columns_biplot"),
               tags$hr(),
               p("This plot may take a few moments to appear when analysing large datasets. You may want to exclude highly correlated variables from the PCA."),
               
               plotOutput("corr_plot"),
               tags$hr(),
               p("Summary of correlations"),
               tableOutput("corr_tables")
      ), # end  tab
      
      tabPanel("Diagnostics",
               
               p("Among SPSS users, these tests are considered to provide some guidelines on the suitability of the data for a principal components analysis. However, they may be safely ignored in favour of common sense. Variables with zero variance are excluded."),
               tags$hr(),
               p("Here is the output of Bartlett's sphericity test. Bartlett's test of sphericity tests whether the data comes from multivariate normal distribution with zero covariances. If p > 0.05 then PCA may not be very informative"),
               verbatimTextOutput("bartlett"),
               tags$hr(),
               p("Here is the output of the Kaiser-Meyer-Olkin (KMO) index test. The overall measure varies between 0 and 1, and values closer to 1 are better. A value of 0.6 is a suggested minimum. "),
               verbatimTextOutput("kmo")
               
               
               
      ), # end  tab
      
      tabPanel("Compute PCA",
               
               p("Choose the columns of your data to include in the PCA."),
               p("Only columns containing numeric data are shown here because PCA doesn't work with non-numeric data."),
               p("The PCA is automatically re-computed each time you change your selection."),
               p("Observations (ie. rows) are automatically removed if they contain any missing values."),
               p("Variables with zero variance have been automatically removed because they're not useful in a PCA."),
               uiOutput("choose_columns_pca"),
               tags$hr(),
               p("Select options for the PCA computation (we are using the prcomp function here)"),
               radioButtons(inputId = 'center',  
                            label = 'Center',
                            choices = c('Shift variables to be zero centered'='Yes',
                                        'Do not shift variables'='No'), 
                            selected = 'Yes'),
               
               radioButtons('scale.', 'Scale',
                            choices = c('Scale variables to have unit variance'='Yes',
                                        'Do not scale variables'='No'), 
                            selected = 'Yes')
               
      ), # end  tab
      
      
      
      tabPanel("PC Plots",
               h2("Scree plot"),
               p("The scree plot shows the variances of each PC, and the cumulative variance explained by each PC (in %) "),
               plotOutput("plot2", height = "300px"),
               tags$hr(),
               h2("PC plot: zoom and select points"),
               p("Select the grouping variable."),
               p("Only variables where the number of unique values is less than 10% of the total number of observations are shown here (because seeing groups with 1-2 observations is usually not very useful)."),
               uiOutput("the_grouping_variable"),
               tags$hr(),
               p("Select the PCs to plot"),
               uiOutput("the_pcs_to_plot_x"),
               uiOutput("the_pcs_to_plot_y"),
               tags$hr(),
               
               p("Click and drag on the first plot below to zoom into a region on the plot. Or you can go directly to the second plot below to select points to get more information about them."),
               p("Then select points on zoomed plot below to get more information about the points."),
               p("You can click on the 'Compute PCA' tab at any time to change the variables included in the PCA, and then come back to this tab and the plots will automatically update."),
               plotOutput ("z_plot1", height = 400,
                           brush = brushOpts(
                             id = "z_plot1Brush",
                             resetOnNew = TRUE)),
               tags$hr(),
               
               p("Click and drag on the plot below to select points, and inspect the table of selected points below"),
               
               plotOutput("z_plot2", height = 400,
                          brush = brushOpts(
                            id = "plot_brush_after_zoom",
                            resetOnNew = TRUE)),
               tags$hr(),
               p("Details of the brushed points"),
               tableOutput("brush_info_after_zoom")
      
      ), # end  tab 
      
      tabPanel('Linear Discriminant Analysis',
               p("In this tab you will be able to apply the Linear Discriminant Analysis. It will help you to assess if there are variables on the dataset that are able to separate the groups of your target variable."),
               p("Since it can be used also for classification, the app will execute a pipeline of training and testing on your data"),
               h2("Target Selection"),
               p("This app is a collection of classification tools, so please select your target variable"),
               p("IMPORTANT! Be sure to select a categorical variable with two groups."),
               uiOutput("target"),
               h2("Columns selection"),
               p("Furthermore, LDA works best with non categorical variables, so to run the algorithm effectively, untick the categorical variables still present in the dataset because of the string.as.factor conversion" ),
               uiOutput('choose_columns_lda'),
               p("Below are displayed the classification performance of the algorithm on the selected columns after having predicted the target variable in a random subsample of 30% of the uploaded dataset"),
               textOutput("lda_performance")
               
      
      ), # end tab
      
      tabPanel('Random Forest with PCA loadings',
               p("In this tab you will be able to apply the Random Forest classifier loaded with PCA outputs. Why this is useful for? It will help you to avoid one of the main issues while dealing with random forest, that is overfitting. In fact, once applied the PCA, the transformed dataset is almost without noise and the Components are computed to appositely describe the most possible variance of the dataset"),
               p("Random Forest can be used with categorical data and also for regression, but for the scope of this app you will able to use it only on dichotomous features"),
               h2("Target Selection"),
               p("This app is a collection of classification tools, so please select your target variable"),
               p("IMPORTANT! Be sure to select a categorical variable with two groups."),
               uiOutput("target_rf_pca"),
               h2("Columns selection"),
               p("Furthermore, PCA works best with non categorical variables, so to run the algorithm effectively, untick the categorical variables still present in the dataset because of the string.as.factor conversion" ),
               uiOutput('choose_columns_rf_pca'),
               p("Below are displayed the classification performance of the algorithm on the selected columns after having predicted the target variable in a random subsample of 30% of the uploaded dataset"),
               textOutput("rf_pca_performance_print"),
               
               
               
      ), #end tab
      
      tabPanel("Logistic Regression",
               sidebarLayout(sidebarPanel(
                 uiOutput("var1_select"),
                 uiOutput("rest_var_select")),
                 mainPanel( helpText("Your results"),
                            verbatimTextOutput("other_val_show")))
      )
      
      
      
    ))) 