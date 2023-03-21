# global items 

# check if pkgs are installed already, if not, install automatically:
# (http://stackoverflow.com/a/4090208/1036500)
list.of.packages <- c("ggplot2", 
                      "DT", 
                      "GGally",
                      "psych",
                      "Hmisc",
                      "MASS",
                      "digest",
                      "dplyr",
                      "randomForest",
                      "rattle")



new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# load all these
lapply(list.of.packages, require, character.only = TRUE)




server <- function(input, output) {
  
  # read in the CSV
  the_data_fn <- reactive({
    inFile <- input$file1
    if (is.null(inFile)) return(the_data=read.csv("GTD_expanded.csv", header = TRUE))
    the_data <-   read.csv(inFile$datapath, header = (input$header == "Yes"),
                           sep = input$sep, quote = input$quote, stringsAsFactors=TRUE)
    the_data <- dplyr::select_if(the_data, is.numeric)
      
    return(the_data)
  })
  
  
  # display a table of the CSV contents
  output$contents <-  DT::renderDataTable({
    #
    the_data_fn()
  })
  
  # display a summary of the CSV contents
  output$summary <-  renderTable({
    the_data <- the_data_fn()
    psych::describe(the_data)
  })
  
  # Check boxes to choose columns
  output$choose_columns_biplot <- renderUI({
    
    the_data <- the_data_fn()
    
    colnames <- names(the_data)
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("columns_biplot", "Choose up to five columns to display on the scatterplot matrix", 
                       choices  = colnames,
                       selected = colnames[1:5])
  })
  
  # corr plot
  output$corr_plot <- renderPlot({
    the_data <- the_data_fn()
    # Keep the selected columns
    columns_biplot <-    input$columns_biplot
    the_data_subset_biplot <- the_data[, columns_biplot, drop = FALSE]
    ggpairs(the_data_subset_biplot)
  })
  
  # corr tables
    output$corr_tables <- renderTable({
    the_data <- the_data_fn()
    # we only want to show numeric cols
    the_data_num <- the_data[,sapply(the_data,is.numeric)]
    # exclude cols with zero variance
    the_data_num <- the_data_num[,!apply(the_data_num, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    
    res <- Hmisc::rcorr(as.matrix(the_data_num))
    cormat <- res$r
    pmat <- res$P
    ut <- upper.tri(cormat)
    df <- data.frame(
      row = rownames(cormat)[row(cormat)[ut]],
      column = rownames(cormat)[col(cormat)[ut]],
      cor  = (cormat)[ut],
      p = pmat[ut]
    )
    with(df, df[order(-cor), ])
    
  })
  
  output$bartlett <- renderPrint({
    the_data <- the_data_fn()
    the_data_num <- na.omit(the_data[,sapply(the_data,is.numeric)])
    # exclude cols with zero variance
    the_data_num <- the_data_num[,!apply(the_data_num, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    cortest.bartlett(cor(the_data_num), n = nrow(the_data_num))
  })  
  
  output$kmo <- renderPrint({
    the_data <- the_data_fn()
    the_data_num <- the_data[,sapply(the_data,is.numeric)]
    # exclude cols with zero variance
    the_data_num <- the_data_num[,!apply(the_data_num, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    # R <- cor(the_data_num)
    # KMO(R)
    
    # http://www.opensubscriber.com/message/r-help@stat.math.ethz.ch/7315408.html
    # KMO Kaiser-Meyer-Olkin Measure of Sampling Adequacy 
    kmo = function( data ){ 
      
      library(MASS) 
      X <- cor(as.matrix(data)) 
      iX <- ginv(X) 
      S2 <- diag(diag((iX^-1))) 
      AIS <- S2%*%iX%*%S2                      # anti-image covariance matrix 
      IS <- X+AIS-2*S2                         # image covariance matrix 
      Dai <- sqrt(diag(diag(AIS))) 
      IR <- ginv(Dai)%*%IS%*%ginv(Dai)         # image correlation matrix 
      AIR <- ginv(Dai)%*%AIS%*%ginv(Dai)       # anti-image correlation matrix 
      a <- apply((AIR - diag(diag(AIR)))^2, 2, sum) 
      AA <- sum(a) 
      b <- apply((X - diag(nrow(X)))^2, 2, sum) 
      BB <- sum(b) 
      MSA <- b/(b+a)                        # indiv. measures of sampling adequacy 
      
      AIR <- AIR-diag(nrow(AIR))+diag(MSA)  # Examine the anti-image of the 
      # correlation matrix. That is the 
      # negative of the partial correlations, 
      # partialling out all other variables. 
      
      kmo <- BB/(AA+BB)                     # overall KMO statistic 
      
      # Reporting the conclusion 
      if (kmo >= 0.00 && kmo < 0.50){ 
        test <- 'The KMO test yields a degree of common variance 
      unacceptable for FA.' 
      } else if (kmo >= 0.50 && kmo < 0.60){ 
        test <- 'The KMO test yields a degree of common variance miserable.' 
      } else if (kmo >= 0.60 && kmo < 0.70){ 
        test <- 'The KMO test yields a degree of common variance mediocre.' 
      } else if (kmo >= 0.70 && kmo < 0.80){ 
        test <- 'The KMO test yields a degree of common variance middling.' 
      } else if (kmo >= 0.80 && kmo < 0.90){ 
        test <- 'The KMO test yields a degree of common variance meritorious.' 
      } else { 
        test <- 'The KMO test yields a degree of common variance marvelous.' 
      } 
      
      ans <- list(  overall = kmo, 
                    report = test, 
                    individual = MSA, 
                    AIS = AIS, 
                    AIR = AIR ) 
      return(ans) 
      
    }    # end of kmo() 
    kmo(na.omit(the_data_num))
    
  }) 
  
  
  
  # Check boxes to choose columns
  output$choose_columns_pca <- renderUI({
    
    the_data <- the_data_fn()
    
    # Get the data set with the appropriate name
    
    # we only want to show numeric cols
    the_data_num <- na.omit(the_data[,sapply(the_data,is.numeric)])
    # exclude cols with zero variance
    the_data_num <- the_data_num[,!apply(the_data_num, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    
    colnames <- names(the_data_num)
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("columns", "Choose columns", 
                       choices  = colnames,
                       selected = colnames)
  })
  
  # choose a grouping variable
  output$the_grouping_variable <- renderUI({
    the_data <- the_data_fn()
    
    
    # for grouping we want to see only cols where the number of unique values are less than 
    # 10% the number of observations
    grouping_cols <- sapply(seq(1, ncol(the_data)), function(i) length(unique(the_data[,i])) < nrow(the_data)/10 )
    
    the_data_group_cols <- the_data[, grouping_cols, drop = FALSE]
    # drop down selection
    selectInput(inputId = "the_grouping_variable", 
                label = "Grouping variable:",
                choices=c("None", names(the_data_group_cols)))
    
  })
  
  
  pca_objects <- reactive({
    # Keep the selected columns
    columns <-    input$columns
    the_data <- na.omit(the_data_fn())
    the_data_subset <- na.omit(the_data[, columns, drop = FALSE])
    
    # from http://rpubs.com/sinhrks/plot_pca
    pca_output <- prcomp(na.omit(the_data_subset), 
                         center = (input$center == 'Yes'), 
                         scale. = (input$scale. == 'Yes'))
    # data.frame of PCs
    pcs_df <- cbind(the_data, pca_output$x)
    
    return(list(the_data = the_data, 
                the_data_subset = the_data_subset,
                pca_output = pca_output, 
                pcs_df = pcs_df))
    
  })
  
  output$the_pcs_to_plot_x <- renderUI({
    pca_output <- pca_objects()$pca_output$x
    
    # drop down selection
    selectInput(inputId = "the_pcs_to_plot_x", 
                label = "X axis:",
                choices= colnames(pca_output), 
                selected = 'PC1')
  })
  
  output$the_pcs_to_plot_y <- renderUI({
    pca_output <- pca_objects()$pca_output$x
    
    # drop down selection
    selectInput(inputId = "the_pcs_to_plot_y", 
                label = "Y axis:",
                choices= colnames(pca_output), 
                selected = 'PC2')
  })
  
  
  
  output$plot2 <- renderPlot({
    pca_output <- pca_objects()$pca_output
    eig = (pca_output$sdev)^2
    variance <- eig*100/sum(eig)
    cumvar <- paste(round(cumsum(variance),1), "%")
    eig_df <- data.frame(eig = eig,
                         PCs = colnames(pca_output$x),
                         cumvar =  cumvar)
    ggplot(eig_df, aes(reorder(PCs, -eig), eig)) +
      geom_bar(stat = "identity", fill = "white", colour = "black") +
      geom_text(label = cumvar, size = 4,
                vjust=-0.4) +
      theme_bw(base_size = 14) +
      xlab("PC") +
      ylab("Variances") +
      ylim(0,(max(eig_df$eig) * 1.1))
  })
  
  
  # PC plot
  pca_biplot <- reactive({
    pcs_df <- pca_objects()$pcs_df
    pca_output <-  pca_objects()$pca_output
    
    var_expl_x <- round(100 * pca_output$sdev[as.numeric(gsub("[^0-9]", "", input$the_pcs_to_plot_x))]^2/sum(pca_output$sdev^2), 1)
    var_expl_y <- round(100 * pca_output$sdev[as.numeric(gsub("[^0-9]", "", input$the_pcs_to_plot_y))]^2/sum(pca_output$sdev^2), 1)
    labels <- rownames(pca_output$x)
    grouping <- input$the_grouping_variable
    
    if(grouping == 'None'){
      # plot without grouping variable
      pc_plot_no_groups  <- ggplot(pcs_df, 
                                   aes_string(input$the_pcs_to_plot_x, 
                                              input$the_pcs_to_plot_y
                                   )) +
        
        
        geom_text(aes(label = labels),  size = 5) +
        theme_bw(base_size = 14) +
        coord_equal() +
        xlab(paste0(input$the_pcs_to_plot_x, " (", var_expl_x, "% explained variance)")) +
        ylab(paste0(input$the_pcs_to_plot_y, " (", var_expl_y, "% explained variance)")) 
      # the plot
      pc_plot_no_groups
      
      
    } else {
      # plot with grouping variable
      
      pcs_df$fill_ <-  as.character(pcs_df[, grouping, drop = TRUE])
      pc_plot_groups  <- ggplot(pcs_df, aes_string(input$the_pcs_to_plot_x, 
                                                   input$the_pcs_to_plot_y, 
                                                   fill = 'fill_', 
                                                   colour = 'fill_'
      )) +
        stat_ellipse(geom = "polygon", alpha = 0.1) +
        
        geom_text(aes(label = labels),  size = 5) +
        theme_bw(base_size = 14) +
        scale_colour_discrete(guide = FALSE) +
        guides(fill = guide_legend(title = "groups")) +
        theme(legend.position="top") +
        coord_equal() +
        xlab(paste0(input$the_pcs_to_plot_x, " (", var_expl_x, "% explained variance)")) +
        ylab(paste0(input$the_pcs_to_plot_y, " (", var_expl_y, "% explained variance)")) 
      # the plot
      pc_plot_groups
    }
    
    
  })
  
  output$brush_info <- renderTable({
    # the brushing function
    brushedPoints(pca_objects()$pcs_df, input$plot_brush)
  })
  
  
  # for zooming
  output$z_plot1 <- renderPlot({
    
    pca_biplot() 
    
  })
  
  # zoom ranges
  zooming <- reactiveValues(x = NULL, y = NULL)
  
  observe({
    brush <- input$z_plot1Brush
    if (!is.null(brush)) {
      zooming$x <- c(brush$xmin, brush$xmax)
      zooming$y <- c(brush$ymin, brush$ymax)
    }
    else {
      zooming$x <- NULL
      zooming$y <- NULL
    }
  })
  
  
  # for zooming
  output$z_plot2 <- renderPlot({
    
    pca_biplot() + coord_cartesian(xlim = zooming$x, ylim = zooming$y) 
    
    
  })
  
  output$brush_info_after_zoom <- renderTable({
    # the brushing function
    brushedPoints(pca_objects()$pcs_df, input$plot_brush_after_zoom)
  })
  
  
  
  
  # target variable selection for LDA
  
  output$target <- renderUI({
    the_data <- the_data_fn()
    
    # drop down selection
    selectInput(inputId = "target", 
                label = "Target variable:",
                choices=names(the_data))
    
  })
  
  # Check boxes to choose columns for LDA
  output$choose_columns_lda <- renderUI({
    
    the_data <- the_data_fn()
    colnames <- names(the_data)
    
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("columns_lda", "Choose columns", 
                       choices  = colnames,
                       selected = colnames)
  })
  
  
  
  
  
  
  lda_objects <- reactive({
    # Keep the selected columns
    columns <-    input$columns_lda
    
    target <- input$target
    
    
    
    the_data <- the_data_fn()
    
    the_data_subset <- the_data[, columns, drop = FALSE]
    
    the_data_subset <- na.omit(the_data_subset)
    
    
    
    sample <- sample(c(TRUE, FALSE), nrow(the_data_subset), replace = T, prob = c(0.8,0.2)) 
    
    # create an index vector
    #to split the dataset in train,test
    
    train <- the_data_subset[sample, ]
    
    test <- the_data_subset[!sample, ]
    
    target_train <-train[,target]
    
    target_test <- test[,target]

    
    
    
    
    lda.fit <- lda(formula = target_train ~. , data = train[,-target_train])# perform the LDA
    df.pred <- predict(lda.fit, test) # test the model over test partition
    
    prior.class0 <- sum(target_test==0)/length(test[,2])
    prior.class1 <- 1-prior.class0
    #prior.class1
    prior.probability <- cbind("0" =prior.class0, "1" = prior.class1)
    #prior.probability
    accuracyLDA <- table(df.pred$class, target_test)
    overall_errorLDA <- sum(accuracyLDA[1,2], accuracyLDA[2,1])/length(test[,2])*100
    errorLDA <- errorMatrix(target_test, df.pred$class, percentage = TRUE)# 
    
    # assess the goodness of the model by 
    #drawing the confusion matrix
    
    error.class0LDA <- errorLDA[1,3] 
    error.class1LDA <- errorLDA[2,3]
    
    return(list(the_data = the_data, 
                the_data_subset = the_data_subset,
                train=train,
                test=test,
                overall_errorLDA=overall_errorLDA,
                prior.probability=prior.probability,
                lda.fit = lda.fit,
                accuracyLDA=accuracyLDA))
    

  })
  
  
  
  output$lda_performance <- renderPrint({
    
    print(paste0("Overall Error: ", lda_objects()$overall_errorLDA))
    print(paste0("Accuracy: ", sum(diag(lda_objects()$accuracyLDA))/length(lda_objects()$test[,3])*100))
    print(paste0("Prior probability of class 0: ", lda_objects()$prior.probability[1]))
    print(paste0("Prior probability of class 1: ", lda_objects()$prior.probability[2]))
  })
  
  
  # target variable selection for Random Forest with pca loadings
  
  output$target_rf_pca <- renderUI({
    the_data <- the_data_fn()
    
    # drop down selection
    selectInput(inputId = "target_rf_pca_1", 
                label = "Target variable:",
                choices=names(the_data))
    
  })
  
  # Check boxes to choose columns for random forest with Pca loadings
  output$choose_columns_rf_pca <- renderUI({
    
    the_data <- the_data_fn()
    
    
    # show numeric cols
    the_data_num <- na.omit(the_data[,sapply(the_data,is.numeric)])
    # exclude cols with zero variance
    the_data_num <- the_data_num[,!apply(the_data_num, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))]
    
    
    colnames <- names(the_data_num)
    
    # Create the checkboxes and select them all by default
    checkboxGroupInput("columns_rf_pca", "Choose columns", 
                       choices  = colnames,
                       selected = colnames)
  })
  
  
  rf_pca_objects <- reactive({
    
    # Create the random forest algorithm that takes as input the PCA rotations
    
    # define as new S3 model class
    
    train_PCA_RF = function(x,y,ncomp=4,...) {
      f.args=as.list(match.call()[-1])
      pca_obj = princomp(x)
      rf_obj = do.call(randomForest,c(alist(x=pca_obj$scores[,1:ncomp]),f.args[-1]))
      out=mget(ls())
      class(out) = "PCA_RF"
      return(out)    
    }
    
    #print method
    print.PCA_RF = function(object) print(object$rf_obj)
    
    #predict method
    
    predict.PCA_RF = function(object,Xtest=NULL,...) {
      print("predicting PCA_RF")
      f.args=as.list(match.call()[-1])
      if(is.null(f.args$Xtest)) stop("cannot predict without newdata parameter")
      sXtest = predict(object$pca_obj,Xtest) #scale Xtest as Xtrain was scaled before
      return(do.call(predict,c(alist(object = object$rf_obj, #class(x)="randomForest" invokes method predict.randomForest
                                     newdata = sXtest),      #newdata input, see help(predict.randomForest)
                               f.args[-1:-2])))  #any other parameters are passed to predict.randomForest
      
    }
    
    #testTrain predict #
    make.component.data = function(
      inter.component.variance = .9,
      n.real.components = 5,
      nVar.per.component = 20,
      nObs=600,
      noise.factor=.2,
      hidden.function = function(x) apply(x,1,mean),
      plot_PCA =F
    ){
      Sigma=matrix(inter.component.variance,
                   ncol=nVar.per.component,
                   nrow=nVar.per.component)
      diag(Sigma)  = 1
      x = do.call(cbind,replicate(n = n.real.components,
                                  expr = {mvrnorm(n=nObs,                                                           mu=rep(0,nVar.per.component),
                                                  Sigma=Sigma)},
                                  simplify = FALSE)
      )
      if(plot_PCA) plot(prcomp(x,center=T,.scale=T))
      y = hidden.function(x)
      ynoised = y + rnorm(nObs,sd=sd(y)) * noise.factor
      out = list(x=x,y=ynoised)
      pars = ls()[!ls() %in% c("x","y","Sigma")]
      attr(out,"pars") = mget(pars) #attach all pars as attributes
      return(out)
    }
    
    # Keep the selected columns
    columns <-    input$columns_rf_pca
    
    target <- input$target_rf_pca_1
    
    
    
    the_data <- the_data_fn()
    
    the_data_subset <- the_data[, columns, drop = FALSE]
    
    the_data_subset <- na.omit(the_data_subset)
    
    
    
    sample <- sample(c(TRUE, FALSE), nrow(the_data_subset), replace = T, prob = c(0.5,0.5)) # then we create an index vector
    
    #to split the dataset in train,test
    train <- the_data_subset[sample, ]
    
    test <- the_data_subset[!sample, ]
    
    target_train <-train$target
    
    target_test <- test$target
    
    #wrapped Random Forest with the loadings of the first three principal 
    #component 
    train <- train[,!(names(train) %in% target)]
    
    rf2 = train_PCA_RF(train, as.factor(target_train), ntree= 15, ncomp=3)
    
    rf2.pred <- predict(rf2, test, type = "class") #we test the model over test set
    
   
    
    errorwrf <- errorMatrix(target_test, rf2.pred)
    error.class0wrf <- errorwrf[1,3] 
    error.class1wrf <- errorwrf[2,3]  
    #errorwrf
    
    
    overall.errorf2 <- mean (rf2.pred!=target_test)*100
    
    
    return(list(the_data = the_data, 
                the_data_subset = the_data_subset,
                train=train,
                test=test,
                rf2.pred=rf2.pred,
                target_test=target_test,
                errorwrf=errorwrf,
                overall.errorf2=overall.errorf2
    ))
  })
  
  
  
  
  
  output$rf_pca_performance_print <- renderPrint({
    
    
    print(paste0("Overall Error wrapped Random Forest: ", rf_pca_objects()$overall.errorf2,"%"))
    
  })
  
  
  
  output$var1_select<-renderUI({
    selectInput("ind_var_select","Select Independent Var", choices =as.list(names(the_data_fn())),multiple = FALSE)
  })
  
  
  output$rest_var_select<-renderUI({
    checkboxGroupInput("other_var_select","Select other Var",choices =as.list(names(the_data_fn())))
  })
  
  
  output$other_val_show<-renderPrint({
    input$other_var_select
    input$ind_var_select
    f<-the_data_fn()
    
    sample <- sample(c(TRUE, FALSE), nrow(f), replace = T, prob = c(0.5,0.5)) # then we create an index vector
    
    #to split the dataset in train,test
    train <- f[sample, ]
    
    test <- f[!sample, ]
    
    
    
    
    
    library(caret)
    form <- sprintf("%s~%s",input$ind_var_select,paste0(input$other_var_select,collapse="+"))
    print(form)
    
    logreg.fit <-glm(as.formula(form),family=binomial(),data=train)
    
    logreg.test<-predict(logreg.fit,newdata = test,type = 'response')
    logreg.class<-ifelse(logreg.test<0.5,0,1)
    table(logreg.class,test[,input$ind_var_select])
    
    errorlogreg <- errorMatrix(test[,input$ind_var_select], logreg.class)
    print(errorlogreg)
    
    
    
    overall.errorlogreg <- mean (logreg.class!=test[,input$ind_var_select])*100
    print(paste0("Overall Error: ", overall.errorlogreg, "%"))
    
    
    print(summary(logreg.fit))
    
  })
}
  
  
  
  
  