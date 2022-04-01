## ----setup1, include = FALSE--------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  echo = TRUE,
  message = FALSE,
  warning = FALSE
)

## ----eval=FALSE---------------------------------------------------------------
#  install.packages("MAFDash")

## ----eval=FALSE---------------------------------------------------------------
#  install.packages(c("dplyr","ensurer","ggplot2","tidyr","DT","rmarkdown","knitr","flexdashboard","htmltools","data.table","ggbeeswarm","plotly","circlize","canvasXpress","crosstalk","bsplus","BiocManager","maftools","ComplexHeatmap"))
#  BiocManager::install(c("TCGAbiolinks"))
#  install.packages(devtools)
#  library(devtools)
#  devtools::install_github("ashishjain1988/MAFDash")

## ----eval=FALSE---------------------------------------------------------------
#  library(MAFDash)
#  maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
#  getMAFDashboard(maf, outputFileName="output", outputFileTitle=paste0("MAF Dashboard - Test"),outputFilePath = tempdir())

## ----eval=FALSE---------------------------------------------------------------
#  library("MAFDash")
#  # Download MAF data from TCGA
#  CancerCode <- c("ACC","UVM")
#  inputFolderPath <- tempdir() ## This folder will be created if it doesn't exist
#  #maf <- getMAFdataTCGA(cancerCode = CancerCode, outputFolder = inputFolderPath)

## ----eval=TRUE----------------------------------------------------------------
library(MAFDash)
library(maftools)
maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
generateOncoPlot(read.maf(maf,verbose = FALSE))

## ----eval=TRUE----------------------------------------------------------------
library(MAFDash)
library(maftools)
maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
generateBurdenPlot(read.maf(maf,verbose = FALSE), plotType="Dotplot")
generateBurdenPlot(read.maf(maf,verbose = FALSE), plotType="Barplot")

## ----eval=TRUE----------------------------------------------------------------
library(MAFDash)
library(maftools)
maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
generateMutationTypePlot(read.maf(maf,verbose = FALSE))

## ----eval=TRUE----------------------------------------------------------------
library(MAFDash)
library(maftools)
library(plotly)
maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
plots<-generateTiTvPlot(read.maf(maf,verbose = FALSE))
plotly::subplot(plotly::subplot(plots$tiTvPatterns,plots$TiTv, nrows = 1, widths = c(0.5, 0.25)),plots$barplot,nrows = 2)

## ----eval=TRUE----------------------------------------------------------------
library(MAFDash)
library(maftools)
maf <- system.file("extdata", "test.mutect2.maf.gz", package = "MAFDash")
maf <- read.maf(maf = maf,verbose = FALSE)
l<-generateTCGAComparePlot(maf = maf, cohortName = "test")
l$tcga_compare_plot

## ----eval=FALSE---------------------------------------------------------------
#  library(ggplot2)
#  library(plotly)
#  library(ComplexHeatmap)
#  
#  data(iris)
#  
#  ## Simple ggplot
#  myplot <- ggplot(iris) + geom_point(aes(x=Sepal.Length, y=Sepal.Width, color=Species))
#  
#  ## Save as PNG (provide absolute file path)
#  mycustomimage_png <- file.path(getwd(),"custom_ggplot.png")
#  ggsave(mycustomimage_png, plot=myplot, width=5, height=4)
#  
#  ## Save as PDF (provide absolute file path)
#  mycustomimage_pdf <- file.path(getwd(),"custom_ggplot.pdf")
#  ggsave(mycustomimage_pdf, plot=myplot, width=5, height=4)
#  
#  ## Convert ggplot to plotly
#  myplotly <- ggplotly(myplot)
#  
#  ## Make heatmap with ComplexHeatmap
#  hmdata <- t(iris[,1:4])
#  hmanno <- HeatmapAnnotation(df=data.frame(Species=iris[,5]))
#  myhm <- Heatmap(hmdata, bottom_annotation = hmanno)
#  
#  ## Customizable plotly from https://github.com/mtandon09/Dynamic_Plotly
#  source("https://raw.githubusercontent.com/mtandon09/Dynamic_Plotly/master/make_cutomizable_plotly.R")
#  custom_plotly <- make_customizable_plotly(iris)
#  
#  ## Put together objects/filepaths into a list
#  toyplotlist <- list("ggplot"= myplot,
#                     "plotly"= myplotly,
#                     "PNG"= mycustomimage_png,
#                     "PDF"= mycustomimage_pdf,
#                     "ComplexHeatmap"= myhm,
#                     "Customizable"= custom_plotly
#  )
#  
#  ## Filename to output to
#  html_filename="toy_dash.html"
#  
#  ## Render dashboard
#  getMAFDashboard(plotList = toyplotlist,
#                  outputFileName = html_filename,
#                  outputFileTitle = "Iris")

## ----eval=FALSE---------------------------------------------------------------
#  library(MAFDash)
#  library(TCGAbiolinks)
#  
#  tcga_code <- c("ACC","UVM")
#  #inputFolderPath <- paste0(tempdir()) ## This folder will be created if it doesn't exist
#  caller = "mutect2"
#  title_label = paste0("TCGA-",tcga_code)
#  
#  #maf_files <- getMAFdataTCGA(tcga_code,outputFolder = tempdir(),variant_caller = caller)

## ----eval=FALSE---------------------------------------------------------------
#  # tcga_clinical <- getTCGAClinicalAnnotation#TCGAbiolinks::GDCquery_clinic(project = paste0("TCGA-",tcga_code), type = "clinical")
#  # tcga_clinical$Tumor_Sample_Barcode <- tcga_clinical$submitter_id
#  defaultW <- getOption("warn")
#  options(warn = -1)
#  tcga_clinical<-getTCGAClinicalAnnotation(cancerCodes = tcga_code)
#  options(warn = defaultW)

## ----eval=FALSE---------------------------------------------------------------
#  #maf_files<- system.file("extdata", "tcga_laml.maf.gz", package = "maftools")
#  filtered_mafdata <- do.call("rbind",lapply(maf_files, function(maf_file){filter_maf_chunked(maf_file)}))

## ----eval=FALSE---------------------------------------------------------------
#  filtered_maf <- read.maf(filtered_mafdata, clinicalData = tcga_clinical$annodata,verbose = FALSE)
#  annotation_colors <- getTCGAClinicalColors(ageRange = range(tcga_clinical$annodata$age_at_diagnosis, na.rm=T))

## ----eval=FALSE---------------------------------------------------------------
#  custom_onco <- generateOncoPlot(filtered_maf,
#                                  add_clinical_annotations = names(annotation_colors),
#                                  clin_data_colors = annotation_colors)
#  custom_onco

## ----eval=FALSE---------------------------------------------------------------
#  tcgaComparePlot<-generateTCGAComparePlot(maf = filtered_maf, cohortName = "test")
#  tcgaComparePlot$tcga_compare_plot

## ----eval=FALSE---------------------------------------------------------------
#  #ribbonplot_file <- file.path(getwd(),"ribbon.pdf")
#  generateRibbonPlot(filtered_maf,save_name = NULL)

## ----eval=FALSE---------------------------------------------------------------
#  customplotlist <- list("summary_plot"=T,
#                         "burden"=T,
#                         "TCGA Comparison"=tcgaComparePlot$tcga_compare_plot,
#                         "oncoplot"=T,
#                         "Annotated Oncoplot"=custom_onco
#                         )
#  
#  ## Filename to output to; if output directory doesn't exist, it will be created
#  html_filename=file.path(paste0(tempdir(),"/TCGA-UVM.custom.mafdash.html"))
#  
#  ## Render dashboard
#  getMAFDashboard(MAFfilePath = filtered_maf,
#                  plotList = customplotlist,
#                  outputFileName = html_filename,
#                  outputFileTitle = "Customized Dashboard")
#  
#  

