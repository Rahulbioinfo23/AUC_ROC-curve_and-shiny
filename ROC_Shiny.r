# Load necessary libraries 
library(shiny) 
library(Biobase) 
library(GEOquery) 
library(pROC) 

# Define UI 
ui <- fluidPage( 
 titlePanel("ROC Curve Plotter"), 
 sidebarLayout( 
 sidebarPanel(
helpText("The AUC (Area Under the Curve) value is a single value that 
represents the overall performance of the 
gene expression levels of the specified gene as a binary classifier across 
all the specified samples. 
In the context of ROC analysis, each point on the ROC curve represents a 
sensitivity/specificity pair corresponding 
to a particular decision threshold. The AUC measures the entire twodimensional area underneath the entire ROC curve 
(from 0 to 1) and provides an aggregate measure of performance across all 
possible classification thresholds. 
So, even though you are plotting the ROC curve based on multiple samples, 
the AUC value is a single number that 
summarizes the overall ability of the gene expression levels of the specified 
gene to distinguish between the two conditions 
(e.g., normal vs. disease) across all these samples. 
The single AUC value seen is indeed the correct and expected output for ROC 
analysis."), 
 textInput("gse_id", "GSE ID"), 
 textInput("gpl_id", "Platform"), 
 textInput("samples", "Sample Conditions (0 for normal, 1 for 
disease)"), 
 textInput("gene", "Gene of Interest"), 
 actionButton("plot", "Plot ROC Curve") 
 ), 
 mainPanel( 
 plotOutput("roc_plot") 
 ) 
 ) 
) 

# Define Server Logic 
server <- function(input, output) { 
 # Define the plot_roc function 
 plot_roc <- function(gse_id, gpl_id, samples, gene) { 

 # Load necessary libraries 
 library(Biobase) 
 library(GEOquery) 
 library(pROC) 

 # Download the dataset from GEO 
 gset <- getGEO(gse_id, GSEMatrix =TRUE, AnnotGPL=TRUE) 

 # Select the platform
if (length(gset) > 1) { 
 idx <- grep(gpl_id, attr(gset, "names")) 
 } else { 
 idx <- 1 
 } 
 gset <- gset[[idx]] 

 # Group membership for all samples 
 gsms <- samples # Use the samples input from the UI 
 sml <- strsplit(gsms, split="")[[1]] 

 # Filter out excluded samples (marked as "X") 
 sel <- which(sml != "X") 
 sml <- sml[sel] 
 gset <- gset[ ,sel] 

 # Check if the gene of interest is in the data 
 if (!(gene %in% featureNames(gset))) { 
 stop("The gene of interest is not in the data.") 
 } 

 # Extract gene expression values 
 gene_expression <- exprs(gset)[gene, ] 

 # Create labels (0 for normal, 1 for disease) 
 labels <- ifelse(sml == "0", 0, 1) 
 # Create a roc object 

 roc_obj <- roc(response = labels, predictor = gene_expression) 
 
 # Plot ROC curve 
 plot(roc_obj, main = paste("ROC Curve for", gene)) 
 
 # Add legend 
 legend("bottomright", legend = paste("Gene:", gene, "AUC:", 
round(auc(roc_obj), 2))) 
 } 
 
 observeEvent(input$plot, { 
 output$roc_plot <- renderPlot({ 

 # Call the function with the inputs from the UI 
 plot_roc(input$gse_id, input$gpl_id, input$samples, input$gene) 
 }) 
 }) 
} 

# Run the app
shinyApp(ui = ui, server = server) 