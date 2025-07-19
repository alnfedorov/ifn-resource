library(iSEE)
library(SummarizedExperiment)

# Load the data
tpm <- read.csv("ld/results/iSEE/tpm.csv.gz", row.names = 1)
zscore <- read.csv("ld/results/iSEE/zscore.csv.gz", row.names = 1)
rld <- read.csv("ld/results/iSEE/rld.csv.gz", row.names = 1)

rowData = read.csv("ld/results/iSEE/row_data.csv.gz", row.names = 1)
colData = read.csv("ld/results/iSEE/col_data.csv.gz", row.names = 1)

stopifnot(rownames(rowData) == rownames(tpm))
stopifnot(rownames(rowData) == rownames(zscore))
stopifnot(rownames(rowData) == rownames(rld))

stopifnot(rownames(colData) == colnames(tpm))
stopifnot(rownames(colData) == colnames(zscore))
stopifnot(rownames(colData) == colnames(rld))

se <- SummarizedExperiment(
  assays = list(
    TPM = as.matrix(tpm),
    "Z-score" = as.matrix(zscore),
    "rlog-transformed" = as.matrix(rld)
  ),
  rowData = rowData,
  colData = colData
)

# Define the initial state of the app
initial <- list(
  FeatureAssayPlot(
    Assay = "TPM",
    YAxisFeatureName = "IFI44-201",
    XAxis = "Column data",
    XAxisColumnData = "Treatment",
    ColorBy = "Column data",
    ColorByColumnData = "Donor",
    PointSize = 3,
    LegendPointSize = 3,
    LegendPosition = "Right",
    PanelWidth = 5L
  ),
  ComplexHeatmapPlot(
    Assay = "Z-score",
    CustomRows = TRUE,
    CustomRowsText = "MIR4432HG-204\nMIR4432HG-203\nMIR4432HG-225\nMIR155HG[201/202]\nMIR3945HG-202\nMIR22HG[201/226]\nMIR223HG-201\nMIR223HG-202\nMIR3142HG-201\nMIR23AHG-201\nMIR600HG-201",
    ClusterRows = FALSE,
    DataBoxOpen = FALSE,
    VisualChoices = "Annotations",
    # Note: The 'ColumnData' parameter is not working as expected.
    # It sorts the columns and swaps Treatment <-> Donor levels.
    # ColumnData = c("Treatment", "Donor"),
    DivergentColormap = "PRGn",
    ShowDimNames = c("Rows", "Columns"),
    LegendPosition = "Bottom",
    LegendDirection = "Horizontal",
    NamesRowFontSize = 14,
    NamesColumnFontSize = 14,
    PanelWidth = 7L
  )
)

# Define the tour steps
element <- intro <- character(0)

element <- c(element, "#Welcome")
intro <-
  c(
    intro,
    "Welcome to the interactive tour of the data used in <strong>Figure 3</strong> of
    <strong>'Single-cell analysis of signalling and transcriptional responses to type I interferons'</strong>.
    <br/><br/>
    You will be shown around the different components of the app and learn the basic usage mechanisms by doing.
    Highlighted elements will respond to the user's actions, while the rest of the user interface will be shaded.
    You will be prompted at particular spots to perform some actions, which will be marked with \"<strong>Action:</strong>\" text.
    Please take care to follow these instructions, since later parts of the tour may assume that all the actions from
    previous steps have been performed.
    <br/><br/>
    If you leave the tour at any point and want to return, please click on the '?' in the right hand corner at the top of the page
    <br/><br/>
    <strong>Action:</strong> click on the 'Next' button to proceed"
  )

element <- c(element, "#allpanels")
intro <-
  c(
    intro,
    "This app provides a Shiny interface that allows you to generate a series of panel types for exploring the data.
    These are described in more detail in the following steps. Please wait for the panels to finish loading before
    beginning your tour."
  )

element <- c(element, "#FeatureAssayPlot1")
intro <-
  c(
    intro,
    "This panel displays expression values (e.g., TPMs) for a single transcript
    across different samples and donors (A, B, and C)."
  )

element <- c(element, "#FeatureAssayPlot1_DataBoxOpen")
intro <-
  c(
    intro,
    "Select the transcript you would like to view from the dropdown list under Y-axis feature in the
    <font color=\"#402ee8\">Data parameters</font> menu<br/><br/>"
  )

element <- c(element, "#ComplexHeatmapPlot1")
intro <-
  c(
    intro,
    "This panel shows the expression of selected transcripts in a heatmap. Here, we are viewing
    transcripts from microRNA host genes (MIR-HG).
    <br/><br/>
    Each column is a sample, and each row is a transcript. Samples are ordered by treatment and donor."
  )

element <- c(element, "#ComplexHeatmapPlot1_DataBoxOpen")
intro <-
  c(
    intro,
    "Rows can be clustered by checking the box under the <font color=\"#402ee8\">Data parameters</font> header.
    <br/><br/> Transcript of interest are entered by clicking the 'Edit feature names' box."
  )

element <- c(element, "#Figures")
intro <- c(intro, "To download any of the figure panels, click on the Download sign in the top right hand corner at the top of the page")

element <- c(element, "#Thanks")
intro <- c(intro, "Thank you for taking the tour and enjoy exploring the data!<br/><br/>To repeat the tour, please click on the '?' in the right hand corner at the top of the page")

tour <- data.frame(
  element = element,
  intro = intro,
  stringsAsFactors = FALSE
)

# Launch the iSEE app
app <- iSEE(se, initial = initial, tour = tour)
shiny::runApp(app, port = 8080, launch.browser = FALSE)
