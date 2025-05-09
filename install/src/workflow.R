# Install packages if needed
#install.packages(c("DiagrammeR", "DiagrammeRsvg", "rsvg"))

# Load libraries
library(DiagrammeR)
library(DiagrammeRsvg)
#library(rsvg)

# Create the diagram

graph <- grViz("
digraph workflow {
  
  graph [layout = dot, rankdir = TB]

  node [shape = rectangle, style = filled, fontname = Helvetica, fontsize = 12, width = 3]

  // Environment Setup
  subgraph cluster_0 {
    label = 'Environment Setup'
    color = white
    style = invis;
    node [fillcolor = lightblue, fontcolor = black]
    create_dirs [label = 'Create database directories\n(DBs/kegg, DBs/cazy, DBs/merops, DBs/iprscan)']
    conda_env [label = 'Create conda environment\n(rbimsenv)']
    install_tools [label = 'Install tools:\nKofamScan, dbCAN,\nInterProScan']
  }

  // Database Preparation
  subgraph cluster_1 {
    label = 'Database Preparation'
    color = white
    style = invis;
    node [fillcolor = gold, fontcolor = black]
    kofam_db [label = 'Download & config Kofam database']
    dbcan_db [label = 'Download, build & config dbCAN database']
    merops_db [label = 'Download MEROPS protease.lib\n+ makeblastdb']
    iprscan_db [label = 'Download & config InterProScan data']
  }

  // Annotation Workflows
  subgraph cluster_2 {
    label = 'Annotation Workflows'
    color = white
    style = invis;
    node [fillcolor = mediumpurple1, fontcolor = black]
    interpro [label = 'Run InterProScan']
    kofam [label = 'Run KofamScan']
    dbcan [label = 'Run dbCAN']
    merops [label = 'Run MEROPS (BLASTp)']
  }

  // Integration & Analysis
  subgraph cluster_3 {
    label = 'Integration & Analysis'
    color = white
    style = invis;
    node [fillcolor = lightpink, fontcolor = black]
    rbims [label = 'Integrate annotations in RbiMs']
    analysis [label = 'Downstream analysis\nand visualization']
  }

  // Connections
  create_dirs -> conda_env -> install_tools
  install_tools -> kofam_db
  install_tools -> dbcan_db
  install_tools -> merops_db
  install_tools -> iprscan_db

  kofam_db -> kofam
  dbcan_db -> dbcan
  merops_db -> merops
  iprscan_db -> interpro

  interpro -> rbims
  kofam -> rbims
  dbcan -> rbims
  merops -> rbims

  rbims -> analysis
}
")

# Visualize in RStudio
graph

# Optional: Export to SVG and PNG
#export_svg(graph) %>% charToRaw %>% rsvg_svg("workflow.svg")
#export_svg(graph) %>% charToRaw %>% rsvg_png("workflow.png")

export_svg(graph) %>% cat(file = "workflow.svg")
