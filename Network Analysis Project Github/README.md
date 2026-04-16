# Brain Network Threshold Analysis

This project investigates how methodological choices, particularly thresholding, influence the interpretation of brain network organisation. It uses graph theoretical analysis to examine how brain functional connectivity networks change across different threshold levels, focusing on small-world organisation, hub structure, community structure, and connectivity patterns.

---
  
## Research Aim
  
To investigate how brain network organisation varies across thresholds by examining:
  
  - Small-world properties (clustering, path length, small-world index)
  - Hub structure (betweenness centrality)
  - Community structure (Louvain modularity)
  - Connectivity patterns (assortativity by degree, hemisphere, and lobe)

---
  
## Data
  
The dataset consists of:
  
  - `AveragedMatrix.txt`: Functional connectivity matrix between brain regions
  - `region_names_full_file.txt`: Brain region labels

Negative correlations were transformed using absolute values prior to analysis.

---
  
## Methods Overview
  
1. Load and preprocess connectivity matrix
2. Apply proportional thresholding (10%–50%)
3. Construct weighted undirected graphs using `igraph`
4. Compute network metrics at each threshold:
  - Small-world index (SWI)
  - Clustering coefficient
  - Characteristic path length
  - Betweenness centrality (hub detection)
  - Louvain modularity
  - Assortativity (degree, hemisphere, lobe)
5. Generate 1000 Erdős–Rényi random networks for baseline comparison
6. Compare real vs random network properties
7. Analyse how metrics change across thresholds

This approach allows evaluation of network organisation across multiple levels, including global structure, community structure, and node-level centrality.

---
  
## File Structure

```
project/
│
├── data/
│   ├── AveragedMatrix.txt
│   ├── region_names_full_file.txt
│   ├── region_names_abbrev_file.txt
│
├── scripts/
│   ├── analysis.R
│
├── results/
│   ├── Final Results.csv
│   ├── Final_Report.RData
│
├── figures/
│   ├── SWI_plot.png
│   ├── hub_plot.png
│
└── README.md
```

---

## How to Run

1. Open R or RStudio
2. Set working directory to project folder
3. Install required package:

```r
install.packages("igraph")
```

4. Run the analysis: 

```r
source("Final_Report_Script.R")
```

---

## Outputs
The script generates:
- `Final Results.csv`: Network metrics across thresholds
- `Final_Report.RData`: Full R workspace snapshot
- Figures showing:
  - Small-world index (SWI)
  - Hub dominance across thresholds
  - Modularity changes
  - Assortativity patterns

---

## Key Findings

- Network density decreases systematically with increasing threshold
- Small-world organisation is present but varies in strength across thresholds
- Hub dominance increases as networks become sparser
- Modularity increases with threshold, indicating increasing segregation
- Connectivity patterns vary across anatomical organisation

## Author

Chua Weng Yan

National University of Singapore
