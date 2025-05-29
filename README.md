# PMconv: Proteometabolomic Network Analysis Tool

<p align="center">
  <img src="https://github.com/Ministreliya131/PMconv/raw/main/data/pmconv.png" alt="PMconv Logo" width="300"/>
</p>

## Overview
PMconv is a web-based application that enables bidirectional conversion between protein and metabolite identifiers based on shared metabolic pathways. The tool facilitates:

- Cross-omics comparisons between proteomic and metabolomic datasets
- Expansion of experimental molecular lists through pathway associations
- Interactive visualization of protein-metabolite interaction networks

## Key Features

### 🔄 Bidirectional Conversion
- Convert between UniProt (proteins) and HMDB (metabolites) identifiers
- Map multiple relationships between metabolites and proteins

### ⚙️ Customizable Analysis
- **Filtering options**:
  - Group by molecular formula or compound class
  - Specify endogenous vs. environmental sources
  - Adjust association thresholds
- **Supported identifiers**:
  - Proteins: UniProt (human-focused)
  - Metabolites: HMDB

### 🌐 Network Visualization
- **Interactive graph**:
  - Proteins: red circles
  - Metabolites: blue squares
- **Navigation controls**: zoom, pan, node selection
- **STRING integration**:
  - High-confidence PPIs (score ≥ 0.9)
  - Overlay of experimental metabolites

### 📁 Data Management
- Load example datasets
- Export results as CSV files

## Quick Start

### Web Access
Available at: [PMconv Web Application](https://pmconv.streamlit.app)
