---
output: html_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

## Binning functions: 
# 1d binning
- Binning by definition (bounds and centers)
- Creating bin definition lists from vectors of bounds and centers
- General rectangular binning
- Standard rectangular binning
- Quantile rectangular binning
- Random rectangular binning
- Frequency binning of reduced binned data

# 2d binning
- Independant 1d binning grid
- Iterative conditional binning
- Reduced binned data - one row per x/y bin

## Binning Loss Functions:
- Spatial loss summary
- Frequency loss summary

## Binned Scatterplot Functions
- Binned scatterplot construction from reduced binned
- Binned scatterplot construction from raw data and bin specs

