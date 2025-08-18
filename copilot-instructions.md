# Copilot Instructions

## Purpose
You are assisting in an R/Tidyverse workflow. Prefer clear, reproducible,
production-ready code with minimal side effects.

## General Style
- Use tidyverse verbs (dplyr, tidyr, stringr, purrr, forcats) and piping.
- Keep objects snake_case; avoid overwriting base names.
- Put all constants and parameters near the top of a script section.
- **Do not include library() calls** in R code snippets — packages are preloaded at startup.

## Commenting & Indentation (MANDATORY)
Use this exact hierarchy in all R code you generate:

- `#!` Issues/bookmarks/important notes
- `#*` Major sections (numbered: `#* 1:`, `#* 2:` …)
- `#+` Subsections (numbered: `#+ 1.1:`, `#+ 2.1:` …)
- `#-` Sub-subsections (numbered: `#- 1.1.1:`, `#- 2.1.1:` …)
- `#_` Checks/verification steps

Indentation rules:
- 2 spaces per hierarchy level, except for Major sections and Subsections (`#*` and `#+`)
- Major sections: 0 spaces; subsections: +0; sub-subsections: +2
- Code blocks match their comment level

**Example**
```r
#* 1: Load data
  df <- readr::read_csv("data.csv")

#+ 1.1: Clean names
  df <- janitor::clean_names(df)

  #- 1.1.1: Verify rows
    #_ Expect > 100 rows
    stopifnot(nrow(df) > 100)
```

## Special Commands
- **CC** (Comment Clean): If I enter 'CC' in the chat, automatically clean and update the entire code/document to conform exactly to the above comment hierarchy and indentation rules.
