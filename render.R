#!/usr/bin/env Rscript
# render.R -- Render the synchronized VTC manuscript and supplement
# Usage: Rscript render.R [main|supplement|both]

args <- commandArgs(trailingOnly = TRUE)
target <- if (length(args) == 0) "both" else tolower(args[1])

if (!target %in% c("main", "pdf", "supplement", "supp", "both")) {
  stop("Unknown render target: ", target,
       ". Use one of: main, pdf, supplement, supp, both.")
}

if (target %in% c("main", "pdf", "both")) {
  rmarkdown::render(
    "VTC_RTG_analysis.Rmd",
    output_file = "VTC_RTG_analysis.pdf"
  )
  message("Main PDF rendered.")
  rmarkdown::render(
    "VTC_RTG_analysis.Rmd",
    output_format = bookdown::word_document2(
      toc = FALSE,
      number_sections = FALSE
    ),
    output_file = "VTC_RTG_analysis.docx"
  )
  message("Main Word manuscript rendered.")
}

if (target %in% c("supplement", "supp", "both")) {
  rmarkdown::render("supplement.Rmd", output_file = "supplement.pdf")
  message("Supplement PDF rendered.")
}
