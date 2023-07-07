# Script to compile all Rmd files

library(rmarkdown)

#part 1
render(input = "./analysis/manuscript_overview.Rmd", output_format = "html_document", output_dir = "results/")

#part 2
render(input = "./analysis/manuscript_drugResponses.Rmd", output_format = "html_document", output_dir = "results/")

#part 3
render(input = "./analysis/manuscript_clinicalAnalysis.Rmd", output_format = "html_document", output_dir = "results/")