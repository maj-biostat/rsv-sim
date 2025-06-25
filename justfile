

site:
  quarto render --to html

cleantmp:
  find tmp -delete
  mkdir tmp

sim01 cfg:
  Rscript --vanilla ./R/sim01.R run_sim01 {{cfg}} 

sap:
  quarto render roadmap-sap.qmd --to pdf

notes:
  quarto render roadmap-notes.qmd --to pdf

sim-report:
  quarto render reports/sim01-report.qmd --to pdf

all: sap notes sim-report



