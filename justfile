

site:
  quarto render --to html

cleantmp:
  find tmp -delete
  mkdir tmp


sim01 cfg:
  Rscript --vanilla ./R/sim01.R run_sim01 {{cfg}} 
  just cleantmp

runsim:
  just sim01 ../etc/sim01/cfg-sim01-v01.yml
  just sim01 ../etc/sim01/cfg-sim01-v02.yml
  just sim01 ../etc/sim01/cfg-sim01-v03.yml
  just sim01 ../etc/sim01/cfg-sim01-v04.yml
  just sim01 ../etc/sim01/cfg-sim01-v05.yml

sap:
  quarto render roadmap-sap.qmd --to pdf

notes:
  quarto render roadmap-notes.qmd --to pdf

sim-report:
  quarto render reports/sim01-report.qmd --to pdf

all: sap notes sim-report



