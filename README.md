
# Stochastic Simulation of Antibiotic Resistance [S.A.R.]

**SAR** is a simulation package originally written to answer questions addressed in the manuscript titled *Genetic basis and patterns of gaining antibiotic resistance of a naïve strain of Staphylococcus aureus*, authored by S., A., and G. Abdel-Azim. While the antibiotic resistance phenomenon is well-documented, the pattern by which bacteria gain resistance remains unknown. The study analyzed the pattern of amoxicillin resistance in Staphylococcus aureus and, through the use of simulation, offered explanations for that pattern. In the lab, we found that bacteria gain resistance to antibiotics in a stepwise pattern. The simulation package was useful in providing evidence of a step-wise pattern caused by varying levels of complexity in the mutations that bacteria undergo to gain resistance. The various cell processes involved in the evolution of bacterial populations, such as population construction, population fitness adjustment, cell division, and mutation, were simulated to ascertain the trend of antibiotic resistance in bacteria. The ability to set unique initial bacterial population parameters provides flexibility in its usage. A Shiny app for SAR is available at https://github.com/ahmadazim/SAR-Shiny-App.

# To install SAR directly from github:
  1. 	Make sure you have devtools and Rtools installed … or install them by:\
    a.  install.packages("devtools")  # in your R prompt\
    b.	install rtools from https://cran.r-project.org/bin/windows/Rtools/  # this is a windows executable program that you install by clicking and following instructions\
  2.	library(devtools)
  3.	install_github("SalmaGAA/SAR")
  4.	library(SAR)

# Using the runSimulation function:
The *runSimulation* function was developed to allow for the investigation of the pattern of fitness development and population growth of bacteria cultured over multiple exposures to antibiotic stress. As mutagenic compounds, antibiotics induce mutations within a bacterial population, which leads to a rapid increase in the fitness of such a population. The complexity of a bacterial genome can be varied, by adding mutations sites to several genes, to investigate its effect on bacterial population behaviors and resistance patterns. The developments of fitness in such simple and complex regions of bacterial genomes can also be individually simulated with the *runSimulation* function.\
Bacteria are simulated to be cultured for 24 hours under environmental stress. A threshold argument is included in the simulation to represent the antibiotic stress placed on a bacterial culture. Bacteria that do not meet this threshold (i.e. bacteria that have not mutated enough in each generation) are barred from proliferation, representing antibiotic stress. This process allows for overall population fitness to increase as only the most mutated bacteria are selected to survive to subsequent generations.
