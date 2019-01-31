
# Stochastic Simulation of Antibiotic Resistance [S.A.R.]

**SAR** is a simulation package originally written to answer questions addressed in the manuscript titled *Genetic basis and patterns of gaining antibiotic resistance of a naïve strain of Staphylococcus aureus*, authored by S., A., and G. Abdel-Azim. While the antibiotic resistance phenomenon is well-documented, the pattern by which bacteria gain resistance remains unknown. The study analyzed the pattern of amoxicillin resistance in Staphylococcus aureus and, through the use of simulation, offered explanations for that pattern. In the lab, we found that bacteria gain resistance to antibiotics in a stepwise pattern. The simulation package was useful in providing evidence of a step-wise pattern caused by varying levels of complexity in the mutations that bacteria undergo to gain resistance. The various cell processes involved in the evolution of bacterial populations, such as population construction, population fitness adjustment, cell division, and mutation, were simulated to ascertain the trend of antibiotic resistance in bacteria. The ability to set unique initial bacterial population parameters provides flexibility in its usage. A Shiny app for SAR is available at https://github.com/ahmadazim/SAR-Shiny-App.

# To install SAR directly from github:
First install devtools and Rtools:
- Install Rtools from https://cran.r-project.org/bin/windows/Rtools/. This is a windows executable program that you install by clicking on the link and following instructions.
- To install devtools: ```install.packages("devtools")```\
Then, run the following in R
```r
library(devtools)
install_github("SalmaGAA/SAR")
library(SAR)
```

# SAR Quickstart: How to Use *runSimulation()* to Study Bacterial Populations Under Antibiotic Stress
The *runSimulation* function was developed to allow for the investigation of the pattern of fitness development and population growth of bacteria cultured over multiple exposures to antibiotic stress. As mutagenic compounds, antibiotics induce mutations within a bacterial population, which leads to a rapid increase in the fitness of such a population. The complexity of a bacterial genome can be varied, by adding mutations sites to several genes, to investigate its effect on bacterial population behaviors and resistance patterns. The developments of fitness in such simple and complex regions of bacterial genomes can also be individually simulated with the *runSimulation* function.\
\
Bacteria are simulated to be cultured for 24 hours under environmental stress. A threshold argument is included in the simulation to represent the antibiotic stress placed on a bacterial culture. Bacteria that do not meet this threshold (i.e. bacteria that have not mutated enough in each generation) are barred from proliferation, representing antibiotic stress. This process allows for overall population fitness to increase as only the most mutated bacteria are selected to survive to subsequent generations.\
\
**Example of *runSimulation* function:**
```r
# Run function with the following arguments,
runSimulation(Ng = c(2, 10, 15),
              Nl = c(1, 6, 10),
              gen.interval = 480,
              Rm = 0.1,
              Psize = 300,
              startingFitness = 0.60,
              thr = 0.51,
              nDays = 2,
              maxPsize = 2000)

# Function Output is,

  Finished gen:  1 of day 1 : 358 0.6028395 0.6333333 0.6013333 0.5997778
  Finished gen:  2 of day 1 : 428 0.6059383 0.674581 0.601676 0.5996276
  Finished gen:  3 of day 1 : 538 0.6089477 0.7219626 0.5950935 0.6031153
  Finished day 1

  Finished gen:  1 of day 2 : 366 0.6104938 0.765 0.5866667 0.6057778
  Finished gen:  2 of day 2 : 440 0.6134386 0.7868852 0.5939891 0.6032787
  Finished gen:  3 of day 2 : 564 0.6159091 0.8056818 0.5920455 0.6065152
  Finished day 2

       [,1]      [,2]      [,3]      [,4]      [,5]
  [1,]  358 0.6028395 0.6333333 0.6013333 0.5997778
  [2,]  428 0.6059383 0.6745810 0.6016760 0.5996276
  [3,]  538 0.6089477 0.7219626 0.5950935 0.6031153
  [4,]  366 0.6104938 0.7650000 0.5866667 0.6057778
  [5,]  440 0.6134386 0.7868852 0.5939891 0.6032787
  [6,]  564 0.6159091 0.8056818 0.5920455 0.6065152
```
As the simulation progresses through bacterial generations, its status will be printed. The population size, the fitness of each region, and the overall fitness of the population will be printed following each generation. After completion of the simulation, a matrix will be outputted, where the first column represents the population size, the second column represents overall population fitness, and subsequent columns will represent section fitness in their order of input. Each row represents a generation simulated.
