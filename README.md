# KeepYourKidsAtHome_IBM
This is a spatially implicit individual-based SIR model for a large range of host and pathogen life-history strategies. We focused specifically on age-dependent between-group interactions, a life history strategy referred as “keep kids at home” (KKH), which occurs in social species rearing their offspring within colonies or communal nurseries, in which juveniles do not engage in between-group contact. We show that such KKH strategy has great influence on the dynamics of diseases, which confer lifelong immunity in group-living species. 

The purpose of our model is to study how disease dynamics in group living species are affected by “keeping their kids at home” (KKH). We are particularly interested in whether such age-dependent patterns in the between-group interactions increase the probability of epidemic fade-out consistently across a range of hosts’ and pathogens’ life-history traits.


##**Entities**
Each individual is characterized by two state variables, which are updated on a weekly basis: age [weeks] and epidemiological status (susceptible S, infected I or recovered R, i.e. immune)

##**Process**
At each time step, the following sequence of processes occurs: (1) pathogen transmission and restructuring of the network, (2) death related to the intrinsic mortality rate and to pathogen infection, (3) reproduction, (4) dispersal of individuals and (5) ageing. The details of the processes are as follows: (1) We considered the two classic forms of pathogen transmission (parameter transmission in Table 1): density-dependent transmission, where the spread of infections depends on the size of the groups, and frequency-dependent transmission where it does not. 
