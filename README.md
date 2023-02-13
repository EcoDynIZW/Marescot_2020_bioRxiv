# Marescot et al. (2020) bioRxiv

> Lucile Marescot, Mathias Franz, Sarah Benhaiem, Heribert Hofer, Marion L. East, Stephanie Kramer-Schadt (2020) ‘Keeping the kids at home’ can limit the persistence of contagious pathogens in social animals. *bioRxiv*. DOI: [10.1101/2020.04.11.036806](https://doi.org/10.1101/2020.04.11.036806)

This spatially implicit individual-based SIR model developped for a large range of host and pathogen life-history strategies. We focused specifically on offspring with restricted between-group contact (ORC), a life history traits, which occurs in social species rearing their offspring within communal nurseries, in which juveniles do not engage in between-group interactions before a certain age. We show that such ORC strategy has great influence on the dynamics of diseases, which confer lifelong immunity in group-living species. 


The purpose of our model is to study how disease dynamics in group living species that keep their kids at home. We are particularly interested in whether such age-dependent patterns in the between-group interactions in the ORC scenario increase the probability of epidemic fade-out consistently across a range of hosts’ and pathogens’ life-history traits.


## **Entities**
Each individual is characterized by two state variables, which are updated on a weekly basis: age [weeks] and epidemiological status (susceptible S, infected I or recovered R, i.e. immune)

## **Process**
At each time step, the following sequence of processes occurs: (1) pathogen transmission and restructuring of the network, (2) death related to the intrinsic mortality rate and to pathogen infection, (3) reproduction, (4) dispersal of individuals and (5) ageing. The details of the processes are as follows: (1) We considered the two classic forms of pathogen transmission (parameter transmission in Table 1): density-dependent transmission, where the spread of infections depends on the size of the groups, and frequency-dependent transmission where it does not. 


In order to assess how consistent the odds of epidemic fadeouts increase due to ORC, we systematically varied the following parameters that can influence pathogen transmission dynamics: (1) host birth rate, (2) host death rate (that occurs independently of infection), (3) contact rates between members of different groups, (4) infection rate at contact, (5) infection length, (6) virulence (i.e. the increase in host death rate due to infection), and (7) whether transmission is frequency-dependent or density-dependent. For each of the specific parameter combinations we investigated two scenarios. First, we set up a baseline scenario in which rates of between-group contact are identical for individuals of all ages. Second, we investigated a scenario of KKH by preventing individuals below an age-threshold (two years) to contact members of other groups, i.e. their respective rates of between-group contact were set to zero. 

Our model uses life-history parameters related to spotted hyena (Crocuta crocuta), such age at first between-group contact (nursing period in the communal den), age of first reproduction and number of offspring (Marescot et al. 2018), but is otherwise kept very simple to yield generalizable results. Generation time and reproductive output as well as time scales can easily be adapted to other species.
