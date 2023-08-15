# Setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
rm(list=ls())

# Load packages
library(deSolve)
library(ggplot2)
library(RColorBrewer)
library(cowplot)
library(writexl)


# Generating colourmap 
colour.aux <- brewer.pal(n = 8, name = 'Dark2')
generate_colourmap <- function(num_colours) {
  # Specify the range of colours
  colour_range <- c("purple", "blue", "green", "orange", "red")
  # Create the colourmap using colorRampPalette
  colourmap <- colorRampPalette(colour_range)(num_colours)
  return(colourmap)
}

# Defining dynamic of vaccine reducing probability of infection model ####
sis <- function(t,Y,p){
  # Parameters #################################################################
  lambda0 = p[1]        # Force of infection at age zero
  c       = p[2]        # Value at birth
  k       = p[3]        # Steepness of the force of infection Increments in early ages
  alpha   = p[4]        # Rate of immunity development after repeated infections
  gamma   = p[5]        # Recovery rate from clinical infection to susceptible 
  gamma2  = p[6]        # Recovery rate from asymptomatic to susceptible 
  p1      = p[7]        # Proportion of low-risk individuals
  p2      = p[8]        # proportion of high-risk individuals
  x1      = p[9]        # Risk factor of low-risk individuals
  x2      = p[10]       # Risk factor of high-risk individuals
  eta     = p[11]       # reduced factor for vaccinated individuals
  
  x        = c(x1,x2)                   # Risk factors for low-risk and high-risk individuals
  lambda   = lambda0*(1-c*exp(-k*t))    # Age-dependent force of infection
  lambda_v = eta*lambda                 # Age-dependent force of infection in the vaccinated group 
    
  # Defining compartments #####################################################
  aux = length(Y)/6
  
  # col1: Low-risk
  # col2: High-risk
  # col3: Low-risk - vaccinated
  # col4: High-risk - vaccinated  
  S = matrix(0,ncol = 2, nrow = aux)      # Susceptible
  I = matrix(0,ncol = 2, nrow = aux)      # Clinical Infection
  A = matrix(0,ncol = 2, nrow = aux)      # Asymptomatic infection
  
  # Low risk group
  S[,1] = Y[(0*aux+1):(1*aux)]      # Susceptible
  I[,1] = Y[(1*aux+1):(2*aux)]      # Clinical Infection
  A[,1] = Y[(2*aux+1):(3*aux)]      # Asymptomatic infection
  
  # High risk group
  S[,2] = Y[(3*aux+1):(4*aux)]      # Susceptible
  I[,2] = Y[(4*aux+1):(5*aux)]      # Clinical Infection
  A[,2] = Y[(5*aux+1):(6*aux)]      # Asymptomatic infection
  
  # dYdt: col1 for Low-risk - col2 for High-risk
  dYdt = matrix(0,ncol = 2, nrow = 3*aux)
  
  # Model Equations ############################################################
  # Low-risk group
  dYdt[0*aux+1,1] = -x[1]*lambda_v*S[1,1]                   # Sv0,1
  dYdt[1*aux+1,1] = +x[1]*lambda_v*S[1,1] - gamma*I[1,1]    # Iv1,1
  dYdt[2*aux+1,1] = 0                                       # Av0,1
  
  # High-risk group
  dYdt[0*aux+1,2] = -x[2]*lambda_v*S[1,2]                   # Sv0,2
  dYdt[1*aux+1,2] = +x[2]*lambda_v*S[1,2] - gamma*I[1,2]    # Iv1,2
  dYdt[2*aux+1,2] = 0                                       # Av0,2
  
  for (i in 2:aux) {
    # Low-risk group
    dYdt[0*aux+i,1] = -x[1]*lambda_v*S[i,1] + gamma*I[i-1,1] + gamma2*A[i,1]                                        # Svi,1
    dYdt[1*aux+i,1] =  x[1]*lambda_v*exp(-alpha*i)*(A[i,1]+S[i,1]) - gamma*I[i,1]                                   # Ivi,1
    dYdt[2*aux+i,1] =  x[1]*lambda_v*(1-exp(-alpha*i))*S[i,1] - x[1]*lambda_v*exp(-alpha*i)*A[i,1] - gamma2*A[i,1]  # Avi,1
    
    # High-risk group
    dYdt[0*aux+i,2] = -x[2]*lambda_v*S[i,2] + gamma*I[i-1,2] + gamma2*A[i,2]                                        # Svi,2
    dYdt[1*aux+i,2] =  x[2]*lambda_v*exp(-alpha*i)*(A[i,2]+S[i,2]) - gamma*I[i,2]                                   # Ivi,2
    dYdt[2*aux+i,2] =  x[2]*lambda_v*(1-exp(-alpha*i))*S[i,2] - x[2]*lambda_v*exp(-alpha*i)*A[i,2] - gamma2*A[i,2]  # Avi,2
  }
  
  return(list(dYdt))
}

# Parameters
p = c(
  lambda0 = 0.6197,     # Force of infection at age zero
  c       = 0.8720,     # Value at birth
  k       = 0.0493,     # Steepness of the force of infection Increments in early ages    
  alpha   = 0.0285,     # Rate of immunity development after repeated infections
  gamma   = 365/28,     # Recovery rate from clinical infection to susceptible
  gamma2  = 365/90,     # Recovery rate from asymptomatic to susceptible
  p1      = 0.8,        # Proportion of low-risk individuals
  p2      = 0.2,        # proportion of high-risk individuals
  x1      = 0.0883,     # Risk factor of low-risk individuals
  x2      = 4.6467,     # Risk factor of high-risk individuals
  eta     = 1.0         # Reduced factor for vaccinated individuals
)

p.vaccine = c(
  lambda0 = 0.6197,     # Force of infection at age zero
  c       = 0.8720,     # Value at birth
  k       = 0.0493,     # Steepness of the force of infection Increments in early ages
  alpha   = 0.0285,     # Rate of immunity development after repeated infections
  gamma   = 365/28,     # Recovery rate from clinical infection to susceptible 
  gamma2  = 365/90,     # Recovery rate from asymptomatic to susceptible
  p1      = 0.8,        # Proportion of low-risk individuals
  p2      = 0.2,        # proportion of high-risk individuals
  x1      = 0.0883,     # Risk factor of low-risk individuals
  x2      = 4.6467,     # Risk factor of high-risk individuals
  eta     = 0           # Reduced factor for vaccinated individuals
)

# Time length
t = seq(0,80,by=1)

# Initial condition (naive individuals are allocated to the compartment S0)
aux = 100                                         # 100 equations for each compartment
init = rep(0,6*aux)
init[c(00*aux+1,01*aux+1,02*aux+1)] = c(1,0,0)
init[c(03*aux+1,04*aux+1,05*aux+1)] = c(1,0,0)

# Solve system using lsoda
sol = lsoda(init,t,sis,p)
sol = as.data.frame(sol)

# Solution for vaccinated individuals
sol_v = sol[,1:ncol(sol)]
sol_v[,] = 0
sol_v[1,] = c(0,init)
for (i in 1:80) {
  t = c((i-1):(i))
  sol_v[i+1,] = lsoda(
    # Initial condition comes from the matrix sol
    y     = as.numeric(sol[i,2:ncol(sol)]),
    # Simulation for one year
    times = t,
    func  = sis,
    # Using p.vaccine parameters which account for reduced factor eta for vaccinated individuals
    p     = p.vaccine
  )[2,] 
}

# Number of cumulative cases by age - with no vaccination
cases.cum = rowSums(
  (
    p["p1"]*(sol[1:81,c((00*aux+1):(01*aux))+1]+sol[1:81,c((01*aux+1):(02*aux))+1]+sol[1:81,c((02*aux+1):(03*aux))+1]) + 
    p["p2"]*(sol[1:81,c((03*aux+1):(04*aux))+1]+sol[1:81,c((04*aux+1):(05*aux))+1]+sol[1:81,c((05*aux+1):(06*aux))+1])
  )*matrix(rep(c(0:(aux-1)), each = 81), nrow = 81, ncol = aux)
)

# Number of cumulative cases by age - with vaccination
cases.cum_v = rowSums(
  (
    p["p1"]*(sol_v[1:81,c((00*aux+1):(01*aux))+1]+sol_v[1:81,c((01*aux+1):(02*aux))+1]+sol_v[1:81,c((02*aux+1):(03*aux))+1]) + 
    p["p2"]*(sol_v[1:81,c((03*aux+1):(04*aux))+1]+sol_v[1:81,c((04*aux+1):(05*aux))+1]+sol_v[1:81,c((05*aux+1):(06*aux))+1])
  )*matrix(rep(c(0:(aux-1)), each = 81), nrow = 81, ncol = aux)
)


# The indirect effect of vaccination (assuming that all ages have the same number of individuals)
# Pop demography
age_distribution = c(
  0.016509165,0.020835705,0.020380280,0.020607993,0.017647729,0.019924855,0.021404987,0.025048389,
  0.025162245,0.020152567,0.022429694,0.023112832,0.023682113,0.025845383,0.025959239,0.024820676,
  0.027439372,0.027439372,0.024137538,0.019697142,0.022315837,0.017647729,0.018672435,0.014915177,
  0.015712171,0.015256746,0.017989298,0.015029033,0.014459752,0.014232039,0.018558579,0.015256746,
  0.017989298,0.017875441,0.014915177,0.017989298,0.015256746,0.011841057,0.010133212,0.011954913,
  0.010816350,0.010247068,0.012410338,0.013435045,0.010133212,0.007969942,0.007856086,0.006262097,
  0.007742229,0.006717522,0.006148241,0.005237390,0.007059091,0.007172948,0.006603666,0.006262097,
  0.005692816,0.005920528,0.006831379,0.004781965,0.006034385,0.005009678,0.003871115,0.004668109,
  0.003984971,0.004326540,0.003984971,0.003757258,0.003415689,0.003757258,0.004554253,0.001935557,
  0.003301833,0.003643402,0.003529546,0.003301833,0.002390983,0.002277126,0.002277126,0.002049414,
  0.002390983
)

# Rho
rho = 0.1 # relative contribution of asymptomatic infection to the population

# Infectious time
infectious.time1 = 4/28
infectious.time2 = 7/28
infectious.time3 = 12/28

# Proportion of individuals in compartments I and A
# The proportion of infectious individuals in the population accounting for the fact that asymptomatic individuals are less infectious
theta = sum(
  rowSums(
    p["p1"]*( # Low-risk group
      infectious.time1*sol[1:81,c((01*aux+1):(02*aux))+1]+   #I1,i
      rho*sol[1:81,c((02*aux+1):(03*aux))+1]                 #A1,i
    ) + 
    p["p2"]*( # High-risk group
      infectious.time1*sol[1:81,c((04*aux+1):(05*aux))+1]+   #I2,i
      rho*sol[1:81,c((05*aux+1):(06*aux))+1]                 #A2,i
    )
  )*age_distribution
)

# Proportion of individuals in compartments I and A assuming that all individuals were vaccinated
theta_v = sum(
  rowSums(
    p["p1"]*( # Low-risk group
      infectious.time1*sol_v[1:81,c((01*aux+1):(02*aux))+1]+   #Iv1,i
      rho*sol_v[1:81,c((02*aux+1):(03*aux))+1]                 #Av1,i
    ) +
    p["p2"]*( # High-risk group
      infectious.time1*sol_v[1:81,c((04*aux+1):(05*aux))+1]+   #Iv2,i
      rho*sol_v[1:81,c((05*aux+1):(06*aux))+1]                 #Av2,i
    )
  )*age_distribution
)

# proportion of people vaccinated
p_v   = 1
sim<-c(0.2, 0.4, 0.6, 0.8)

# Specify the number of colors in the coluormap
num_colours <- length(sim)

# Generate the colourmap using the function
line_colours <- generate_colourmap(num_colours)
iter_line = 1

for(k in sim) {
eta_p_fullv = ( theta_v*k + theta*(1-k) ) / theta
eta_p_fullv

# Age min and max for vaccinated individuals
age.vac.min = 18 - k*10/2
age.vac.max = 35 + k*10/2
aux.age.vac = (age.vac.min+1):(age.vac.max+1)

theta_v_age.vac = sum(
  rowSums(
    p["p1"]*( # Low-risk group
      infectious.time1*
        rbind(k*sol_v[ aux.age.vac,   c((01*aux+1):(02*aux))]+
          (1-k)*sol  [ aux.age.vac,   c((01*aux+1):(02*aux))],
                  sol  [-aux.age.vac,   c((01*aux+1):(02*aux))])+   #c(Iv1,i;I1,i)
    rho*rbind(k*sol_v[ aux.age.vac,   c((02*aux+1):(03*aux))]+
          (1-k)*sol  [ aux.age.vac,   c((02*aux+1):(03*aux))],      
                  sol  [-aux.age.vac,   c((02*aux+1):(03*aux))])    #c(Av1,i;A1,i)
    ) +
      p["p2"]*( # High-risk group
        infectious.time1*
          rbind(k*sol_v[ aux.age.vac, c((04*aux+1):(05*aux))]+
            (1-k)*sol  [ aux.age.vac, c((04*aux+1):(05*aux))],
                    sol  [-aux.age.vac, c((04*aux+1):(05*aux))])+   #c(Iv2,i;I2,i)
      rho*rbind(k*sol_v[ aux.age.vac, c((05*aux+1):(06*aux))]+
            (1-k)*sol  [ aux.age.vac, c((05*aux+1):(06*aux))],      
                    sol  [-aux.age.vac, c((05*aux+1):(06*aux))])    #c(Av2,i;A2,i)
      )
  )*age_distribution[c(aux.age.vac,setdiff(1:81, aux.age.vac))]
)

# Indirect effect of the vaccination in the unvaccinated population
eta_p = ( theta_v_age.vac ) / theta
eta_p

# simulating the indirect effect of vaccination on the non vaccinated population
p.nvac.ind = c(
  lambda0 = 0.6197*eta_p, # Force of infection at age zero
  c       = 0.8720,       # Value at birth
  k       = 0.0493,       # Steepness of the force of infection Increments in early ages
  alpha   = 0.0285,       # Rate of immunity development after repeated infections
  gamma   = 365/28,       # Recovery rate from clinical infection to susceptible 
  gamma2  = 365/90,       # Recovery rate from asymptomatic to susceptible 
  p1      = 0.8,          # Proportion of low-risk individuals
  p2      = 0.2,          # proportion of high-risk individuals
  x1      = 0.0883,       # Risk factor of low-risk individuals
  x2      = 4.6467,       # Risk factor of high-risk individuals
  eta     = 1.0           # Reduced factor for vaccinated individuals
)

# Time length
t = seq(0,80,by=1)

# Solve system
sol_p = sol[,1:ncol(sol)]
sol_p[,] = 0
sol_p[1,] = c(0,init)
for (i in 1:80) {
  t = c((i-1):(i))
  sol_p[i+1,] = lsoda(
    # Initial condition comes from the matrix sol
    y     = as.numeric(sol[i,2:ncol(sol)]),
    # Simulation for one year
    times = t,
    func  = sis,
    # Using p.vaccine parameters which account for reduced factor eta for vaccinated individuals
    p     = p.nvac.ind
  )[2,] 
}

# Number of cumulative cases by age - with indirect effect of vaccination on the non vaccinated population
cases.cum_p = rowSums(
  (
    # Vaccinated
    p["p1"]*(sol_p[1:81,c((00*aux+1):(01*aux))+1]+sol_p[1:81,c((01*aux+1):(02*aux))+1]+sol_p[1:81,c((02*aux+1):(03*aux))+1]) + 
    p["p2"]*(sol_p[1:81,c((03*aux+1):(04*aux))+1]+sol_p[1:81,c((04*aux+1):(05*aux))+1]+sol_p[1:81,c((05*aux+1):(06*aux))+1])
  )*matrix(rep(c(0:(aux-1)), each = 81), nrow = 81, ncol = aux)
)

# Final matrix
sol_f = cbind( 
  time=c(aux.age.vac,setdiff(1:81, aux.age.vac)),
  # Low-risk group
  rbind(k*sol_v[ aux.age.vac,c((00*aux+1):(03*aux))+1]+
    (1-k)*sol_p[ aux.age.vac,c((00*aux+1):(03*aux))+1],
            sol_p[-aux.age.vac,c((00*aux+1):(03*aux))+1]),   #c(Sv1,i:A1,i)
  # High-risk group
  rbind(k*sol_v[ aux.age.vac,c((03*aux+1):(06*aux))+1]+
    (1-k)*sol_p[ aux.age.vac,c((03*aux+1):(06*aux))+1],
            sol_p[-aux.age.vac,c((03*aux+1):(06*aux))+1])   #c(Sv1,i:A1,i)  
)

sol_f = sol_f[order(sol_f$time),]

#Final incidence
cases.cum_f = rowSums(
  (
    # Indirect effect of non vaccinated
    p["p1"]*(sol_f[1:81,c((00*aux+1):(01*aux))+1]+sol_f[1:81,c((01*aux+1):(02*aux))+1]+sol_f[1:81,c((02*aux+1):(03*aux))+1]) + 
    p["p2"]*(sol_f[1:81,c((03*aux+1):(04*aux))+1]+sol_f[1:81,c((04*aux+1):(05*aux))+1]+sol_f[1:81,c((05*aux+1):(06*aux))+1])
  )*matrix(rep(c(0:(aux-1)), each = 81), nrow = 81, ncol = aux)
)

#Cumulative cases
cum.sum.aux = cbind(
  cases.cum,
  cases.cum_v,
  cases.cum_p,
  cases.cum_f
)

# Incidence stratified by age
incidence = data.frame(
  age            = c(1:80),
  no.vacination  = diff(cases.cum),
  vacination     = abs(cases.cum_v[2:81] - cases.cum[1:80]),
  no.vacination.indirect = (cases.cum_p[2:81] - cases.cum[1:80]),
  final.case     = abs(cases.cum_f[2:81] - cases.cum[1:80])
)

# Plots
if (k==sim[1]){
# par(mfrow=c(2,1))

FILENAME = 'out10/model1_adult35_eta_00'
FILENAME_IMG = paste(FILENAME, '.png', sep="")
FILENAME_EXCL = paste(FILENAME, '.xlsx', sep="")

png(file=FILENAME_IMG, width = 20, height = 15, units = 'cm', res = 300)
plot(
   incidence$age,
   incidence$vacination,
   xlab ="Age (years)",
   ylab = "Incidence",
   ylim = c(0,0.3),
   type = "l",
   main = " ",
   lwd = 1.75,
   lty = 3
 )
  
lines(
  incidence$age,
  incidence$no.vacination,
  lwd = 1.75,
  lty = 2
)
}

lines(
  incidence$age,
  incidence$final.case,
  col = line_colours[iter_line]
)

iter_line = iter_line + 1

}

# Setting plot attribute
legend_name = c( 'No Vaccination',  '20%' , '40%', '60%', '80%', 'Full Vaccination')
line_colours = append(line_colours, 'black', after=0)
line_colours = append(line_colours, 'black', after=length(line_colours))
legend("topright", legend = legend_name, col = line_colours, lty = c(3,1,1,1,1,2), cex = 1, bty = "n")
mtext("B", side = 3, line =-2, adj = 0, cex=1.5, outer= TRUE)
dev.off()

 
# CHILDREN
# incidence_baseline = sum(incidence$no.vacination * age_distribution[1:80]) / sum(age_distribution[1:80])
# incidence_final_case = sum(incidence$no.vacination.indirect[19:80] * age_distribution[19:80]) +  sum(incidence$vacination[1:18] * age_distribution[1:18])
# incidence_vacination = (sum(incidence$vacination[1:18] * age_distribution[1:18])) / sum(age_distribution[1:18])
# incidence_vacination_indirect = (sum(incidence$no.vacination.indirect[19:80] * age_distribution[19:80])) / (sum(age_distribution[19:80]))
# incidence_full_vac = (sum(incidence$vacination[1:80] * age_distribution[1:80])) / sum(age_distribution[1:80])

# df_out = data.frame(
#   inc_vac_indirect = incidence_vacination_indirect,
#   inc_vac = incidence_vacination,
#   inc_final_case = incidence_final_case,
#   reduction_final_case = (1-(incidence_final_case/incidence_baseline)) * 100,
#   reduction_direct = (1-(incidence_vacination/incidence_baseline)) * 100,
#   reduction_indirect = (1-(incidence_vacination_indirect/incidence_baseline)) * 100,
#   reduction_full_vac = (1-(incidence_full_vac/incidence_baseline)) * 100
# )

# ADULT 
incidence_baseline = sum(incidence$no.vacination * age_distribution[1:80]) / sum(age_distribution[1:80])
incidence_final_case = (sum(incidence$no.vacination.indirect[1:18] * age_distribution[1:18]) + sum(incidence$no.vacination.indirect[37:80] * age_distribution[37:80])) + (sum(incidence$vacination[19:36] * age_distribution[19:36]))
incidence_vacination = (sum(incidence$vacination[19:36] * age_distribution[19:36])) / sum(age_distribution[19:36])
incidence_vacination_indirect = (sum(incidence$no.vacination.indirect[1:18] * age_distribution[1:18]) + sum(incidence$no.vacination.indirect[37:80] * age_distribution[37:80]))/(sum(age_distribution[1:18]) + sum(age_distribution[37:80]))
incidence_full_vac = (sum(incidence$vacination[1:80] * age_distribution[1:80])) / sum(age_distribution[1:80])

df_out = data.frame(
  inc_vac_indirect= incidence_vacination_indirect,
  inc_vac = incidence_vacination ,
  inc_final_case =incidence_final_case,
  reduction_final_case = (1-(incidence_final_case/incidence_baseline)) * 100,
  reduction_direct = (1-(incidence_vacination/incidence_baseline)) * 100,
  reduction_indirect = (1-(incidence_vacination_indirect/incidence_baseline)) * 100,
  reduction_full_vac = (1-(incidence_full_vac/incidence_baseline)) * 100
)

# Saving output file
colnames(df_out) = c('inc_vac_indirect', 'inc_vac', 'inc_final_case', 'reduction_final_case', 'reduction_direct', 'reduction_indirect', 'reduction_full_vac')
write_xlsx(df_out, FILENAME_EXCL)