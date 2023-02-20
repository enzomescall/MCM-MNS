# Utility function
gamma = 1     # speed of poaching adjustment
r_p =  0.05   # revenue per poach
sigma_p = 0.01# risk of poaching
c = 0.375     # culling rate per predator
E_r = 0.5     # baseline recovery rate
G = 0.4       # impact of global warming on degradation
b = 0.3       # birth rate of prey
k = 120       # carrying capacity of prey
eta = 1       # sensitivity of degradation on carrying cap
h_0 = 0.05    # effect of hunting on prey population
alpha = 0.01  # poaching effectiveness
m_0 = 1.5     # baseline migration inwards
m = 0.5       # effect of degradation on outwards migration
h_1 = 0.05    # conversion of prey biomass into predator biomass
d = 0.17      # death rate of predators
beta = 0.1    # effectiveness of retaliatory killings
S = 0.01      # sensitivity of predators to tourist distractions
delta = 0.01  # effect of degradation on Tour

# Changeable parameters

# utility function for policy generation
utility = function(x) {
  
  # weights for each of the variables
  a0 = 10000 # one unit in prey represents 100,000 animals
  a1 = 100 # one unit in predator represents 1,000 animals
  a2 = 10000 # one unit in tourism represents 100,000 yearly visitors
  a3 = 0 # chose to measure poaching through prey population
  a4 = 10000 # high weight on retaliatory killings to keep maasai happy 
  a5 = 1000 # deg
  
  cost_multiplier = 100000 # weight of the cost, assuming millions of USD
  
  B = x[1]
  Omega =  x[2]  
  lambda_0 = x[3]
  lambda = x[4] 
  I_0 = x[5]
  PI_t = x[6]
  A = x[7]
  
  Tour = (E_r + I_0 - G)/(PI_t)
  
  if (is.infinite(Tour)) {print(c(E_r, I_0, G, PI_t))}
  
  Pred = Tour * B/c
  
  Deg = A/(delta*Tour)
  
  Prey = (Omega + sigma_p * (lambda_0 + lambda * Tour))/r_p
  
  RK =  (h_1 * Prey - d - S * Tour)/beta
  
  Poach = ((b * (1 - Prey / (k - eta * Deg)) - h_0 * Pred) * Prey + m_0 - m * Deg)/(alpha * Prey)
  
  cost = B**2 + Omega**2 + lambda_0**2/5 + lambda**2/35 + PI_t**2 + A**2
  
ret = a0 * Prey + a1 * Pred + a2 * Tour - a3 * Poach - a4 * RK - a5 * Deg - cost_multiplier * cost
  
  if (is.nan(ret)) {
    print(c(Prey, Pred, Tour, Deg, Poach, RK))
  }
  
  return(-ret)
}

initial_vals = c(B = 0.5,   # B
                 Omeega = 0.6,   # Omega
                 lambda_0 = 5,     # lambda_0
                 lambda = 35,    # lambda
                 I_0 = 0.5,   # I_0
                 P_tminusI_t = 0.2,   # P_t - I_t
                 A = 0.01)  # A

upper = c(1,10,10,20,5,0.99,1)
lower = c(0.01,0.01,0.01,0.01,0.001,0.01)


utility(initial_vals)

result = optim(par = initial_vals, fn = utility, method = "L-BFGS-B", lower = lower, upper = upper)

result$par
