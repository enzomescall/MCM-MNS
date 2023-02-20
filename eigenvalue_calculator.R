# Eigenvalue calculation
gamma = 1     # speed of poaching adjustment
r_p =  0.05   # revenue per poach
Omega =  0.6  # opportunity cost
sigma_p = 0.01# risk of poaching
lambda_0 = 5  # base law enforcement level
lambda = 35   # law enforcement per investment
c = 0.375     # culling rate per predator
B = 0.5       # boma construction per unit invested
E_r = 0.5     # baseline recovery rate
I_0 = 0.5     # baseline government investment
P_t = 0.3     # average pollution per tourist
I_t = 0.1     # investment into ecology per tourism income
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
delta = 0.01  # effect of degradation on tourism
A = 0.01      # tourism destination attractiveness saturation

tau = (E_r + I_0 - G)/(P_t - I_t)
mu = tau*B/c
omega = A/(delta*tau)
pi = (Omega + sigma_p * (lambda_0 + lambda * tau))/r_p
upsilon = (h_1 * pi - d - S * tau)/beta
iota = ((b * (1 - pi / (k - eta * omega)) - h_0 * mu) * pi + m_0 - m * omega)/(alpha * pi)

c(Prey = pi,Pred = mu,Tour = tau, Deg =  omega, Poach = iota, RK = upsilon)

TAU = b - 2*pi/(k-eta*omega) - h_0*mu - alpha*iota
PI = h_1*pi - d - beta*upsilon - S*tau
MU = gamma*r_p*pi - gamma*(Omega+sigma_p*lambda_0) - gamma*sigma_p*lambda*tau


J = matrix(c(
            TAU, -h_0*pi, 0, eta*pi**2/((k-eta*omega)**2)-m, -alpha*pi, 0,
            h_1*mu, PI, -S*mu, 0, 0, -beta*mu,
            0, 0, -A/(tau**2), -delta, 0, 0,
            0, 0, P_t - I_t, 0, 0, 0,
            gamma*r_p*tau, 0, -gamma*r_p*lambda*iota, 0, MU, 0,
            0, c, -B, 0, 0, 0),
            6, 6, byrow = TRUE)

Re(eigen(J)$values)


vars = c(gamma, r_p, Omega, sigma_p, lambda_0, lambda, c, B, E_r,
         I_0, P_t, I_t, G, b, k, eta, h_0, alpha, m_0, m, h_1, d, beta, S)

for (m in seq(0.01, 100, 0.01)) {
  tau = (E_r + I_0 - G)/(P_t - I_t)
  mu = tau*B/c
  omega = A/(delta*tau)
  pi = (Omega + sigma_p * (lambda_0 + lambda * tau))/r_p
  upsilon = (h_1 * pi - d - S * tau)/beta
  iota = ((b * (1 - pi / (k - eta * omega)) - h_0 * mu) * pi + m_0 - m * omega)/(alpha * pi)
  
  c(Prey = pi,Pred = mu,Tour = tau, Deg =  omega, Poach = iota, RK = upsilon)
  
  TAU = b - 2*pi/(k-eta*omega) - h_0*mu - alpha*iota
  PI = h_1*pi - d - beta*upsilon - S*tau
  MU = gamma*r_p*pi - gamma*(Omega+sigma_p*lambda_0) - gamma*sigma_p*lambda*tau
  
  
  J = matrix(c(
    TAU, -h_0*pi, 0, eta*pi**2/((k-eta*omega)**2)-m, -alpha*pi, 0,
    h_1*mu, PI, -S*mu, 0, 0, -beta*mu,
    0, 0, -A/(tau**2), -delta, 0, 0,
    0, 0, P_t - I_t, 0, 0, 0,
    gamma*r_p*tau, 0, -gamma*r_p*lambda*iota, 0, MU, 0,
    0, c, -B, 0, 0, 0),
    6, 6, byrow = TRUE)
  
  if (any(Re(eigen(J)$values) > 0)) {
    print(c(m = m, Eigenvalues = which(Re(eigen(J)$values) > 0)))
  }
}


