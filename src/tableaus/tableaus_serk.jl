"Tableau for the explicit 1-stage stochastic Euler method"
function TableauStochasticEuler()
  a = zeros(Float64, 1, 1)
  b = [1.0]
  c = [0.0]

  TableauSERK(:Stochastic_Euler_explicit_method, 1, a, b, c, 1, a, b, c)
end


"Tableau for the explicit 2-stage stochastic Heun method"
function TableauStochasticHeun()

    a = [[0.0 0.0]
         [1.0 0.0]]
    b = [0.5, 0.5]
    c = [0.,  1.]

    TableauSERK(:Stochastic_Heun_explicit_method, 2, a, b, c, 2, a, b, c)
end


"""
Tableau for the explicit Platen method
  Platen's method cited in Eq. (52) in
   K. Burrage, P. Burrage (1996) "High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations".
  According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion.
  Appears to have a rather poor long-time performance.
"""
function TableauPlaten()

    a_drift = [[0.0 0.0]
               [1.0 0.0]]
    b_drift = [1.0, 0.0]
    c_drift = [0.0, 1.0]

    a_diff = [[0.0 0.0]
              [1.0 0.0]]
    b_diff = [0.5, 0.5]
    c_diff = [0.0, 1.0]

    TableauSERK(:Platen_explicit_method, 1, a_drift, b_drift, c_drift, 2, a_diff, b_diff, c_diff)
end


"""
Tableau for the explicit 2-stage R2 method due to K. Burrage and P. Burrage
  Method cited in Eq. (51) in
   K. Burrage, P. Burrage (1996) "High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations".
  According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion
"""
function TableauBurrageR2()

    a = [[0.0 0.0]
         [2. / 3. 0.0]]
    b = [0.25, 0.75]
    c = [0., 2. / 3.]

    TableauSERK(:BurrageR2_explicit_method, 2, a, b, c, 2, a, b, c)
end


"""
Tableau for the explicit 4-stage CL method due to K. Burrage and P. Burrage
  Method cited in Eq. (56) in
   K. Burrage, P. Burrage (1996) "High strong order explicit Runge-Kutta methods for stochastic ordinary differential equations".
  According to the paper, the method has strong order 1.5 for one-dimensional Brownian motion. Reduces to the classical R-K method
  of order 4 when noise is zero.
"""
function TableauBurrageCL()

    a_drift = [[0.0 0.0 0.0 0.0]
               [0.5 0.0 0.0 0.0]
               [0.0 0.5 0.0 0.0]
               [0.0 0.0 1.0 0.0]]
    b_drift =  [1. / 6., 1. / 3., 1. / 3., 1. / 6.]
    c_drift =  [0.0, 0.5, 0.5, 1.0]

    a_diff =  [[ 0.0           0.0         0.0         0.0]
               [-0.7242916356  0.0         0.0         0.0]
               [ 0.4237353406 -0.199443705 0.0         0.0]
               [-1.578475506   0.840100343 1.738375163 0.0]]
    b_diff =   [-0.7800788474, 0.0736376824, 1.486520013, 0.2199211524]
    c_diff =   [0.0, -0.7242916356, 0.2242916356, 1.0]

    a_diff2=  [[ 0.0           0.0         0.0         0.0]
               [ 2.70200041    0.0         0.0         0.0]
               [ 1.757261649   0.0         0.0         0.0]
               [-2.918524118   0.0         0.0         0.0]]
    b_diff2=   [ 1.693950844, 1.636107882, -3.024009558, -0.3060491602]
    c_diff2=   [ 0.0, 2.70200041, 1.757261649, -2.918524118]

    TableauSERK(:BurrageCL_explicit_method, 4, a_drift, b_drift, c_drift, 4, a_diff, b_diff, c_diff, 0, a_diff2, b_diff2, c_diff2)
end


"""
Tableau for the explicit 4-stage E1 method due to K. Burrage and P. Burrage
  Method cited in Eq. (4.2)-(4.3) in
   K. Burrage, P. Burrage (2000) "Order conditions for stochastic Runge-Kutta methods by B-series".
  According to the paper, the method has strong order 1.0 for one-dimensional Brownian motion.
"""
function TableauBurrageE1()

    a_drift = [[0.0      0.0     0.0 0.0]
               [2. / 3.  0.0     0.0 0.0]
               [3. / 2. -1. / 3. 0.0 0.0]
               [7. / 6.  0.0     0.0 0.0]]
    b_drift =  [0.25, 0.75, -0.75, 0.75]
    c_drift =  [0.0, 2. / 3., 7. / 6., 7. / 6.]

    a_diff =  [[ 0.0      0.0      0.0    0.0]
               [ 2. / 3.  0.0      0.0    0.0]
               [ 0.5      1. / 6.  0.0    0.0]
               [-0.5      0.0      0.5    0.0]]
    b_diff =   [-0.5, 1.5, -0.75, 0.75]
    c_diff =   [0.0, 2. / 3., 2. / 3., 0.0]

    a_diff2=  [[ 0.0    0.0    0.0    0.0]
               [ 0.0    0.0    0.0    0.0]
               [ -2. / 3. 0.0    0.0    0.0]
               [ 1. / 6.  0.5    0.0    0.0]]
    b_diff2=   [ 1.5, -1.5, 0.0, 0.0]
    c_diff2=   [ 0.0, 0.0, -2. / 3., 2. / 3.]

    TableauSERK(:BurrageE1_explicit_method, 4, a_drift, b_drift, c_drift, 4, a_diff, b_diff, c_diff, 0, a_diff2, b_diff2, c_diff2)
end


"""
Tableau for the explicit 5-stage G5 method due to K. Burrage and P. Burrage
  Method cited in Section 4 of
   K. Burrage, P. Burrage (2000) "Order conditions for stochastic Runge-Kutta methods by B-series".
  According to the paper, the method has strong order 1.5 for one-dimensional Brownian motion.
"""
function TableauBurrageG5()

    a_drift = [[ 0.0               0.0               0.0               0.0              0.0]
               [ 0.52494822322232  0.0               0.0               0.0              0.0]
               [ 0.07167584568902  0.27192330512685  0.0               0.0              0.0]
               [ 0.13408162649312  0.24489042208103 -0.02150276857782  0.0              0.0]
               [-0.07483338680171 -0.07276896351874  0.55202897082453 -0.50752343840006 0.0]]

    b_drift =  [-5.60958180689351, -0.67641638321828, -5.44025143434789, 8.76396506407891, 3.96228456038077]
    c_drift =  [ 0.0, 0.52494822322232, 0.34359915081587, 0.35746927999633, -0.10309681789598002]

    a_diff =  [[ 0.0               0.0               0.0              0.0              0.0]
               [ 0.52494822322232  0.0               0.0              0.0              0.0]
               [ 0.49977623528582 -0.14576793502675  0.0              0.0              0.0]
               [ 0.60871134749146  0.58291821365556 -0.94596532788804 0.0              0.0]
               [-0.04005606091567 -0.22719654397712 -0.1292628422212  0.42881625288868 0.0]]

    b_diff =   [ 6.68050246229861, 0.0, 4.28273528343281, -3.25408735237225, -6.7091503933593]
    c_diff =   [ 0.0, 0.52494822322232, 0.35400830025907, 0.24566423325898, 0.03230080577469002]

    a_diff2=  [[ 0.0               0.0              0.0               0.0              0.0]
               [ 0.0               0.0              0.0               0.0              0.0]
               [-0.23101439602069  0.59278042710702 0.0               0.0              0.0]
               [-0.54946055077234  0.86811263829203 0.06772607159055  0.0              0.0]
               [ 0.03847082280344 -0.16953882944054 0.88387761274601 -0.85833118389518 0.0]]

    b_diff2=   [ 1.90494977554482, -1.90494977554482, 0.0, 0.0, 0.0]
    c_diff2=   [ 0.0, 0.0, 0.36176603108633, 0.38637815911024, -0.10552157778627003]

    TableauSERK(:BurrageG5_explicit_method, 4, a_drift, b_drift, c_drift, 4, a_diff, b_diff, c_diff, 0, a_diff2, b_diff2, c_diff2)
end
