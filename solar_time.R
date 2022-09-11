### Parameters of Earth Orbit
ecc <- .0167                                                                  ## eccentricity
lambda_1 <- (185 * 86400 + 22 * 3600 + 27 * 60) / (365.2422 * 86400) * 2 * pi ## eclipse longitude when aphelion occurs in 2021
lambda_2 <- (78 * 86400 + 9 * 3600 + 37 * 60) / (365.2422 * 86400) * 2 * pi   ## spring equinox in 2021
T_rev <- 365.2422 * 86400                                                     ## solar year
Omega_rot <- 2 * pi / 86400 * 366.2422 / 365.2422                             ## Earth rotation angular speed 
delta_max <- 23.4 / 180 * pi                                                  ## Eclipse Plane Declination

### Sampling in the lambda and wave number domain
### Lambda is the eclipse longitude
res <- 365 * 24
lambda <- seq(0, 2 * pi, length.out = res + 1)
lambda <- lambda[-length(lambda)]
k_lambda <- 0:(res / 2 - 1)
k_lambda <- c(k_lambda, -rev(1:(res / 2)))

### Time - Lambda Relationship
diff_t <- 1 / (1 - ecc * cos(lambda - lambda_1))^2
DIFF_t <- fft(diff_t)
K_rev <- T_rev / 2 / pi / (Re(DIFF_t[1]) / res)
Re(DIFF_t[1]) / res

INT_t <- DIFF_t / complex(imaginary = k_lambda)
INT_t[1] <- 0
t <- K_rev * fft(INT_t, inverse = T) / res
t <- Re(t) + K_rev * Re(DIFF_t[1]) / res * lambda

### Alpha - Lambda Relationship
### Alpha is the equatorial longitude
diff_alpha <- cos(delta_max) * (1 + tan(lambda - lambda_2)^2) / (1 + cos(delta_max)^2 * tan(lambda - lambda_2)^2)
DIFF_alpha <- fft(diff_alpha)
INT_alpha <- DIFF_alpha / complex(imaginary = k_lambda)
INT_alpha[1] <- 0
alpha <- fft(INT_alpha, inverse = T) / res
alpha <- Re(alpha) + Re(DIFF_alpha[1]) / res * lambda

### h - Lambda Relationship
### h is the local solar time
diff_h <- Omega_rot * K_rev * diff_t - diff_alpha
DIFF_h <- fft(diff_h)
INT_h <- DIFF_h / complex(imaginary = k_lambda)
INT_h[1] <- 0
h <- fft(INT_h, inverse = T) / res
h <- Re(h) + Re(DIFF_h[1]) / res * lambda

## Mean Solar Time and Error
h_bar <- t / 86400 * 2 * pi
h_error <- h - h_bar

plot(t / T_rev * 365.2422, h_error / 2 / pi * 1440, type = "l", xlab = "Day of Year", ylab = "Equation of Time (Minutes)")
lines(c(0, max(t / T_rev * 365.2422)), c(0, 0), lty = 2)