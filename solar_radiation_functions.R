#### Solar Radiation Model with Atmospheric Mass
### Function for determining daytime or not at a site
daytime_get <- function(lon, lat, lon_0){
  time_del <- lon - lon_0
  
  # Normal Vector of the Sun
  n_Sun = matrix(0, nrow = res, ncol = 3)
  # n_Sun[, 1] = sin(delta_Sun) * cos(lat) - cos(delta_Sun) * sin(lat) * cos(h + time_del)
  # n_Sun[, 2] = cos(delta_Sun) * sin(h + time_del)
  n_Sun[, 3] = sin(delta_Sun) * sin(lat) + cos(delta_Sun) * cos(lat) * cos(h + time_del)
  
  day = n_Sun[, 3] >= 0
  return(day)
}

### Function for finding maximum full load hours of a specific site and PV angle
FLH_PV <- function(phi_PV, azi_PV, zen_PV){
  # Normal vector of the panel
  n_PV = c(sin(zen_PV) * cos(azi_PV), sin(zen_PV) * sin(azi_PV), cos(zen_PV))
  
  # Normal Vector of the Sun
  n_Sun = matrix(0, nrow = res, ncol = 3)
  
  n_Sun[, 1] = sin(delta_Sun) * cos(phi_PV) - cos(delta_Sun) * sin(phi_PV) * cos(h)
  n_Sun[, 2] = cos(delta_Sun) * sin(h)
  n_Sun[, 3] = sin(delta_Sun) * sin(phi_PV) + cos(delta_Sun) * cos(phi_PV) * cos(h)
  
  # Atmospheric Mass and Diffusion Coefficient
  AM <- sqrt(ratio_atmo^2 * n_Sun[, 3]^2 + 2 * ratio_atmo + 1) - ratio_atmo * n_Sun[, 3]
  coeff_diff <- 1 - 0.75^AM
  
  # Coefficient of Direct Radiation Loss 
  coeff_Sun_PV = n_Sun %*% n_PV * (n_Sun[, 3] >= 0)
  coeff_Sun_PV = coeff_Sun_PV * (coeff_Sun_PV >= 0)
  coeff_Sun_PV = coeff_Sun_PV * (1 - coeff_diff)
  FLH_PV = (coeff_Sun_PV[1:(length(t) - 1)] + coeff_Sun_PV[2:length(t)]) / 2 * diff(t)
  FLH_PV = sum(FLH_PV) / 3600
  
  return(FLH_PV)
}

## Loacation and Site Dependent Variables
phi_PV = 25 * pi / 180                                           ## Latitude of location
azi_PV = 180 * pi / 180                                          ## Azimuth angle of the panel
zen_PV = 5 * pi / 180                                            ## Zenith angle of the panel

zen_angle = seq(0, 90, .2) * pi / 180
FLH_diff_zen_65 = c()
FLH_diff_zen_45 = c()
FLH_diff_zen_25 = c()

for(angle in 1:length(zen_angle)){
  FLH_diff_zen_65[angle] = FLH_PV(65 * pi / 180, azi_PV, zen_angle[angle])
  FLH_diff_zen_45[angle] = FLH_PV(45 * pi / 180, azi_PV, zen_angle[angle])
  FLH_diff_zen_25[angle] = FLH_PV(25 * pi / 180, azi_PV, zen_angle[angle])
}

plot(zen_angle * 180 / pi, FLH_diff_zen_65, type = "l", ylim = range(c(FLH_diff_zen_25, FLH_diff_zen_45, FLH_diff_zen_65)), ylab = "Maximum Full Load Hours of Direct Solar Radiation (hr)", xlab = "Zenith Angle of PV (Degrees)")
lines(zen_angle * 180 / pi, FLH_diff_zen_45, col = "blue")
lines(zen_angle * 180 / pi, FLH_diff_zen_25, col = "red")