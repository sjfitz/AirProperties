# AirProperties
 Calculates a variety of air thermodynamic properties from the measured temperature, pressure and humidity.

Variables able to calculated:
 * rho - [kg m^-3] - Density
 * mu - [N s m^-2] - Dynamic viscosity
 * k - [W m^-1 K^-1] - Thermal conductivity
 * c_p - [J kg^-1 K^-1] - Specific heat capacity (constant pressure)
 * c_v - [J kg^-1 K^-1] - Specific heat capacity (constant volume)
 * gamma - [1] - Ratio of specific heats
 * c - [m s^-1] - Speed of sound: c = (gamma*R*T/M)^0.5
 * nu - [m^2 s^-1] - Kinematic viscosity: nu = mu/rho
 * alpha - [m^2 s^-1] - Thermal diffusivity: alpha = k/(rho*c_p)
 * Pr - [1] - Prandtl number: Pr = mu*c_p/k
 * M - [kg mol^-1] - Molar mass of humid air
 * R - [J kg^-1 K^-1] - Specific gas constant
 * h - [%] - Relative humidity (if dew point inputted)

Calculations are based around atmospheric temperatures and pressures, i.e. not for the use of high temperature combustion. See references for details and limitations.

References:
Picard, A, Davis, RS, Glaser, M, Fujii, K, 2008, 'Revised formula for the density of moist air (CIPM-2007)', Metrologia, vol. 45, no. 2, pp. 149-155. DOI: http://dx.doi.org/10.1088/0026-1394/45/2/004

Tsilingiris, P, 2008, 'Thermophysical and transport properties of humid air at temperature range between 0 and 100°C', Energy Conversion and Management, vol. 49, no. 5, pp.1098-1110. DOI: https://doi.org/10.1016/j.enconman.2007.09.015
