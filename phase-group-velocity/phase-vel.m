# Phase Velocity of Plane Monochromatic Wave along Z-Axis
# Author: Arvin Sharma
# E(z,t)=E_0*cos(kz-wt)
# - k = wave number 2pi/lambda [proportionality constant to convert between radians and span (meters, micrometers, etc)
# - z = instantaneous distance along span being analyzed
# - omega = angular frequency, frequency expressed in rad/sec
# - time = instantaneous time in seconds
# We want to determine the movement of a particular phase, therefore
# - kz-wt = K' where K' is a constant value for the phase to analyze
#   - k dz/dt - w = 0 provides the variation with respect to time of the constant phase, defined as the phase velocity
#       - This is how the point of constant phase in the cycle varies in time.
#           - We see that this point of constant phase varies in space, or span along the z-axis in time
#   - dz/dt = w/k = vp 
#   - w/k = 2pif/(2pi/lambda)=lambda*f=c/n
#   - vp=c/n
# In this demonstration, we take a fixed point along the span of transmission line being analyzed
#   - this is arbitrarily defined at 0.75 micro-meters
# We then use the wave-function varied in time and plot the result, showing the variation in phase by plotting amplitude versus time.
#
# References: 
#   - M. S. Wartak, Computational Photonics: An Introduction with MATLAB. Cambridge: Cambridge University Press, 2013. 

E_0=10; # Amplitude
n=3.4; # Refractive Index
c=3e14; # Velocity of Light in micro-meters/s
v_p=c/n; # Phase Velocity
lambda=1.0; # Micro-meters
k=2*pi/lambda; # Wave Number rad/micro-meters
frequency=v_p/lambda; # Cycles/sec
z=0.75; # Fixed Distance to Analyze in Micro-meters
omega=2*pi*frequency;
t=linspace(0,2*1/frequency,101); # Over 2 Periods
A = E_0*cos(k*z - omega*t);
plot(t/(1e-14),A);
title("Phase Velocity");
xlabel("Time (10^-14)");
ylabel("Amplitude");
axis("tight");
waitfor(gcf);

