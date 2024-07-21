# Rectangular Pulses
# Author: Arvin Sharma
#
# Theory:
# A Simple rectangular pulse is defined in terms of the time interval T:
# s_rect(t) = s_0 if 0<t<T else 0      
# Next we consider the synthesis of Square waves through the addition of odd harmonics.
# Sum_j=1^infinity (s_0/n)*sin(2*pi*n*f*t) where n=2*j-1
#
# We then consider more general Rectangular pulses formed using a Fourier series of approximate Sinc functions.
# Sum_j=1^infinity (s_0*(t_max-t_min)/(2*n*pi))*sin((n*pi*(t-t_min) / (t_max-t_min)) where n=2*j-1
# - We break this down as follows:
# - (n*pi/T) is equivalent to n*pi*f where n is the order of the harmonic, 1/T is the frequency
# - (t-t_min) translates the pulse to t_min for the purposes of comparison to the pure rectangular pulse.
# - s_0 is the amplitude
# - (t_max-t_min) is the period
# - 1/(t_max-t_min) is the frequency
# - (t_max-t_min)/(2*n*pi) forms the approximation of the sinc function with sin((n*pi/T)*(t-t_min)) where sincx=sinx/x
#
# Analytically we evalute the Fourier Transform for the rectangular pulse as follows:
# - s_rect,sym(t) = 1 if |t|<tau else 0 if |t|>tau
# - Fourier transform:
#   - S_rect,sym(w)=integral_-inf^+inf (x(t)e^(-iwt)dt)=int_-tau^+tau (e^(-iwt)dt) = 1/(-iw) * e^(-iwt) |_-tau^+tau
#       => 1/(-iw) * (e^(-iwt)-e^(iwt))
#       => 1/(-iw) * ((cos(-wt) + isin(-wt)) - (coswt + isinwt))
#       => 1/(-iw) * ((cos(wt) - isin(wt)) - coswt + isinwt)
#       => 1/(-iw) * (-2isinwt)
#       => 2sinwt/w
#
s_0=10.0; # Pulse Amplitude
t_min=2; # Pulse Start (sec)
t_max=10.0; # Pulse End (sec)
T=t_max-t_min; # Period (sec)
frequency=1/T; # Cycles/sec
N=101; # Number of points
dt=1/N;
t=0:dt:1.5*t_max; # Varied Time plotted over the boundary tmax to provide extra space.
x=s_0*0.5*(sign(t-t_min)-sign(t-t_max)); # We determine if t > t_min and t < t_max by taking the sign of the difference of t with t_min and t_max
                                     # which determines whether t is greator or less than the quantity. Since t must be > t_min and < t_max, and 
                                     # t_min < t_max then t is either < t_min and < t_max, > t_min and < t_max, or > t_min and > t_max.
                                     # We are interested in the second case, using the sign operator we get (-1) and (-1) in the first case,
                                     # (-1) and (1) in the second case, and (1) and (1) in the third case. Subtracting these quantites eliminates
                                     # the first and third case, leaving only the second case. Multiplying by one half provides us with unity.
#
h=plot(t,x); # Plot the Pure Rectangular Pulse
# Plot the Rectangular pulse formed by Odd Harmonics.
#
hold on
#
s=0.0;
for n = [1,3,5,7,9,11,13]
    s = s + s_0*(t_max-t_min) / (n*pi*2) * (sin(n*pi*(t-t_min)./(t_max-t_min)));
    plot(t,s);
endfor
waitfor(gcf);

