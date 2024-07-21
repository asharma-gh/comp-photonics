# The Super-Gaussian Pulses
# Author: Arvin Sharma
#  
# Theory:
#
# A Gaussian pulse is defined by the Gaussian profile:
# - s_G(t)=A/sigma*sqrt(2pi) * e^(-0.5)*(t/sigma)^2
#   - the pulse has area A=P_0*tau_0 where tau_0=sigma*(sqrt(2pi)) and P_0 is maximum at t=0
#   - The spectral function S_G(w) obtained by the Fourier transform is also a Gaussian function
#       - S_G(w)=Ae^((-1/2)*(w*sigma)^2)
#
# The derivation of the Gaussian function is explained as follows:
#
# The probability density function p(x)=1/(sigma*sqrt(2*pi))*e^((-1/2)*(x-mu/sigma)^2 is used to compute probability values for a given mu and sigma for a continuous random variable x.
# This function p(x) represents the probability of a particular value having the value x.
# We have two distributions in this example, the X and Y distributions for throws along the X and Y axis. 
# The probability of a particular value in the span from x to x+dX is defined as p(x)dX. 
# This choice of differentials allows for the probability distribution to be constructed by integrating over the infinitessimal areas at each point in the domain. 
# We similarly define the distribution p(y)dY as the probability of a particular value in the span between y and y+dY. 
# We assume the values of dX and dY to be infinitessimal such that the value of P(x) and P(y) is constant over their span for all points x and y.
# Now, we define the probability of a set of paricular random variables having the value (x,y) as p(x)dX * p(y)dY which we define as a new function g(r)dXdY where g(r)=p(x)p(y)
# The independent variable (r) is in polar coordinates where (r,theta) indicate the point (x,y)
# We differentiate both sides of g(r)=p(x)p(y) with respect to theta,
# 0 = p(x) * dp(y) / dtheta + p(y) * dp(x) / dtheta, where dg(r)/dtheta = 0 since g is not dependent on theta.
# Then using x=rcos(theta) and y=rsin(theta), we expand dp(y)/dtheta to dp(y)/dtheta*dy/dtheta and dp(x)/dtheta to dp(x)/dtheta * dx/dtheta:
# 0 = p(x)*dp(y)/dtheta * (rcos(theta)) + p(y) * dp(x)/dtheta * (-rsin(theta))
# where dx/dtheta = -rsin(theta) and dy/dtheta = rcos(theta)
# Rewriting this as 0 = p(x)p'(y)x - p(y)p'(x)y, we solve this differential equation using the separation of variables
# p'(x)/xp(x) = p'(y)/yp(y)
# We set both sides of this equation to a constant C, since both x and y are independent variables, where
# p'(x)/xp(x) = C, such that p'(x)/p(x) = Cx, Integrating both sides
# Int(dp/dx / p * dx) = Int(Cxdx) = Int(dp(x) / p(x) = Int(Cxdx) = Int(1/p(x) dpx) = Int(Cxdx)
# => ln(p(x)) = Cx^2/2 + C
# => p(x) = e^(Cx^2/2 + C) = Ke^(Cx^2/2).
# The assumption is then made for the value C to be negative based on larger values of X being less likely to occur, forming the characteristic shape of the Gaussian distribution.
# => p(x) = Ke^(kx^2/2) where K=1/(sigma*sqrt(2*pi)) and k=(-1/sigma^2) forms the equation for the Gaussian.
# The Super Gaussian generalizes the Gaussian replacing the exponent square on the variable x/sigma with 2m, forming the equation:
# S_SG= A/(sigma*sqrt(2*pi))*e^((-1/2)*(t/sigma)^2m)
# Where m controls the degree of edge sharpness:
#       - m=1 creates the ordinary Gaussian Pulse
#       - as m becomes large the pulse more closely approximates the rectangular pulse.
# 
# References: 
#   - M. S. Wartak, Computational Photonics: An Introduction with MATLAB. Cambridge: Cambridge University Press, 2013. 
#   - 
A=1;
sigma=30;
N=303;
t=linspace(-50,50,N);
hold on
for m=[1,2,3,4]
    p=A/(sigma*sqrt(2*pi))*exp(-t.^((2*m))/(sigma^(2*m)));
    h=plot(p);
endfor
xlabel("Time");
ylabel("Amplitude");
axis("tight")
waitfor(gcf);

