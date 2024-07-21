# Group Velocity of Superposition of Two Plane Waves
# Author: Arvin Sharma
# 
# A group is formed through the superposition of waves, E_1 and E_2:
# - E_1(z,t)=E_0*cos(k_1z - w_1t)
# - E_2(z,t)=E_0*cos(k_2z - w_2t)
# The difference between these two waves is defined as 2dW and 2dK
# where 2dW is the difference in angular frequency and 2dK is the difference in wave number.
# - E_1=E_0*cos((k+dK)z - (w+dW)t)
# - E_2=E_0*cos((k-dK)z - (w-dW)t)
# Through Superposition we find:
# - E(z,t)=E_1(z,t)+E_2(z,t)
#   => E(z,t)=E_0(cos((k+dK)z-(w+dW)t) + cos((k-dK)z - (w-dW)t))
# In order to simplify this further, we consider two points along the unit circle, B=(cosB,sinB) and A=(cosa,sina) where a > B
# The difference between these two points is:
# - D = sqrt((x_1-x_2)^2 + (y_1 - y_2)^2)
#   => AB = sqrt((cosa-cosB)^2 + (sina-sinB)^2)
#   => (AB)^2 = (cosa-cosB)^2 + (sina-sinB)^2
#   => (AB)^2 = cos^2a - 2cosacosB+cos^2B + sin^2a-2sinasinB+sin^2B
#   [Identity: sin^2t+cos^2t=1]
#   => (AB)^2 = 2 - 2cosacosB - 2sinasinB
#   => (AB)^2 = 2 - 2(cosacosB+sinasinB)
# If these two points are rotated by angle B, such that point B is now on the X-axis, we get two new points C:(cost,sint) and D:(1,0)
# Note that the distance from C to D is equal to the distance from A to B, since both points are equally rotated about the origin, where t=a-B and a > B as mentioned above.
# Therefore, CD = AB, and we find CD to simplify the expression found in calculating AB.
# - CD = sqrt((cost - 1)^2 + (sint - 0)^2)
#   => (CD)^2 = (cost - 1)^2 + (sint - 0)^2
#   => (CD)^2 = cos^2 - 2cost + 1 + sin^2 t
#   [Identity: sin^2t + cos^2t = 1]
#   => (CD)^2 = 2 - 2cos(t)
#   [Note: t = a - B]
#   => (CD)^2 = 2 - 2cos(a -B)
# Since CD=AB, (CD)^2=(AB)^2
#   - (CD)^2=(AB)^2
#   => (2 - 2cos(a - B)) = (2 - 2(cosacosB+sinasinB))
#   => -2cos(a - B) = -2(cosacosB+sinasinB)
#   => cos(a - B) = cosacosB+sinasinB
# We need one more identity to simplify this expression, specifically for cos(a+B)
# - cos(a+B) = cos(a - (-B))
#   => cosacos(-B) + sinasin(-B)
#   => cosacosB+sina(-sinB) // cos(-X)=cos(X), sin(-X)=-sin(X)
#   => cosacosB-sinasinB
# Therefore, cos(a+B)=cosacosB-sinasinB
# Combining these two:
# - cos(a-B) + cos(a+B)=(cosacosB+sinasinB) + (cosacosB-sinasinB)
#   => cos(a-B) + cos(a+B) = 2cosacosB
# Recall, our expression for the superposition of plane waves E_1 and E_2:
# - E(z,t)=E_0(cos((k+dK)z-(w+dW)t) + cos((k-dK)z - (w-dW)t))
#   => ((k+dK)z-(w+dW)t) = (kz+dKz-wt-dWt)=((kz-wt)+(dKz-dWt))
#   => ((k-dK)z-(w-dW)t) = (kz-dKz-wt+dWt)=((kz-wt)-(dKz-dWt))
#   Notice the two expressions on the write are of the form (a+B) and (a-B)
#   where a = (kz-wt) and B=(dKz-dWt)
#   Therefore we apply the trig identity: cos(a-B)+cos(a+B)=2cosacosB,
#   => E(z,t)=E_0(2cos(kz-wt)*cos(dKz-dWt)) which is the simplified expression representing the modulated wave. 
# Note that the phase velocity is w/k while the envelope, group velocity is dw/dk.
#
# References: 
#   - M. S. Wartak, Computational Photonics: An Introduction with MATLAB. Cambridge: Cambridge University Press, 2013. 
#
E_0=10;
c=3e14; # Velocity of light in micro-meter/sec
n=3.4; # Refractive Index
v_p=c/n; # Phase Velocity
lambda=1.0; # Wavelength in micro-meters
frequency=v_p/lambda; # v=flambda, v/lambda=f
z=0.75; # Arbitrary point to analyze
omega=2*pi*frequency; # Angular frequency
k=2*pi/lambda; # Wave number
dOmega=omega/15.0; # Differential of angular frequency
dK = k/15.0; # Differential of wave number 
#
omega_1=omega+dOmega;
omega_2=omega-dOmega;
k_1=k+dK;
k_2=k-dK;
t=linspace(0,300e-15,101);
A=E_0*(cos(k*z - omega*t).*cos(dK*z - dOmega*t));
plot(t/(1e-15),A);
xlabel("Time (10^-15)");
ylabel("Amplitude");
axis("tight");
waitfor(gcf);

