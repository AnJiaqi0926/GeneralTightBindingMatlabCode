function t=Hopping(R)

Rz = R(3);
Rd = sqrt(R(1)^2+R(2)^2+R(3)^2);

Vpppi0 = -2.7;
Vppsigma0 = 0.48;

a = 2.438977765130281;
a0 = a/sqrt(3);
d0 = 3.35;
r0 = 0.184*a;

Vpppi = Vpppi0*exp(-(Rd-a0)/r0);
Vppsigma = Vppsigma0*exp(-(Rd-d0)/r0);

t = Vpppi*(1-(Rz/Rd)^2)+Vppsigma*((Rz/Rd)^2);
