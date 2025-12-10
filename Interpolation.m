function Phi = Interpolation(x,y,A,B,C,ForG)

if abs(B-C) < 0.0000001
    fprintf('we have a singularity\n')
    D = (A-B)/(10^(-16));
else
    D = (A-B)/(B-C);
end

x1 = 0.0;
y1 = 1/sqrt(3);
x2 = 0.0;
y2 = 0.0;
G1 = 4.0*pi/(sqrt(3.0));
m1 = cos(sqrt(3.0)*G1*x1/2.0);
m2 = cos(sqrt(3.0)*G1*x2/2.0);
a1 = cos(G1*y1) - cos(G1*y2);
a3 = 2.0*cos(G1*y1/2.0) * m1 - 2.0*cos(G1*y2/2.0)* m2;
alpha = a1 + a3;

a2 = sin(G1*y1) - sin(G1*y2);
a4 = 2.0*sin(G1*y1/2.0) * m1 - 2.0*sin(G1*y2/2.0)* m2;
beta = a2 - a4;

x1 = 0.0;
y1 = 0.0;
x2 = 0.0;
y2 = 2.0/sqrt(3);
G1 = 4.0*pi/(sqrt(3.0));
m1 = cos(sqrt(3.0)*G1*x1/2.0);
m2 = cos(sqrt(3.0)*G1*x2/2.0);
a1 = cos(G1*y1) - cos(G1*y2);
a3 = 2.0*cos(G1*y1/2.0) * m1 - 2.0*cos(G1*y2/2.0)* m2;
gamma = a1 + a3;

a2 = sin(G1*y1) - sin(G1*y2);
a4 = 2.0*sin(G1*y1/2.0) * m1 - 2.0*sin(G1*y2/2.0)* m2;
delta = a2 - a4;

phi = atan((1.0/(delta/beta*D-1.0)*((delta*alpha-beta*gamma)/(beta*delta)))-(gamma)/(delta));

if abs(B-C) < 0.0000001
    c1 = (10^(-16))/(2.0*(gamma*cos(phi)+delta*sin(phi)));
else
    c1= (B-C)/(2.0*(gamma*cos(phi)+delta*sin(phi)));
end

Ax = 0.0;
Ay = 1/sqrt(3);
G1 = 4.0*pi/(sqrt(3.0));
c0 = A - 2.0 * c1 * cos(phi - G1 * Ay) - 4.0 * c1 * cos(G1 * Ay / 2.0 + phi) * cos(sqrt(3.0) * G1 * Ax / 2.0);

acc = 1.0;
G1 = 4.0*pi/(sqrt(3.0)*acc);
if ForG == 1
    f1 = 2.0*c1*cos(phi-G1*y) + 4.0*c1*cos(G1*y/2.0 + phi)*cos(sqrt(3.0)*G1*x/2.0);
    Phi = c0 + f1;
elseif ForG == 2
    Phi = c0*cos(sqrt(3)*G1*x/2)*cos(G1*y/2-phi)-2*c1*cos(G1*y+phi)-1j*2*sqrt(3)*c1*sin(sqrt(3)*G1*x/2)*sin(G1*y/2-phi);
end
