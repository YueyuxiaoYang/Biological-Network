function dx = metabolicModelGeneRegulationGrowthControlNewBiomass(t,x,p)
   % Defines ODE model describing simple metabolic network used in
   % CompSysBio course.
   % The model includes regulation on genetic level: inhibition of synthesis of
   % e2 by x1 or v1.
   % Moreover, it includes the synthesis of enzymes from precursors and defines growth rate as the relative accumulation of the total amount of cellular constituents.
   
   % HdJ 16/3/17

   
   % Define variables
   
   s1 = x(1);
   s2 = x(2);
   x1 = x(3);
   x2 = x(4);
   m = x(5);
   b = x(6);
   e1 = x(7);
   e2 = x(8);
   e3 = x(9);
   e4 = x(10);

   
   
   % Define parameters
   
   K1 = p(1);
   K2 = p(2);
   k1 = p(3);
   k2 = p(4);
   k3 = p(5);
   k4 = p(6);
   c1 = p(7);
   c2 = p(8);
   c3 = p(9);
   c4 = p(10);
   L2 = p(11);

   
   % Define reaction rates
   
   v1 = k1*e1*s1/(s1+K1);
   v2 = k2*e2*x1;
   v3 = k3*e3*s2/(s2+K2);
   v4 = k4*e4*x2;
   r1 = c1*m;
   r2 = c2*m;
   r3 = c3*m*L2^2/(L2^2 + x1^2);
   r4 = c4*m;
   mu = (4*v1 + v3)/(4*x1+x2+m+2.5*e1+2.5*e2+2.5*e3+2.5*e4);
   
   % Define dynamics
   
   ds1 = -v1*b;
   ds2 = -v3*b;
   dx1 = v1 - v2;
   dx2 = v3 - v4; 
   dm = 4*v2 + v4 - 2.5*r1 - 2.5*r2 - 2.5*r3 - 2.5*r4;
   db = mu*b;
   de1 = r1 - mu*e1;
   de2 = r2 - mu*e2;
   de3 = r3 - mu*e3;
   de4 = r4 - mu*e4;
   
   dx = [ds1; ds2; dx1; dx2; dm; db; de1; de2; de3; de4];
   
end
