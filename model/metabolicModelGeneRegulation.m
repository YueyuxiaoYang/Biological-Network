function dx = metabolicModelGeneRegulation(t,x,p)
   % Defines ODE model describing simple metabolic network used in
   % CompSysBio course.
   % The model includes regulation on genetic level: inhibition of synthesis of
   % e2 by x1 or v1.

   % HdJ 16/3/17

    
   % Define variables
   
   s1 = x(1);
   s2 = x(2);
   x1 = x(3);
   x2 = x(4);
   m = x(5);
   b = x(6);
   e1 = x(7);
   e3 = x(8);
   
   
   % Define parameters
   
   K1 = p(1);
   K2 = p(2);
   k1 = p(3);
   k2 = p(4);
   k3 = p(5);
   k4 = p(6);
   k5 = p(7);
   beta = p(8); 
   c1 = p(9);
   c3 = p(10);
   g1 = p(11);
   g3 = p(12);
   L2 = p(13);

   
   % Define reaction rates
   
   v1 = k1*e1*s1/(s1+K1);
   v2 = k2*x1;
   v3 = k3*e3*s2/(s2+K2);
   v4 = k4*x2;
   v5 = k5*m;
   r1 = c1;
   r2 = g1*e1;
   r5 = c3*L2^2/(L2^2 + x1^2);
%   r5 = c3*L2^2/(L2^2 + v1^2);
   r6 = g3*e3;
   
   
   % Define dynamics
   
   ds1 = -v1*b;
   ds2 = -v3*b;
   dx1 = v1 - v2;
   dx2 = v3 - v4; 
   dm = 4*v2 + v4 - 10*v5;
   db = v5 * beta * b;
   de1 = r1 - r2;
   de3 = r5 - r6;
   
   dx = [ds1; ds2; dx1; dx2; dm; db; de1; de3];
   
end
