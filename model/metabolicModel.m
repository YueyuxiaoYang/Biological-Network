function dx = metabolicModel(t,x,p)
   % Defines ODE model describing simple metabolic network used in
   % CompSysBio course.
   % The model does not include regulation
   
   % HdJ 16/3/17


   % Define variables
   
   s1 = x(1);
   s2 = x(2);
   x1 = x(3);
   x2 = x(4);
   m = x(5);
   b = x(6);
   
   
   % Define parameters
   
   K1 = p(1);
   K2 = p(2);
   k1 = p(3);
   k2 = p(4);
   k3 = p(5);
   k4 = p(6);
   k5 = p(7);
   beta = p(8); 

   
   % Define reaction rates
   
   v1 = k1*s1/(s1+K1);
   v2 = k2*x1;
   v3 = k3*s2/(s2+K2);
   v4 = k4*x2;
   v5 = k5*m;
   
   
   % Define dynamics
   
   ds1 = -v1*b;
   ds2 = -v3*b;
   dx1 = v1 - v2;
   dx2 = v3 - v4; 
   dm = 4*v2 + v4 - 10*v5;
   db = v5 * beta * b;

   
   dx = [ds1; ds2; dx1; dx2; dm; db];
   
end
