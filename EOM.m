function dY = EOM(T, Y, m, g, lambda1, lambda2, lambda3, type)

dY = zeros(20,1);
% define state variables
xdot= Y(1);
ydot= Y(2);
zdot= Y(3);
x=Y(4);
y=Y(5);
z=Y(6);
omega1 = Y(7);
omega2 = Y(8);
omega3 = Y(9);
e0= Y(10);
e1= Y(11);
e2= Y(12);
e3= Y(13);
omega1dot = Y(14);
omega2dot = Y(15);
omega3dot = Y(16);
e0dot=Y(17);
e1dot=Y(18);
e2dot=Y(19);
e3dot=Y(20);

L = [-e1, e0, e3, -e2;
  -e2, -e3, e0, e1;
  -e3, e2, -e1, e0];

p = [ e0;
      e1;
      e2; 
      e3];
  
omega =[omega1;
        omega2;
        omega3];
    
omegadot =[omega1dot;
    omega2dot;
    omega3dot];

dY=[0 ;
    0;
    -g;
    xdot;
    ydot;
    zdot;
    (lambda2-lambda3)/lambda1*omega2*omega3;
    (lambda3-lambda1)/lambda2*omega1*omega3;
    (lambda1-lambda2)/lambda3*omega2*omega1;
    1/2*L'*omega;
    (lambda2-lambda3)/lambda1*(omega2dot*omega3+omega2*omega3dot);
    (lambda3-lambda1)/lambda2*(omega1dot*omega3+omega1*omega3dot);
    (lambda1-lambda2)/lambda3*(omega2dot*omega1+omega2*omega1dot);
    1/2*L'*omegadot-1/2*(omegadot'*omegadot)*p];

end


