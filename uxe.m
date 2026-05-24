function [uxv,uyv]=uxe(z) %u_x,u_y
x=z(1); y=z(2);
uxv=x*y*(y - 1) + y*(x - 1)*(y - 1);
uyv=x*y*(x - 1) + x*(x - 1)*(y - 1);