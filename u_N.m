function unv=u_N(P1,P2) % normal derivative of u

N=-[-(P2(2)-P1(2)) P2(1)-P1(1)]/norm(P2-P1);

MP=(P1+P2)/2;
[ux uy]=uxe(MP);
unv=[ux uy]*N';