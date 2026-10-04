% Runs the convergence study (ocfem.m), prints a results table and saves a log-log plot.
% Works in MATLAB and in GNU Octave.
ocfem;                       % defines h, L2e, H1e, ocl2, och1 (h is proportional to the mesh size)
hh = 1 ./ 2.^(1:numel(L2e)); % mesh size (leg length) at levels 1..6

fprintf('Level  h        L2 error     H1 error     L2 rate  H1 rate\n');
for k = 1:numel(L2e)
    if k == 1
        fprintf('%d      1/%-5d  %.3e    %.3e    -        -\n', k, round(1/hh(k)), L2e(k), H1e(k));
    else
        fprintf('%d      1/%-5d  %.3e    %.3e    %.2f     %.2f\n', k, round(1/hh(k)), L2e(k), H1e(k), ocl2(k-1), och1(k-1));
    end
end

figure; loglog(hh, L2e, 'o-', hh, H1e, 's-', hh, L2e(1)*(hh/hh(1)).^2, 'k--', hh, H1e(1)*(hh/hh(1)), 'k:');
xlabel('mesh size h'); ylabel('error'); grid on;
legend('L2 error', 'H1 seminorm error', 'slope 2', 'slope 1', 'Location', 'southeast');
title('P1 FEM convergence');
print('-dpng', 'convergence.png');
