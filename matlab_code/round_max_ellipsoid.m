function [A, b, center, radius, T, T_shift, t_r] = round_max_ellipsoid(A, b, x0)

    t_r=1;
    dim = size(A,2);
    T = eye(dim);
    T_shift = zeros(dim,1);
    if dim<200

        [T_sh, Tmve, converged] = mve_run_cobra(A, b, x0, 1e-8);
        %Tmve
        %p_shift = p_shift + N_total*T_sh;
        %N_total = N_total * Tmve;
        T_shift = T_shift + T_sh;
        T = T * Tmve;
        b = b - A*T_sh;
        A = A * Tmve;
        %[N_total, p_shift, T] = shiftPolytope(N_total, p_shift, T, Tmve, T_sh);
        if converged~=1
            fprintf('There was a problem with finding the maximum volume ellipsoid.\n');
        end
        %x0 = zeros(dim, 1);
        [x0, ~]=get_cheb(A, b);
    else
        %if dimension is large, be a little more careful
        %the below loop looks silly for an interior point method, but is actually
        %quite important for numerical stability. while normally you'd only call an optimization routine
        %once, we call it iteratively--where in each iteration, we map a large
        %ellipsoid to the unit ball. the idea is that the "iterative roundings"
        %will make each subsequent problem easier to solve numerically.
        max_its = 20;
        its = 0;
        reg=1e-3;
        Tmve = eye(dim);
        converged = 0;
        while (max(eig(Tmve))>6*min(eig(Tmve)) && converged~=1) || reg>1e-6 || converged==2
            tic;
            its = its+1
            %check if we can use the good lp solver
            %if exist('solveCobraLP')==2
            %    [x0,dist] = getCCcenter(P.A,P.b);
            %else
                %let mve_run use matlab's lp solver to select a starting point
            %    [~,x0] = mve_presolve_cobra(P.A,P.b,150,1e-6);
            %end

            reg = max(reg/10,1e-10);
            [T_sh, Tmve, converged] = mve_run_cobra(A, b, x0, reg);
            %p_shift = p_shift + N_total * T_sh;
            %N_total = N_total * Tmve;
            T_shift = T_shift + T_sh;
            T = T * Tmve;
            b = b - A * T_sh;
            A = A * Tmve;
            t_r = t_r * det(Tmve');
            %[P, N_total, p_shift, T] = shiftPolytope(P, N_total, p_shift, T, Tmve, T_shift);
            row_norms = sqrt(sum(A.^2,2));
            A = diag(1./row_norms)*A;
            b = diag(1./row_norms)*b;
            if its==max_its
                break;
            end
            rr=toc;
            rr
            [x0, ~]=get_cheb(A, b);

            %fprintf('Iteration %d: reg=%.1e, ellipsoid vol=%.1e, longest axis=%.1e, shortest axis=%.1e, x0 dist to bdry=%.1e, time=%.1e seconds\n', its, reg, det(Tmve), max(eig(Tmve)), min(eig(Tmve)), dist, toc);
        end
        

        if its==max_its
            fprintf('Reached the maximum number of iterations, rounding may not be ideal.\n');
        end
    end
    [center, radius]=get_cheb(A, b);
end