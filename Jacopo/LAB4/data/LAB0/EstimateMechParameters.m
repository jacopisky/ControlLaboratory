
[Beq_a, tausf_a] = estimateBeq_tausf("data/least_square_positive.mat");
[Beq_b, tausf_b] = estimateBeq_tausf("data/least_square_negative.mat");

est_par.B_eq = (Beq_a + Beq_b)/2;
est_par.tau_sf = (tausf_a + tausf_b)/2;
est_par.J_eq = estimateJeq("data/accel_decel.mat", est_par.B_eq, est_par.tau_sf);


% clearing useless variables
clear tausf_a;
clear tausf_b;
clear Beq_a;
clear Beq_b;

function Jeq = estimateJeq(file, Beq, tausf)
    run("LoadParameters.m");

    load(file);
    t = ia_best.time;
    n_sample = length(t);
    i = ia_best.signals.values;
    vs = sens.curr.Rs * i;
    wm = deg2rad(al_best_meas1.signals.values * gbox.N);
    am = deg2rad(al_best_meas.signals.values * gbox.N);
    torque = torquem_best.signals.values - (Beq * wm + tausf * sign(wm) / gbox.N);
    am = abs(am);
    torque = abs(torque);
    sum = 0;
    n = 0;
    for i = 1:n_sample
        if am(i) ~= 0
            sum = sum + torque(i)/am(i);
            n = n + 1;
        end
    end
    Jeq = sum / n;
end

function [Beq, tausf] = estimateBeq_tausf(file)
    run("LoadParameters.m");
    Tchange = 5;      % every 5s
    margin  = 100e-3; % 100ms
    n_remove = margin/sens.enc.T_s;

    % importing data
    load(file);
    t  = wm_best_meas.time;
    n_sample = length(t);
    wm = deg2rad(wm_best_meas.signals.values);
    ia = ia_best.signals.values;
    vs = sens.curr.Rs * ia;

    % suppressing peaks
    filtered_vs = [];
    filtered_wm = [];
    filtered_t  = [];
    invalid_v   = 0;
    for i=1:n_sample
       if(mod(t(i),Tchange)==0)
           start = i - n_remove;
           stop  = i + n_remove;
           if(start<1)
               start = 1;
           end
           if(stop > n_sample)
               stop = n_sample;
           end
           for j = start:1:stop
               vs(j) = invalid_v;
               wm(j) = invalid_v;
           end
       end
    end
    j = 1;
    for i=1:n_sample
        if(vs(i)~=invalid_v)
            filtered_t(j)  = t(i);
            filtered_vs(j) = vs(i);
            filtered_wm(j) = wm(i);
            j = j+1;
        end
    end
    filtered_t = transpose(filtered_t);
    filtered_wm = transpose(filtered_wm);
    filtered_vs = transpose(filtered_vs);

    % least squares algorithm
    Xt = [filtered_wm, sign(filtered_wm)];
    Y = filtered_vs;
    X = transpose(Xt);
    theta_star = (X*Xt)\X*Y;
    a = theta_star(1,1);
    b = theta_star(2,1);
    Beq   = a*mot.Kt/sens.curr.Rs;
    tausf = b*mot.Kt*gbox.N/sens.curr.Rs;
end