% known parameters
Kt = 7.68e-3;
Rs = 0.5;
N  = 14;
Tsample = 1e-3;
deg2rad    = pi/180;

Tchange = 5; % every 5s
margin  = 100e-3; % 100ms
n_remove = margin/Tsample;

% importing data
load("data/least_square_positive.mat");
t  = wm_best_meas.time;
n_sample = length(t);
wm = wm_best_meas.signals.values * deg2rad;
ia = ia_best.signals.values;
vs = Rs * ia;

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

% plotting results
nexttile
plot(filtered_t, filtered_vs)
title('Voltage')
nexttile
plot(filtered_t, filtered_wm)
title('Angular velocity')

% least squares algorithm
Xt = [filtered_wm, sign(filtered_wm)];
Y = filtered_vs;
X = transpose(Xt);
theta_star = (X*Xt)\X*Y;
a = theta_star(1,1);
b = theta_star(2,1);
Beq = a*Kt/Rs;
tausf = b*Kt*N/Rs;

% found values
sprintf(num2str(Beq))
sprintf(num2str(tausf))
