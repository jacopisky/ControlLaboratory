%run("LoadHubParameters.m");
load("hub_test.mat");
t = out.hub_test.time;
theta_b = abs(out.hub_test.signals(1).values());
prev = theta_b(1);
dir = 0;
count = 1;
peaks = [];
visual = zeros(size(t));
for i=2:size(t)
    now = theta_b(i);
    if dir == 1 && now < prev
        % record peak
        peaks(count) = theta_b(i-1);
        count = count + 1;
        dir = -1;
        visual(i) = theta_b(i);
    elseif dir == -1 && now > prev
        dir = 1;
    elseif dir == 0
        if prev < now
            dir = 1;
        else
            dir = -1;
        end
    end
    prev = now;
end
figure
plot(t, abs(theta_b))
hold on
plot(t, visual)
Y = transpose(log(peaks));
PHI = zeros(count-1,2);
for i=1:count-1
    PHI(i, 1) = i;
    PHI(i, 2) = 1;
end
THETA_HAT = (transpose(PHI)*PHI)\transpose(PHI)*Y;
a_hat = THETA_HAT(1,1);