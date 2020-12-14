
function sim_ecrv(N, z0, tau, dt, tF, mu, sig)

t = 0:dt:tF;

figure();
for i = 1:N
dt = dt;
tau = tau;
z(1) = z0;
for j = 1:length(t)-1
    z(j+1) = z(j) * exp(-dt/tau) + normrnd(mu, sig);
end
plot(t,z); hold on
end
end