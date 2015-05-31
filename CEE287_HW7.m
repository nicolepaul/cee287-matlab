%% CEE 287: Homework 7

%% Pushover Curve

load HW7PushoverData
figure;
plot(d_npd,v_npd,'b-',d_pd,v_pd,'r--','LineWidth',2);
grid on;
xlabel('Roof Displacement [in]');
ylabel('Base Shear [kips]');
title('Nonlinear Static Pushover Curve');
legend('No P-\Delta','With P-\Delta','Location','best');

%% Drift curves

load HW7DriftData;
figure;
plot(Drift_NPD, Elev, 'b*-', Drift_PD, Elev, 'r*--', 'LineWidth', 2);
grid on;
xlabel('Interstory Drift Ratio');
ylabel('Elevation [ft]');
title('Interstory Drift Ratios');
legend('No P-\Delta','With P-\Delta','Location','best');
ylim([0 max(Elev)]);

