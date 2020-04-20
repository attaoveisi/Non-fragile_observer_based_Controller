figure;
plot(controlleduncontrolled(:,1),controlleduncontrolled(:,2)/100)
hold on
plot(controlleduncontrolled(:,1),controlleduncontrolled(:,3)/100)
hold on
plot(controlleduncontrolled(:,1),controlleduncontrolled(:,4)/100)
xlabel('time (sec)')
ylabel(['measured voltage (1 V=125 m/sec)']);
legend('uncontrolled','Lien \it{et al}. (2007)','present')

figure;
subplot(2,1,1)
plot(effort1(:,1),effort1(:,3))
hold on
plot(effort1(:,1),effort1(:,2))
hold on
plot(disturbance(:,1),disturbance(:,2),'LineWidth',2)
xlabel('time (sec)')
ylabel('control signal (V)')
ylim([-60 50])
legend('present','Lien \it{et al}. (2007)','disturbance signal','Orientation','horizontal','Location','south')
title('first actuator')
subplot(2,1,2)
plot(effort2(:,1),effort2(:,3))
hold on
plot(effort2(:,1),effort2(:,2))
xlabel('time (sec)')
ylabel('control signal (V)')
title('second actuator')

figure;
plot(observation_error_Lien(:,1),observation_error_Lien(:,2))
hold on
plot(observation_error_present(:,1),observation_error_present(:,2))
xlabel('time (sec)')
ylabel('observation error (mm/sec)')
legend('Lien \it{et al}. (2007)','present')