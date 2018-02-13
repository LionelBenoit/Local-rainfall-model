function [] = plot_result_metropolis( V_M, num_fig )

[sx,sy]=size(V_M);

figure(num_fig)
clf;

subplot(4,6,1)
hist(V_M(:,1),50)

subplot(4,6,7)
plot(V_M(:,1),1:length(V_M(:,1)))

subplot(4,6,2)
hist(V_M(:,2),50)

subplot(4,6,8)
plot(V_M(:,2),1:length(V_M(:,2)))

subplot(4,6,3)
hist(V_M(:,3),50)

subplot(4,6,9)
plot(V_M(:,3),1:length(V_M(:,3)))

subplot(4,6,4)
hist(V_M(:,4),50)

subplot(4,6,10)
plot(V_M(:,4),1:length(V_M(:,4)))

subplot(4,6,5)
hist(V_M(:,5),50)

subplot(4,6,11)
plot(V_M(:,5),1:length(V_M(:,5)))

subplot(4,6,6)
hist(V_M(:,6),50)

subplot(4,6,12)
plot(V_M(:,6),1:length(V_M(:,6)))

%---

subplot(4,6,13)
hist(V_M(:,7),50)

subplot(4,6,19)
plot(V_M(:,7),1:length(V_M(:,7)))

subplot(4,6,14)
hist(V_M(:,8),50)

subplot(4,6,20)
plot(V_M(:,8),1:length(V_M(:,8)))

subplot(4,6,15)
hist(V_M(:,9),50)

subplot(4,6,21)
plot(V_M(:,9),1:length(V_M(:,9)))

subplot(4,6,16)
hist(V_M(:,10),50)

subplot(4,6,22)
plot(V_M(:,10),1:length(V_M(:,10)))

subplot(4,6,17)
hist(V_M(:,11),50)

subplot(4,6,23)
plot(V_M(:,11),1:length(V_M(:,11)))

end