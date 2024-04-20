clear all
close all
clc

load distributed_localization_data.mat


figure(1)
F1=digraph(Q_8');
plot(F1)

figure(2)
F2=digraph(Q_12'); 
plot(F2)

figure(3)
F3=digraph(Q_18'); 
plot(F3)

figure(4)
F4=digraph(Q_4')
plot(F4)