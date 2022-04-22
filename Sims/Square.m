clear all
close all
clc

d = 0.4;
a = linspace(0,2*pi,100);

A = linspace(-d/2,d/2,100);
B = linspace(-d/2,d/2,100);

dist = [A; B];

%points = [(A./cos(a)),(B./cos(a))];
