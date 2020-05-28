%% (a)
clc; clear;

T = [3000 4000 5000];
loopsT = 3;

lambdaMax = []; %on the x axis

options = optimset( 'TolX' , 1e-10);

upperLimitLambda = (4.*10.^(-6)); 
lambdaVal = 0:upperLimitLambda./100:upperLimitLambda;

for i = 1:loopsT
   
    fun = @(lambda) -planck(lambda, T(1,i));
    plot(lambdaVal, -fun(lambdaVal));
    xlabel('');
    ylabel('Monochromatic Emittance')
    hold on
    max = fminbnd(fun, 0, 3.*10.^(-6), options);
    lambdaMax = [lambdaMax, max];
end

legend('Temperature: 3000 °C', 'Temperature: 4000 °C', 'Temperature: 5000 °C', 'Location', 'NorthEast');

disp("Lambda at which the monochromatic emittance is maximal")
disp(lambdaMax) %on the x axis
disp("Value of the monochromatic emittance at the maximal lambda")
disp(planck(lambdaMax,T)) % on the y axis, describes the monochromatic emittance

disp("Wiens displacement law: ")
b = [];
for i = 1:loopsT
   currentB = lambdaMax(1,i).*T(1,i);
   b = [b, currentB];
end
disp(b)
disp("Values in b should be the same if Wiens displacement law is true.")
%% (b)

clc; clear;

upperLimitX = (5.3); 
X = linspace(0,upperLimitX,300);
fOfX = @(x) (5-x).*exp(x) - 5; % the derivative of f, its max values

plot(X, fOfX(X));
hold on

yline(0, 'k'); % y = 0
xline(5, '--r'); %guess
xline(4.9651, '--g');
legend('f(x)','y = 0', 'Guess', 'Solution', 'Location', 'NorthWest');
xlabel('x = {hc}/{k?T}');
ylabel('{d(me_e)}/{d?}');

loops = 100;
h = 5;
for i = 1:loops
   
    h = h - (((5-h).*exp(h) - 5)./((4-h).*exp(h)));
    
end

xStar = h;
disp('The numerical solution to f(x) = 0 via Newtons method is; ')
disp(xStar)

h=6.6256e-34; c=2.9979e8; k=1.3805e-23;
b = (h.*c)./(k.*xStar);
disp('Where Wiens constant is equal to: ')
disp(b)
disp('The result for Wiens constant in (a) is 0.0029.')

%% (c)

clear; clc;

sigma = 5.67.*10.^(-8);

from = 4.*10.^(-7);
to = 7.*10.^(-7);

loopFrom = 100;
loopTo = 10000;
points = 1000;
step = (loopTo - loopFrom)./(points);

tic
T = [];
q = [];
for i = loopFrom:step:loopTo
    T = [T, i];
    fun = @(lambda) planck(lambda,i);
    s = integral(fun, from, to);
    e = planckC(sigma, i);
    quota = (s)./(e);
    q = [q, quota];
end
toc

set(0, 'defaultfigureposition', [800 500 900 400])
plot(T,q);

MaxQ = max(q);
index = find(q == MaxQ);
MaxT = T(index);

yline(MaxQ, '--y', 'color', [0.9294 0.6902 0.1294]);
xline(MaxT, '--y', 'color', [0.9294 0.6902 0.1294]);

hold on
plot(MaxT,MaxQ,'r*')

xlabel('Temperature in celcius (°C)');
ylabel('Ratio between the visiable and total emittance');
legend('{M_s(T)}/_{M_e(T)} function of °C','Max quota', 'Temperatur of max quota', 'Max point', 'Location', 'NorthWest');

temprM = num2str(MaxT); 
text(7200, 0.35, ['T_M = ' temprM ' °C']);

%% (d)

clear; clc;


sigma = 5.67.*10.^(-8);

lambdaSpan = [4.*10.^(-7) 7.*10.^(-7)];
y0 = [0];

loopFrom = 100;
loopTo = 10000;
points = 1000;
step = (loopTo - loopFrom)./(points);

tic
T = [];
q = [];
for i = loopFrom:step:loopTo
T = [T, i];    
[lambda, s] = ode45(@(lambda,m) planck(lambda,i), lambdaSpan, y0);
e = planckC(sigma, i);

last = size(s);
quota = (s(last(1,1)))./(e);
q = [q, quota];
end
toc

set(0, 'defaultfigureposition', [800 500 900 400])
plot(T,q);

MaxQ = max(q);
index = find(q == MaxQ);
MaxT = T(index);

yline(MaxQ, '--y', 'color', [0.9294 0.6902 0.1294]);
xline(MaxT, '--y', 'color', [0.9294 0.6902 0.1294]);

hold on
plot(MaxT,MaxQ,'r*')

xlabel('Temperature in celcius (°C)');
ylabel('Ratio between the visiable and total emittance');
legend('{M_s(T)}/_{M_e(T)} function of °C','Max quota', 'Temperatur of max quota', 'Max point', 'Location', 'NorthWest');

temprM = num2str(MaxT); 
text(7200, 0.35, ['T_M = ' temprM ' °C']);
