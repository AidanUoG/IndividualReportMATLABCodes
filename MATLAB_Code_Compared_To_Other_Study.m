a1 = 0.18; K1 = 420; v1 = 18;
a2 = 0.25; b2 = 0.02;
a3 = 0.25; b3 = 0.024; gamma = 0.9;
Initial_Tumour_Size = 0.26;
% All values taken from Shabanisamghabady (2016)

% Murine Data (from Shabanisamghabady (2016))
TVI = 0.26;
TV21 = 1.6338;
TV27 = 27.6923;
TV35 = 144.6923;

f = @(t,y) [a1*y(1)*(1 - (y(1)/K1)^v1);
            a2*y(1)*exp(-b2*t);
            a3*y(1)^gamma - b3*y(1)]; % y(1) is the x variable (tumour size)
% If exact solution not known, set function below to return a zero vector
exact = @(t) [0;0;0];

h=0.25; % Step size
t0 = 1; tN = 35; % range of t

y0 = [Initial_Tumour_Size; Initial_Tumour_Size; Initial_Tumour_Size];

n = 3; % Number of unknown functions of y_i

if length(y0) ~= n
    disp('Error, you entered an incorrect number of equations. Try again ')
    return;
end

% Currently assumes no more than 10-dimensional problem for plotting purposes. Change if needed.
if length(y0) > 10
    disp('Error, for plotting purposes only up to 10 equations supported.')
    disp('Edit the file to change this and try again.')
    return;
end

t = t0:h:tN;
sizet = length(t);

y = zeros(n,sizet);
yexact = zeros(n, sizet);
ytemp = zeros(n,1);

% Initial conditions set in following loop
for k = 1:n
    y(k,1) = y0(k);
end

% Main Heun's method loop
for k = 2:sizet
    ytemp = y(:,k-1) + h*f(t(k-1), y(:,k-1));
    y(:,k) = y(:,k-1) + (h/2)*(f(t(k-1), y(:,k-1)) + f(t(k), ytemp));
end

for k = 1:sizet
    yexact(:,k) = exact(t(k));
end

for mm = 1: n
    fprintf('PRINTING INFORMATION FOR FUNCTION %d\n\n',mm);
    fprintf(' i    TIME       Yi (HEUN)        y(ti) (EXACT)      ABS. ERROR\n')
    for k = 1:sizet
        fprintf('%3d   %8f   %10f      %10f        %10f\n',k-1,t(k),y(mm,k),yexact(mm,k), abs(y(mm,k)-yexact(mm,k)))
        if mm == 1
            if t(k) == 1
                AbsError_1_1 = TVI - y(mm,k);
                RelError_1_1 = AbsError_1_1/y(mm,k);
                fprintf('\nAt t = 1, the absolute error is %f\10 and the relative error is %f\10 for logistic model. \n', AbsError_1_1, RelError_1_1);
            elseif t(k) == 21
                AbsError_1_21 = TV21 - y(mm,k);
                RelError_1_21 = AbsError_1_21/y(mm,k);
                fprintf('\nAt t = 21, the absolute error is %f\10 and the relative error is %f\10 for logistic model. \n', AbsError_1_21, RelError_1_21);
            elseif t(k) == 27
                AbsError_1_27 = TV27 - y(mm,k);
                RelError_1_27 = AbsError_1_27/y(mm,k);
                fprintf('\nAt t = 27, the absolute error is %f\10 and the relative error is %f\10 for logistic model. \n', AbsError_1_27, RelError_1_27);
            elseif t(k) == 35
                AbsError_1_35 = TV35 - y(mm,k);
                RelError_1_35 = AbsError_1_35/y(mm,k);
                fprintf('\nAt t = 35, the absolute error is %f\10 and the relative error is %f\10 for logistic model. \n', AbsError_1_35, RelError_1_35);
            end
        elseif mm == 2
            if t(k) == 1
                AbsError_2_1 = TVI - y(mm,k);
                RelError_2_1 = AbsError_2_1/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for gompertz model. \n', AbsError_2_1, RelError_2_1);
            elseif t(k) == 21
                AbsError_2_21 = TV21 - y(mm,k);
                RelError_2_21 = AbsError_2_21/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for gompertz model. \n', AbsError_2_21, RelError_2_21);
            elseif t(k) == 27
                AbsError_2_27 = TV27 - y(mm,k);
                RelError_2_27 = AbsError_2_27/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for gompertz model. \n', AbsError_2_27, RelError_2_27);
            elseif t(k) == 35
                AbsError_2_35 = TV35 - y(mm,k);
                RelError_2_35 = AbsError_2_35/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for gompertz model. \n', AbsError_2_35, RelError_2_35);
            end
        elseif mm == 3
            if t(k) == 1
                AbsError_3_1 = TVI - y(mm,k);
                RelError_3_1 = AbsError_3_1/t(k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for von bertalanffy model. \n\n', AbsError_3_1, RelError_3_1);
            elseif t(k) == 21
                AbsError_3_21 = TV21 - y(mm,k);
                RelError_3_21 = AbsError_3_21/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for von bertalanffy model. \n\n', AbsError_3_21, RelError_3_21);
            elseif t(k) == 27
                AbsError_3_27 = TV27 - y(mm,k);
                RelError_3_27 = AbsError_3_27/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for von bertalanffy model. \n\n', AbsError_3_27, RelError_3_27);
            elseif t(k) == 35
                AbsError_3_35 = TV35 - y(mm,k);
                RelError_3_35 = AbsError_3_35/y(mm,k);
                fprintf('The absolute error is %f\10 and the relative error is %f\10 for von bertalanffy model. \n\n', AbsError_3_35, RelError_3_35);
            end
        end
    end
end

clf %clear all figures
cmap = zeros(10,3); % Currently assumes no more than 10-dimensional problem. Change if needed.
cmap(1:3,1:3) = eye(3);
cmap(4,:) = [1 1 0]; cmap(5,:) = [1 0 1]; cmap(6,:) = [0 1 1];
cmap(7,:) = [0.5 1 0]; cmap(8,:) = [1 0 0.5]; cmap(9,:) = [0.5 0.5 0.5];
set(gca,'fontsize',14)
hold on

for k = 1:n
    if k == 1
        legentry = sprintf('Generalised Logistic Model');
    elseif k == 2
        legentry = sprintf('Gompertz Model');
    elseif k == 3
        legentry = sprintf('Von Bertalanffy Model');
    end
    plot(t,y(k,:),'Color',cmap(k,:),'linewidth',2,'DisplayName',legentry)
    xlabel('Time (days)'); ylabel('Tumour Volume (mm^3)');
    fprintf('PLOTTING INFORMATION FOR FUNCTION %d\n\n',k);
end
legend(gca,'show');
print('Tumour Models 1.jpeg','-djpeg');
hold off

LogModAbs = [AbsError_1_1; AbsError_1_21; AbsError_1_27; AbsError_1_35];
GomModAbs = [AbsError_2_1; AbsError_2_21; AbsError_2_27; AbsError_2_35];
VonModAbs = [AbsError_3_1; AbsError_3_21; AbsError_3_27; AbsError_3_35];
Days = [1, 21, 27, 35];

x = Days;
y = [LogModAbs, GomModAbs, VonModAbs];
plot(x,y);
title('Plot of absolute error of the three models');
xlabel('Time (days)'); ylabel('Absolute Error');
legend('Generalised Logistic Model', 'Gompertz Model', 'Von Bertalanffy Model');
print('Absolute Error','-djpeg');

LogModRel = [RelError_1_1; RelError_1_21; RelError_1_27; RelError_1_35];
GomModRel = [RelError_2_1; RelError_2_21; RelError_2_27; RelError_2_35];
VonModRel = [RelError_3_1; RelError_3_21; RelError_3_27; RelError_3_35];
Days = [1, 21, 27, 35];

x = Days;
y = [LogModRel, GomModRel, VonModRel];
plot(x,y);
title('Plot of relative error of the three models');
xlabel('Time (days)'); ylabel('Relative Error');
legend('Generalised Logistic Model', 'Gompertz Model', 'Von Bertalanffy Model');
print('Relative Error','-djpeg');