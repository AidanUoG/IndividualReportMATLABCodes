a2 = 0.25; b2 = 0.02;
Initial_Tumour_Size = 0.26;
% All values taken from Shabanisamghabady (2016)

f = @(t,y)  a2*y(1)*exp(-b2*t);
% y(1) is the x variable (tumour size)
% If exact solution not known, set function below to return a zero vector
exact = @(t) 0;

h=0.25; % Step size
t0 = 1; tN = 500; % range of t

y0 = Initial_Tumour_Size;

n = 1; % Number of unknown functions of y_i

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
        legentry = sprintf('Gompertz Model');
    end
    plot(t,y(k,:),'Color',cmap(k,:),'linewidth',2,'DisplayName',legentry)
    xlabel('Time (days)'); ylabel('Tumour Volume (mm^3)');
    fprintf('PLOTTING INFORMATION FOR FUNCTION %d\n\n',k);
end
legend(gca,'show');
print('Tumour Model GM.jpeg','-djpeg');
hold off