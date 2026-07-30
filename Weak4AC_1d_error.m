clc;
clear;
close all;
rng(100,'twister');

rounds = 20;

N = 2^12;
Mref = 2^16;
T = 10;

tau = T/Mref;
Ms = 2.^(9:14);
stepsizes = T./Ms;

sigma = 1;
alpha = 1/4;
beta1 = 1;
beta2 = 1;
theta = 1;
rho = 1;
PhiN = N^2;

modeIndex = (1:N)';
Lambda = -pi^2*modeIndex.^2;

Eref = exp(Lambda*tau);
Ecoarse = exp(Lambda*stepsizes);

U0 = zeros(N,1);
U0(1) = 1/sqrt(2);

Diff1 = zeros(1,length(Ms));
Diff2 = zeros(1,length(Ms));
Diff3 = zeros(1,length(Ms));

Rfinest = Mref/Ms(end);

for sample = 1:rounds
    sampleTimer = tic;

    X = U0;
    Y = repmat(U0,1,length(Ms));

    accumulatorFinest = zeros(N,1);
    accumulatorCoarse = zeros(N,length(Ms)-1);

    for m = 1:Mref
        dWfine = sqrt(tau)*randn(N,1);

        x = sqrt(2)*dst(X);

        denominatorX = ...
            (1 + (beta1*tau^theta + beta2*PhiN^(-rho)) ...
            .*abs(x).^(2/alpha)).^alpha;

        fx = (sigma*x-x.^3)./denominatorX;

        fX = idst(fx)/sqrt(2);

        X = Eref.*(X + tau*fX + dWfine);

        accumulatorFinest = accumulatorFinest + dWfine;

        if mod(m,Rfinest)==0
            finestStep = m/Rfinest;
            deltaW = accumulatorFinest;
            accumulatorFinest(:) = 0;

            j = length(Ms);
            y = sqrt(2)*dst(Y(:,j));

            denominatorY = ...
                (1 + (beta1*stepsizes(j)^theta + beta2*PhiN^(-rho)) ...
                .*abs(y).^(2/alpha)).^alpha;

            fy = (sigma*y-y.^3)./denominatorY;

            fY = idst(fy)/sqrt(2);

            Y(:,j) = Ecoarse(:,j).*(Y(:,j) + stepsizes(j)*fY + deltaW);

            for j = length(Ms)-1:-1:1
                accumulatorCoarse(:,j) = accumulatorCoarse(:,j) + deltaW;

                if mod(finestStep,2^(length(Ms)-j))==0
                    deltaW = accumulatorCoarse(:,j);
                    accumulatorCoarse(:,j) = 0;

                    y = sqrt(2)*dst(Y(:,j));

                    denominatorY = ...
                        (1 + (beta1*stepsizes(j)^theta + beta2*PhiN^(-rho)) ...
                        .*abs(y).^(2/alpha)).^alpha;

                    fy = (sigma*y-y.^3)./denominatorY;

                    fY = idst(fy)/sqrt(2);

                    Y(:,j) = Ecoarse(:,j).*(Y(:,j) + stepsizes(j)*fY + deltaW);
                else
                    break;
                end
            end
        end
    end

    normX = norm(X);
    phiX1 = sin(normX);
    phiX2 = cos(normX^2-pi/2);
    phiX3 = exp(-normX^2);

    for j = 1:length(Ms)
        normY = norm(Y(:,j));
        phiY1 = sin(normY);
        phiY2 = cos(normY^2-pi/2);
        phiY3 = exp(-normY^2);

        Diff1(j) = Diff1(j) + (phiX1-phiY1);
        Diff2(j) = Diff2(j) + (phiX2-phiY2);
        Diff3(j) = Diff3(j) + (phiX3-phiY3);
    end

    if mod(sample,10)==0
        fprintf('Round: %d / %d, elapsed: %.2f s\n', ...
            sample,rounds,toc(sampleTimer));
    end
end

Error1 = abs(Diff1/rounds);
Error2 = abs(Diff2/rounds);
Error3 = abs(Diff3/rounds);

disp('Error1 = ');
disp(Error1);
disp('Error2 = ');
disp(Error2);
disp('Error3 = ');
disp(Error3);

figure;
loglog(stepsizes,Error1,'o-r','LineWidth',1.5);
hold on;
loglog(stepsizes,Error2,'^-g','LineWidth',1.5);
loglog(stepsizes,Error3,'x-m','LineWidth',1.5);

loglog(stepsizes,stepsizes.^0.25/4,'c--','LineWidth',1);
loglog(stepsizes,stepsizes.^0.5/2,'b--','LineWidth',1);
loglog(stepsizes,5*stepsizes,'k--','LineWidth',1);

axis tight;
grid on;

h = legend( ...
    '$\varphi(\cdot)=\sin(\|\cdot\|)$', ...
    '$\varphi(\cdot)=\cos(\|\cdot\|^2-\frac{\pi}{2})$', ...
    '$\varphi(\cdot)=\exp(-\|\cdot\|^2)$', ...
    'Order 0.25', ...
    'Order 0.5', ...
    'Order 1.0', ...
    'Interpreter','latex', ...
    'Location','best');
set(h,'Interpreter','latex');

xlabel('Time stepsizes');
ylabel('Weak errors');
title('Temporal weak errors for the one-dimensional scheme','FontSize',13);
hold off;

Design = [ones(length(Ms),1),log(stepsizes)'];

rhs1 = log(max(Error1,realmin))';
rhs2 = log(max(Error2,realmin))';
rhs3 = log(max(Error3,realmin))';

sol1 = Design\rhs1;
sol2 = Design\rhs2;
sol3 = Design\rhs3;

q1 = sol1(2)
q2 = sol2(2)
q3 = sol3(2)

resid1 = norm(Design*sol1-rhs1)
resid2 = norm(Design*sol2-rhs2)
resid3 = norm(Design*sol3-rhs3)
