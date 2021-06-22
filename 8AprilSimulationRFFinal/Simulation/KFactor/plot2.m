clc

phaseD1_prime = exp(1i*phaseD1);

%%
figure

for i = 1: 1000

    plot(real(phaseD1_prime(i)),imag(phaseD1_prime),'o','LineWidth',2);
    
    pause(2)
end
