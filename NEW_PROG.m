    clear
    close all
    clc

    %NAKAGAMI-M
    %com o rho maior que 1 todas as variaveis vao para infinito
    rho = [0.0001:0.01:0.9999];
    lambda = 1  ;            %atraso do feedback
    omega = 1;
    Rb = 10e3;
    B = 10e3; 
    %M = 1;
    M = [1 2 3];             %ordem de desvanecimento   
    NodB = -204;             %dB
    No = 10^(NodB./10);      %No linear
    Ml = 20;                 %Ml = 20dB
    Ml_n = 10^(Ml./10);
    Gl = 5;                  %Gl = 5dB
    Gl_n = 10^(Gl./10);
    fc = 2.5e9;              %fc=2.5GHz
    alpha = 2.5;               %path-loss exponent
    laco = 50:10:500;
    SNR = -50:0.1:150;       %SNR em dB
    SNR_n = 10.^(SNR./10);
    Ptx = 97.9e-3;           %100mW   
    Prx = 112e-3;            %100mW
    n = 0.35;                %eficiencia do amplificador de potencia
    new_l = 2;               %new_l = 2, pois n�o repetimos transmissao de pacotes
    tax = 1;                 %taxa unica
    n_out = 1:6;             %numero de outages consecutivas                    
    
for t = 1:length(rho)
    t
    for l = 1:length(M)
        for d = 1:length(laco)
            for j = 1:length(tax)
            R(j) = tax(j)*Rb/B;
            z1(j) = (2^(R(j)) - 1);
            m = M(l);
            Pr = ((Gl_n.*(3e8/fc)^2).*SNR_n.*(No*B))./(((4*pi).^2).*Ml_n.*(laco(d).^alpha));
            L_grande = new_l;
                for h = 1:length(SNR)
                    for b = 1:length(L_grande)
                    c = zeros(1,length(L_grande));
                    f = ones(1,length(L_grande));
                        for a = 1:L_grande(b)
                        c_new(b) = ((rho(t))^(2*(a+lambda-1)))/(1-((rho(t))^(2*(a+lambda-1))));
                        c(b) = c(b) + c_new(b);
                        w(b) = 1+c(b);
                        f_new(b) = (1-((rho(t))^(2*(a+lambda-1))));
                        f(b) = f(b)*f_new(b);
                        l_new(b) = (w(b)*f(b))^(-m);
                        up(b) = (m^(m*L_grande(b)))*l_new(b)*(z1(j)^(L_grande(b)*m));
                        %down(b) = (gamma(m+1)^L_grande(b))*(((Pr(h)/(No*B))^L_grande(b)^m)); 
                        down(b) = (gamma(m+1)^L_grande(b))*(((SNR_n(h))^L_grande(b)^m));
                        pout_NAKA(h) = up(b)/down(b);
                        end
                    end
                    value = NaN;
                    pout(h) = (pout_NAKA(h));
                        for r = 1:length(n_out)
                            if pout(h) < ((1e-5)^n_out(r))
                                pout_new(r,h) = pout(h);
                                SNR_all(r,h) = SNR_n(h);
                                else
                                pout_new(r,h) = value;
                                SNR_all(r,h) = value;
                            end
                        end
                    down_pt =(((4*pi).^2).*Ml_n.*(laco(d).^alpha));
                    up_pt(:,h) = down_pt.*(SNR_all(:,h).*(No*B));
                    Pt_all(:,h) = up_pt(:,h)./(Gl_n.*(3e8/fc)^2);
                    E_up(:,h) = (Pt_all(:,h)./n) + Ptx + Prx;
                end
                E = E_up./R(j);
                E_new(:,:,j) = E;
            end
            for q = 1:length(n_out)
                [E_col index_row]  = min(E_new(q,:,:));
                [E_min(d,l,t), col(d,l,t)] = min(E_col);
                row(d,l,t) = index_row(col(d,l,t));
                E_end(q,d,t,l) = E_new(q,row(d,l,t),col(d,l,t));
            end
        end
    end
end

for w = 1:length(rho)
    E1(w) = E_end(2,1,w,1);
    E2(w) = E_end(2,1,w,2);
    E3(w) = E_end(2,1,w,3);
end

%FAZENDO D = 50, TAX = 1
figure(1) 
semilogy(n_out,E_end(:,1,100,1),'LineWidth',2)  
hold all
semilogy(n_out,E_end(:,6,100,1),'-*','LineWidth',2)  
semilogy(n_out,E_end(:,16,100,1),'-x','LineWidth',2)  
semilogy(n_out,E_end(:,26,100,1),'-s','LineWidth',2)  
semilogy(n_out,E_end(:,36,100,1),'-.','LineWidth',2)  
semilogy(n_out,E_end(:,46,100,1),'--','LineWidth',2)  
grid on
hold off 
leg1 = xlabel('Consecutive Outage')
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15); 
leg2 = ylabel('Energy Consumption Per Bit ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);  
leg3 = legend('d = 50m','d = 100m','d = 200m','d = 300m','d = 400m','d = 500m','Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10);  


figure(2) 
semilogy(laco,E_end(2,:,1,1),'LineWidth',2)  
hold all
semilogy(laco,E_end(2,:,1,2),'--','LineWidth',1.5)  
semilogy(laco,E_end(2,:,1,3),'-.','LineWidth',1.5)
semilogy(laco,E_end(2,:,100,1),'-*','LineWidth',2)  
grid on
hold off
leg1 = xlabel('Dist\^ancia ($d_\mathrm{ca}$) [m]')
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15); 
leg2 = ylabel('Energia Consumida por Bit ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);  
leg3 = legend('m = 1 e $\rho$ = 0','m = 2 e $\rho$ = 0','m = 3 e $\rho$ = 0','m = 1 e $\rho$ = 1','Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10);  


figure(3) 
loglog(rho,E1,'LineWidth',2)  
hold all
loglog(rho,E2,'--','LineWidth',2)  
loglog(rho,E3,'-.','LineWidth',2)  
axis([0.7 0.99 0.2 0.5])
grid on
hold off
leg1 = xlabel('Correlac\~ao no tempo ($\rho$)')
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15); 
leg2 = ylabel('Energia Consumida por Bit ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);  
leg3 = legend('m = 1','m = 2','m = 3','Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10);  

figure(4) 
loglog(rho,E2,'--','LineWidth',2)
hold on
loglog(rho,E1,'LineWidth',2)  
loglog(rho,E3,'-.','LineWidth',2)
axis([0.7 0.99 0.2 0.5])
grid on
leg1 = xlabel('Time Correlation ($\rho$)')
set(leg1,'Interpreter','latex');
set(leg1,'FontSize',15); 
leg2 = ylabel('Energy Consumption Per Bit ($E_\mathrm{b}$) [W]')
set(leg2,'Interpreter','latex');
set(leg2,'FontSize',15);  
leg3 = legend('m = 2','m = 1','m = 3','Location','best')
set(leg3,'Interpreter','latex');
set(leg3,'FontSize',10);
hold off
    
axes('position',[0.18 0.3 0.8 0.8])
box on

your_index= 0.9<rho & rho<0.999;

ax = gca;
ax.ColorOrderIndex = 2;
loglog(rho(your_index),E2(your_index),'--','LineWidth',2)  
hold on
ax.ColorOrderIndex = 3;
loglog(rho(your_index),E3(your_index),'-.','LineWidth',2)
axis tight
grid on
hold off
  
