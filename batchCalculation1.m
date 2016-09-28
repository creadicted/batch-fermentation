function [ t_save , cx ] = batchCalculation1

    mu_max      = 0.8;      % 1/h           Maximale Wachstumsrate
    mu_death    = 0.1;      % 1/h           Biomasse/Substratausbeute
    Y_xs        = 0.44;     % gx 1/gs       Maintenance Bedarf
    M_s         = 0.0057;   % gs gx 1/h     Accetat/Biomasseausbeute
    Y_xp        = 5;    % gs 1/L        Accetat/Biomassenausbeute
    K_s         = 0.004;    % gs 1/L        Substrataffinität
    K_02        = 0.0001;   % mmol 1/L      Sauerstoffaffinität
    K_la        = 1000;     % 1/h           Gasübergang
    Y_x02       = 0.10;     % gx 1/mmol     Biomasse/O_2 Aubeute
    y_02        = 0.2;      % 1/V_Luft      O2-Gehalt
    k_hpc       = 1057;     % L atm 1/mol   Henrykonstante O_2 
    
    % Aufgabenparameter
    c_x0 = 0.15; % initale Biomasse
    c_s0 = 7.5;  % initale Substratmenge
    c_P0 = 0.0;  %
    c_O2 = 0.0;  %
            
  % mu      = .8; % µ als Wachstumsrate    
    h       = .5; % Laufweite als Delta
    
    t_span(1)   =  0; % Anfangszeitpunkt
    t_span(2)   = 5; % Endzeitpunkt
       
    cx = c_x0;
    
    %****** Startwerte *******    
    % Biomasse
    start(1) = c_x0;
    % Substrat
    start(2) = c_s0;
    % Produkt
    start(3) = c_x0;
    % Sauerstoff
        % Berrechnung von c_02*
        c_O2s = (((mu_max * (c_s0/c_s0 + K_s)))/Y_x02)/k_hpc;
    start(4) = c_O2s;
    
    [t_save, cx] = ode15s(@RHS, t_span, start);
    
    function dx = RHS( t , x )
        % µ nach Monod 
        mu = mu_max * (x(2)/(x(2) + K_s));
        
        % Biomasse mit absterbephase
        dx(1) = (mu - mu_death) * x(1);
        
        % Substrat mit Maintenance
        dx(2) = -((mu/Y_xs)+M_s) * x(1);
        
        % Produkt / Acetat
        dx(3) = (mu/Y_xp) * x(1);              
        
        % Sauerstofverbrauch
        q_O2 = mu/Y_x02;
        c_02s = q_O2/k_hpc;
        dx(4) = K_la * (c_02s - x(4)) - (q_O2 * x(1));
        
        % Transponieren
        dx = transpose(dx);
    end
    %****** Plots *******    
    % Biomasse
    subplot(2,3,2);
    plotCx = plot(t_save,cx(:,1),'Color','g','DisplayName','C_x');    
    title('C_x');
    % Biomasse
    subplot(2,3,3);
    title('C_s');
    plotCs = plot(t_save,cx(:,2),'-.','Color','b','DisplayName','C_s');
    subplot(2,3,5)
    % Produkt
    plotCp = plot(t_save,cx(:,3),'--','Color','m','DisplayName','C_p');
    subplot(2,3,6)
    % Sauerstoff
    plotO2 = plot(t_save,cx(:,4),'--','Color','r','DisplayName','O^2')
       
    hL = legend([plotCx,plotCs,plotCp,plotO2],{'C_x','C_s','C_p','O^2'});
  
    newPosition = [0.1 0.65 0.2 0.2];
    newUnits = 'normalized';
    set(hL,'Position', newPosition,'Units', newUnits);

    
   % subplot(t_save,cx(:,1),'Color','g','DisplayName','C_x')
    
end
