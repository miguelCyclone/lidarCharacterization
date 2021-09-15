% Attribution-NonCommercial-ShareAlike 3.0 Unported (CC BY-NC-SA 3.0)
% This is the main script
% It is possible to run multiple simulations with different configurations by changing the code number

clear all;
close all;

%Boot parametric file
code = 2;
[Energy, PRFVar, primaryLens, eMF, eENF, rIF_BW, rMF, anodeDarkCurrent, anodeRadiantSens, vsm, LidarRatio, boundaryLayer, opertaionTimeVar] = parameters(code);
opertaionDayTime = opertaionTimeVar; % PARAMETER Operation time, night time or day time

% compute until end of troposhpere [km]
R_end = 15;
iterations = 1000;

% constants
c=2.99793*10^8; %speed of light [m/s]
h=6.6262*10^-34; % plancks constant [J*s]
q = 1.602*10^-19; % electron charge [C]

%laser
em_wL = 532*10^-9; % emission wavelength from [nm] to [m]
E=Energy*10^-3; %PARAMETER, laser energy in J=W*s from [mJ] to [J]
PRF=PRFVar;%PARAMETER pulse repetition frequency [Hz]

%telescope
d=primaryLens; %PARAMETER primary lens diamater [m]
dsh=0.06858; %shade diamater [m]
fl=2; %focal length [m]
Lt=0.6;   % losses for the telescope [.]


% ***********************************
% START ELASTIC receiving channel
% ***********************************

lambda_0 = 532*10^-9; % reception wavelength, from [nm] to [m]

%Interferance filter
if_bw=10; % interferance BW [nm]
Lif=0.65; % losses for the interferance filter [.]

%Photodiode
Dapd=3*10^-3; %Active area diameter from [mm] to [m]
I_ds = 7.64*10^-8; %dark signal that runs on the surface of teh lattice [A]
I_db = 3.10*10^-10; %dark signal that runs inside the lattice [A]
R_io = 240*10^-3; %intrinsic responsivity from [mA/W] to [A/W]
M = eMF; %PARAMETER multiplication factor [.]
F = eENF; %PARAMETER excess noise factor [.]

%Signal conditionning stages
G_t=5750; % transimpedance Gain, 1st stage, [ohms]=[V/I]
G_acn=20.3; % Voltage conditioning Gain (2nd Stage) [V/V]
B = 10*10^6; % Noise equivalent BW, from [MHz] to [Hz]
sigma_th_i = 5*10^-12; % equivalent input noise from [pA/sqrt(Hz)] to [A/sqrt(Hz)]

% Acqusition system
%analog to digital recorder ADC
f_s = 20*10^6; % sampling frequency from [M*sps] to [sps]

% *****************************
% End Elastic receiving channel
% *****************************


% *****************************
% START RAMEN receiving channel
% *****************************

lambda_r = 607.4*10^-9; % Reception wavelenght from [nm] to [m]

% Interferance filter ramen
if_bw_ramen= rIF_BW; % PARAMETER interferance BW in [nm], do not change! the formula for background is in nm
Lif_ramen=0.65; % losses for the interferance filter [.]

% Photo multiplier
D=3*10^-3; %equivalent area diameter from [mm] to [m]
I_da = anodeDarkCurrent*10^-9; % PARAMETER anode dark current from [nA] to [A]
R_i_ramen = anodeRadiantSens; % Anode radiant sensitivity [A/W]
M_ramen = rMF; % PARAMETER multiplication factor [.]
F_ramen = 1.8; %excess noise factor [.]

%Signal conditionning stages
PMT_load_resistance = 50; % G_t [ohms]
B_ramen = 10*10^6; % Noise equivalent BW from [MHz] to [Hz]

% Acqusition system
% Type: photon counter
t_pc = 1; % temporal resolution [bin]
binning = 50*10^-9; %sampling period from [ns/bin] to [s/bin]

% *****************************
% END RAMEN receiving channel
% *****************************

%Atmosphere aerosol (lambda_0)
Vm = vsm; %PARAMETER visibility margin at lambda_0 [km]
Sm = LidarRatio; % PARAMETER Lidar ratio [sr]

%Atmosphere aerosol (lambda_R)
k_wdc = 1.8; %wavelength dependant coefficient [.]

%molecular components at lambda_0 and lambda_R
alpha_mol = 0;%obtained from loop, (1.2593*10^-2)-(7.7599*10^-4)*R; % Rayleigh extinction [km^-1]
alpha_mol_lamdaR = 0;%obtained from loop, (7.3219*10^-3)-(4.5204*10^-4)*R; % Rayleigh extinction [km^-1]
S_r = (8*pi)/3; %Raylegh ratio

backscatterCrosssection = 3.71*10^-41; %N_2: Raman backscattering cross-section [km^2*sr^-1]
molecNumbDensity = 0;%obtained from loop, 2.1145*10^34 - 2.0022*10^33 * R + 5.4585*10^31 * R^2; %molecule number density [km^-3]

%atmosphere backround radiance component
Lradiation(1) = 3*10^-11; % MOON [W * cm^-2 * nm^-1 * sr^-1]
Lradiation(2) = 3*10^-6;  %SUN [W * cm^-2 * nm^-1 * sr^-1]

%bounder layer
Rpbl = boundaryLayer; % PARAMETER bounder layer [km]
%other parameters
Rovf = 0.2; %full-overlap factor [km]
SNR_min = 1; % Maximum-range criteration [dbV]
t_obs_simulate = 1/PRF; %[s]

% *****************************
% *****************************
% *****************************
%    - END OF PARAMETERS -
% *****************************
% *****************************
% *****************************

[K,area_ef] = systemConstant(E,d,dsh,c);

%Elastic receiver voltage chain
G_T =G_t*G_acn; % ELASTIC equivalent gain [V/I]
L_elastic = Lt * Lif; % overall losses not an acurate term, better optical transmissions
[Rv,Rv_net, Ri] = receiverChainVoltage(R_io , M, G_T, L_elastic);

%Ramen receiver voltage chain
G_T_ramen =PMT_load_resistance; % ELASTIC equivalent gain [V/I]
L_ramen = Lt * Lif_ramen; % overall losses not an acurate term, better optical transmissions
[Rv_ramen,Rv_net_ramen, Rio_ramen] = receiverChainVoltageRamen(R_i_ramen, M_ramen, G_T_ramen, L_ramen);

%aerosol components
[alpha_aer,beta_aer] = aerosolCoeffElastic(Vm,Sm);%
alpha_aer_lambdaR = aerosolCoeffRamen(alpha_aer,em_wL,lambda_r);%

% elastic sigma shot dark and sigma thermal are constant values
sigma_sh_d = sqrt(2*q*(G_T^2)*(I_ds+F*(M^2)*I_db));
sigma_th = sqrt((sigma_th_i^2)*(G_T^2));
N_sh_d_aux = (sigma_sh_d^2)*B; %auxiliar variables so we dont calculate it all the time inside the SNR
N_th_aux=(sigma_th^2)*B;
% Raman
x= 1;
I_db_Ramen = I_da / M_ramen;
%ramen_sigma_sh_d = sqrt(2*q*(G_T_ramen^2)*(0+F_ramen*(M_ramen^2)*I_db_Ramen)); %Ids=0
N_d = I_db/q; %Dark current humber of photona
N_d_prime = N_d*x; %dark current in counts/bin

% Plot
R = linspace(Rovf,R_end,iterations);% Full overlap from 200m onwards
P_elastic = zeros(iterations,0); %init vector for faster computations
P_ramen = zeros(iterations,0); % only molecule backscatter component
P_back_elastic = zeros(iterations,0);
P_back_ramen = zeros(iterations,0);
SNR_0 = zeros(iterations,0);
SNR_R = zeros(iterations,0);

for i = 1:length(R)
    %molecular components
    [alpha_mol, beta_mol] = molecCoeffElasic(S_r,R(i));
    [beta_mol_lambdaR, alpha_mol_lamdaR] = molecCoeffRamen(backscatterCrosssection,R(i));
    
    %Only molecs can answer for Ramen, so it has constant drop, no backscatter aerosol
    alphatot_lamnda0 = alpha_mol + alpha_aer;
    alphatot_lamndaR = alpha_mol_lamdaR + alpha_aer_lambdaR;
    P_ramen(i)= (K/R(i)^2)*1*(beta_mol_lambdaR)*exp(-R(i)*(alphatot_lamnda0 + alphatot_lamndaR));
    
    %Elastic returned power
    if R(i) <= Rpbl
        P_elastic(i) = (K/R(i)^2)*(beta_aer + beta_mol)*exp(-2*R(i)*(alpha_aer + alpha_mol));
       else 
        P_elastic (i)= (K/R(i)^2)*(beta_mol)*exp(-2*( Rpbl*(alpha_aer + alpha_mol) + (R(i) - Rpbl)*(alpha_mol) ));
    end
    
    % Background component
    P_back_elastic(i) = Pback(Lradiation(opertaionDayTime),area_ef,fl,Dapd,if_bw);
    P_back_ramen(i) = Pback(Lradiation(opertaionDayTime),area_ef,fl,D,if_bw_ramen);
    
    %SNR Elastic
    sigma_sh_s =  sqrt(2*q*(G_T^2)*F*(M^2)*R_io*(P_elastic(i)+P_back_elastic(i))*L_elastic); %[V^2/Hz]
    N_sh_s(i) = (sigma_sh_s^2)*B;
    N_sh_d(i) = N_sh_d_aux;
    N_th(i) = N_th_aux;
    N_total(i) = N_sh_s(i) + N_sh_d(i) + N_th(i);   
    SNR_0(i) = (Rv_net*P_elastic(i)) / (sqrt(N_total(i)));    
    
    %SNR Raman  
    %measuredtime= 1/PRF;%[s]
    %measuredtime=(1/(2*B_ramen));
    %measuredtime=binning*t_pc;
    %measuredtime=1/(binning*t_pc);
    
    N_back = (P_back_ramen(i)*L_ramen*Rio_ramen)/q; % number of photos from Pback component
    N_back_prime(i) = N_back*x; % background component [counts/s]
    
    N_d_prime_aux(i)=N_d_prime; % to create the plot
    
    N_ph = (P_ramen(i)*L_ramen*Rio_ramen)/q; %number of photons from returned power
    N_ph_prime(i) = N_ph*x; % returned power [counts/s]
    
    SNR_R(i) = N_ph_prime(i)*(sqrt(binning))/(sqrt(N_ph_prime(i)+2*(N_back_prime(i)+N_d_prime_aux(i))));
   
    %from linear to dB    
    P_elastic(i) = 10.*log10(P_elastic(i)); %pass p to dBW
    P_ramen(i) = 10.*log10(P_ramen(i));    
    SNR_R(i)=20.*log10(abs(SNR_R(i))); % Pass value to dbV
    SNR_0(i)=20.*log10(abs(SNR_0(i)));
    
end 

%Q7
SNRN = 2*10^3; %linear unit 
SNRN=20.*log10(abs(SNRN));
maxObserbationTime = 10^4; %s
[t_obs, maxRangeObsTime]=pulseIntegration(SNR_R, SNRN, PRF, maxObserbationTime, 20);

%Q8
sigma_photodiode_sh_d = sqrt(2*q*(I_ds+F*(M^2)*I_db)); %elastic
[NEP,NEP_s,quantum_eficiency] = nep(sigma_photodiode_sh_d,sigma_sh_d,sigma_th, Ri,Rv_net, M,q,h,c,em_wL);%elastic
R_sigma_photodiode_sh_d = sqrt(2*q*(0+F_ramen*(M_ramen^2)*I_db_Ramen)); %Raman
[R_NEP,R_NEP_s,R_quantum_eficiency] = nep(R_sigma_photodiode_sh_d,N_d_prime,0, R_i_ramen,Rv_net_ramen, M_ramen,q,h,c,em_wL);%raman

%Q9
tObs_simul = linspace(1/PRF,maxObserbationTime,iterations);
[SNR_E_aux]=question9(SNR_min, PRF, tObs_simul, Rv_net, P_elastic);
[SNR_R_aux]=question9(SNR_min, PRF, tObs_simul, Rv_net_ramen, P_ramen);
maxRange_E=maxRange(SNR_E_aux, R, SNR_min);
maxRange_R=maxRange(SNR_R_aux, R, SNR_min);

% to plot the SNR Y limit
if(SNR_0(1) > SNR_R(1))
    limitYMax = SNR_0(1);
else
    limitYMax = SNR_R(1);
end

titleLabel = "";
if opertaionDayTime == 1
    titleLabel = "Night";
end
if opertaionDayTime == 2
    titleLabel = "Day";
end



figure(1)
hold on
plot(R, P_elastic, 'b')
plot(R, P_ramen,'m')
plot(R, 10.*log10(P_back_elastic), '--b')
plot(R, 10.*log10(P_back_ramen),'--m')
xlim([Rovf R_end]) %make axis start at 200m instead of 0m
set(gca, 'XTick', sort([Rovf, get(gca, 'XTick')])); % add extra tick at 200m
set(gca, 'XTick', sort([get(gca, 'XTick'), R_end])); % add extra tick at 15km
xlabel('Range in [km]');
ylabel('Power in [dBW]');
legend('Returned power Elastic','Returned power Ramen', "Background Elastic", "Background Ramen");
title("Returned power of Elastic and Raman - "+titleLabel);

figure(2)
hold on
plot(R, SNR_0,'b');
plot(R, SNR_R,'m');
xlim([Rovf R_end]); %make axis start at 200m instead of 0m
ylim([-inf limitYMax+5]);
set(gca, 'XTick', sort([Rovf, get(gca, 'XTick')])); % add extra tick at 200m
set(gca, 'XTick', sort([get(gca, 'XTick'), R_end])); % add extra tick at Rend
xlabel('Range in [km]');
ylabel('SNR in [dBV]');
legend('SNR_elastic', "SNR_Raman");
title("SNR raman and Elastic - "+titleLabel);

figure(3)
hold on
plot(R, SNR_0,'b');
plot(R, 20.*log10(abs(N_sh_s)),'--g');
plot(R, 20.*log10(abs(N_sh_d)),'--black');
plot(R, 20.*log10(abs(N_th)),'--r');
xlim([Rovf R_end]); %make axis start at 200m instead of 0m
set(gca, 'XTick', sort([Rovf, get(gca, 'XTick')])); % add extra tick at 200m
set(gca, 'XTick', sort([get(gca, 'XTick'), R_end])); % add extra tick at Rend
xlabel('Range in [km]');
ylabel('SNR in [dBV]');
legend('SNR_elastic','N_sh_s','N_sh_d',"N_th");
title("SNR elastic");

figure(4)
hold on
plot(R, SNR_R,'m');
plot(R, 20.*log10(abs(N_ph_prime)),'--g');
plot(R, 20.*log10(abs(N_d_prime_aux)),'--black');
plot(R, 20.*log10(abs(N_back_prime)),'--r');
xlim([Rovf R_end]); %make axis start at 200m instead of 0m
set(gca, 'XTick', sort([Rovf, get(gca, 'XTick')])); % add extra tick at 200m
set(gca, 'XTick', sort([get(gca, 'XTick'), R_end])); % add extra tick at Rend
xlabel('Range in [km]');
ylabel('SNR in [dBV]');
legend('SNR_raman','N_sh_s','N_sh_d',"N_back");
title("SNR raman");

figure(5)
hold on
plot(R, t_obs,'m');
plot(R(maxRangeObsTime),maxObserbationTime,'r*')
xlim([Rovf R_end]); %make axis start at 200m instead of 0m
set(gca, 'XTick', sort([Rovf, get(gca, 'XTick')])); % add extra tick at 200m
set(gca, 'XTick', sort([get(gca, 'XTick'), R_end])); % add extra tick at Rend
xlabel('Range [km]');
ylabel('Observation time [s]');
legend('Ramen_observationTime', "Maximum distance");
title("Obersvation time Vs range to ensure SNR min - "+titleLabel);


figure(6)
hold on
plot(tObs_simul, maxRange_E,'b');
plot(tObs_simul, maxRange_R,'m');
ylim([Rovf R_end]); %make axis start at 200m instead of 0m
xlim([tObs_simul(1) tObs_simul(length(tObs_simul))]);
set(gca, 'YTick', sort([Rovf, get(gca, 'YTick')])); % add extra tick at 200m
set(gca, 'YTick', sort([get(gca, 'YTick'), R_end])); % add extra tick at Rend
xlabel('Observation Time [s]');
ylabel('Range [km]');
legend('Elastic', "Raman");
title("Question 9");

disp("********************");
disp("Code: "+code);
disp("Start Params: Energy: "+Energy+", PRFVar: "+PRFVar+", primaryLens: "+primaryLens+", eMF: "+eMF+", eENF: "+eENF+", rIF_BW: "+rIF_BW+", rMF: "+rMF+", anodeDarkCurrent: "+anodeDarkCurrent+", anodeRadiantSens: "+anodeRadiantSens+", vsm: "+vsm+", LidarRatio: "+LidarRatio+", boundaryLayer: "+boundaryLayer+", opertaionTimeVar: "+opertaionTimeVar);
disp("End init params");
disp(" ");
disp("System constant: "+K+"[W*km^3]");
disp("Effective area: "+area_ef+"[m^2]");
disp("Elastic, PbackComponent: "+P_back_elastic(1)+"[W]");
disp("Ramen, PbackComponent: "+P_back_ramen(1)+"[W]");
disp("Elastic, alpha_aerosol: "+alpha_aer+"[km^-1], backscatter_aerosol: "+beta_aer+"[km^-1*sr^-1]");
disp("Ramen, alpha_aerosol: "+alpha_aer_lambdaR+"[km^-1]");
disp("Elastic, Rv: "+Rv+"[V/W]"+", Rv_net:"+Rv_net+"[V/W]"+", Ri:"+Ri+"[A/W]");
disp("Ramen, Rv: "+Rv_ramen+"[V/W]"+", Rv_net:"+Rv_net_ramen+"[V/W]"+", Rio:"+Rio_ramen+"[A/W]");
disp("Ramen, Max range for time observation max of "+maxObserbationTime+" [s] is: "+R(maxRangeObsTime)+" [km]")
disp("Elastic: NEP in the photodiode: "+NEP+" [femto_W/sqrt(Hz)], with a quantum efficiency of: "+quantum_eficiency+"%, "+"NEP_system: "+NEP_s+" [femto_W/sqrt(Hz)]");
disp("Raman : NEP in the PMT: "+R_NEP+" [femto_W/sqrt(Hz)], with a quantum efficiency of: "+R_quantum_eficiency+"%");
disp("********************");
disp(" ");
