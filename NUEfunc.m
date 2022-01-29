function [NUE, leach] = NUEfunc(soils, precip, apptime1, mass, simlength)
%%
%   Corn N-NUE Simulation Model, 2022 
%   Patrick J. Dunn, Leanne M. Gilbertson
%   University of Pittsburgh
%   email: pjd66@pitt.edu
%   
%   References: 
%   

%   See README for more information on variables and model function.

%   Inputs: soils: a 7x1 vector conatining the following soil parameters in
%           this order [Ksat(cm/day), FC(cm^3/cm^3), WP(cm^3/cm^3),
%           SAT(cm^3/cm^3), N, a(cm^-1), SATinit(cm^3/cm^3)]
%           
%           precip: a vector containing the hourly precipitation in cm^3.
%           Length must be equivalent to (simulation time * 24 hours)
        
%           apptime1: integer that is the number of days before or after 
%           germination that application of fertilizer occurs

%           mass: the amount of fertilizer applied in (mg/cm^2)
          
%           simlength: the length of the simulation in days    
%  
%           EX) soils = [25, 0.3, 0.15, 0.45, 1.15, 0.01, 0.2]
%               precip = 0.1*ones(1,24*90)
%               apptime1 = [1, 14]                
%               mass = 1.5
%               simlength = 90
%
%   Outputs: NUE: the cumulative amount of nitrate recovered by the root
%           system divided by the amount applied. Can obtain the cumulative
%           amount recovered by multiplying this value by NUE*mass*X*Y,
%           where X and Y are the surface dimensions.

%           leach: the cumulative amount of nitrate that has leached
%           from the study system at the end of simulation (mg)
%
%           Model simulates the transport, uptake, and leaching of nitrate
%           in 3D due to subsurface drainage of water from precipitation.
%           Root files must be in folder with this file for
%           function to operate. Root files give information
%           on the dynamic root architecture and physiological root uptake
%           of nitrate.
%           
%           Study plot is 18X24X50cm
%%
application = mass;
apptime = apptime1;

%% Initial Conditions
load("ROOTfiles");
%%
SATly = soils(4);
KS = soils(1);
FC = soils(2);
WP = soils(3);
TT = (SATly - FC)/(KS/24);
P = zeros(1,24*simlength);
P(1:24*simlength) = precip(:,1);
%%
n = 50;
r1 = 18;
r2 = 24;

timelength = simlength;
NitLay(2:n-1,1) = ones(1,n-2)*0; %mass of N in layer
NitMob = zeros(n-1,r1-1,r2-1,timelength,24); %mass of mobile N in layer
nitratel(:,:,:,1,:) = (repmat(NitLay(:,1), 1, r1-1, r2-1,24));
nitratel(:,:,:,1,:) = nitratel(:,:,:,1,:);

SWlyex = zeros(n-1,r1-1, r2-1,timelength);
%WATER UPTAKE
Kx =  0.0432/24; % cm^3/hour from Javaux 2008
Kr = 1.728E-4/24; %per hour from Javaux 2008
prw = 10197.2; %cm
N = soils(5);
a = soils(6);
M = 1-1/N;
mapsurf0 = zeros(n-1, r1-1, r2-1, timelength);
mapsurf00 = zeros(n-1, r1-1, r2-1, timelength);
map0 = zeros(n-1, r1-1, r2-1, timelength);
map00 = zeros(n-1, r1-1, r2-1, timelength);
maprad0 = zeros(n-1, r1-1, r2-1, timelength);
maprad00 = zeros(n-1, r1-1, r2-1, timelength);
mapvol0 = zeros(n-1, r1-1, r2-1, timelength);
mapvol00 = zeros(n-1, r1-1, r2-1, timelength);
SWly = soils(7)*ones(n-1, r1-1, r2-1, timelength,24);
Wperc = zeros(n-1, r1-1, r2-1, timelength,24);
uptake = zeros(n-1, r1-1, r2-1, timelength,24);
S = zeros(n-1, r1-1, r2-1, timelength,24);
agefact = zeros(n-1, r1-1, r2-1, timelength);
wuptake = zeros(n-1, r1-1, r2-1, timelength,24);
h = zeros(n-1, r1-1, r2-1, timelength,24);
psw = zeros(n-1, r1-1, r2-1, timelength,24);
A = zeros(n-1, r1-1, r2-1, timelength,24);
before = 0;
%% SIMULATION
for t = 2:timelength
    % SWAT INFLITRATION
    if FC > SATly
        continue
    end

    
    mapsurf0(:,:,:,t) = permute(mapsurf(:,4:20,:,t), [3,2,1,4]);
    mapsurf00(:,:,:,t) = flip(mapsurf0(:,:,:,t),1);
    map0(:,:,:,t) = permute(map(:,4:20,:,t), [3,2,1,4]);
    map00(:,:,:,t) = flip(map0(:,:,:,t),1);
    maprad0(:,:,:,t) = permute(maprad(:,4:20,:,t), [3,2,1,4]);
    maprad00(:,:,:,t) = flip(maprad0(:,:,:,t),1);
    mapvol0(:,:,:,t) = permute(mapvol(:,4:20,:,t), [3,2,1,4]);
    mapvol00(:,:,:,t) = flip(mapvol0(:,:,:,t),1);
    ageprim0(:,:,:,t) = ageprim(:,4:20,:,t);
    agelat0(:,:,:,t) = agelat(:,4:20,:,t);
    mapage00(:,:,:,t) = agelat0(:,:,:,t)+ageprim0(:,:,:,t);

    simtime = t;
    %Soil Hydraulic Properties
    %Soil Water Pressure
    for hour = 1:24
        totalhour = (t-1)*24 + hour;
        dayuptake(t) = sum(uptake(:,:,:,t,:),'all');
        alluptake = sum(dayuptake(:),'all');
        Cint = alluptake/(totalrootvol(t));
        Cgoal = (massdemand(t))/(totalrootvol(t));
        for x = 1:n-1
            for y = 1:r1-1
                for z = 1:r2-1
                    if x == 1
                        if hour == 1
                            SWly(1,y,z,t,hour) = P(totalhour)  - Wperc(1,y,z,t-1,24) + SWly(1,y,z,t-1,24)-wuptake(1,y,z,t-1,24);
                            if SWly(1,y,z,t,hour) < WP
                                SWly(1,y,z,t,hour)= WP;
                            else
                            end
                            if SWly(1,y,z,t,hour) <= FC
                                SWlyex(1,y,z,t,hour) = 0;
                            else
                                SWlyex(1,y,z,t,hour) = SWly(1,y,z,t,hour) - FC;
                                if SWlyex(1,y,z,t,hour) > SATly - FC
                                    SWlyex(1,y,z,t,hour) = SATly - FC;
                                end
                            end
                            Wperc(1,y,z,t,hour) =  SWlyex(1,y,z,t,hour)*(1-exp(-1/TT));
                            
                        else
                            SWly(1,y,z,t,hour) =  P(totalhour) - Wperc(1,y,z,t,hour-1) + SWly(1,y,z,t,hour-1) - wuptake(1,y,z,t,hour-1);
                            if SWly(1,y,z,t,hour) < WP
                                SWly(1,y,z,t,hour)= WP;
                            else
                            end
                            if SWly(1,y,z,t,hour) <= FC
                                SWlyex(1,y,z,t,hour) = 0;
                            else
                                SWlyex(1,y,z,t,hour) = SWly(1,y,z,t,hour) - FC;
                                if SWlyex(1,y,z,t,hour) > SATly - FC
                                    SWlyex(1,y,z,t,hour) = SATly - FC;
                                end
                            end
                            Wperc(1,y,z,t,hour) =  SWlyex(1,y,z,t,hour)*(1-exp(-1/TT));
                        end
                        %Deeper Layer Water
                    elseif x > 1
                        if hour == 1
                            SWly(x,y,z,t,hour) = Wperc(x-1,y,z,t-1,24) - Wperc(x,y,z,t-1,24) + SWly(x,y,z,t-1,24) - wuptake(x,y,z,t-1, 24);
                            if SWly(x,y,z,t,hour)>SATly
                                SWly(x,y,z,t,hour) = SATly;
                            end
                            if SWly(x,y,z,t,hour) < WP
                                SWly(x,y,z,t,hour) = WP;
                            else
                            end
                            if SWly(x,y,z,t,hour) <= FC
                                SWlyex(x,y,z,t,hour) = 0;
                            else
                                SWlyex(x,y,z,t,hour) = SWly(x,y,z,t,hour) - FC;
                            end
                            Wperc(x,y,z,t,hour) =  SWlyex(x,y,z,t,hour)*(1-exp(-1/TT));
                        else
                            SWly(x,y,z,t,hour) = Wperc(x-1,y,z,t,hour-1) - Wperc(x,y,z,t,hour-1) + SWly(x,y,z,t,hour-1) - wuptake(x,y,z,t, hour -1) ;
                            if SWly(x,y,z,t,hour)>SATly
                                SWly(x,y,z,t,hour) = SATly;
                            end
                            if SWly(x,y,z,t,hour) <= WP
                                SWly(x,y,z,t,hour)= WP;
                            else
                            end
                            if SWly(x,y,z,t,hour) <= FC
                                SWlyex(x,y,z,t,hour) = 0;
                            else
                                SWlyex(x,y,z,t,hour) = SWly(x,y,z,t,hour) - FC;
                            end
                            Wperc(x,y,z,t,hour) =  SWlyex(x,y,z,t,hour)*(1-exp(-1/TT));
                        end
                    end
                        if hour > 1
                                S(x,y,z,t,hour) = (SWly(x,y,z,t,hour) - WP)/(SATly - WP);
                                if S(x,y,z,t,hour) < 0
                                    S(x,y,z,t,hour) = 0;
                                elseif S(x,y,z,t,hour) > 1
                                    S(x,y,z,t,hour) = 1;
                                end
                                h(x,y,z,t,hour) = ((((S(x,y,z,t,hour)).^(-1/M)-1).^(1/N))/a); %ml/d^2*cm
                                psw(x,y,z,t,hour) = h(x,y,z,t,hour);
                            if map00(x,y,z,t) ~= 0
                                wuptake(x,y,z,t,hour) = -Kr*(psw(x,y,z,t,hour)-prw)*((t^0.6)*mapsurf00(x,y,z,t));
                                if wuptake(x,y,z,t,hour) < 0
                                    wuptake(x,y,z,t,hour) = 0;
                                elseif wuptake(x,y,z,t,hour) >= SWly(x,y,z,t,hour)
                                    
                                    wuptake(x,y,z,t,hour) = SWly(x,y,z,t,hour) -WP;
                                end
                            else
                                A(x,y,z,t,hour) = 0;
                                wuptake(x,y,z,t,hour) = 0;
                            end
                        else
                                S(x,y,z,t,hour) = (SWly(x,y,z,t,hour) - WP)/(SATly - WP);
                                if S(x,y,z,t,hour) < 0
                                    S(x,y,z,t,hour) = 0;
                                 elseif S(x,y,z,t,hour) > 1
                                    S(x,y,z,t,hour) = 1;
                                end
                                h(x,y,z,t,hour) = ((((S(x,y,z,t,hour)).^(-1/M)-1).^(1/N))/a); %cm
                                psw(x,y,z,t,hour) =  h(x,y,z,t,hour);
                                if psw(x,y,z,t,hour) < 0
                                    psw(x,y,z,t,hour) = 0;
                                end
                            if map00(x,y,z,t) ~= 0
                                wuptake(x,y,z,t,hour) = -Kr*(psw(x,y,z,t,hour)-prw)*((t^0.6)*mapsurf00(x,y,z,t));
                                if wuptake(x,y,z,t,hour) < 0
                                    wuptake(x,y,z,t,hour)= 0;
                                elseif wuptake(x,y,z,t,hour) >= SWly(x,y,z,t,hour)
                                    
                                    wuptake(x,y,z,t,hour) = SWly(x,y,z,t,hour) - WP;
                                end
                            else
                                wuptake(x,y,z,t,hour) = 0;
                            end
                        end
                    % SWAT NITRATE
                    if x == 1
                        if hour == 1
                            if before == 0 && t == apptime
                                NitMob(1,y,z,t,hour) = nitratel(1,y,z,t-1,24)*(1-exp((-Wperc(1,y,z,t-1,24))/(SATly)));
                                nitratel(1,y,z,t,hour) = nitratel(1,y,z,t-1,24) - NitMob(1,y,z,t-1,24)+application;
                                if isnan(nitratel(1,y,z,t,hour))
                                    disp("stop1, error: N concentration is NaN value due to error in soil moisture value")
                                    break
                                end
                            else
                                NitMob(1,y,z,t,hour) = nitratel(1,y,z,t-1,24)*(1-exp((-Wperc(1,y,z,t-1,24))/(SATly)));
                                nitratel(1,y,z,t,hour) = nitratel(1,y,z,t-1,24) - NitMob(1,y,z,t-1,24);
                                if isnan(nitratel(1,y,z,t,hour))
                                    disp("stop2, error: N concentration is NaN value due to error in soil moisture value")
                                    break
                                end
                            end
                        else
                            NitMob(1,y,z,t,hour) = nitratel(1,y,z,t,hour-1)*(1-exp((-Wperc(1,y,z,t,hour-1))/(SATly)));
                            nitratel(1,y,z,t,hour) = nitratel(1,y,z,t,hour-1) - NitMob(1,y,z,t,hour-1);
                            if isnan(nitratel(1,y,z,t,hour))
                                disp("stop3, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                        end
                    elseif x > 1
                        if hour == 1
                            NitMob(x,y,z,t,hour) = nitratel(x,y,z,t-1,24)*(1-exp((-Wperc(x,y,z,t-1,24))/(SATly)));
                            nitratel(x,y,z,t,hour) = nitratel(x,y,z,t-1,24) - NitMob(x,y,z,t-1,24) + NitMob(x-1,y,z,t-1,24);
                            if isnan(nitratel(1,y,z,t,hour))
                                disp("stop4, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                        else
                            NitMob(x,y,z,t,hour) = nitratel(x,y,z,t,hour-1)*(1-exp((-Wperc(x,y,z,t,hour-1))/(SATly)));
                            nitratel(x,y,z,t,hour) = nitratel(x,y,z,t,hour-1) - NitMob(x,y,z,t,hour-1) + NitMob(x-1,y,z,t,hour-1);
                            if isnan(nitratel(1,y,z,t,hour))
                                disp("stop5, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                        end
                    end
                    if before == 1 && t<apptime
                        uptake(x,y,z,t,hour) = 0;
                    elseif before == 1 && t>apptime
                        if round(map00(x,y,z,t)) == 0
                            uptake(x,y,z,t, hour) = 0;
                            continue
                        elseif round(map00(x,y,z,t)) == 1
                            age = mapage00(x,y,z,t);
                            agefact(x,y,z,t) = 1-((t - age)/120);
                            Imax = (30.03 * 86400 * 10^-12 * 62.05*1000); % mg/cm^2*d
                            Km = (6.72*0.001*62.05*10^-6*1000); %mg/cm^3
                            C = (nitratel(x,y,z,t,hour))/(SWly(x,y,z,t,hour)); %mg/cm^3
                            if isnan(C)
                                disp("stop6, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                            if Cgoal > Cint
                                uptake(x,y,z,t,hour) = agefact(x,y,z,t).*devfact(t).*((t^0.6)*(mapsurf00(x,y,z,t))*(Imax * (C)*(Cgoal-Cint))./(Km + (C)));
                            elseif Cgoal < Cint
                                uptake(x,y,z,t,hour) = 0;
                                %mg/cm^2*hour
                            end
                            if nitratel(x,y,z,t,hour) <= uptake(x,y,z,hour)
                                uptake(x,y,z,t,hour) = nitratel(x,y,z,t,hour);
                            end
                            if uptake(x,y,z,t,hour) < 0
                                uptake(x,y,z,t,hour) = 0;
                            end
                        elseif round(map00(x,y,z,t)) == 2
                            age = mapage00(x,y,z,t);
                            agefact(x,y,z,t) = 1-((t - age)/120);
                            Imax = (35.81 * 86400 * 10^-12 * 62.05*1000); % mg/cm^2*d
                            Km = 17.25*0.001*62.05*10^-6 * 1000 ; %mg/cm^3
                            C = (nitratel(x,y,z,t,hour))/(SWly(x,y,z,t,hour)); %mg/cm^3
                            if isnan(C)
                                disp("stop7, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                            if Cgoal > Cint
                                uptake(x,y,z,t,hour) = agefact(x,y,z,t).*devfact(t).*((t^0.6)*(mapsurf00(x,y,z,t))*(Imax * (C)*(Cgoal-Cint))./(Km + (C)));
                            elseif Cgoal < Cint
                                uptake(x,y,z,t) = 0;
                                %mg/cm^2*d
                            end
                            if nitratel(x,y,z,t,hour) <= uptake(x,y,z,t,hour)
                                uptake(x,y,z,t,hour) = nitratel(x,y,z,t,hour);
                            end
                            if uptake(x,y,z,t,hour) < 0
                                uptake(x,y,z,t,hour) = 0;
                            end
                        end
                    elseif before == 0
                        if round(map00(x,y,z,t)) == 1
                            age = mapage00(x,y,z,t);
                            agefact(x,y,z,t) = 1-((t - age)/120);
                            Imax = (30.03 * 86400 * 10^-12 * 62.05*1000); % mg/cm^2*d
                            Km = (6.72*0.001*62.05*10^-6*1000); %mg/cm^3
                            C = (nitratel(x,y,z,t,hour))/(SWly(x,y,z,t,hour)); %mg/cm^3
                            if isnan(C)
                                disp("stop8, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                            if Cgoal > Cint
                                uptake(x,y,z,t,hour) = agefact(x,y,z,t).*devfact(t).*((t^0.6)*(mapsurf00(x,y,z,t))*(Imax * (C)*(Cgoal-Cint))./(Km + (C)));
                            elseif Cgoal < Cint
                                uptake(x,y,z,t,hour) = 0;
                                %mg/cm^2*hour
                            end
                            if nitratel(x,y,z,t,hour) <= uptake(x,y,z,t,hour)
                                uptake(x,y,z,t,hour) = nitratel(x,y,z,t,hour);
                            end
                            if uptake(x,y,z,t,hour) < 0
                                uptake(x,y,z,t,hour) = 0;
                            end
                        elseif round(map00(x,y,z,t)) == 2
                            age = mapage00(x,y,z,t);
                            agefact(x,y,z,t) = 1-((t - age)/120);
                            Imax = (35.81 * 86400 * 10^-12 * 62.05*1000); % mg/cm^2*h
                            Km = 17.25*0.001*62.05*10^-6 * 1000 ; %mg/cm^3
                            C = (nitratel(x,y,z,t,hour))/(SWly(x,y,z,t,hour)); %mg/cm^3
                            if isnan(C)
                                disp("stop9, error: N concentration is NaN value due to error in soil moisture value")
                                break
                            end
                            if Cgoal > Cint
                                uptake(x,y,z,t,hour) = agefact(x,y,z,t).*devfact(t).*((t^0.6)*(mapsurf00(x,y,z,t))*(Imax * (C)*(Cgoal-Cint))./(Km + (C)));
                            elseif Cgoal < Cint
                                uptake(x,y,z,t,hour) = 0;
                                %mg/cm^2*d
                            end
                        end
                        if nitratel(x,y,z,t,hour) <= uptake(x,y,z,t,hour)
                            uptake(x,y,z,t,hour) = nitratel(x,y,z,t,hour);
                            nitratel(x,y,z,t,hour) = 0;
                        end
                        if uptake(x,y,z,t,hour) < 0
                            uptake(x,y,z,t,hour) = 0;
                        end
                    end
                end
            end
        end
        nitratel(:,:,:,t,hour) = nitratel(:,:,:,t,hour) - uptake(:,:,:,t,hour);
    end
    
end


%%
total1 = zeros(1,120);
totaluptake = zeros(1,120);
summeduptake = zeros(1,120);
for t = 1:simlength
if FC >SATly
    continue
end
nitrate1(:,:,t) = sum(nitratel(:,:,:,t),3);
SWly1(:,:,t) = sum(SWly(:,:,:,t),3);
conc1(:,:,t) = nitrate1(:,:,t)./SWly1(:,:,t);
nitratedepther(:,t) = sum(nitrate1(:,:,t),2);
upt1(:,:,t) = sum(uptake(:,:,:,t),3);
if t >= apptime
    nitratedepth(t) = find(nitratedepther(:,t) == max(nitratedepther(:,t)),1);
else
    nitratedepth(t) = 0;
end
total1(:,t) = sum(nitratel(:,:,:,t,:),'all');
totaluptake(t) = sum(uptake(:,:,:,t,:),'all');
nittotal(:,t) = sum(nitratel(:,:,:,t,24),'all') ;
waterup = sum(wuptake,'all');
summeduptake(:,t) = sum(totaluptake(:,1:t));
SWlysum(:,t) = sum(SWly(2:49,:,:,t),'all');
end
leached = zeros(1,simlength);
if apptime == 0
    for t = 1:simlength
        leached(t) = (mass*17*23)-summeduptake(:,t)-nittotal(:,t);
    end
else
    for t = apptime:simlength
        leached(t) = (mass*17*23)-summeduptake(:,t)-nittotal(:,t);
    end
end
if FC > SATly
    NUE = NaN;
    lost = NaN;
else
NUE = summeduptake(simlength)/(mass*(r1-1)*(r2-1));
disp(NUE)
leftinsystem = sum(nitrate1(:,:,simlength),'all')/total1(1+apptime);
moisture = mean(SWlysum)/((r1-1)*(r2-1)*(n-1));
for t = 2:simlength
    if round(nittotal(:,t),2) < round(nittotal(:,t-1),2)
        loss = t;
        break
    end
end
end
if NUE == 0 
    RLDnit = 0;
elseif max(nitratedepth) >= 5
    RLDnit = find(nitratedepth > 5,1); 
else
    RLDnit = 0;
end
%%
for t = 1:simlength
   a(t) = sum(wuptake(:,:,:,t,:), 'all');
end

%%
for i = 25:simlength*24
    day = ceil(i/24);
    Ho = (i-(24*(day-1)));
    nitrate2(:,:,i) = sum(nitratel(:,:,:,day,Ho),3);
    SWly2(:,:,i) = sum(SWly(:,:,:,day,Ho),3);
    upt2(:,:,i) = sum(uptake(:,:,:,day,Ho),3);
    mapsurf2(:,:,i) = sum(mapsurf00(:,:,:,day),3);
    nittotal2(:,i) = sum(nitratel(:,:,:,day, Ho),'all') + sum(uptake(:,:,:,day,Ho),'all');
    for a = 1:49
        nitrate3(a,i) = mean(nitrate2(a,:,i));
        SWly3(a,i) = mean(SWly2(a,:,i));
        upt3(a,i) = sum(nonzeros(upt2(a,:,i)));
        mapsurf3(a,i) = sum(nonzeros(mapsurf2(a,:,i)));
    end
end
for t = 25:simlength*24
    if round(nittotal2(:,t)) < round(nittotal2(:,t-1))
        loss = t;
        break
    end
end

for i = 2:simlength
        I = 24*(i-1);
        upt4(:,:,i) = sum(upt2(:,:,I:(I+24)),3);
end
leach = leached(:,simlength);
leachtime = find(leached >= 0.5*(leach),1);
leachrate = 50/(leachtime-apptime);

