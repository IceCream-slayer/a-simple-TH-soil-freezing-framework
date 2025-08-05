% Solves a TH coupling problem using the 
% finite difference method. 
% use two equation system, consider heat of solid and heat of fluid
%consider frost heave will grows even the energy flux larger than the
%needed energy to converting incoming water...
clear all

tic
fprintf('-------Calculation start---------\n')
fprintf('The default parameters replicate the results of the NS-5 experiment on Devon silt published by Konrad & Morgenstern (1980) \n')
fprintf(['Total simulated freezing duration is 50 h  \n' ...
    'The simulation will take several hours \n'])
fprintf(['For a quick look at the results, we can set the default time step to 0.01 s, elements1 to 40. \n' ...
    'Then, the calculation time will around 2 hours. However, the accuracy of this calculation is relatively a bit lower compared to the default setting \n'])



%% Input parameters %%

%% default maximum time step [s]
default = 0.0025; % %% can be adjusted if choose different elements number per cm, For example in the replications of Devon silts:elements=20 - 0.02s, elements=40 - 0.01s, elements=80 - 0.0025s

%% switch on or off for the convection during calculation
Convection_on = 1;  % 1 is on, 0 is off

%% define the space interval %%
column_length= 10; % the length of the soil column [cm]
elements1=80; % elements number per cm

%% Density, gravity accleration, inital porosity %%
RhoW = 998.2;    % reference water density at 20C
Rhoio = 917.0;   % reference ice density at 0C [kg/m3]
RhoS = 1700.0;    % solid particles density [kg/m3]
gAccel = 9.806;  % earth gravity acceleration [m/s2]
n_ini = 0.38;    %initial porosity

%% hydraulic conductivity %%
Tref = 273.16 + 0.0; % reference temperature
mu = (243.18 * 1e-7)*10^(247.8/(Tref - 140)); %Liquid dynamic viscosity [pa.s]
Ksat_1 = 1E-9;   %hydraulic conductivity m/s
%%-- if kref is given, please use the equations below --%%
%kref = 1.5e-16 % referenced intrinsic permeability [m2]
%Ksat_1 = (gAccel * RhoW * kref)/mu  %Hydraulic conductivitiy (permeability) at full saturation

%% thermal properties %%
Cs = 900;                             %heat capacity of solid component [J/kg/K]
Cw = 4180;                            %heat capacity of water component [J/kg/K]
Ci = 2108;                            %heat capacity of ice component   [J/kg/K]
lemda_soil_grain = 3;                 % thermal conductivity of soil grain material [W/m/K]
lemda_water = 0.56;                   % thermal conductivity of water at 0°C [W/m/K]
lemda_ice = 2.18;                     % thermal conductivity of ice at 0°C [W/m/K]
Latentiw = 3.34*10^5;                 % latent heat of water freezing [J/kg]

%% Parameters for heat transfer coefficient %%
dp = 0.03*0.001; % soil particle diameter [m]

%% Time settings %%
TotalTime=3600 *50*1; %total simulation time in seconds
Time_step_checked = 60*10; %Time interval for the data storing [s]

%% Capillary bundle model %%
bundle_number=5;
r1_orignal=1.42E-6;r2_orignal=2.13E-07;r3_orignal=4.6E-8;r4_orignal=1E-8;r5_orignal=4E-9;  %% pore radius [m]

%% initial Temperature and boundary condition %%
FreezingPoint = 273.16;
T_top = -6.2; 
T_bottom = 1.1;
T_body = 1.1;

%% set the initial state %%
%-----length of each grid [m]----------
L=column_length/100;  %length of the soil column in [m]
NumberDivisions = column_length*elements1; %number of divisions
NumberPoints = NumberDivisions+1;%number of points
delx =  zeros (1,NumberPoints-1);
delx(1,1:elements1*column_length)=0.01/elements1;
delx_ini = delx;
delx2 =delx(1)*delx(1);

%thermal properties-------------------
Ch_water = Cw*RhoW;
Ch_ice = Ci * Rhoio;
Ch_soil = Cs*RhoS;
Ch_smallest = min(Ch_soil,Ch_water);
lemda_biggest = max(lemda_soil_grain, lemda_ice);

%%--------calculate Pr,Re,Nu and h_ht---------
vs = Ksat_1; % the velocity of the fluid relate to the particle,setting velocity [m/s]
Pr = mu*Cw/lemda_water;
Re = (RhoW*dp*vs/n_ini)/mu;
Nu_fx = 2.4E-5 + 285.6*Pr^(1/3)*Re^2.7;
h_ht = Nu_fx*lemda_water/dp;

%--------Grid settings for plot--------
x(1) = 0.0;  %first grid point
for i=2:NumberPoints 
x(i) = x(i-1)+delx(1,i-1);  %grid points - for plotting
end 

%----calculate total time step------------
TotalSteps=ceil((TotalTime)/Time_step_checked);
%%--- Capillary bundle model ---%%
temporary_size = 2;
bundle_radi_orignal=[r1_orignal,r2_orignal,r3_orignal,r4_orignal,r5_orignal];
radius_1 =zeros (temporary_size,NumberPoints) + r1_orignal;
radius_2 =zeros (temporary_size,NumberPoints) + r2_orignal;
radius_3 =zeros (temporary_size,NumberPoints) + r3_orignal;
radius_4 =zeros (temporary_size,NumberPoints) + r4_orignal;
radius_5 =zeros (temporary_size,NumberPoints) + r5_orignal;
bundle_radi=zeros (NumberPoints,5);
bundle_tubes=zeros (NumberPoints,5);
SSA_bundle = zeros(NumberPoints,5);
Si_bundle =zeros (NumberPoints,5);
Volume_all_bundles = n_ini;
L_t = 1; % length of capillary tube. %%this value doesn't affect the calculation of specific surface area...
for i=1:NumberPoints
    bundle_radi(i,:)=[radius_1(1,i); radius_2(1,i); radius_3(1,i); radius_4(1,i); radius_5(1,i)];
    for j =1:bundle_number
        bundle_tubes(i,j) = n_ini*0.2/(3.14*bundle_radi(i,j)*bundle_radi(i,j)*L_t); % number of tubes in each bundle
    end
    for j =1:bundle_number
        SSA_bundle(i,j) =2*3.14*n_ini*bundle_tubes(i,j)*bundle_radi(i,j)*L_t/Volume_all_bundles;
    end
end
Number_of_capillarise=bundle_tubes(1,:);
h_heat_transfer = h_ht*sum(SSA_bundle(1,:));

SSA_bundle_orignal = SSA_bundle;
%%%--weighting factor for the contribution of hydraulic conductivity--%%%
w =zeros(1,5);
Q_wT_total = sum(bundle_radi_orignal.^4.*Number_of_capillarise);
for i =1:5
    w(i)=bundle_radi_orignal(i).^4*Number_of_capillarise(i)/Q_wT_total;
end

%%--------------------------------nodal variable storaged---------------------------
n =zeros (TotalSteps+1,NumberPoints) + n_ini;
n_fluid = zeros (TotalSteps+1,NumberPoints) + n_ini;       %porosity
Si = zeros (TotalSteps+1,NumberPoints);            %Ice saturation
%Sl = zeros (TotalSteps+1,NumberPoints) + 1;        %Water saturation
lemda_solid = zeros(TotalSteps+1,NumberPoints) +  lemda_soil_grain;           % thermal conductivity of the material
EnergyFlux_fluid = zeros (TotalSteps+1,NumberPoints);
EnergyFlux_solid = zeros (TotalSteps+1,NumberPoints);
EnergyFluxrate_fluid = zeros (TotalSteps+1,NumberPoints);
EnergyFluxrate_solid = zeros (TotalSteps+1,NumberPoints);
H_f = zeros (TotalSteps+1,NumberPoints);
C_fluid = zeros (TotalSteps+1,NumberPoints) + Cw;
Rho_Fluid= zeros (TotalSteps+1,NumberPoints) + RhoW;
n_ice = zeros(TotalSteps+1,NumberPoints) + 0;
convection_1 = zeros(TotalSteps+1,NumberPoints);

%heat transfer 
Heat_transfer_storaged =  zeros(TotalSteps+1,NumberPoints);

%related to water transfer equations
WaterHead = zeros (TotalSteps+1,NumberPoints);
n_total_water = n_fluid;                                %porosity computed with increased water volume.
WaterFlux = zeros (TotalSteps+1,NumberPoints);
WaterFluxRate = zeros (TotalSteps+1,NumberPoints);
Ksat = zeros (TotalSteps+1,NumberPoints)+Ksat_1;
n_WaterIntake = zeros(TotalSteps+1,NumberPoints);
%%--------------------------------nodal variable temporary---------------------------
temporary_size = 2;
n_temporary = zeros (temporary_size,NumberPoints) + n_ini;
n_fluid_temporary = zeros (temporary_size,NumberPoints) + n_ini;       %porosity
Si_temporary = zeros (temporary_size,NumberPoints);            %Ice saturation
Si_conduction_temporary = zeros (temporary_size,NumberPoints); 
Si_convection_temporary = zeros (temporary_size,NumberPoints); 
Si_local_thermal_transfer_temporary = zeros (temporary_size,NumberPoints);
lemda_solid_temporary = zeros(temporary_size,NumberPoints) +  lemda_soil_grain;           % thermal conductivity of the material
EnergyFlux_fluid_temporary = zeros (temporary_size,NumberPoints);
EnergyFlux_solid_temporary = zeros (temporary_size,NumberPoints);
EnergyFluxrate_fluid_temporary = zeros (temporary_size,NumberPoints);
EnergyFluxrate_solid_temporary = zeros (temporary_size,NumberPoints);
Heat_fluid_temporary = zeros(temporary_size,NumberPoints);
Heat_solid_temporary = zeros(temporary_size,NumberPoints);
H_f_temporary = zeros (temporary_size,NumberPoints);
dH_f_local_thermal_transfer = zeros (temporary_size,NumberPoints);
dH_f_convection = zeros (temporary_size,NumberPoints);
dH_f_conduction = zeros (temporary_size,NumberPoints);
n_ice_temporary = zeros(temporary_size,NumberPoints) + 0;                       % factors I used to avoid numerical error when porosity closed to zero
convection_temporary = zeros(temporary_size,NumberPoints);

%heat transfer 
Heat_transfer_temporary = zeros(temporary_size,NumberPoints);

%related to water transfer equations
WaterHead_temporary = zeros (temporary_size,NumberPoints);
n_WaterIntake_temporary = zeros (temporary_size,NumberPoints);
WaterFlux_temporary = zeros (temporary_size,NumberPoints);
WaterFluxRate_temporary = zeros (temporary_size,NumberPoints);
Ksat_temporary = zeros (temporary_size,NumberPoints)+Ksat_1;
n_WaterIntake1 = 0;
%%-------------------boundary and inital conditions for heat equation, calculate Enthalpy at inital time step----------
Heat_water =zeros (TotalSteps+1,NumberPoints); %defining the matrix containing Temperature
Heat_water(:,1) = 273.16 + T_top;     %Initial Temperature at top layer [°C]
Heat_water(:, NumberPoints) =273.16 + T_bottom;       %Initial Temperature at bottom layer [°C]
for i=2:NumberPoints-1
    Heat_water(1, i) = 273.16 + T_body;
    H_f(1,i) = Heat_water(1, i)*(RhoW*Cw);
end

if Heat_water(1:1) > FreezingPoint
    H_f(:,1) = Heat_water(1:1)*(RhoW*Cw);
else
    H_f(:,1) = Heat_water(1:1)*(Rhoio*Ci);
end

Heat_soil = Heat_water;
Heat_solid_temporary(1:2,:) = Heat_soil(1:2,:);
Heat_fluid_temporary(1:2,:) = Heat_water(1:2,:);
H_f_temporary(1:2,:) = H_f(1:2,:);
%%-----------------------------------------------------------------------------
% add in the initial condition of ice saturation
for i=1:NumberPoints
    if Heat_water(1, i) <= FreezingPoint
        Si(:, i) = 1;
    elseif Heat_water(1, i) >FreezingPoint
           Si(:, i) = 0;
    end
end

n_ice(:,1) = Si(1,1)*n_ini;
n_ice_temporary(1:2,:) = n_ice(1:2,:);
Si_temporary(1:2,:) = Si(1:2,:);
% assign porosity and lemda at the boundary
 n_fluid(:,1) = (1-Si(1,1))*n_ini;
 if n_fluid(1,1) == 0
     n_fluid(:,1) = 1E-05;
 end
 n_fluid_temporary(:,1) = n_fluid(1,1);
 lemda_solid(:,1) = lemda_soil_grain;
 lemda_solid_temporary(1:2,:) = lemda_soil_grain;


%%-----------------------------water head and Pressure-----------------------
SuctionHead = 2*0.032/(r1_orignal*gAccel*RhoW);
SampleHeight = abs(L);
WaterHead(:,1) = 0; %upper boundary condition
WaterHead_temporary(1:2,:) = WaterHead(1:2,:);

%----------------calculation start from here-------------------
A=0;
delx_fixed =delx_ini(1);
Cumulative_step = 2;
Cumulative_TimeStep = 0;
Total_step_Time = 0;
i = 1;
checkpoint = 1; % check wether we have partily frozen soil

%% Define Time step %%
safety = 2; %safety factor for the stability requirement
Time_step1 = delx2/(safety*(lemda_biggest/Ch_smallest)); %stability requirement for Thermal - safety is safety factor,Ch is Soil heat capacity(Rhos * c)
Time_step2 = ((n_ini*Ch_smallest/h_heat_transfer));
Max_Time_step = default;
array =[Time_step1,Time_step2,Max_Time_step];
dt = round(min(array),4); %Time step

fprintf('Calculating Progress: %d / %d  \n', (Cumulative_step-1),TotalSteps)

while 1
    Total_step_Time = Total_step_Time + dt;
    Cumulative_TimeStep = Cumulative_TimeStep + dt;
    if Total_step_Time >= TotalTime
        break
    end
    for j = 2:(NumberPoints-1)
        if Si_bundle(j,1) < 0
            Ksat_temporary(i,j) = Ksat_1;
        else
            Ksat_temporary(i,j) =  Ksat_1*sum(w.*((n_ini*(1-Si_bundle(j,:)*5)).^3./(1-(n_ini*(1-Si_bundle(j,:)*5))).^2) * ((1-n_ini)^2/n_ini^3));
        end
        if Si_temporary(i,j) == 1
             Ksat_temporary(i,j) = 0;
        end

%%%-----water head gradient-----%%% 
        if  Si_temporary(i,j) < 1 &&  Si_temporary(i,j-1) == 1  
            WaterHead_temporary(i,j) = min(2*0.032/(bundle_radi(j,1)*gAccel*RhoW),10); % maximum capillary rise value should not exceed 10m
            if  WaterHead_temporary(i,j+1) > WaterHead_temporary(i,j)
                for l=j+1:elements1*column_length
                    WaterHead_temporary(i,l)= WaterHead_temporary(i,j)-WaterHead_temporary(i,j)*delx_ini(1)/(delx_ini(1)*(elements1*column_length+1-j))*(l-j);
                end
                WaterHead_temporary(i,elements1*column_length+1)=0;
            end
            if  WaterHead_temporary(i,j) > WaterHead_temporary(i+1,j)
                for l=j+1:elements1*column_length
                    WaterHead_temporary(i,l)= WaterHead_temporary(i,j)-WaterHead_temporary(i,j)*delx_ini(1)/(delx_ini(1)*(elements1*column_length+1-j))*(l-j);
                end
                WaterHead_temporary(i,elements1*column_length+1)=0;
            end
        elseif Si_temporary(i,j) == 1 &&  Si_temporary(i,j-1) == 1
            WaterHead_temporary(i,j) = 0;
        end


%%%-----------------------------------------------
total_SSA = sum(SSA_bundle(j,:));
Heat_transfer_temporary(i,j) = total_SSA*h_ht;


%%% calculate energy and water flux according to the finite difference method------------------------------------

        Value = [0,0,0];
        Value(1)=(2/(delx(j)+delx(j-1)))*((WaterHead_temporary(i,j+1)-WaterHead_temporary(i,j))/delx(j)-(WaterHead_temporary(i,j)-WaterHead_temporary(i,j-1))/delx(j-1));
        
        if Si_temporary(i,j) == 1 &&  Si_temporary(i,j-1) == 1
            Value(1)=0;
        elseif Si_temporary(i,j) < 1 &&  Si_temporary(i,j-1) == 1 
            Value(1)=0;
        end

        Value(2)=(2/(delx(j)+delx(j-1)))*((Heat_fluid_temporary(i,j+1)-Heat_fluid_temporary(i,j))/delx(j)-(Heat_fluid_temporary(i,j)-Heat_fluid_temporary(i,j-1))/delx(j-1));
        Value(3)=(2/(delx(j)+delx(j-1)))*((Heat_solid_temporary(i,j+1)-Heat_solid_temporary(i,j))/delx(j)-(Heat_solid_temporary(i,j)-Heat_solid_temporary(i,j-1))/delx(j-1));


        EnergyFlux_solid_temporary(i,j-1) = -lemda_solid_temporary(i,j)*(Heat_solid_temporary(i,j) - Heat_solid_temporary(i,j-1))/abs(delx(j-1)); % Energy Flux [W/m2 or J/(m2 * s)] 
        EnergyFlux_solid_temporary(i,j) = -lemda_solid_temporary(i,j)*(Heat_solid_temporary(i,j+1) - Heat_solid_temporary(i,j))/abs(delx(j)); % Energy Flux [W/m2 or J/(m2 * s)] 
        EnergyFluxrate_solid_temporary(i,j) = -lemda_solid_temporary(i,j) * Value(3);  % Energy flux rate [W/m3]

        lemda_pc = lemda_water*(1-Si_temporary(i,j))+lemda_ice*Si_temporary(i,j);
        EnergyFlux_fluid_temporary(i,j-1) = -lemda_pc*(Heat_fluid_temporary(i,j) - Heat_fluid_temporary(i,j-1))/abs(delx(j-1)); % Energy Flux [W/m2 or J/(m2 * s)] 
        EnergyFlux_fluid_temporary(i,j) = -lemda_pc*(Heat_fluid_temporary(i,j+1) - Heat_fluid_temporary(i,j))/abs(delx(j)); % Energy Flux [W/m2 or J/(m2 * s)] 
        EnergyFluxrate_fluid_temporary(i,j) = -lemda_pc * Value(2);  % Energy flux rate [W/m3]


        WaterFluxRate_temporary(i,j) = (WaterFlux_temporary(i+1,j+1) - WaterFlux_temporary(i+1,j))/abs(delx(j)); %WaterFluxRate, actually this calculated the waterfluxrate from previous steps
        WaterFlux_temporary(i,j) = -Ksat_temporary(i,j) *(WaterHead_temporary(i,j+1) - WaterHead_temporary(i,j))/abs(delx(j));  % darcy's law - water velocity, water Flux [m2*m/s]
       

        
%%% soil freezing and thawing %%% 
        heat_convection =Convection_on*RhoW*Cw*(WaterFlux_temporary(i,j)*(Heat_fluid_temporary(i,j+1) - Heat_fluid_temporary(i,j))/abs(delx(j)));  % W/m3
        convection_temporary(i,j)=heat_convection;

        if  H_f_temporary(i,j) > RhoW*Cw*FreezingPoint %Check Whether it is before a phase transition occurs

            Heat_solid_temporary(i+1,j) = Heat_solid_temporary(i,j) + dt*((lemda_solid_temporary(i,j)*Value(3))/Ch_soil + (Heat_transfer_temporary(i,j)*(Heat_fluid_temporary(i,j) - Heat_solid_temporary(i,j)))/((1-n_ini) * Ch_soil));
            H_f_temporary(i+1,j) = H_f_temporary(i,j) + (lemda_water*Value(2) + (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j)))*dt + heat_convection*dt;
            Heat_fluid_temporary(i+1,j) = H_f_temporary(i,j)/(Ch_water);
            Si_temporary(i+1,j) = 0;


        elseif H_f_temporary(i,j) <= RhoW*Cw*FreezingPoint && H_f_temporary(i,j) >= RhoW*Cw*FreezingPoint - RhoW*Latentiw %Check Whether it is in a phase transition state
           Energy_phase_change_waterintake = (EnergyFluxrate_fluid_temporary(i,j))/(Latentiw * RhoW);% how much intake_water can be converted into ice based on energy flux we had...
            if Si_temporary(i,j-1) == 1 && abs(EnergyFluxrate_fluid_temporary(i,j)) <= abs(WaterFluxRate_temporary(i,j-1)*Latentiw * RhoW) %Check Whether or not we're on the freezing front and there is not enough energy to convert all the incoming water into ice...

                 n_WaterIntake_temporary(i+1,j-1) =  n_WaterIntake_temporary(i,j-1) + dt*Energy_phase_change_waterintake; % we stored the increased volume due to water migration
                 Heat_solid_temporary(i+1,j) = Heat_solid_temporary(i,j) + dt*((lemda_solid_temporary(i,j)*Value(3))/Ch_soil + (Heat_transfer_temporary(i,j)*(Heat_fluid_temporary(i,j) - Heat_solid_temporary(i,j)))/((1-n_ini) * Ch_soil));
                 H_f_temporary(i+1,j) = H_f_temporary(i,j) + (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j))*dt + heat_convection/n_ini*dt;
                 Heat_fluid_temporary(i+1,j) = FreezingPoint;
                 dH_f_local_thermal_transfer(i+1,j) = (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j))*dt;
                 dH_f_convection(i+1,j) = convection_temporary(i,j)/n_ini*dt;
                 dH_f_conduction(i+1,j) = 0;
                 %delx variation
                 delx(j-1) = delx_ini(j-1)*(1+Si_temporary(i+1,j-1)*(RhoW/Rhoio-1))+n_WaterIntake_temporary(i+1,j-1)*delx_ini(j-1)*(RhoW/Rhoio);

            elseif Si_temporary(i,j-1) == 1                                                                                                                             %If we are on the freezing front and had enough energy flux to freeze downwards
                Heat_solid_temporary(i+1,j) = Heat_solid_temporary(i,j) + dt*((lemda_solid_temporary(i,j)*Value(3))/Ch_soil + (Heat_transfer_temporary(i,j)*(Heat_fluid_temporary(i,j) - Heat_solid_temporary(i,j)))/((1-n_ini) * Ch_soil));
                n_WaterIntake_temporary(i+1,j-1) =  n_WaterIntake_temporary(i,j-1) + dt*(WaterFlux_temporary(i,j))/delx(j-1);
                H_f_temporary(i+1,j) = H_f_temporary(i,j) + ((-EnergyFluxrate_fluid_temporary(i,j)+(Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j)) +heat_convection/n_ini+abs(WaterFluxRate_temporary(i,j-1)*Latentiw * RhoW)))*dt;
                Heat_fluid_temporary(i+1,j) = FreezingPoint;
                dH_f_local_thermal_transfer(i+1,j) = (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j))*dt;
                dH_f_convection(i+1,j) = convection_temporary(i,j)/n_ini*dt;
                dH_f_conduction(i+1,j) =(-EnergyFluxrate_fluid_temporary(i,j)+abs(WaterFluxRate_temporary(i,j-1)*Latentiw * RhoW))*dt;
                 %delx variation
                 delx(j-1) = delx_ini(j-1)*(1+Si_temporary(i+1,j-1)*(RhoW/Rhoio-1))+n_WaterIntake_temporary(i+1,j-1)*delx_ini(j-1)*(RhoW/Rhoio);
            
            else %not on the freezing front
                Heat_solid_temporary(i+1,j) = Heat_solid_temporary(i,j) + dt*((lemda_solid_temporary(i,j)*Value(3))/Ch_soil + (Heat_transfer_temporary(i,j)*(Heat_fluid_temporary(i,j) - Heat_solid_temporary(i,j)))/((1-n_ini) * Ch_soil));
                H_f_temporary(i+1,j) = H_f_temporary(i,j) + (-EnergyFluxrate_fluid_temporary(i,j) + (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j)))*dt + heat_convection/n_ini*dt;
                Heat_fluid_temporary(i+1,j) = FreezingPoint;
                dH_f_local_thermal_transfer(i+1,j) = (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j))*dt;
                dH_f_convection(i+1,j) = convection_temporary(i,j)/n_ini*dt;
                dH_f_conduction(i+1,j) =-EnergyFluxrate_fluid_temporary(i,j)*dt;

                 %delx variation
                 n_WaterIntake_temporary(i+1,j-1) =  n_WaterIntake_temporary(i,j-1);
                 delx(j-1) = delx_ini(j-1)*(1+Si_temporary(i+1,j-1)*(RhoW/Rhoio-1))+n_WaterIntake_temporary(i+1,j-1)*delx_ini(j-1)*(RhoW/Rhoio);
            end
           
           Si_local_thermal_transfer_temporary(i+1,j) =Si_local_thermal_transfer_temporary(i,j) + (-dH_f_local_thermal_transfer(i,j))/(Latentiw*RhoW);
           Si_convection_temporary(i+1,j) = Si_convection_temporary(i,j) + (-dH_f_convection(i,j))/(Latentiw*RhoW);
           Si_temporary(i+1,j) =  (RhoW*Cw*FreezingPoint-H_f_temporary(i,j))/(Latentiw*RhoW);
           Si_conduction_temporary(i+1,j) = Si_temporary(i+1,j)- Si_local_thermal_transfer_temporary(i+1,j)-Si_convection_temporary(i+1,j);
           bundles = bundle_number;
           for k=1:bundle_number
               if Si_bundle(j,6-k) >= 0.2
                     Si_bundle(j,6-k)= 0.2;
                     bundles = bundles-1;
               else
                   if -dH_f_conduction(i,j)/(Latentiw*RhoW) + (-dH_f_convection(i,j))/(Latentiw*RhoW) > 0
                       Si_bundle(j,6-k)=Si_bundle(j,6-k) + ( -dH_f_conduction(i,j)/(Latentiw*RhoW) + (-dH_f_convection(i,j))/(Latentiw*RhoW))/bundles + ((-dH_f_local_thermal_transfer(i,j))/(Latentiw*RhoW))*SSA_bundle(j,6-k)/total_SSA;
                   else
                       Si_bundle(j,6-k)=Si_bundle(j,6-k) +((-dH_f_local_thermal_transfer(i,j))/(Latentiw*RhoW)+(-dH_f_convection(i,j))/(Latentiw*RhoW) + -dH_f_conduction(i,j)/(Latentiw*RhoW))*SSA_bundle(j,6-k)/total_SSA;
                   end
                if Si_bundle(j,6-k) >= 0.2
                     Si_bundle(j,6-k)= 0.2;
                end
               end
               bundle_radi(j,6-k) = sqrt(1-Si_bundle(j,6-k)*5)*bundle_radi_orignal(1,6-k);
               if bundle_radi(j,6-k) < 1e-20
                   bundle_radi(j,6-k)=1e-20;
               end
               SSA_bundle(j,6-k) = SSA_bundle_orignal(j,6-k)*(1-Si_bundle(j,6-k)*5);  %update SSA due to pore radius change during freezing
           end

        else  % it is fully frozen
            Heat_transfer_temporary(i,j)=h_heat_transfer;
            Heat_solid_temporary(i+1,j) = Heat_solid_temporary(i,j) + dt*((lemda_soil_grain*Value(3))/Ch_soil + (Heat_transfer_temporary(i,j)*(Heat_fluid_temporary(i,j) - Heat_solid_temporary(i,j)))/((1-n_ini) * Ch_soil));
            H_f_temporary(i+1,j) = H_f_temporary(i,j) + (-EnergyFluxrate_fluid_temporary(i,j) + (Heat_transfer_temporary(i,j)/n_ini)*(Heat_solid_temporary(i,j) - Heat_fluid_temporary(i,j)))*dt + heat_convection*dt;
            Heat_fluid_temporary(i+1,j) = (H_f_temporary(i,j)-(RhoW*Cw*FreezingPoint - RhoW*Latentiw))/(Ch_ice) + 273.16;
            Si_temporary(i+1,j) = 1;
    end



  %%% %update soil volume fraction  %%% 
       n_fluid_temporary(i+1,j) = (1-Si_temporary(i+1,j))*n_ini;
       n_ice_temporary(i+1,j) = Si_temporary(i+1,j)*n_ini;
       if n_fluid_temporary(i+1,j) <= 1E-3
           n_fluid_temporary(i+1,j) = 1E-8;
       end
       soil_volume_fraction = (1-n_ini)/(1-n_fluid_temporary(i+1,j));


    end

%%% storage check, storing the temporary data into datasets %%%
 if Cumulative_TimeStep >= Time_step_checked
    Cumulative_TimeStep = Cumulative_TimeStep-Time_step_checked;
    n_ice(Cumulative_step,:) = n_ice_temporary(i+1,:);
    n_fluid(Cumulative_step,:) = n_fluid_temporary(i+1,:);
    Si(Cumulative_step,:) = Si_temporary(i+1,:);
    Heat_soil(Cumulative_step,:) = Heat_solid_temporary(i+1,:);
    Heat_water(Cumulative_step,:) = Heat_fluid_temporary(i+1,:); 
    lemda_solid(Cumulative_step,:) = lemda_solid_temporary(i+1,:);
    EnergyFlux_fluid(Cumulative_step,:) = EnergyFlux_fluid_temporary(i,:);
    EnergyFlux_solid(Cumulative_step,:) = EnergyFlux_solid_temporary(i,:);
    EnergyFluxrate_fluid(Cumulative_step,:) = EnergyFluxrate_fluid_temporary(i,:);
    EnergyFluxrate_solid(Cumulative_step,:) = EnergyFluxrate_solid_temporary(i,:);
    WaterHead(Cumulative_step,:) = WaterHead_temporary(i,:);
    n_WaterIntake(Cumulative_step,:) = n_WaterIntake_temporary(i+1,:);
    WaterFlux(Cumulative_step,:) = WaterFlux_temporary(i,:);
    WaterFluxRate(Cumulative_step,:) = WaterFluxRate_temporary(i,:);
    Ksat(Cumulative_step,:) = Ksat_temporary(i,:);
    H_f(Cumulative_step,:) = H_f_temporary(i,:);
    convection_1(Cumulative_step,:) = convection_temporary(i,:);
    Heat_transfer_storaged(Cumulative_step,:) = Heat_transfer_temporary(i,:);

    fprintf('Calculated freezing duration: %.1f h / %d h \n', (Cumulative_step-1)*Time_step_checked/3600,TotalTime/3600 )
    Cumulative_step = Cumulative_step + 1;
    fprintf('Calculating Progress: %d / %d  \n', (Cumulative_step-1),TotalSteps)
 end

%%% update the temporary data %%%
%--update data of next step into current step
n_fluid_temporary(i,:) = n_fluid_temporary(i+1,:);
n_ice_temporary(i,:) = n_ice_temporary(i+1,:);
Si_temporary(i,:) = Si_temporary(i+1,:);
Si_conduction_temporary(i,:) = Si_conduction_temporary(i+1,:); 
Si_convection_temporary(i,:) = Si_convection_temporary(i+1,:);
Si_local_thermal_transfer_temporary(i,:) = Si_local_thermal_transfer_temporary(i+1,:);
H_f_temporary(i,:) = H_f_temporary(i+1,:);
dH_f_local_thermal_transfer(i,:) = dH_f_local_thermal_transfer(i+1,:);
dH_f_convection(i,:)= dH_f_convection(i+1,:);
dH_f_conduction(i,:)= dH_f_conduction(i+1,:);
lemda_solid_temporary(i,:) = lemda_solid_temporary(i+1,:);
Heat_fluid_temporary(i,:) = Heat_fluid_temporary(i+1,:);
Heat_solid_temporary(i,:) = Heat_solid_temporary(i+1,:);

%--store the data from previous steps into i+1
EnergyFlux_fluid_temporary(i+1,:) = EnergyFlux_fluid_temporary(i,:);
EnergyFluxrate_fluid_temporary(i+1,:) = EnergyFluxrate_fluid_temporary(i,:);
EnergyFlux_solid_temporary(i+1,:) = EnergyFlux_solid_temporary(i,:);
EnergyFluxrate_solid_temporary(i+1,:) = EnergyFluxrate_solid_temporary(i,:);
WaterHead_temporary(i,:) = WaterHead_temporary(i+1,:);
n_WaterIntake_temporary(i,:) = n_WaterIntake_temporary(i+1,:);
WaterFlux_temporary(i+1,:) = WaterFlux_temporary(i,:);
WaterFluxRate_temporary(i+1,:) = WaterFluxRate_temporary(i,:);
Ksat_temporary(i+1,:) = Ksat_temporary(i,:);


end

%% Calculation done %%

%% Ploting %%
Heat_soil = Heat_soil-273.16;    %Changing the unit of temperature from kelvins to degrees Celsius
Heat_water = Heat_water-273.16;  %Changing the unit of temperature from kelvins to degrees Celsius
%%--------------
calculated_time = zeros(TotalSteps,1);
for i=1:TotalSteps+1
    calculated_time(i) = 0+(i-1)*Time_step_checked;
end

%%%%ploting total heave%%%
Heave_insitu = zeros(TotalSteps+1,1);
Heave_waterintake = zeros(TotalSteps+1,1);
Heave_total = zeros(TotalSteps+1,1);
Heave_time = zeros(TotalSteps+1,1);
for i=2:TotalSteps-1
    Sum_heave_WaterIntake = 0;
    Sum_heave_insitu = 0;
    for j =1:NumberPoints-1
        Sum_heave_WaterIntake = Sum_heave_WaterIntake + n_WaterIntake(i,j)*delx_ini(j)*(RhoW/Rhoio);
        Sum_heave_insitu = Sum_heave_insitu + n_ini*Si(i,j)*delx_ini(j)*(RhoW/Rhoio-1);
    end
    Heave_waterintake(i) = Sum_heave_WaterIntake*1000;  % volume expansion of waterintake in [mm]
    Heave_insitu(i)= Sum_heave_insitu*1000; % volume expansion of in situ water in [mm]
    Heave_total(i) = Heave_insitu(i) + Heave_waterintake(i);
    Heave_time(i) = i*Time_step_checked; 
end
plot(Heave_time(1:Cumulative_step-2)/3600,Heave_total(1:Cumulative_step-2),'r--')
xlabel('Heave time h') 
ylabel('Total heave mm')
data_heave_experimetal =[2.5, 3.8, 5.5, 7.5, 8.5, 9.4,10.9,11.5];
data__heave_hour_experimetal =[5,10,15,20,25,30,40,50];
hold on
plot(data__heave_hour_experimetal(1:7),data_heave_experimetal(1:7),'r*')
legend('total heave(simulated)','total heave(experimental data, Konrad 1980)')
fprintf('-------calculation done---------\n')
toc

%print('converge1', '-dpng', '-r1200') % Print high DPI(1200) figures in png 