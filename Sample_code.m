clc;
clear;
%% Load data for sample station 4033000
region = load('region4_Smax100.txt');
stationID = region(:,1); % station ID
drainArea = region(:,2); % station drainage area
elevation = region(:,3); % station elevation
Sest = region(:,4); % estimated soil moisture

% Select sample station
sample = 5;

ID = stationID(sample); 
DA = drainArea(sample); 
Elev = elevation(sample); 
b = Sest(sample);

file_pre = strcat(num2str(ID),'.pre');
file_pet = strcat(num2str(ID),'.pet');
file_q = strcat(num2str(ID),'.AMM');
file_tmax = strcat(num2str(ID),'.TMAX');
file_tmin = strcat(num2str(ID),'.TMIN');

pre_dat = load(file_pre); %Precipitation
pet_dat = load(file_pet); %PET
q_obs_dat = load(file_q); %Discharge
tmax_dat = load(file_tmax); %Tmax
tmin_dat = load(file_tmin); %Tmin

%% Data Preprocessing
nrow_pre = size(pre_dat,1); %Number of years
count = 0;

% Match data timeframes
for i = 1:nrow_pre
    year = pre_dat(i,2); %Year
    pet_coord = find(pet_dat(:,2) == year);
    q_coord = find(q_obs_dat(:,2) == year);
    tmax_coord = find(tmax_dat(:,2) == year);
    tmin_coord = find(tmin_dat(:,2) == year);
    if ~isempty(pet_coord) && ~isempty(q_coord) && ~isempty(tmax_coord) && ~isempty(tmin_coord)
        count = count + 1;
        P(count,:) = pre_dat(i,:);
        PET(count,:) = pet_dat(pet_coord,:);
        Q(count,:) = q_obs_dat(q_coord,:);
        TMAX(count,:) = tmax_dat(tmax_coord,:);
        TMIN(count,:) = tmin_dat(tmin_coord,:);
    end
end

%% Estimating snowfall based on temperature
m  = size(P,1); % Number of effective years
S = P(:,1:2); % Snowfall Matrix

for i = 1:m
    for j = 1:12
        if ((TMIN(i,j+2) + TMAX(i,j+2))/2)*0.1 > 0.5
            S(i,j+2) = 0;
        elseif ((TMIN(i,j+2) + TMAX(i,j+2))/2)*0.1 < -0.5
            S(i,j+2) = P(i,j+2);
        else
            S(i,j+2) = (1 - (((TMIN(i,j+2) + TMAX(i,j+2))/2)*0.1 + 0.5)) * P(i,j+2);
        end
    end
end

%% Organize data
count = 1;
max_count = 1;
ly = P(1,2);
last = P(1,2);
for i = 2 : size(P,1)
   ny = P(i,2);
   if ny == ly+1
       count = count + 1;
   else
       count = 0;
   end
   
   if count > max_count
       max_count = count;
       last = ny;
   end
   ly = ny;
end

coord = find(P(:,2) == last);
P = P(coord-max_count+1:coord,:);
PET = PET(coord-max_count+1:coord,:);
Q = Q(coord-max_count+1:coord,:);
S = S(coord-max_count+1:coord,:);
TMAX = TMAX(coord-max_count+1:coord,:);
TMIN = TMIN(coord-max_count+1:coord,:);

%% Convert matrix data to 1-d array
months = 12;
pre_data = P(:,3:months+2);
pet_data = PET(:,3:months+2);
q_data = Q(:,3:months+2);
s_data = S(:,3:months+2);
tmax_data = TMAX(:,3:months+2);
tmin_data = TMIN(:,3:months+2);

pre_col = pre_data(1,:)';
pet_col = pet_data(1,:)';
q_col = q_data(1,:)';
s_col = s_data(1,:)';
tmax_col = tmax_data(1,:)';
tmin_col = tmin_data(1,:)';

for i = 2 : size(pre_data,1)
    pre_col = vertcat(pre_col,pre_data(i,:)');
    pet_col = vertcat(pet_col,pet_data(i,:)');
    q_col = vertcat(q_col,q_data(i,:)');
    s_col = vertcat(s_col,s_data(i,:)');
    tmax_col = vertcat(tmax_col,tmax_data(i,:)');
    tmin_col = vertcat(tmin_col,tmin_data(i,:)');
end

%% Convert cfs to mm/month
q_col = q_col*(74.47/DA);

%% Account for days in month
 for i=0:(size(q_col,1)/12)-1
     q_col(1+12*i) = q_col(1+12*i)*(31/30); %OCT
     q_col(3+12*i) = q_col(3+12*i)*(31/30); %DEC
     q_col(4+12*i) = q_col(4+12*i)*(31/30); %JAN
     q_col(5+12*i) = q_col(5+12*i)*(28/30); %FEB
     q_col(6+12*i) = q_col(6+12*i)*(31/30); %MAR
     q_col(8+12*i) = q_col(8+12*i)*(31/30); %MAY
     q_col(10+12*i) = q_col(10+12*i)*(31/30); %JUL
     q_col(11+12*i) = q_col(11+12*i)*(31/30); %AUG
 end

%% Convert calendar year to water year
n = size(q_col,1);
q = q_col(13:n,1);
rf = pre_col(10:n-3,1) - s_col(10:n-3,1);
pet = pet_col(10:n-3,1);
sf = s_col(10:n-3,1);
tmax = tmax_col(10:n-3,1);
tmin = tmin_col(10:n-3,1);
nmonths = size(q,1);

% End of data preprocessing %

%% Optimization

% Replace 0 discharge to 0.001
[ii,~] = find(q==0);
q(ii,1) = 0.001;

% Initialize parameters
initial_beta(1:12,1) = 0.3;
initial_alpha1(1:12,1) = 0.4;
initial_alpha2(1:12,1) = 0.3;
initial_d(1:12,1) = 0.8;
initial_Smax(1:12,1) = b;
paramGuess = [initial_beta, initial_alpha1, initial_alpha2, initial_d, initial_Smax];

% Optimization settings
options = optimset('Algorithm','interior-point', 'MaxFunEvals', 8e6, 'MaxIter', 8e6,'Display','notify','DerivativeCheck','on');

% Set parameter boundaries
lowerbound_beta(1:12,1) = 0;
lowerbound_alpha1(1:12,1) = 0;
lowerbound_alpha2(1:12,1) = 0;
lowerbound_d(1:12,1) = 0;
lowerbound_Smax(1:12,1) = 0.2*b;
upperbound_beta(1:12,1) = 1;
upperbound_alpha1(1:12,1) = 1;
upperbound_alpha2(1:12,1) = 1;
upperbound_d(1:12,1) = 1;
upperbound_Smax(1:12,1) = 5*b;
    
lowerBound = [lowerbound_beta,lowerbound_alpha1,lowerbound_alpha2,lowerbound_d,lowerbound_Smax];
upperBound = [upperbound_beta,upperbound_alpha1,upperbound_alpha2,upperbound_d,upperbound_Smax];

% Minimize zhangModelError
[paramCalib] = fmincon(@(paramCalib)zhangModelError(paramCalib, q, rf, pet, sf, nmonths, tmax, tmin),...
    paramGuess,[],[],[],[],lowerBound,upperBound,[],options);

beta = paramCalib(:,1);
alpha1 = paramCalib(:,2);
alpha2 = paramCalib(:,3);
d = paramCalib(:,4);
s_max = paramCalib(:,5);

% Adjust beta values for no snow conditions
[qEst evap stor storFinal1 groundFinal1 snowFinal1 Snowpack] = zhang_model_snow(rf,pet,sf,nmonths,beta,alpha1,alpha2,...
    s_max,d,0,0,0,tmax,tmin);

    yr = nmonths/12;
      for k = 1:yr
       for j = 1:12
            monthly_snowpack(k,j) = Snowpack(12*(k-1)+j,1);
       end
      end
      
      for k = 1:yr
       for j = 1:12
            monthly_snowfall(k,j) = sf(12*(k-1)+j,1);
       end
      end
    annual_snowpack = mean(monthly_snowpack,1);
    [~,ii] = find(annual_snowpack < 0.1);
    beta(ii) = 1;
    
% Estimate Q with adjusted beta values
[qEst evap stor snowFinal2 storFinal2 groundFinal2 Snowpack] = zhang_model_snow(rf,pet,sf,nmonths,beta,alpha1,alpha2,...
    s_max,d,0,0,0, tmax, tmin); % Full-cycle spin-up

[qEst evap stor snowFinal3 storFinal3 groundFinal3 Snowpack] = zhang_model_snow(rf,pet,sf,nmonths,beta,alpha1,alpha2,...
    s_max,d,storFinal2,groundFinal2,snowFinal2, tmax, tmin);

betas = beta;
alpha1s = alpha1;
alpha2s = alpha2;
ds = d;
Smaxs = s_max;

%% Evaluate NSE
nashutcliffe(qEst,q)