function[m0,V_step_walk,M_prior]=set_Metropolis_hastings(MeasuredRainTS)

% Infer a0 by inversion of proportion of zero measurements
data=[];
nb_pts=0;
nb_zeros=0;
for i=1:length(MeasuredRainTS)
    for j=1:length(MeasuredRainTS(i).RainRate)
        if ~isnan(MeasuredRainTS(i).RainRate)
            nb_pts=nb_pts+1;
            if MeasuredRainTS(i).RainRate(j)>0
                data=[data;MeasuredRainTS(i).RainRate(j)];
            else
                nb_zeros=nb_zeros+1;
            end
        end
    end
end
prop_zeros=nb_zeros/nb_pts;

%initial values of transform fonction by quantile / quantile mapping
data_sort=sort(data);
M_anamorphose=zeros(10+1,2);
M_anamorphose(1,1)=0;
M_anamorphose(1,2)=norminv(prop_zeros);
for i=1:10
    ind_datasort=floor(length(data_sort)/10)*i;
    if isinf(norminv((ind_datasort+nb_zeros)/nb_pts))==0
        M_anamorphose(i+1,1)=data_sort(ind_datasort);
        M_anamorphose(i+1,2)=norminv((ind_datasort+nb_zeros)/nb_pts);
    else
         M_anamorphose(i+1,1)= M_anamorphose(i,1);
         M_anamorphose(i+1,2)= M_anamorphose(i,2);
    end  
end
par0(1)=1;
par0(2)=1;
[par]=fminsearch(@(par) calc_rms(par,M_anamorphose), par0);

function [rms]=calc_rms(par, M_anamorphose)
    obs=M_anamorphose(1:end,2);
    mod=M_anamorphose(1,2)+par(1)*(M_anamorphose(1:end,1).^par(2));
    rms=mean((obs-mod).^2);
end
 
%set initial parameters
m0=zeros(1,11);
m0(1)=5000;                 %spatial range
m0(2)=0.5;                  %space regularity parameter
m0(3)=800;                  %temporal range
m0(4)=0.5;                  %time regularity parameter
m0(5)=0.5;                  %space-time interaction parameter 
m0(6)=0.05;                 %measurement noise / sig_e
m0(7)=norminv(prop_zeros);  %Transform function; parameter a0
m0(8)=par(1);               %Transform function; parameter a1
m0(9)=par(2);               %Transform function; parameter a2
m0(10)=5.0;                 %advection velocity
m0(11)=0;                   %advection direction

%set the size of random walk steps
V_step_walk=zeros(11,1);
V_step_walk(1)=1000;        %spatial range
V_step_walk(2)=0.01;        %space regularity parameter
V_step_walk(3)=100;         %temporal range
V_step_walk(4)=0.01;        %time regularity parameter
V_step_walk(5)=0.01;        %space-time interaction parameter 
V_step_walk(6)=0.01;        %measurement noise / sig_e
%m(7) doesn't change. 
V_step_walk(8)=0.0001;       %Transform function; parameter a1
V_step_walk(9)=0.0001;       %Transform function; parameter a2
V_step_walk(10)=0.5;        %advection velocity
V_step_walk(11)=5.0;        %advection direction

%set priors
M_prior=zeros(11,2);
M_prior(1,1)=500;%spatial range - min
M_prior(1,2)=100000;%spatial range - max
M_prior(2,1)=0.1;%space regularity parameter- min
M_prior(2,2)=0.9;%space regularity parameter- max
M_prior(3,1)=100;%temporal range- min
M_prior(3,2)=15000;%temporal range- max
M_prior(4,1)=0.1;%time regularity parameter- min
M_prior(4,2)=0.9;%time regularity parameter- max
M_prior(5,1)=0.0;%space-time interaction parameter- min 
M_prior(5,2)=1.0;%space-time interaction parameter - max
M_prior(6,1)=0.01; %measurement noise / sig_e- min
M_prior(6,2)=0.5;%1.0; %measurement noise / sig_e- max
M_prior(8,1)=0.01;%a1- min
M_prior(8,2)=2.0;%a1- max
M_prior(9,1)=0.1;%a2- min
M_prior(9,2)=1.5;%a2- max
M_prior(10,1)=0.5;%0.5;%advection velocity- min
M_prior(10,2)=20.0;%advection velocity- max
M_prior(11,1)=-180;%advection direction- min
M_prior(11,2)=180;%advection direction- max


end