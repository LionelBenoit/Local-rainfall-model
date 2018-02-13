%%
close all
clear all
%%
%--------------------------------------------load and visualize data-------------------------------------------------------------------
%****************parameters to set********************
%Rq: Data files are in the Data folder               %*
input_file='metadata_Nant_Event2.csv';           %*
%**************end parameters to set******************
[MeasuredRainTS]=read_raw_data(input_file);

%for this dataset, two rain gauges did not work. Remove the corresponding time series
if strcmp(input_file,'metadata_Nant_Event2.csv')==1
    MeasuredRainTS(1)=[];
    MeasuredRainTS(5)=[];
end
%clean workspace
clear input_file;

%%
%----------------------------------Parameter inference - Metropolis sampling----------------------------------
%****************parameters to set********************
%sampling parameters                                %*
nb_ep_warmup=1000;                                  %*
nb_sample=100;                                      %*
step=10;                                            %*
%likelihood parameters                              %*
size_block=40; %composite likelihood - block size    %*
nb_iter_Gibbs=5; %nb iterations gibbs sampler       %*
%**************end parameters to set******************


%--set parameters within the function set_Metropolis_hastings is needed - but default parameters should give correct results--
[m0,V_step_walk,M_prior]=set_Metropolis_hastings(MeasuredRainTS);
%run Metropolis sampler
[V_tot, V_sample]=Metropolis_block_likelihood_gibbs(m0,V_step_walk,M_prior,MeasuredRainTS,size_block,nb_ep_warmup,nb_sample,step, nb_iter_Gibbs);
%plot results
plot_result_metropolis( V_tot, 3 )
%clean workspace
clear nb_ep_warmup nb_sample step nb_iter_Gibbs size_block m0 M_prior V_step_walk;


%%
%-------------------------------Unconditional simulation (1 single realization)------------------------------------------
%*************parameters to set*********************************************************
%index of model parameter vector used for simulation                                  %*
ind_sim=1;                                                                            %*
%points to simulate (here points were rain gauges were located)                       %*
Coord_Euler=[[872;653;488;652;1696;1402;523;50],[2912;51;1916;50;2713;2111;716;125]]; %*
%length of the simulated time series (here same length as observations)               %*
length_simulation=length(MeasuredRainTS(1).t);                                        %*
%temporal resolution (here same as data)                                              %*
time_resolution=MeasuredRainTS(1).t(2)-MeasuredRainTS(1).t(1);                        %*
%conditioning neighborhood                                                            %*
nb_conditioning_time_steps=4;                                                         %*
%***********end parameters to set*******************************************************


%simulation of latent field
V_m=V_sample(ind_sim,1:11);
[SimulatedRainTS_gaussian]=simul_multigrid(Coord_Euler,[],V_m,length_simulation, time_resolution,nb_conditioning_time_steps);
%add noise and back transform (because simul_multigrid generates latent and non noisy random fields)
SimulatedRainTS=SimulatedRainTS_gaussian;
for i=1:length(SimulatedRainTS_gaussian)
    for j=1:length(SimulatedRainTS_gaussian(i).RainRate)
        temp=SimulatedRainTS_gaussian(i).RainRate(j)+randn*V_m(6);
        if temp>V_m(7)
            SimulatedRainTS(i).RainRate(j)=((temp-V_m(7))/V_m(8))^(1/V_m(9));
        else
            SimulatedRainTS(i).RainRate(j)=0;
        end
    end
end
%plot the unconditional realization
figure(4)
clf
hold on
for i=1:length(SimulatedRainTS)
    plot(SimulatedRainTS(i).t,SimulatedRainTS(i).RainRate)
end
%clear workspace
clear Coord_Euler i ind_sim j length_simulation nb_conditioning_time_steps SimulatedRainTS_gaussian temp time_resolution V_m;


%%
%---------------------------------------Conditional simulation of dense space-time rain fileds------------------------------
%*******************parameters to set****************
%index of model parameter vector used for simulation*
ind_sim=1;                                         %*
%conditioning neighborhood                         %*
nb_conditioning_time_steps=4;                      %*
%area to simulate                                  %*
min_X=0;                                           %*
step_X=200;                                        %*
max_X=2000;                                        %*
min_Y=0;                                           %*
step_Y=200;                                        %*
max_Y=3000;                                        %*
%temporal resolution                               %*
dt_condi=60;                                       %*
%parameters for movie export                       %*                             
output_file='test_gif.gif';                        %*
speed=10;% in frames per second                    %*
max_color_scale=40; %max rain rate in color scale  %* 
%*************end parameters to set******************


%use observation dataset as conditioning dataset
V_sample_simul=V_sample(ind_sim,:);
X_condi=[];
Y_condi=[];
for i=1:length(MeasuredRainTS)
    X_condi=[X_condi;MeasuredRainTS(i).X];
    Y_condi=[Y_condi;MeasuredRainTS(i).Y];
end
%perform conditional simulation
[M_Simul]=simulate_space_time_rain_fields(min_X, step_X, max_X, min_Y, step_Y, max_Y, V_sample_simul, dt_condi , X_condi, Y_condi,nb_conditioning_time_steps);
%Create a movie (gif file) to visualize the realization
M_for_gif=zeros((max_Y-min_Y)/step_Y+1,(max_X-min_X)/step_X+1,length(M_Simul(1,1).t));
for t=1:length(M_Simul(1,1).t)
    my_ind=0;
    temp=[];
    for i=min_X:step_X:max_X
        V_Y=[];
        for j=min_Y:step_Y:max_Y
            my_ind=my_ind+1;
            V_Y=[V_Y;M_Simul(1,my_ind).RainRate(t)];
        end
        V_Y=flipud(V_Y);
        temp=[temp,V_Y];
    end
    M_for_gif(:,:,t)=temp;
end
make_gif( M_for_gif,output_file,speed,max_color_scale)
%clean workspace
clear dt_condi i j ind_sim M_Simul max_color_scale max_X max_Y min_X min_Y my_ind nb_conditioning_time_steps output_file speed step_X step_Y t temp V_sample_simul V_Y X_condi Y_condi;