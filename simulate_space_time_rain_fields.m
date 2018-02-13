function[M_Simul]=simulate_space_time_rain_fields(min_X, step_X, max_X, min_Y, step_Y, max_Y,V_sample,dt ,X_condi, Y_condi, nb_conditioning_times)

%check compatibility of V_t, X_condi and Y_condi with V_sample
%check V_t has constant steps

V_t=(0:dt:(length(V_sample(1,12:end))-1)/length(X_condi)*dt-1)';%add -1

size(V_t)
Coord_Euler=[];
for i=min_X:step_X:max_X
    for j=min_Y:step_Y:max_Y
        Coord_Euler=[Coord_Euler;[i,j]];
    end
end

[nb_realizations,~]=size(V_sample);

M_Simul=[];
%simulation
for ind_simul=1:nb_realizations
    display(strcat('simul:',num2str(ind_simul),'/',num2str(nb_realizations)))
    %hard data
    if ~isempty(X_condi)
        Condi_RainTS_gaussian=struct();
        for i=1:length(X_condi)
            Condi_RainTS_gaussian(i).X=X_condi(i);
            Condi_RainTS_gaussian(i).Y=Y_condi(i);
            Condi_RainTS_gaussian(i).t=V_t;
        
            for j=1:length(Condi_RainTS_gaussian(i).t)
                Condi_RainTS_gaussian(i).RainRate(j)=V_sample(ind_simul,11+(j-1)*length(X_condi)+i);
            end
            Condi_RainTS_gaussian(i).RainRate=Condi_RainTS_gaussian(i).RainRate';
        end
    
    else
        Condi_RainTS_gaussian=[];
    end

    %simulate latent field
    nb_ep=length(V_t);
    step_t=V_t(2)-V_t(1);
    %nb_conditioning_times=2;
    [SimulatedRainTS_gaussian]=simul_multigrid(Coord_Euler,Condi_RainTS_gaussian,V_sample(ind_simul,1:11), nb_ep, step_t,nb_conditioning_times);
    
    %transform to obtain simulated rainfall
    SimulatedRainTS_skewed=SimulatedRainTS_gaussian;
    for i=1:length(SimulatedRainTS_skewed)
        for j=1:length(SimulatedRainTS_skewed(i).RainRate)
            if SimulatedRainTS_gaussian(i).RainRate(j)<=V_sample(ind_simul,7)
                SimulatedRainTS_skewed(i).RainRate(j)=0;
            else
                SimulatedRainTS_skewed(i).RainRate(j)=((SimulatedRainTS_gaussian(i).RainRate(j)-V_sample(ind_simul,7))/V_sample(ind_simul,8))^(1/V_sample(ind_simul,9));
            end    
        end
    end
    M_Simul=[M_Simul;SimulatedRainTS_skewed];

    
end

end