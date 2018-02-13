function[MeasuredRainTS]=read_raw_data(input_file)
cd('Data')
M = readtable(input_file);
[sx,~]=size(M);


figure(1)
clf
figure(2)
clf

%load data and vizualize them
M_time=[];
M_rain=[];
M_coord=[];


for i=1:sx
    file=M.data_file{i};
    Mt=load(file);

    t = datetime(Mt(:,1),Mt(:,2),Mt(:,3),Mt(:,4),Mt(:,5),Mt(:,6));
    figure(1)
    hold on
 
    plot(t,Mt(:,7))
    M_time=[M_time; datenum(t')*86400];
    M_rain=[M_rain; Mt(:,7)'];
    
    X=M.Easting(i);
    Y=M.Northing(i);
    name=M.rain_gauge_description{i};
    txt2 = strcat('\leftarrow ',name);
    figure(2)
    hold on
    plot(X,Y,'k+')
    axis square
    text(X,Y,txt2)
    M_coord=[M_coord;[X, Y]];
end

%fill MeasuredRainTS structure
%!!in MeasuredRainTS everything is relativ (time, coordinates).
MeasuredRainTS=struct();
for i=1:sx
    MeasuredRainTS(i).t=M_time(i,:)'-M_time(i,1);
    MeasuredRainTS(i).RainRate=M_rain(i,:)';
    MeasuredRainTS(i).X=M_coord(i,1)-min(M_coord(:,1))+50;
    MeasuredRainTS(i).Y=M_coord(i,2)-min(M_coord(:,2))+50;
end
cd('..')
end