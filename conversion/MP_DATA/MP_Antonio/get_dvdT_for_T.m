function [dvdT_geo]=get_dvdT_for_T(model,T_geo,type);


T=load([ 'MP_' model '/T.dat' ]);
P=load([ 'MP_' model '/P_G.dat' ]);

load([ 'MP_' model '/d' type 'dT_MP_' model '.mat' ]);

if type=='vs'
    dvdT=dvsdT;
elseif type=='vp'
    dvdT=dvpdT;
end

size(dvdT);

for i=1:length(T_geo)

dvdT_geo(i)=interp1(T(1:80),dvdT(i,:)',T_geo(i));

end

%figure;

%plot(dvdT_geo,P);

%set(gca,'YDir','reverse');