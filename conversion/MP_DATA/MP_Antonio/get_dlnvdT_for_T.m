function [dlnvdT_geo]=get_dlnvdT_for_T(model,T_geo,type);


T=load([ 'MP_' model '/T.dat' ]);
P=load([ 'MP_' model '/P_G.dat' ]);

load([ 'MP_' model '/dln' type 'dT_MP_' model '.mat' ]);

if type=='vs'
    dlnvdT=dvsdT;
elseif type=='vp'
    dlnvdT=dvpdT;
elseif type=='vc'
    dlnvdT=dvcdT;
end

size(dlnvdT);

for i=1:length(T_geo)

dlnvdT_geo(i)=interp1(T(1:80),dlnvdT(i,:)',T_geo(i));

end

%figure;

%plot(dvdT_geo,P);

%set(gca,'YDir','reverse');