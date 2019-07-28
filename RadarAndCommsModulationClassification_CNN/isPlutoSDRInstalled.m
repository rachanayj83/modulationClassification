function flag = isPlutoSDRInstalled
%isPlutoSDRInstalled Check if ADALM-PLUTO is installed

spkg = matlabshared.supportpkg.getInstalled;
flag = ~isempty(spkg) && any(contains({spkg.Name},'ADALM-PLUTO','IgnoreCase',true));
end