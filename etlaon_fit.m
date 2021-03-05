function [C,mean_signal,FWHM, FSR, s, xx, locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  wave_comb = vertcat(WavelengthTime{ind_start:ind_end});   
  signal_comb = vertcat(DetectorSignal{ind_start:ind_end}); 

 [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
 mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths
 xx = min(C):step:max(C);
 s = pchip(C,mean_signal,xx);


 [vals,locs] = findpeaks(s);
 [val,loc] = maxk(vals,3);
 FSR_1 = abs((xx(locs(loc(2)))-xx(locs(loc(1))))*1000);  % in pm  
 FSR_2 = abs((xx(locs(loc(3))))-xx(locs(loc(1)))*1000);  % in pm  
 FSR = min(FSR_1, FSR_2);
 
 max(s)/2;
 idx = find(s > max(s)/2);
 idx_diff = find(diff(idx)>2);
 FWHM = (xx(idx(1)+idx_diff)- xx(idx(1)))*1000;  % in pm  
 
 
end

