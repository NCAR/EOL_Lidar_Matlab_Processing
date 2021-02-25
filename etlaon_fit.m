function [C,mean_signal,FWHM, FSR, s, xx, locs, loc, val] = etlaon_fit(ind_start,ind_end, WavelengthTime,DetectorSignal, step)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  wave_comb = vertcat(WavelengthTime{ind_start:ind_end});   
  signal_comb = vertcat(DetectorSignal{ind_start:ind_end}); 

 [C,ia,idx] = unique(wave_comb );  %find the unique wavelengths and their indices 
 mean_signal = accumarray(idx,signal_comb,[],@mean); % calc mean at unique wavelengths
 xx = min(C):step:max(C);
 s = spline(C,mean_signal,xx);
 max(s)/2;
 idx = find(s > max(s)/2);
 idx_diff = find(diff(idx)>2);
 FWHM = (xx(idx(1)+idx_diff)- xx(idx(1)))*1000;  % in pm  
 [vals,locs] = findpeaks(s);
 [val,loc] = maxk(vals,2);
 FSR = xx(locs(loc(2)))-xx(locs(loc(1)));

end

