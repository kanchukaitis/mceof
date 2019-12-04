function [otime,odata,otrough] = agemodel(itime, idepth, depth, data,otime) 

% AGEMODEL  function [otime,odata,otrough] = agemodel(itime, idepth, depth, data, otime) 
%
% This function takes data with given tie-point dates, specified at certain depths, and makes a linear timeline
% between them.  It then interpolates to a specified constant time step (i.e monthly).  
%
%  [otime,odata] = agemodel(itime, idepth, depth, data, otime)         
%
% otime  = output time scale, evenly spaced at a given constant time step (i.e monthly)
% odata  = output data that corresponds to that time vector
% otrought = output data age modeled but no consistent time step (debugging)
%
% itime  = input time markers (i.e. tie points) (1 x p)
% idepth = corresponding input depth marks for those tie points (1 x p)
% depth  = the depth that every measurement in 'data' was made (1 x m)
% data   = the data we are age modeling (1 x m)
% otime   = the constant time steps between each otime.  (i.e. 1/12) (1 x 1)
 
% KJA 07/2006 1st version, hacked from Peter Huybers' original code


%% First, interpolate fully between tie points so that every data point also has a time 
for place = 1:(length(idepth)-1)
  t1 = itime(place);
  t2 = itime(place+1);
  d1 = idepth(place);
  d2 = idepth(place+1);
  [d,pd1] =min(abs(depth(:,1)-d1));
  [d,pd2] =min(abs(depth(:,1)-d2));
  if place+1 == length(idepth), pd2 = length(depth); end; % deal with end points
   
 % Now, Linear Interpolation, simply: (percent of depth) * (time interval) + (begin time)
 for row = pd1:pd2
     %otrough(row) =  (depth(row,1)-d1)/(d2-d1) * (t2-t1) + t1;
     % otrough =  interp1(idepth,itime,depth,'linear','extrap');
     otrough =  interp1(idepth,itime,depth,'pchip','extrap');
   end; % end the interpolating loop
   
end; % end the time/depth loop

%otrough(otrough == 0) = NaN; % otrough = otrough'

%% Now, use interpolate again to make time steps which are at a constant spacing
%if otrough(max(find(~isnan(otrough)))) > otrough(min(find(~isnan(otrough)))) % if the larger date is later
%  otime = min(itime):step:max(itime);
% else
%  otime = max(itime):-step:min(itime);
%end

%otime = otime';
odata = interp1(otrough,data,otime,'linear','extrap');
