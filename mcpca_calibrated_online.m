% prepare the workspace and start clean
clear; close all; clc;

%% HARD CODE
% Hard code site number, iterations, and time span here
nrecords   = 7;               % How many lake records?
iterations = 100;             % How many Monte Carlo iterations? Start small until you know how long it will take!
iyear      = [1040:5:1950]';  % Over what common time period will we do the analysis? At what resolution?

% Hard code data paths here - you'll have to change these to match your system
rootPath = '/Users/kja/projects/africa/distribute/';

% where we are now ...
functionpath      = [rootPath]; % cd(functionpath);

% where the age model information is ...
agemodelPath      = [rootPath,'data/agemodels'];

% where the CALIB B00 files are ...
radiocarbonPath   = [rootPath,'data/14c'];

% where the actual proxy and depth data are ...
timeseriesPath    = [rootPath,'data/timeseries'];

addpath(functionpath);

%% GET AGE, DEPTH, AND DATA
cd(agemodelPath);   ageFiles  = wildfiles('*.txt');
cd(timeseriesPath); dataFiles = wildfiles('*.txt');
cd(functionpath);

%% PREALLOCATE FOR SPEED
published  = NaN(length(iyear),nrecords);
alliters   = NaN(length(iyear),nrecords,iterations);
member     = NaN(length(iyear),nrecords);     

tic % Check how long this takes ... (now would be a good time to develop new hobbies, interests; meet new people! Get out more!)

%% MCEOF PROCEDURE
for j = 1:iterations % Loop over the Monte Carlo iterations ...
 disp(['Interation ',num2str(j),' of ',num2str(iterations)])   
 ensemble = NaN(length(iyear),nrecords);
  for i = 1:nrecords % for each lake record ... 
   radiocarbonfound = 0;

   sitename = dataFiles{i}(1:6); 
   % disp(sitename)
   
   cd(timeseriesPath);
   data     = textread(dataFiles{i},'','headerlines',1); % get the time series data
     
     % Log transform the _highly_ non Gaussian charcoal data from Tanganyika
     if i == 6; data(:,3) = -log(-data(:,3)); end
     
   cd(agemodelPath);
   ages     = textread(ageFiles{i},'','headerlines',1);   
   cd(functionpath);
   
   if i == 2 || i > 3 % if the lake is NOT varved ...
    fullradiocarbonpath = [radiocarbonPath,'/',sitename];  % Get the directory with the CALIB B00 files ...
    cd(fullradiocarbonpath) % Go to the directory with the CALIB B00 files ...5
    c14Files = wildfiles('*.B00'); % get a list of all the CALIB B00 files there ...
   end
   
   % make a 'published' time series to compare to the MCPCA results later ..
    published(:,i) = interp1(data(:,2),data(:,3),iyear);
  
   original(i).Depth = ages(:,1);
   
   resampledDepth = [];
    if i == 5  % Naivasha has uncertainty in the depths, as well as the ages
     resampledDepth = normrnd(ages(:,1),ages(:,2));
     ages(:,2) = []; 
    else% the other time series are assumed to have no depth error
      resampledDepth = ages(:,1);
    end
    
   %% Age remodeling loop 
   original(i).Ages = ages(:,2);
   
   resampledAges = []; % preallocate for speed    
  for k = 1:size(ages,1)
      
   %% CALIBRATED AGE MODEL OUTPUT
    if i == 2 || i > 3  % for the other records, we'll resample from the actual calibration density function
          
     parfor c = 1:size(c14Files,1)
       calibrationAge(c) = str2double(c14Files{c}(end-7:end-4))
     end % end (c) loop over individual radiocarbon files
       cd(fullradiocarbonpath)
       if any(calibrationAge == ages(k,2))
        ic14 = find(calibrationAge==ages(k,2));
        radiocarbonfound = radiocarbonfound + 1;
        
        [c14pdf,~] = readB00(c14Files{ic14});
        
        if k > 1
         c14pdf(c14pdf(:,1)>=resampledAges(k-1,1),2) = 0;
        end
        
        if  nansum(c14pdf(:,2))==0      
          resampledAges(k,1) = resampledAges(k-1,1)-1;
        else
          resampledAges(k,1) = c14pdf(randsample(length(c14pdf),1,true,c14pdf(:,2)),1);
        end
        
       else % if there is no radiocarbon information, go ahead back to Gaussian resampling ... 
           
         if k == 1 % free reign to sample if it is the top date ... 
          resampledAges(k,1)  = normrnd(ages(k,2),ages(k,3));    
         else
             
            % Large Random Draw with age mean and measurement standard deviation
             R = normrnd(ages(k,2),ages(k,3),10e5,1);

             % Truncate this 'empirical distribution' based on resampled age upcore
             R(R>resampledAges(k-1,1)) = [];
      
             % Now, redraw from this distribution
              if isempty(R)
               resampledAges(k,1) = resampledAges(k-1,1)-1; % badDistributionFlag(k) = 2;
              elseif length(R) == 1
               resampledAges(k,1) = R;  % badDistributionFlag(k) = 1;
              else
               resampledAges(k,1)  = emprand(R); 
              end
         end   % end loop for CALIBRATED but no radiocarbon
       end   % end resampling loop in CALIBRATED  
     end % end loop with age model resampling from CALIBRATED 
    cd(functionpath);   
  end % end age model (k) loop ...
   

if i == 3 % for the Malawi varved record, we'll do cumulative varve counting error accounting
    alldepth = [0:708]';   
     for wx = 0:708
      windowedbias(wx+1,1) = normrnd(0,0.5);          
     end
     
      cumulativebias  = cumsum(windowedbias);
      [~,bindx,cindx] = intersect(ages(:,1),alldepth);
      resampledAges(:,1) = ages(:,2)+cumulativebias(cindx);

      % Double check to make sure we don't create a freak reversal
      if any(diff(resampledAges)>0)
         for wx = 0:708
          windowedbias(wx+1,1) = normrnd(0,0.5);          
         end
     
      cumulativebias  = cumsum(windowedbias);
      [~,bindx,cindx] = intersect(ages(:,1),alldepth);
      resampledAges(:,1) = ages(:,2)+cumulativebias(cindx); 
      end
       
end % if loop for Challa or Malawi Gaussian resampled ages ...


if i == 1  % for Challa's tightly packed varve dating
      
     uncertainty = 0.3;
     alldepth = ages(:,1);   
     
     
     for wx = 1:length(alldepth)
      windowedbias(wx,1) = normrnd(0,uncertainty);        
     end
     
      cumulativebias  = cumsum(windowedbias);
      [commonDepth,bindx,cindx] = intersect(data(:,1),alldepth);
      resampledAges(:,1) = data(bindx,2)+cumulativebias(cindx);;

      % Double check to make sure we don't create a freak reversal
      if any(diff(resampledAges)>0)
       for wx = 1:length(alldepth)
        windowedbias(wx,1) = normrnd(0,uncertainty);          
       end
     
      cumulativebias  = cumsum(windowedbias);
      % [~,bindx,cindx] = intersect(ages(:,1),alldepth);
      resampledAges(:,1) = data(bindx,2)+cumulativebias(cindx);;
      end
      
      resampledDepth = data(bindx,1);
      
end % if loop for Challa or Malawi Gaussian resampled ages ...


%% deal with Victoria's linear trend final date
if i == 7
 stats = regstats(resampledAges(1:end-1),resampledDepth(1:end-1),'linear');
 yhat = floor(resampledDepth(end) * stats.beta(2) + stats.beta(1));
 resampledAges(end) = yhat;
 % ci = ceil(sqrt(stats.mse)); resampledAges(end) = floor(normrnd(yhat,ci)); % Not used, but could add additional error to bottom from regression statistics
 victoriaend(j) = yhat;
end
  
   [otime,member]  = agemodel2(resampledAges, resampledDepth, data(:,1), data(:,3), iyear);
   
   if i == 3
    clip = find(iyear==1270); % cut off for Malawi short cores
    member(1:clip) = NaN;   
   end
   
   if i == 7
    clip = findnearest(resampledAges(end),iyear,1); % cut off Victoria
    member(1:clip) = NaN;   
   end
   
   alliters(:,i,j) = member;      
   ensemble(:,i)   = member;
   % disp(['Number of radiocarbon ages found and used: ', num2str(radiocarbonfound)])
  end % end over-record (i) loop 
  
  ensemble = standardize(ensemble); % make [0,1] 
  ensemble(iyear<1270,3) = NaN; % clip the (non-existant) part of Malawi
  ensemble(isnan(ensemble))=0;  % set missing values to zero
  [eofs,pc,expvar(j,:)] = caleof(standardize(ensemble),3,3);                            
  
  pc1.ensemble(j,:) = pc(1,:);
  pc2.ensemble(j,:) = pc(2,:);

  pc1.loadings(j,:) = eofs(1,:);
  pc2.loadings(j,:) = eofs(2,:);  

% For very long runs, requested by the Reviewers ... when do we converge?   
%   if j == 20000
%       save eastAfrica_mcpca_calibrated_2varve_20000.mat                                 
%
%   elseif j == 30000
%       save eastAfrica_mcpca_calibrated_2varve_30000.mat     
%
%   elseif j == 40000
%       save eastAfrica_mcpca_calibrated_2varve_40000.mat     
%       
%   elseif j == 50000
%       save eastAfrica_mcpca_calibrated_2varve_50000.mat     
%
%   elseif j == 60000
%       save eastAfrica_mcpca_calibrated_2varve_60000.mat     
%       
%   elseif j == 70000
%       save eastAfrica_mcpca_calibrated_2varve_70000.mat     
%       
%   elseif j == 80000
%       save eastAfrica_mcpca_calibrated_2varve_80000.mat     
%       
%   elseif j == 90000
%       save eastAfrica_mcpca_calibrated_2varve_90000.mat     
%   end

end % end Monte Carlo iteration loop

toc % Wow, that took a long time ...

% set raw, age modeled, interpolated, and PC of published data
publishedRaw = published;
published = standardize(published);
published(isnan(published)) = 0;
[peofs,ppc,pexpvar] = caleof(standardize(published),3,3);

save eastAfrica_mcpca_calibrated.mat                                 