function [c14pdf,wmean] = readB00(filename,version)

% READb00 Read in CALIB calibrated radiocarbon ages
%
% Reads output from CALIB 5 for the probability density function.  
% This information is output in a file with the suffix with the format available here:
% http://radiocarbon.pa.qub.ac.uk/calib/manual/chapter3.html
%
%
% Inputs:
% filename = string; filename, possibly including full path, of the calibration file, including suffix 
% version = scalar; either 4 or 5 depending on CALIB version used (optional)
% 
% Outputs:
% c14pdf = N x 2 column with calibrated date and probability density
% wmean = weighted mean, rounded to nearest year
%

% References:
% 
% Reimer, P. J., Baillie, M. G. L., Bard, E., Bayliss, A., Beck, J. W., Bertrand, C. J. H., 
% Blackwell, P. G., Buck, C. E., Burr, G. S., Cutler, K. B., Damon, P. E., Edwards, R. L., 
% Fairbanks, R. G., Friedrich, M., Guilderson, T. P., Hogg, A. G., Hughen, K. A., Kromer, B., 
% McCormac, F. G., Manning, S. W., Ramsey, C. B., Reimer, R. W., Remmele, S., Southon, J. R., 
% Stuiver, M., Talamo, S., Taylor, F. W., van der Plicht, J., and Weyhenmeyer, C. E. 2004a. 
% IntCal04 Terrestrial radiocarbon age calibration, 26 - 0 ka BP. Radiocarbon 46(3), 1029-1058
%
% Stuiver, M., Reimer, P. J., and Reimer, R. W. 2005. CALIB 5.0. [WWW
% program and documentation].
%
% Stuiver, M., Reimer, P.J., 1993. Extended C-14 database and revised CALIB 3.0 C-14 age
% calibration program. Radiocarbon 35, 215-230. 
%
% Stuiver, M., Reimer, P.J., Bard, E., Beck, J.W., Burr, G.S., Hughen, K.A., Kromer, B., 
% McCormac, G., van der Plicht, J., Spurk, M., 1998. INTCAL98 radiocarbon age 6
% calibration, 24,000 -- 0 cal B.P.. Radiocarbon 40 (3), 1041-1083.
%
% Telford, R.J., Heegaard, E., Birks, H.J.B., 2004. The intercept is a poor
% estimate of a calibrated radiocarbon age. The Holocene 14(2), 296-298.
%
 
% error handling
if nargin < 2; version = 5;end;
if nargin == 0; error('Function Requires Valid Filename or Path'); end;

% Open the file, allow it to be read in text mode ('rt'). Need this ('t') for the PC, but not for Linux
fl=fopen(filename,'rt');

% initialize the count
count=0;

% read the file line by line, as long as the line is not yet the data, advance the count, got to the next line
if version == 4;
    while ~strcmp(' -99999.0',fgetl(fl))
        count=count+1;
    end   
 count = count+2; 
else
    while strcmp('#',strtok(fgetl(fl)))
        count=count+1;
    end
end

fclose(fl); % close the file

% open the file again, skip the header lines, read in two columns of numbers
[c14pdf]=textread(filename,'','headerlines',count,'delimiter',' '); 

% calculated the weighted mean; the pdf sums to unity
wmean = round(sum(c14pdf(:,1) .* c14pdf(:,2)));

% Type of date returned will be a function of options selected in CALIB, so
% no effort is made here to assign BP, AD, BC suffix, nor to convert
% between them.

