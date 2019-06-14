function [outputfinal]=ITCFitnlin2(fname,V0,SyC,T,npd,startka,startn,corrected,ptype)
%% ITC Data fitting %%
% ITC fitting of data using one set of binding sites based on Origin fitting function
% fname: filename obtained from ITC instrument in text format [xyz.DAT]
% V0: Volume of Cell in ml [1.42747]
% SyC: Syringe concentration in mM [0.3]
% T: Temperature of the ITC experiment in C [15]
% npd: Number of initial points to discard [2]
% startka: starting values for fit for Ka [8e5]
% startn: starting values for fit for n [1]
% corrected: flag for correcting the isotherm using last 4 injection points ['y']
% ptype: plot type in figure for the data points ['ko']

%% Read the input file
fileID = fopen(fname,'r');
dataArray = textscan(fileID, '%s%s%s%s%s%s%s%s%[^\n\r]' , 'Delimiter', '\t', 'HeaderLines' ,1, 'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2,3,4,5,6,7,8]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {0}; % Replace non-numeric cells

% Allocate imported array to column variable names
injv = cell2mat(raw(:, 2));
xt = cell2mat(raw(:, 3));
mt = cell2mat(raw(:, 4));
xmt = cell2mat(raw(:, 5));
ndh = cell2mat(raw(:, 6));

clearvars filename delimiter startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;
%% Correct the values of cell vol. and syringe conc. to L and M respectively
V0=V0*1e-3;% V0=Active Cell Volume (ml)
SyC=SyC*1e-3; %SyC=input('Concentration of Xt in Syringe (mM)

%% First injection point stored in I0
I0=[injv(npd+1)*1e-6 mt(npd+2)*1e-3 xmt(npd)];

%% Remove data points according to npd
injv=injv(npd+1:length(injv)-1);
xt=xt(npd+2:length(xt));
mt=mt(npd+2:length(mt));
xmt=xmt(npd+1:length(xmt)-1);
ndh=ndh(npd+1:length(ndh)-1);

%% Correct the baseline using the last 4 injection points if corrected is 'y' 
if corrected=='y'
ndh = ndh - mean(ndh(length(ndh)-4:length(ndh)));
end

% If the user has deleted points during ITC analysis, zero values are read
% and dealt with
foundzerondh=find(ndh==0);
injv(foundzerondh+1)=injv(foundzerondh+1)+injv(foundzerondh);
xt(foundzerondh)=[];mt(foundzerondh)=[];xmt(foundzerondh)=[];ndh(foundzerondh)=[];

%% Set up fittype and options.
opts=statset('nlinfit');
opts.UseParallel='true';
opts.MaxFunEvals=1e6;
opts.MaxIter=1e6;

s0=[startn startka ndh(1)-ndh(length(ndh))];
[beta,R,J,CovB]=nlinfit( [xmt, mt], ndh, @(par,x)onesitemodelf2(par,x(:,1),x(:,2),injv,V0,SyC,I0), s0);
stdcv=sqrt(diag(CovB))';

fitresult.n=beta(1); fitresult.Ka=beta(2); fitresult.dH=beta(3);
std.n=stdcv(1); std.Ka=stdcv(2); std.dH=stdcv(3);

Kd=1/fitresult.Ka;
erKd=std.Ka/(fitresult.Ka)^2; % Error Propogation of Ka to get errors of Kd
dS=-(-1.9858775*(273.15+T)*log(fitresult.Ka)-fitresult.dH)/(273.15+T); %dS = -(dG-dH)/T; dG = -RTln(Ka);

outputfinal=[Kd erKd fitresult.Ka std.Ka fitresult.n std.n fitresult.dH std.dH];
%% Plot fit with data.

figure( 'Name', ['ITC Fit of ' fname]);
title(['One Site Fitting of ' strrep(fname, '_', ' ')]);
% Create fake points in molar ratio and injection values
injf=2*ones(150,1); % 150 points with 2ul injection
[xtf,mtf,xmtf]=simulatedxtmt(injf.*1e-6,I0(1),I0(2),V0,SyC); % get simulated xt,mt and xmt
yfitf=onesitemodelf2([fitresult.n fitresult.Ka fitresult.dH],...
    xmtf,mtf.*1e3,injf,V0,SyC, ...
    [injf(1)*1e-6 mtf(npd)*1e-3 xmtf(npd)]); % created the fitted line
plot(xmtf,yfitf,'-','Color',[0.8 0.8 0.8],'LineWidth',3); % Plot the line

hold on;
plot(xmt,ndh,ptype); % plot the points of style ptype

% Comment/Uncomment to print fitted values on figure
x=0.95*mean(xmt);
y1=1.75*mean(ndh);
text(x,y1,{['K_{d} = (' sprintf('%5.2f',Kd*1e6) ' \pm ' sprintf('%5.2f',erKd*1e6) ' ) x 10^{-6} M'],['n = ' sprintf('%5.2f',fitresult.n) ' \pm ' sprintf('%5.2f',std.n) ],['\Delta' 'H = ' sprintf('%5.2f',fitresult.dH*1e-3) ' \pm ' sprintf('%5.2f',std.dH*1e-3) ' kcal/mol'],['\Delta' 'S = ' sprintf('%5.2f',dS) ' cal/mol']},'EdgeColor','red');
hold off;
end