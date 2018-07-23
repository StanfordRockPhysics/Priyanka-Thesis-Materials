function [h1,slimits,cmap]=scatterqq(x,y,s,slimits,cmap)
%        [h1,slimits,cmap]=scatterqq(x,y,s,slimits,cmap)
%
% Makes a "quick" scatter plot.  Similar to Matlab's scatter, except that
% this version takes the data, divides it into 64 data ranges and makes
% regular plots with different colors.  Therefore, it is very much faster
% than scatter.
% This version will not work with LSELECT, because each color is a line
% with a different length.
% Alternate version SCATTERQ plots each color as a vector of length equal
% to the original data, with lots of NANs for the invisible points. 
% Inputs
%    x,y     arrays to plot along x,y axes
%    s       array to plot in color
%    slimits (optional) vector to specify range for color scale [smin, smax]
% Output
%    h1      handle to colorbar
%    slimits final color range [smin, smax]

% G. Mavko, April 2004

n2 = max([size(x,2) size(y,2) size(s,2)]);
if size(x,2)<n2, x=repmat(x,1,n2); end;
if size(y,2)<n2, y=repmat(y,1,n2); end;
if size(s,2)<n2, s=repmat(s,1,n2); end;

if nargin<4 | (nargin>3 & isempty(slimits)),
    smin = min(min(min(s))); 
    smax = max(max(max(s)));
    smin = smin - .01*(abs(smax-smin));
    smax = smax + .01*(abs(smax-smin));
    answer=inputdlg({'Min value of color range','Maximum value:'},'Quick Scatter Plot',1,{num2str(smin) num2str(smax)});
    smin = str2num(answer{1});
    smax = str2num(answer{2});
else,
    smin = slimits(1);
    smax = slimits(2);
end;
slimits = [smin smax];

if nargin==5,
    colormap(cmap); 
else,
    answer=inputdlg({'Colormap'},'Quick Scatter Plot',1,{'jet'});
    eval(['colormap(' answer{1} ')']);
    cmap = answer{1};
end;
linecolor = colormap; 
n = size(linecolor,1);
srange = linspace(smin, smax, n+1);

% this piece tricks the colorbar by making a dummy imagesc plot, then
% deleting it.
%h=imagesc(smin+(smax-smin)*rand(200,1));hold on;h1=colorbar;p=get(gca,'position');set(h1,'tag','colorbar');delete(h);

% now make the separate plots
inds = srange(2) > s;
h=plot(x(inds),y(inds),'o'); set(h,'markeredgecolor',linecolor(1,:),'markerfacecolor',linecolor(1,:),'color',linecolor(1,:),'markersize',4); 
hold on;
for k = 2: n-1
    inds = srange(k) <= s & srange(k+1) > s;
%    xp = nan*x; xp(inds) = x(inds);
%    yp = nan*y; yp(inds) = y(inds);
    h=plot(x(inds),y(inds),'o'); set(h,'markeredgecolor',linecolor(k,:),'markerfacecolor',linecolor(k,:),'color',linecolor(k,:),'markersize',4); 
    hold on;
end;
inds = srange(n) <= s; 
h=plot(x(inds),y(inds),'o'); set(h,'markeredgecolor',linecolor(n,:),'markerfacecolor',linecolor(n,:),'color',linecolor(n,:),'markersize',4); 
hold on;
set(gca,'tag','scatterplot');
caxis([smin smax]);
h1=colorbar;set(h1,'tag','colorbar');%delete(h);

% fix the axes
%p=get(gca,'position');
%set(gca,'position',p);
%axis xy
axis auto

