function cbtitle(titstr,fs)
%   cbtitle(titstr,fs)
%
% Puts title on a colorbar in the current figure
%
% titstr     is the desired title string
% fs         is the fontsize

% G Mavko

if nargin<2, fs=10; end;

% set(get(findobj(gcf,'tag','Colorbar'),'title'),'string',titstr,'fontsize',fs);
% ht = get(findobj(gcf,'tag','colorbar'),'title');
hcb = [findobj(gcf,'tag','colorbar') findobj(gcf,'tag','Colorbar')];

for k=1:length(hcb),
%     ht = findobj(hcb(k),'type','title');
%     set(ht,'string',titstr,'fontsize',fs);
    set(get(hcb(k),'title'),'string',titstr,'fontsize',fs);
end;
