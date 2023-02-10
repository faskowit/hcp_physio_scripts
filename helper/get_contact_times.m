function [ictlens] = get_contact_times(binseq)
% given binary sequence, get the length of the consecutive 1's sequences

if ~all(ismember(binseq,[0 1]))
    error('needs to be binary sequence')
end

if isempty(binseq) || sum(binseq) == 0
    ictlens = nan ;
    return
end

n = length(binseq) ;
f1 = find(diff(binseq)==1) ; % step "up"
f2 = find(diff(binseq)==-1) ; % step "down"

% find out what the sequnce begins as
if binseq(1) % if its already up
    % have to add a 0 for already up
    f1 = [ 0 ; f1 ] ;
end

% find the shorter one...
if binseq(end) % ended on up
    % pretend like there is a step down at the end
    f2 = [ f2 ; n ] ;
end

ictlens = f2 - f1 ; 
