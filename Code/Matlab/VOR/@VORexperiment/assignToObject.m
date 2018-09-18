function [s,x] = assignToObject(s, x)
%function [s,x] = assignToObject(s, x)
%
%if x{i} is a property of s, 
%s.(x{i}) = x{i+1}
%
%removes any parameter, value pairs from x and returns
%needs to be made member function is there are properties with private set
%access
% re-modified by Subhaneil Lahiri
% modified by marc gershow from pvpmod by
% (c) U. Egert 1998

%############################################
% this loop is assigns the parameter/value pairs in x to the calling
% workspace.
used = [];
p=properties(s);
if ~isempty(x)
    skipnext = false;
   for i = 1:size(x,2)
       if skipnext
           skipnext = false;
           continue;
       end
       if (ischar(x{i}) && ismember(x{i},p))         
          [s.(x{i})] = deal(x{i+1});
          used = [used i];
          skipnext = true;
       end
   end;
end;
if (~isempty(used))
    used = [used used+1];
    inds = setdiff(1:length(x), used);
    x = x(inds);
end

%############################################

