function [ toobj ] = CopyFields( fromstruct,toobj )
%COPYFIELDS utility for copying properties from struct fields. Compatible with arrays of
%objects
%needs to be made a private member function:
%     methods (Access=private)
%         toobj=CopyFields(fromstruct,toobj)
%     end %methods
% 
%     methods
%         %constructor
%         function obj=Foo(varargin)
%             switch nargin
%                 case 0
%                     %do nothing
%                 case 1
%                     [Unmatched,varargin]=extractArgOfType(varargin,'struct');
%                     if ~isempty(Unmatched)
%                         tempobj=CopyFields(Unmatched,tempobj);
%                     end
%                 otherwise
%                     error('Unknown inputs');
%             end %switch
%         end %constructor
%     end %methods

if ~isscalar(fromstruct)
    error('can only work with scalar originals');
end


fns=fields(fromstruct);
props=properties(toobj);
for i=1:length(fns)
   try
       if ismember(fns{i},props)
           [toobj.(fns{i})]=deal(fromstruct.(fns{i}));
       end
   catch exception
       if ~strcmp(exception.identifier,'MATLAB:class:SetProhibited');
           rethrow(exception);
       end%if
  end%try
end%for


end

