function [ copy ] = CopyProps( original,copy )
%COPYPROPS utility for writing copy constructors. Compatible with arrays of
%objects
%needs to be made a private member function:
%     methods (Access=private)
%         copy=CopyProps(original,copy)
%     end %methods
% 
%     methods
%         %constructor
%         function obj=Foo(varargin)
%             switch nargin
%                 case 0
%                     %do nothing
%                 case 1
%                     if isa(varargin{1},'Foo')
%                         obj=CopyProps(varargin{1},obj); %copy constructor
%                     else
%                         error('Unknown inputs');
%                     end %if
%                 otherwise
%                     error('Unknown inputs');
%             end %switch
%         end %constructor
%     end %methods

if ~isscalar(original)
    error('can only work with scalar originals');
end


fns=properties(original);
for i=1:length(fns)
   try
       [copy.(fns{i})]=deal(original.(fns{i}));
   catch exception
       if ~strcmp(exception.identifier,'MATLAB:class:SetProhibited');
           rethrow(exception);
       end%if
  end%try
end%for


end

