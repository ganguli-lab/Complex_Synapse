function [ bin_num ] = WhichBin( bin_edges, data )
%WHICHBIN(BIN_EDGES,DATA) In which bin (edges in BIN_EDGES) is DATA?
%   Detailed explanation goes here

bin_num=floor( interp1(bin_edges, 1:length(bin_edges), data) );
bin_num(bin_num==length(bin_edges))=length(bin_edges)-1;

end

