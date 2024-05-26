function [ output_args ] = indexOf( input_args )
%INDEXOF Summary of this function goes here
%   Detailed explanation goes here

TEs = evalin('base','TEs');
output_args = find(TEs==input_args);

end

