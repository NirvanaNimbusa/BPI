function [A]=find_A(pre_fname, num_traces, T_amb)
% Estimates natural response matrix
%
% This function estimates the natural response matrix A in T[k]=AT[k-1]+P[k} 
%
% Authors:  S.Reda and A.Belouchrani
% Supported by US NAS Grant 2016, Brown
%
    Tk_1 = [];    
    Tk = [];

    for i=1:num_traces
        fname = sprintf('%s%d', pre_fname, i);
        D = dlmread(fname, '\t');
        Tk_1 = [Tk_1; D(1:end-1,:)];
        Tk = [Tk; D(2:end, :)];
    end

    Tk = Tk-T_amb;
    Tk_1 = Tk_1-T_amb;    
    Tk = Tk';
    Tk_1 = Tk_1';
    n = size(Tk, 1);
    A = zeros(n, n);
    for i=1:n    
        A(i,:)=lsqnonneg(Tk_1',Tk(i,:)');
    end
end
