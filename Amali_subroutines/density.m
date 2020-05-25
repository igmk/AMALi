function [Dns,DnsAE] = density(p,pAE,T,TAE)
%DENSITY  calculate molecular density from pressure and temperature
%
%           [Dns,DnsAE] = density(p,pAE,T,TAE)

% (c) GeBe 1996

  Dns = p ./ (T * boltzman);
  DnsAE = abserdiv(p,pAE,T,TAE) / boltzman;

