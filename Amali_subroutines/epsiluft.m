function Res = epsiluft()
%EPSILUFT  Depolarisationskonstante Epsilon von Luft ist 0.22
%          Quelle: A. T. Young, Physics Today, p42, Jan 82.
%                  A. T. Young, Appl. Opt, (1980) 19 p3427.

%	Gebe 27-03-92

  rho_to = 0.0279;  % f"ur trockene Luft

  Res = 45 * rho_to / (6 - 7 * rho_to);


