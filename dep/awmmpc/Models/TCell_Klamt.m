function [dx,p] = TCell_Klamt(t,x,p,u,drugs)
%TCELL_KLAMT Simulates ODEs for T cell model developed by Klamt (2006).
%   The following differential equations were automatically generated from
%   the original boolean model by the Odefy toolbox (using non-normalized
%   HillCubes as the continuous homologue for Boolean functions). Contains
%   all differential equations to run the complete model and can be called
%   by any of MATLAB's ODE solvers. If this file is called without initial
%   conditions, then the initial conditions are returned.
%   
%   SYNTAX:
%   [y0,p] = TCell_Klamt()
%   [T,Y] = odexxx(@TCell_Klamt,tspan,y0,options,p,u,drugs)
%   
%   INPUTS:
%   t: [scalar] time point at current step.
%   x: [vector] model states at current time point.
%   p: [structure] model parameter names and values.
%   u: [vector] administered dose of control drugs.
%   drugs: [structure] indicates drug characteristics (e.g.
%       response time, duration of action).
%   
%   OUTPUTS:
%   dx: [vector] reaction rates used to compute model states.
%   	OR
%       [vector] initial conditions (if nargin == 0).
%   p: [structure] model parameter names and values.
%   
%   Revision by: Jeffrey Perley (jperley@purdue.edu)
%   Last revision: 7/19/2012


% Determine if model parameters are needed and not provided
if ~exist('p','var') || isempty(p)
    %% Model Parameters
    
    % Define initial conditions (at equilibrium when not stimulated)
    x0 = zeros(40,1);
    x0([14,26]) = [1,1];
    x0(15) = 0.026290165530672;
    % Stimulation
    x0([18,36]) = [1,1];      % Receptor stimulation
    p.x0 = x0;
    
    % Define reaction parameters
    p.t_zap70 = 2.77;
    p.n_tcrphos = 2;
    p.k_tcrphos = 0.29;
    p.n_lck = 1;
    p.k_lck = 0.33;
    p.n_ccbl = 1;
    p.k_ccbl = 0.3;
    p.t_tcrphos = 0.03;
    p.n_tcr = 1;
    p.k_tcr = 0.34;
    p.n_lck_11 = 1;
    p.k_lck_12 = 0.34;
    p.n_fyn = 1;
    p.k_fyn = 0.13;
    p.t_tcr = 0.59;
    p.n_TCRlig = 1;
    p.k_TCRlig = 0.2;
    p.n_ccbl_18 = 1;
    p.k_ccbl_19 = 0.41;
    p.t_slp76 = 0.24;
    p.n_gads = 1;
    p.k_gads = 0.35;
    p.t_sek = 0.01;
    p.n_pkcth = 1;
    p.k_pkcth = 0.26;
    p.t_rsk = 1;
    p.n_erk = 3;
    p.k_erk = 0.3;
    p.t_rlk = 1.02;
    p.n_lck_30 = 1;
    p.k_lck_31 = 0.36;
    p.t_rasgrp = 2.76;
    p.n_pkcth_33 = 1;
    p.k_pkcth_34 = 0.4;
    p.n_dag = 1;
    p.k_dag = 0.39;
    p.t_ras = 7.5;
    p.n_rasgrp = 3;
    p.k_rasgrp = 0.33;
    p.n_grb2sos = 17;
    p.k_grb2sos = 0.39;
    p.t_raf = 11.8;
    p.n_ras = 3;
    p.k_ras = 0.31;
    p.t_plcgbind = 0.01;
    p.n_lat = 5;
    p.k_lat = 0.36;
    p.t_plcgact = 0.01;
    p.n_zap70 = 1;
    p.k_zap70 = 0.1;
    p.n_slp76 = 1;
    p.k_slp76 = 0.29;
    p.n_rlk = 1;
    p.k_rlk = 0.24;
    p.n_plcgbind = 1;
    p.k_plcgbind = 0.31;
    p.n_itk = 2;
    p.k_itk = 0.34;
    p.t_pkcth = 8.66;
    p.n_dag_60 = 2;
    p.k_dag_61 = 0.31;
    p.t_pagcsk = 1.33;
    p.n_tcr_63 = 1;
    p.k_tcr_64 = 0.32;
    p.n_fyn_65 = 1;
    p.k_fyn_66 = 0.69;
    p.t_nfkb = 1;
    p.n_ikb = 3;
    p.k_ikb = 0.3;
    p.t_nfat = 5.03;
    p.n_calcen = 1;
    p.k_calcen = 0.37;
    p.t_mek = 0.01;
    p.n_raf = 1;
    p.k_raf = 0.3324;%0.28;
    p.t_lck = 0.15;
    p.n_pagcsk = 2;
    p.k_pagcsk = 0.3;
    p.n_cd4 = 2;
    p.k_cd4 = 0.37;
    p.n_cd45 = 1;
    p.k_cd45 = 0.33;
    p.t_lat = 0.01;
    p.n_zap70_84 = 3;
    p.k_zap70_85 = 0.14;
    p.t_jun = 1;
    p.n_jnk = 3;
    p.k_jnk = 0.3;
    p.t_jnk = 0.08;
    p.n_sek = 4;
    p.k_sek = 0.57;
    p.t_itk = 0.01;
    p.n_zap70_93 = 2;
    p.k_zap70_94 = 0.1;
    p.n_slp76_95 = 1;
    p.k_slp76_96 = 0.31;
    p.t_ip3 = 0.01;
    p.n_plcgact = 2;
    p.k_plcgact = 0.38;
    p.t_ikkbeta = 10.55;
    p.n_pkcth_101 = 18;
    p.k_pkcth_102 = 0.33;
    p.t_ikb = 1;
    p.n_ikkbeta_104 = 3;
    p.k_ikkbeta_104 = 0.3;
    p.t_grb2sos = 5.26;
    p.n_lat_107 = 8;
    p.k_lat_108 = 0.46;
    p.t_gads = 0.12;
    p.n_lat_110 = 1;
    p.k_lat_111 = 0.36;
    p.t_fyn = 0.46;
    p.n_tcr_113 = 6;
    p.k_tcr_114 = 0.27;
    p.n_lck_115 = 2;
    p.k_lck_116 = 0.3;
    p.n_cd45_117 = 4;
    p.k_cd45_118 = 0.37;
    p.t_fos = 1;
    p.n_erk_120 = 3;
    p.k_erk_121 = 0.3;
    p.t_erk = 0.1;
    p.n_mek = 2;
    p.k_mek = 0.28;
    p.t_dag = 0.01;
    p.n_plcgact_126 = 3;
    p.k_plcgact_127 = 0.28;
    p.t_creb = 1;
    p.n_rsk = 3;
    p.k_rsk = 0.3;
    p.t_cre = 1;
    p.n_creb = 3;
    p.k_creb = 0.3;
    p.t_ccbl = 0.29;
    p.n_zap70_135 = 1;
    p.k_zap70_136 = 0.7;
    p.t_calcen = 4.53;
    p.n_ca = 11;
    p.k_ca = 0.3;
    p.t_ca = 0.2;
    p.n_ip3 = 3;
    p.k_ip3 = 0.33;
    p.t_ap1 = 1;
    p.n_jun = 3;
    p.k_jun = 0.3;
    p.n_fos = 3;
    p.k_fos = 0.3;
    
    % Additional parameters to fit experimental data in Jurkat cell line
    p.t_dx = 11.0503;%4.6704;           % dx = dx*p.t_dx % Scales all rates
    p.n_cd45_stim = 19.4718;%8.7008;    % Rate parameter for cd45
    p.k_cd45_stim = 0.2866;%0.4141;     % Rate parameter for cd45
    p.t_cd45_stim = 1;                  % Time constant for cd45
    
    % Parameters of control reagents
    p.drugs.sanguinarine.w = 0.9;
    p.drugs.sanguinarine.k = 1;
    p.drugs.u0126.w = 9;
    p.drugs.u0126.k = 1;
    %{
    p = zeros(1,147);
    p(1) = 2.77;        % t_zap70_1
    p(2) = 2;           % n_tcrphos_2
    p(3) = 0.29;        % k_tcrphos_3
    p(4) = 1;           % n_lck_4
    p(5) = 0.33;        % k_lck_5
    p(6) = 1;           % n_ccbl_6
    p(7) = 0.3;         % k_ccbl_7
    p(8) = 0.03;        % t_tcrphos_8
    p(9) = 1;           % n_tcr_9
    p(10) = 0.34;       % k_tcr_10
    p(11) = 1;          % n_lck_11
    p(12) = 0.34;       % k_lck_12
    p(13) = 1;          % n_fyn_13
    p(14) = 0.13;       % k_fyn_14
    p(15) = 0.59;       % t_tcr_15
    p(16) = 1;          % n_TCRlig_16
    p(17) = 0.2;        % k_TCRlig_17
    p(18) = 1;          % n_ccbl_18
    p(19) = 0.41;       % k_ccbl_19
    p(20) = 0.24;       % t_slp76_20
    p(21) = 1;          % n_gads_21
    p(22) = 0.35;       % k_gads_22
    p(23) = 0.01;       % t_sek_23
    p(24) = 1;          % n_pkcth_24
    p(25) = 0.26;       % k_pkcth_25
    p(26) = 1;          % t_rsk_26
    p(27) = 3;          % n_erk_27
    p(28) = 0.3;        % k_erk_28
    p(29) = 1.02;       % t_rlk_29
    p(30) = 1;          % n_lck_30
    p(31) = 0.36;       % k_lck_31
    p(32) = 2.76;       % t_rasgrp_32
    p(33) = 1;          % n_pkcth_33
    p(34) = 0.4;        % k_pkcth_34
    p(35) = 1;          % n_dag_35
    p(36) = 0.39;       % k_dag_36
    p(37) = 7.5;        % t_ras_37
    p(38) = 3;          % n_rasgrp_38
    p(39) = 0.33;       % k_rasgrp_39
    p(40) = 17;         % n_grb2sos_40
    p(41) = 0.39;       % k_grb2sos_41
    p(42) = 11.8;       % t_raf_42
    p(43) = 3;          % n_ras_43
    p(44) = 0.31;       % k_ras_44
    p(45) = 0.01;       % t_plcgbind_45
    p(46) = 5;          % n_lat_46
    p(47) = 0.36;       % k_lat_47
    p(48) = 0.01;       % t_plcgact_48
    p(49) = 1;          % n_zap70_49
    p(50) = 0.1;        % k_zap70_50
    p(51) = 1;          % n_slp76_51
    p(52) = 0.29;       % k_slp76_52
    p(53) = 1;          % n_rlk_53
    p(54) = 0.24;       % k_rlk_54
    p(55) = 1;          % n_plcgbind_55
    p(56) = 0.31;       % k_plcgbind_56
    p(57) = 2;          % n_itk_57
    p(58) = 0.34;       % k_itk_58
    p(59) = 8.66;       % t_pkcth_59
    p(60) = 2;          % n_dag_60
    p(61) = 0.31;       % k_dag_61
    p(62) = 1.33;       % t_pagcsk_62
    p(63) = 1;          % n_tcr_63
    p(64) = 0.32;       % k_tcr_64
    p(65) = 1;          % n_fyn_65
    p(66) = 0.69;       % k_fyn_66
    p(67) = 1;          % t_nfkb_67
    p(68) = 3;          % n_ikb_68
    p(69) = 0.3;        % k_ikb_69
    p(70) = 5.03;       % t_nfat_70
    p(71) = 1;          % n_calcen_71
    p(72) = 0.37;       % k_calcen_72
    p(73) = 0.01;       % t_mek_73
    p(74) = 1;          % n_raf_74
    p(75) = 0.28;       % k_raf_75
    p(76) = 0.15;       % t_lck_76
    p(77) = 2;          % n_pagcsk_77
    p(78) = 0.3;        % k_pagcsk_78
    p(79) = 2;          % n_cd4_79
    p(80) = 0.37;       % k_cd4_80
    p(81) = 1;          % n_cd45_81
    p(82) = 0.33;       % k_cd45_82
    p(83) = 0.01;       % t_lat_83
    p(84) = 3;          % n_zap70_84
    p(85) = 0.14;       % k_zap70_85
    p(86) = 1;          % t_jun_86
    p(87) = 3;          % n_jnk_87
    p(88) = 0.3;        % k_jnk_88
    p(89) = 0.08;       % t_jnk_89
    p(90) = 4;          % n_sek_90
    p(91) = 0.57;       % k_sek_91
    p(92) = 0.01;       % t_itk_92
    p(93) = 2;          % n_zap70_93
    p(94) = 0.1;        % k_zap70_94
    p(95) = 1;          % n_slp76_95
    p(96) = 0.31;       % k_slp76_96
    p(97) = 0.01;       % t_ip3_97
    p(98) = 2;          % n_plcgact_98
    p(99) = 0.38;       % k_plcgact_99
    p(100) = 10.55;     % t_ikkbeta_100
    p(101) = 18;        % n_pkcth_101
    p(102) = 0.33;      % k_pkcth_102
    p(103) = 1;         % t_ikb_103
    p(104) = 3;         % n_ikkbeta_104
    p(105) = 0.3;       % k_ikkbeta_105
    p(106) = 5.26;      % t_grb2sos_106
    p(107) = 8;         % n_lat_107
    p(108) = 0.46;      % k_lat_108
    p(109) = 0.12;      % t_gads_109
    p(110) = 1;         % n_lat_110
    p(111) = 0.36;      % k_lat_111
    p(112) = 0.46;      % t_fyn_112
    p(113) = 6;         % n_tcr_113
    p(114) = 0.27;      % k_tcr_114
    p(115) = 2;         % n_lck_115
    p(116) = 0.3;       % k_lck_116
    p(117) = 4;         % n_cd45_117
    p(118) = 0.37;      % k_cd45_118
    p(119) = 1;         % t_fos_119
    p(120) = 3;         % n_erk_120
    p(121) = 0.3;       % k_erk_121
    p(122) = 0.1;       % t_erk_122
    p(123) = 2;         % n_mek_123
    p(124) = 0.28;      % k_mek_124
    p(125) = 0.01;      % t_dag_125
    p(126) = 3;         % n_plcgact_126
    p(127) = 0.28;      % k_plcgact_127
    p(128) = 1;         % t_creb_128
    p(129) = 3;         % n_rsk_129
    p(130) = 0.3;       % k_rsk_130
    p(131) = 1;         % t_cre_131
    p(132) = 3;         % n_creb_132
    p(133) = 0.3;       % k_creb_133
    p(134) = 0.29;      % t_ccbl_134
    p(135) = 1;         % n_zap70_135
    p(136) = 0.7;       % k_zap70_136
    p(137) = 4.53;      % t_calcen_137
    p(138) = 11;        % n_ca_138
    p(139) = 0.3;       % k_ca_139
    p(140) = 0.2;       % t_ca_140
    p(141) = 3;         % n_ip3_141
    p(142) = 0.33;      % k_ip3_142
    p(143) = 1;         % t_ap1_143
    p(144) = 3;         % n_jun_144
    p(145) = 0.3;       % k_jun_145
    p(146) = 3;         % n_fos_146
    p(147) = 0.3;       % k_fos_147
    %}
    
end


%% Initial Conditions, Reaction Rates & Ordinary Differential Equations

% Set initial conditions if called with no current states
n = 40;     % Number of substrates
if nargin == 0
    %% Set initial conditions
    
    % Return initial conditions
    dx = p.x0;
    
    
else
    %% Evaluate System of ODEs
    
    % Simulate drug activity
    if ~isempty(u) && ~isempty(drugs)
        %% Simulate drug effects on kinetic parameters
        for j = 1:length(drugs)
            switch lower(drugs(j).name)
                case 'sanguinarine' % MAPKPH inhibitor (Sanguinarine)
                    drugs(j).w = p.drugs.sanguinarine.w;
                    p.drugs.sanguinarine.k = sim_Drug(p.drugs.sanguinarine.k,u(:,j),t,drugs(j));
                case 'u0126'        % MEK inhibitor (U0126)
                    drugs(j).w = p.drugs.u0126.w;
                    p.drugs.u0126.k = sim_Drug(p.drugs.u0126.k,u(:,j),t,drugs(j));
                case 'pma'          % PKC stimulant (PMA)
                case 'pp1'          % LCK inhibitor (PP1)
            end
        end
    end
    
    % Ensure states are within interval [0 1]
%     x(x<0) = 0;
%     x(x>1) = 1;
    
    % Set variable indices
    zap70 = 1;
    tcrphosp = 2;
    tcr = 3;
    slp76 = 4;
    sek = 5;
    rsk = 6;
    rlk = 7;
    rasgrp = 8;
    ras = 9;
    raf = 10;
    plcgbind = 11;
    plcgact = 12;
    pkcth = 13;
    pagcsk = 14;
    nfkb = 15;
    nfat = 16;
    mek = 17;
    TCRlig = 18;
    lck = 19;
    lat = 20;
    jun = 21;
    jnk = 22;
    itk = 23;
    ip3 = 24;
    ikkbeta = 25;
    ikb = 26;
    grb2sos = 27;
    gads = 28;
    fyn = 29;
    fos = 30;
    erk = 31;
    dag = 32;
    creb = 33;
    cre = 34;
    cd4 = 35;
    cd45 = 36;
    ccbl = 37;
    calcen = 38;
    ca = 39;
    ap1 = 40;
    
    % System of ordinary differential equations
    dx = zeros(n,1);
    dx(zap70) = (x(tcrphosp)^p.n_tcrphos/(x(tcrphosp)^p.n_tcrphos+p.k_tcrphos^p.n_tcrphos)*x(lck)^p.n_lck/(x(lck)^p.n_lck+p.k_lck^p.n_lck)*(1-x(ccbl)^p.n_ccbl/(x(ccbl)^p.n_ccbl+p.k_ccbl^p.n_ccbl))-x(zap70)) / p.t_zap70;
    dx(tcrphosp) = (x(tcr)^p.n_tcr/(x(tcr)^p.n_tcr+p.k_tcr^p.n_tcr)*x(lck)^p.n_lck_11/(x(lck)^p.n_lck_11+p.k_lck_12^p.n_lck_11)*(1-x(fyn)^p.n_fyn/(x(fyn)^p.n_fyn+p.k_fyn^p.n_fyn))+(1-x(tcr)^p.n_tcr/(x(tcr)^p.n_tcr+p.k_tcr^p.n_tcr))*(1-x(lck)^p.n_lck_11/(x(lck)^p.n_lck_11+p.k_lck_12^p.n_lck_11))*x(fyn)^p.n_fyn/(x(fyn)^p.n_fyn+p.k_fyn^p.n_fyn)+x(tcr)^p.n_tcr/(x(tcr)^p.n_tcr+p.k_tcr^p.n_tcr)*(1-x(lck)^p.n_lck_11/(x(lck)^p.n_lck_11+p.k_lck_12^p.n_lck_11))*x(fyn)^p.n_fyn/(x(fyn)^p.n_fyn+p.k_fyn^p.n_fyn)+(1-x(tcr)^p.n_tcr/(x(tcr)^p.n_tcr+p.k_tcr^p.n_tcr))*x(lck)^p.n_lck_11/(x(lck)^p.n_lck_11+p.k_lck_12^p.n_lck_11)*x(fyn)^p.n_fyn/(x(fyn)^p.n_fyn+p.k_fyn^p.n_fyn)+x(tcr)^p.n_tcr/(x(tcr)^p.n_tcr+p.k_tcr^p.n_tcr)*x(lck)^p.n_lck_11/(x(lck)^p.n_lck_11+p.k_lck_12^p.n_lck_11)*x(fyn)^p.n_fyn/(x(fyn)^p.n_fyn+p.k_fyn^p.n_fyn)-x(tcrphosp)) / p.t_tcrphos;
    dx(tcr) = (x(TCRlig)^p.n_TCRlig/(x(TCRlig)^p.n_TCRlig+p.k_TCRlig^p.n_TCRlig)*(1-x(ccbl)^p.n_ccbl_18/(x(ccbl)^p.n_ccbl_18+p.k_ccbl_19^p.n_ccbl_18))-x(tcr)) / p.t_tcr;
    dx(slp76) = (x(gads)^p.n_gads/(x(gads)^p.n_gads+p.k_gads^p.n_gads)-x(slp76)) / p.t_slp76;
    dx(sek) = (x(pkcth)^p.n_pkcth/(x(pkcth)^p.n_pkcth+p.k_pkcth^p.n_pkcth)-x(sek)) / p.t_sek;
    dx(rsk) = (x(erk)^p.n_erk/(x(erk)^p.n_erk+p.k_erk^p.n_erk)-x(rsk)) / p.t_rsk;
    dx(rlk) = (x(lck)^p.n_lck_30/(x(lck)^p.n_lck_30+p.k_lck_31^p.n_lck_30)-x(rlk)) / p.t_rlk;
    dx(rasgrp) = (x(pkcth)^p.n_pkcth_33/(x(pkcth)^p.n_pkcth_33+p.k_pkcth_34^p.n_pkcth_33)*x(dag)^p.n_dag/(x(dag)^p.n_dag+p.k_dag^p.n_dag)-x(rasgrp)) / p.t_rasgrp;
    dx(ras) = (x(rasgrp)^p.n_rasgrp/(x(rasgrp)^p.n_rasgrp+p.k_rasgrp^p.n_rasgrp)*(1-x(grb2sos)^p.n_grb2sos/(x(grb2sos)^p.n_grb2sos+p.k_grb2sos^p.n_grb2sos))+(1-x(rasgrp)^p.n_rasgrp/(x(rasgrp)^p.n_rasgrp+p.k_rasgrp^p.n_rasgrp))*x(grb2sos)^p.n_grb2sos/(x(grb2sos)^p.n_grb2sos+p.k_grb2sos^p.n_grb2sos)+x(rasgrp)^p.n_rasgrp/(x(rasgrp)^p.n_rasgrp+p.k_rasgrp^p.n_rasgrp)*x(grb2sos)^p.n_grb2sos/(x(grb2sos)^p.n_grb2sos+p.k_grb2sos^p.n_grb2sos)-x(ras)) / p.t_ras;
    dx(raf) = (x(ras)^p.n_ras/(x(ras)^p.n_ras+p.k_ras^p.n_ras)-x(raf)) / p.t_raf;
    dx(plcgbind) = (x(lat)^p.n_lat/(x(lat)^p.n_lat+p.k_lat^p.n_lat)-x(plcgbind)) / p.t_plcgbind;
    dx(plcgact) = (x(zap70)^p.n_zap70/(x(zap70)^p.n_zap70+p.k_zap70^p.n_zap70)*x(slp76)^p.n_slp76/(x(slp76)^p.n_slp76+p.k_slp76^p.n_slp76)*x(rlk)^p.n_rlk/(x(rlk)^p.n_rlk+p.k_rlk^p.n_rlk)*x(plcgbind)^p.n_plcgbind/(x(plcgbind)^p.n_plcgbind+p.k_plcgbind^p.n_plcgbind)*(1-x(itk)^p.n_itk/(x(itk)^p.n_itk+p.k_itk^p.n_itk))+x(zap70)^p.n_zap70/(x(zap70)^p.n_zap70+p.k_zap70^p.n_zap70)*x(slp76)^p.n_slp76/(x(slp76)^p.n_slp76+p.k_slp76^p.n_slp76)*(1-x(rlk)^p.n_rlk/(x(rlk)^p.n_rlk+p.k_rlk^p.n_rlk))*x(plcgbind)^p.n_plcgbind/(x(plcgbind)^p.n_plcgbind+p.k_plcgbind^p.n_plcgbind)*x(itk)^p.n_itk/(x(itk)^p.n_itk+p.k_itk^p.n_itk)+x(zap70)^p.n_zap70/(x(zap70)^p.n_zap70+p.k_zap70^p.n_zap70)*x(slp76)^p.n_slp76/(x(slp76)^p.n_slp76+p.k_slp76^p.n_slp76)*x(rlk)^p.n_rlk/(x(rlk)^p.n_rlk+p.k_rlk^p.n_rlk)*x(plcgbind)^p.n_plcgbind/(x(plcgbind)^p.n_plcgbind+p.k_plcgbind^p.n_plcgbind)*x(itk)^p.n_itk/(x(itk)^p.n_itk+p.k_itk^p.n_itk)-x(plcgact)) / p.t_plcgact;
    dx(pkcth) = (x(dag)^p.n_dag_60/(x(dag)^p.n_dag_60+p.k_dag_61^p.n_dag_60)-x(pkcth)) / p.t_pkcth;
    dx(pagcsk) = ((1-x(tcr)^p.n_tcr_63/(x(tcr)^p.n_tcr_63+p.k_tcr_64^p.n_tcr_63))*(1-x(fyn)^p.n_fyn_65/(x(fyn)^p.n_fyn_65+p.k_fyn_66^p.n_fyn_65))+(1-x(tcr)^p.n_tcr_63/(x(tcr)^p.n_tcr_63+p.k_tcr_64^p.n_tcr_63))*x(fyn)^p.n_fyn_65/(x(fyn)^p.n_fyn_65+p.k_fyn_66^p.n_fyn_65)+x(tcr)^p.n_tcr_63/(x(tcr)^p.n_tcr_63+p.k_tcr_64^p.n_tcr_63)*x(fyn)^p.n_fyn_65/(x(fyn)^p.n_fyn_65+p.k_fyn_66^p.n_fyn_65)-x(pagcsk)) / p.t_pagcsk;
    dx(nfkb) = ((1-x(ikb)^p.n_ikb/(x(ikb)^p.n_ikb+p.k_ikb^p.n_ikb))-x(nfkb)) / p.t_nfkb;
    dx(nfat) = (x(calcen)^p.n_calcen/(x(calcen)^p.n_calcen+p.k_calcen^p.n_calcen)-x(nfat)) / p.t_nfat;
    dx(mek) = (x(raf)^p.n_raf/(x(raf)^p.n_raf+p.k_raf^p.n_raf)-x(mek)) / p.t_mek;
    dx(TCRlig) = 0;
%     dx(lck) = ((1-x(pagcsk)^p.n_pagcsk/(x(pagcsk)^p.n_pagcsk+p.k_pagcsk^p.n_pagcsk))*x(cd4)^p.n_cd4/(x(cd4)^p.n_cd4+p.k_cd4^p.n_cd4)*x(cd45)^p.n_cd45/(x(cd45)^p.n_cd45+p.k_cd45^p.n_cd45)-x(lck)) / p.t_lck;% Original
    dx(lck) = ((1-x(pagcsk)^p.n_pagcsk/(x(pagcsk)^p.n_pagcsk+p.k_pagcsk^p.n_pagcsk))*x(cd45)^p.n_cd45/(x(cd45)^p.n_cd45+p.k_cd45^p.n_cd45)-x(lck)) / p.t_lck;% cd4 removed
    dx(lat) = (x(zap70)^p.n_zap70_84/(x(zap70)^p.n_zap70_84+p.k_zap70_85^p.n_zap70_84)-x(lat)) / p.t_lat;
    dx(jun) = (x(jnk)^p.n_jnk/(x(jnk)^p.n_jnk+p.k_jnk^p.n_jnk)-x(jun)) / p.t_jun;
    dx(jnk) = (x(sek)^p.n_sek/(x(sek)^p.n_sek+p.k_sek^p.n_sek)-x(jnk)) / p.t_jnk;
    dx(itk) = (x(zap70)^p.n_zap70_93/(x(zap70)^p.n_zap70_93+p.k_zap70_94^p.n_zap70_93)*x(slp76)^p.n_slp76_95/(x(slp76)^p.n_slp76_95+p.k_slp76_96^p.n_slp76_95)-x(itk)) / p.t_itk;
    dx(ip3) = (x(plcgact)^p.n_plcgact/(x(plcgact)^p.n_plcgact+p.k_plcgact^p.n_plcgact)-x(ip3)) / p.t_ip3;
    dx(ikkbeta) = (x(pkcth)^p.n_pkcth_101/(x(pkcth)^p.n_pkcth_101+p.k_pkcth_102^p.n_pkcth_101)-x(ikkbeta)) / p.t_ikkbeta;
    dx(ikb) = ((1-x(ikkbeta)^p.n_ikkbeta_104/(x(ikkbeta)^p.n_ikkbeta_104+p.k_ikkbeta_104^p.n_ikkbeta_104))-x(ikb)) / p.t_ikb;
    dx(grb2sos) = (x(lat)^p.n_lat_107/(x(lat)^p.n_lat_107+p.k_lat_108^p.n_lat_107)-x(grb2sos)) / p.t_grb2sos;
    dx(gads) = (x(lat)^p.n_lat_110/(x(lat)^p.n_lat_110+p.k_lat_111^p.n_lat_110)-x(gads)) / p.t_gads;
    dx(fyn) = (x(tcr)^p.n_tcr_113/(x(tcr)^p.n_tcr_113+p.k_tcr_114^p.n_tcr_113)*(1-x(lck)^p.n_lck_115/(x(lck)^p.n_lck_115+p.k_lck_116^p.n_lck_115))*x(cd45)^p.n_cd45_117/(x(cd45)^p.n_cd45_117+p.k_cd45_118^p.n_cd45_117)+(1-x(tcr)^p.n_tcr_113/(x(tcr)^p.n_tcr_113+p.k_tcr_114^p.n_tcr_113))*x(lck)^p.n_lck_115/(x(lck)^p.n_lck_115+p.k_lck_116^p.n_lck_115)*x(cd45)^p.n_cd45_117/(x(cd45)^p.n_cd45_117+p.k_cd45_118^p.n_cd45_117)+x(tcr)^p.n_tcr_113/(x(tcr)^p.n_tcr_113+p.k_tcr_114^p.n_tcr_113)*x(lck)^p.n_lck_115/(x(lck)^p.n_lck_115+p.k_lck_116^p.n_lck_115)*x(cd45)^p.n_cd45_117/(x(cd45)^p.n_cd45_117+p.k_cd45_118^p.n_cd45_117)-x(fyn)) / p.t_fyn;
    dx(fos) = (x(erk)^p.n_erk_120/(x(erk)^p.n_erk_120+p.k_erk_121^p.n_erk_120)-x(fos)) / p.t_fos;
%     dx(erk) = (x(mek)^p.n_mek/(x(mek)^p.n_mek+p.k_mek^p.n_mek)-x(erk)) / p.t_erk;
    dx(erk) = (p.drugs.u0126.k*x(mek)^p.n_mek/(x(mek)^p.n_mek+p.k_mek^p.n_mek)-p.drugs.sanguinarine.k*x(erk)) / p.t_erk;
    dx(dag) = (x(plcgact)^p.n_plcgact_126/(x(plcgact)^p.n_plcgact_126+p.k_plcgact_127^p.n_plcgact_126)-x(dag)) / p.t_dag;
    dx(creb) = (x(rsk)^p.n_rsk/(x(rsk)^p.n_rsk+p.k_rsk^p.n_rsk)-x(creb)) / p.t_creb;
    dx(cre) = (x(creb)^p.n_creb/(x(creb)^p.n_creb+p.k_creb^p.n_creb)-x(cre)) / p.t_cre;
    dx(cd4) = 0;
    dx(cd45) = -(x(cd45)^p.n_cd45_stim/(x(cd45)^p.n_cd45_stim+p.k_cd45_stim^p.n_cd45_stim)) / p.t_cd45_stim;
    dx(ccbl) = (x(zap70)^p.n_zap70_135/(x(zap70)^p.n_zap70_135+p.k_zap70_136^p.n_zap70_135)-x(ccbl)) / p.t_ccbl;
    dx(calcen) = (x(ca)^p.n_ca/(x(ca)^p.n_ca+p.k_ca^p.n_ca)-x(calcen)) / p.t_calcen;
    dx(ca) = (x(ip3)^p.n_ip3/(x(ip3)^p.n_ip3+p.k_ip3^p.n_ip3)-x(ca)) / p.t_ca;
    dx(ap1) = (x(jun)^p.n_jun/(x(jun)^p.n_jun+p.k_jun^p.n_jun)*x(fos)^p.n_fos/(x(fos)^p.n_fos+p.k_fos^p.n_fos)-x(ap1)) / p.t_ap1;
    
    dx = dx*p.t_dx;
    
end

end