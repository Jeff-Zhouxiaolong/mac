// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2008 - INRIA - Allan CORNET
// Copyright (C) 2011 - DIGITEO - Allan CORNET
//
// This file is released under the 3-clause BSD license. See COPYING-BSD.

function subdemolist = demo_gateway()

demopath = get_absolute_file_path("emd_toolbox.dem.gateway.sce");

subdemolist=['EMD for the sum of 2 sinusoidal FMs + 1 Gaussian logon', 'emd_fmsin.sce';..
'error measure in the EMD estimation of a single tone', 'emd_sampling.sce';...
'error measure in the EMD separation of two tones', 'emd_separation.sce';...
'EMD for the sum of 2 triangular  waveforms + 1 tone' , 'emd_triang.sce';...
'on-line EMD',  'ex_online.sce';...
'Hilbert-Huang spectrum of an EMD',  'emd_hhspectrum.sce';...
'bivariate EMD extension on a real-world oceanographic signal', 'bivariate_EMD_principle.sce';...
'shows principle of the bivariate EMD extension', 'bivariate_EMD_illustration.sce';...
'shows how the center of the tube envelope is defined', 'bivariate_EMD_mean_definitions.sce';...
];


subdemolist(:,2) = demopath + subdemolist(:,2);
// ====================================================================

endfunction

subdemolist = demo_gateway();
clear demo_gateway; // remove demo_gateway on stack
