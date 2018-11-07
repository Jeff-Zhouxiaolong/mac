// Copyright (C) 2010 - DIGITEO - Michael Baudin
//
// This file must be used under the terms of the CeCILL.
// This source file is licensed as described in the file COPYING, which
// you should have received as part of this distribution.  The terms
// are also available at
// http://www.cecill.info/licences/Licence_CeCILL_V2-en.txt

// Updates the .xml files by deleting existing files and
// creating them again from the .sci with help_from_sci.


//
cwd = get_absolute_file_path("update_help.sce");
mprintf("Working dir = %s\n",cwd);
//
// Generate the ambiguity help
mprintf("Updating emd_toolbox/emd\n");
helpdir = fullfile(cwd,"emd");
funmat = [
"emd"
"memd"
"emd_online"
"emd_local"
  ];
macrosdir = cwd +"../../macros";
//demosdir = cwd +"../../demos";
demosdir = [];
modulename = "emd_toolbox";
helptbx_helpupdate ( funmat , helpdir , macrosdir , demosdir , modulename , %t );
//
// Generate the amplitude help
mprintf("Updating emd_toolbox/emd\n");
helpdir = fullfile(cwd,"emd");
funmat = [
  "emdc"
  "emdc_fix"
  "cemdc"
  "cemdc_fix"
  "cemdc2"
  "cemdc2_fix"
  ];
macrosdir = cwd +"../../macros/help_building";
//demosdir = cwd +"../../demos";
demosdir = [];
modulename = "emd_toolbox";
helptbx_helpupdate ( funmat , helpdir , macrosdir , demosdir , modulename , %t );
//
// Generate the datasets help
mprintf("Updating emd/help_functions\n");
helpdir = fullfile(cwd,"help_functions");
funmat = [
"emd_zero_crossings"
"emd_peaks"
"emd_local_peaks"
"emd_io"
"boundary_conditions_emd"
"emd_dirstretch"
"emd_rmtag"
"emd_hastag"
"emd_findtag"
"emd_addtag"
  ];
macrosdir = cwd +"../../macros";
//demosdir = cwd +"../../demos";
demosdir = [];
modulename = "emd_toolbox";
helptbx_helpupdate ( funmat , helpdir , macrosdir , demosdir , modulename , %t );
//
// Generate the datasets help
mprintf("Updating emd/visualization\n");
helpdir = fullfile(cwd,"visualization");
funmat = [
"disp_hhs"
"emd_visu"
"toimage"
"hhspectrum"
"cemd_visu"
"cemd_disp"
"plotc"
"plot3c"
"cenvelope"
"emd_reconstruct"
"emd_instfreq"
  ];
macrosdir = cwd +"../../macros";
//demosdir = cwd +"../../demos";
demosdir = [];
modulename = "emd_toolbox";
helptbx_helpupdate ( funmat , helpdir , macrosdir , demosdir , modulename , %t );
//
