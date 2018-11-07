// This file is released under the 3-clause BSD license. See COPYING-BSD.

function builder_gw_c()
CURRENT_PATH = strsubst(get_absolute_file_path("builder_gateway_c.sce"), "\", "/");

if getos() <> "Windows" then
  INCLUDES_PATHS = "-I" + CURRENT_PATH + " -DC99_OK ";;
else
  INCLUDES_PATHS = "";
end

//INCLUDES_PATHS = INCLUDES_PATHS + " -DC99_OK ";

// PutLhsVar managed by user in sci_sum and in sci_sub
// if you do not this variable, PutLhsVar is added
// in gateway generated (default mode in scilab 4.x and 5.x)
WITHOUT_AUTO_PUTLHSVAR = %T;

FUNCTIONS_GATEWAY = ["emdc","int_emdc";"emdc_fix","int_emdc_fix";"cemdc","int_cemdc";...
"cemdc_fix","int_cemdc_fix";"cemdc2","int_cemdc2";"cemdc2_fix","int_cemdc2_fix";]

FILES_GATEWAY = [ "extr.c",  "interpolation.c",  "io.c", "io_fix.c", "local_mean.c","emdc_fix.c", "emdc.c",...
"emd_complex.c", "cextr.c", "clocal_mean.c", "clocal_mean2.c","cio.c", "cio_fix.c","cemdc.c","cemdc_fix.c","cemdc2.c","cemdc2_fix.c"];

tbx_build_gateway("emd_c", FUNCTIONS_GATEWAY, FILES_GATEWAY, CURRENT_PATH, [],"",INCLUDES_PATHS);

endfunction

builder_gw_c();
clear builder_gw_c; // remove builder_gw_c on stack
