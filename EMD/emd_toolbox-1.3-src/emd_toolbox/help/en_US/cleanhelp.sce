// ====================================================================
// 2010 - DIGITEO - Michael Baudin
// This file is released into the public domain
// ====================================================================
libpath = get_absolute_file_path('cleanhelp.sce');

mdelete(fullfile(libpath,"master_help.xml"));
rmdir(fullfile(libpath,"scilab_en_US_help"),"s");


clear libpath;

// ====================================================================
