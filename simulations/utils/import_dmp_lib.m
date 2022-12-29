
dmp_lib_path = strrep(mfilename('fullpath'), 'import_dmp_lib','');

dmp_lib_path = [dmp_lib_path '/dmp_lib/'];

addpath(dmp_lib_path);
addpath([dmp_lib_path '/DMP/']);
addpath([dmp_lib_path '/QPMP/']);
addpath([dmp_lib_path '/CanonicalClock/']);
addpath([dmp_lib_path '/GatingFunction/']);
addpath([dmp_lib_path '/trainMethods/']);
addpath([dmp_lib_path '/enums/']);
addpath([dmp_lib_path '/math/']);
