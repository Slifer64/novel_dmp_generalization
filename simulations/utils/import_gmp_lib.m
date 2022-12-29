
gmp_lib_path = strrep(mfilename('fullpath'), 'import_gmp_lib','');

addpath([gmp_lib_path '/gmp_lib/']);
addpath([gmp_lib_path '/gmp_lib/GMP/']);
addpath([gmp_lib_path '/gmp_lib/GMPo/']);
addpath([gmp_lib_path '/gmp_lib/TrajScale/']);
addpath([gmp_lib_path '/gmp_lib/CanonicalSystem/']);

addpath([gmp_lib_path '/gmp_lib/deps/QP_Goldfarb_Idnani/']);
addpath([gmp_lib_path '/gmp_lib/deps/OSQP/']);
