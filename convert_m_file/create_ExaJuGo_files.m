%==========================================================================
% Generate Output Files (.raw, .rop, .con, .inl) for ExaJuGo from MATPOWER Case
%==========================================================================

% Define the output directory for the generated files
example_dir = './examples';            % directory to save output files
casefile_dir = 'California';           % case directory
% casefile_dir = '9bus';                 % <- alternate example: 9-bus system

output_dir = fullfile(example_dir, casefile_dir);

% Create the output directory if it does not exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% Specify the MATPOWER case filename (without extension)
filename = "CaliforniaTestSystem";       % <- name of .m case file in ./example_m_files/
% filename = "case9_mod";                % <- alternate example: 9-bus system

% Change to directory containing the MATPOWER .m file
cd ./example_m_files/

% Load the MATPOWER case
mpc = loadcase(filename + ".m");

% Return to the base directory (where the output directory exists)
cd ..

% Define types of output files to be created
types = [".raw", ".rop", ".con", ".inl"];

% Define common prefix for generated files
parser_filename = "case";

%---------------------------------------------------------------------------
% Save MATPOWER case as PSS/E RAW file
raw_file = save2psse(char(output_dir + "/" + parser_filename + types(1)), mpc);
disp(['file saved: ', char(parser_filename + types(1))]);

% Save MATPOWER case as PSS/E ROP file
rop_file = save2psse_rop(char(output_dir + "/" + parser_filename + types(2)), mpc);
disp(['file saved: ', char(parser_filename + types(2))]);

% Create empty .con and .inl files (required by ExaJuGo)
fclose(fopen(char(output_dir + "/" + parser_filename + types(3)), 'w'));
disp(['empty file created: ', char(parser_filename + types(3))]);

fclose(fopen(char(output_dir + "/" + parser_filename + types(4)), 'w'));
disp(['empty file created: ', char(parser_filename + types(4))]);
%---------------------------------------------------------------------------
