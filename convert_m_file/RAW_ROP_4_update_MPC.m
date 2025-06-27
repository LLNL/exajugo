%==========================================================================
% Generate ExaJuGo Files with Time-Varying Loads and Optional .m Export 
% from MATPOWER Case
%==========================================================================

% List of JSON file(s) containing time-dependent load data
filenames = {"load_CATS_2018-04-04.json", "load_CATS_2020-04-04.json"};
% filenames = "load_CATS_2018-04-04.json";  % Alternative single-file format

% Directory where load JSON files are stored
json_dir = './seasonal_data/load_data';

% Hour of interest (0–23, following 24-hour format from midnight to 11pm)
hour = 5; % This coresponse to 5am PST

% Convert PST to UCT
UTC_hour = mod(hour - 17, 24);

% MATPOWER case file with modified California system
Cal_filename = "CaliforniaTestSystem_fixed_imports_AM.m";

% Directory to save the output files
example_dir = './../examples';            % directory to save output files
casefile_dir = "Mod_Cal_System";          % case directory
output_dir = fullfile(example_dir, casefile_dir);

% Flag to control saving of the original case with no load modification
save_original = true;

% Create the output directory if it does not exist
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

% boolian to know if we want to save the modified mpc as m file
make_m_file = true;
if make_m_file
    new_mpc_dir = fullfile(output_dir, "Updated_MPC");
    if ~exist(new_mpc_dir, 'dir')
        mkdir(new_mpc_dir);
    end
end

% Load the MATPOWER case (note: changes working directory temporarily)
cd seasonal_data/updated_California_system/
mpc = loadcase(Cal_filename);
cd ../..

% Save the base MATPOWER case with no load modification (optional)
if save_original
    store_files(output_dir, mpc, '')
end

% Ensure filenames is a cell array of strings
if isstring(filenames) || ischar(filenames)
    filenames = {filenames};  % Corrected variable name from `filename` to `filenames`
end

% Loop through each JSON load file
for i = 1:length(filenames)
    fname = filenames{i};  % Current file name

    % Full path to the load JSON file
    load_path = fullfile(json_dir, fname);

    % Check if the JSON file exists
    if ~isfile(load_path)
        fprintf('The file "%s" was not found in the directory "%s".\n', fname, json_dir);
        continue;  % Skip to the next file
    end

    % Extract date string from the filename (format: YYYY-MM-DD)
    dateStr = regexp(fname, '\d{4}-\d{2}-\d{2}', 'match');
    if isempty(dateStr)
        warning('Could not extract date from filename: %s', fname);
        continue;
    end
    dateStr = dateStr{1};  % Extracted date

    % Construct a descriptive raw file name: case_YYYY-MM-DD_HH
    raw_name = "case" + "_" + dateStr + "_" + sprintf('%02d', hour);

    % Extract load IDs and corresponding demand at the specified hour
    hour_index = UTC_hour + 1;
    [load_ids, load_demand] = extract_load_data(load_path, hour_index);

    % Update the bus load values in the MATPOWER case
    for j = 1:length(load_demand)
        bus_loc = find(mpc.bus(:,1) == load_ids(j));
        mpc.bus(bus_loc, 3) = load_demand(j);  % Pd (real power demand)
    end

    % save as MPC as m file
    if make_m_file
        mpc_new_name = "CaliforniaTestSystem"+ "_" + strrep(dateStr, '-', '_') + "_" + sprintf('%02d', hour) + ".m";
        m_path = fullfile(new_mpc_dir, mpc_new_name);
        savecase(char(m_path), mpc);
    end    

    % Save the updated case using the new raw file name
    store_files(output_dir, mpc, raw_name)
end



function store_files(output_dir, mpc, diff_load)
    % STORE_FILES generates power system data files in various formats.
    % 
    % Inputs:
    %   output_dir - string or char, the directory where files will be saved
    %   mpc        - MATPOWER case struct
    %   diff_load  - (optional) string, alternative base name for output files;
    %                if empty, use default "case" and create full file set
    %
    % This function generates one or more of the following files:
    %   .raw - RAW file (PSS/E format)
    %   .rop - ROP file (PSS/E format)
    %   .con - empty file required by ExaJuGo
    %   .inl - empty file required by ExaJuGo

    % Define types of output files to be created
    types = [".raw", ".rop", ".con", ".inl"];
    
    % Default file prefix if diff_load is not provided
    parser_filename = "case";
    
    % Check if diff_load is empty (no alternative name given)
    if isempty(diff_load)
        %---------------------------------------------------------------
        % Save RAW file using default name
        save2psse(char(output_dir + "/" + parser_filename + types(1)), mpc);
        disp(['file saved: ', char(parser_filename + types(1))]);
        
        % Save ROP file using default name
        save2psse_rop(char(output_dir + "/" + parser_filename + types(2)), mpc);
        disp(['file saved: ', char(parser_filename + types(2))]);
        
        % Create empty .con file (needed by ExaJuGo even if unused)
        fclose(fopen(char(output_dir + "/" + parser_filename + types(3)), 'w'));
        disp(['empty file created: ', char(parser_filename + types(3))]);
        
        % Create empty .inl file (needed by ExaJuGo even if unused)
        fclose(fopen(char(output_dir + "/" + parser_filename + types(4)), 'w'));
        disp(['empty file created: ', char(parser_filename + types(4))]);
        %---------------------------------------------------------------
    else
        % If an alternative file prefix is given, only generate the RAW file
        save2psse(char(output_dir + "/" + diff_load + types(1)), mpc);
        disp(['file saved: ', char(diff_load + types(1))]);
    end
end


