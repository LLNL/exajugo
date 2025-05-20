function [load_ids, load_demand] = extract_load_data(filename, hour)
    % EXTRACT_LOAD_DATA extracts load IDs and corresponding demand values
    % from a JSON file at a specified hour.
    %
    % Inputs:
    %   filename - string, path to the JSON file containing load data
    %   hour     - integer, the hour index to extract load demand
    %
    % Outputs:
    %   load_ids     - vector of numeric load bus IDs
    %   load_demand  - vector of load demand values (in MW) at the given hour

    % Read the file contents into a character vector
    jsonText = fileread(filename);
    
    % Decode the JSON text into a MATLAB data structure (a struct)
    data = jsondecode(jsonText);
    
    % Get all field names (load IDs in the format 'x#') from the load_MW struct
    keys = fieldnames(data.load_MW);  
    n = length(keys);  % Number of loads
    
    % Preallocate arrays for efficiency
    load_ids = zeros(n, 1);
    load_demand = zeros(n, 1);  
    
    % Loop through each load entry
    for i = 1:n
        keyStr = keys{i};                       % Get the current field name (e.g., 'x123')
        load_demand(i) = data.load_MW.(keyStr)(hour);  % Extract load value at the specified hour
        numStr = extractAfter(keyStr, 1);       % Remove 'x' prefix to isolate the number
        load_ids(i) = str2double(numStr);       % Convert the string number to a numeric ID
    end
end
