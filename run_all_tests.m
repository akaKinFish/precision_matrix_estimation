%% Run all unit tests

% Add paths
run('startup.m');

% Get all test files
test_dir = fullfile(pwd, 'tests', 'unit');
test_files = dir(fullfile(test_dir, 'test_*.m'));

fprintf('\\n=== Running all unit tests ===\\n');
fprintf('Found %d test files\\n\\n', length(test_files));

% Run each test
passed = 0;
failed = 0;

for i = 1:length(test_files)
    [~, test_name, ~] = fileparts(test_files(i).name);
    fprintf('Running %s...\\n', test_name);
    
    try
        eval(test_name);
        passed = passed + 1;
    catch ME
        fprintf('FAILED: %s\\n', ME.message);
        failed = failed + 1;
    end
    fprintf('\\n');
end

fprintf('\\n=== Test Summary ===\\n');
fprintf('Passed: %d\\n', passed);
fprintf('Failed: %d\\n', failed);
