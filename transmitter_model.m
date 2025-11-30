clc; clear; close all;
% PRE PROCESSING DATA
file = fopen("Data/Transmitter/data.csv", "r");
sanizited_data = fopen("Data/Transmitter/sanitized_data.csv", "w");
while ~feof(file)
    line = fgetl(file);
    new_line = strrep(line, ',','.');
    new_line = strrep(new_line, ';',',');
    new_line = strcat(new_line,'\n');
    fprintf(sanizited_data, new_line);
    disp(new_line);
end
fclose(file);
fclose(sanizited_data);

data = readtable("Data/Transmitter/sanitized_data.csv",'PreserveVariableNames',true);
figure;
scatter(data.x, data.larvae_2);
grid on;
figure;
scatter(data.x, data.larvae_4);
grid on;
figure;
scatter(data.x, data.larvae_8);
grid on;