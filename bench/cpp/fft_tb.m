% Read the file
fid = fopen('fft_tb.dbl','r');
ftlen = fread(fid, 1, 'int32');
raw = fread(fid, [2 inf], 'double');
fclose(fid);

% Convert the raw doubles into complex values
datc = raw(1,:)+j*raw(2,:);
% Reshape the matrix into one line per FFT
ndat = reshape(datc, ftlen*2, floor(length(datc)/(ftlen*2)));

truth  = ndat((ftlen+1):(2*ftlen), :);
output = ndat(1:ftlen,:);

% Create a time axis, for use in plotting if desired
tm = 0:(ftlen-1);

% Now, the data from the test is ready for inspection
k = 1; plottst; k

