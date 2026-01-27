% TD = Read_Ndataset(filename)
% returns the Temporal Difference (TD) events from binary file for the
% N-MNIST and N-Caltech101 datasets. See garrickorchard.com\datasets for
% more info

% modified for 31 bit timestamp (up to 2147 seconds or ~30minutes)
function TD = Read_Ndataset_longTs(filename)
eventData = fopen(filename);
evtStream = fread(eventData);
fclose(eventData);

TD.x    = evtStream(1:6:end); %+1; %pixel x address, with first pixel having index 1
TD.y    = evtStream(2:6:end); %+1; %pixel y address, with first pixel having index 1
TD.p    = bitshift(evtStream(3:6:end), -7)+1; %polarity, 1 means off, 2 means on
TD.ts   = bitshift(bitand(evtStream(3:6:end), 127), 24); %time in microseconds
TD.ts   = TD.ts + bitshift(evtStream(4:6:end), 16);
TD.ts   = TD.ts + bitshift(evtStream(5:6:end), 8);
TD.ts   = TD.ts + evtStream(6:6:end);
return