signal1 = load('../data/MLII/1 NSR/100m (0).mat');
signal2 = load('../data/MLII/1 NSR/100m (1).mat');
signal3 = load('../data/MLII/1 NSR/100m (2).mat');

ecg = [signal1.val signal2.val];

Fs = 360;