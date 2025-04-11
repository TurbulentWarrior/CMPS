%% 1. Generate a Structured Magnetic Disk Configuration

% Disk layout parameters
numTracks   = 10;         % Number of concentric tracks (radial divisions)
bitsPerTrack = 60;        % Number of bits (angular sectors) per track
Rmin        = 20;         % Inner radius of the disk (inactive region inside)
Rmax        = 100;        % Outer radius of the disk
N           = 500;        % Resolution of the simulation grid (NxN)

% Create a Cartesian grid covering the full disk
[x, y] = meshgrid(linspace(-Rmax, Rmax, N));
R     = sqrt(x.^2 + y.^2);      % Radial coordinate
Theta = atan2(y, x);            % Angular coordinate (-pi to pi)

% Define the valid (active) region: points between Rmin and Rmax
valid = (R >= Rmin) & (R <= Rmax);

% Normalize the radial coordinate (only for valid points)
Rnorm = (R - Rmin) / (Rmax - Rmin);
Rnorm(~valid) = NaN;  % Mark outside the disk with NaN

% Map each valid point to a track (radial division)
trackIdx = floor(Rnorm * numTracks) + 1;

% Normalize the angle to [0, 1) and divide into sectors (bits)
angleNorm = mod(Theta, 2*pi) / (2*pi);
bitIdx = floor(angleNorm * bitsPerTrack) + 1;

% Create a deterministic bit pattern.
patternBase = [1 0 1 1 0 0 1 1 0 1];
bitPattern = repmat(patternBase, 1, ceil(bitsPerTrack/length(patternBase)));
bitPattern = bitPattern(1:bitsPerTrack);  % truncate to exactly bitsPerTrack

% Initialize the disk configuration matrix.
diskConfig = NaN(N,N);
for t = 1:numTracks
    for b = 1:bitsPerTrack
        domainValue = 2 * bitPattern(mod(t + b, bitsPerTrack) + 1) - 1;
        mask = (trackIdx == t) & (bitIdx == b);
        diskConfig(mask) = domainValue;
    end
end

diskConfig(~valid) = NaN;

% Show the initial disk configuration (active region only)
figure;
imagesc(diskConfig);
axis equal off;
colormap(gray);
title('Initial Structured Magnetic Disk Configuration');

%% 2. Apply the Metropolis Algorithm to the Disk (Interactions Only Inside the Disk)

% Simulation parameters
J      = 1;        % Coupling constant
kb     = 1;        % Boltzmann constant
T      = 2;     % Temperature (adjust this value as desired)
nSteps = 1e4;      % Total number of attempted spin flips

% Pre-calculate the list of valid indices (sites with spin ±1)
[valid_i, valid_j] = find(~isnan(diskConfig));
numValid = numel(valid_i);

% Set up the VideoWriter to save the simulation as a video
videoFileName = 'magnetic_disk_simulation.avi';
v = VideoWriter(videoFileName, 'Motion JPEG AVI');
v.FrameRate = 20;  % Adjust frame rate (frames per second)
open(v);  % Open the video file for writing

% Create a figure for the simulation animation.
figure;

% Initialize tracking for corruption rate over time
corruptionRateHistory = [];  % To store corruption rate at each step

% Store the initial disk configuration for comparison
diskConfig_initial = diskConfig;  % This is the starting point for comparison

for m=1:250
    for step = 1:nSteps
        % Choose a random valid site (i,j) from the list
        idx = randi(numValid);
        i   = valid_i(idx);
        j   = valid_j(idx);
        
        cur = diskConfig(i, j);  % Current spin at (i,j)
        
        % Calculate the sum of spins on the four neighbors.
        neighborSum = 0;
        if i > 1 && ~isnan(diskConfig(i-1, j))  % Up neighbor
           neighborSum = neighborSum + diskConfig(i-1, j);
        end
        if i < N && ~isnan(diskConfig(i+1, j))  % Down neighbor
           neighborSum = neighborSum + diskConfig(i+1, j);
        end
        if j > 1 && ~isnan(diskConfig(i, j-1))  % Left neighbor
           neighborSum = neighborSum + diskConfig(i, j-1);
        end
        if j < N && ~isnan(diskConfig(i, j+1))  % Right neighbor
           neighborSum = neighborSum + diskConfig(i, j+1);
        end
        
        % Compute the change in energy if the spin is flipped.
        H     = -J * cur * neighborSum;
        % Metropolis acceptance criterion:
        if -2*H < 0 || rand() < exp(H/(kb*T))
            diskConfig(i, j) = -cur;
        end
    end

    % Optionally update 
    if mod(m, 2) == 0
        imagesc(diskConfig);
        axis equal off;
        colormap(gray);
        drawnow;
        pause(0.2)

        % Capture the current frame for the video
        frame = getframe(gcf);
        writeVideo(v, frame);  % Add the current frame to the video

        % Calculate the number of flipped bits (change from initial config)
        flippedBits = sum(diskConfig ~= diskConfig_initial & ~isnan(diskConfig));

        % Calculate total active bits (valid bits inside the disk)
        totalActiveBits = sum(~isnan(diskConfig));

        % Calculate the corruption rate (percentage of flipped bits)
        corruptionRate = (flippedBits / totalActiveBits) * 100;

        % Track the corruption rate over time
        corruptionRateHistory = [corruptionRateHistory; corruptionRate];
    end
end

% Close the video writer
close(v);

% Show final configuration
figure;
imagesc(diskConfig);
axis equal off;
colormap(gray);
title('Final Disk Configuration after Metropolis Simulation');

%% 3. Plot the Corruption Rate Over Time

figure;
plot(corruptionRateHistory, 'LineWidth', 2);
xlabel('Time Step');
ylabel('Corruption Rate (%)');
title('Corruption Rate Over Time');
grid on;
