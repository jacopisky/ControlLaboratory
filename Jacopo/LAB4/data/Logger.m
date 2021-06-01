packetSpec = {'single','single','single','single','single','single','single'};

% interactive data−logging
data = serial_datalog( 'COM3', ... % serial port name
packetSpec, ... % packet specifier
'BaudRate', 38400, ... % serial port baud rate
'TxSampleTime', 0.01, ... % TX sample time
'TxStartCmd', 65, ... % TX start cmd (default)
'TxStopCmd', 66, ... % TX stop cmd (default)
'BufferSize', 10000, ... % RX buffer size
'PlotWidth', 5, ... % time−axis width (default)
'PlotRefreshRatio', 10, ... % plot refresh ratio
'PlotList', [1,2,3,4,5,6]); % plot list

t = data.time;
gamma = data.out{1};
theta = data.out{2};
dot_gamma = data.out{3};
dot_theta = data.out{4};
duty = data.out{5};
ref = data.out{6};
disturb = data.out{7};