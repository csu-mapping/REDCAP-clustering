function test()
data = load('data.txt');
connect = load('data1.txt');

warning('off');

k = size(data,1) - 1;
% First-slk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'slk','first');
xlswrite('RegionIDX.xls',RegionIDX, 'FisrtSLK');
xlswrite('Z.xls',Z, 'FisrtSLK');

% First-alk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'alk','first');
xlswrite('RegionIDX.xls',RegionIDX, 'FisrtALK');
xlswrite('Z.xls',Z, 'FisrtALK');

% First-clk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'clk','first');
xlswrite('RegionIDX.xls',RegionIDX, 'FisrtCLK');
xlswrite('Z.xls',Z, 'FisrtCLK');

% Full-slk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'slk','full');
xlswrite('RegionIDX.xls',RegionIDX, 'FullSLK');
xlswrite('Z.xls',Z, 'FullSLK');

% Full-alk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'alk','full');
xlswrite('RegionIDX.xls',RegionIDX, 'FullALK');
xlswrite('Z.xls',Z, 'FullALK');

% Full-clk
[RegionIDX,Z,LinkageIDX] = redcap(data,connect,k,'clk','full');
xlswrite('RegionIDX.xls',RegionIDX, 'FullCLK');
xlswrite('Z.xls',Z, 'FullCLK');
%save('full-clk.mat', 'RegionIDX','Z');


clc;

end





