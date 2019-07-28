function im = readModulationsDatabaseForAlexNet(filename)
%Read Image
im = imread(filename);

%Repeat RGB dimension
im = repmat(im,[1 1 3]);