function optionsTree = getDefaultOptionsTreeDaniel_dec
optionsTree.zscore             = 1;
optionsTree.featuresAll        = 0;
optionsTree.featuresCross      = 0;
optionsTree.nFold              = 5;
optionsTree.storeTree          = 0;
optionsTree.verbose            = 1;
optionsTree.minLeafClassifier  = 1;
optionsTree.minLeafRegression  = 5;
optionsTree.N_bagging          = 128;
optionsTree.minTrials          = 4;
optionsTree.save               = 'allTree';
optionsTree.path               = 'F:\Diego\Matlab\temp\Daniel\Results\decoderRF_deconv\';
optionsTree.rootPath             = optionsTree.path;
optionsTree.pathPopulation       = 'population\';
optionsTree.pathTrees            = 'onlyTrees\';
optionsTree.pathSingleNeuron     = 'singleNeurons\';
optionsTree.baseFilename         = '';
optionsTree.decoder.calciumShifts = [-2 -1 0 1 2]; 
warning ('off');
mkdir (optionsTree.rootPath,optionsTree.pathPopulation);
mkdir (optionsTree.rootPath,optionsTree.pathSingleNeuron);
mkdir (sprintf('%s%s',optionsTree.rootPath,optionsTree.pathPopulation),optionsTree.pathTrees);
mkdir (sprintf('%s%s',optionsTree.rootPath,optionsTree.pathSingleNeuron),optionsTree.pathTrees);
warning ('on');