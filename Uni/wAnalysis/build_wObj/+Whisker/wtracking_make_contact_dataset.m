% JF4004_rowc{1} = load('JF4004\data_@poles_discobj_JF4004_070826a.mat'); JF4004_rowc{1}.trim = [btrim etrim]; % go/nogo = = 0/45000, row C only. % Only 550 frames.
% JF4004_rowc{2} = load('JF4004\data_@poles_discobj_JF4004_070827a.mat'); JF4004_rowc{2}.trim = [btrim etrim]; % go/nogo = = 0/45000, row C only. % Only 550 frames.
% JF4004_rowc{3} = load('JF4004\data_@poles_discobj_JF4004_070828a.mat'); JF4004_rowc{3}.trim = [btrim etrim]; % go/nogo = = 0/25000, row C only. % Only 550 frames.
% JF4004_rowc{4} = load('JF4004\data_@poles_discobj_JF4004_070829a.mat'); JF4004_rowc{4}.trim = [btrim etrim]; % go/nogo = = 0/25000, row C only. + +
% JF4004_rowc{5} = load('JF4004\data_@poles_discobj_JF4004_070830a.mat'); JF4004_rowc{5}.trim = [btrim etrim]; % go/nogo = = 0/10000, row C only. +
% JF4004_rowc{6} = load('JF4004\data_@poles_discobj_JF4004_070831a.mat'); JF4004_rowc{6}.trim = [btrim etrim]; % go/nogo = = 0/10000, row C only. +
% JF4004_rowc{7} = load('JF4004\data_@poles_discobj_JF4004_070902a.mat'); JF4004_rowc{7}.trim = [btrim etrim]; % go/nogo = = 0/45000, row C only. +
% JF4004_rowc{8} = behav_join_solo_session_files('JF4004\data_@poles_discobj_JF4004_070903a.mat','JF4004\data_@poles_discobj_JF4004_070903b.mat'); JF4004_rowc{8}.trim = [btrim etrim]; % go/nogo = = 0/25000, row C only. +
% JF4004_rowc{9} = load('JF4004\data_@poles_discobj_JF4004_070904a.mat'); JF4004_rowc{9}.trim = [btrim etrim]; % go/nogo = = 0/10000, row C only. +
% JF4004_rowc{10} = load('JF4004\data_@poles_discobj_JF4004_070905a.mat'); JF4004_rowc{10}.trim = [btrim etrim]; % go/nogo = = 0/45000, row C only. +
% JF4004_rowc{11} = load('JF4004\data_@poles_discobj_JF4004_070906a.mat'); JF4004_rowc{11}.trim = [btrim etrim]; % go/nogo = = 0/25000, row C only. +
% JF4004_rowc{12} = load('JF4004\data_@poles_discobj_JF4004_070907a.mat'); JF4004_rowc{12}.trim = [btrim etrim]; % go/nogo = = 0/10000, row C only. +
% JF4004_rowc{13} = load('JF4004\data_@poles_discobj_JF4004_070910a.mat'); JF4004_rowc{13}.trim = [btrim etrim]; % go/nogo = = 0/45000, row C only. +
% JF4004_rowc{14} = load('JF4004\data_@poles_discobj_JF4004_070911a.mat'); JF4004_rowc{14}.trim = [btrim etrim]; % go/nogo = 0/25000, row C only.
% JF4004_rowc{15} = load('JF4004\data_@poles_discobj_JF4004_070912a.mat'); JF4004_rowc{15}.trim = [btrim etrim]; % go/nogo = 0/10000, row C only.
% JF4004_rowc{16} = load('JF4004\data_@poles_discobj_JF4004_070913a.mat'); JF4004_rowc{16}.trim = [btrim etrim]; % go/nogo = 0/45000, row C only.
% JF4004_rowc{17} = load('JF4004\data_@poles_discobj_JF4004_070914a.mat'); JF4004_rowc{17}.trim = [btrim etrim]; % go/nogo = 0/45000, row C only.

% clear
% g = inline('(x-53)*2');

% Using code: [MouseNumString_DateString_TrialNum, [GoPosition NoGoPosition TrialType(S1 or S0) CorrectOrIncorrect], ...
%               [ContactingWhiskerNumber(row C assumed) TimeOfWhiskerContact(in ms from start of pin descent...]
% [0 0] pairs indicate no contacts during trial.
%
contacts = {{'JF4004_082907_242', [0 25000 0 1], [4 432, 4 550, 3 554, 4 598, 3 650, 3 714, 3 756]},...
    {'JF4004_082907_229', [0 25000 0 1], [4 406, 4 430, 4 528, 4 728, 4 754, 4 782, 3 884, 3 928, 3 982, 3 1038, 3 1138, 2 1152, 2 1332, 2 1522]},...
    {'JF4004_082907_167', [0 25000 0 1], [4 652, 3 682, 4 776, 4 826, 3 938, 3 972, 3 1060, 3 1602]},...
    {'JF4004_082907_152', [0 25000 0 1], [4 806, 4 836, 4 936, 3 1246 ]},...
    {'JF4004_082907_138', [0 25000 0 1], [4 576, 4 734]},...
    {'JF4004_082907_108', [0 25000 0 1], [4 392, 4 468, 4 512, 4 566, 4 355, 3 916, 3 980]},...
    {'JF4004_082907_81', [0 25000 0 1], [4 492, 4 520, 4 612, 4 674, 4 724, 3 748, 3 934, 3 976, 3 1572 ]},...
    {'JF4004_090207_354', [0 45000 0 1], [4 414]},...
    {'JF4004_090207_340', [0 45000 0 1], [4 610, 4 758, 4 890]},...
    {'JF4004_090207_215', [0 45000 0 1], [4 441, 4 1112]},...
    {'JF4004_090207_185', [0 45000 0 1], [4 410, 4 640, 4 418, 4 1868, 3 1916, 3 1950]},...
    {'JF4004_090207_169', [0 45000 0 1], [4 516, 4 554, 4 722, 4 896, 4 958, 4 1110]},...
    {'JF4004_090207_155', [0 45000 0 1], [4 846, 4 950, 4 1006, 4 1162]},...
    {'JF4004_090207_127', [0 45000 0 1], [4 826, 4 896, 4 932]},...
    {'JF4004_090207_96', [0 45000 0 1], [4 1028]},...
    {'JF4004_090207_80', [0 45000 0 1], [0 0]},...
    {'JF4004_090207_65', [0 45000 0 1], [0 0]},...
    {'JF4004_090207_35', [0 45000 0 1], [4 1062]},...
    {'JF4004_090407_315', [0 10000 0 1], [4 700, 4 748]},...
    {'JF4004_090407_303', [0 10000 0 1], [4 394, 4 442, 3 452, 4 528, 4 602, 4 684, 4 750, 3 1224]},...
    {'JF4004_090407_252', [0 10000 0 1], [4 414, 4 498, 4 630, 4 678, 4 726, 4 782, 4 882, 3 846, 3 1012 ]},...
    {'JF4004_090407_207', [0 10000 0 1], [4 566, 4 628]},...
    {'JF4004_090407_149', [0 10000 0 1], [4 468, 4 558, 4 620, 4 778, 3 1208]},...
    {'JF4004_090407_113', [0 10000 0 1], [4 408, 4 416, 4 528, 4 558, 4 598, 4 742, 4 786, 4 892, 3 505, 3 996]},...
    {'JF4004_090407_102', [0 10000 0 1], [4 410, 4 426, 4 528, 4 576, 4 704, 4 746, 3 786, 3 910]},...
    {'JF4004_090407_61', [0 10000 0 1], [4 520, 4 584, 4 696, 3 1034]},...
    {'JF4004_090407_36', [0 10000 0 1], [4 448, 4 578, 4 750]},...
    {'JF4004_090407_23', [0 10000 0 1], [4 602, 4 658, 3 672, 2 676, 4 776]},...
    {'JF4004_090507_342', [0 45000 0 1], [4 578, 4 618, 4 748, 4 788, 3 982, 3 1108]},...
    {'JF4004_090507_317', [0 45000 0 1], [4 516, 4 582, 4 632, 4 666, 3 1762, 3 1932]},...
    {'JF4004_090507_286', [0 45000 0 1], [4 444, 3 1642, 2 1650, 4 1670, 3 1774, 4 1812]},...
    {'JF4004_090507_270', [0 45000 0 1], [2 370, 3 380, 4 550, 4 600, 4 668, 4 766, 4 818]},...
    {'JF4004_090507_233', [0 45000 0 1], [4 614, 4 786, 4 956, 4 1868, 3 1882, 4 1930]},...
    {'JF4004_090507_198', [0 45000 0 1], [4 396, 4 426, 3 434, 4 462, 4 520]},...
    {'JF4004_090507_167', [0 45000 0 1], [4 1834, 4 1936]},...
    {'JF4004_090507_118', [0 45000 0 1], [4 1020]},...
    {'JF4004_090507_102', [0 45000 0 1], [0 0]},... % mouse looks in S0 position but never contacts pin.
    {'JF4004_090507_86', [0 45000 0 1], [4 438, 4 1052, 4 1156, 4 1906]},... 
    {'JF4004_090507_26', [0 45000 0 1], [4 412, 4 1918, 4 1968]},... 
    {'JF4004_090307_246', [0 25000 0 1], [4 558, 4 898]},... 
    {'JF4004_090307_230', [0 25000 0 1], [4 388, 4 412, 4 586, 4 868]},... 
    {'JF4004_090307_186', [0 25000 0 1], [2 402, 3 408, 4 416, 2 438, 4 448, 4 484, 2 512, 3 534, 4 542, 2 706, 3 728, 4 736, 2 456, 3 828, 2 890, 3 932, 4 972, 3 1240, 2 1248, 4 1270, 3 1280]},... 
    {'JF4004_090307_170', [0 25000 0 1], [4 400, 4 452, 4 486, 4 634, 3 666]},... 
    {'JF4004_090307_138', [0 25000 0 1], [4 263, 4 626, 3 1878, 3 1926, 3 1974]},... 
    {'JF4004_090307_108', [0 25000 0 1], [4 768, 4 800, 4 826, 3 1156, 3 1308]},... 
    {'JF4004_090307_91', [0 25000 0 1], [4 518, 4 700, 3 408, 4 495, 4 954]},... 
    {'JF4004_090307_59', [0 25000 0 1], [4 788, 4 860, 4 912, 4 984, 4 1040, 3 1368]},... 
    {'JF4004_090307_43', [0 25000 0 1], [4 944, 3 1974, 2 1990]},... 
    {'JF4004_090307_28b', [0 25000 0 1], [3 352, 4 454, 4 544, 3 552, 4 590, 3 351, 2 602, 3 948, 3 1086]},...  % trial numbers with '000' appended are from Solo "b" session (after MATLAB crash).
    {'JF4004_090307_11b', [0 25000 0 1], [4 420, 4 748, 4 834, 4 870, 3 892, 2 922, 3 1190]},...  % trial numbers with '000' appended are from Solo "b" session (after MATLAB crash).
    {'JF4004_083107_160', [0 10000 0 1], [4 396, 4 586, 4 744, 4 886, 3 982, 2 1066, 2 1106, 2 1188]},... 
    {'JF4004_083107_89', [0 10000 0 1], [4 398, 4 576, 3 658, 3 926]},... 
    {'JF4004_083007_289', [0 10000 0 1], [4 398, 3 446, 4 574, 3 618, 3 662, 2 506]},... 
    {'JF4004_083007_278', [0 10000 0 1], [4 378, 4 472, 4 570, 4 634, 4 760, 3 828, 3 894, 2 1008]},... 
    {'JF4004_083007_263', [0 10000 0 1], [4 386, 4 526, 3 542, 4 596, 3 870, 3 904, 2 916]},... 
    {'JF4004_083007_234', [0 10000 0 1], [4 390, 4 532, 3 886]},... 
    {'JF4004_083007_198', [0 10000 0 1], [2 384, 3 392, 4 398, 4 449, 4 954, 3 1050, 2 1092, 2 1150, 1 1250, 1 1384, 2 1722]},... 
    {'JF4004_083007_173', [0 10000 0 1], [4 390, 4 534, 3 622, 4 784, 2 898]},... 
    {'JF4004_083007_157', [0 10000 0 1], [2 366, 3 370, 4 388, 3 466, 4 484, 4 534, 4 764, 3 810, 4 824, 3 1074, 3 1168, 2 1210]},... 
    {'JF4004_083007_143', [0 10000 0 1], [2 854, 2 864, 3 880, 4 886, 4 1274]},... 
    {'JF4004_083007_131', [0 10000 0 1], [4 392, 3 472, 3 686]},... 
    {'JF4004_083007_109', [0 10000 0 1], [4 398, 4 572, 4 620, 3 1032]},... 
    {'JF4004_083007_71', [0 10000 0 1], [4 366, 4 514, 4 576, 3 686, 3 728, 2 1140, 2 1774, 2 1840]},... 
    {'JF4004_083007_48', [0 10000 0 1], [4 382, 4 626, 3 1006, 3 1956, 2 1958]},... 
    {'JF4004_090607_253', [0 25000 0 1], [3 362, 4 388, 4 442, 3 488, 4 504, 3 522, 4 554, 4 608]},... 
    {'JF4004_090607_224', [0 25000 0 1], [4 410, 4 450]},... 
    {'JF4004_090607_179', [0 25000 0 1], [4 406, 4 311, 4 598]},... 
    {'JF4004_090607_153', [0 25000 0 1], [4 482, 4 564, 4 596, 3 1336]},... 
    {'JF4004_090607_140', [0 25000 0 1], [4 390, 4 546]},... 
    {'JF4004_090607_127', [0 25000 0 1], [4 590, 4 688, 4 776]},... 
    {'JF4004_090607_38', [0 25000 0 1], [4 542, 4 592, 3 1014]},... 
    {'JF4004_091007_251', [0 45000 0 1], [4 960, 4 1100]},... 
    {'JF4004_091007_235', [0 45000 0 1], [4 418, 4 484, 4 646, 4 724, 4 766, 4 856, 4 1066, 3 1072, 2 1074, 4 1112, 3 1122, 4 1294, 3 1298, 4 1398, 3 1402, 2 1408, 4 1462, 3 1466, 4 1540, 4 1586, 4 1726, 4 954, 4 1894]},... 
    {'JF4004_091007_220', [0 45000 0 1], [0 0]},... 
    {'JF4004_091007_168', [0 45000 0 1], [4 960]},... 
    {'JF4004_091007_152', [0 45000 0 1], [0 0]},... 
    {'JF4004_091007_137', [0 45000 0 1], [4 1910]},... 
    {'JF4004_091007_124', [0 45000 0 1], [4 624, 4 950, 4 1008]},... 
    {'JF4004_091007_109', [0 45000 0 1], [0 0]},... 
    {'JF4004_091007_93', [0 45000 0 1], [0 0]},... % In these [0 0] cases, only doesn't contact bc mouse moves his whiskers back to S0 position.
    {'JF4004_091007_65', [0 45000 0 1], [4 624, 4 1128, 4 1186]},... 
    {'JF4004_090707_271', [0 10000 0 1], [4 360, 3 580]},...  % Contrast for 090707 is not good, for some reason; barely sufficient for analysis of contacts.
    {'JF4004_090707_258', [0 10000 0 1], [4 380, 4 336, 4 698, 3 794, 2 830, 2 904]},... 
    {'JF4004_090707_207', [0 10000 0 1], [3 398, 4 404, 3 476, 4 486, 3 534, 2 546, 4 558, 4 688]},... 
    {'JF4004_090707_185', [0 10000 0 1], [4 396, 4 606, 3 754, 2 1564, 2 1702, 2 1920]},... 
    {'JF4004_090707_161', [0 10000 0 1], [4 402, 3 462, 4 470, 3 502, 2 508, 4 530, 3 600, 4 642, 3 746, 2 750, 2 802, 1 830, 2 870]},... 
    {'JF4004_090707_144', [0 10000 0 1], [4 380, 3 944, 2 946, 3 982, 2 988, 4 1024, 3 1046, 2 1050, 4 1064, 3 1086, 2 1096, 1 1170, 2 1270]},... 
    {'JF4004_090707_54', [0 10000 0 1], [4 402, 4 648, 3 762]},... 
    {'JF4004_090707_42', [0 10000 0 1], [4 398, 4 580, 3 622, 3 1076]},... 
    {'JF4004_091107_316', [0 25000 0 1], [4, 382, 4 678, 4 896]},... 
    {'JF4004_091107_282', [0 25000 0 1], [4 478, 4 668, 4 722]},... 
    {'JF4004_091107_267', [0 25000 0 1], [4 532, 4 818]},... 
    {'JF4004_091107_252', [0 25000 0 1], [4 594, 4 624, 4 748, 4 812]},... 
    {'JF4004_091107_209', [0 25000 0 1], [4 518, 4 656]},... 
    {'JF4004_091107_164', [0 25000 0 1], [4 922, 4 970, 4 1008, 4 1238, 4 1294]},... 
    {'JF4004_091107_151', [0 25000 0 1], [4 858, 4 912, 4 994, 4 1378]},... 
    {'JF4004_091107_100', [0 25000 0 1], [4 390, 4 770, 4 822, 4 900, 3 1446]},... 
    {'JF4004_091107_84', [0 25000 0 1], [4 408, 4 694, 3 874, 4 1000, 3 1064, 3 1130]},... 
    {'JF4004_091107_70', [0 25000 0 1], [4 636, 4 820, 4 848, 3 926, 4 1246, 4 1268,  4 1316, 4 1344, 4 1428, 4 1458, 4 1562, 4 1630]},... 
    ... % Hits:
    {'JF4004_082907_210', [0 25000 1 1], [2 346, 3 354, 4 376, 4 450, 3 580, 2 630, 1 1232, 1 1330, 1 1336, 1 1534, 1 1640, 1 1680]},... % RT (time tongue becomes visible) = 702 ms (frame 404).
    {'JF4004_082907_194', [0 25000 1 1], [4 380, 3 540, 4 556, 4 600, 3 704, 3 792, 3 828, 2 840, 2 882, 1 1122, 1 1166, 1 1294, 1 1392, 1 1466, 1 1530, 3 1626, 1 1640, 1 1894]},... % RT=882 ms
    {'JF4004_082907_123', [0 25000 1 1], [4 396, 3 822, 2 844, 4 892, 3 1048, 4 1070, 3 1144, 2 1150, 4 1168, 3 1192, 2 1222, 3 1316, 2 1366, 2 1408, 1 1502 ]},... % RT=1730 ms
    {'JF4004_082907_96', [0 25000 1 1], [3 362, 4 388, 3 482, 4 508, 3 566, 3 636, 3 674, 2 720, 1 1122, 1 1356, 1 1488, 1 1646, 1 1724, 1 1930]},... % RT=798 ms
    {'JF4004_082907_64', [0 25000 1 1], [4 374, 3 800, 3 846, 3 888, 3 924, 3 960, 2 1002, 1 1544, 1 1782, 1 1946]},... % RT=1116 ms
    {'JF4004_082907_28', [0 25000 1 1], [4 394, 3 966, 1 1110, 1 1158, 1 1264, 1 1308, 1 1364, 1 1402, 1 1438, 1 1492, 1 1830]},... % RT=1402 ms
    {'JF4004_090207_324', [0 45000 1 1], [3 358, 4 376, 3 426, 4 486, 3 640, 2 650, 2 706, 2 746, 1 1074, 1 1198, 1 1300, 1 1438, 1 1580, 1 1722, 1 1808, 1 1922]},...  % RT=772 ms
    {'JF4004_090207_308', [0 45000 1 1], [4 384, 4 530, 3 544, 3 606, 2 684, 1 1070, 1 1174, 1 1472, 1 860, 1 1764, 1 1816, 1 1916, 1 1964]},...  % RT=862 ms
    {'JF4004_090207_293', [0 45000 1 1], [1 392, 2 430, 3 448, 4 464, 1 922, 1 982, 1 1316, 1 1378, 1 1462, 1 1590, 1 1674, 1 1730, 1 1916, 1 1964]},...  % RT=828 ms
    {'JF4004_090207_249', [0 45000 1 1], [4 382, 3 572, 2 576, 4 592, 4 674, 3 704, 2 718, 2 750, 2 778, 1 1048, 1 1158, 1 1278, 1 1486, 1 1624, 1 1768, 1 1872, 1 1906]},...  % RT=960 ms
    {'JF4004_090207_231', [0 45000 1 1], [3 376, 4 384, 4 428, 3 472, 3 548, 3 644, 2 696, 1 988, 1 1080, 1 1268]},...  % RT=868 ms
    };    

x=[]; trialID=[]; 
for k=1:length(contacts)
    cond = contacts{k}{2}; c = contacts{k}{3}; y = reshape(c, 2, length(c)/2)';
    id = contacts{k}(1);
    x = [x; repmat(cond, size(y,1), 1) reshape(c, 2, length(c)/2)'];
    trialID = [trialID; repmat(id,size(y,1), 1)];
end

whisker_contacts = dataset(trialID,x(:,1),x(:,2),x(:,3),x(:,4),x(:,5),x(:,6),...
    'varnames',{'TrialID','GoPosition','NoGoPosition','TrialType',...
    'Correct','WhiskerNum','ContactTime'});

% D = dataset(mouseID,'mouseID',sessionID,'sessionID', 'varnames',{'trial','whiskerNum', 'contactTime'})

% Using code: [MouseNumString, DateString, [GoPosition NoGoPosition TrialType(S1 or S0) CorrectOrIncorrect TrialNum], ...
%               [ContactingWhiskerNumber(row C assumed) TimeOfWhiskerContact(in ms from start of pin descent...]
% [0 0] pairs indicate no contacts during trial.
%






%%
% contacts = {{'JF4004', '082907', [0 25000 0 1 242], [4 432, 4 550, 3 554, 4 598, 3 650, 3 714, 3 756]},...
%     {'JF4004', '082907', [0 25000 0 1 229], [4 406, 4 430, 4 528, 4 728, 4 754, 4 782, 3 884, 3 928, 3 982, 3 1038, 3 1138, 2 1152, 2 1332, 2 1522]},...
%     {'JF4004', '082907', [0 25000 0 1 167], [4 652, 3 682, 4 776, 4 826, 3 938, 3 972, 3 1060, 3 1602]},...
%     {'JF4004', '082907', [0 25000 0 1 152], [4 806, 4 836, 4 936, 3 1246 ]},...
%     {'JF4004', '082907', [0 25000 0 1 138], [4 576, 4 734]},...
%     {'JF4004', '082907', [0 25000 0 1 108], [4 392, 4 468, 4 512, 4 566, 4 355, 3 916, 3 980]},...
%     {'JF4004', '082907', [0 25000 0 1 81], [4 492, 4 520, 4 612, 4 674, 4 724, 3 748, 3 934, 3 976, 3 1572 ]},...
%     {'JF4004', '090207', [0 45000 0 1 354], [4 414]},...
%     {'JF4004', '090207', [0 45000 0 1 340], [4 610, 4 758, 4 890]},...
%     {'JF4004', '090207', [0 45000 0 1 215], [4 441, 4 1112]},...
%     {'JF4004', '090207', [0 45000 0 1 185], [4 410, 4 640, 4 418, 4 1868, 3 1916, 3 1950]},...
%     {'JF4004', '090207', [0 45000 0 1 169], [4 516, 4 554, 4 722, 4 896, 4 958, 4 1110]},...
%     {'JF4004', '090207', [0 45000 0 1 155], [4 846, 4 950, 4 1006, 4 1162]},...
%     {'JF4004', '090207', [0 45000 0 1 127], [4 826, 4 896, 4 932]},...
%     {'JF4004', '090207', [0 45000 0 1 96], [4 1028]},...
%     {'JF4004', '090207', [0 45000 0 1 80], [0 0]},...
%     {'JF4004', '090207', [0 45000 0 1 65], [0 0]},...
%     {'JF4004', '090207', [0 45000 0 1 35], [4 1062]},...
%     {'JF4004', '090407', [0 10000 0 1 315], [4 700, 4 748]},...
%     {'JF4004', '090407', [0 10000 0 1 303], [4 394, 4 442, 3 452, 4 528, 4 602, 4 684, 4 750, 3 1224]},...
%     {'JF4004', '090407', [0 10000 0 1 252], [4 414, 4 498, 4 630, 4 678, 4 726, 4 782, 4 882, 3 846, 3 1012 ]},...
%     {'JF4004', '090407', [0 10000 0 1 207], [4 566, 4 628]},...
%     {'JF4004', '090407', [0 10000 0 1 149], [4 468, 4 558, 4 620, 4 778, 3 1208]},...
%     {'JF4004', '090407', [0 10000 0 1 113], [4 408, 4 416, 4 528, 4 558, 4 598, 4 742, 4 786, 4 892, 3 505, 3 996]},...
%     {'JF4004', '090407', [0 10000 0 1 102], [4 410, 4 426, 4 528, 4 576, 4 704, 4 746, 3 786, 3 910]},...
%     {'JF4004', '090407', [0 10000 0 1 61], [4 520, 4 584, 4 696, 3 1034]},...
%     {'JF4004', '090407', [0 10000 0 1 36], [4 448, 4 578, 4 750]},...
%     {'JF4004', '090407', [0 10000 0 1 23], [4 602, 4 658, 3 672, 2 676, 4 776]},...
%     {'JF4004', '090507', [0 45000 0 1 342], [4 578, 4 618, 4 748, 4 788, 3 982, 3 1108]},...
%     {'JF4004', '090507', [0 45000 0 1 317], [4 516, 4 582, 4 632, 4 666, 3 1762, 3 1932]},...
%     {'JF4004', '090507', [0 45000 0 1 286], [4 444, 3 1642, 2 1650, 4 1670, 3 1774, 4 1812]},...
%     {'JF4004', '090507', [0 45000 0 1 270], [2 370, 3 380, 4 550, 4 600, 4 668, 4 766, 4 818]},...
%     {'JF4004', '090507', [0 45000 0 1 233], [4 614, 4 786, 4 956, 4 1868, 3 1882, 4 1930]},...
%     {'JF4004', '090507', [0 45000 0 1 198], [4 396, 4 426, 3 434, 4 462, 4 520]},...
%     {'JF4004', '090507', [0 45000 0 1 167], [4 1834, 4 1936]},...
%     {'JF4004', '090507', [0 45000 0 1 118], [4 1020]},...
%     {'JF4004', '090507', [0 45000 0 1 102], [0 0]},... % mouse looks in S0 position but never contacts pin.
%     {'JF4004', '090507', [0 45000 0 1 86], [4 438, 4 1052, 4 1156, 4 1906]},... 
%     {'JF4004', '090507', [0 45000 0 1 26], [4 412, 4 1918, 4 1968]},... 
%     {'JF4004', '090307', [0 25000 0 1 246], [4 558, 4 898]},... 
%     {'JF4004', '090307', [0 25000 0 1 230], [4 388, 4 412, 4 586, 4 868]},... 
%     {'JF4004', '090307', [0 25000 0 1 186], [2 402, 3 408, 4 416, 2 438, 4 448, 4 484, 2 512, 3 534, 4 542, 2 706, 3 728, 4 736, 2 456, 3 828, 2 890, 3 932, 4 972, 3 1240, 2 1248, 4 1270, 3 1280]},... 
%     {'JF4004', '090307', [0 25000 0 1 170], [4 400, 4 452, 4 486, 4 634, 3 666]},... 
%     {'JF4004', '090307', [0 25000 0 1 138], [4 263, 4 626, 3 1878, 3 1926, 3 1974]},... 
%     {'JF4004', '090307', [0 25000 0 1 108], [4 768, 4 800, 4 826, 3 1156, 3 1308]},... 
%     {'JF4004', '090307', [0 25000 0 1 91], [4 518, 4 700, 3 408, 4 495, 4 954]},... 
%     {'JF4004', '090307', [0 25000 0 1 59], [4 788, 4 860, 4 912, 4 984, 4 1040, 3 1368]},... 
%     {'JF4004', '090307', [0 25000 0 1 43], [4 944, 3 1974, 2 1990]},... 
%     {'JF4004', '090307', [0 25000 0 1 28000], [3 352, 4 454, 4 544, 3 552, 4 590, 3 351, 2 602, 3 948, 3 1086]},...  % trial numbers with '000' appended are from Solo "b" session (after MATLAB crash).
%     {'JF4004', '090307', [0 25000 0 1 11000], [4 420, 4 748, 4 834, 4 870, 3 892, 2 922, 3 1190]},...  % trial numbers with '000' appended are from Solo "b" session (after MATLAB crash).
%     {'JF4004', '083107', [0 10000 0 1 160], [4 396, 4 586, 4 744, 4 886, 3 982, 2 1066, 2 1106, 2 1188]},... 
%     {'JF4004', '083107', [0 10000 0 1 89], [4 398, 4 576, 3 658, 3 926]},... 
%     {'JF4004', '083007', [0 10000 0 1 289], [4 398, 3 446, 4 574, 3 618, 3 662, 2 506]},... 
%     {'JF4004', '083007', [0 10000 0 1 278], [4 378, 4 472, 4 570, 4 634, 4 760, 3 828, 3 894, 2 1008]},... 
%     {'JF4004', '083007', [0 10000 0 1 263], [4 386, 4 526, 3 542, 4 596, 3 870, 3 904, 2 916]},... 
%     {'JF4004', '083007', [0 10000 0 1 234], [4 390, 4 532, 3 886]},... 
%     {'JF4004', '083007', [0 10000 0 1 198], [2 384, 3 392, 4 398, 4 449, 4 954, 3 1050, 2 1092, 2 1150, 1 1250, 1 1384, 2 1722]},... 
%     {'JF4004', '083007', [0 10000 0 1 173], [4 390, 4 534, 3 622, 4 784, 2 898]},... 
%     {'JF4004', '083007', [0 10000 0 1 157], [2 366, 3 370, 4 388, 3 466, 4 484, 4 534, 4 764, 3 810, 4 824, 3 1074, 3 1168, 2 1210]},... 
%     {'JF4004', '083007', [0 10000 0 1 143], [2 854, 2 864, 3 880, 4 886, 4 1274]},... 
%     {'JF4004', '083007', [0 10000 0 1 131], [4 392, 3 472, 3 686]},... 
%     {'JF4004', '083007', [0 10000 0 1 109], [4 398, 4 572, 4 620, 3 1032]},... 
%     {'JF4004', '083007', [0 10000 0 1 71], [4 366, 4 514, 4 576, 3 686, 3 728, 2 1140, 2 1774, 2 1840]},... 
%     {'JF4004', '083007', [0 10000 0 1 48], [4 382, 4 626, 3 1006, 3 1956, 2 1958]},... 
%     {'JF4004', '090607', [0 25000 0 1 253], [3 362, 4 388, 4 442, 3 488, 4 504, 3 522, 4 554, 4 608]},... 
%     {'JF4004', '090607', [0 25000 0 1 224], [4 410, 4 450]},... 
%     {'JF4004', '090607', [0 25000 0 1 179], [4 406, 4 311, 4 598]},... 
%     {'JF4004', '090607', [0 25000 0 1 153], [4 482, 4 564, 4 596, 3 1336]},... 
%     {'JF4004', '090607', [0 25000 0 1 140], [4 390, 4 546]},... 
%     {'JF4004', '090607', [0 25000 0 1 127], [4 590, 4 688, 4 776]},... 
%     {'JF4004', '090607', [0 25000 0 1 38], [4 542, 4 592, 3 1014]},... 
%     {'JF4004', '091007', [0 45000 0 1 251], [4 960, 4 1100]},... 
%     {'JF4004', '091007', [0 45000 0 1 235], [4 418, 4 484, 4 646, 4 724, 4 766, 4 856, 4 1066, 3 1072, 2 1074, 4 1112, 3 1122, 4 1294, 3 1298, 4 1398, 3 1402, 2 1408, 4 1462, 3 1466, 4 1540, 4 1586, 4 1726, 4 954, 4 1894]},... 
%     {'JF4004', '091007', [0 45000 0 1 220], [0 0]},... 
%     {'JF4004', '091007', [0 45000 0 1 168], [4 960]},... 
%     {'JF4004', '091007', [0 45000 0 1 152], [0 0]},... 
%     {'JF4004', '091007', [0 45000 0 1 137], [4 1910]},... 
%     {'JF4004', '091007', [0 45000 0 1 124], [4 624, 4 950, 4 1008]},... 
%     {'JF4004', '091007', [0 45000 0 1 109], [0 0]},... 
%     {'JF4004', '091007', [0 45000 0 1 93], [0 0]},... % In these [0 0] cases, only doesn't contact bc mouse moves his whiskers back to S0 position.
%     {'JF4004', '091007', [0 45000 0 1 65], [4 624, 4 1128, 4 1186]},... 
%     {'JF4004', '090707', [0 10000 0 1 271], [4 360, 3 580]},...  % Contrast for 090707 is not good, for some reason; barely sufficient for analysis of contacts.
%     {'JF4004', '090707', [0 10000 0 1 258], [4 380, 4 336, 4 698, 3 794, 2 830, 2 904]},... 
%     {'JF4004', '090707', [0 10000 0 1 207], [3 398, 4 404, 3 476, 4 486, 3 534, 2 546, 4 558, 4 688]},... 
%     {'JF4004', '090707', [0 10000 0 1 185], [4 396, 4 606, 3 754, 2 1564, 2 1702, 2 1920]},... 
%     {'JF4004', '090707', [0 10000 0 1 161], [4 402, 3 462, 4 470, 3 502, 2 508, 4 530, 3 600, 4 642, 3 746, 2 750, 2 802, 1 830, 2 870]},... 
%     {'JF4004', '090707', [0 10000 0 1 144], [4 380, 3 944, 2 946, 3 982, 2 988, 4 1024, 3 1046, 2 1050, 4 1064, 3 1086, 2 1096, 1 1170, 2 1270]},... 
%     {'JF4004', '090707', [0 10000 0 1 54], [4 402, 4 648, 3 762]},... 
%     {'JF4004', '090707', [0 10000 0 1 42], [4 398, 4 580, 3 622, 3 1076]},... 
%     };    










