addpath(genpath('/Users/christopherpenfold/Desktop/Thorsten/spacemonkey/code/gpml-matlab-v3.6-2015-07-07'))
addpath(genpath('/Volumes/GoogleDrive/My\ Drive/Marmoset_shared_folder/Shiny_plots_Erin/SpatialModelling/'))



%Load in the 3D scaffold
[OBJ2,section] = LoadCS6('3D');
%[OBJ2A,a2,b2] = transformCS6(OBJ1,'notall')
%Load shot locations
[D,Locations,XYZ,CellType,Shots] = LoadShots('CS6');
%Process the shots onto scaffold
[Output] = loadCS6Scaffold(D,Locations,Shots);
[Output] = MarmosetGP_CS6_v3(D,Output,'DPPA2');
[Output] = MarmosetGPInfer_CS6_v3(Output,OBJ2);
%f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceVertexCData',Output.m1,'FaceColor','interp','LineStyle','none');
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',Output.m2,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',Output.m3,'FaceColor','interp','LineStyle','none');%
axis equal
axis off
material dull 
colormap(parula);
%view([44.1000,-39.2638])
%view([41,90])
view([253.1703,-28.4480])
camlight('right')
axis equal
axis off
colorbar
%ax = gca
%ax.CLim = [min([Output.cLim1,Output.cLim3]),max([Output.cLim1,Output.cLim3])];

%view([44.1000,-39.2638])
%iew([41,90])



g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'TFAP2C');
g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'TFAP2A');
g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'VTCN1');
g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'T');
g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'MIXL1');
g1 =[g1;Output.Ytrain];
[Output] = MarmosetGP_CS6_v3(D,Output,'OTX2');
g1 =[g1;Output.Ytrain];




%Load CS6
[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set1.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set2.csv'
[newX2,newY2,newZ2,X2,Y2,Z2,W2,CellType2] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set3.csv'
[newX3,newY3,newZ3,X3,Y3,Z3,W3,CellType3] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set4.csv'
[newX4,newY4,newZ4,X4,Y4,Z4,W4,CellType4] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set5.csv'
[newX5,newY5,newZ5,X5,Y5,Z5,W5,CellType5] = ProjectData(File1,File2,ShotsCS6,10);

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set1.csv'
[newXs1,newYs1,newZs1,Xs1,Ys1,Zs1,Ws1,CellTypes1] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set2.csv'
[newXs2,newYs2,newZs2,Xs2,Ys2,Zs2,Ws2,CellTypes2] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set3.csv'
[newXs3,newYs3,newZs3,Xs3,Ys3,Zs3,Ws3,CellTypes3] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set4.csv'
[newXs4,newYs4,newZs4,Xs4,Ys4,Zs4,Ws4,CellTypes4] = ProjectData(File1,File2,ShotsCS6,10);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set5.csv'
[newXs5,newYs5,newZs5,Xs5,Ys5,Zs5,Ws5,CellTypes5] = ProjectData(File1,File2,ShotsCS6,10);

S1 = sqrt( (newX1-X1).^2 + (newY1-Y1).^2 + (newZ1-Z1).^2 );
S2 = sqrt( (newX2-X2).^2 + (newY2-Y2).^2 + (newZ2-Z2).^2 );
S3 = sqrt( (newX3-X3).^2 + (newY3-Y3).^2 + (newZ3-Z3).^2 );
S4 = sqrt( (newX4-X4).^2 + (newY4-Y4).^2 + (newZ4-Z4).^2 );
S5 = sqrt( (newX5-X5).^2 + (newY5-Y5).^2 + (newZ5-Z5).^2 );

S1_ = sqrt( (newXs1-Xs1).^2 + (newYs1-Ys1).^2 + (newZs1-Zs1).^2 );
S2_ = sqrt( (newXs2-Xs2).^2 + (newYs2-Ys2).^2 + (newZs2-Zs2).^2 );
S3_ = sqrt( (newXs3-Xs3).^2 + (newYs3-Ys3).^2 + (newZs3-Zs3).^2 );
S4_ = sqrt( (newXs4-Xs4).^2 + (newYs4-Ys4).^2 + (newZs4-Zs4).^2 );
S5_ = sqrt( (newXs5-Xs5).^2 + (newYs5-Ys5).^2 + (newZs5-Zs5).^2 );

RME = [S1,S1_,S2,S2_,S3,S3_,S4,S4_,S5,S5_];

csvwrite('MSE_N=10.csv',RME)
T = cell2table([CellType1,CellType2,CellType3,CellType4,CellType5,CellTypes1,CellTypes2,CellTypes3,CellTypes4,CellTypes5])
writetable(T,'CellType_N=10.csv')


[OBJ2,section] = LoadCS6('3D');
[D,Locations,XYZ,CellType,ShotsCS6] = LoadShots('CS6');

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set1.csv'
[newX1,newY1,newZ1,X1,Y1,Z1,W1,CellType1,RegCellType1] = ProjectData(File1,File2,ShotsCS6,5);

[newX1c,newY1c,newZ1c] = CorrectProjectData(newX1,newY1,newZ1,W1,CellType1,RegCellType1,OBJ2);

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set2.csv'
[newX2,newY2,newZ2,X2,Y2,Z2,W2,CellType2] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set3.csv'
[newX3,newY3,newZ3,X3,Y3,Z3,W3,CellType3] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set4.csv'
[newX4,newY4,newZ4,X4,Y4,Z4,W4,CellType4] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1_set5.csv'
[newX5,newY5,newZ5,X5,Y5,Z5,W5,CellType5] = ProjectData(File1,File2,ShotsCS6,5);

File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set1.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set1.csv'
[newXs1,newYs1,newZs1,Xs1,Ys1,Zs1,Ws1,CellTypes1] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set2.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set2.csv'
[newXs2,newYs2,newZs2,Xs2,Ys2,Zs2,Ws2,CellTypes2] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set3.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set3.csv'
[newXs3,newYs3,newZs3,Xs3,Ys3,Zs3,Ws3,CellTypes3] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set4.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set4.csv'
[newXs4,newYs4,newZs4,Xs4,Ys4,Zs4,Ws4,CellTypes4] = ProjectData(File1,File2,ShotsCS6,5);
File1 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/S1_set5.csv'
File2 = '/Users/christopherpenfold/Desktop/Thorsten/FINAL/MappingwMissingData/C1perm_set5.csv'
[newXs5,newYs5,newZs5,Xs5,Ys5,Zs5,Ws5,CellTypes5] = ProjectData(File1,File2,ShotsCS6,5);

S1 = sqrt( (newX1-X1).^2 + (newY1-Y1).^2 + (newZ1-Z1).^2 );
S2 = sqrt( (newX2-X2).^2 + (newY2-Y2).^2 + (newZ2-Z2).^2 );
S3 = sqrt( (newX3-X3).^2 + (newY3-Y3).^2 + (newZ3-Z3).^2 );
S4 = sqrt( (newX4-X4).^2 + (newY4-Y4).^2 + (newZ4-Z4).^2 );
S5 = sqrt( (newX5-X5).^2 + (newY5-Y5).^2 + (newZ5-Z5).^2 );

S1_ = sqrt( (newXs1-Xs1).^2 + (newYs1-Ys1).^2 + (newZs1-Zs1).^2 );
S2_ = sqrt( (newXs2-Xs2).^2 + (newYs2-Ys2).^2 + (newZs2-Zs2).^2 );
S3_ = sqrt( (newXs3-Xs3).^2 + (newYs3-Ys3).^2 + (newZs3-Zs3).^2 );
S4_ = sqrt( (newXs4-Xs4).^2 + (newYs4-Ys4).^2 + (newZs4-Zs4).^2 );
S5_ = sqrt( (newXs5-Xs5).^2 + (newYs5-Ys5).^2 + (newZs5-Zs5).^2 );

RME = [S1,S1_,S2,S2_,S3,S3_,S4,S4_,S5,S5_];

csvwrite('MSE_N=5.csv',RME)
T = cell2table([CellType1,CellType2,CellType3,CellType4,CellType5,CellTypes1,CellTypes2,CellTypes3,CellTypes4,CellTypes5])
writetable(T,'CellType_N=5.csv')


ind1_1 = find(strcmp(CellType1,'EmDisc')==1);
ind1_2 = find(strcmp(CellType2,'EmDisc')==1);
ind1_3 = find(strcmp(CellType3,'EmDisc')==1);
ind1_4 = find(strcmp(CellType4,'EmDisc')==1);
ind1_5 = find(strcmp(CellType5,'EmDisc')==1);

inds1_1 = find(strcmp(CellTypes1,'EmDisc')==1);
inds1_2 = find(strcmp(CellTypes2,'EmDisc')==1);
inds1_3 = find(strcmp(CellTypes3,'EmDisc')==1);
inds1_4 = find(strcmp(CellTypes4,'EmDisc')==1);
inds1_5 = find(strcmp(CellTypes5,'EmDisc')==1);

ind2_1 = find(strcmp(CellType1,'VE')==1);
ind2_2 = find(strcmp(CellType2,'VE')==1);
ind2_3 = find(strcmp(CellType3,'VE')==1);
ind2_4 = find(strcmp(CellType4,'VE')==1);
ind2_5 = find(strcmp(CellType5,'VE')==1);

inds2_1 = find(strcmp(CellTypes1,'VE')==1);
inds2_2 = find(strcmp(CellTypes2,'VE')==1);
inds2_3 = find(strcmp(CellTypes3,'VE')==1);
inds2_4 = find(strcmp(CellTypes4,'VE')==1);
inds2_5 = find(strcmp(CellTypes5,'VE')==1);


ind3_1 = find(strcmp(CellType1,'Am')==1);
ind3_2 = find(strcmp(CellType2,'Am')==1);
ind3_3 = find(strcmp(CellType3,'Am')==1);
ind3_4 = find(strcmp(CellType4,'Am')==1);
ind3_5 = find(strcmp(CellType5,'Am')==1);

inds3_1 = find(strcmp(CellTypes1,'Am')==1);
inds3_2 = find(strcmp(CellTypes2,'Am')==1);
inds3_3 = find(strcmp(CellTypes3,'Am')==1);
inds3_4 = find(strcmp(CellTypes4,'Am')==1);
inds3_5 = find(strcmp(CellTypes5,'Am')==1);


ind4_1 = find(strcmp(CellType1,'Tb')==1);
ind4_2 = find(strcmp(CellType2,'Tb')==1);
ind4_3 = find(strcmp(CellType3,'Tb')==1);
ind4_4 = find(strcmp(CellType4,'Tb')==1);
ind4_5 = find(strcmp(CellType5,'Tb')==1);

inds4_1 = find(strcmp(CellTypes1,'Tb')==1);
inds4_2 = find(strcmp(CellTypes2,'Tb')==1);
inds4_3 = find(strcmp(CellTypes3,'Tb')==1);
inds4_4 = find(strcmp(CellTypes4,'Tb')==1);
inds4_5 = find(strcmp(CellTypes5,'Tb')==1);

%EmDisc
h=figure(1)
subplot(2,4,1);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])
camlight('left')
material dull 
hold on
scatter3(newX1(ind1_1,1),newY1(ind1_1,1),newZ1(ind1_1,1),W1(ind1_1)*60,'k','filled')
scatter3(newX2(ind1_2,1),newY1(ind1_2,1),newZ1(ind1_2,1),W2(ind1_2)*60,'k','filled')
scatter3(newX3(ind1_3,1),newY1(ind1_3,1),newZ1(ind1_3,1),W3(ind1_3)*60,'k','filled')
scatter3(newX4(ind1_4,1),newY1(ind1_4,1),newZ1(ind1_4,1),W4(ind1_4)*60,'k','filled')
scatter3(newX5(ind1_5,1),newY1(ind1_5,1),newZ1(ind1_5,1),W5(ind1_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,5);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);


subplot(2,4,2);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind2_1,1),newY1(ind2_1,1),newZ1(ind2_1,1),W1(ind2_1)*60,'k','filled')
scatter3(newX2(ind2_2,1),newY1(ind2_2,1),newZ1(ind2_2,1),W2(ind2_2)*60,'k','filled')
scatter3(newX3(ind2_3,1),newY1(ind2_3,1),newZ1(ind2_3,1),W3(ind2_3)*60,'k','filled')
scatter3(newX4(ind2_4,1),newY1(ind2_4,1),newZ1(ind2_4,1),W4(ind2_4)*60,'k','filled')
scatter3(newX5(ind2_5,1),newY1(ind2_5,1),newZ1(ind2_5,1),W5(ind2_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,6);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds2_1,1),newYs1(inds2_1,1),newZs1(inds2_1,1),Ws1(inds2_1)*60,'r','filled')
scatter3(newXs2(inds2_2,1),newYs1(inds2_2,1),newZs1(inds2_2,1),Ws2(inds2_2)*60,'r','filled')
scatter3(newXs3(inds2_3,1),newYs1(inds2_3,1),newZs1(inds2_3,1),Ws3(inds2_3)*60,'r','filled')
scatter3(newXs4(inds2_4,1),newYs1(inds2_4,1),newZs1(inds2_4,1),Ws4(inds2_4)*60,'r','filled')
scatter3(newXs5(inds2_5,1),newYs1(inds2_5,1),newZs1(inds2_5,1),Ws5(inds2_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);








subplot(2,4,3);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind3_1,1),newY1(ind3_1,1),newZ1(ind3_1,1),W1(ind3_1)*60,'k','filled')
scatter3(newX2(ind3_2,1),newY1(ind3_2,1),newZ1(ind3_2,1),W2(ind3_2)*60,'k','filled')
scatter3(newX3(ind3_3,1),newY1(ind3_3,1),newZ1(ind3_3,1),W3(ind3_3)*60,'k','filled')
scatter3(newX4(ind3_4,1),newY1(ind3_4,1),newZ1(ind3_4,1),W4(ind3_4)*60,'k','filled')
scatter3(newX5(ind3_5,1),newY1(ind3_5,1),newZ1(ind3_5,1),W5(ind3_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,7);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds3_1,1),newYs1(inds3_1,1),newZs1(inds3_1,1),Ws1(inds3_1)*60,'r','filled')
scatter3(newXs2(inds3_2,1),newYs1(inds3_2,1),newZs1(inds3_2,1),Ws2(inds3_2)*60,'r','filled')
scatter3(newXs3(inds3_3,1),newYs1(inds3_3,1),newZs1(inds3_3,1),Ws3(inds3_3)*60,'r','filled')
scatter3(newXs4(inds3_4,1),newYs1(inds3_4,1),newZs1(inds3_4,1),Ws4(inds3_4)*60,'r','filled')
scatter3(newXs5(inds3_5,1),newYs1(inds3_5,1),newZs1(inds3_5,1),Ws5(inds3_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,4);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
%view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newX1(ind4_1,1),newY1(ind4_1,1),newZ1(ind4_1,1),W1(ind4_1)*60,'k','filled')
scatter3(newX2(ind4_2,1),newY1(ind4_2,1),newZ1(ind4_2,1),W2(ind4_2)*60,'k','filled')
scatter3(newX3(ind4_3,1),newY1(ind4_3,1),newZ1(ind4_3,1),W3(ind4_3)*60,'k','filled')
scatter3(newX4(ind4_4,1),newY1(ind4_4,1),newZ1(ind4_4,1),W4(ind4_4)*60,'k','filled')
scatter3(newX5(ind4_5,1),newY1(ind4_5,1),newZ1(ind4_5,1),W5(ind4_5)*60,'k','filled')
%scatter3(newXs1(inds1_1,1),newYs1(inds1_1,1),newZs1(inds1_1,1),Ws1(inds1_1)*60,'r','filled')
%scatter3(newXs2(inds1_2,1),newYs1(inds1_2,1),newZs1(inds1_2,1),Ws2(inds1_2)*60,'r','filled')
%scatter3(newXs3(inds1_3,1),newYs1(inds1_3,1),newZs1(inds1_3,1),Ws3(inds1_3)*60,'r','filled')
%scatter3(newXs4(inds1_4,1),newYs1(inds1_4,1),newZs1(inds1_4,1),Ws4(inds1_4)*60,'r','filled')
%scatter3(newXs5(inds1_5,1),newYs1(inds1_5,1),newZs1(inds1_5,1),Ws5(inds1_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);

subplot(2,4,8);
f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
view([41,90])

%view([c1,d1])
camlight('left')
material dull 
hold on
scatter3(newXs1(inds4_1,1),newYs1(inds4_1,1),newZs1(inds4_1,1),Ws1(inds4_1)*60,'r','filled')
scatter3(newXs2(inds4_2,1),newYs1(inds4_2,1),newZs1(inds4_2,1),Ws2(inds4_2)*60,'r','filled')
scatter3(newXs3(inds4_3,1),newYs1(inds4_3,1),newZs1(inds4_3,1),Ws3(inds4_3)*60,'r','filled')
scatter3(newXs4(inds4_4,1),newYs1(inds4_4,1),newZs1(inds4_4,1),Ws4(inds4_4)*60,'r','filled')
scatter3(newXs5(inds4_5,1),newYs1(inds4_5,1),newZs1(inds4_5,1),Ws5(inds4_5)*60,'r','filled')
%scatter3(newX(ind8(i),1),newY(ind8(i),1),newZ(ind8(i),1),W(ind8(i))*60,'r','filled')
%scatter3(newX(ind7(i),1),newY(ind7(i),1),newZ(ind7(i),1),W(ind7(i))*60,'b','filled')
colormap(parula);


%EmDisc
%h=figure(1)
subplot(2,4,2);
%f1 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',  OBJ2.vertices,'FaceColor',[95, 84, 199]/255,'LineStyle','none','FaceAlpha',.1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[33, 105, 192]/255,'LineStyle','none','FaceAlpha',.1);
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[191,191,4]/255,'LineStyle','none','FaceAlpha',.1);
f4 = patch('Faces',OBJ2.objects(24).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[117,76,36]/255,'LineStyle','none','FaceAlpha',.1);
%f5 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[215, 68, 4]/255,'LineStyle','none','FaceAlpha',.1);
%f6 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[209,118,0]/255,'LineStyle','none','FaceAlpha',.1);
%f7 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[128, 23, 194]/255,'LineStyle','none','FaceAlpha',.1);
%f8 = patch('Faces',OBJ2.objects(32).data.vertices,'Vertices',OBJ2.vertices,'FaceColor',[196,154,0]/255,'LineStyle','none','FaceAlpha',.1);
axis equal
axis off
view([44.1000,-39.2638])
%view([c1,d1])
camlight('left')
material dull 
hold on
for i = 1:length(ind1)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind1(i),1),Y(ind1(i),1),Z(ind1(i),1),W(ind1(i))*60,'k','filled')
end
for i = 1:length(ind8)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind8(i),1),Y(ind8(i),1),Z(ind8(i),1),W(ind8(i))*60,'r','filled')
end
for i = 1:length(ind7)
%plot3([newX(ind1(i),1);X(ind1(i),1)],[newY(ind1(i),1);X(ind1(i),1)],[newZ(ind1(i),1);Z(ind1(i),1)],'k.-')
scatter3(X(ind7(i),1),Y(ind7(i),1),Z(ind7(i),1),W(ind7(i))*60,'b','filled')
end
colormap(parula);


%Density mapping
allemdisc=[newX1(ind1_1,1),newY1(ind1_1,1),newZ1(ind1_1,1);
newX2(ind1_2,1),newY1(ind1_2,1),newZ1(ind1_2,1);
newX3(ind1_3,1),newY1(ind1_3,1),newZ1(ind1_3,1);
newX4(ind1_4,1),newY1(ind1_4,1),newZ1(ind1_4,1);
newX5(ind1_5,1),newY1(ind1_5,1),newZ1(ind1_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Tb')
ax = gca
ax.CLim = nClim;


allemdisc=[newXs1(ind1_1,1),newYs1(ind1_1,1),newZs1(ind1_1,1);
newXs2(ind1_2,1),newYs1(ind1_2,1),newZs1(ind1_2,1);
newXs3(ind1_3,1),newYs1(ind1_3,1),newZs1(ind1_3,1);
newXs4(ind1_4,1),newYs1(ind1_4,1),newZs1(ind1_4,1);
newXs5(ind1_5,1),newYs1(ind1_5,1),newZs1(ind1_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('EmDisc to Tb')
ax = gca
ax.CLim = nClim;

print(h,['~/Desktop/EmDisc.png'],'-dpng','-r600')
close all


%Density mapping
allemdisc=[newX1(ind2_1,1),newY1(ind2_1,1),newZ1(ind2_1,1);
newX2(ind2_2,1),newY1(ind2_2,1),newZ1(ind2_2,1);
newX3(ind2_3,1),newY1(ind2_3,1),newZ1(ind2_3,1);
newX4(ind2_4,1),newY1(ind2_4,1),newZ1(ind2_4,1);
newX5(ind2_5,1),newY1(ind2_5,1),newZ1(ind2_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind2_1,1),newYs1(ind2_1,1),newZs1(ind2_1,1);
newXs2(ind2_2,1),newYs1(ind2_2,1),newZs1(ind2_2,1);
newXs3(ind2_3,1),newYs1(ind2_3,1),newZs1(ind2_3,1);
newXs4(ind2_4,1),newYs1(ind2_4,1),newZs1(ind2_4,1);
newXs5(ind2_5,1),newYs1(ind2_5,1),newZs1(ind2_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('VE to Tb')
ax = gca
ax.CLim = nClim;


%Now Amnion
%Density mapping
allemdisc=[newX1(ind3_1,1),newY1(ind3_1,1),newZ1(ind3_1,1);
newX2(ind3_2,1),newY1(ind3_2,1),newZ1(ind3_2,1);
newX3(ind3_3,1),newY1(ind3_3,1),newZ1(ind3_3,1);
newX4(ind3_4,1),newY1(ind3_4,1),newZ1(ind3_4,1);
newX5(ind3_5,1),newY1(ind3_5,1),newZ1(ind3_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind3_1,1),newYs1(ind3_1,1),newZs1(ind3_1,1);
newXs2(ind3_2,1),newYs1(ind3_2,1),newZs1(ind3_2,1);
newXs3(ind3_3,1),newYs1(ind3_3,1),newZs1(ind3_3,1);
newXs4(ind3_4,1),newYs1(ind3_4,1),newZs1(ind3_4,1);
newXs5(ind3_5,1),newYs1(ind3_5,1),newZs1(ind3_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Am to Tb')
ax = gca
ax.CLim = nClim;




%Now Tb
%Density mapping
allemdisc=[newX1(ind4_1,1),newY1(ind4_1,1),newZ1(ind4_1,1);
newX2(ind4_2,1),newY1(ind4_2,1),newZ1(ind4_2,1);
newX3(ind4_3,1),newY1(ind4_3,1),newZ1(ind4_3,1);
newX4(ind4_4,1),newY1(ind4_4,1),newZ1(ind4_4,1);
newX5(ind4_5,1),newY1(ind4_5,1),newZ1(ind4_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,1);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to EmDisc')
ax = gca
nClim = ax.CLim;
ax.CLim = nClim;

subplot(2,4,2);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,3);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Am')
ax = gca
ax.CLim = nClim;


subplot(2,4,4);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
%f3 = patch('Faces',OBJ2.objects(12).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Tb')
ax = gca
ax.CLim = nClim;



allemdisc=[newXs1(ind4_1,1),newYs1(ind4_1,1),newZs1(ind4_1,1);
newXs2(ind4_2,1),newYs1(ind4_2,1),newZs1(ind4_2,1);
newXs3(ind4_3,1),newYs1(ind4_3,1),newZs1(ind4_3,1);
newXs4(ind4_4,1),newYs1(ind4_4,1),newZs1(ind4_4,1);
newXs5(ind4_5,1),newYs1(ind4_5,1),newZs1(ind4_5,1)];
f = mvksdensity(allemdisc/400,OBJ2.vertices/400,'Bandwidth',.2,'Kernel','normal','Function','pdf'); % PDF

subplot(2,4,5);
f2 = patch('Faces',OBJ2.objects(4).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
f3 = patch('Faces',OBJ2.objects(8).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to EmDisc')
ax = gca
ax.CLim = nClim;

subplot(2,4,6);
f3 = patch('Faces',OBJ2.objects(16).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to VE')
ax = gca
ax.CLim = nClim;


subplot(2,4,7);
f2 = patch('Faces',OBJ2.objects(20).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Am')
ax = gca
ax.CLim = nClim;

subplot(2,4,8);
f2 = patch('Faces',OBJ2.objects(28).data.vertices,'Vertices',OBJ2.vertices,'FaceVertexCData',f,'FaceColor','interp','LineStyle','none');
axis equal
axis off
view([253.1703,-28.4480])
camlight('left')
material dull 
colormap(parula);
title('Tb to Tb')
ax = gca
ax.CLim = nClim;
