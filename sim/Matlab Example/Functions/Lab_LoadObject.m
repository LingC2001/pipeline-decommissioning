function [] = Lab_LoadObject(iO,type,pos,properties,ColourF,ColourE,AlphaF,AlphaE,WidthE,ax1)
global DTL


if type == 0
    T = pos;
    name = properties;

    DTL.Object{iO}.stl = stlread(append("STL/Objects/",name,".stl"));
    DTL.Object{iO}.Transform = hgtransform(ax1); %gg
    DTL.Object{iO}.fig=trimesh(DTL.Object{iO}.stl,'Parent',ax1,'Parent',DTL.Object{iO}.Transform); %%hh
    DTL.Object{iO}.fig.FaceColor = ColourF;
    DTL.Object{iO}.fig.EdgeColor = ColourE;
    DTL.Object{iO}.fig.FaceAlpha = AlphaF;
    DTL.Object{iO}.fig.EdgeAlpha = AlphaE;
    DTL.Object{iO}.fig.LineWidth = WidthE;
    DTL.Object{iO}.fig.FaceLighting = "flat";
    
    DTL.Object{iO}.Transform.Matrix = T;
end


if type == 1
    dx = properties(1);
    dy = properties(2);
    dz = properties(3);

    x = pos(1);
    y = pos(2);
    z = pos(3);

    pdx = x+dx;
    pdy = y+dy;
    pdz = z+dz;

    Xc = [x pdx pdx x x x; pdx pdx x x pdx pdx; pdx pdx x x pdx pdx; x pdx pdx x x x];
    Yc = [y y pdy pdy y y; y pdy pdy y y y; y pdy pdy y pdy pdy; y y pdy pdy pdy pdy];
    Zc = [z z z z z pdz; z z z z z pdz; pdz pdz pdz pdz z pdz; pdz pdz pdz pdz z pdz];

    DTL.Object{iO}.Transform = hgtransform(ax1); %gg
    DTL.Object{iO}.fig=patch('Xdata',Xc,'YData',Yc,'ZData',Zc,'Parent',ax1,'Parent',DTL.Object{iO}.Transform); %%hh
    DTL.Object{iO}.fig.FaceColor = ColourF;
    DTL.Object{iO}.fig.EdgeColor = ColourE;
    DTL.Object{iO}.fig.FaceAlpha = AlphaF;
    DTL.Object{iO}.fig.EdgeAlpha = AlphaE;
    DTL.Object{iO}.fig.LineWidth = WidthE;
    DTL.Object{iO}.fig.FaceLighting = "flat";
end
end