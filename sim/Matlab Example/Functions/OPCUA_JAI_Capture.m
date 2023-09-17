function [img] = OPCUA_JAI_Capture(save,name)

%% Connect to OPCUA
uaClient = opcua('172.31.1.236', 4840);
connect(uaClient)

%% Define All Relevant Nodes
NODE_Printer1 = findNodeByName(uaClient.Namespace,'Camera1','-once');
NODE_C1f_Ready = findNodeByName(NODE_Printer1,'C1f_Ready','-once');
NODE_C1f_End = findNodeByName(NODE_Printer1,'C1f_End','-once');
NODE_C1c_Snap = findNodeByName(NODE_Printer1,'C1c_Snap','-once');
NODE_C1c_Name = findNodeByName(NODE_Printer1,'C1c_Name','-once');
NODE_C1d_Img = findNodeByName(NODE_Printer1,'C1d_Img','-once');
NODE_C1d_Img_h = findNodeByName(NODE_Printer1,'C1d_Img_h','-once');
NODE_C1d_Img_w = findNodeByName(NODE_Printer1,'C1d_Img_w','-once');

%% Communication

[C1f_Ready,~,~] = readValue(uaClient,NODE_C1f_Ready);

if C1f_Ready==true
    writeValue(uaClient, NODE_C1c_Snap, true);
    [C1f_End,~,~] = readValue(uaClient,NODE_C1f_End);
    while C1f_End==false
        [C1f_End,~,~] = readValue(uaClient,NODE_C1f_End);
    end
    %pause(2)
    [C1d_Img,~,~] = readValue(uaClient,NODE_C1d_Img);
    [C1d_Img_h,~,~] = readValue(uaClient,NODE_C1d_Img_h);
    [C1d_Img_w,~,~] = readValue(uaClient,NODE_C1d_Img_w);
    img = reshape(C1d_Img,C1d_Img_h,C1d_Img_w);
    if save==1
        writeValue(uaClient,NODE_C1c_Name,name);
        imwrite(img,append(name,'.png'));
    end
    writeValue(uaClient, NODE_C1c_Snap, false);
else
    img = 0;
    disp('CamSys not ready')
end

end

