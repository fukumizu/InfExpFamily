% graphs
%  used for the graphs in the following paper
%  "Density Estimation in Infinite Dimensional Exponential Families" 
%   by Bharath Sriperumbudur, Kenji Fukumizu, Arthur Gretton, Aapo
%   Hyvarinen and Revant Kumar, Journal of Machine Learning Research 2017 
%--------------------------------------------------------------------------

close all;

PRN=1;


dim=[2 4 6 8 10 12 14 16 18 20];
n=[100 200 300 400 500];

%  Gaussian 



kde_gauss_d=[
500	2	0.027008	0.015778	0.991078	0.003904
500	4	0.110518	0.010424	0.980256	0.005765
500	6	0.289641	0.018392	0.963430	0.007449
500	8	0.564937	0.031854	0.940379	0.008713
500	10	0.929083	0.042940	0.911176	0.010306
500	12	1.376057	0.052803	0.875832	0.011363
500	14	1.886278	0.069615	0.838066	0.012321
500	16	2.451274	0.081074	0.796901	0.010858
500	18	3.059019	0.089798	0.750554	0.014849
500	20	3.712199	0.102881	0.701854	0.024851
];

kde_gauss_n=[
100	7	0.430135	0.055803	0.881111	0.024101
200	7	0.416414	0.042877	0.925579	0.017328
300	7	0.419343	0.031322	0.943024	0.009971
400	7	0.416585	0.023343	0.946368	0.008726
500	7	0.413765	0.024421	0.953360	0.007319
];

score_gauss_d = [
500	2	0.029271	0.088894	9.967109e-01	7.058296e-03
500	4	0.025286	0.013200	9.974033e-01	1.956023e-03
500	6	0.047184	0.014078	9.957782e-01	1.486380e-03
500	8	0.080765	0.022615	9.927796e-01	4.308311e-03
500	10	0.116766	0.021746	9.904412e-01	3.600105e-03
500	12	0.162620	0.025723	9.872208e-01	3.266548e-03
500	14	0.220234	0.031539	9.836682e-01	3.460585e-03
500	16	0.281739	0.038204	9.796329e-01	4.685804e-03
500	18	0.357014	0.046828	9.762506e-01	3.751798e-03
500	20	0.447096	0.053868	9.718895e-01	5.225880e-03
];

score_gauss_n=[
100	7	0.297618	0.100766	9.776322e-01	6.002034e-03
200	7	0.146044	0.052354	9.874892e-01	4.742516e-03
300	7	0.096406	0.029817	9.908737e-01	4.860536e-03
400	7	0.078637	0.021061	9.925466e-01	3.483150e-03
500	7	0.063273	0.016809	9.943615e-01	1.978836e-03
];



figure
errorbar(score_gauss_d(:,2),score_gauss_d(:,3),score_gauss_d(:,4),...
    'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([1 21 0 4]);
xlabel('Dimension','FontName','Arial','FontSize',18);
ylabel('Socre objective function','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gauss_d(:,2),kde_gauss_d(:,3),kde_gauss_d(:,4),...
    'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Northwest');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian distribution: n = 500');
hold off;
if PRN
    fgname=sprintf('Gauss_dimension_score_new.eps');
    print('-depsc', fgname);
end

figure
errorbar(score_gauss_n(:,1),score_gauss_n(:,3),score_gauss_n(:,4),...
    'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([50 550 0 0.6]);
xlabel('Sample Size','FontName','Arial','FontSize',18);
ylabel('Socre objective function','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gauss_n(:,1),kde_gauss_n(:,3),kde_gauss_n(:,4),...
    'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian distribution: d = 7');
hold off;
if PRN
    fgname=sprintf('Gauss_datasize_score_new.eps');
    print('-depsc', fgname);
end

%correlation
figure
errorbar(score_gauss_d(:,2),score_gauss_d(:,5),score_gauss_d(:,6),...
    'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([1 21 0.65 1.0]);
xlabel('Dimension','FontName','Arial','FontSize',18);
ylabel('Correlation','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gauss_d(:,2),kde_gauss_d(:,5),kde_gauss_d(:,6),...
    'bs:','LineWidth', 1.5,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Southwest');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian distribution: n = 500');
hold off;
if PRN
    fgname=sprintf('Gauss_dimensoin_corr_new.eps');
    print('-depsc', fgname);
end

figure
errorbar(score_gauss_n(:,1),score_gauss_n(:,5),score_gauss_n(:,6),...
    'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([50 550 0.8 1]);
xlabel('Sample Size','FontName','Arial','FontSize',18);
ylabel('Correlation','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gauss_n(:,1),kde_gauss_n(:,5),kde_gauss_n(:,6),...
    'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Southeast');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian distribution: d = 7');
hold off;
if PRN
    fgname=sprintf('Gauss_datasize_corr_new.eps');
    print('-depsc', fgname);
end



%  Gaussian mixture

kde_gm_d=[
300 2 0.044640 0.016120 0.983560 0.006628
300 4 0.207280 0.055036 0.962107 0.011940
300 6 0.325189 0.097650 0.917177 0.018596
300 8 0.554308 0.043383 0.877444 0.033603
300 10 0.892611 0.094547 0.838866 0.050134
300 12 1.302024 0.177183 0.797976 0.064770
300 14 1.768008 0.279702 0.747725 0.077537
300 16 2.278649 0.410303 0.695050 0.088735
300 18 2.845143 0.560605 0.652686 0.089329
300 20 3.419712 0.726264 0.604779 0.096062
];

kde_gm_n=[
100    7    0.683516    0.331054    0.836487    0.053165
200    7    0.461754    0.146331    0.865775    0.038885
300    7    0.418918    0.031846    0.895939    0.024422
400    7    0.420578    0.032010    0.918667    0.011914
500    7    0.412799    0.021251    0.929500    0.013148
];

score_gm_d = [
300	2	0.149466	0.072206	9.295429e-01	2.318078e-02
300	4	0.197329	0.016823	9.245088e-01	2.112997e-02
300	6	0.302679	0.022612	9.138348e-01	3.150892e-02
300	8	0.410021	0.033810	9.112762e-01	2.409625e-02
300	10	0.502006	0.047139	9.027640e-01	3.273237e-02
300	12	0.603765	0.061444	9.005135e-01	2.665879e-02
300	14	0.721420	0.071264	8.933337e-01	3.097286e-02
300	16	0.846234	0.085084	8.900202e-01	3.598744e-02
300	18	0.989568	0.090535	8.838647e-01	3.272485e-02
300	20	1.165564	0.107757	8.727604e-01	3.923815e-02
];

score_gm_n=[
100	7	0.639280	0.102966	8.947028e-01	3.657254e-02
200	7	0.430522	0.042310	9.019187e-01	4.383203e-02
300	7	0.360312	0.035215	9.112326e-01	3.314864e-02
400	7	0.326811	0.029468	9.038860e-01	4.269999e-02
500	7	0.306099	0.032819	9.075245e-01	3.021784e-02
];

figure
errorbar(score_gm_d(:,2),score_gm_d(:,3),2.*score_gm_d(:,4),'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([0 21 0 5]);
xlabel('Dimension','FontName','Arial','FontSize',18);
ylabel('Score objective function','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gm_d(:,2),kde_gm_d(:,3),2.*kde_gm_d(:,4),'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Northwest');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian mixture: n = 300');
hold off;
if PRN
    fgname=sprintf('GM_dimension_score_new.eps');
    print('-depsc', fgname);
end

figure
errorbar(score_gm_n(:,1),score_gm_n(:,3),2.*score_gm_n(:,4),'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([50 550 0 1]);
xlabel('Sample Size','FontName','Arial','FontSize',18);
ylabel('Score objective function','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gm_n(:,1),kde_gm_n(:,3),2.*kde_gm_n(:,4),'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian mixture: d = 7');
hold off;
if PRN
    fgname=sprintf('GM_datasize_score_new.eps');
    print('-depsc', fgname);
end


%correlation
figure
errorbar(score_gm_d(:,2),score_gm_d(:,5),2.*score_gm_d(:,6),'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([0 21 0.5 1]);
xlabel('Dimension','FontName','Arial','FontSize',18);
ylabel('Correlation','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gm_d(:,2),kde_gm_d(:,5),2.*kde_gm_d(:,6),'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Southwest');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian mixture: n = 300');
hold off;
if PRN
    fgname=sprintf('GM_dimension_corr_new.eps');
    print('-depsc', fgname);
end

figure
errorbar(score_gm_n(:,1),score_gm_n(:,5),2.*score_gm_n(:,6),'ro-','LineWidth', 2,'MarkerFaceColor','r');
axis([50 550 0.7 1]);
xlabel('Sample Size','FontName','Arial','FontSize',18);
ylabel('Correlation','FontName','Arial','FontSize',18);
hold on;
errorbar(kde_gm_n(:,1),kde_gm_n(:,5),2.*kde_gm_n(:,6),'bs:','LineWidth', 2,'MarkerFaceColor','b');
legend('Score match', 'KDE','Location','Southeast');
set(gca,'FontName','Arial','FontSize',18);
title('Gaussian mixture: d = 7');
hold off;
if PRN
    fgname=sprintf('GM_datasize_corr_new.eps');
    print('-depsc', fgname);
end






