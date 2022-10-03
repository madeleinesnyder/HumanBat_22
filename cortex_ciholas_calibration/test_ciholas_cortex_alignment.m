% This test uses test data from 220103

ciholas = load('/home/batlab/Desktop/HumanBat/data/tests/processed/220103/b149f/ciholas/extracted_220103_cdp_1.mat')

cortex = load('/home/batlab/Desktop/HumanBat/data/tests/processed/220103/b149f/cortex/220103_14592_tracking_1_track.mat')

ciholas2cortex = load('ciholas2cortex_scaling_factors.mat').ciholas2cortex;

alignedCiholasCortex = alignCiholasCortex(cortex,ciholas,ciholas2cortex); % align ciholas with cortex

posK = alignedCiholasCortex.ciholas.kq.pos;
posM = alignedCiholasCortex.ciholas.ms.pos;
posB = alignedCiholasCortex.cortex.avgMarkerPos;

figure;
scatter(posK(:,1), posK(:,2), 3,'red','filled');
hold on
scatter(posM(:,1), posM(:,2), 3,'green','filled');
hold on
scatter(posB(:,1), posB(:,2), 3,'blue','filled');
legend('KQ', 'MS', 'Bat');