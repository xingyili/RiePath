function res = RiePath(GeneExprsseionForAllNetwork)

% RiePath for training cohort

% INPUT
% GeneExprsseionForAllNetwork: a cell array with 1*294 dimension, corresponding to the number of pathways.
% Each element in the array is (1+m)*n, the first row is the label of samples (0:Normal, -1: DiseaseA, 1:DiseaseB)
% m is the differetial expression features belonging to each pathway, n is the samples.

for i = 1:length(GeneExprsseionForAllNetwork)
    EachNetwork = GeneExprsseionForAllNetwork{1,i};
    EachNetworkData = EachNetwork(2:end,2:end);
    EachNetworkLabel  = EachNetwork(1,2:end);
    NormalData = EachNetworkData(:,EachNetworkLabel==0);
    DiseaseA = EachNetworkData(:,EachNetworkLabel==-1);
    DiseaseB = EachNetworkData(:,EachNetworkLabel==1);
    NormalCov = covariances(NormalData);
    DiseaseA_Cov = zeros(size(DiseaseA,1),size(DiseaseA,1),size(DiseaseA,2));
    DiseaseB_Cov = zeros(size(DiseaseB,1),size(DiseaseB,1),size(DiseaseB,2));
    DiseaseA_Distance_Single = [];
    DiseaseB_Distance_Single = [];
    for j = 1:size(DiseaseA,2)
        DiseaseA_Cov(:,:,j) = covariances([NormalData,DiseaseA(:,j)]);
    end
    [Feat,~] = Tangent_space(DiseaseA_Cov,NormalCov);  
    DiseaseA_Tangent = real(Feat);
    DiseaseA_Distance_Single = sqrt(sum(DiseaseA_Tangent.^2,1));
    for j = 1:size(DiseaseB,2)
        DiseaseB_Cov(:,:,j) = covariances([NormalData,DiseaseB(:,j)]);
    end
    [Feat,~] = Tangent_space(DiseaseB_Cov,NormalCov);  
    DiseaseB_Tangent = real(Feat);
    DiseaseB_Distance_Single = sqrt(sum(DiseaseB_Tangent.^2,1));
    DiseaseA_Distance = [DiseaseA_Distance;DiseaseA_Distance_Single];
    DiseaseB_Distance = [DiseaseB_Distance;DiseaseB_Distance_Single];
end

res = [DiseaseA_Distance,DiseaseB_Distance];

end