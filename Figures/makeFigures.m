function makeFigures()
    close all;
    visualiseABCResultsNEW( '../Data/TumorToyData2D_0.001merr_100pop_GCKI67ECM');
    
     close all;
     visualiseComparison({'../Data/Tumor2dGC', '../Data/Tumor2dGCKI67', '../Data/Tumor2dGCECM', '../Data/Tumor2dGCKI67ECM'});
    
    close all;
    visualiseComparison({'../Data/TumorToyData2D_0.001merr_200pop_GCKI67ECM','../Data/TumorToyData2D_0.001merr_100pop_GCKI67ECM','../Data/TumorToyData2D_0.001merr_10pop_GCKI67ECM'});
    
    close all;
    visualiseABCResultsNEW( '../Data/Tumor3dGCKI67ECM');
    
%     i=1;
%     M = [-1.4 2 1.1 -0.2 -2.3 -3.3 -2.3];
%     for p=1:7
%         m=M(p);
%         
%         for e=1:4
%             s=e;
%             
%             X(:,i) = normrnd(m,s, [100,1]);
%             L(1,i) = p;
%             G(1,i) = e;
%             i=i+1; 
%         end
%     end
%     
%     figure(1)
%     boxplot(X, {L,G})
%     
%     figure(2)
%     [m,s, idTab] = identifiabilityTab(X, {L,G})
%     errorbar(m,s)
%     
%     matrix2latex(idTab, unique(G), unique(L))
end
