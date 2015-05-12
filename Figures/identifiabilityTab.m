function [m,s, idTab] = identifiabilityTab( X, GroupData)
    %{L, G}
    disp('RawData and GroupData')
    size(X), size(GroupData)
    size(GroupData{1})
    size(GroupData{2})
    L=GroupData{1};
    G=GroupData{2};
    
    %% count different groups
    Groups = unique(G); % different data sets
    Labels = unique(L); % different parameters
    
    for p=1:length(L)
        m(p) = mean(X(:,p));
        s(p) = std( X(:,p));
    end
    
    % min s for each parameter
    for p=1:length(Labels)
        s_min(p) = min( s(L == Labels(p)) );
    end
    s_min
    
    % build identifiability table
    %for i=1:length(L)
    %    idTab()
    %end
    for p=1:length(Labels)
        for e=1:length(Groups)
            L==Labels(p)
            G==Groups(e)
            idTab(p,e) = s_min(p) / s( (L==Labels(p)) & (G==Groups(e)));
        end
    end
    idTab
end