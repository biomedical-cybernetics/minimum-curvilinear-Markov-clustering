function [comm,nc] = mc_mcl(x,C,cygwin_path,dist,factor,max_nc,mod)


%mcl function that works with the MCL algorithm tool installed with cygwin
%on a windows environment.
% Inputs:
%
% x:  matrix data where in the rows are the samples and in the
% columns the features.
%
% C: number of clusters that the user is trying to find.
%
% dist:  number that represent the distance applied on "x" to construct the graph
%       [1] apply the MC matrix created with euclidean distance
%       [2] apply the MC matrix created with correlation distance
%       [3] apply the MC matrix created with correlation distance
% default: 1
%
% cygwin_path: cygwin path to run MCL
%
% factor: the option for the MC distance:
%       [1] nothing
%       [2] sqrt
%       [3] log(1+x)
% default: 1
%
%
% cutoff:  threshold in the similarity of the samples used to compute
% the graph. The higher the cutoff it is, the more number of components the graph may have.
% default: the program calculate for itselfs the higher cutoff where the
% graph contains 1 component.
%
% max_nc: maximum number of components of the graph for choosing the cutoff
% default : 1
%
% mod: perform modality:
%          [1] choosing cutoff automatically, matrix-dataset goes from 0
%          until 1-cutoff.
%          [2] sparcifying network using cutoff selected by the user
% default: 1


%% Initialization of data

if nargin<3, error('The number of Inputs must be at least 3');  end
if nargin<4, dist=1; factor=1; max_nc = 1; mod = 1; end
if nargin<5, factor=1; max_nc = 1; mod = 1; end
if nargin<6, max_nc = 1; mod = 1; end
if nargin<7, mod = 1; end

if max_nc > C-1
    error('Max number of components can not be greater than the number of cluster are beeing looking for');
end

switch dist
    case 1
        %%%%%%%%%%%%%%%For MC-matrix euclidean norm%%%%%%%%%%%%%%%
        if factor==1
            matr=squareform(pdist(x,'euclidean')); matr=tril(matr);
        elseif factor==2
            matr=sqrt(squareform(pdist(x,'euclidean'))); matr=tril(matr);
        elseif factor==3
            matr=log(1+squareform(pdist(x,'euclidean'))); matr=tril(matr);
        end
        
        xx= graphallshortestpaths(graphminspantree(sparse(matr),'method','kruskal'),'directed','false');
        
        if ~issymmetric(xx) % xx return as a nonsymetric matrix (not detected by eye)
            xx = triu(xx,0) + triu(xx,1)';
        end
        
        x1 = xx./max(max(xx));
        x1= 1-x1;
    case 2
        %%%%%%%%%%%%%%%For MC-matrix correlation norm%%%%%%%%%%%%%%%
        if factor==1
            matr=squareform(pdist(x,'correlation')); matr=tril(matr);
        elseif factor==2
            matr=sqrt(squareform(pdist(x,'correlation'))); matr=tril(matr);
        elseif factor==3
            matr=log(1+squareform(pdist(x,'correlation'))); matr=tril(matr);
        end
        
        xx= graphallshortestpaths(graphminspantree(sparse(matr),'method','kruskal'),'directed','false');
        x1= 1-xx;
    case 3
        %%%%%%%%%%%%%%%For MC-matrix spearman norm%%%%%%%%%%%%%%%
        if factor==1
            matr=squareform(pdist(x,'spearman')); matr=tril(matr);
        elseif factor==2
            matr=sqrt(squareform(pdist(x,'spearman'))); matr=tril(matr);
        elseif factor==3
            matr=log(1+squareform(pdist(x,'spearman'))); matr=tril(matr);
        end
        
        xx= graphallshortestpaths(graphminspantree(sparse(matr),'method','kruskal'),'directed','false');
        x1= 1-xx;
end


nf = size(x,2);

if mod == 1
    
    [x1,nc]=choose_cut(x1,max_nc);
    
    [comm]=compute_mcl(x1,C,cygwin_path,nf);
    
elseif mod == 2
    %user put the cutoff
    cutoff = input('Insert the cutoff -> ');
    
    %%%% counting the number of components inside the graph %%%%
    x1(x1<cutoff) = 0;
    x1(x1>=cutoff) = x1(x1>=cutoff)-cutoff;
    
    S=sparse(x1);
    [nc,~]=graphconncomp(S,'Directed', false);
    
    [comm]=compute_mcl(x1,C,cygwin_path,nf);
end

end

function [x1,nc]=choose_cut(x1,max_nc)

uniq_weigths = unique(round(x1,2));
idxs = find(uniq_weigths>0);

for i=1:length(idxs)
    
    cutoff = uniq_weigths(idxs(i));
    
    tmp_x1 = x1;
    tmp_x1(tmp_x1<cutoff) = 0;
    tmp_x1(tmp_x1>=cutoff) = tmp_x1(tmp_x1>=cutoff)-cutoff;
    
    S=sparse(tmp_x1);
    [nc,~]=graphconncomp(S,'Directed', false);
    if nc > max_nc
        
        if i == 1
            warning('For the first cutoff the number of components %d are already greater than %d',nc,max_nc);
            break
        end
        cutoff = uniq_weigths(idxs(i-1));
        tmp_x1 = x1;
        tmp_x1(tmp_x1<cutoff) = 0;
        tmp_x1(tmp_x1>=cutoff) = tmp_x1(tmp_x1>=cutoff)-cutoff;
        
        S=sparse(tmp_x1);
        [nc,~]=graphconncomp(S,'Directed', false);
        break
    end
end

x1(x1<cutoff) = 0;
x1(x1>=cutoff) = x1(x1>=cutoff)-cutoff;

end

function [comm]=compute_mcl(x1,C,cygwin_path,nf)
%% Writing .mci file
currdir = strrep([pwd '/'],'\','/');
fileprefix = 'MCL';
fileinput = [currdir fileprefix '.in'];
fileoutput = [currdir fileprefix '.out'];

fileID = fopen(fileinput,'w');
fprintf(fileID,'(mclheader\n');
fprintf(fileID,'mcltype matrix\n');
fprintf(fileID,'dimensions %dx%d\n',size(x1));
fprintf(fileID,')\n');
fprintf(fileID,'(mclmatrix\n');
fprintf(fileID,'begin\n');

for i = 1:size(x1,1)
    fprintf(fileID,'%d:%d ',i-1,nf);
    for j = 1:size(x1,2)
        if(x1(j,i)>0)
            fprintf(fileID,'%d:%f ',j-1,x1(j,i));
        end
    end
    fprintf(fileID,'$\n');
end

fprintf(fileID,')\n');
fclose(fileID);

%% Computing MCL

minInflation = 1.1;
maxInflation = 20;
resolution = [0.1,0.01,0.001];

left = 1;
I = minInflation:resolution(1):maxInflation;
right = length(I);
i = left + round((right-left)/2);
c = 1;

for res = resolution
    if left > right
        if lastL == lastR
            if c-1 > C
                lastL = lastL - 1;
            else
                lastR = lastR + 1;
            end
            if lastR > length(I) || lastL < 1
                break
            end
        end
        I = I(lastL):res:I(lastR);
        left = 1;
        right = length(I);
        i = left + round((right-left)/2);
    end
    
    while c-1 ~= C
        
        % Run MCL
        command = [cygwin_path '\bash --login -c "mcl ' fileinput ' -I ' num2str(I(i)) ' -V all -o ' fileoutput '"'];
        system(command);
        
        % Reading output file
        fileID = fopen(fileoutput,'r');
        content = textscan(fileID,'%s','delimiter','\n','collectoutput',true);
        comm = NaN(size(x1,1),1);
        c = 1;
        last = '$';
        for j=8:length(content{1})-1
            ln = content{1}(j);
            sp = strsplit(ln{1});
            if last == '$'
                sp = sp(2:end);
            end
            last = sp{end};
            if last == '$'
                sp = sp(1:end-1);
            end
            for k=1:length(sp)
                idx = str2double(sp{k})+1;
                comm(idx) = c;
            end
            if last == '$'
                c = c + 1;
            end
        end
        
        % Check number of clusters
        lastL = left;
        lastR = right;
        if c-1 > C
            right = i-1;
        elseif c-1 < C
            left = i+1;
        end
        if left > right
            fclose(fileID);
            delete(fileoutput);
            break
        end
        fclose(fileID);
        delete(fileoutput);
        i = left + round((right-left)/2);
    end
    
    if c-1 == C
        break
    end
end
delete(fileinput);
end

