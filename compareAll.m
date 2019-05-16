function rotamerPlot(res_1, res_2, xRange, resName, binSize)
    % compare the side-chain distance distribution of proteins from X-Ray
    % and EM
    % For example:
    %       compareAll('../X-RayR_1.5/ASN.csv', '../EM_3.5_4.0/ASN.csv', 0:0.05:10, 'ASN', 0.05)
    % read matrix
    l = getMatrix(res_1); %get X-Ray data
    l_2 = getMatrix(res_2); %get EM data
    
    % phi
    dist_1 = l(:, 1);
    dist_2 = l_2(:, 1);
    process(dist_1, dist_2, 0:5:360, resName, 5, 'phi');
    
    % psi
    figure;
    dist_1 = l(:, 2);
    dist_2 = l_2(:, 2);
    process(dist_1, dist_2, 0:5:360, resName, 5, 'psi');
    
    % get the column representing dist_1, which is the distance between CA and the side-chain mass center
    figure;
    dist_1 = l(:, 3); %d_sc of X-Ray
    dist_2 = l_2(:, 3); %d_sc of EM    
    process(dist_1, dist_2, xRange, resName, binSize, 'SC');
    
    % plot histogram for block
    figure;
    dist_1 = l(:, 4); % get the column representing dist_1, which is the distance between CA and the side-chain mass center
    dist_2 = l_2(:, 4);
    process(dist_1, dist_2, xRange, resName, binSize, 'Block');
    
    % chi_1
    figure;
    dist_1 = l(:, 5);
    dist_2 = l_2(:, 5);
    process(dist_1, dist_2, 0:5:360, resName, 5, 'chi1');
end

function process(dist_1, dist_2, xRange, name, binSize, label)
    % dist_1, distance of X-Ray data
    % dist_2, distance of EM data
    % xRange, resolution in x-axis, for instance, 2:0.1:7
    % name, residue name
    % binSize, size of bin
    % label, SC or Block
    
    % process X-Ray data
    %histogram(dist_1, xRange, 'Normalization', 'pdf', 'FaceColor', 'r'); % red color for the first residue
    [N_1_pdf, edges_1] = histcounts(dist_1, xRange, 'Normalization', 'pdf');
    peak = max(N_1_pdf);
    N_1_npdf = N_1_pdf/peak;
    plotHistcounts(N_1_npdf, edges_1, 'r', 2); % red for X-Ray
    [N_1_prob, edges_1] = histcounts(dist_1, xRange, 'Normalization', 'probability');
    %outputFile(strcat(name, '_X-Ray_', label, '_', num2str(binSize), '.csv'), N_1_npdf, N_1_pdf, N_1_prob, edges_1);
    %histogram(dist_1, xRange, 'Normalization', 'probability', 'FaceColor', 'r'); % red color for the first residue
    
    % process EM data
    hold on;
    %histogram(dist_2, xRange, 'Normalization', 'pdf', 'FaceColor', 'b'); % blue color for the second residue
    [N_2_pdf, edges_2] = histcounts(dist_2, xRange, 'Normalization', 'pdf');
    N_2_npdf = N_2_pdf/peak;
    plotHistcounts(N_2_npdf, edges_2, 'b', 1); % blue for EM
    [N_2_prob, edges_2] = histcounts(dist_2, xRange, 'Normalization', 'probability');
    %outputFile(strcat(name, '_EM_', label, '_', num2str(binSize), '.csv'), N_2_npdf, N_2_pdf, N_2_prob, edges_2);
    %histogram(dist_2, xRange, 'Normalization', 'probability', 'FaceColor', 'b'); % blue color for the second residue
    title(strcat(label, ' Histogram of name'));
end

function outputFile(fileName, npdf, pdf, prob, edges)
    low = edges(1:end-1);
    high = edges(2:end);
    bins = (low+high)/2;
    csvwrite(fileName, [bins' npdf' pdf' prob']);
end

function plotHistcounts(N, edges, color, width)
    low = edges(1:end-1);
    high = edges(2:end);
    bins = (low+high)/2;
    %bar(bins, N);
    plot(bins, N, color, 'LineWidth',width);
end
    
function l = getMatrix(fileName)
    l = csvread(fileName, 0, 4);
    a = l(:, 1:2);%d_sc, d_block
    b = l(:, 3:4);%phi, psi
    c = l(:, 5); % chi_1
    l = [b a c]; % use phi, psi, d_sc, d_block, chi_1
    %l = l(:, 1:4); % use phi, psi, d_sc, d_block
end