function distancePlot(fileName, bin, cutoff)
    %fileName, file name of raw data
    %bin, interval of angles
    %cutoff
    %For example:
    %       distancePlot('X-RayR_1.5/ARG.csv', 60, 10)
    l = csvread(fileName, 0, 4);
    l = l(:, 1:4); % use Dist_SC, Dist_Block, Phi, Psi
    l = [l(:, 3:4), l(:, 2)]; % use Phi, Psi, Dist_Block
    phiDist = [l(:, 1), l(:, 3)];
    psiDist = [l(:, 2), l(:, 3)];
    %ndhist(phiDist);
    %phiDist = phiDist/max(phiDist(:));
    ndhist(phiDist, 'axis', [0 360 3.55 6.5], '3d', 'columns');
    xticks([0 120 240 360])
%     xlabel('Phi');
%     ylabel('Block Length');
    set(gca,'fontsize',20)
    xlabel('Phi', 'FontWeight','bold')
    ylabel('Block Length', 'FontWeight','bold')
    %return
    figure;
    ndhist(psiDist);
    %ndhist(psiDist, 'axis', [0 360 6.1 6.7]);
    xlabel('Psi');
    ylabel('Block Length');
    binNum = 360/bin+1;
    bins = linspace(0, 360, binNum);
    bins = bins(2:end);
    index = 1;
    figure;
    for i = bins
        for j = bins
            t = l( (l(:, 1) >= i-bin) & (l(:, 1) < i) & (l(:, 2) >= j-bin) & (l(:, 2) < j), :); 
            chi_1 = t(:, 3); % get the column representing chi_1
            ax = subplot(size(bins, 2), size(bins, 2), index);
            index = index + 1;
            %histogram(chi_1);
            %figure
            [bincounts] = histc(chi_1, 0:0.1:10);
%             if max(bincounts) < cutoff
%                 uplimit = cutoff;
%             else
%                 uplimit = max(bincounts);
%             end
%             if uplimit > 2*cutoff
%                 bar(0:10:360, bincounts, 'r');
%             else
%                 bar(0:10:360, bincounts, 'histc');
%             end
            bar(0:0.1:10, bincounts);
            axis(ax, [0 10 0 cutoff]);
            set(ax, 'XTick', [], 'XTickLabel', []);
            xlabel(ax, strcat('(', int2str(i), ',', int2str(j), ')'));
        end
    end
    figure;
    plot3k(l, 'Labels', {'Rotamer Density', 'Phi', 'Psi', 'Block Length'});
end