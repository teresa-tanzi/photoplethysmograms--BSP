files = dir('TBME2013-PPGRR-Benchmark_R3/data/*');

all_y = [];
all_y_pred = [];

size = 144002;

for f = 1 : length(files)
	if files(f).isdir == 0
		data = load(strcat('TBME2013-PPGRR-Benchmark_R3/data/', files(f).name));
		labels = data.labels.pleth.artif.x;
				
		%{
		if length(labels) > 0
			fprintf("%s\n", files(f).name);
			display(labels);
		end
		%}
		
		% salvo in y tutte le coordinate di punti considerati 
		y = zeros(size,1);
		
		for i = 1 : length(labels)
			% considero solo gli indici pari, perché indicano la fine
			% dell'intervallo
			if mod(i, 2) == 0
				%y = [y, labels(i-1) : labels(i)];
				y(labels(i-1) + 1 : labels(i) + 1) = ones(length(labels(i-1) + 1 : labels(i) + 1), 1);
			end 
		end
		
		detected = load(strcat('data/', files(f).name(1:end-4), '_art.mat'));
		artifacts = detected.artifacts;
		
		%all_y = [all_y, y];
		all_y = cat(1, all_y, y);
		
		% salvo tutti i punti considerati artefatto dalla mia analisi
		y_pred = zeros(size,1);
		
		for i = 1 : length(artifacts)
			%y_pred = [y_pred, (artifacts(i).start : artifacts(i).end)];
			y_pred(artifacts(i).start : artifacts(i).end) = ones(length(artifacts(i).start : artifacts(i).end), 1);
		end
		
		%all_y_pred = [all_y_pred, y_pred];
		all_y_pred = cat(1, all_y_pred, y_pred);
		
		%figure(1)
		%cm = confusionchart(y, y_pred);
		%cm.RowSummary = 'row-normalized';
		%cm.ColumnSummary = 'column-normalized';
		%savefig(strcat('conf/', files(f).name(1:end-4),'.fig'));
	end
end

classLabels = ["Correct", "Artifact"];

figure(2)
cm_all = confusionchart(all_y, all_y_pred);
%cm_all.ColumnSummary = 'column-normalized';
%cm_all.RowSummary = 'row-normalized';
%plotconfusion(all_y,all_y_pred)
%savefig(strcat('conf/all.fig'));