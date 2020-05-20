% an algorithm for real-time pulse waveform segmentation and artifact
% detection in photoplethysmograms

%% dataset: http://www.capnobase.org/database/pulse-oximeter-ieee-tbme-benchmark/

% The recordings obtained included PPG (100 Hz) signals
% raw PPG signal for 42 cases of 8 min duration 
sampling_rate = 300;    % 144000 / 480 (8 * 60) (c'è scritto nei file)

files = dir('TBME2013-PPGRR-Benchmark_R3/data/*');

% memorizzo in questo array il numero di volte in cui triggero un
% check
artifact_triggered = zeros(1,13);

% array per analizzare la distribuzione degli errori
rise_time = zeros(1,0);
s_d_ratio = zeros(1,0);
n_peaks = zeros(1,0);
systolic_peak = zeros(1,0);

for f = 1 : length(files)
	if files(f).isdir == 0
 		%if files(f).name == '0147_8min.mat' 
		%if files(f).name == '0150_8min.mat' 
		%if files(f).name == '0122_8min.mat'
		%if files(f).name == '0029_8min.mat'
		
		data = load(strcat('TBME2013-PPGRR-Benchmark_R3/data/', files(f).name));
			
		raw_signal = data.signal.pleth.y;

		N = length(raw_signal); % 144001
		time_axis = (0:1/sampling_rate:(N-1)*1/sampling_rate);

		figure(1)
		plot(time_axis, raw_signal)
		xlabel('Time [s]') 
		ylabel('Ampitude') 

		%% 1. Raw Signal Artifact Detection

		% detects clipping values in the raw signal
		% the thresholds for top and bottom clipping must be determined only once
		% for each type of recording system: these values where derived fro the
		% signal channel minimum and maximum

		% check 1: clipping bottom detected
		% RawSignalValue <= 0
		% il mio segnale però è centrato in 0, quindi il minimo (dove fa clipping)
		% è -10.1600
		%display(min(raw_signal));

		% check 2: clipping top detection
		% RawSignalValue >= Sleep/Ergometry Lab Uml: 16777215
		%                   ICU PhysioBank: 1.0
		%                   PICU CSL: 1023
		%display(max(raw_signal));
		%marks = strings(1,size(raw_signal,1));
		disturbed = zeros(1,size(raw_signal,1));

		for i = 1:size(raw_signal)
			if raw_signal(i) <= -10.1600
				%fprintf("%u\n",i);          % %u: unsigned int base 10
				%marks(i) = "disturbed";
				disturbed(i) = -1;
			end
			if raw_signal(i) >= 10.2400
				disturbed(i) = 1;
			end
		end

		%display(disturbed);

		x = find(disturbed == -1);
		y = -10.1600 * ones(length(x),1);

		hold all;
		scatter((x-1)/sampling_rate, y, 'rx');

		x = find(disturbed == 1);
		y = 10.2400 * ones(length(x),1);

		scatter((x-1)/sampling_rate, y, 'm+');
		title('Raw Signal Artifact Detection');
		savefig(strcat('img/', files(f).name(1:end-4),'_raw.fig'));
		hold off;

		%% 2. Low-pass Filter: 4th order Butterworth

		% the power spectrum of the PPG consistes of frequency components up to 15
		% Hz with a dominant power range between 6 and 8 Hz
		% The cutoff frequency of the low-pass filter was set to 15 Hz in order to
		% remove higher artifact frequencies (e.g.: power line interference)

		filter_order = 4;
		cutoff_frequency = 15/(sampling_rate/2);
		[b,a] = butter(filter_order, cutoff_frequency, 'low');
		%freqz(b,a);

		low_pass_signal = filtfilt(b,a,raw_signal);
		
		figure(2)
		plot(time_axis, low_pass_signal);
		xlabel('Time [s]'); 
		ylabel('Ampitude');
		title('Low-pass Filter');
		savefig(strcat('img/', files(f).name(1:end-4),'_lp.fig'));

		%% 3. High-pass Filter: 4th order Butterworth

		% to suppress the dc part without removing any information about the
		% autonomic nervous system control of the cardiovascular system, the cutoff
		% frequency of the high-pass filter in the third stage was set to 0.01 Hz

		filter_order = 4;
		%cutoff_frequency = .01/(sampling_rate/2);
		cutoff_frequency = .4/(sampling_rate/2);
		[b,a] = butter(filter_order, cutoff_frequency, 'high');
		%freqz(b,a);

		high_pass_signal = filtfilt(b,a,low_pass_signal);
		
		%{
		figure(4)
		plot(time_axis, raw_signal);
		hold on;
		%plot(time_axis, low_pass_signal);
		plot(time_axis, high_pass_signal);
		xlabel('Time [s]') 
		ylabel('Ampitude') 
		title('Low pass and high pass filter');
		hold off;
		%}
		
		figure(3)
		plot(time_axis, high_pass_signal);
		xlabel('Time [s]'); 
		ylabel('Ampitude');
		title('High-pass Filter');
		hold on
		
		%% 4. Pulse Waves Valleys and Peak Detection

		% the filtered PPG and the raw signal annotation of the PPG are stored in
		% two ring buffers the size of two times the longest permitted duration of
		% a pulse wave (4.8 s each)

		max_pwd = 2.4 * sampling_rate;  % 2.4 è il limite teorico in secondi
		window_size = 2 * max_pwd;

		% array vuoto che conterrà i valori della moving average man mano che la
		% calcolo
		moving_average = zeros(length(high_pass_signal), 1);

		% array di zeri che avrà 1 in corrispondenza dei potenziali picchi
		%peaks = zeros(length(high_pass_signal), 1);
		peaks = zeros(0,1);

		% array di zeri che avrà 1 in corrispondenza delle potenziali valli
		%valleys = zeros(length(high_pass_signal), 1);
		valleys = zeros(0,1);

		% array di zeri che conterr 1 in corrispondenza degli artefatti
		artifact = zeros(length(high_pass_signal), 1);

		% informazioni sugli artefatti trovati
		artifacts = zeros(0,1);

		% informazioni sulle pulse wave passate
		waves = zeros(0,1);

		% booleano che determina quando uscire dal ciclo
		cont = true;
		% valore che mi dice da che sample cominciare la finestra
		start = 1;
		%start = 78600;

		% 1) scegliere una PWD teorica iniziale
		pwd = max_pwd;
		% inizializzo pwa a 0, in questo modo la prima pwa verrà sempre considerata
		% corretta
		pwa = 0;

		% 2) cominciare a filtrare il segnale dall'inizio con una filtro a media 
		% mobile lungo il 75% della PWD teorica. Filtra 4.8s di segnale e quindi 
		% cerca picchi e valli e localizza la nuova PWD.
		while cont
			% filtro 4.8 secondi di segnale
			window = high_pass_signal(start : min(start + window_size -1, length(high_pass_signal)));
			%fprintf("dimensione finestra: %u\n", length(window));

			% an adaptive threshold is calculated by applying a moving average filter
			% with a span size of 75% of the last valid PWD (pulse wave duration)
			b = (1/(pwd * .75)) * ones(1, fix(pwd * .75));
			a = 1;
			ma = filter(b, a, window);
			%fprintf("ma: %u\n", length(ma));

			%{
			hold on
			plot(time_axis(1 : start + window_size - 1), ma);
			hold off
			%}

			% the absolute maximum is identified as a potential pulse wave peak in 
			% every part of the signal above this threshold the absolute minimum is 
			% selectes as a potential pulse wave valley in every signal part below
			% this allow the detection of potential valleys and peaks even in signals
			% with dc drift

			% ciclo sui sample della finestra per trovare massimi e minimi assoluti
			% sopra e sotto la soglia 

			% memorizzo in un array lungo quanto la finestra le distanze tra il
			% segnale e la moving average

			distance = zeros(length(window), 1);

			for i = 1 : length(window)
				distance(i) = (window(i) - ma(i));
			end

			% scorro le distanze e cerco massimi e minimi

			x_max = zeros(0,1);
			x_min = zeros(0,1);

			% salvo se mi trovo sopra o sotto alla moving average per capire se sto
			% cercando un massimo o un minimo
			above = 1;

			maximum = 1;
			minimum = 1;

			if distance(1) < 0
				above = -1;
			end

			for i = 1 : length(distance)
				if above * distance(i) < 0	% sto cambiando segno della distanza
					%fprintf("intersezione: %u\n", i);
					above = -1 * above;
					%fprintf("above: %u\n", above);

					% salvo il massimo o il minimo che ho trovato al giro
					% precedente
					if above < 0	% vuol dire che prima cercavo il massimo
						x_max = [x_max, maximum];
					else
						%fprintf("min: %u\n", minimum);
						x_min = [x_min, minimum];
					end

					% azzero nuovamente massimo e minimo perché devo cercarne di
					% nuovi
					maximum = i;
					minimum = i;
				end

				% se mi trovo sopra alla moving average, cerco un massimo
				if above > 0
					if window(i) > window(maximum)
					%if distance(i) > distance(maximum)
						maximum = i;
					end
				else
					%if distance(i) < distance(minimum)
					if window(i) < window(minimum)
						minimum = i;
					end
				end
			end

			%{
			hold on
			scatter((x_min + start)/sampling_rate, high_pass_signal(x_min+start), 'g*');
			hold off
			%}

			%{
			display(x_max);
			display(x_min);
			%}
			
			% wrongly detected potential peaks and corresponding valleys belonging to
			% diastolic peaks are discarded by demanding that the vertical distance
			% between them multiplied by a factor (3) must be bigger than the previous
			% valid PWA (pulse wave ampitude)

			% cerco valli e picchi immediatamente successivi e calcolo la loro
			% distanza

			remove_max = zeros(0,1);
			remove_min = zeros(0,1);

			%fprintf("len min: %u, len max: %u\n", length(x_max), length(x_min));

			% TODO: faccio partire da 2, ma non so se è ok
			for i = 1 : min(length(x_max), length(x_min))
			%for i = 1 : length(x_max)
				if x_min(1) < x_max(1)		% nel mio intervallo c'è prima una valle
					%fprintf("pwa: %u, dist: %u\n", pwa, abs(x_max(i) - x_min(i))*3);
					if abs(window(x_max(i)) - window(x_min(i))) * 3 < pwa
						%fprintf("x min: %u\n", x_min(i));
						remove_max = [remove_max, x_max(i)];
						remove_min = [remove_min, x_min(i)];
					end
				else						% c'è prima un picco
					if i < min(length(x_max), length(x_min))
						if abs(window(x_max(i+1)) - window(x_min(i))) * 3 < pwa
							%fprintf("x min: %u\n", x_min(i));
							remove_max = [remove_max, x_max(i+1)];
							remove_min = [remove_min, x_min(i)];
						end
					end
				end
			end	
			
			%{
			if start == 78459
				fprintf("start: %u\n", start);
				fprintf("pwa: %u", pwa);

				display(x_max);
				display(x_min);
				
				fprintf("x_max: %u, x_min: %u\n", (start+337)/300, (start+132)/300);
				fprintf("y_max: %u, y_min: %u, y_diff: %u\n", window(337), window(132), window(337)-window(132));
				
				%{
				for j = 1 : length(x_max)
					y_max = window(x_max);
					y_min = window(x_min);
					fprintf("pwa: %u\n", (y_max(j) - y_min(j)));
					fprintf("diff > pwa: %d\n", (((y_max(j) - y_min(j))*3)>pwa));
				end
				%}

				display(remove_max);
				display(remove_min);
			end
			%}
			
			% se scarto la prima valle, allora la accorpo all'artefatto o
			% all'onda trovata al giro precedente
			%{
			if length(remove_min) > 0
				if remove_min(1) == x_min(1)
					if waves(length(waves)).end == start - 1
						%fprintf("start: %u\n", start);
						x_min_prov = setdiff(x_min, remove_min);
						waves(length(waves)).end = x_min_prov(1) - 1;
					end
					
					if artifacts(length(artifacts)).end == start - 1
						%fprintf("start: %u\n", start);
						x_min_prov = setdiff(x_min, remove_min);
						artifacts(length(artifacts)).end = x_min_prov(1) - 1;
					end
				end
			end
			%}
			
			% rimuovo picchi e valli scartati
			x_max = setdiff(x_max, remove_max);
			x_min = setdiff(x_min, remove_min);
			
			%{
			display(x_max);
			display(x_min);
			%}
			
			% 3) a partire da punto in cui ha localizzato la nuova PWD, ricomincia a
			% filtrare 4.8 s ma facendo la media mobile su un numero di campioni pari
			% al 75% della nuova PWD (il filtro è adattativo) e cerca un nuovo
			% picco/valle validi per aggiornare la PWD e quindi la lunghezza del filtro.
			% È probabile che 4.8s sia anche troppo lunga come finestra, ma gli autori 
			% l'hanno scelta il modo che corrisponda ad una lunghezza fisiologica.

			% the results of these detections are stored synchronously to the filtered
			% PPG ring buffer in the PPG annotation buffer

			if length(x_min) > 2	% ho individuato una nuova pwd
				% aggiorno la pwd
				if (x_min(2) - x_min(1)) <= max_pwd
					pwd = abs(x_min(2) - x_min(1));
					pwa = abs(window(x_max(1)) - window(x_min(1)));
				end

				% memorizzo la moving average nella finestra
				moving_average(start : start + x_min(2) - 1) = ma(1 : x_min(2));

				% memorizzo picchi e valli nella finestra
				peaks = [peaks, (start + x_max(1))];
				%valleys = [valleys, (start + x_min(1)), (start + x_min(2))];
				valleys = [valleys, (start + x_min(1))];

				%% 5. Absolute and Relative Artifact Detection in N-1

				% se diventa true, segno tutta la pulse vawe individuata come un
				% artefatto
				artifact_detected = false;
				artifact_type = 0;

				% single raw sample in N-1 pulse wave area marked as disturbed
				disturbed_window = disturbed(start : start + x_min(2) - 1);

				if sum(abs(disturbed_window)) > 0
					%display(unique(disturbed_window));
					artifact_detected = true;
					artifact_type = 1;					% clipped top o clipped bottom
					artifact_triggered(1) = artifact_triggered(1) + 1;
				end

				% check 3: too small relative ampitude detected
				% pwa(n-1) <= 2 * mean(abs(diff(RingBufferSig)))
				if (pwa <= 2* mean(abs(diff(window))))
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 3;
						artifact_triggered(3) = artifact_triggered(3) + 1;
					end
				end

				% check 4: rise time outside absolute range detected
				% 0.08s > pwrt(n-1) >0.49s
				pwrt = (x_max(1) - x_min(1)) / sampling_rate;
				
				rise_time = [rise_time, pwrt];

				if pwrt < 0.08 || pwrt > 0.49
					%fprintf("pwrt: %u\n", pwrt);
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 4;
						artifact_triggered(4) = artifact_triggered(4) + 1;
					end
				end

				% check 5: s/d duration ratio outside absolute range detected
				% pwsdratio(n-1) > 1.1

				% syastolic phase: fase che va dall'inizio dell'onda al picco
				% diastolic phase: fase che va dal picco alla fine dell'onda
				s = (x_max(1) - x_min(1)) / sampling_rate;
				d = (x_min(2) - x_max(1)) / sampling_rate;
				s_to_d = s/d;
				
				s_d_ratio = [s_d_ratio, s_to_d];

				if s_to_d > 1.1
					%fprintf("s/d duration rate: %u\n", s_to_d);
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 5;
						artifact_triggered(5) = artifact_triggered(5) + 1;
					end
				end

				% check 6: duration outside absolute range detected
				% 0.27s > pwd(n-1) > 2.4s
				if (pwd / sampling_rate) < 0.27 || (pwd / sampling_rate) > 2.4
					%fprintf("pwd: %u\n", (pwd/sampling_rate));
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 6;
						artifact_triggered(6) = artifact_triggered(6) + 1;
					end
				end

				% check 7: number of diastolic peaks outside absolute range detected
				% numberofdiastolicpeaks(n-1) > 2

				% cerco i picchi diastolici prendendo una porzione della finestra
				% che va dal picco alla fine dell'onda e cercando al suo interno i
				% massimi relativi
				diastolic_phase = window(x_max(1) : x_min(2));
				
				n_peaks = [n_peaks, sum(islocalmax(diastolic_phase))];

				if sum(islocalmax(diastolic_phase)) > 2
					%fprintf("# massimi relativi: %u\n", sum(islocalmax(diastolic_phase)));
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 7;
						artifact_triggered(7) = artifact_triggered(7) + 1;
					end
				end

				% check 8: not monotonically increasing systolic phase detected
				% in n-1 systolic phase, one sample detected which is not bigger or
				% equal to the sample before
				systolic_phase = window(x_min(1) : x_max(1));
				
				% valuto la percentuale del picco rispetto all'intera altezza dell'onda (pwa)
				pks_x = islocalmax(systolic_phase);
				pks = systolic_phase(pks_x);
				vls_x = islocalmin(systolic_phase);
				vls = systolic_phase(vls_x);
				
				loc_dist = 0;
				
				if length(pks) > 0
					for j = 1 : length(pks)
						s_dist = abs(pks(j) - vls(j));
						
						if s_dist > loc_dist
							loc_dist = s_dist;
						end
					end
				end
				
				systolic_peak = [systolic_peak, loc_dist/pwa];

				for j = 2 : length(systolic_phase)
					if (systolic_phase(j) - systolic_phase(j-1)) < 0
						%fprintf("diff: %u\n", systolic_phase(j) - systolic_phase(j-1));
						artifact_detected = true;
						
						if artifact_type == 0
							artifact_type = 8;
							artifact_triggered(8) = artifact_triggered(8) + 1;
						end
					end
				end

				% check 9: negative valley in diastolic phase detected
				% in n-1 diastolic phase one sample detected wich is smaller than
				% pulse wave n-1 begin or end
				for j = 1 : length(diastolic_phase)
					if (diastolic_phase(j) - window(x_min(1))) < 0 && (diastolic_phase(j) - window(x_min(2))) < 0
						%fprintf("start: %u\n", start);
						%fprintf("diff1: %u, diff2: %u\n", diastolic_phase(j) - window(x_min(1)), diastolic_phase(j) - window(x_min(2)));
						artifact_detected = true;

						if artifact_type == 0
							artifact_type = 9;
							artifact_triggered(9) = artifact_triggered(9) + 1;
						end
					end
				end

				% check 10: pulse distortion detected
				% pwa(n-1 left) / pwa(n-1 right) > 0.4
				% or
				% 0.4 < pwa(n-1 right) / pwa(n-1 left)
				pwa_left = window(x_max(1)) - window(x_min(1));
				pwa_right = window(x_max(1)) - window(x_min(2));

				% TODO: > o < ?
				if (pwa_left / pwa_right) < 0.4 || (pwa_right / pwa_left) < 0.4
					%fprintf("left/right: %u, right/left: %u\n", pwa_left / pwa_right, pwa_right / pwa_left);
					artifact_detected = true;

					if artifact_type == 0
						artifact_type = 10;
						artifact_triggered(10) = artifact_triggered(10) + 1;
					end
				end

				if artifact_detected
					% in pulse wave N-1 marks area from first potential valley to
					% sample before second potential valley as artifact
					artifact(start : start + x_min(2) - 2) = 1;

					% salvo le informazioni dell'artefatto
					art = struct;
					art.start = start;
					art.end = start + x_min(2) - 3;
					art.type = artifact_type;
					
					%display(art);

					artifacts = [artifacts, art];
				else
					% considero l'onda come corretta e salvo le sue informazioni			
					wave = struct;
					wave.start = start;
					wave.end = start + x_min(2) - 3;
					wave.rise_time = pwrt;
					wave.pwd = pwd / sampling_rate;
					wave.pwa = pwa_left;
					%display(wave);

					waves = [waves, wave];

					%% 6. Relative Artifact Detection between N-1 and N-2

					if length(waves) > 2
						% compares relative changes of the last complete pulse wave N-1 with the
						% previous pulse wave N-2
						old_wave = waves(length(waves)-1);

						% pulse wave N-1 and N-2 are direct neighbor
						if (wave.start - old_wave.end) == 1
							% check 11: rise time variation outside relative range detected
							% 33% > pwrt(n-1 to n-2) > 300%
							if wave.rise_time < (old_wave.rise_time / 3) || wave.rise_time > (old_wave.rise_time * 3)
								%fprintf("rise time: %u, rise time old: %u\n", wave.rise_time, old_wave.rise_time);
								artifact_detected = true;

								if artifact_type == 0
									artifact_type = 11;
									artifact_triggered(11) = artifact_triggered(11) + 1;
								end
							end

							% check 12: duration variation outside relative range detected
							% 33% > pwd(n-1 to n-2) > 300%
							if wave.pwd < (old_wave.pwd / 3) || wave.pwd > (old_wave.pwd * 3)
								%fprintf("pwd: %u, pwd old: %u\n", wave.pwd, old_wave.pwd);
								artifact_detected = true;

								if artifact_type == 0
									artifact_type = 12;
									artifact_triggered(12) = artifact_triggered(12) + 1;
								end
							end

							% check 13: ampitude variation outside relative range detected
							% 25% > pwa(n-1 to n-2) > 400%
							if wave.pwa < (old_wave.pwa / 4) || wave.pwa > (old_wave.pwa * 4)
								%fprintf("pwa: %u, èwa old: %u\n", wave.pwa, old_wave.pwa);
								artifact_detected = true;

								if artifact_type == 0
									artifact_type = 13;
									artifact_triggered(13) = artifact_triggered(13) + 1;
								end
							end

							if artifact_detected
								% in pulse wave N-1 marks area from first potential valley to
								% sample before second potential valley as artifact
								artifact(start : start + x_min(2) - 2) = 1;

								% salvo le informazioni dell'artefatto
								art = struct;
								art.start = start;
								art.end = start + x_min(2) - 3;
								art.type = artifact_type;

								artifacts = [artifacts, art];
							end
						end 
					end
				end

				% aggiorno lo start
				start = start + x_min(2) - 2;
				%fprintf("start: %u\n", start);start = start + x_min(2) - 1;
				window_size = 2 * max_pwd;
			else
				if (length(moving_average) - start) <= length(ma)-1
					moving_average(start : (length(moving_average) - 1)) = ma(1 : (length(ma) - 1));

					for j = 1 : length(x_max)
						peaks = [peaks, (start + x_max(j))];
					end

					for j = 1 : length(x_min)
						valleys = [valleys, (start + x_min(j))];
					end

					cont = false;
				else
					%{
					moving_average(start : start + length(ma) - 1) = ma(1 : length(ma));

					display(x_min);
					display(x_max);
					%}

					%start = start + length(ma);
					window_size = window_size + 2 * max_pwd;
				end

				% 4) continua in questo modo fino alla fine
				if moving_average(size(moving_average)) == 0
					cont = false;
				end
			end
		end
		
		save(strcat('data/', files(f).name(1:end-4),'_art.mat'), 'artifacts');
		save(strcat('data/', files(f).name(1:end-4),'_wav.mat'), 'waves');
		
		hold on
		plot(time_axis, moving_average, '-.');
		scatter((peaks - 2)/sampling_rate, high_pass_signal(peaks - 1), 'm+');
		scatter((valleys - 2)/sampling_rate, high_pass_signal(valleys - 1), 'rx');
		% plotto gli artefatti
		%scatter(time_axis(find(artifact == 1)), zeros(length(time_axis(find(artifact == 1))), 1), 'g*');
		x = find(artifact == 1);
		y = 0 * ones(length(x),1);
		scatter((x-1)/sampling_rate, y, 'g*');
		hold off

		base = min(high_pass_signal);
		top = max(high_pass_signal);
		amp = top - base;

		for i = 1 : length(artifacts)
			beg = artifacts(i).start;
			fin = artifacts(i).end;
			type = artifacts(i).type;

			% plotto gli artefatti
			hold on
			rectangle('Position',[(beg/sampling_rate) base ((fin-beg)/sampling_rate) amp]', ...
				%{
				'EdgeColor', [1, 0, 0, 0.5]);
				%}
				'FaceColor', [1, 0, 0, 0.2], ...
				'EdgeColor', [1, 0, 0, 0.5]);
				
			text((beg + 1)/sampling_rate, top - .2, int2str(type));
			%text(475.2867,10,'\leftarrow sin(\pi)')
			hold off	
		end
		
		savefig(strcat('img/', files(f).name(1:end-4),'_hp.fig'));
		%end	
	end		
end

figure(5)
histogram(s_d_ratio, 500);
xlabel('S/D ratio') 
ylabel('Frequency') 
title('Absolute frequency of S/D ratio');
hold on
xline(1.1,'--r');
hold off

figure(6)
histogram(rise_time,300);
xlabel('Rise time [s]') 
ylabel('Frequency') 
title('Absolute frequency of rise time duration');
hold on
xline(0.08,'--r');
hold on
xline(0.49,'--r');
hold off

figure(7)
histogram(systolic_peak);
xlabel('Ampitude of local peak [%]') 
ylabel('Frequency') 
title('Absolute frequency of the ampitude of local peak in systolic phase');
hold off