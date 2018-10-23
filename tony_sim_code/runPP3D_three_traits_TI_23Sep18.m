% PP-evolution, 3D

% -----------------
% Adaptive dynamics version
% -----------------


q = 3; % number of traits
mutProb = 0.01; % mutation rate per time step, per density unit
SDmut = .03; %.05; %.1; %0.07;%0.05; % SD of mutational steps

timeT = 6e6;


for iter = 1
	figure(300+iter), clf reset
	figure(400+iter), clf reset
	drawnow
	if iter < 3
		name = sprintf('PP3Dsim1_q3_%d',iter)
		nameStart = sprintf('Start_PP3Dsim1_q3_%d',iter)

		%initV = 0.1;%.517;
		%initU = 0.9;%0.85;
		V = rand(1,q);%initV*ones(1,q);
		U = rand(1,q);%initU*ones(1,q);

		initP = 1;%1/B*(1-A*M/cc/B); % predator equilibrium population size
		initN = 1;%M/cc/B; % prey equilibrium

		N = initN + zeros(size(V,1),1);
		P = initP; % a single predator

		time_samples = 0:1000:timeT;
		nt = length(time_samples);
		% create an empty strucure array with fields 'V' and 'N', one element per sample:
		Vsamples = struct('V',cell(1,nt),'N',cell(1,nt));
		Usamples = struct('U',cell(1,nt),'P',cell(1,nt));

		%sample time 0:
		Vsamples(1).V = V;
		Vsamples(1).N = N;
		Usamples(1).U = U;
		Usamples(1).P = P;

		tic
		save_interval = 20; % save (and plot) every 20 secs
		nextsave = save_interval;
		last_plotted_sample = 0;
		% run dynamics with mutations:
		for sample = 2:length(time_samples)
			[sample size(N,1) size(V,1) size(P,1) size(U,1)]

			[V,N,U,P,WN,WP] = doManyTimeSteps_TI_23Sep18(V,N,U,P,diff(time_samples(sample-1:sample)),q,mutProb,SDmut);
			Vsamples(sample).V = V;
			Vsamples(sample).N = N;
			Usamples(sample).U = U;
			Usamples(sample).P = P;

			if toc>nextsave || sample==length(time_samples)
				save(name,'Vsamples','Usamples','time_samples','sample');
				nextsave = toc + save_interval;
				set(0,'currentfigure',300+iter)

				if (last_plotted_sample + 1) < sample

					ssi = (last_plotted_sample + 1):10:sample;
					for di=1:q
						subplot(q,1,di)
						for si = ssi
							plot(time_samples(si+zeros(size(Vsamples(si).V,1),1)), Vsamples(si).V(:,di), '.b')
							hold on
						end
						xlim([0 timeT])
						ylim([0 2])
						xlabel('time')
						ylabel(['V_' num2str(di)])
						drawnow
					end
					set(0,'currentfigure',400+iter)
					for di=1:q
						subplot(q,1,di)
						for si=ssi
							plot(time_samples(si+zeros(size(Usamples(si).U,1),1)), Usamples(si).U(:,di), '.b')
							hold on
						end
						xlim([0 timeT])
						ylim([0 2])
						xlabel('time')
						ylabel(['U_' num2str(di)])
						hold on
						drawnow
					end
				end
				last_plotted_sample = ssi(end)+9;
			end
		end
		save(nameStart,'V','U','N','P');
	end

	if iter == 3

		timeT = .54e6

		name = sprintf('PP3Dsim1_q3_%d',iter)
		nameStart = sprintf('Start_PP3Dsim1_q3_%d',iter)

		load('Start_PP3Dsim1_q3_1')
		NStart = N;
		PStart = P;
		VStart = V;
		UStart = U;

		load('Start_PP3Dsim1_q3_2')
		N = [NStart;N];
		P = [PStart;P];
		V = [VStart;V];
		U = [UStart;U];

% 		[xN, orderN]=sort(N);
% 		orderN=flip(orderN);
% 		pick=orderN(1:round(length(N)/10));
% 		N=N(pick);
% 		V=V(pick,:);
%
		[xP, orderP]=sort(P);
		orderP=flip(orderP);
		pick=orderP(1:round(length(P)/2));
		initP=P(pick);
		initU=U(pick,:);

 		N = .5*N;
 		P = .5*P;

		time_samples = 0:1000:timeT;
		nt = length(time_samples);
		% create an empty strucure array with fields 'V' and 'N', one element per sample:
		Vsamples = struct('V',cell(1,nt),'N',cell(1,nt));
		Usamples = struct('U',cell(1,nt),'P',cell(1,nt));

		%sample time 0:
		Vsamples(1).V = V;
		Vsamples(1).N = N;
		Usamples(1).U = U;
		Usamples(1).P = P;

		tic
		save_interval = 20; % save (and plot) every 20 secs
		nextsave = save_interval;
		last_plotted_sample = 0;
		% run dynamics with mutations:
		for sample = 2:length(time_samples)
			[sample size(N,1) size(V,1) size(P,1) size(U,1)]

			[V,N,U,P,WN,WP] = doManyTimeSteps_TI_21Sep18(V,N,U,P,diff(time_samples(sample-1:sample)),q,mutProb,SDmut);

			%%%%%%%%%%%%%%%%%%%
			if length(P) < length(initP)/10 & sample < 40
				P = initP;
				U = initU;
			end


			Vsamples(sample).V = V;
			Vsamples(sample).N = N;
			Usamples(sample).U = U;
			Usamples(sample).P = P;

			if toc>nextsave || sample==length(time_samples)
				save(name,'Vsamples','Usamples','time_samples','sample');
				nextsave = toc + save_interval;
				set(0,'currentfigure',300+iter)

				if (last_plotted_sample + 1) < sample

					ssi = (last_plotted_sample + 1):10:sample;
					for di=1:q
						subplot(q,1,di)
						for si = ssi
							plot(time_samples(si+zeros(size(Vsamples(si).V,1),1)), Vsamples(si).V(:,di), '.b')
							hold on
						end
						xlim([0 timeT])
						ylim([0 2])
						xlabel('time')
						ylabel(['V_' num2str(di)])
						drawnow
					end
					set(0,'currentfigure',400+iter)
					for di=1:q
						subplot(q,1,di)
						for si=ssi
							plot(time_samples(si+zeros(size(Usamples(si).U,1),1)), Usamples(si).U(:,di), '.b')
							hold on
						end
						xlim([0 timeT])
						ylim([0 2])
						xlabel('time')
						ylabel(['U_' num2str(di)])
						hold on
						drawnow
					end
				end
				last_plotted_sample = ssi(end)+9;
			end
		end
		save(nameStart,'V','U','N','P');
	end
end
figure(500+iter)
plot3(V(:,1), V(:,2),V(:,3),'.b')
