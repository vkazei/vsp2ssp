% Creates data at the surface from VSP data ------------------------------
% (c) Vladimir Kazei, Oleg Ovcharenko and Yan Yang, KAUST 2019

%% This program demonstrates how to redatum the DAS VSP data to virtual SSP data by crosscorrelation interferometry
close all; clear all;


zsrc = 100; %% z coordinate index of sources 
xsrc = 70; %% x coordinate index of sources

nsrc = 25; dsrc = 4; %% nsrc--number of sources; dsrc--sources interval in grid points;
recordVideoFlag = 1;

% receivers on a horizontal line at depth
rec_depth = 3;
model_ghost_flag = 1;

%% MODEL
% Model dimensions
nx = 370; nz = 200; %% nx-- numbers of x coordinate index, nz--numbers of z coordinate index
dx = 10;    % [m] dx--interval of x coordinate index

% Velocity
vp = 2000.0 * ones(nz, nx); %constant velocity 

% velocity of compressional waves, [m/s]

%% TIME STEPPING
t_total = 4.0; % [sec] recording duration
dt = 2e-3; %[sec] time interval

%% SOURCE
f0 = 10.0;                          % dominant frequency of the wavelet
t0 = 1.2 / f0;                     % half Ricker wavelet excitation time
source_ampl = 1e10;                      % amplitude coefficient

x = xsrc;                 % source location along OX

fsFlag = 0;
dt2 = dt^2;

upscale = 4; % decimates data for correlations
min_wavelengh = 0.5*min(vp(vp>330))/f0;     % shortest wavelength bounded by velocity in the air


%% ABSORBING BOUNDARY (ABS)
abs_thick = 50;%min(floor(0.15*nx), floor(0.15*nz)); % thicknes of the layer



abs_rate = 0.3/abs_thick;                           % decay rate
lmargin = [abs_thick abs_thick];
rmargin = [abs_thick abs_thick];
weights = ones(nz+2,nx+2);
IT_DISPLAY = 10; % display every .. snapshot for test source location
zsrcArr = dsrc*[1:nsrc] + abs_thick;

opts = v2struct();

%% test modeling code
imagesc(single_shot(opts, zsrc))
caxis(caxis()/50);

%% generate VSP data
opts.IT_DISPLAY = 100000; % don't display anything
opts.fsFlag = 1;


tic;
parfor i=1:nsrc
    data_direct_t_src_rec(:,i,:) = imresize(single_shot(opts, zsrcArr(i)),1/upscale);
end
disp('DIRECT WAVE DATA GENERATED');
toc;

opts.vp(abs_thick+100:end,:) = 2500;
opts.fsFlag = 0;
tic;
parfor i=1:nsrc
    data_ghost_t_src_rec(:,i,:) = imresize(single_shot(opts, zsrcArr(i)),1/upscale);
end
disp('GHOST NO FS DATA GENERATED');
toc;

opts.fsFlag = 1;
tic;
parfor i=1:nsrc
    data_t_src_rec(:,i,:) = imresize(single_shot(opts, zsrcArr(i)),1/upscale);
end
disp('FULL DATA GENERATED');
toc;

save('TRUE_DATA');

%% correlating
figure();
newRecord = 1;
load('TRUE_DATA');
nt = size(data_t_src_rec,1);
t = linspace(0,t_total,nt);

data_new_t_s_r = zeros(2*size(data_t_src_rec,1)-1,size(data_t_src_rec,3),size(data_t_src_rec,3));
last_source =0;

%if ~exist('vido','var')
    
    filename = 'virtualData.gif';
   
%end

for new_sou = round(abs_thick/upscale):size(data_t_src_rec,3)-round(abs_thick/upscale)
    for i = [1:size(data_t_src_rec,3) * dx * upscale]
        t_src(i) = timeRefl(dx*new_sou*upscale, ...
            i, 1000, 2000);
    end
    for new_rec = max(new_sou,round(abs_thick/upscale)):size(data_t_src_rec,3)-round(abs_thick/upscale)
        traces_direct = squeeze(data_direct_t_src_rec(:,:,new_sou));
        traces_full = squeeze(data_t_src_rec(:,:,new_rec));
        traces_refl = traces_full-data_direct_t_src_rec(:,:,new_rec);
        
        for real_rec = 1:size(data_t_src_rec,2)
            data_new_t_s_r(:, new_sou, new_rec) = data_new_t_s_r(:, new_sou, new_rec) + ...
                xcorr(nmz(traces_refl(:,real_rec)),nmz(traces_direct(:,real_rec)));
        end
        if mod(new_rec, 10) == 0 && mod(new_sou, 10) == 0
            h = figure(555);
            subplot 131
            imagesc(t,[],[traces_direct 10*traces_refl ]');% ...
            title('Direct from */ no direct from v')
            caxis([-100 100])
            hold on 
            plot(t, nsrc+0.7*nsrc* ...
                nmz(data_new_t_s_r(size(data_t_src_rec,1):end, new_sou, new_rec)), ...
                'Color','r','LineWidth',2);
            legend('virtual trace')
            ylabel('Receiver depth (m)')
            yticks([5 : 5 : 2*nsrc])
            yticklabels(dx*dsrc*[[5:5:nsrc],[5:5:nsrc]])
            xlabel('time (sec)')
            
            
            subplot 132
            imagesc(upscale*dx*[-1:size(data_new_t_s_r,3)-1],t, ...
                squeeze(data_new_t_s_r(size(data_t_src_rec,1):end, ...
                new_sou, :)));
            %caxis(cax);
            caxis auto
            caxis(caxis()/5);
            title('New shot gather for *')
            hold on
            
            plot(t_src,'LineWidth',2,'Color','r');
            
            
            
            SP_well = opts.xsrc/upscale;
            line(upscale*dx*[SP_well SP_well],get(gca,'YLim'),'Color','black','LineWidth',3)
            
            SP=SP_well+3*(new_sou-SP_well);
            line(upscale*dx*[SP SP],get(gca,'YLim'),'Color','g','LineWidth',3)
            
            line(upscale*dx*[new_sou new_sou],get(gca,'YLim'),'Color','r','LineWidth',1)
            
            line(upscale*dx*[new_rec new_rec],get(gca,'YLim'),'Color','blue','LineWidth',1)
            
            legend('target reflection time', 'well location', 'ray coverage zone', 'shot', 'virtual receiver');
            xlabel('x (m)')
            ylabel('time (sec)')
            
            subplot 133
            
            imagesc(imresize(opts.vp,10,'nearest'));
            caxis([1000 3000])
            hold on
    
            scatter(dx*xsrc*ones(size(zsrcArr)), dx*zsrcArr, 20, 'black', 'v');
            scatter(dx*new_sou*upscale, dx*(opts.rec_depth + abs_thick), 40, 'r', '*');            
            scatter(dx*new_rec*upscale, dx*(opts.rec_depth + abs_thick), 40, 'b', 'v', 'filled');
            
            legend('real receivers', 'real shot', 'virtual receiver')
            
            %axis equal tight
            title('Velocity')
            xlabel('x (m)')
            ylabel('z (m)')
            
            colorbar
            
            drawnow
            set (gcf,'units','pixels','position',[1 1 1400 450])
            pause(0.01);
            
            if recordVideoFlag
                
                % Capture the plot as an image
                frame = getframe(gcf);
                im = frame2im(frame);
                [imind,cm] = rgb2ind(im,256);
                % Write to the GIF File
                if newRecord == 1
                    imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                    newRecord = 0;
                else
                    imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',1);
                end
                
            end
            
        end
    end
    
    
end


%%
% normalizes the data and corrects geometrical spreading
function nmzz = nmz(a)
nmzz = a/(max(abs(a(:)))+10^-30);% geometrical spreading -> .* (1:length(a))'.^0.5 ;
end

function tr = timeRefl(xs, xr, refl_depth, vel)

tr = sqrt((xs-xr).^2 + 4*refl_depth^2)/vel;

end

%% puts source to depth zsrc (kept se and returns seismogram at opts.z
function seismogram = single_shot(opts, zsrc)
% --------------------------------------------------------------
% based on Oleg Ovcharenko and Vladimir Kazei, 2018
% https://github.com/ovcharenkoo/WaveProp_in_MATLAB

% Output every ... time steps
opts.zsrc = zsrc;
v2struct(opts);
rec_depth = rec_depth + abs_thick;

dz = dx;    % [m]



%dt = min(dt, 0.5 * min(dx,dz)/max(vp(:)));   % min grid space / max velocity
nt = ceil(t_total/dt);             % number of time steps

CFL = max(vp(:)) * dt / min(dx,dz);

t = [0:nt]*dt;
t_n = pi*f0*(t-t0);   % normalized time for Ricker
source_term = (1.0 - 2.0*(t_n.^2)).*exp(-t_n.^2);        % Ricker source time function (second derivative of a Gaussian):
source_term = source_ampl * source_term * dt2 / (dx * dz);


for iz = 1:nz+2
    for ix = 1:nx+2
        i = 0;
        j = 0;
        k = 0;
        if (ix < lmargin(1) + 1)
            i = lmargin(1) + 1 - ix;
        end
        if (iz < lmargin(2) + 1)
            k = lmargin(2) + 1 - iz;
            if fsFlag
                k = 10^20;
            end
            
        end
        if (nx - rmargin(1) < ix)
            i = ix - nx + rmargin(1);
        end
        if (nz - rmargin(2) < iz)
            k = iz - nz + rmargin(2);
        end
        if (i == 0 && j == 0 && k == 0)
            continue
        end
        rr = abs_rate * abs_rate * double(i*i + j*j + k*k );
        weights(iz,ix) = exp(-rr);
    end
end

% %% SUMMARY
% fprintf('#################################################\n');
% fprintf('2D acoustic FDTD wave propagation in isotripic \nmedium in displacement formulation with \nCerjan(1985) boundary conditions\n');
% fprintf('#################################################\n');
% fprintf('Model:\n\t%d x %d\tgrid nz x nx\n\t%.1e x %.1e\t[m] dz x dx\n',nz, nx, dz,dx);
% fprintf('\t%.1e x %.1e\t[m] model size\n',nx*dx, nz*dz);
% fprintf('\t%.1e...%.1e\t[m/s] vp\n', min(vp(:)), max(vp(:)));
% fprintf('Time:\n\t%.1e\t[sec] total\n\t%.1e\tdt\n\t%d\ttime steps\n',t_total,dt,nt);
% fprintf('Source:\n\t%.1e\t[Hz] dominant frequency\n\t%.1f\t[sec] index time\n',f0,t0);
% fprintf('Other:\n\t%.1f\tCFL number\n', CFL);
% fprintf('\t%.2f\t[m] shortest wavelength\n\t%d, %d\t points-per-wavelength OX, OZ\n', min_wavelengh, floor(min_wavelengh/dx), floor(min_wavelengh/dz));
% fprintf('#################################################\n');

%% ALLOCATE MEMORY FOR WAVEFIELD
p3 = zeros(nz+2,nx+2);            % Wavefields at t
p2 = zeros(nz+2,nx+2);            % Wavefields at t-1
p1 = zeros(nz+2,nx+2);            % Wavefields at t-2
% Coefficients for derivatives
co_dxx = 1/dx^2;
co_dzz = 1/dz^2;

%% Loop over TIME
%tic;
for it = 1:nt
    p3 = zeros(size(p2));
    % Second-order derivatives
    dp_dxx = co_dxx * (p2(2:end-1,1:end-2) - 2*p2(2:end-1,2:end-1) + p2(2:end-1,3:end));
    dp_dzz = co_dzz * (p2(1:end-2,2:end-1) - 2*p2(2:end-1,2:end-1) + p2(3:end,2:end-1));
    % U(t) = 2*U(t-1) - U(t-2) + G dt2/rho;
    p3(2:end-1,2:end-1) = 2.0*p2(2:end-1,2:end-1) - p1(2:end-1,2:end-1) + (vp.^2).*(dp_dxx + dp_dzz).*dt2;
    % Add source term
    p3(zsrc, xsrc) = p3(zsrc, xsrc) + source_term(it);
    % Exchange data between t-2 (1), t-1 (2) and t (3) and apply ABS
    p1 = p2 .* weights;
    p2 = p3 .* weights;
    % Output
    if mod(it,IT_DISPLAY) == 0
        fprintf('Time step: %d \t %.4f s\n',it, single(t(it)));
        imagesc(p3); colorbar;
        axis equal tight; colormap jet;
        caxis([-100 100])
        title(['Step = ',num2str(it),'/',num2str(nt),', Time: ',sprintf('%.4f',t(it)),' sec']);
        drawnow;
        pause(0.1);
    end
    seismogram(it,:) = p3(rec_depth,:);
    if model_ghost_flag
        seismogram(it,:) = seismogram(it,:) - p3(rec_depth-2*opts.rec_depth,:);
    end
end
%toc; disp('End');

end
