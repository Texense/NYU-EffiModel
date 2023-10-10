function [spectrogram]=spectrogram(sd,param)

%unit of timewindow and bin is milisecond while unit of frequency is times/second.
tw=param.spectrogram_timewindow; % Integration Time Win
bin=param.sdbin;                 % Binsize??
fr=param.frequency_range;        % frequency range [min max]
%duration = param.duration;
initial_index = 1;
N=tw/bin;
kernel=(exp(-2*pi*1i*bin/1000)*ones(1,N)).^(1:N)';
num_spec= length(sd) - N+1 - initial_index+1;

% Build up kernal mat
kernelMat = kernel.^(fr(1):1:fr(2));
% Build up sliding signal matrix
SdMat = zeros(num_spec,N);
for SlideInd=1:num_spec
    SdMat(SlideInd,:) = sd(SlideInd + initial_index-1: SlideInd+initial_index+N-2);
end
spectrogram = (abs(bin/1000 * (SdMat * kernelMat)/sqrt(tw/1000)).^2)';
% for j=1:num_spec
%     for k=1:fr(2)-fr(1)+1
%         spectrogram(k,j)=abs(bin/1000 * sd(initial_index+j -1: initial_index+j+N-2)*(kernel.^(k+fr(1)-1))...
%             /sqrt(tw/1000))^2;
%     end
% end

%smoothing
norm(spectrogram,'fro');
spectrogram1=conv2(spectrogram,ones(1,8)/8, 'valid');
frs   = linspace(fr(1), fr(2), fr(2)-fr(1)+1);
times = 0: num_spec;
times = times *bin;
imagesc(times, frs, spectrogram1);
set(gca,'YDir','normal')%, 'Position',[0.1, 0.41, 0.77, 0.24]);
%set(gca,'ColorScale','log')
xlabel('Time(ms)','fontsize',11);
ylabel('Freq(Hz)','fontsize',11);
hcb=colorbar;
%set(hcb, 'YScale','log')%,'Position',[0.89,0.41, 0.03, 0.24]);
ylabel(hcb,'(spikes/sec)^2/Hz','fontsize',11);
%yticks([20 40 60 80]);
%xticks([2000, 2200, 2400, 2600, 2800, 3000]);
set(gca,'fontsize',11);
%set(gca,'xtick',[]);
colormap turbo;
%caxis([0 300]);
%title(name);
end




         
    
    