# matlab
Tasks 1 in Matlab

close all
type='Sinc';
fr=10;
fsampl=100;
sig=0.01;
x=-5:1/fsampl:5;%1 fsampl+1 data points as middle point at 0 
L=length(x);
xmid=round(length(x)/2,0);

switch type
   case 'Cos'
      Signal=(0*exp(-fr*1i*2*pi*x)+1*exp(fr*1i*2*pi*x)); 
   case 'Gauss'
      Signal=1/(4*sqrt(2*pi*sig))*(exp(-2*x.^2/(2*sig)));
   case 'Sinc'
      Signal=((-1*exp(-fr*pi*1i*x)+1*exp(fr*1i*pi*x))./(pi*2*x*1i*fr));
      Signal(1,xmid)=1;%correction for Sinc(t=0) is NaN!
   case 'SincSQ'
      tempSignal=((-1*exp(-fr*pi*1i*x)+1*exp(fr*1i*pi*x))./(pi*2*x*1i*fr));
      Signal=(tempSignal).^2;
      Signal(1,xmid)=1;%correction for Sinc(t=0) is NaN!
end


%Do the Fast Fourier Transform and Discrete Fourier Transform
FFT=fft(Signal)/L;
DFT=dft(Signal)/L;

%Shifting x-Axis
freq=(0:(L-1))*(fsampl/L);
freqshift=(-(L-1)/2:(L)/2)*(fsampl/L);

%Rearranging the data for x-Axis
shiftedFFT=fftshift(FFT);
shiftedDFT=fftshift(DFT);

%FWHM and TBP
fwhmIn=fwhm(abs(Signal),x);
fwhmDFT=fwhm(abs(shiftedDFT),freqshift);
fwhmFFT=fwhm(abs(shiftedFFT),freqshift);
TBPdft=fwhmIn*fwhmDFT
TBPfft=fwhmIn*fwhmFFT;

%Parsevals Theorem
absSignal=Signal.*conj(Signal);
absDFT=DFT.*conj(FFT);
ParsevalcompareSignaltoDFT=[sum(absSignal)/L,sum(absDFT)]




%Plot
figure;
subplot(2,1,1)
p0=plot(x,Signal);
title(type)
xlabel('"time"')
ylabel('Amplitude')

subplot(2,1,2)
p1=plot(freqshift,abs(shiftedFFT));
hold on
p2=plot(freqshift,abs(shiftedDFT),'--');
hold off
set([p0,p1,p2],'LineWidth',1.2)
title(['FT of ' type])
xlabel('"frequency"')
ylabel('Amplitude')
legend(['FFT' TBPfft],['DFT' TBPdft])

function mat = createTransfMat(L)
%Creates the matrix needed for diskrete fourier transformation
%mat(k,l)=LxL with elements W^k*l
mat=zeros(L,L);
    for k=1:L
        for l=1:L
            mat(k,l)=exp(-1i*2*pi*(k-1)*(l-1)/L);
        end
    end
end

function vec = dft(signal)
%Applies Discrete Fourier Transformation
%F[f]=M*f
L=length(signal);
vec=signal*createTransfMat(L);
end


function fnc=fwhm(data,axis)
%Finds FWHM
[pks,locs,widths,prom] = findpeaks(data,axis,'NPeaks',1,'SortStr','descend','WidthReference','halfheight');
fnc=widths;
 end

