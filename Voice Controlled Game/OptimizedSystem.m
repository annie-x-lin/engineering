clear; clc; clf; format compact;
disp("Optmized voice controlled game. Press any key to start")
pause

% initialization
tic % reset clock
recObj = audiorecorder(8000,8,1);
str = 'Optimized system ';
Nx=400;
transient=800;
recTime=0.15;
% load filter coefs from Lec 19 HW for optimized system
load('a_filter','FilterCoefs')
a_coefs = FilterCoefs;
load('i_filter','FilterCoefs')
i_coefs = FilterCoefs;

% load Hmag_dB for /a/ and /i/ from Lec 17 HW
Ha=getHmag(Nx,a_coefs);
Hi=getHmag(Nx,i_coefs);

      % length of signal x and H_mag
% define Hamming window of length Nx
hw = 0.54 - 0.46*cos(2*pi*(0:Nx-1)/(Nx-1));

% Add code to define and normalize condition templates Ta and Ti. Normalize
% with respect to mean and STD. 
Ta = Ha(1:Nx/2+1);
Ti = Hi(1:Nx/2+1);
Ta = Ta - mean(Ta); Ta = Ta/std(Ta);
Ti = Ti - mean(Ti); Ti = Ti/std(Ti);

nGames = 10; % increase to 10 after debugged
nWins = 0;
dx = 2;  % me_x and targetx increment, increasing dx makes game faster (bigger steps)  dx >= 1
dy = 3;  % me_y control increment
initTime = toc  % get initialization time
tic % reset clock
for game = 1:nGames
    % initial locations
    targetx = 100;%right side this is the x axis, how fast it moves
    targety = randi([1,70]);  % on right side at random height
    me_x = 1;
    me_y = randi([1,70]); % on left side at random height
    
    for n = 2:dx:100 %dx is how fast things move
        targetx = 100-n;
        me_x = n;
        me_y = OptimizedSys(me_y,dy,recObj,recTime,Ta,Ti,hw, transient);
        [done, win, result] = displayGame(targetx,targety,me_x,me_y);
        if done == 1
            nWins = nWins + win;
            fprintf('Game %d / %d %s\n',game,nGames,result)
            break
        end
    end
end
fprintf('%s  %d wins in %d games\n',str,nWins,nGames)
gTime = toc;
fprintf('Total time= %.2f seconds\n',gTime)

function x = getSpeech(recObj,recTime,Nx, transient)
    tau = 0.001;
    recordblocking(recObj, recTime);    %recTime is the variable we have to figure out. Its in the code 
    signal = getaudiodata(recObj)';
    x = signal(transient:end); % eliminates start transient for laptop mic
   
    if max(abs(x(1:100))) < tau % no speech at beginning
        x = zeros(1,Nx);
        return % exit function
    end
    nx = length(x);
    if nx < Nx
        x = [x zeros(1,Nx-nx)];
    else
        x = x(1:Nx);
    end
end

function me_y = OptimizedSys(me_y,dy,recObj,recTime,Ta,Ti,hw, transient)
myData = getSpeech(recObj,recTime,length(hw), transient);
x = myData;
nx = length(fft(hw.*x));
Xmag = abs(fft(hw.*x));
S_dB = 20*log10(Xmag);
S_dB = S_dB - max(S_dB);
S = S_dB(1:nx/2+1);
S = S - mean(S);
S = S/std(S);
Va = sum(Ta.*S);
Vi = sum(Ti.*S);
if Vi > Va
    me_y = me_y + dy;
elseif abs(Vi-Va)<45
    me_y = me_y + dy;    
elseif Va-Vi>45
    me_y = me_y - dy;
end
end

function [done,win,result] = displayGame(tx,ty,mx,my)
done = 0;
win = 0;
result = 'Playing';
figure(1)
plot(tx,ty,'ro','markersize',25)
axis([1 100 1 70])
axis off
hold on
plot(mx,my,'ko')
hold off
if mx >= tx % reached target
    done = 1;
    result = 'LOSS';
    if abs(my-ty) < 4
        result = 'WIN';
        win = 1;
    end
    text(10,35,result)
    pause(1)
end
end

function Hmag = getHmag(Nx,a_coefs)
hw = 0.54 - 0.46*cos(2*pi*(0:Nx-1)/(Nx-1));
d= [1 zeros(1,Nx-1)];
a_F1=a_coefs(1,:);
a_F2=a_coefs(2,:);
a_F3=a_coefs(3,:);

%F1
h = AR(a_F1,d,Nx);
%F-1
h = AR(a_F2,h,Nx);
%F3
h = AR(a_F3,h,Nx);
x = hw.*h(1:Nx);
X = fft(x);
Xmag = abs(X);
Hmag_dB = 20*log10(Xmag);
Hmag_dB = Hmag_dB - max(Hmag_dB); %normalize to 0 dB
Hmag = Hmag_dB;
end

function y = AR(a,x,ny)
ybuf = zeros(size(a));  % output buffer memory
y = zeros(1,ny); % length of AR output is infinite, finite number must be specified.
for i = 1:length(y)  
    if i > length(x)
        xin = 0; 
    else
        xin = x(i);
    end
    y(i) = sum(a.*ybuf) + xin;  % sum of products
    % shift values in xbuf by one unit
    if length(ybuf) >= 2
        ybuf(2:length(ybuf)) = ybuf(1:length(ybuf)-1);
    end
    ybuf(1) = y(i);    
end
end
