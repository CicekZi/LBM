%% inputting from user
Uinitial = 0.1;	% Uinitial in lattice Units
Vinitial = 0.0 ;	% Vinitial in lattice Units
Rhoinitial = 1.0;	% Rhoinitial in lattice Units

fprintf('choose your method... (or press enter for defualt: 1)\n');
fprintf('1.bounceback\n');
fprintf('2.flippova\n');
fprintf('3.mei\n');
fprintf('4.bozidi\n');
fprintf('5.yu\n');


fprintf('\n>>>');
replymethode = input('');
if isempty(replymethode)
    replymethode = 1 ;
end
methode = replymethode;

Tau = 0.8 ;		% Tau - Relaxation time
fprintf('\nTau = ? (default: %d)\n',Tau);
replyTau = input('');
if isempty(replyTau)
    replyTau = Tau ;
end
Tau = replyTau;

R = 11;
fprintf('\nR = ? (default: %d)\n',R);
replyR = input('');
if isempty(replyR)
    replyR = R ;
end
R = replyR;

if methode==2 || methode==3 || methode==5 %%wall rotation available only in two methods
    
    wall_rotation=0;
    
    fprintf('wall rotation = ? (default: %d)\n', wall_rotation);
    replywall_rotation = input('');
    
    if isempty(replywall_rotation)
        replywall_rotation = wall_rotation ;
    end
    
    wall_rotation = replywall_rotation;
end

nx = 600;
fprintf('nx = ? (default: %d)\n',nx);
replynx = input('');
if isempty(replynx)
    replynx = nx ;
end
nx = replynx;

ny = 150;
fprintf('ny = ? (default: %d)\n',ny);
replyny = input('');
if isempty(replyny)
    replyny = ny ;
end
ny = replyny;


Nu_physical= 1e-3;  % Nu = 0.003 
fprintf('Nu_physical = ? (default: %f)\n',Nu_physical);
replyNu_physical = input('');
if isempty(replyNu_physical)
    replyNu_physical = Nu_physical ;
end
Nu_physical = replyNu_physical;


channel_height=0.01; %0.01 m
fprintf('channel_height = ? (default: %f)\n',channel_height);
replychannel_height = input('');
if isempty(replychannel_height)
    replychannel_height = Uinitial ;
end
channel_height = replychannel_height;


fprintf('U initial = ? (default: %f)\n',Uinitial);
replyu = input('');
if isempty(replyu)
    replyu = Uinitial ;
end
Uinitial = replyu;



Nu = (Tau-0.5)/3 ;  %lattice kinematic viscosity
Re=Uinitial*2*ny/Nu;
Re_cylinder=Uinitial*2*R/Nu;
FD=R*Rhoinitial*Uinitial^2/105.6430/Re_cylinder;
t_lattice=channel_height^2/ny^2/3*(Tau-1/2)/Nu_physical;
dt=t_lattice;
dx=channel_height/ny;
fprintf('*****************************************\n')
fprintf('               Re of flow  = %f              \n',Re)
fprintf('               Re Cylinder = %f              \n',Re_cylinder)
fprintf('               FD          = %f              \n',FD)
fprintf('               Tau         = %f              \n',Tau)
fprintf('               time lattice= %f              \n',t_lattice)
fprintf('*****************************************\n')
reply = input('Please review flow paramters, continue? Y/N [Y]', 's');
if isempty(reply)
    reply = 'Y';
end
if reply == 'N' || reply == 'n'
    return
end