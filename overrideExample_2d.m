%% Simulate model with gravity override
% We consider a model consisting of two zones of different permeability. A
% light fluid of high mobility is injected into the lower zone by a
% vertical well placed near the east side. Fluids are produced from a well
% placed near the west side and perforated in the lower zone only.
mrstModule add incomp coarsegrid ad-core

% time step size from 1 day to 0.1 day, nothing change
% vertical refinement from 40 mesh to 100 mesh, nothing change
% dx from 15m to 5m, nothing change


%% Make grid and assign petrophysical properties
% The grid has two zones with different permeabilities: high permeability
% on top and low below, or opposite if inequality sign is reversed in the
% definition of layer function
G = cartGrid([50,1,40],[2000 100 100]);
%G = cartGrid([60,1,10],[1500 10 200]);
G.nodes.coords(:,3) = G.nodes.coords(:,3)+2050;
G = computeGeometry(G);

%figure(1); clf
%plotGrid(G); view(3); axis tight

[K1,K2,p1,p2] = deal(10,100,.25,.25);
%layer = @(c) (c(:,3)-2150)>0; % <0: high perm on top, >0: low on top
layer = @(c) (c(:,3)-2100)>0; % <0: high perm on top, >0: low on top

rock = makeRock(G, K1*milli*darcy, p1);
rock.poro(layer(G.cells.centroids)) = p2;
rock.perm(layer(G.cells.centroids)) = K2*milli*darcy;

hT = computeTrans(G,rock);

clf
pargs = {'EdgeAlpha',.1,'EdgeColor','k'};
hs = plotCellData(G,rock.perm,pargs{:});
view(3), axis tight
set(hs,'FaceAlpha',.35);
zoom(1.4); set(gca,'dataasp',[2 2 1]); view(27,-12);

%% Setup wells
% Both wells are perforated in the lower zone only.
% No wells are included in this test
T  = (365*2)*day();
% x = G.cells.centroids(:,1:2);
% W = addWell([], G, rock,...
%             find( sum(bsxfun(@minus,x,[67.5 487.5]).^2,2)<320 ...
%                   & G.cells.centroids(:,3)>2150),  ...
%             'InnerProduct', 'ip_tpf', ...
%             'Type', 'bhp', 'Val', 0*barsa, ...
%             'Comp_i', [0 1], 'Name', 'P', 'Dir','z');

% x = G.cells.centroids(:,1:2);
% 
% W = addWell([],  G, rock, ...
%             find( sum(bsxfun(@minus,x,[1437.5 487.5]).^2,2)<320 ...
%                   & G.cells.centroids(:,3)>2150),  ...
%             'InnerProduct', 'ip_tpf',...
%             'Type', 'rate', 'Val', 0*.8*sum(poreVolume(G,rock))/T, ...
%             'Comp_i', [1 0], 'Name', 'I', 'Dir','z');

%plotWell(G,W,'height',75,'radius',.01);
% plotWell(G,W,'height',75,'radius',10000);

% src = addSource([], 1, .8*sum(poreVolume(G,rock))/T,'sat',[1 0]);
% add the source of CO2 at the bottom layer
ci = linspace(1000,2000,21);
% this injection rate is for every cell/ in c++ code is in the whole domain
src = addSource([], ci, repmat(42.33*meter^3/day,numel(ci),1), 'sat', [1,0]);


CG = generateCoarseGrid(G,(layer(G.cells.centroids)>0)+1);
plotFaces(CG,1:CG.faces.num,'FaceColor','none','LineWidth',1);
plotFaces(CG,11,'FaceColor','y','FaceAlpha',.3);


%% Fluid model
% Inject a light and mobilt fluid into a denser and less mobile fluid
fluid = initSimpleFluid('mu' , [  0.0726,   0.72] .* centi*poise     , ...
                        'rho', [ 810,1000] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);               %    ,...
%                        'sr' , [0, 0.3]);
                    
% boundary condition  
bc = psideh([], G, 'xmax', fluid, 'sat', [0 1]);

%% Simulation loop
N  = 365*2;
dT = T/N*ones(N,1);
dT = [dT(1)*sort(2.^-[1:4 4])'; dT(2:end)];

gravity reset on
rSol = initState(G, [], 0, [0, 1]);
rSol = incompTPFA(rSol, G, hT, fluid, 'src', src);

t = 0; 
colormap(flipud(winter))
%wellSols = cell(numel(dT),1);
set(gca,'XTick',[],'Ytick',[],'ZTick',[]);
%
for i=1:numel(dT)
   rSol = implicitTransport(rSol, G, dT(i), rock, fluid, 'src', src);

   % Check for inconsistent saturations
   assert(max(rSol.s(:,1)) < 1+eps && min(rSol.s(:,1)) > -eps);

   % Update solution of pressure equation.
   rSol  = incompTPFA(rSol , G, hT, fluid, 'src', src);

   % Measure water saturation in production cells in saturation
   %wellSols{i} = getWellSol(W, rSol, fluid);

   % Increase time
   t = t + dT(i);

   % show the current status
   X = ['the current time step is', i];
   disp(i);
%    % Plot saturation
%    delete(hs);
%    hs = plotCellData(G, rSol.s(:,1), (rSol.s(:,1)>.01), pargs{:});
%    view(45, 35);
%    %hs = plotCellData(G, rSol.s(:,1));
%    title([num2str(convertTo(t,day)), 'days']),
%    caxis([0 1]); drawnow 

end

% plot saturation
Sco2Matrix = reshape(rSol.s(:,1),[50, 40])' ;
Sco2Matrix = fliplr(Sco2Matrix);
%SbrineMatrix = reshape(rSol.s(:,1),[40, 50]);
SbrineMatrix = reshape(rSol.s(:,2),[50, 40])' ;

dx = 5;
dz = 0.5;

x = 20:40:1980;
z = 1.25:2.5:(100-1.25);

figure(1);
Sco2Matrix = flipud(Sco2Matrix);
h1 = uimagesc(x,z,Sco2Matrix);
colormap(jet);
% shading flat; 
% shading interp;
% caxis([ 700 850]); colorbar;
caxis([0.0 1.0]);colorbar('ticks',0:0.1:1.0);
% colorbar('ytick',0:0.1:1);
% xlim([0 2000]);ylim([0 50]);
xlabel('x [m]','fontsize',18);
ylabel('z [m]','fontsize',18);
set(gca, 'fontsize',18,'ylim',[0 100],'xlim',[0 2000],'ytick',0:20:100);
set(gca,'YDir','normal');
fig = gcf;
fig.PaperPosition = [0 2 7.5 3];
%saveas(gcf,'MRST2D_reference_solution.pdf');

% Plot saturation
%    Sco2Matrix = reshape(rSol.s(:,1),[50, 40])' ;
%    Sco2Matrix = fliplr(Sco2Matrix);
%    %SbrineMatrix = reshape(rSol.s(:,1),[40, 50]);
%    SbrineMatrix = reshape(rSol.s(:,2),[50, 40])' ;
%    
%    imagesc(Sco2Matrix);
   colorbar;
   %colormap jet;
   
%% Oil rate, with peak production indicated by red line
% figure,
% [Ym,Tm] = meshgrid(G.cells.centroids(W(1).cells,3),cumsum(dT)/year);
% p       = cellfun(@(x) abs(x(1).qO)', wellSols,'UniformOutput',false);
% po      = vertcat(p{:})*day;
% [~,j]   = max(po);
% m       = sub2ind(size(po),j,1:numel(W(1).cells));
% surf(Ym,Tm,po); shading interp; colormap(parula(20));
% hold on; plot3(Ym(m),Tm(m),po(m), '-r','LineWidth',2); hold off
% axis tight, view(100,35)

%% Water production, with breakthrough indicated by red line
% figure,
% p  = cellfun(@(x) abs(x(1).qW)', wellSols,'UniformOutput',false);
% pw = vertcat(p{:})*day;
% j  = sum(~(pw>1e-3));
% m  = sub2ind(size(pw),j,1:numel(W(1).cells));
% surf(Ym,Tm,pw); shading interp; colormap(parula(20));
% hold on; plot3(Ym(m),Tm(m),pw(m), '-r','LineWidth',2); hold off
% axis tight, view(45,35)

%% Total rate, with breakthrough indicated by red line
% figure,
% p = po + pw;
% surf(Ym,Tm,p); shading interp; colormap(parula(20));
% hold on; plot3(Ym(m),Tm(m),p(m), '-r','LineWidth',2); hold off
% axis tight, view(45,35)

%% Plot surface rates etc using GUI
% plotWellSols(wellSols, cumsum(dT));