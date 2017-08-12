%% Simulate model with gravity override
% We consider a model consisting of two zones of different permeability. A
% light fluid of high mobility is injected into the lower zone by a
% vertical well placed near the east side. Fluids are produced from a well
% placed near the west side and perforated in the lower zone only.
mrstModule add incomp coarsegrid ad-core


%% Make grid and assign petrophysical properties
% The grid has two zones with different permeabilities: high permeability
% on top and low below, or opposite if inequality sign is reversed in the
% definition of layer function
G = cartGrid([50,1,40],[1500 1000 200]);
%G = cartGrid([60,1,10],[1500 10 200]);
G.nodes.coords(:,3) = G.nodes.coords(:,3)+2050;
G = computeGeometry(G);

%figure(1); clf
%plotGrid(G); view(3); axis tight

[K1,K2,p1,p2] = deal(10,100,.3,.3);
layer = @(c) (c(:,3)-2150)>0; % <0: high perm on top, >0: low on top

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
T  = 10*day();
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

%src = addSource([], 1, .8*sum(poreVolume(G,rock))/T,'sat',[1 0]);
ci = linspace(1000,2000,21);
src = addSource([], ci, repmat(0.3,numel(ci),1), 'sat', [1,0]);


CG = generateCoarseGrid(G,(layer(G.cells.centroids)>0)+1);
plotFaces(CG,1:CG.faces.num,'FaceColor','none','LineWidth',1);
plotFaces(CG,11,'FaceColor','y','FaceAlpha',.3);


%% Fluid model
% Inject a light and mobilt fluid into a denser and less mobile fluid
fluid = initSimpleFluid('mu' , [  .1,   1] .* centi*poise     , ...
                        'rho', [ 700,1000] .* kilogram/meter^3, ...
                        'n'  , [   2,   2]);               %    ,...
%                        'sr' , [0, 0.3]);
                    
% boundary condition  
bc = psideh([], G, 'xmax', fluid, 'sat', [0 1]);

%% Simulation loop
N  = 10;
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

% Plot saturation
   Sco2Matrix = reshape(rSol.s(:,1),[50, 40])' ;
   Sco2Matrix = fliplr(Sco2Matrix);
   %SbrineMatrix = reshape(rSol.s(:,1),[40, 50]);
   SbrineMatrix = reshape(rSol.s(:,2),[50, 40])' ;
   
   imagesc(Sco2Matrix);
   colorbar;
   %colormap jet;
   %delete(hs);
   %hs = plotCellData(G, rSol.s(:,1), (rSol.s(:,1)>.01), pargs{:});
   %view(45, 35);
   %hs = plotCellData(G, rSol.s(:,1));
   %title([num2str(convertTo(t,day)), 'days']),
   %caxis([0 1]); drawnow 

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