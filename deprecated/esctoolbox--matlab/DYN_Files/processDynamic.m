% --------------------------------------------------------------------
% function processDynamic
%
% Technical note: PROCESSDYNAMIC assumes that specific Arbin test 
% scripts have been executed to generate the input files. 
% "makeMATfiles.m" converts the raw Excel data files into "MAT" format
% where the MAT files have fields for time, step, current, voltage, 
% chgAh, and disAh for each script run.
%
% The results from three scripts are required at every temperature.  
% The steps in each script file are assumed to be:
%   Script 1 (thermal chamber set to test temperature):
%     Step 1: Rest @ 100% SOC to acclimatize to test temperature
%     Step 2: Discharge @ 1C to reach ca. 90% SOC
%     Step 3: Repeatedly execute dynamic profiles (and possibly
%             intermediate rests) until SOC is around 10%
%   Script 2 (thermal chamber set to 25 degC):
%     Step 1: Rest ca. 10% SOC to acclimatize to 25 degC
%     Step 2: Discharge to min voltage (ca. C/3)
%     Step 3: Rest
%     Step 4: Constant voltage at vmin until current small (ca. C/30)
%     Steps 5-7: Dither around vmin
%     Step 8: Rest
%   Script 3 (thermal chamber set to 25 degC):
%     Step 2: Charge @ 1C to max voltage
%     Step 3: Rest
%     Step 4: Constant voltage at vmax until current small (ca. C/30)
%     Steps 5-7: Dither around vmax
%     Step 8: Rest
% 
% All other steps (if present) are ignored by PROCESSDYNAMIC. The time 
% step between data samples must be uniform -- we assume a 1s sample
% period in this code
%
% The inputs:
% - data: An array, with one entry per temperature to be processed. 
%         One of the array entries must be at 25 degC. The fields of 
%         "data" are: temp (the test temperature), script1, 
%         script 2, and script 3, where the latter comprise data 
%         collected from each script.  The sub-fields of these script 
%         structures that are used by PROCESSDYNAMIC are the vectors: 
%         current, voltage, chgAh, and disAh
% - model: The output from processOCV, comprising the OCV model
% - numpoles: The number of R-C pairs in the model
% - doHyst: 0 if no hysteresis model desired; 1 if hysteresis desired
%
% The output:
% - model: A modified model, which now contains the dynamic fields
%         filled in.

% Copyright (c) 2015 by Gregory L. Plett of the University of Colorado 
% Colorado Springs (UCCS). This work is licensed under a Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 Intl. License, v. 1.0.
% It is provided "as is", without express or implied warranty, for 
% educational and informational purposes only.
%
% This file is provided as a supplement to: Plett, Gregory L., "Battery
% Management Systems, Volume I, Battery Modeling," Artech House, 2015.

function model = processDynamic(data,model,numpoles,doHyst)
  global bestcost
  
  % used by fminbnd later on
  if exist('fminbnd.m','file') % part of optimization toolbox
    options=optimset('TolX',1e-8,'TolFun',1e-8,'MaxFunEval',100000, ...
      'MaxIter',1e6); % for later optimization
  end
    
  % ------------------------------------------------------------------
  % Step 1: Compute capacity and coulombic efficiency for every test
  % ------------------------------------------------------------------
  alltemps = [data(:).temp];
  alletas  = 0*alltemps;
  allQs    = 0*alltemps;
  
  ind25 = find(alltemps == 25); 
  if isempty(ind25),
    error('Must have a test at 25degC');
  end
  not25 = find(alltemps ~= 25);
  
  for k = ind25,    
    totDisAh = data(k).script1.disAh(end) + ...
               data(k).script2.disAh(end) + ...
               data(k).script3.disAh(end);
    totChgAh = data(k).script1.chgAh(end) + ...
               data(k).script2.chgAh(end) + ...
               data(k).script3.chgAh(end);
    eta25 = totDisAh/totChgAh; 
    data(k).eta = eta25; alletas(k) = eta25;
    data(k).script1.chgAh = data(k).script1.chgAh*eta25;
    data(k).script2.chgAh = data(k).script2.chgAh*eta25;
    data(k).script3.chgAh = data(k).script3.chgAh*eta25;    

    Q25 = data(k).script1.disAh(end) + data(k).script2.disAh(end) -...
          data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    data(k).Q = Q25; allQs(k) = Q25;
  end
  eta25 = mean(alletas(ind25));
  
  for k = not25,    
    data(k).script2.chgAh = data(k).script2.chgAh*eta25;
    data(k).script3.chgAh = data(k).script3.chgAh*eta25;
    eta = (data(k).script1.disAh(end) + data(k).script2.disAh(end)+...
           data(k).script3.disAh(end) - data(k).script2.chgAh(end)-...
           data(k).script3.chgAh(end))/data(k).script1.chgAh(end);
    data(k).script1.chgAh = eta*data(k).script1.chgAh;
    data(k).eta = eta; alletas(k) = eta;
    
    Q = data(k).script1.disAh(end) + data(k).script2.disAh(end) - ...
          data(k).script1.chgAh(end) - data(k).script2.chgAh(end);
    data(k).Q = Q; allQs(k) = Q;
  end
  
  model.temps    = unique(alltemps); numTemps = length(model.temps);
  model.etaParam = NaN(1,numTemps);
  model.QParam   = NaN(1,numTemps);
  for k = 1:numTemps,
    model.etaParam(k) = mean(alletas(alltemps == model.temps(k)));
    model.QParam(k)   = mean(allQs(alltemps == model.temps(k)));
  end
  
  % ------------------------------------------------------------------
  % Step 2: Compute OCV for "discharge portion" of test
  % ------------------------------------------------------------------
  for k = 1:length(data),
    etaParam = model.etaParam(k);
    etaik = data(k).script1.current; 
    etaik(etaik<0)= etaParam*etaik(etaik<0);
    data(k).Z = 1 - cumsum([0,etaik(1:end-1)])*1/(data(k).Q*3600); 
    data(k).OCV = OCVfromSOCtemp(data(k).Z(:),alltemps(k),model);
  end
  
  % ------------------------------------------------------------------
  % Step 3: Now, optimize!
  % ------------------------------------------------------------------
  model.GParam  = NaN(1,numTemps); % "gamma" hysteresis parameter
  model.M0Param = NaN(1,numTemps); % "M0" hysteresis parameter
  model.MParam  = NaN(1,numTemps); % "M" hysteresis parameter
  model.R0Param = NaN(1,numTemps); % "R0" ohmic resistance parameter
  model.RCParam = NaN(numTemps,numpoles); % time const.
  model.RParam  = NaN(numTemps,numpoles); % Rk

  for theTemp = 1:numTemps, 
    fprintf('Processing temperature %d\n',model.temps(theTemp));
    bestcost = Inf;
%     theGammas = 1:5000;
%     theFit = zeros(size(theGammas));
%     for indGamma = 1:length(theGammas),
%       theFit(indGamma) = optfn(theGammas(indGamma),data,model,...
%                                model.temps(theTemp),doHyst);
%     end
%     figure(4); plot(theFit); stop
    if doHyst,
      if exist('fminbnd.m','file'),
        model.GParam(theTemp) = abs(fminbnd(@(x) optfn(x,data,...
                                  model,model.temps(theTemp),...
                                  doHyst),1,250,options));
      else
        model.GParam(theTemp) = abs(gss(@(x) optfn(x,data,...
                                  model,model.temps(theTemp),...
                                  doHyst),1,250,1e-8));
      end
    else
      model.GParam(theTemp) = 0;
      optfn(theGParam,data,model,model.temps(theTemp),doHyst);
    end
    [~,model] = minfn(data,model,model.temps(theTemp),doHyst);                          
  end
return

% --------------------------------------------------------------------
% This minfn works for the enhanced self-correcting cell model
% --------------------------------------------------------------------
function cost=optfn(theGParam,data,model,theTemp,doHyst)
  global bestcost 
  
  model.GParam(model.temps == theTemp) = abs(theGParam);
  [cost,model] = minfn(data,model,theTemp,doHyst);
  if cost<bestcost, % update plot of model params for every improvement
    bestcost = cost;
    disp('Best ESC model values yet!');
    figure(3); theXlim = [min(model.temps) max(model.temps)];
    subplot(2,4,1); plot(model.temps,model.QParam); 
                    title('Capacity (Ah)'); 
                    xlim(theXlim);
    subplot(2,4,2); plot(model.temps,1000*model.R0Param); 
                    title('Resistance (m\Omega)');
                    xlim(theXlim);
    subplot(2,4,3); plot(model.temps,1000*model.M0Param); 
                    title('Hyst Magnitude M0 (mV)');
                    xlim(theXlim);
    subplot(2,4,4); plot(model.temps,1000*model.MParam); 
                    title('Hyst Magnitude M (mV)');
                    xlim(theXlim);
    subplot(2,4,5); plot(model.temps,getParamESC('RCParam',...
                    model.temps,model));
                    title('RC Time Constant (tau)');
                    xlim(theXlim);
    subplot(2,4,6); plot(model.temps,1000*getParamESC('RParam',...
                    model.temps,model));
                    title('R in RC (m\Omega)');
                    xlim(theXlim);
    subplot(2,4,7); plot(model.temps,abs(model.GParam)); 
                    title('Gamma');    
                    xlim(theXlim);
  end
return

% --------------------------------------------------------------------
% Using an assumed value for gamma (already stored in the model), find 
% optimum values for remaining cell parameters, and compute the RMS 
% error between true and predicted cell voltage
% --------------------------------------------------------------------
function [cost,model]=minfn(data,model,theTemp,doHyst)
  alltemps = [data(:).temp];
  ind = find(alltemps == theTemp); numfiles = length(ind);

  xplots = ceil(sqrt(numfiles));
  yplots = ceil(numfiles/xplots);
  rmserr = zeros(1,xplots*yplots);
  
  G = abs(getParamESC('GParam',theTemp,model));
  Q = abs(getParamESC('QParam',theTemp,model));
  eta = abs(getParamESC('etaParam',theTemp,model));
  RC = getParamESC('RCParam',theTemp,model);
  numpoles = length(RC);
  
  for thefile = 1:numfiles;
    ik = data(ind(thefile)).script1.current(:);
    vk = data(ind(thefile)).script1.voltage(:);
    tk = (1:length(vk))-1;
    etaik = ik; etaik(ik<0) = etaik(ik<0)*eta;

    h=0*ik; sik = 0*ik;
    fac=exp(-abs(G*etaik/(3600*Q)));
    for k=2:length(ik),
      h(k)=fac(k-1)*h(k-1)-(1-fac(k-1))*sign(ik(k-1));
      sik(k) = sign(ik(k));
      if abs(ik(k))<Q/100, sik(k) = sik(k-1); end
    end
    
    % First modeling step: Compute error with model = OCV only
    vest1 = data(ind(thefile)).OCV;
    verr = vk - vest1;
    
    % Second modeling step: Compute time constants in "A" matrix
    np = numpoles; 
    while 1,
      A = SISOsubid(-diff(verr),diff(etaik),np);
      eigA = eig(A); 
      eigA = eigA(eigA == conj(eigA));  % make sure real
      eigA = eigA(eigA > 0 & eigA < 1); % make sure in range
      okpoles = length(eigA); np = np+1;
      if okpoles >= numpoles, break; end
      fprintf('Trying np = %d\n',np);
    end    
    RCfact = sort(eigA); RCfact = RCfact(end-numpoles+1:end);
    RC = -1./log(RCfact);
    % Simulate the R-C filters to find R-C currents
    if exist('dlsim.m','file') % in the control-system toolbox
      vrcRaw = dlsim(diag(RCfact),1-RCfact,...
                   eye(numpoles),zeros(numpoles,1),etaik);
    else % a somewhat slower workaround if no control-system toolbox
      vrcRaw = zeros(length(RCfact),length(etaik));
      for vrcK = 1:length(etaik)-1,
        vrcRaw(:,vrcK+1) = diag(RCfact)*vrcRaw(:,vrcK)+(1-RCfact)*etaik(vrcK);
      end
      vrcRaw = vrcRaw';
    end

    % Third modeling step: Hysteresis parameters
    if doHyst,
      H = [h,sik,-etaik,-vrcRaw]; 
      if exist('lsqnonneg.m','file'), % in optimization toolbox
        W = lsqnonneg(H,verr); %  W = H\verr;    
      else
        W = nnls(H,verr); %  W = H\verr;    
      end
      M = W(1); M0 = W(2); R0 = W(3); Rfact = W(4:end)';
    else
      H = [-etaik,-vrcRaw]; 
      W = H\verr;    
      M=0; M0=0; R0 = W(1); Rfact = W(2:end)';
    end
    ind = find(model.temps == data(ind(thefile)).temp,1);
    model.R0Param(ind) = R0;
    model.M0Param(ind) = M0;
    model.MParam(ind) = M;
    model.RCParam(ind,:) = RC';
    model.RParam(ind,:) = Rfact';
    
    vest2 = vest1 + M*h + M0*sik - R0*etaik - vrcRaw*Rfact';
    verr = vk - vest2;
    
    % plot voltages: decimate to make faster
    figure(1); subplot(yplots,xplots,thefile); 
    plot(tk(1:10:end)/60,vk(1:10:end),tk(1:10:end)/60,...
         vest1(1:10:end),tk(1:10:end)/60,vest2(1:10:end));  
    xlabel('Time (min)'); ylabel('Voltage (V)'); 
    title(sprintf('Voltage and estimates at T=%d',...
                  data(ind(thefile)).temp));
    legend('voltage','vest1 (OCV)','vest2 (DYN)','location','southwest');

    % plot modeling errors: decimate to make faster
    figure(2); subplot(yplots,xplots,thefile); 
    thetitle=sprintf('Modeling error at T = %d',data(ind(thefile)).temp);
    plot(tk(1:10:end)/60,verr(1:10:end)); title(thetitle);
    xlabel('Time (min)'); ylabel('Error (V)');
    ylim([-0.1 0.1]); 
    drawnow
    
    % Compute RMS error only on data roughly in 5% to 95% SOC
    v1 = OCVfromSOCtemp(0.95,data(ind(thefile)).temp,model);
    v2 = OCVfromSOCtemp(0.05,data(ind(thefile)).temp,model);
    N1 = find(vk<v1,1,'first'); N2 = find(vk<v2,1,'first');
    if isempty(N1), N1=1; end; if isempty(N2), N2=length(verr); end
    rmserr(thefile)=sqrt(mean(verr(N1:N2).^2));    
  end 

  cost=sum(rmserr); 
  fprintf('RMS error = %0.2f (mV)\n',cost*1000);
  if isnan(cost), stop, end
return

% A = SISOsubid(y,u,n);
%  Identifies state-space "A" matrix from input-output data.
%     y: vector of measured outputs
%     u: vector of measured inputs 
%     n: number of poles in solution
%           
%     A: discrete-time state-space state-transition matrix.
%                 
%  Theory from "Subspace Identification for Linear Systems
%               Theory - Implementation - Applications" 
%               Peter Van Overschee / Bart De Moor (VODM)
%               Kluwer Academic Publishers, 1996
%               Combined algorithm: Figure 4.8 page 131 (robust)
%               Robust implementation: Figure 6.1 page 169
%
%  Code adapted from "subid.m" in "Subspace Identification for 
%               Linear Systems" toolbox on MATLAB CENTRAL file 
%               exchange, originally by Peter Van Overschee, Dec. 1995
function A = SISOsubid(y,u,n)
  y = y(:); y = y'; ny = length(y); % turn y into row vector
  u = u(:); u = u'; nu = length(u); % turn u into row vector
  i = 2*n; % #rows in Hankel matrices. Typically: i = 2 * (max order)
  twoi = 4*n;           

  if ny ~= nu, error('y and u must be same size'); end
  if ((ny-twoi+1) < twoi); error('Not enough data points'); end

  % Determine the number of columns in the Hankel matrices
  j = ny-twoi+1;

  % Make Hankel matrices Y and U
  Y=zeros(twoi,j); U=zeros(twoi,j);
  for k=1:2*i
    Y(k,:)=y(k:k+j-1); U(k,:)=u(k:k+j-1);
  end
  % Compute the R factor
  R = triu(qr([U;Y]'))'; % R factor
  R = R(1:4*i,1:4*i); 	 % Truncate

  % ------------------------------------------------------------------
  % STEP 1: Calculate oblique and orthogonal projections
  % ------------------------------------------------------------------
  Rf = R(3*i+1:4*i,:);              % Future outputs
  Rp = [R(1:1*i,:);R(2*i+1:3*i,:)]; % Past inputs and outputs
  Ru  = R(1*i+1:2*i,1:twoi); 	      % Future inputs
  % Perpendicular future outputs 
  Rfp = [Rf(:,1:twoi) - (Rf(:,1:twoi)/Ru)*Ru,Rf(:,twoi+1:4*i)]; 
  % Perpendicular past inputs and outputs
  Rpp = [Rp(:,1:twoi) - (Rp(:,1:twoi)/Ru)*Ru,Rp(:,twoi+1:4*i)]; 

  % The oblique projection is computed as (6.1) in VODM, page 166.
  % obl/Ufp = Yf/Ufp * pinv(Wp/Ufp) * (Wp/Ufp)
  % The extra projection on Ufp (Uf perpendicular) tends to give 
  % better numerical conditioning (see algo on VODM page 131)

  % Funny rank check (SVD takes too long)
  % This check is needed to avoid rank deficiency warnings
  if (norm(Rpp(:,3*i-2:3*i),'fro')) < 1e-10
    Ob = (Rfp*pinv(Rpp')')*Rp; 	% Oblique projection
  else
    Ob = (Rfp/Rpp)*Rp;
  end

  % ------------------------------------------------------------------
  % STEP 2: Compute weighted oblique projection and its SVD
  %         Extra projection of Ob on Uf perpendicular
  % ------------------------------------------------------------------
  WOW = [Ob(:,1:twoi) - (Ob(:,1:twoi)/Ru)*Ru,Ob(:,twoi+1:4*i)];
  [U,S,~] = svd(WOW);
  ss = diag(S);

  % ------------------------------------------------------------------
  % STEP 3: Partitioning U into U1 and U2 (the latter is not used)
  % ------------------------------------------------------------------
  U1 = U(:,1:n); % Determine U1

  % ------------------------------------------------------------------
  % STEP 4: Determine gam = Gamma(i) and gamm = Gamma(i-1) 
  % ------------------------------------------------------------------
  gam  = U1*diag(sqrt(ss(1:n)));
  gamm = gam(1:(i-1),:);
  gam_inv  = pinv(gam); 			% Pseudo inverse of gam
  gamm_inv = pinv(gamm); 			% Pseudo inverse of gamm

  % ------------------------------------------------------------------
  % STEP 5: Determine A matrix (also C, which is not used) 
  % ------------------------------------------------------------------
  Rhs = [gam_inv*R(3*i+1:4*i,1:3*i),zeros(n,1); R(i+1:twoi,1:3*i+1)];
  Lhs = [gamm_inv*R(3*i+1+1:4*i,1:3*i+1); R(3*i+1:3*i+1,1:3*i+1)];
  sol = Lhs/Rhs;    % Solve least squares for [A;C]
  A = sol(1:n,1:n); % Extract A
return

function X = gss(f,a,b,tol)
  % golden section search to find the minimum of f on [a,b]
  % based on code: https://en.wikipedia.org/wiki/Golden-section_search
  gr = (sqrt(5)+1)/2; % golden ratio used in search

  c = b - (b - a) / gr;
  d = a + (b - a) / gr;
  while abs(c - d) > tol
    if f(c) < f(d),
      b = d;
    else
      a = c;
    end

    % we recompute both c and d here to avoid loss of precision which 
    % may lead to incorrect results or infinite loop
    c = b - (b - a) / gr;
    d = a + (b - a) / gr;
  end
  X = (b+a)/2;
return

function [x,w,info]=nnls(C,d,opts)
  % nnls  Non negative least squares Cx=d x>=0 w=C'(d-Cx)<=0
  %  2012-08-21  Matlab8  W.Whiten
  %  2013-02-17  Line 52 added
  %  Copyright (C) 2012, W.Whiten (personal W.Whiten@uq.edu.au) BSD license
  %  (http://opensource.org/licenses/BSD-3-Clause)
  %
  % [x,w,info]=nnls(C,d,opts)
  %  C    Coefficient matrix
  %  d    Rhs vector
  %  opts Struct containing options: (optional)
  %        .Accy  0 fast version, 1 refines final value (default), 
  %                 2 uses accurate steps but very slow on large cases, 
  %                 faster on small cases, result usually identical to 1
  %        .Order True or [], or order to initially include positive terms
  %                 if included will supply info.Order, if x0 available use 
  %                 find(x0>0), but best saved from previous run of nnls
  %        .Tol   Tolerance test value, default zero, use multiple of eps
  %        .Iter  Maximum number of iterations, should not be needed.
  %
  %  x    Positive solution vector x>=0
  %  w    Lagrange multiplier vector w(x==0)<= approx zero
  %  info Struct with extra information: 
  %        .iter  Number of iterations used
  %        .wsc0  Estimated size of errors in w
  %        .wsc   Maximum of test values for w
  %        .Order Order variables used, use to restart nnls with opts.Order
  %
  % Exits with x>=0 and w<= zero or slightly above 0 due to
  %  rounding and to ensure for convergence
  % Using faster matrix operations then refines answer as default (Accy 1).
  % Accy 0 is more robust in singular cases.
  %
  % Follows Lawson & Hanson, Solving Least Squares Problems, Ch 23.

  [~,n]=size(C);
  maxiter=4*n;

  % inital values
  P=false(n,1);
  x=zeros(n,1);
  z=x;

  w=C'*d;

  % wsc_ are scales for errors
  wsc0=sqrt(sum(w.^2));
  wsc=zeros(n,1);
  tol=3*eps;
  accy=1;
  pn1=0;
  pn2=0;
  pn=zeros(1,n);

  % see if option values have been given
  ind=true;
  if(nargin>2)
    if(isfield(opts,'Tol'))
      tol=opts.Tol;
      wsc(:)=wsc0*tol;
    end
    if(isfield(opts,'Accy'))
      accy=opts.Accy;
    end
    if(isfield(opts,'Iter'))
      maxiter=opts.Iter;
    end
  end

  % test if to use normal matrix for speed
  if(accy<2)
    A=C'*C;
    b=C'*d;
    %L=zeros(n,n);
    LL=zeros(0,0);
    lowtri=struct('LT',true);
    uptri=struct('UT',true);
  end

  % test if initial information given
  if(nargin>2)
    if(isfield(opts,'Order') && ~islogical(opts.Order))
      pn1=length(opts.Order);
      pn(1:pn1)=opts.Order;
      P(pn(1:pn1))=true;
      ind=false;
    end
    if(~ind && accy<2)
      %L(1:pn1,1:pn1)=chol(A(pn(1:pn1),pn(1:pn1)),'lower');
      UU(1:pn1,1:pn1)=chol(A(pn(1:pn1),pn(1:pn1)));
      LL=UU';
    end
    pn2=pn1;
  end

  % loop until all positive variables added
  iter=0;
  while(true)
    % Check if no more terms to be added
    if(ind && (all(P==true) || all(w(~P)<=wsc(~P))))
      if(accy~=1)
        break
      end
      accy=2;
      ind=false;
    end

    % skip if first time and initial Order given
    if(ind)
      % select best term to add
      ind1=find(~P);
      [~,ind2]=max(w(ind1)-wsc(ind1));
      ind1=ind1(ind2);
      P(ind1)=true;
      pn2=pn1+1;
      pn(pn2)=ind1;
    end

    % loop until all negative terms are removed
    while(true)

      % check for divergence
      iter=iter+1;
      if(iter>=2*n)
        if(iter>maxiter)
          error(['nnls Failed to converge in ' num2str(iter)  ...
              ' iterations'])
        elseif(mod(iter,n)==0)
          wsc=(wsc+wsc0*tol)*2;
        end
      end

      % solve using suspected positive terms
      z(:)=0;
      if(accy>=2)
        z(P)=C(:,P)\d;
      else
        % add row to the lower triangular factor
        for i=pn1+1:pn2
          i1=i-1;
          %LL=L(1:i1,1:i1);
          %LL=LL(1:i1,1:i1);
          t=linsolve(LL,A(pn(1:i1),pn(i)),lowtri);
          %t=LL\A(pn(1:i1),pn(i));
          %L(i,1:i1)=t;
          %LL(i,1:i1)=t;
          AA=A(pn(i),pn(i));
          tt=AA-t'*t;
          if(tt<=AA*tol)
              tt=1e300;
          else
              tt=sqrt(tt);
          end
          %L(i,i)=sqrt(tt);
          %LL(i,i)=sqrt(tt);
          LL(i,1:i)=[t',tt];
          UU(1:i,i)=[t;tt];
        end

        % solve using lower triangular factor
        %LL=L(1:pn2,1:pn2);
        t=linsolve(LL,b(pn(1:pn2)),lowtri);
        %t=LL\b(pn(1:pn2));
        %UU=LL';
        %z(pn(1:pn2))=linsolve(UU,t,uptri);
        z(pn(1:pn2))=linsolve(UU,t,uptri);
        %z(pn(1:pn2))=LL'\t;
        % or could use this to solve without updating factors
        %z(pn(1:pn2))=A(pn(1:pn2),pn(1:pn2))\b(pn(1:pn2));
      end
      pn1=pn2;

      % check terms are positive
      if(all(z(P)>=0))
        x=z;
        if(accy<2)
          w=b-A*x;
        else
          w=C'*(d-C*x);
        end
        wsc(P)=max(wsc(P),2*abs(w(P)));
        ind=true;
        break
      end

      % select and remove worst negative term
      ind1=find(z<0);
      [alpha,ind2]=min(x(ind1)./(x(ind1)-z(ind1)+realmin));
      ind1=ind1(ind2);

      % test if removing last added, increase wsc to avoid loop
      if(x(ind1)==0 && ind)
        w=C'*(d-C*z);
        wsc(ind1)=(abs(w(ind1))+wsc(ind1))*2;
      end
      P(ind1)=false;
      x=x-alpha*(x-z);
      pn1=find(pn==ind1);
      pn(pn1:end)=[pn(pn1+1:end),0];
      pn1=pn1-1;
      pn2=pn2-1;
      if(accy<2)
        LL=LL(1:pn1,1:pn1);
        UU=UU(1:pn1,1:pn1);
      end
      ind=true;
    end
  end

  % info result required
  if(nargout>2)
    info.iter=iter;
    info.wsc0=wsc0*eps;
    info.wsc=max(wsc);
    if(nargin>2 && isfield(opts,'Order'))
      info.Order=pn(1:pn1);
    end
  end

return

