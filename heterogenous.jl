Pkg.add("JuMP")
Pkg.add("Gurobi")
using JuMP;
using Gurobi;

N = 12;
# functions
function getChildren(par)
	child = Array{Int64, 1}[];
	N = length(par);
	for i=1:N
		childi = Int64[];
		for j=1:N
			if par[j] == i
				push!(childi, j);
			end
		end
		push!(child, childi);
	end

	child
end

function countEdgesOnPath(par)
	nEdgesOnPath = Int64[];
	for i = 1:length(par)
		node = i;
		j = 1;
		while par[node] !=0
			node = par[node];
			j+=1;
		end
		push!(nEdgesOnPath, j);
	end
	nEdgesOnPath
end

function countCommonEdges(par)
	paths = Array{Int64,2}[];
	for i = 1:N
		path = Int64[];
		node = i;
		while node != 0;
			push!(path,node);
			node = par[node];
		end
		push!(paths, path');
	end

	commonEdgesList = Array{Int64,2}[];
	for i = 1:N
		commonEdges = zeros(Int64, N);
		p1 = paths[i];


		for j = 1:N
			k = 0;
			p2 = paths[j];
            l1 = length(p1);
			l2 = length(p2);
			while l1 > 0 && l2 > 0 && p1[l1] == p2[l2]
				k += 1;
				l1 -= 1;
				l2 -= 1;
			end
            commonEdges[j] = k;
		end
        push!(commonEdgesList, commonEdges');
	end

    commonEdgesList
end

function countCommonResistance(par)
	paths = Array{Int64,2}[];
  commonResistances = Array(Float64,N,N);
  commonReactances = Array(Float64,N,N);

	for i = 1:N
		path = Int64[];
		node = i;
		while node != 0;
			push!(path,node);
			node = par[node];
		end
		push!(paths, path');
	end

	for i = 1:N
		p1 = paths[i];
		for j = 1:N
			resistance = 0;
      reactance = 0;
			p2 = paths[j];
      l1 = length(p1);
			l2 = length(p2);
			while l1 > 0 && l2 > 0 && p1[l1] == p2[l2]
				resistance += r[i];
        reactance += x[i];
				l1 -= 1;
				l2 -= 1;
			end
      commonResistances[i,j] = resistance;
      commonReactances[i,j] = reactance;
		end
	end

  commonResistances, commonReactances
end
a = Array(Float64,N,N);

function getBasicModel(gammaMin)
  m = Model(solver=GurobiSolver());
  @defVar(m, l >= 0);
  @defVar(m, v >= 0);
  # Binary Variables
  @defVar(m, deltas[1:N], Bin)
  # continous variables
  @defVar(m, gammaMin <= gammas[1:N] <= 1)
  # voltage and current magnitude squares
  @defVar(m, nus[1:N] >= 0)
  @defVar(m, lij[1:N] >= 0)
  # real and reactive and complex loads and DGs
  @defVar(m, pc[1:N] >= 0)
  @defVar(m, qc[1:N] )
  @defVar(m, pg[1:N] >= 0)
  @defVar(m, qg[1:N] )
  # real and reactive power flows
  @defVar(m, P[1:N] )
  @defVar(m, Q[1:N] )
  #auxiliary variable
  @defVar(m, t >= 0)

  m, l, v, deltas, gammas, nus, lij, pc, qc, pg, qg, P, Q, t
end

function getLinearModel(Clovr, gammaMin)
  m, l, v, deltas, gammas, nus, lij, pc, qc, pg, qg, P, Q, t = getBasicModel(gammaMin);
	for i=1:N

		# load constraints
		addConstraint(m, pc[i] == gammas[i]*PCD[i]);
		addConstraint(m, qc[i] == gammas[i]*QCD[i]);

		# PV constraints
    addConstraint(m, pg[i] <= deltas[i]*SGmax[i]);
    addConstraint(m, qg[i] <= -SGmax[i] + 2 * SGmax[i] * deltas[i]);
    addConstraint(m, pg[i]^2  + qg[i]^2 <= SGmax[i]^2);
    addConstraint(m, qg[i] >= -SGmax[i]);

		# power conservation
		PSUM = pc[i] - pg[i];
		QSUM = qc[i] - qg[i];

		for j in child[i]
			PSUM += P[j];
			QSUM += Q[j];
		end

		addConstraint(m, P[i] == PSUM);
		addConstraint(m, Q[i] == QSUM);

		#ohms law
		eqSum = -2*(r[i]*P[i] + x[i]*Q[i]) ;

		if i == 1
			eqSum += nus0;
		else
			eqSum += nus[par[i]];
		end

		addConstraint(m, nus[i] == eqSum);
		addConstraint(m, t >= nusL - nus[i]);

	end

	# define objective
	addConstraint(m, l == Clovr*t);
  obj = l;
  setObjective(m, :Min, l);
	m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg
end
function getNonLinearModel(Cll, Clovr, gammaMin)
  m, l, v, deltas, gammas, nus, lij, pc, qc, pg, qg, P, Q, t = getBasicModel(gammaMin);

	for i=1:N
		# load constraints
		addConstraint(m, pc[i] == gammas[i]*PCD[i]);
		addConstraint(m, qc[i] == gammas[i]*QCD[i]);

    # PV constraints
    addConstraint(m, pg[i] <= deltas[i]*SGmax[i]);
    addConstraint(m, qg[i] <= -SGmax[i] + 2 * SGmax[i] * deltas[i]);
    addConstraint(m, qg[i] >= -SGmax[i]);
    addConstraint(m, pg[i]^2 + qg[i]^2 <= SGmax[i]^2);

		# power conservation
		PSUM = pc[i] - pg[i] + r[i] * lij[i];
		QSUM = qc[i] - qg[i] + x[i] * lij[i];

		for j in child[i]
			PSUM += P[j];
			QSUM += Q[j];
		end

		addConstraint(m, P[i] == PSUM);
		addConstraint(m, Q[i] == QSUM);

		#ohms law
		eqSum = -2*(r[i]*P[i] + x[i]*Q[i]) + (r[i]^2 + x[i]^2) * lij[i];

		if i == 1
			eqSum += nus0;
		else
			eqSum += nus[par[i]];
		end

		addConstraint(m, nus[i] == eqSum);
		addConstraint(m, t >= nusL - nus[i]);

		#current constraint
		addConstraint(m, lij[i] * nus[i] >= P[i]^2 + Q[i]^2);

	end

	# define objective
	voll = 0;
  obj = 0;
	for i=1:N
    voll += Cll * (1-gammas[i]) * PCD[i];
		obj += Cll * (1-gammas[i]) * PCD[i] + r[i]*lij[i];
	end
  addConstraint(m, v == voll);
  addConstraint(m, l == Clovr*t);

	obj += Clovr*t;

  setObjective(m, :Min, obj);
	m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg
end

ny = 8N+1;
function getRBar(Cll, Clovr)
  ny = 8N+1;
	diameter = maximum(nEdgesOnPath);
	R = zeros(Float64, ny);
  for i = 1:ny
    if i <= N
      # pg
      R[i] = 2 * nEdgesOnPath[i] * r * Clovr  +  Cll;
    elseif i <= 2N
      # qg
      R[i] = 2 * nEdgesOnPath[i-N] * x * Clovr +  Cll; # qg
    elseif i <= 3N
      # pc
      #R[i] = -2 * nEdgesOnPath[i-2N] * r * Clovr;
      R[i] = 0;
    elseif i <= 4N
      # qc
      # R[i] = -2 * nEdgesOnPath[i-3N] * x * Clovr;
      R[i] = 0;
    elseif i <= 5N
      # P
      # R[i] = -2 * diameter * r * Clovr;
      R[i] = 0;
    elseif i <= 6N
      # Q
      # R[i] = -2 * diameter * x * Clovr;
      R[i] = 0;
    elseif i <= 7N
      # nu
      # R[i] = Clovr;
      R[i] = 0;
    elseif i <= 8N
      # gamma
      R[i] = Cll*PCD[i-7N] - 2 * (r * PCD[i-7N] + x * QCD[i-7N]) * Clovr;
      R[i] = Cll;
    else
      # t
      #R[i] = -Clovr;
      R[i] = 0;
    end
  end
	R
end


function getADLP2Model(Cll, Clovr, ndeltas, rbar)
  ny = 8N+1;
	m = Model(solver=GurobiSolver());

	# Binary Variables
	@defVar(m, deltas[1:ny], Bin)

	# continous variables
  @defVar(m, gammaMin <= gammas[1:N] <= 1)

	# voltage and current magnitude squares
	@defVar(m, nus[1:N] >= 0)
	#@defVar(m, lij[1:N] >= 0)

	# real and reactive and complex loads and DGs
	@defVar(m, pc[1:N] >= 0)
	@defVar(m, qc[1:N] >= 0)
	@defVar(m, pg[1:N] >= 0)
	@defVar(m, qg[1:N] >= 0)
  #@defVar(m, qg[1:N])

	# real and reactive power flows
	@defVar(m, P[1:N] )
	@defVar(m, Q[1:N] )

	#auxiliary variable
	@defVar(m, t >= 0)

	for i=1:N
		addConstraint(m, deltas[i] == deltas[i+N] );

		# load constraints
    addConstraint(m, pc[i] == gammas[i]*PCD[i]);
    addConstraint(m, qc[i] == gammas[i]*QCD[i]);

		# PV constraints
    addConstraint(m, pg[i] <= PGSP[i]);
    #addConstraint(m, qg[i] <= QGSP[i]);
    #addConstraint(m, qg[i] >= -SGmax[i]);
    addConstraint(m, qg[i] <= SGmax[i] + QGSP[i]);

		# power conservation
		PSUM = pc[i] - pg[i];
    QSUM = qc[i] - qg[i] + SGmax[i];

		for j in child[i]
			PSUM += P[j];
			QSUM += Q[j];
		end

		addConstraint(m, P[i] == PSUM);
		addConstraint(m, Q[i] == QSUM);

		#ohms law
		eqSum = -2*(r*P[i] + x*Q[i]);

		if i == 1
			eqSum += nus0;
		else
			eqSum += nus[par[i]];
		end

		addConstraint(m, nus[i] == eqSum);
		addConstraint(m, t >= nusL - nus[i]);

	end

	# deltas constraints
	for i = 1:N
		addConstraint(m, deltas[i] == ndeltas[i]);
	end

	# define cost vector c and defense y vector, and attack vector d
	c = zeros(Float64,ny);
	for i = 7N+1:8N
        c[i] = Cll*PCD[i-7N];
	end
	c[ny] = -Clovr;

	y = typeof(t)[];
	for i = 1:ny
		if i <= N
			push!(y, pg[i]);
		elseif i <= 2N
			push!(y, qg[i-N]);
		elseif i <= 3N
			push!(y, pc[i-2N]);
		elseif i <= 4N
			push!(y, qc[i-3N]);
		elseif i <= 5N
			push!(y, P[i-4N]);
		elseif i <= 6N
			push!(y, Q[i-5N]);
		elseif i <= 7N
			push!(y, nus[i-6N]);
		elseif i <= 8N
			push!(y, gammas[i-7N]);
		else
			push!(y, t);
		end
	end

	# define objective
	cc = zeros(Float64, ny);
	for i = 1:ny
		cc[i] = (c[i] - ndeltas[i] * rbar[i]) ;
	end
  @show ndeltas[1:2N]
  @show cc[1:2N]

	obj = sum (cc' * y);
  # if master problem, then minimize the cost for defender
  @setObjective(m, :Max, obj);

	m, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg, c, y, cc
end

function getDeltaEffects(comb, gammaMin, ratio, K)
  #Cll = K/(ratio^2+1)^0.5;
  Cll = K;
  Clovr = ratio*Cll;
  m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getNonLinearModel(Cll, Clovr, gammaMin);

  for i in 1:N
    if i in comb
      addConstraint(m, deltas[i] == 0);
    else
      addConstraint(m, deltas[i] == 1);
    end
  end

  @show comb
  targetNodes = copy(comb);

  solve(m);
  ngammas = copy(getValue(gammas[1:N]));
  nnus = copy(getValue(nus[1:N]));
  ndeltas = copy(getValue(deltas[1:N]));
  npg = copy(getValue(pg[1:N]));
  nqg = copy(getValue(qg[1:N]));

  lovr = getValue(l);
  voll = getValue(v);

  lovr, voll, targetNodes, ngammas, nnus, ndeltas, npg, nqg
end

function bruteForce(gammaMin, ratio, M, K)
  #Cll = K/(ratio^2+1)^0.5;
  #Clovr = ratio*Cll;
  Cll = K;
  Clovr = ratio*K;

   # initialize ngammas
  targetNodesStar = zeros(Int64, M);
  ndeltas = ones(Int64, M);
  ngammas = ones(Float64, N);
  nnus = ones(Float64, N);
  currentBestObjValue = -1000(Clovr+N*Cll);
  finalObjValue = currentBestObjValue;
  lovrValue = 0;
  vollValue = 0;
  minpgqgRatio = 0;
  maxpgqgRatio = 0;
  for comb in combinations(pvNodes, M)
    lovr, voll, targetNodes, ncgammas, ncnus, ncdeltas, npg, nqg =  getDeltaEffects(comb, gammaMin, ratio, K);

    finalObjValue = lovr + voll;
    if finalObjValue >= currentBestObjValue
      # update values
      currentBestObjValue = finalObjValue;
      lovrValue = lovr;
      vollValue = voll;
      ngammas = copy(ncgammas);
      nnus = copy(ncnus);
      ndeltas = copy(ncdeltas);
      targetNodesStar = copy(targetNodes);

      minpgqgRatio = 3;
      maxpgqgRatio = 0;
      for i in 1:N
        if !(i in comb)
          minpgqgRatio = min(minpgqgRatio,npg[i]/nqg[i]);
          maxpgqgRatio = max(maxpgqgRatio,npg[i]/nqg[i]);
        end
      end
    end
    #break;
  end
  lovrValue, vollValue, targetNodesStar, ngammas, nnus, ndeltas, minpgqgRatio, maxpgqgRatio
end
@show targetNodesStar
gammaMin = 0.5;
M = 2;
ratio = 6;
bfl, bfv, bftargetNodesStar, bfngammas, bfnnus, ndeltas  = bruteForce(gammaMin, ratio, M, K);
@show bfl bfv bftargetNodesStar bfngammas bfnnus

function greedyApproach(gammaMin, ratio, M, K)
  #Cll = K/(ratio^2+1)^0.5;
  #Clovr = ratio*Cll;
  Cll = K;
  Clovr = ratio*K;

  @show gammaMin ratio M Cll Clovr
  # initialize ngammas
  targetNodes = zeros(Int64, M);
  ngammas = ones(Float64, N);
  nnus = ones(Float64, N);
  ndeltas = ones(Float64, N);
  currentBestObjValue = -1000(Clovr  +N*Cll);
  finalObjValue = currentBestObjValue;
  ndeltas = zeros(Int64,N);
  maxIterations = 5;
  lovr = 0;
  voll = 0;
  tempTargetNodesList = Array{Int64,1}[];
  maxMaxIterations = -1;
  for iter = 1:maxIterations
    @show iter
    if iter > maxMaxIterations
      maxMaxIterations = iter
      #@show maxMaxIterations
      #@show gammaMin ratio M
    end
    ################################################################
    # Master problem
    ################################################################
    optimalAttackValues = zeros(Float64, N);
    optimalDeltasForPivotNode = Array{Int64,1}[];

    m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getLinearModel(Clovr, gammaMin);

    for i = 1:N
      addConstraint(m, gammas[i] == min(ngammas[i],1));
      addConstraint(m, deltas[i] == 1);
    end

    solve(m);

    for pivotNode = 1:N

      objValues = zeros(Float64, N);
      for targetNode = 1:N
        changeInNu = commonResistances[pivotNode,targetNode] * PGSP[targetNode] + commonResistances[pivotNode,targetNode] * (QGSP[targetNode] + SGmax[targetNode]);
        objValues[targetNode] = 2*changeInNu;
      end

      # find optimal target nodes for the pivot node
      optimalDeltas = Int64[];
      optimalAttackValue = getValue(nus[pivotNode]);
      for i = 1:M
        objValue = -Clovr;
        maxJ = -1;
        for j = 1:N
          if !(j in optimalDeltas) && objValues[j] >= objValue
            objValue = objValues[j];
            maxJ = j;
          end
        end

        optimalAttackValue -= objValues[maxJ];
        push!(optimalDeltas, maxJ);
      end

      optimalAttackValues[pivotNode] = optimalAttackValue;
      push!(optimalDeltasForPivotNode, optimalDeltas);
    end

    # find the best pivot node to target
    minValue = N*(Clovr+Cll);
    minI = -1;
    for i = 1:N
      if optimalAttackValues[i] <= minValue
        minValue = optimalAttackValues[i];
        minI = i;
      end
    end
  @show optimalAttackValues minValue
    optimalTargetNodes = optimalDeltasForPivotNode[minI];
    # initialize optimal attack vector for the best pivot node to target
    for i = 1:N
      if i in optimalTargetNodes
        ndeltas[i] = 0;
      else
        ndeltas[i] = 1;
      end
    end

    ################################################################
    # Slave problem
    ################################################################
    m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getNonLinearModel(Cll, Clovr, gammaMin);

    for i = 1:N
      addConstraint(m, deltas[i] == ndeltas[i]);
      #addConstraint(m, deltas[i] == 1);
      #addConstraint(m, gammas[i] == 1);
    end

    solve(m);

    finalObjValue = getObjectiveValue(m);
    @show finalObjValue currentBestObjValue
    @show getValue(gammas[1:N])

    if finalObjValue >= currentBestObjValue
      # update values
      currentBestObjValue = finalObjValue;
      nnus = copy(getValue(nus[1:N]));
      ngammas = copy(getValue(gammas[1:N]));
      ndeltaStar = copy(ndeltas);

      lovr = getValue(l);
      voll = getValue(v);

      j = 1;
      for i=1:N
        if ndeltas[i] == 0
          targetNodes[j] = i;
          j += 1;
        end
      end

      if targetNodes in tempTargetNodesList
        #@show targetNodes
        break;
      else
        push!(tempTargetNodesList, targetNodes);
      end
    end
  end
  lovr, voll, targetNodes, ngammas, nnus, ndeltas, Cll, Clovr
end

gammaMin = 0.5;
M = 12;
ratio = 6;
gal, gav, gatargetNodesStar, gangammas, gannus, gandeltas, Cll, Clovr  = greedyApproach(gammaMin, ratio, M, K);
@show gal gav gatargetNodesStar gangammas gannus
@show Cll, Clovr, ratio

@show min(2,3)
Clovr = ratio*Cll;
m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getLinearModel(Clovr, gammaMin);
ngammas = ones(Int64,N);
for i = 1:N
  addConstraint(m, gammas[i] == ngammas[i]);
  addConstraint(m, deltas[i] == 1);
end
solve(m);

@show PCD

ratio = 12;
comb = {1,2,3,4,5,6,7,8,9,10,11,12}
  lovr, voll, targetNodes, ngammas, nnus, ndeltas  =  getDeltaEffects(comb, gammaMin, ratio, K);
@show   lovr voll targetNodes ngammas nnus



#######################################################
########## Generate Topology ##############################
#######################################################
par = zeros(Int64, N);
for i = 1:N
    par[i] = i-1;
end
par[6] = 3;
par[8] = 2;
par[11] = 8;

pvNodes = [1:N]
child = getChildren(par);
nEdgesOnPath = countEdgesOnPath(par);
commonEdgesList = countCommonEdges(par);

# initialize the power requirements
#Base Values
gammaMin = 0.6
PDmax = 25000;
Pbase = N * PDmax;
Qbase = 0.3Pbase;
Sbase = Pbase * (1+0.3^2)^0.5;
Vbase = 4000;
Zbase = Vbase^2/Sbase;

#average length of lines
lineLength = 3.0; # miles

# line impedance parameters
rxRatios = {0.5,1,1.5,2};
r = Array(Float64,N);
x = Array(Float64,N);
for i = 1:N
  rr = rand();
  r[i] = 0.33*lineLength/Zbase;
  if rr <= 0.25
    rxRatio = rxRatios[1];
  elseif rr <= 0.5
    rxRatio = rxRatios[2];
  elseif rr <= 0.75
    rxRatio = rxRatios[3];
  else
    rxRatio = rxRatios[4];
  end
  x[i] = 0.38*lineLength/Zbase * rxRatio;
end

commonResistances, commonReactances = countCommonResistance(par);

# Constraints
nusL = 0.95^2;
nusU = 1.05^2;
nus0 = 1;

# load parameters
pcdRatios = Array(Float64, N);
qcdRatios = Array(Float64, N);
for i = 1:N
  pcdRatios[i] = rand();
  qcdRatios[i] = rand();
end
pcdRatioSum = sum(pcdRatios);
qcdRatioSum = sum(qcdRatios);

PCDmax = PDmax/Sbase;
PCDS = zeros(Float64, N);
QCDS = zeros(Float64, N);
for i = 1:N
	PCDS[i] = Pbase*pcdRatios[i]/pcdRatioSum/Sbase;
  QCDS[i] = Qbase*qcdRatios[i]/qcdRatioSum/Sbase;
end

# PV parameters
SGmaxS = zeros(Float64, N);
PGmaxS = zeros(Float64, N);
PGSPS = zeros(Float64, N);
QGSPS = zeros(Float64, N);

ravg = mean(r);
xavg = mean(x);
pvvvNodes = shuffle(pvNodes);
for i = 1:length(pvvvNodes)
  if i <= length(pvvvNodes)/3
    PGmaxS[i] = 0.3PDmax/Sbase;
  elseif i <= 2*length(pvvvNodes)/3
    PGmaxS[i] = 0.5PDmax/Sbase;
  else
    PGmaxS[i] = 0.7PDmax/Sbase;
  end

  SGmaxS[i] = 1.1PGmaxS[i];
  PGSPS[i] = ravg/(ravg^2+xavg^2)^0.5 * SGmaxS[i];
  QGSPS[i] = xavg/(ravg^2+xavg^2)^0.5 * SGmaxS[i];
end


############################# Load Values ############################

PCD = zeros(Float64, N);
QCD = zeros(Float64, N);
PGSP = zeros(Float64, N);
SGmax = zeros(Float64, N);
QGSP = zeros(Float64, N);
PGmax = zeros(Float64, N);

for i = 1:N
		PCD[i] = PCDS[i];
		QCD[i] = QCDS[i];
		SGmax[i] = SGmaxS[i];
		PGmax[i] = PGmaxS[i];
		PGSP[i] = PGSPS[i];
		QGSP[i] = QGSPS[i];
end

Cll = 7*Pbase/1000
K = Cll;

######################## Test ratios ########################
gammaMin = 0;
M = 0;
ratio = 6;
gal, gav, gatargetNodesStar, gangammas, gannus, gandeltas, Cll, Clovr  = greedyApproach(gammaMin, ratio, M, K);
@show gal gav gatargetNodesStar gangammas gannus
@show Cll, Clovr, ratio
######################## Test ratios ########################

######## ######## ######## ########
######## Brute Force #########
######## ######## ######## ########
gammaMins = {0,0.5,0.7,0.9,1};
ratios = {3,6,12};

minpgqgRatios = Float64[];
maxpgqgRatios = Float64[];
ls = Float64[];
vs = Float64[];
targetNodesList = Array{Int64,1}[];
for gammaMinIndex = 1:length(gammaMins)
  gammaMin = gammaMins[gammaMinIndex]
  for ratioIndex = 1:length(ratios)
    ratio = ratios[ratioIndex]
    for M = 0:N
      l, v, targetNodes, ngammas, nnus, ndeltas, minpgqgRatio, maxpgqgRatio = bruteForce(gammaMin, ratio, M, K);

      push!(ls, l);
      push!(vs, v);
      push!(targetNodesList, targetNodes);
      push!(minpgqgRatios,minpgqgRatio);
      push!(maxpgqgRatios,maxpgqgRatio);
      #writecsv(file, ["BF" gammaMin ratio M maxMaxIterations l v targetNodes'])
    end
  end
end
#close(file)
bfls = copy(ls);
bfvs = copy(vs);
bftargetNodesList = copy(targetNodesList);
@show bfvs
#############################################
######### Greedy Approach ################
#############################################

ls = Float64[];
vs = Float64[];
targetNodesList = Array{Int64,1}[];
for gammaMinIndex = 1:length(gammaMins)
    gammaMin = gammaMins[gammaMinIndex]
    for ratioIndex = 1:length(ratios)
        ratio = ratios[ratioIndex]
        for M = 0:N
            l, v, targetNodes, ngammas, nnus, ndeltas = greedyApproach(gammaMin, ratio, M, K);
            push!(ls, l);
            push!(vs, v);
            push!(targetNodesList, targetNodes);
        end
    end
end
gals = copy(ls);
gavs = copy(vs);
gatargetNodesList = copy(targetNodesList);
@show gavs

############################ write ###########################

pvl = length(pvNodes);
filename = "12hetero0.csv";
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BFL2" "BFV2" "GAL2" "GAV2" "BFL3" "BFV3" "GAL3" "GAV3" "PGQGMin1" "PGQGMax1" "PGQGMin2" "PGQGMax2"  "PGQGMin2" "PGQGMax2"]);
for M = 0:pvl
  i = 0*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] minpgqgRatios[i] maxpgqgRatios[i]  minpgqgRatios[i+pvl+1] maxpgqgRatios[i+pvl+1]  minpgqgRatios[i+2pvl+2] maxpgqgRatios[i+2pvl+2] ]);
end
close(file)

filename = "12hetero50.csv"
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BFL2" "BFV2" "GAL2" "GAV2" "BFL3" "BFV3" "GAL3" "GAV3" "PGQGMin1" "PGQGMax1" "PGQGMin2" "PGQGMax2"  "PGQGMin2" "PGQGMax2"]);
for M = 0:pvl
  i = 1*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] minpgqgRatios[i] maxpgqgRatios[i]  minpgqgRatios[i+pvl+1] maxpgqgRatios[i+pvl+1]  minpgqgRatios[i+2pvl+2] maxpgqgRatios[i+2pvl+2] ]);
end
close(file)

filename = "12hetero70.csv"
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BFL2" "BFV2" "GAL2" "GAV2" "BFL3" "BFV3" "GAL3" "GAV3" "PGQGMin1" "PGQGMax1" "PGQGMin2" "PGQGMax2"  "PGQGMin2" "PGQGMax2"]);
for M = 0:pvl
  i = 2*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] minpgqgRatios[i] maxpgqgRatios[i]  minpgqgRatios[i+pvl+1] maxpgqgRatios[i+pvl+1]  minpgqgRatios[i+2pvl+2] maxpgqgRatios[i+2pvl+2] ]);
end
close(file)

filename = "12hetero90.csv";
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BFL2" "BFV2" "GAL2" "GAV2" "BFL3" "BFV3" "GAL3" "GAV3" "PGQGMin1" "PGQGMax1" "PGQGMin2" "PGQGMax2"  "PGQGMin2" "PGQGMax2"]);
for M = 0:pvl
  i = 3*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] minpgqgRatios[i] maxpgqgRatios[i]  minpgqgRatios[i+pvl+1] maxpgqgRatios[i+pvl+1]  minpgqgRatios[i+2pvl+2] maxpgqgRatios[i+2pvl+2] ]);
end
close(file)

filename = "12hetero100.csv";
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BFL2" "BFV2" "GAL2" "GAV2" "BFL3" "BFV3" "GAL3" "GAV3" "PGQGMin1" "PGQGMax1" "PGQGMin2" "PGQGMax2"  "PGQGMin2" "PGQGMax2"]);
for M = 0:pvl
  i = 4*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] minpgqgRatios[i] maxpgqgRatios[i]  minpgqgRatios[i+pvl+1] maxpgqgRatios[i+pvl+1]  minpgqgRatios[i+2pvl+2] maxpgqgRatios[i+2pvl+2] ]);
end
close(file)
