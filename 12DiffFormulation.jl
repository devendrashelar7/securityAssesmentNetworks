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
  @defVar(m, t[1:N] >= 0)

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
		eqSum = -2*(r*P[i] + x*Q[i]) ;

		if i == 1
			eqSum += nus0;
		else
			eqSum += nus[par[i]];
		end

		addConstraint(m, nus[i] == eqSum);
		addConstraint(m, t[i] >= nusL - nus[i]);

	end

	# define objective
	addConstraint(m, l == Clovr*sum(t));
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
		PSUM = pc[i] - pg[i] + r * lij[i];
		QSUM = qc[i] - qg[i] + x * lij[i];

		for j in child[i]
			PSUM += P[j];
			QSUM += Q[j];
		end

		addConstraint(m, P[i] == PSUM);
		addConstraint(m, Q[i] == QSUM);

		#ohms law
		eqSum = -2*(r*P[i] + x*Q[i]) + (r^2 + x^2) * lij[i];

		if i == 1
			eqSum += nus0;
		else
			eqSum += nus[par[i]];
		end

		addConstraint(m, nus[i] == eqSum);
		addConstraint(m, t[i] >= nusL - nus[i]);

		#current constraint
		addConstraint(m, lij[i] * nus[i] >= P[i]^2 + Q[i]^2);

	end

	# define objective
	voll = 0;
  obj = 0;
	for i=1:N
    voll += Cll * (1-gammas[i]) * PCD[i];
		obj += Cll * (1-gammas[i]) * PCD[i] + r*lij[i];
	end
  addConstraint(m, v == voll);
  addConstraint(m, l == Clovr*sum(t));

	obj += l;

  setObjective(m, :Min, obj);
	m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg
end

ny = 9N;
function getRBar(Cll, Clovr)
  ny = 9N;
	diameter = maximum(nEdgesOnPath);
	R = zeros(Float64, ny);
  for i = 1:ny
    if i <= N
      # pg
      R[i] = 2 * nEdgesOnPath[i] * r * Clovr * N  +  Cll;
    elseif i <= 2N
      # qg
      R[i] = 2 * nEdgesOnPath[i-N] * x * Clovr * N +  Cll; # qg
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
  ny = 9N;
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
	@defVar(m, t[1:N] >= 0)
  @defVar(m, l >= 0)
  @defVar(m, v >= 0)

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
		addConstraint(m, t[i] >= nusL - nus[i]);

	end

  addConstraint(m, l == Clovr*sum(t));

	# deltas constraints
	for i = 1:N
		addConstraint(m, deltas[i] == ndeltas[i]);
	end

	# define cost vector c and defense y vector, and attack vector d
	c = zeros(Float64,ny);
	for i = 7N+1:8N
        c[i] = Cll*PCD[i-7N];
	end
  for i = 8N+1:9N
	  c[i] = -Clovr;
  end

	y = typeof(l)[];
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
			push!(y, t[i-8N]);
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

	m, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg, c, y, cc, l, v
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

  #@show comb
  targetNodes = copy(comb);

  solve(m);
  ngammas = copy(getValue(gammas[1:N]));
  nnus = copy(getValue(nus[1:N]));
  ndeltas = copy(getValue(deltas[1:N]));
  nP = copy(getValue(P[1:N]));
  nQ = copy(getValue(Q[1:N]));
  npc = copy(getValue(pc[1:N]))- copy(getValue(pg[1:N]));
  nqc = copy(getValue(qc[1:N])) - copy(getValue(qg[1:N]));
  nqg = copy(getValue(qg[1:N]));

  #@show getValue(t[1:N]) sum(getValue(t[1:N])) Clovr*sum(getValue(t[1:N]))
  lovr = getValue(l);
  voll = getValue(v);

  lovr, voll, targetNodes, ngammas, nnus, ndeltas, nP, nQ, npc, nqc, nqg
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
  for comb in combinations(pvNodes, M)
    lovr, voll, targetNodes, ncgammas, ncnus, ncdeltas, nP, nQ, npc, nqc, nqg =  getDeltaEffects(comb, gammaMin, ratio, K);

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
    end
    #break;
  end
  lovrValue, vollValue, targetNodesStar, ngammas, nnus, ndeltas
end
@show targetNodesStar
gammaMin = 0.5;
M = 10;
ratio = 3;
bfl, bfv, bftargetNodesStar, bfngammas, bfnnus, ndeltas  = bruteForce(gammaMin, ratio, M, K);
@show bfl bfv bftargetNodesStar bfngammas bfnnus

comb = {3,4,5,6,7,8,9,10,11,12}
lovr, voll, targetNodes, ngammas, nnus, ndeltas, nP, nQ, npc, nqc = getDeltaEffects(comb, gammaMin, ratio, K)
@show lovr voll targetNodes ngammas nnus ndeltas

a = {2,3,1};
ind = sortperm(a, rev = true);
a

function greedyApproach(gammaMin, ratio, M, K)
  #Cll = K/(ratio^2+1)^0.5;
  #Clovr = ratio*Cll;
  Cll = K;
  Clovr = ratio*K;

  #@show gammaMin ratio M Cll Clovr
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
    #optimalAttackValues = zeros(Float64, N);
    #optimalDeltasForPivotNode = Array{Int64,1}[];

    #m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getLinearModel(Clovr, gammaMin);

    #for i = 1:N
      #addConstraint(m, gammas[i] == ngammas[i]);
      #addConstraint(m, deltas[i] == 1);
    #end

    #solve(m);
    objValues = zeros(Float64, N);
    for targetNode = 1:N
      changeInNu = 2*(r * PGSP[targetNode] + x * (QGSP[targetNode] + SGmax[targetNode]));
      for pivotNode = 1:N
        objValues[targetNode] += commonEdgesList[pivotNode][targetNode] * changeInNu;
      end
    end

    # find optimal target nodes
    ind = sortperm(objValues, rev = true);
    @show ind

    optimalTargetNodes = zeros(Int64,M);
    # initialize optimal attack vector for the best pivot node to target
    for i = 1:M
      optimalTargetNodes[i] = ind[i];
    end
    #@show optimalTargetNodes
    for i = 1:N
      if i in optimalTargetNodes
        ndeltas[i] = 0;
      else
        ndeltas[i] = 1;
      end
    end
    #@show 1 in optimalTargetNodes
    #@show ndeltas
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
    #@show finalObjValue currentBestObjValue
    #@show getValue(gammas[1:N])

    if finalObjValue >= currentBestObjValue
      # update values
      currentBestObjValue = finalObjValue;
      nnus = copy(getValue(nus[1:N]));
      ngammas = copy(getValue(gammas[1:N]));
      ndeltaStar = copy(ndeltas);

      lovr = Clovr*sum(getValue(t[1:N]));
      voll = getValue(v);

      targetNodes = copy(optimalTargetNodes);

      if targetNodes in tempTargetNodesList
        #@show targetNodes
        break;
      else
        push!(tempTargetNodesList, targetNodes);
      end
    end
  end

  #lovr1, voll1, targetNodes1, ngammas1, nnus1, ndeltas1 = getDeltaEffects(targetNodes, gammaMin, ratio, K)

  lovr, voll, targetNodes, ngammas, nnus, ndeltas, Cll, Clovr
end

gammaMin = 0.5;
M = 12;
ratio = 0.35;
gal, gav, gatargetNodesStar, gangammas, gannus, gandeltas, Cll, Clovr  = greedyApproach(gammaMin, ratio, M, K);
@show gal gav gatargetNodesStar gangammas gannus
@show Cll, Clovr, ratio
@show gandeltas
comb = {3,4,5,6,7,8,9,10,11,12}
lovr, voll, targetNodes, ngammas, nnus, ndeltas = getDeltaEffects(comb, gammaMin, ratio, K)
@show lovr voll targetNodes ngammas nnus ndeltas

@show (12 in [11 12])


function bendersCut(gammaMin, ratio, M, K, maxIterations)
  #Cll = K/(ratio^2+1)^0.5;
  #Clovr = ratio*Cll;
  Cll = K;
  Clovr = ratio*K;
  ny = 9N;
  rbar = getRBar(Cll, Clovr);
  @show gammaMin
  # initialize ngammas
  targetNodes = zeros(Int64, M);
  targetNodesStar = zeros(Int64, M);
  ndeltas = zeros(Int64, ny);
  ngammas = ones(Float64,N);
  nnus = ones(Float64,N);
  ndeltaStar = zeros(Int64, ny);
  alovr = 0;
  avoll = 0;
  nP = 0;
  nQ = 0;
  npc  =0;
  nqc = 0;
  nqg = 0;

  zL = -(Clovr + N * Cll);
  zU = -zL;
  lhs = 0;
  Y = Array{Float64, 1}[];

  epsilon = 1;

  # master problem ADMP2
  mp = Model();
  @defVar(mp, z );
  @defVar(mp, mdeltas[1:ny], Bin);
  d = typeof(z)[];
  for i = 1:N
    push!(d, mdeltas[i]);
  end

  addConstraint(mp, sum(d) == M);
  # delta constraints
  for i = 1:N
    addConstraint(mp, mdeltas[i] == mdeltas[i+N]);
  end
  for i = 2N + 1:ny
    addConstraint(mp, mdeltas[i] == 0);
  end

  @setObjective(mp, :Min, z);

  for iter = 1:maxIterations
    #@show iter ndeltas[1:N]

    # slave problem
    m, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg, c, y, cc,l, v = getADLP2Model(Cll, Clovr,ndeltas, rbar);
    #@show obj

    solve(m);

    y = Float64[];
    for i = 1:ny
      if i <= N
        push!(y, getValue(pg[i]));
      elseif i <= 2N
        push!(y, getValue(qg[i-N]));
      elseif i <= 3N
        push!(y, getValue(pc[i-2N]));
      elseif i <= 4N
        push!(y, getValue(qc[i-3N]));
      elseif i <= 5N
        push!(y, getValue(P[i-4N]));
      elseif i <= 6N
        push!(y, getValue(Q[i-5N]));
      elseif i <= 7N
        push!(y, getValue(nus[i-6N]));
      elseif i <= 8N
        push!(y, getValue(gammas[i-7N]));
      else
        push!(y, getValue(t[i-8N]));
      end
    end
    #@show y[1:N]'
    #@show y[N+1:2N]'
    zhat = getObjectiveValue(m);

    #@show zhat
    #@show zU
    #@show getValue(t[1:N])
    #@show sum(getValue(t[1:N])) sum(getValue(t[1:N]))*Clovr
    #@show getValue(pg[1:N]) getValue(qg[1:N])


    push!(Y, y);
    if zhat < zU
      zU = zhat;
      j = 1;
      ndeltaStar = copy(ndeltas);
      targetNodesStar = copy(targetNodes);

      nnus = copy(getValue(nus[1:N]));
      ngammas = copy(getValue(gammas[1:N]));
      nP = copy(getValue(P[1:N]))
      nQ = copy(getValue(Q[1:N]))
      #@show ngammas nP nQ
      npc = copy(getValue(pc[1:N])) - copy(getValue(pg[1:N])) ;
      nqc = copy(getValue(qc[1:N])) - copy(getValue(qg[1:N])) + copy(SGmax[1:N]);
      nqg = copy(getValue(qg[1:N])) - copy(SGmax[1:N]);

      alovr = getValue(l);
      avoll = 0;
      for i = 1:N
        avoll += Cll*(1 - ngammas[i])*PCD[i];
      end
      #@show alovr avoll getValue(nus[1:N])

    end
    #@show targetNodesStar
    if zU - zL <= epsilon
      break;
    end

    # master problem
    lhs = z;
    for i = 1:ny
      lhs += y[i] * rbar[i] * mdeltas[i];
    end
    addConstraint(mp, lhs >= sum(c' * y));
    #@show lhs
    #@show sum(c' * y)
    solve(mp);

    z2 = getObjectiveValue(mp);
    zL = z2;
    j = 1;
    for i=1:ny
      if getValue(mdeltas[i]) < 0.5
        ndeltas[i] = 0
      else
        ndeltas[i] = 1
      end
      if i <= N && ndeltas[i] == 1
        targetNodes[j] = i;
        j += 1;
      end
    end
    #@show targetNodes

    #@show iter ngammas avoll alovr targetNodesStar
    #@show targetNodes
    #@show ndeltas[1:N]
    if zU - zL > epsilon
      continue;
    end
    #@show obj

  end
  comb = targetNodesStar;
  lovr, voll, targetNodes, ngammas, nnus, ndeltas, nPP, nQQ, npcc, nqcc, nqgg = getDeltaEffects(comb, gammaMin, ratio, K);

  lovr, voll, targetNodes, ngammas, nnus, ndeltas, Cll, Clovr, alovr, avoll, nP, nQ, nPP, nQQ, npc, nqc, npcc, nqcc, nqg, nqgg
end

maxIterations = 100;
M = 4;
ratio = 1;
gammaMin = 0.5;
bcl, bcv, bctargetNodesStar, bcngammas, bcnnus, ndeltas, Cll, Clovr, bcal, bcav, nP, nQ, nPP, nQQ, npc, nqc, npcc, nqcc, nqg, nqgg = bendersCut(gammaMin, ratio, M, K, maxIterations);
@show bcl bcv bcal bcav bctargetNodesStar bcngammas
@show nP-nPP nQ-nQQ
@show npc-npcc nqc nqcc nqg nqgg
comb = {1,2,3,4,5,6,7,8,9,10,11,12}
lovr, voll, targetNodes, ngammas, nnus, ndeltas, nPP, nQQ, npcc, nqcc  =  getDeltaEffects(comb, gammaMin, ratio, K);
@show lovr voll ngammas nnus npcc nqcc

@show ratio,K, Cll, Clovr
@show bcavs[55]

@show QGSP
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
Sbase = Pbase * (1+0.3^2)^0.5;
Vbase = 4000;
Zbase = Vbase^2/Sbase;

#average length of lines
lineLength = 4.0; # miles

# line impedance parameters
r = 0.33*lineLength/Zbase;
x = 0.38*lineLength/Zbase;


# Constraints
nusL = 0.95^2;
nusU = 1.05^2;
nus0 = 1;

# load parameters
PCDmax = PDmax/Sbase;
PCDS = zeros(Float64, N);
QCDS = zeros(Float64, N);
for i = 1:N
	PCDS[i] = PCDmax;
	QCDS[i] = 0.3PCDS[i];
end

# PV parameters
SGmaxS = zeros(Float64, N);
PGmaxS = zeros(Float64, N);
PGSPS = zeros(Float64, N);
QGSPS = zeros(Float64, N);

for i = 1:N
	if i in pvNodes
	    PGmaxS[i] = 0.5PDmax/Sbase;

	    SGmaxS[i] = 1.1PGmaxS[i];
	    PGSPS[i] = r/(r^2+x^2)^0.5 * SGmaxS[i];
	    QGSPS[i] = r/(r^2+x^2)^0.5 * SGmaxS[i];
	end
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

######## ######## ######## ########
######## Brute Force #########
######## ######## ######## ########

Clovr = 5;
Cll=5;

rbar = getRBar();
objValue = -100Clovr;
ls = Float64[];
vs = Float64[];
targetNodesList = Array{Int64,1}[];
targetNodes = 0;
nnus = zeros(Float64, N);
ngammas = ones(Float64, N);
l = 0;
v = 0;

ALL = [1:N]
Cll = 7*Pbase/1000
Clovrp = 5Cll;
#K = (Cllp^2 + Clovrp^2)^0.5
K = Cll;


gammaMin = 0.5;
gammaMins = {0,0.5,0.7,0.9,1};
ratios = {0.35,0.45,1.4}
filename = "12diff50.csv"
#file = open(filename, "w")
maxMaxIterations = -1

ls = Float64[];
vs = Float64[];
targetNodesList = Array{Int64,1}[];
for gammaMinIndex = 1:length(gammaMins)
    gammaMin = gammaMins[gammaMinIndex]
    for ratioIndex = 1:length(ratios)
        ratio = ratios[ratioIndex]
        for M = 0:N
            l, v, targetNodes, ngammas, nnus, ndeltas = bruteForce(gammaMin, ratio, M, K);

            push!(ls, l);
            push!(vs, v);
            push!(targetNodesList, targetNodes);
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

gammaMins = {0,0.5,0.7,0.9,1};
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

@show gals[1:N+1] gals[N+2:2N+2] gals[2N+3:3N+3] gals[3N+4:4N+4] gals[4N+5:5N+5] gals[5N+6:6N+6] gals[6N+7:7N+7] gals[7N+8:8N+8] gals[8N+9:9N+9]
############################################################
############ Benders Decomposition ##################
############################################################
ls = Float64[];
vs = Float64[];
bcals = Float64[];
bcavs = Float64[];
gammaMins = {0,0.5,0.7,0.9,1};
targetNodesList = Array{Int64,1}[];
maxIterations = 100;
for gammaMinIndex = 1:length(gammaMins)
  gammaMin = gammaMins[gammaMinIndex]
  for ratioIndex = 1:length(ratios)
    ratio = ratios[ratioIndex]
    for M = 0:N
      bcl, bcv, bctargetNodesStar, bcngammas, bcnnus, ndeltas, Cll, Clovr, bcal, bcav = bendersCut(gammaMin, ratio, M, K, maxIterations);

      push!(ls, bcl);
      push!(vs, bcv);
      push!(bcals, bcal);
      push!(bcavs, bcav);
      push!(targetNodesList, bctargetNodesStar);
    end
  end
end
bcls = copy(ls);
bcvs = copy(vs);
bctargetNodesList = copy(targetNodesList);
@show bcvs[1:39] bcvs[40:78] bcls[1:39] bcls[40:78]
@show bcavs[1:39] bcavs[40:78]
@show gammaMins*10
length(bctargetNodesList)
@show bcls[27:39]


filename = "12diff0.csv";
file = open(filename,"w");
pvl = length(pvNodes);
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BCL1" "BCV1" "BFL2" "BFV2" "GAL2" "GAV2" "BCL2" "BCV2" "BFL3" "BFV3" "GAL3" "GAV3" "BCL3" "BCV3" "BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3" ]);
for M = 0:pvl
  i = 0*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bcls[i] bcvs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bcls[i+pvl+1] bcvs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] bcls[i+2pvl+2] bcvs[i+2pvl+2]  bcals[i] bcavs[i]  bcals[i+pvl+1] bcavs[i+pvl+1]  bcals[i+2pvl+2] bcavs[i+2pvl+2] ]);
end
close(file)


filename = "12diff50.csv"
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BCL1" "BCV1" "BFL2" "BFV2" "GAL2" "GAV2" "BCL2" "BCV2" "BFL3" "BFV3" "GAL3" "GAV3" "BCL3" "BCV3" "BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3" ]);
for M = 0:pvl
  i = 1*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bcls[i] bcvs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bcls[i+pvl+1] bcvs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] bcls[i+2pvl+2] bcvs[i+2pvl+2]  bcals[i] bcavs[i]  bcals[i+pvl+1] bcavs[i+pvl+1]  bcals[i+2pvl+2] bcavs[i+2pvl+2] ]);
end
close(file)

filename = "12diff70.csv"
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BCL1" "BCV1" "BFL2" "BFV2" "GAL2" "GAV2" "BCL2" "BCV2" "BFL3" "BFV3" "GAL3" "GAV3" "BCL3" "BCV3" "BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3" ]);
for M = 0:pvl
  i = 2*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bcls[i] bcvs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bcls[i+pvl+1] bcvs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] bcls[i+2pvl+2] bcvs[i+2pvl+2]  bcals[i] bcavs[i]  bcals[i+pvl+1] bcavs[i+pvl+1]  bcals[i+2pvl+2] bcavs[i+2pvl+2] ]);
end
close(file)

filename = "12diff90.csv";
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BCL1" "BCV1" "BFL2" "BFV2" "GAL2" "GAV2" "BCL2" "BCV2" "BFL3" "BFV3" "GAL3" "GAV3" "BCL3" "BCV3" "BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3" ]);
for M = 0:pvl
  i = 3*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bcls[i] bcvs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bcls[i+pvl+1] bcvs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] bcls[i+2pvl+2] bcvs[i+2pvl+2]  bcals[i] bcavs[i]  bcals[i+pvl+1] bcavs[i+pvl+1]  bcals[i+2pvl+2] bcavs[i+2pvl+2] ]);
end
close(file)

filename = "12diff100.csv";
file = open(filename,"w")
writecsv(file, ["M" "BFL1" "BFV1" "GAL1" "GAV1" "BCL1" "BCV1" "BFL2" "BFV2" "GAL2" "GAV2" "BCL2" "BCV2" "BFL3" "BFV3" "GAL3" "GAV3" "BCL3" "BCV3" "BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3" ]);
for M = 0:pvl
  i = 4*(3pvl+3) + M+1
  writecsv(file, [M bfls[i] bfvs[i] gals[i] gavs[i] bcls[i] bcvs[i] bfls[i+pvl+1] bfvs[i+pvl+1] gals[i+pvl+1] gavs[i+pvl+1] bcls[i+pvl+1] bcvs[i+pvl+1] bfls[i+2pvl+2] bfvs[i+2pvl+2] gals[i+2pvl+2] gavs[i+2pvl+2] bcls[i+2pvl+2] bcvs[i+2pvl+2]  bcals[i] bcavs[i]  bcals[i+pvl+1] bcavs[i+pvl+1]  bcals[i+2pvl+2] bcavs[i+2pvl+2] ]);
end
close(file)


gammaMin = 0.7;
M = 6;
ratio = 6;
Clovr = ratio*Cll;
gal, gav, gatargetNodesStar, gangammas, gannus, ndeltas  = greedyApproach(gammaMin, ratio, M, K);

maxIterations = 100;

m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg = getNonLinearModel(Cll, Clovr, gammaMin);
for i = 1:N
  addConstraint(m, gammas[i] == 1);
  addConstraint(m, deltas[i] == 1);
end
solve(m);
nonus = copy(getValue(nus[1:N]));
m, l, v, deltas, gammas, t, obj, nus, P, Q, pc, qc, pg, qg= getNonLinearModel(Cll, Clovr, gammaMin);
for i = 1:N
  addConstraint(m, gammas[i] == 1);
  addConstraint(m, deltas[i] == ndeltas[i]);
end
solve(m);
prelcnus = copy(getValue(nus[1:N]));

bfl, bfv, bftargetNodesStar, bfngammas, bfnnus, ndeltas  = bruteForce(gammaMin, ratio, M, K);
gal, gav, gatargetNodesStar, gangammas, gannus, gandeltas, Cll, Clovr  = greedyApproach(gammaMin, ratio, M, K);
bcl, bcv, bctargetNodesStar, bcngammas, bcnnus, ndeltas,  Cll, Clovr, bcal, bcav  = bendersCut(gammaMin, ratio, M, K, maxIterations);
@show bfnnus gannus bcnnus bftargetNodesStar gatargetNodesStar bctargetNodesStar bfngammas gngammas bcngammas bfl gal bcl bfv gav bcv

path = {1,2,8,11,12};
levels = Array{Int64,1}[];
push!(levels,{1});
push!(levels,{2});
push!(levels,{3,8});
push!(levels,{4,6,9,11});
push!(levels,{5,7,10,12});

nAttacks = 100;
filename = "12nus.csv";
file = open(filename,"w");
VV = Array(Float64, length(levels)+1, 6+nAttacks);
for j = 2:7
  VV[1,j]=1
end
VV[1,1] = 0;
for i = 2:length(levels)+1
  nonv = 2;preclv = 2; bfv = 2; gav = 2; bcv=2;
  for j in levels[i-1]
    nonv = min(nonv,nonus[j]);
    preclv = min(preclv,prelcnus[j]);
    bfv = min(bfv,bfnnus[j]);
    gav = min(gav,gannus[j]);
    bcv = min(bcv,bcnnus[j]);
  end
  VV[i,1] = i-1; VV[i,2] = nonv; VV[i,3] = preclv; VV[i,4] = bfv; VV[i,5] = gav; VV[i,6] = bcv;
end
VV
newLabels = ["N" "Before Attack" "Before LC" "BF" "GA" "BC"];
for rr = 7:nAttacks+6
  pvvNodes = shuffle(pvNodes);
  comb = Array(Int64,M);
  nv = ones(Float64,length(levels)+1);
  newLabels = hcat(newLabels,string("rr",rr));
  VV[1,rr] = 1;
  for mj = 1:M
    comb[mj] = pvvNodes[mj];
    lovr, voll, targetNodes, ncgammas, ncnus, ncdeltas =  getDeltaEffects(comb, gammaMin, ratio, K);
    for i = 2:length(levels)+1
      for j in levels[i-1]
        nv[i-1] = min(nv[i-1],ncnus[j]);
      end
      VV[i,rr] = nv[i-1];
    end
  end
end
VVV = vcat(newLabels,VV);
writedlm("12nus.csv",VVV,",");

close(file);
min(2,3)
@show bfnnus gannus bcnnus
@show levels[4]
@show bctargetNodesList


gammaMin = 0.5;
M = 4;
ratio = 5;
Cll = K/(ratio^2+1)^0.5;
Clovr = ratio*Cll;
rbar = getRBar();
bfl, bfv, bftargetNodesStar, bfgammas, bfnnus, ndeltas  = bruteForce(gammaMin, ratio, M, K);
gl, gv, gtargetNodesStar, gngammas, gannus, ndeltas  = greedyApproach(gammaMin, ratio, M, K);
@show bfl gl bfv gv bftargetNodesStar gtargetNodesStar bfgammas gngammas
maxIterations = 200;
bcl, bcv, bctargetNodesStar, bcngammas, bcnnus, ndeltas  = bendersCut(gammaMin, ratio, M, K, maxIterations);
@show bcl bcv bctargetNodesStar bcngammas
@show bfnnus gannus bcnnus bftargetNodesStar gtargetNodesStar bctargetNodesStar bfngammas gngammas bcngammas bfl gl bcl bfv gv bcv

@show bfls bfvs bftargetNodesList


@show bfls-bcls

@show bctargetNodesList[2]

ls = Float64[];
vs = Float64[];
ratios = {100000}
gals = Array(Float64,N+1,length(ratios));
gavs = Array(Float64,N+1,length(ratios));
targetNodesList = Array{Int64,1}[];
gammaMins = {0.5};

for gammaMinIndex = 1:length(gammaMins)
    gammaMin = gammaMins[gammaMinIndex]
    for ratioIndex = 1:length(ratios)
        ratio = ratios[ratioIndex]
        for M = 0:length(pvNodes)
            l, v, targetNodes, ngammas, nnus, ndeltas = greedyApproach(gammaMin, ratio, M, K);
            gals[M+1,ratioIndex] = l;
            gavs[M+1,ratioIndex] = v;
        end
    end
end
@show Cll Clovr gals[:,1];

@show gavs[:,1] gavs[:,2] gavs[:,3] gavs[:,4] gavs[:,5] gavs[:,6]
@show gals[:,1] gals[:,2] gals[:,3] gals[:,4] gals[:,5] gals[:,6]


gals = copy(ls);
gavs = copy(vs);
filename = "ratios.csv"
file = open(filename, "w");
for i = 1:N
end


maximum(gavs[40:78])


filename = "12lateralNNew50.csv"
file = open(filename,"w");
i = 1;
rr = @show string("r", i);
newLabels = ["BCAL1" "BCAV1" "BCAL2" "BCAV2" "BCAL3" "BCAV3"];
hcat(newLabels,rr)
data = readcsv("12lateralNew100.csv");
nd = copy(data);
gi = 4;
filename = "12lateralNNew100.csv";
bcaValues = Array(Any,N+1,0);
bcValues = Array(Any,N+1,0);
for ri = 0:2
  id = gi*(3*(N+1)) + ri*(N+1)+1;
  bcaValues = hcat(bcaValues,bcals[id:id+N], bcavs[id:id+N]);
  bcValues = hcat(bcValues,bcls[id:id+N], bcvs[id:id+N]);
end
bcaValues = vcat(newLabels, bcaValues);
nd = hcat(nd, bcaValues);
writedlm(filename,nd,",");

@show nd
@show bcaValues
@show bcls-gals

