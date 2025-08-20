morpho = readDiscreteCharacterData("Softshell Dataset.nex")
taxa <- morpho.names()
num_taxa <- taxa.size()
num_branches <- 2 * num_taxa - 2

moves    = VectorMoves()
monitors = VectorMonitors()

phylogeny ~ dnUniformTopologyBranchLength(taxa, branchLengthDistribution=dnExponential(10.0))

tree_length := phylogeny.treeLength()

moves.append( mvNNI(phylogeny, weight=num_branches) )
moves.append( mvSPR(phylogeny, weight=num_branches/5.0) )
moves.append( mvBranchLengthScale(phylogeny, weight=num_branches) )

Q_morpho <- fnJC(7)

alpha_morpho ~ dnUniform( 0.0, 1E8 )
alpha_morpho.setValue( 1.0 )
rates_morpho := fnDiscretizeGamma( alpha_morpho, alpha_morpho, 4 )

moves.append( mvScale(alpha_morpho,lambda=1, weight=2.0) )

phyMorpho ~ dnPhyloCTMC(tree=phylogeny, siteRates=rates_morpho, Q=Q_morpho, type="Standard")
phyMorpho.clamp(morpho)

mymodel = model(phylogeny)

monitors.append( mnModel(filename="output/mk.log", printgen=10) )

monitors.append( mnFile(filename="output/mk.trees", printgen=10, phylogeny) )

monitors.append( mnScreen(printgen=100) )

mymcmc = mcmc(mymodel, monitors, moves, nruns=2, combine="mixed")

mymcmc.run(generations=250000, tuningInterval=200)

trace = readTreeTrace("output/mk.trees", treetype="non-clock")
trace.setBurnin(0.25)

mapTree(trace, file="output/mk.map.tre")
consensusTree(trace, file="output/mk.majrule.tre")

q()