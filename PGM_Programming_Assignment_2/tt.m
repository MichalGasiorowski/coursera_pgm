pedigree = struct('parents', [0,0;1,3;0,0;1,3;2,6;0,0;2,6;4,9;0,0]);
pedigree.names = {'Ira','James','Robin','Eva','Jason','Rene','Benjamin','Sandra','Aaron'};
phenotypeList = {'CysticFibrosis', 'NoCysticFibrosis'};
alleleFreqsThree = [0.1; 0.7; 0.2];
alleleListThree = {'F', 'f', 'n'};
alphaListThree = [0.8; 0.6; 0.1; 0.5; 0.05; 0.01];
factorListDecoupled = constructDecoupledGeneticNetwork(pedigree, alleleFreqsThree, alphaListThree);


ret = length(factorListDecoupled);
for j = 1:18,
	cc = factorListDecoupled(j).card;
	ret(j) = (cc(1) - 1) * prod(cc(2:end));
end;
ret
sum(ret)