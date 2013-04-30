function genotypeFactor = genotypeGivenParentsGenotypesFactor(numAlleles, genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo)
% This function computes a factor representing the CPD for the genotype of
% a child given the parents' genotypes.

% THE VARIABLE TO THE LEFT OF THE CONDITIONING BAR MUST BE THE FIRST
% VARIABLE IN THE .var FIELD FOR GRADING PURPOSES

% When writing this function, make sure to consider all possible genotypes 
% from both parents and all possible genotypes for the child.

% Input:
%   numAlleles: int that is the number of alleles
%   genotypeVarChild: Variable number corresponding to the variable for the
%   child's genotype (goes in the .var part of the factor)
%   genotypeVarParentOne: Variable number corresponding to the variable for
%   the first parent's genotype (goes in the .var part of the factor)
%   genotypeVarParentTwo: Variable number corresponding to the variable for
%   the second parent's genotype (goes in the .var part of the factor)
%
% Output:
%   genotypeFactor: Factor in which val is probability of the child having 
%   each genotype (note that this is the FULL CPD with no evidence 
%   observed)

% The number of genotypes is (number of alleles choose 2) + number of 
% alleles -- need to add number of alleles at the end to account for homozygotes

genotypeFactor = struct('var', [], 'card', [], 'val', []);

% Each allele has an ID.  Each genotype also has an ID.  We need allele and
% genotype IDs so that we know what genotype and alleles correspond to each
% probability in the .val part of the factor.  For example, the first entry
% in .val corresponds to the probability of having the genotype with
% genotype ID 1, which consists of having two copies of the allele with
% allele ID 1, given that both parents also have the genotype with genotype
% ID 1.  There is a mapping from a pair of allele IDs to genotype IDs and 
% from genotype IDs to a pair of allele IDs below; we compute this mapping 
% using generateAlleleGenotypeMappers(numAlleles). (A genotype consists of 
% 2 alleles.)

[allelesToGenotypes, genotypesToAlleles] = generateAlleleGenotypeMappers(numAlleles);

% One or both of these matrices might be useful.
%
%   1.  allelesToGenotypes: n x n matrix that maps pairs of allele IDs to 
%   genotype IDs, where n is the number of alleles -- if 
%   allelesToGenotypes(i, j) = k, then the genotype with ID k comprises of 
%   the alleles with IDs i and j
%
%   2.  genotypesToAlleles: m x 2 matrix of allele IDs, where m is the 
%   number of genotypes -- if genotypesToAlleles(k, :) = [i, j], then the 
%   genotype with ID k is comprised of the allele with ID i and the allele 
%   with ID j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INSERT YOUR CODE HERE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

% Fill in genotypeFactor.var.  This should be a 1-D row vector.
% Fill in genotypeFactor.card.  This should be a 1-D row vector.
genotypeNum = size(genotypesToAlleles, 1);

genotypeFactor.var = [genotypeVarChild, genotypeVarParentOne, genotypeVarParentTwo];
genotypeFactor.card = [genotypeNum, genotypeNum, genotypeNum];

genotypeFactor.val = zeros(1, prod(genotypeFactor.card));
% Replace the zeros in genotypeFactor.val with the correct values.

for j = 1:prod(genotypeFactor.card),
	assigment = IndexToAssignment(j, genotypeFactor.card);
	child_ass = genotypesToAlleles(assigment(1),:);
	child_par1 = genotypesToAlleles(assigment(2),:);
	child_par2 = genotypesToAlleles(assigment(3),:);
	possible_child = zeros(length(child_par1) * length(child_par1), 2);
	possible_child_genotype = zeros(length(child_par1) * length(child_par1), 1);
	c = 1;
	for p1 = 1:length(child_par1),
		for p2 = 1:length(child_par2),
			possible_child(c,:) = sort([child_par1(p1), child_par2(p2)]);
			possible_child_genotype(c) = allelesToGenotypes(possible_child(c,1), possible_child(c,2));
			c = c + 1; 
		end;
	end;
	%child_par1
	%child_par2
	pp = struct('c1',child_par1 , 'c2',child_par2, 'possible_child', possible_child,'possible_child_genotype', possible_child_genotype);
	%allelesToGenotypes%(possible_child)
	allelesToGenotypes(child_ass(1),child_ass(2));
	mismatchGeno = possible_child_genotype - allelesToGenotypes(child_ass(1),child_ass(2));
	mismatchNum = length(find(mismatchGeno));
	all_possibilty = length(possible_child_genotype);
	genotypeFactor.val(j) = (all_possibilty - mismatchNum)/all_possibilty;
	%genotypesToAlleles(mod(j, genotypeNum))
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%