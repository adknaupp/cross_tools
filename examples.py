from library import *

print("FRUITFLY TESTCROSS EXAMPLE 1:")
print("Cross between a white eyed, vestigial winged parent and a red eyed, wild type winged parent.")
print("The first parent is homozygous recessive for both traits, the other is homozygous dominant for both.")

# White Eyed, Vestigial winged
parent2 = Two_factor_unlinked(Assortment('w','w', dominant_phenotype='Red Eyes', recessive_phenotype='White Eyes'), Assortment('b','b',dominant_phenotype='WT wings', recessive_phenotype='Vestigial wings'))
# Red Eyed, WT winged
parent1 = Two_factor_unlinked(Assortment('W','W', dominant_phenotype='Red Eyes', recessive_phenotype='White Eyes'), Assortment('B','B',dominant_phenotype='WT wings', recessive_phenotype='Vestigial wings'))

f1_generation = parent1.cross(parent2)

print("")
print("Progeny display only the dominant phenotype and only one genotype is present:")
print(f1_generation.phenotypic_ratio)
for i, genotype in enumerate(f1_generation.genotypes):
  print("Genotype ", str(i+1), ": ", genotype.genotype)

print("")
print("Crossing two individuals of the f1 generation yields the following results:")
f2_generation = f1_generation.genotypes[0].cross(f1_generation.genotypes[0])
print(f2_generation.phenotypic_ratio)
print(str(len(f2_generation.genotypes)), " genotypes are found in the f2 generation:")
for i, genotype in enumerate(f2_generation.genotypes):
  print("Genotype ", str(i+1), ": ", genotype.genotype)

print("")
print("FRUITFLY TESTCROSS EXAMPLE 2:")
print("Cross between red eyed, wild type winged parent and a white eyed, vestigial winged parent.")
print("This time, the second parent is heterozygous for both traits. The first parent remains homozygous recessive for both.")

# Red Eyed, WT winged but heterozygous
parent3 = Two_factor_unlinked(Assortment('W','w', dominant_phenotype='Red Eyes', recessive_phenotype='White Eyes'), Assortment('B','b',dominant_phenotype='WT wings', recessive_phenotype='Vestigial wings'))

f1_generation = parent1.cross(parent3)

# Progeny display only dominant phenotype
print("")
print("Progeny display only the dominant phenotype, but this time, four genotypes exist:")
print(f1_generation.phenotypic_ratio)
for i, genotype in enumerate(f1_generation.genotypes):
  print("Genotype ", str(i+1), ": ", genotype.genotype)
