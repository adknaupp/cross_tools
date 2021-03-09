class Assortment:
  def __init__(self, a1, a2, **kwargs):
    self.allele1=a1
    self.allele2=a2
    self.gene_name='(default)'
    self.assortment=a1+a2
    #self.dominant_phenotype='(default)'
    #self.recessive_phenotype='(default)'
    # if provided by user, assign dominance
    if 'allele1_is_dominant' in kwargs: self.allele1_is_dominant=kwargs['allele1_is_dominant']
    else:
      if self.allele1.islower(): self.allele1_is_dominant=False # otherwise, if allele is all lowercase, assign dominance as recessive
      else: self.allele1_is_dominant=True
    if 'allele2_is_dominant' in kwargs: self.allele2_is_dominant=kwargs['allele2_is_dominant']
    else:
      if self.allele2.islower(): self.allele2_is_dominant=False # otherwise, if allele is all lowercase, assign dominance as recessive
      else: self.allele2_is_dominant=True
    # initialize zygosity
    if self.allele1_is_dominant and self.allele2_is_dominant: self.zygosity='Homozygous Dominant'
    elif not self.allele1_is_dominant and not self.allele2_is_dominant: self.zygosity='Homozygous Recessive'
    else: self.zygosity='Heterozygous'
    # if provided by user, assign phenotypes
    if 'dominant_phenotype' in kwargs: 
      self.dominant_phenotype=kwargs['dominant_phenotype']
    if 'recessive_phenotype' in kwargs: 
      self.recessive_phenotype=kwargs['recessive_phenotype']
  def _print(self):
    output_string=self.gene_name+' '+self.allele1
    if self.allele1_is_dominant:
      output_string+=' (dom), '
    else: output_string+=' (rec), '
    output_string+=self.allele2
    if self.allele2_is_dominant:
      output_string+=' (dom)'
    else: output_string+=' (rec)'
    return output_string
  def __eq__(self, other): # declares equivalence between two genes with the same alleles in any order.
    return ( self.allele1+self.allele2 == other.allele1+other.allele2 or self.allele1+self.allele2 == other.allele2+other.allele1 )

class Three_factor_linkage:
  def __init__(self, assortment1, assortment2, assortment3):
    self.gene1 = assortment1
    self.gene2 = assortment2
    self.gene3 = assortment3
    self.non_recomb=[[self.gene1.allele1, self.gene2.allele1, self.gene3.allele1],[self.gene1.allele2, self.gene2.allele2, self.gene3.allele2]]
    self.single_crossover=[[self.gene1.allele1, self.gene2.allele1, self.gene3.allele2],
                           [self.gene1.allele1, self.gene2.allele2, self.gene3.allele2],
                           [self.gene1.allele2, self.gene2.allele2, self.gene3.allele1],
                           [self.gene1.allele2, self.gene2.allele1, self.gene3.allele1]]
    self.double_crossover=[[self.gene1.allele1, self.gene2.allele2, self.gene3.allele1],[self.gene1.allele2, self.gene2.allele1, self.gene3.allele2]]
    self.gametes=self.non_recomb + self.single_crossover + self.double_crossover
    self.genotype='Gene 1: {}\nGene 2: {}\nGene 3: {}\n'.format(self.gene1.assortment,self.gene2.assortment,self.gene3.assortment)
    self.chr1=self.gene1.allele1 + self.gene2.allele1 + self.gene3.allele1
    self.chr2=self.gene1.allele2 + self.gene2.allele2 + self.gene3.allele2
  def print_gametes(self):
    a='non-recomb1:\t\t{}-{}-{}'.format(self.gene1.allele1, self.gene2.allele1, self.gene3.allele1)
    b='non-recomb2:\t\t{}-{}-{}'.format(self.gene1.allele2, self.gene2.allele2, self.gene3.allele2)
    c='single cross-over1:\t{}-{}-{}'.format(self.gene1.allele1, self.gene2.allele1, self.gene3.allele2)
    d='single cross-over2:\t{}-{}-{}'.format(self.gene1.allele1, self.gene2.allele2, self.gene3.allele2)
    e='single cross-over3:\t{}-{}-{}'.format(self.gene1.allele2, self.gene2.allele2, self.gene3.allele1)
    f='single cross-over4:\t{}-{}-{}'.format(self.gene1.allele2, self.gene2.allele1, self.gene3.allele1)
    g='double cross-over1:\t{}-{}-{}'.format(self.gene1.allele1, self.gene2.allele2, self.gene3.allele1)
    h='double cross-over2:\t{}-{}-{}'.format(self.gene1.allele2, self.gene2.allele1, self.gene3.allele2)
    return '{}\n{}\n{}\n{}\n{}\n{}\n{}\n{}\n'.format(a,b,c,d,e,f,g,h)
  def reorder(self,order):
    tmp1=self.gene1
    tmp2=self.gene2
    tmp3=self.gene3
    if order[0] == 2: self.gene2=tmp1
    if order[0] == 3: self.gene3=tmp1
    if order[1] == 1: self.gene1=tmp2
    if order[1] == 3: self.gene3=tmp2
    if order[2] == 1: self.gene1=tmp3
    if order[2] == 2: self.gene2=tmp3

class Two_factor_unlinked:
  def __init__(self, assortment1, assortment2):
    self.gene1 = assortment1
    self.gene2 = assortment2
    self.non_recombs=[[self.gene1.allele1, self.gene2.allele1],[self.gene1.allele2, self.gene2.allele2]]
    self.recombinants=[[self.gene1.allele1, self.gene2.allele2],[self.gene1.allele2, self.gene2.allele1]]
    self.gametes=self.non_recombs + self.recombinants
    self.genotype='{}, {}'.format(self.gene1.assortment, self.gene2.assortment)
    self.chr1=self.gene1.allele1 + self.gene2.allele1
    self.chr2=self.gene1.allele2 + self.gene2.allele2
  def print_gametes(self):
    a='non-recomb1:\t{}-{}'.format(self.gene1.allele1, self.gene2.allele1)
    b='non-recomb2:\t{}-{}'.format(self.gene1.allele2, self.gene2.allele2)
    c='recombinant1:\t{}-{}'.format(self.gene1.allele1, self.gene2.allele2)
    d='recombinant2:\t{}-{}'.format(self.gene1.allele2, self.gene2.allele1)
    return '{}\n{}\n{}\n{}\n'.format(a,b,c,d)
  def __eq__(self, other):
    return (self.gene1 == other.gene1 and self.gene2 == other.gene2)
  def cross(self, mate):
    zygotes=[]
    products=[]
    for i in range(0,4):
      for j in range(0,4):
        products.append([[self.gametes[i][0], mate.gametes[j][0]], [self.gametes[i][1], mate.gametes[j][1]]])
    for product in products: ## create Two_factor_unlinked objects
      new_gene1 = Assortment(product[0][0], product[0][1], dominant_phenotype=self.gene1.dominant_phenotype, recessive_phenotype=self.gene1.recessive_phenotype)
      new_gene2 = Assortment(product[1][0], product[1][1], dominant_phenotype=self.gene2.dominant_phenotype, recessive_phenotype=self.gene2.recessive_phenotype)
      zygote = Two_factor_unlinked(new_gene1, new_gene2)
      zygotes.append(zygote)
    return Two_factor_cross_progeny(zygotes)

class Two_factor_cross_progeny:
  def __init__(self, zygotes):
    self.zygotes = zygotes
    self.genotype_counts = []
    self.genotypes = []
    for zygote in zygotes: # populate genotypes with only unique zygotes and find the genotype ratio
      alreadyThere = False
      for i,genotype in enumerate(self.genotypes): # check genotypes to see if current zygote should be added
        if zygote == genotype:
          alreadyThere = True
      if alreadyThere: # if this zygote is already in genotypes
        self.genotype_counts[i]+=1 # increment count of this zygote
      else:
        self.genotypes.append(zygote) # this zygote isn't represented in genotypes. Add it and initialize its count in the ratio to 1.
        self.genotype_counts.append(1)
    self.zygosities = {}
    for zygote in zygotes:
      zygosity = zygote.gene1.zygosity + ', ' + zygote.gene2.zygosity
      if zygosity in self.zygosities.keys(): # if zygosity is already in dictionary
        self.zygosities[zygosity]+=1 # increment count of zygosity
      else:
        self.zygosities[zygosity]=1
    self.phenotypic_ratio = {}
    for zygote in zygotes:
      phenotype = ''
      if zygote.gene1.zygosity == 'Homozygous Recessive':
        phenotype += zygote.gene1.recessive_phenotype + ', '
      else: phenotype = zygote.gene1.dominant_phenotype + ', '
      if zygote.gene2.zygosity == 'Homozygous Recessive':
        phenotype += zygote.gene2.recessive_phenotype
      else: phenotype += zygote.gene2.dominant_phenotype
      if phenotype in self.phenotypic_ratio.keys(): # if phenotype is already in dictionary
        self.phenotypic_ratio[phenotype]+=1 # increment count of phenotype
      else:
        self.phenotypic_ratio[phenotype]=1
  def print_ratio(self): # genotypic ratio
    output_string=''
    for i in range(len(self.genotype_counts)):
      output_string+=self.genotypes[i].genotype + ': ' + str(self.genotype_counts[i]) + '\n'
    return output_string

class Three_factor_cross_progeny:
  def __init__(self, zygotes):
    self.zygotes = zygotes
    self.ratio = []
    self.genotypes = []
    for zygote in zygotes: # populate genotypes with only unique zygotes and find the genotype ratio
      alreadyThere = False
      for i,genotype in enumerate(self.genotypes): # check genotypes to see if current zygote should be added
        if zygote == genotype:
          alreadyThere = True
      if alreadyThere: # if this zygote is already in genotypes
        self.ratio[i]+=1 # increment count of this zygote in the ratio
      else:
        self.genotypes.append(zygote) # this zygote isn't represented in genotypes. Add it and initialize its count in the ratio to 1.
        self.ratio.append(1)
