cProperties for AIMs_generator.py should be configured as key=value pairs in a 
file called aims_properties.txt, located in the same directory as the script 
is run from. Alternatively, run the script as
   > python AIMs_generator.py <full_path_to_properties_file>

The following properties must/can be configured:

PROPERTY    | DESCRIPTION                                | MANDATORY? | DEFAULT
-------------------------------------------------------------------------------
ldfile      | Full path name of PLINK .ld file           | Yes        | None

ldthresh    | The r^2 LD threshold. If the LD r^2 value  | No         | 0.1
            | between 2 SNPs exceeds this value, the     |            |
            | second SNP (the least informative SNP)     |            |
            | will be excluded as an AIM                 |            |

posfile     | Full path name of PLINK .bim file          | Yes        | None

nraims      | A comma-delimited list of the number of    | No         | 500 
            | SNPs required in each AIM set. A value of  |            |
            | 100,200,300 would for example result in    |            | 
            | three sets, of sizes 100, 200 and 300      |            | 
            | respectively

distances   | A comma-delimited list of the minimum      | No         | 100000
            | distance required between SNPs selected    |            |
            | as AIMs, specified in number of base pair  |            |
            | positions. A value of 100000,1000000 would |            |
            | for example result in two sets being       |            |
            | generated for each AIM set specified by    |            |
            | nraims; the first set containing no SNPs   |            |
            | that are less than 100 000 base pairs      |            |
            | apart and the second set containing no     |            |
            | SNPs that are less than 1 000 000 base     |            |
            | pairs apart.                               |            |

populations | A comma-delimited list of the names of     | Yes        | None
            | the ancestral populations, e.g.            |            |
            | european,african                           |            |

<pop>.frq   | A PLINK frequency (.frq) file is required  | No         | None
            | for each of the ancestral populations      |            |
            | specified in the populations property. By  |            |
            | default it is assumed that the frequency   |            |
            | file is located in the script directory    |            |
            | with the same name as the population with  |            |
            | .frq as file extension, e.g. european.frq. |            |
            | Different frequency file names may be      |            |
            | specified using this sub-property, e.g.    |            |
            | european.frq=eur.frq                       |            |

strategy    | The strategy to use for AIMs selection.    | No         | lsbl
            | Valid values are lsbl and In. The lsbl     |            |
            | strategy can only be used with three       |            |
            | ancestral populations. The In strategy     |            |
            | uses Rosenberg's In statistic calculated   |            |
            | accross all the ancestral populations      |            |
            | to select the most informative SNPs.       |            |

propmultiIn | If the In strategy is used, the proportion | No         | 0.5
            | of SNPs that should be selected according  |            |
            | to the In statistic calculated across all  |            |
            | ancestral populations. A value of 0.5 for  |            |
            | a set of size 100 would result in the      |            |
            | first 50 SNPs being selected according to  |            |
            | this criteria. Valid values range from 0   |            |
            | to 1.                                      |            |

<pop>.sub-  | If SNPs should be excluded if they are not | No         | None
populations | homogeneously distributed in the sub-      |            |
            | populations of an ancestral population,    |            |
            | the sub-populations should be specified    |            |
            | for that population as a comma-delimited   |            |
            | list. For example, if 'european' has been  |            |
            | specified in the populations property and  |            |
            | it is comprised of subpopulations 'tsi'    |            |
            | and 'ceu', the following property should   |            |
            | be set: european.subpopulations=tsi.ceu    |            |

hetthresh   | The P value threshold used to test whether | No         | 0.01
            | a SNP is homogenously distributed          |            |
            | between the subpopulations of an ancestral |            |
            | population.                                |            |

<subpop>.frq| A PLINK frequency (.frq) file is required  | No         | None
            | for each of the populations specified in   |            |
            | each <pop>.subpopulations property. By     |            |
            | default it is assumed that the frequency   |            |
            | file is located in the script directory    |            |
            | with the same name as the subpopulation    |            |
            | with .frq as file extension, e.g. tsi.frq. |            |
            | Different frequency file names may be      |            |
            | specified using this sub-property, e.g.    |            |
            | tsi.frq=hapmap_tsi.frq                     |            |

outstem     | The file prefix to add to each .aims       | No         | aims_
            | output file that will be generated         |            |
