Steps to identify non-random enrichment in association between two SNPs.
1. Ancestry-search
Ancestry search - Program executed to determine the ancestry from each SNP or
variant available in our dataset from each Caribbean Hispanic patient in treatment with clopidogrel. The resulting ancestry fractions were used for further steps in script development. Output: A cvs file was obtained with the following columns: rsnumber, chr (chromosome), pos (position), cM (centimorgan), average ancestry, and samples.

2. Bracket search
Bracket search - Program executed to determine the matching positions between variants associated to cardiovascular conditions and variants from clopidogrel resistance in Caribbean Hispanics. In addition, it determined if the matching positions were  0.5 cM from each other. The ancestry was determined for the resulting variants. Output: three cvs files were generated. The first one indicated the positions that were on our list of available variants with ancestry. The second csv specified the variants that were +- 0.5 cM from each other. The third csv had information about the ancestry from matching variants.

3. Program executed to determine random SNPs, excluding the variants matching positions described previously in step 2.Output: A csv file is generated containing N random positions with its corresponding ancestry.

