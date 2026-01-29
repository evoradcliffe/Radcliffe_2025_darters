
#Vector describing the frequencies of the five female genotypes in orangethroat
# XXaa XXa'a XXa'a' XYa'a' YYa'a'
# The first one should be heavily weighted (XX on chromosome 9, aa on chromosome 23)
# the frequencies of the others could depend on gene flow, etc. 
# E.g. if there's gene flow, there might be individuals from the second genotype
female_genos_ORANGE <- c(0.95, 0.05, 0.0, 0, 0)


#Vector describing the frequencies of the three male genotypes in orangethroat
# XYaa XYa'a YYa'a
## The first one should be heavily weighted (XX on chromosome 9, aa on chromosome 23), 
## but gene flow could put some weight on the second
male_genos_ORANGE <- c(0.95, 0.05, 0)

#Vector describing the frequencies of the five female genotypes in rainbow
# XXaa XXa'a XXa'a' XYa'a' YYa'a' 
# the last genotype will be the most common, but genotype 4 might be around due to gene flow
female_genos_RAINBOW <- c(0, 0, 0, 0.1, 0.9)

#Vector describing the frequencies of the three male genotypes in rainbow
# XYaa XYa'a YYa'a
#the last one will be most common
male_genos_RAINBOW <- c(0, 0.1, 0.9)


female_genos_ORANGE <- c(0.8851, 0.0629, 0.0500, 0.0010, 0.0010)
male_genos_ORANGE   <- c(0.9371, 0.0629, 0.0000)

female_genos_RAINBOW <- c(0.0000, 0.0000, 0.0000, 0.14, 0.86)
male_genos_RAINBOW   <- c(0.01, 0.15, 0.84)

#the cross matrix describes the sex ratio arising from all the possible matings 
#of the different genotypes

cross_matrix <- matrix(NA, nrow = 5, ncol = 3)
cross_matrix[1,1] = 0.5
cross_matrix[1,2] = 0.5
cross_matrix[1,3] = 1

cross_matrix[2,1] = 0.5
cross_matrix[2,2] = 3/8
cross_matrix[2,3] = 3/4

cross_matrix[3,1] = 0.5
cross_matrix[3,2] = 1/4
cross_matrix[3,3] = 3/4

cross_matrix[4,1] = 3/4
cross_matrix[4,2] = 3/8
cross_matrix[4,3] = 0.5

cross_matrix[5,1] = 1
cross_matrix[5,2] = 0.5
cross_matrix[5,3] = 0.5


sex_ratio_female_ORANGE_male_RAINBOW = sum(outer(female_genos_ORANGE, male_genos_RAINBOW) * cross_matrix)
sex_ratio_female_RAINBOW_male_ORANGE = sum(outer(female_genos_RAINBOW, male_genos_ORANGE) * cross_matrix)


#sex_ratio_female_ORANGE_male_RAINBOW is ~0.85 
#and sex_ratio_female_RAINBOW_male_ORANGE is ~0.95

sex_ratio_female_ORANGE_male_RAINBOW
sex_ratio_female_RAINBOW_male_ORANGE

