################
### LIBRARIES ###
#################

library("data.table")

# Statistical Wilcoxon test. We use this one because the distributions are not normal
plot_marco = fread("_tmp_XXXrefdirXXX_XXXdirXXX")

plot_split = split(plot_marco, factor(plot_marco$V1))
ind1=1
ind2=2

sample1 = plot_split[[ind1]]$V2
sample2 = plot_split[[ind2]]$V2
aa = wilcox.test(sample1, sample2, alternative = "two.sided")

s=sprintf("%d %d %f %f %d %d %f %f %f %f",ind1,length(sample1),mean(sample1),sd(sample1),ind2,length(sample2),mean(sample2),sd(sample2),sqrt((mean(sample1)-mean(sample2))*(mean(sample1)-mean(sample2))),aa$p.value)
print(s)
