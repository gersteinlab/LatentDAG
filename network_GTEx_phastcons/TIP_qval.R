library("qvalue")

args = commandArgs(trailingOnly=TRUE)

df = read.table(file = args[1], sep="\t", header = T)
qresult = qvalue(df$p.value, fdr.level=0.05, lambda=0.1)
df$qvalue = qresult$qvalues
df2 = df[qresult$significant, ]
write.table(df2, file = args[2], sep="\t", quote=F, row.names = F, col.names = T)

