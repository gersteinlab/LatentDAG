peak2weight<-function(peakfilename, annotationfilename, width=10000, smooth=T, outfilename)
{
    chr.len = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
                145138636,138394717,133797422,135086622,133275309,114364328,107043718,
                101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
                156040895,57227415,16569)
    chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
    chr.len = chr.len[1:(length(chr.len)-2)]
    chr.nam = chr.nam[1:(length(chr.nam)-2)]

    # read in gene info
    mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
    mygene = mygene[, 1:5]
    colnames(mygene) = c("name", "chr", "str", "sta", "end")
    mygene = unique(mygene)
    mygene$str[mygene$str == "-1"] = "-"
    mygene$str[mygene$str == "1"] = "+"
    mygene$chr = paste("chr", mygene$chr, sep="")
    mygene$chr[mygene$chr == "chrMT"] = "chrM"

    tmp = mygene$sta
    tmp1 = tmp-width
    tmp2 = tmp+width

    tag = mygene[,3]=="-"
    tmp[tag==1] = mygene[tag==1,5]
    tmp1[tag==1] = tmp[tag==1]+width 
    tmp2[tag==1] = tmp[tag==1]-width 
    mygene$sta = tmp1
    mygene$end = tmp2
    for(k in 1:length(chr.nam))
    {
      tag = mygene$chr==chr.nam[k]
      tmp = mygene$sta[tag==1]
      tmp[tmp<1] = 1
      tmp[tmp>chr.len[k]] = chr.len[k]
      mygene$sta[tag==1] = tmp
      tmp = mygene$end[tag==1]
      tmp[tmp<1] = 1
      tmp[tmp>chr.len[k]] = chr.len[k]
      mygene$end[tag==1] = tmp
    }

    # read in peak file
    data = read.table(peakfilename, sep="\t", header=F)
    colnames(data) = c("chr", "start", "end", "c4", "c5", "c6", "height", "c8", "c9", "max_ht_offset")

    ## create weight file
    myw = rep(0, width*2+1)
    for(k in 1:length(chr.nam))
    {
      cat("\r", chr.nam[k])
      # signal vector
      read.cov = rep(0, chr.len[k])
      tag = data$chr==chr.nam[k]
      if(sum(tag)==0) next
      mysig = data[tag==1,]
      for(i in 1:nrow(mysig))
      {
        tmp1 = as.numeric(mysig$start[i])+1   # start line
        tmp2 = as.numeric(mysig$end[i])   # end line
        read.cov[tmp1:tmp2] = read.cov[tmp1:tmp2] + as.numeric(mysig$height[i])
      }
      # refseq
      curgene = mygene[mygene[,2]==chr.nam[k],]
      for(i in 1:nrow(curgene))
      {
        myw = myw+read.cov[curgene[i,4]:curgene[i,5]]
      }
    }

    myw = myw/nrow(mygene)
    tmp = myw
    if(smooth==T)
    {
      for(i in 1:length(myw))
      {
        myw[i] = mean(tmp[max(i-250,1):min(i+250, length(myw))])
      }
    }
    myw = myw/sum(myw)
    write.table(myw, outfilename, sep="\t", row.names=F, col.names=F, quote=F)
}

calscoreonpeak<-function(peakfilename, annotationfilename, weightfilename, outfilename)
{
    chr.len = c(248956422,242193529,198295559,190214555,181538259,170805979,159345973,
                145138636,138394717,133797422,135086622,133275309,114364328,107043718,
                101991189,90338345,83257441,80373285,58617616,64444167,46709983,50818468,
                156040895,57227415,16569)
    chr.nam = paste("chr", c(1:22, "X", "Y", "M"), sep="")
    chr.len = chr.len[1:(length(chr.len)-2)]
    chr.nam = chr.nam[1:(length(chr.nam)-2)]

    # read in gene info
    mygene = read.table(annotationfilename, sep="\t", header=T, colClasses=c(rep("character", 3), rep("numeric", 5)))
    mygene = mygene[, 1:5]
    colnames(mygene) = c("name", "chr", "str", "sta", "end")
    mygene = unique(mygene)
    mygene$str[mygene$str == "-1"] = "-"
    mygene$str[mygene$str == "1"] = "+"
    mygene$chr = paste("chr", mygene$chr, sep="")
    mygene$chr[mygene$chr == "chrMT"] = "chrM"

    # read in weight
    conIn = file(weightfilename, "r")
    myw = readLines(conIn)
    close(conIn)
    myw = as.numeric(myw)
    width = (length(myw)-1)/2

    # read in peak file
    data = read.table(peakfilename, sep="\t", header=F)
    colnames(data) = c("chr", "start", "end", "c4", "c5", "c6", "height", "c8", "c9", "max_ht_offset")

    ## calculate score
    mysco = rep(0, nrow(mygene))
    gname = rep("",nrow(mygene))
    count = 0
    for(k in 1:length(chr.nam))
    {
      cat("\r", chr.nam[k])
      # signal vector
      read.cov = rep(0, chr.len[k])
      tag = data$chr==chr.nam[k]
      if(sum(tag)==0) next
      mysig = data[tag==1,]
      for(i in 1:nrow(mysig))
      {
        tmp1 = as.numeric(mysig$start[i])+1   # start line
        tmp2 = as.numeric(mysig$end[i])   # end line
        read.cov[tmp1:tmp2] = read.cov[tmp1:tmp2] + as.numeric(data$height[i])
      }
      # refseq
      curgene = mygene[mygene[,2]==chr.nam[k],]
      for(i in 1:nrow(curgene))
      {
        count = count + 1
        if(curgene$str[i]=="+")
        {
          tmp= read.cov[max(1,curgene$sta[i]-width):min(curgene$sta[i]+width, chr.len[k])]
          b1 = 1-(curgene$sta[i]-width)
          b2 = curgene$sta[i]+width-chr.len[k]
          if(b1>0)
          {
            tmp=c(rep(0, b1), tmp)
          }
          if(b2>0)
          {
            tmp=c(tmp, rep(0, b2))
          }
          mysco[count] = sum(tmp*myw)
        }
        if(curgene$str[i]=="-")
        {
          tmp= read.cov[min(chr.len[k],curgene$end[i]+width):max(curgene$end[i]-width, 1)]
          b1 = curgene$end[i]+width-chr.len[k]
          b2 = 1-(curgene$end[i]-width)
          if(b1>0)
          {
            tmp=c(rep(0, b1), tmp)
          }
          if(b2>0)
          {
            tmp=c(tmp, rep(0, b2))
          }
          mysco[count] = sum(tmp*myw)
        }
        gname[count] = curgene[i,1]
      }
    }
    zscore = (mysco-mean(mysco))/sd(mysco)
    pvalue = pnorm(-zscore)
    res = cbind(gname, mysco, zscore, pvalue)
    colnames(res) = c("name", "raw.score", "zscore", "p.value")
    res = res[order(res[,1], decreasing=T), ]
    write.table(res, outfilename, sep="\t", row.names=F, quote=F)
}

args = commandArgs(trailingOnly=TRUE)
annotationfilename = args[1]
peakfilename = args[2]
weightfilename = args[3]
outfilename = args[4]

peak2weight(peakfilename, annotationfilename, 10000, T, weightfilename)
calscoreonpeak(peakfilename, annotationfilename, weightfilename, outfilename)