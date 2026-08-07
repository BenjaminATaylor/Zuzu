quasis = "[/Users/benjamin/Repositories/Zuzu/work/db/3d9529d0234c11da5fca650d53402c/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/91/6a22f7aaf04ef54b63e513d6b589b8/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/03/03bfc3b0fbf3692096077d8257c3cb/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/a7/201e264b1522f33caf245fd002d01b/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/e6/ac5d38c26dd9ca53a327daa3f111f8/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/89/f823f879da0aacbc546df4b1f87652/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/0c/e998761031ad2c63387fba1c5a28ae/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/57/a904310dd89340003d24ab32d7b66f/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/5b/f7329b77058a4e69c719a4916d7e02/outframe.csv, /Users/benjamin/Repositories/Zuzu/work/38/14fc9c3a349a65f4829e9112d83222/outframe.csv]"

deseq.inlist = str_remove_all($deseq_quasis,"[\\[\\] ]") %>% 
  strsplit(split = ",") %>% unlist()

deseq.quasi.input = 
  sapply(deseq.inlist, read.csv) %>% 
  t() %>% 
  `row.names<-`(NULL) %>% 
  data.frame() %>%
  mutate_all(as.numeric) %>%
  mutate(method = "DESeq2") %>% 
  melt(id.vars = "method")

deseq.inlist = str_remove_all($edger_quasis,"[\\[\\] ]") %>% 
  strsplit(split = ",") %>% unlist()

edger.quasi.input = 
  sapply(deseq.inlist, read.csv) %>% 
  t() %>% 
  `row.names<-`(NULL) %>% 
  data.frame() %>%
  mutate_all(as.numeric) %>%
  mutate(method = "edgeR") %>% 
  melt(id.vars = "method")

gg.quasi.input = rbind(deseq.quasi.input, edger.quasi.input)


ggplot(permhist.in, aes(x = value)) +
  geom_histogram(color = "black", fill = "white", bins = 100) +
  scale_x_continuous(limits = c(1,(nreps+1)), breaks = cutbreaks) +
  labs(x = "Number of permutations in which gene was\nincorrectly identified as a DEG",
       y = "Number of genes") +
  theme_bw() +
  facet_grid(~Var2) +
  theme(strip.background = element_blank())


ggplot(gg.quasi.input, aes(x = method, y = value)) +
  geom_point(size = 3, alpha = 0.7) +
  scale_y_continuous(limits = c(0,1)) +
  labs(x = "Method",y = "Value") +
  facet_grid(~variable) +
  theme_bw() +
  theme(strip.background = element_blank())

  
