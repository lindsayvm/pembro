#######FUNCTIONS for plots
my_barplot <- function(df, inputCol, sort=F) {
  #make barplot with default that it does not sort unless specified
  
  #sort on frequence
  freq.t = table(inputCol)
  freq_sorted.t = sort(freq.t, decreasing=TRUE)
  if(sort == T){
    inputCol = factor(inputCol, 
                      levels = names(freq_sorted.t))
  }
  
  #factorize
  inputCol = factor(inputCol, levels=names(table(inputCol)))
  
  #freq
  t = plyr::count(inputCol)
  
  #colors
  if(length(unique(inputCol)) <= 3){
    values = c("brown",  "orange", "skyblue")
    
    #values = c("#2172b6","#EBA924")
    
  } else {
    values = randomColor(length(unique(inputCol)), luminosity = "light")
  }
  
  #barplot
  p2 =  df %>% ggplot(aes(x = inputCol, fill = inputCol)) + 
    
    geom_bar() + 
    
    theme_bw() +
    
    theme(axis.text.x = element_text(size = 16, 
                                     angle = 45, 
                                     vjust = 1, 
                                     hjust = 1),
          axis.text.y = element_text(size = 12)) +
    
    theme(axis.title.x=element_blank(), 
          axis.title.y=element_blank(), 
          legend.position = "none") +
    
    scale_fill_manual(values = values) 
  
  return(p2)
}



my_heatmap <- function(clin.df, var.df, var.m, survivalCol) {
  (clin.df$patientID)
  colnames(var.df)[-1]
  my_group = as.factor(clin.df[[survivalCol]]) #blue is clinical benefit
  colSide = brewer.pal(length(unique(my_group)), "Set1")[my_group]
  heatmap(var.m,
          Colv = TRUE,
          Rowv = TRUE,
          labRow = var.df$gene,
          labCol = NA,
          ColSideColors=colSide)
}



my_heatmap2 <- function(clin.df, var.m, survivalCol) {
  my_group = as.factor(clin.df[[survivalCol]]) #blue is clinical benefit
  colSide = brewer.pal(length(unique(my_group)), "Set1")[my_group]
  heatmap(var.m,
          Colv = TRUE,
          Rowv = TRUE,
          labCol = NA,
          ColSideColors=colSide)
}





















#########PLAYGROUND







# 
# 
# 
# 
# #distribution
# my_barplot <- function(df, inputCol, sort=F) {
#   #make barplot with default that it does not sort unless specified
#   
#   #sort on frequence
#   freq.t = table(inputCol)
#   freq_sorted.t = sort(freq.t, decreasing=TRUE)
#   if(sort == T){
#     inputCol = factor(inputCol, 
#                       levels = names(freq_sorted.t))
#   }
#   
#   #factorize
#   inputCol = factor(inputCol, levels=names(table(inputCol)))
#   
#   #freq
#   t = plyr::count(inputCol)
#   
#   #colors
#   if(length(unique(inputCol)) <= 3){
#     #values = c("#2c7a58", "salmon", "grey") #seagreen 
#     values = c("#2172b6","#EBA924")
#     
#   } else {
#     values = randomColor(length(unique(inputCol)), luminosity = "light")
#   }
#   
#   #barplot
#   p2 =  df %>% ggplot(aes(x = inputCol, fill = inputCol)) + 
#     
#     geom_bar() + 
#     
#     theme_bw() +
#     
#     theme(axis.text.x = element_text(size = 16, 
#                                      angle = 45, 
#                                      vjust = 1, 
#                                      hjust = 1),
#           axis.text.y = element_text(size = 12)) +
#     
#     theme(axis.title.x=element_blank(), 
#           axis.title.y=element_blank(), 
#           legend.position = "none") +
#     
#     scale_fill_manual(values = values) +
#     ylim(0, 125)
#   
#   return(p2)
# }
# 
# 
# 
# my_cb_stackedbarplot <- function(df, xlab, pval=NA){
#   
#   df$cb = ifelse(df$cb == "CB", "DCB", "NCB")
#   df$cb = factor(df$cb, levels = c("NCB", "DCB"))
#   df$n = as.numeric(df$n)
#   df$perc = as.numeric(df$perc)
#   df$ypos = (df$perc/100) + 0.05 
#   perc= df %>% filter(cb == "DCB")
#   
#   total.df = df %>% 
#     dplyr::group_by(status) %>% 
#     dplyr::summarise(total = sum(n)) %>% 
#     dplyr::left_join(perc, by = "status")
#   
#   xlab = toupper(gsub("Lab_BL_|_clinical|_bool|wgs_|rna_|sig_|_RESPONSE|_","",xlab))
#   xlab = gsub("INTERFERONGAMMA",paste0("IFN-","\U03B3"),xlab)
#   xlab = gsub("ALBUMIN",paste0("Albumin"),xlab)
#   
#   
#   total.df$pval = c(paste0("P=",   format(pval, digits = 3)),
#                     rep("", nrow(total.df)-1))
#   if(nrow(total.df) == 2){
#     hjust= 0.1
#   }else{
#     hjust= -0.8
#   }
#   
#   ggplot(df, aes(fill=cb, 
#                  y=n, 
#                  x=status)) + 
#     
#     geom_bar(position="fill", 
#              stat="identity") +
#     
#     scale_fill_manual(values = c("salmon", "#2c7a58"))+
#     
#     theme_bw(base_size = 15)+
#     
#     theme(legend.position = "None",
#           #        axis.text.x = element_text(size = 16, angle = 45, vjust = 1.5, hjust=0.5),
#           axis.text.x = element_text(size = 17, angle = 45, vjust = 1, hjust=1),
#           #    axis.text.y = element_text(size = 10, angle = 0, vjust = 1, hjust=1),
#           axis.text.y = element_blank(),
#           axis.title.y = element_blank(),
#           axis.title.x = element_text(size = 20),
#           panel.border = element_blank(), 
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           axis.ticks = element_blank()) +
#     
#     scale_y_continuous(expand=expansion(mult = c(0.1, 0.1))) +
#     
#     xlab(xlab) + 
#     
#     #count
#     geom_text(inherit.aes = FALSE, data = total.df, 
#               aes(x = status, y = rep(1, nrow(total.df)), 
#                   label = paste0(total)), 
#               vjust = 19.5,#-0.3,
#               size = 6) +
#     
#     #perc
#     geom_text(inherit.aes = T, data = total.df,
#               aes(x = status, y = ypos, group = status,
#                   label = paste0(perc,"%")), size = 4) +
#     
#     #pval
#     geom_text(inherit.aes = T, data = total.df,
#               aes(x = status, y = 1.05, hjust = hjust,  
#                   label = pval), size = 6) 
#   #    geom_signif(annotations = 0.03, y_position = 0.5 ,xmin="High", xmax="Low")
#   
#   
# }
# 
# 
# my_wilcoxon <- function(df, VAR){
#   
#   df$Clinical_benefit = as.factor(ifelse(df$PFS_immuno_6mo == 0,"nonCB","CB"))
#   
#   df = df[!is.na(df[[VAR]]), ]
#   
#   p = df %>%  ggplot( aes(x=Clinical_benefit, 
#                           y=df[[VAR]], 
#                           fill=Clinical_benefit)) +
#     geom_boxplot(outlier.shape = 18) +
#     scale_fill_manual(values = c("skyblue", "orange")) +
#     geom_jitter(color="black", size=0.5, alpha=0.3) +
#     theme_minimal(base_size = 20)+
#     theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.4, hjust=0.2)) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 20),
#           legend.position = "none") +
#     ylab(VAR)
#   
#   return(p)
# }  
# 
# my_BOR <- function(df, VAR){
#   
#   df$BOR_treatment1_naCPCT = factor(df$BOR_treatment1_naCPCT, 
#                                     levels=c("PD","SD","PR","CR")) ###!!"NE"
#   df %>%  ggplot( aes(x=BOR_treatment1_naCPCT, y=df[[VAR]], fill=BOR_treatment1_naCPCT)) +
#     geom_boxplot(outlier.shape = 18) +
#     scale_fill_manual(values = c("darkred","orange","seagreen", "skyblue")) + ##!!!"grey"
#     geom_jitter(color="black", size=0.5, alpha=0.3) +
#     theme_minimal(base_size = 12)+
#     theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.4, hjust=1)) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 10),
#           legend.position = "none") +
#     ylab(VAR)
# }
# my_violin <- function(df, VAR){
#   
#   df$Clinical_benefit = as.factor(ifelse(df$PFS_immuno_6mo == 0,"nonCB","CB"))
#   
#   df = df[!is.na(df[[VAR]]), ]
#   
#   p = df %>%  ggplot( aes(x=Clinical_benefit, 
#                           y=df[[VAR]], 
#                           fill=Clinical_benefit)) +
#     geom_violin(outlier.shape = 18) +
#     scale_fill_manual(values = c("skyblue", "orange")) +
#     geom_jitter(color="black", size=0.5, alpha=0.3) +
#     theme_minimal(base_size = 20)+
#     theme(axis.text.x = element_text(size = 20, angle = 45, vjust = 0.4, hjust=0.2)) +
#     theme(axis.title.x=element_blank(), axis.title.y=element_text(size = 20),
#           legend.position = "none") +
#     ylab(VAR)
#   
#   return(p)
# }  
# 
# 
# 
# my_pieplot <- function(df, inputCol) {
#   #factorize
#   inputCol = factor(inputCol, 
#                     levels = c("NE","PD","SD","PR","CR"))
#   #Calc freq
#   t = plyr::count(inputCol)
#   prob = t$freq/sum(t$freq)
#   print(round(prob,2))
#   #colors = randomColor(length(unique(inputCol)), luminosity = "light")
#   colors = c("lightgrey", "#4292C6", "#2171B5", "#08519C", "#08306B")
#   pie(prob, t$x, border="white", col= colors,init.angle=180 )
# }
# 
# my_ringplot <- function(inputCol) {
#   #Freq table different hospital centers
#   final.df = as.data.frame(plyr::count(inputCol))
#   
#   #Order freq table
#   final.df = final.df[order(final.df$freq, decreasing = F), ]
#   
#   #Change col and row names
#   colnames(final.df) = c("groups","count")
#   rownames(final.df) = c()
#   
#   # Compute percentages
#   final.df$fraction = final.df$count/sum(final.df$count)
#   
#   # Compute the cumulative percentages (top of each rectangle)
#   final.df$ymax = cumsum(final.df$fraction)
#   
#   # Compute the bottom of each rectangle
#   final.df$ymin = c(0, head(final.df$ymax, n=-1))
#   
#   # Compute label position
#   final.df$labelPosition = c(0.25,0.75)
#   
#   # Compute a good label
#   final.df$label = paste0(final.df$groups, "\n n = ", final.df$count)
#   
#   # Make the circl plot
#   final.df = final.df %>%
#     mutate(groups = forcats::fct_inorder(groups))
#   
#   ggplot(final.df, aes(ymax=ymax, ymin=ymin, 
#                        xmax=4, xmin=3, 
#                        fill=groups)) +
#     geom_rect() +
#     geom_label( x=3.5, aes(y=labelPosition, label=label), size=7) +
#     theme_void() +
#     coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#     xlim(c(2.5, 4)) + # Try to remove that to see how to make a pie chart
#     scale_fill_manual(values = c("salmon","#2c7a58")) +
#     theme(legend.position = "none")
# }
# 
# my_ringplot2 <- function(inputCol) {
#   #Freq table different hospital centers
#   final.df = as.data.frame(plyr::count(inputCol))
#   
#   #Change col and row names
#   colnames(final.df) = c("groups","count")
#   rownames(final.df) = c()
#   
#   # Compute percentages
#   final.df$fraction = final.df$count/sum(final.df$count)
#   
#   # Compute the cumulative percentages (top of each rectangle)
#   final.df$ymax = cumsum(final.df$fraction)
#   
#   # Compute the bottom of each rectangle
#   final.df$ymin = c(0, head(final.df$ymax, n=-1))
#   
#   # Compute label position
#   final.df$labelPosition = (final.df$ymax + final.df$ymin) / 2
#   
#   # Compute a good label
#   final.df$label = paste0(final.df$groups, "\n n = ", final.df$count)
#   
#   # Make the circl plot
#   final.df = final.df %>%
#     mutate(groups = forcats::fct_inorder(groups))
#   
#   ggplot(final.df, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=groups)) +
#     geom_rect() +
#     geom_text( x=2.3, aes(y=labelPosition, label=label, color=groups), size=4) + # x here controls label position (inner / outer)
#     scale_fill_brewer(palette=3) +
#     scale_color_brewer(palette=3) +
#     coord_polar(theta="y") +
#     xlim(c(1.2, 4)) +
#     theme_void() +
#     theme(legend.position = "none")
# }
# 
# 
# #Fit the model
# my_survplot <- function(df, fit) {
#   surv.p = ggsurvplot(
#     fit,
#     data = df,
#     size = 1,                 # change line size
#     palette = c("#EBA924","#2172b6" ,"salmon","seagreen"),# custom color palettes
#     conf.int = F,          # Add confidence interval
#     pval = TRUE,              # Add p-value
#     risk.table = TRUE,        
#     risk.table.col = "strata",# Risk table color by groups
#     # legend.labs =  c("Male", "Female"),    # Change legend labels
#     risk.table.height = 0.25, 
#     ggtheme = theme_light(),  #theme_ipsum    
#     xlab = xlab,
#     xlim = c(0, x_axis),
#     break.time.by = breaktime,
#     surv.median.line = "hv",
#     risk.table.y.text = F,
#     legend.title=""
#   )
#   return(surv.p)
# }
# x_axis = 60
# breaktime = 6
# xlab = "Time in months"
# 
# 
# 


# 
# 
# 
# 
# makeProfilePlot <- function(mylist,names)
# {
#   require(RColorBrewer)
#   # find out how many variables we want to include
#   numvariables <- length(mylist)
#   # choose 'numvariables' random colours
#   colours <- brewer.pal(numvariables,"Set1")
#   # find out the minimum and maximum values of the variables:
#   mymin <- 1e+20
#   mymax <- 1e-20
#   for (i in 1:numvariables)
#   {
#     vectori <- mylist[[i]]
#     mini <- min(vectori)
#     maxi <- max(vectori)
#     if (mini < mymin) { mymin <- mini }
#     if (maxi > mymax) { mymax <- maxi }
#   }
#   # plot the variables
#   for (i in 1:numvariables)
#   {
#     vectori <- mylist[[i]]
#     namei <- names[i]
#     colouri <- colours[i]
#     print(c(namei, colouri))
#     if (i == 1) { plot(vectori,col=colouri,type="l",ylim=c(mymin,mymax)) }
#     else         { points(vectori, col=colouri,type="l")                                     }
#     lastxval <- length(vectori)
#     lastyval <- vectori[length(vectori)]
#     text((lastxval-10),(lastyval),namei,col="black",cex=0.6)
#   }
# }
# 
# 
# my_forest <- function(df, survival, method) {
#   
#   #plot
#   ggplot(data = df, aes(y = Biomarkers, x = OR)) + 
#     geom_point() +
#     # xlim(0,10) +
#     # coord_trans(x = "log")+
#     scale_x_log10() + 
#     # scale_x_continuous(trans = "log2")+
#     geom_errorbarh(aes(xmax = conf.high, xmin = conf.low)) +
#     geom_vline(xintercept = 1, 
#                color = "lightblue", size=1, alpha = 0.6) +
#     facet_grid(groups ~ ., scales = "free_y", space = "free") +
#     # scale_y_discrete(breaks=df$Biomarkers, #https://stackoverflow.com/questions/44013753/conditional-formatting-of-axis-text-in-faceted-ggplot
#     #                   labels=df$Biomarkers) + #or: https://stackoverflow.com/questions/44013753/conditional-formatting-of-axis-text-in-faceted-ggplot
#     theme(plot.title = element_text(hjust = 0.5, size=16,face="bold"),
#           axis.title.y = element_blank(),
#           axis.text = element_text(size=16),
#           axis.title=element_text(size=16),
#           strip.text.y = element_blank(),
#           strip.background = element_rect(colour="white", fill="white"),
#           panel.background = element_rect(fill = "white",
#                                           colour = "black",
#                                           size = 0.5, linetype = "solid"),
#           panel.grid.major = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgrey"),
#           panel.grid.minor = element_line(size = 0.25, linetype = 'solid',
#                                           colour = "lightgrey")) +
#     ggtitle(survival) + 
#     xlab(paste0("Log10(",method,")")) 
# }
# 
# 
