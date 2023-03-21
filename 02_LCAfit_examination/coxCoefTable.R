###
### Function to extract data a cox ph model.
###


### Custom function to dynamically generate cox coefs table and forest plot data
coxCoefTable<- function(cox_surv_obj,data,cutpoints=NULL,truncatecutpoints=NULL){
  
  #Generate summary of the coefficients table
  sum<-summary(cox_surv_obj) 
  
  fullsum<- merge(as.data.frame(sum$coefficients), 
                  as.data.frame(sum$conf.int[,-1]), 
                  by='row.names') %>%
    mutate(Variables=case_when(grepl("\\:strata\\(tpoint\\)",Row.names) ~ gsub("\\:strata\\(tpoint\\)","",Row.names),
                               grepl("strata\\(tpoint\\)tpoint=[0-9]+\\:",Row.names) ~ gsub("strata\\(tpoint\\)(tpoint=[0-9]+)\\:(.*$)","\\2\\1",Row.names),
                               TRUE ~ Row.names),
           Variables=if_else(grepl("`",Variables),gsub("`","",Variables),Variables) #Remove any backticks in variable names
    ) %>%
    select(-Row.names)
  
  ### Generate a summary table which includes reference levels for cat vars
  xlevs <- unlist(cox_surv_obj$xlevels) #unlist the variable levels from the model
  coxres<- data.frame(vars= gsub("[0-9]+$","",names(xlevs)), #Construct data frame with variables and levels
                      var_levels= unname(xlevs),
                      stringsAsFactors = FALSE
  )
  
  nonCat<- unique(gsub("tpoint=[0-9]+","",fullsum$Variables[!is.na(fullsum$coef)]))
  nonCat<-nonCat[!nonCat %in% paste0(coxres[[1]],coxres[[2]])]
  if(length(nonCat)>0){
    coxres[(nrow(coxres)+1):(nrow(coxres)+length(nonCat)),1] <- nonCat
  }
  
  #If there is time-dependence in the dataset, subset to non-time columns and and represent all levels and times in a table
  incTime<- grepl("tstart",attributes(cox_surv_obj$terms)$variables)
  if(any(incTime)){
    
    data <- na.omit(data) #Filter to complete cases, as used in analysis
    
    coxres<- coxres %>%
      mutate(Var1 = case_when(!is.na(var_levels) ~ paste0(vars,var_levels),
                              is.na(var_levels) ~ vars)
      ) %>%
      filter(!grepl("tpoint",Var1)) %>%
      distinct() %>%
      full_join(as.data.frame(table(gsub("tpoint=[0-9]+$","",fullsum$Variables))),by="Var1") %>%
      rename(Variables=Var1,Time=Freq) %>%
      mutate(Time = if_else(!is.na(Time), Time,1L),
             Variables = if_else(Time>1,paste0(Variables,"tpoint="),Variables)
      ) %>%
      slice(rep(row_number(), Time)) %>%
      group_by(Variables) %>%
      mutate(Time=if_else(grepl("tpoint=",Variables),row_number(),NA_integer_)) %>%
      rowwise() %>%
      mutate(N=case_when(!is.na(var_levels) & !is.na(Time) ~ sum(data[data$tpoint==Time,][[vars]]==var_levels), ### N's determined by N at each relevant level at each time
                         is.na(var_levels) & !is.na(Time) ~ sum(!is.na(data[data$tpoint==Time,][[vars]])),
                         !is.na(var_levels) & is.na(Time) ~ sum(data[data$tpoint==1,][[vars]]==var_levels),
                         TRUE ~ sum(!is.na(data[data$tpoint==1,][[vars]]))
      )
      ) %>%
      ungroup() %>%
      mutate(Variables=if_else(grepl("tpoint=",Variables),paste0(Variables,Time),Variables))
    
  } else {
    
    ### Filter to complete cases, as per analysis
    cols<- gsub("`","",rownames(attributes(cox_surv_obj$terms)$factors)) #identify variables used in analysis, drop any backticks
    cols[1] <- gsub("Surv\\(|\\)","",cols[1]) #Remove the survival variable, and extract the strings
    if(grepl("\\,",cols[1])){
      strBreak<- strsplit(cols[1],split=", ")[[1]]
      cols[(length(cols)+1):(length(cols)+length(strBreak))] <- strBreak 
      cols <- cols[-1]
    }
    cols<- cols[which(cols %in% colnames(data))]  #Drop to columns which have been correctly labelled
    if(length(cols)>0){data<- na.omit(data[,cols])} #Subset dataset to those in the analysis
    
    #Mutate to add Variables and N columns
    coxres <- coxres %>%
      rowwise() %>%
      mutate(N=if_else(!is.na(var_levels),sum(data[[vars]]==var_levels),
                       sum(!is.na(data[[vars]]))
      )) %>%
      ungroup() %>%
      mutate(Variables =if_else(!is.na(var_levels),paste0(vars,var_levels),vars))
    
  }
  
  #Combine the main tables and assign pretty names for variables
  coxres <- full_join(coxres,fullsum,by="Variables") %>%
    mutate(`Pr(>|z|)` = if_else(`Pr(>|z|)`< 2e-16,"<2e-16",paste0(signif(`Pr(>|z|)`,3))) #Tidy P-value summaries
           )
  
  if(any(incTime)){ #For a time-dependent analysis, drop any reference-level factors not at time 1
    
    timeRefgroup<- which(coxres$Time>1 & !is.na(coxres$var_levels) & is.na(coxres$coef))
    if(length(timeRefgroup)>0){
      coxres <- coxres[-timeRefgroup,] # Remove the T>1 reference groups
      
      timeRefgroup_t1<- which(coxres$Time==1 & !is.na(coxres$var_levels) & is.na(coxres$coef))
      coxres$Time[timeRefgroup_t1] <- NA
    }
    
    ## Generate meaningful time interval labels
    if(length(cutpoints)==1){
      TimeLabels <- c(paste0("0<:",cutpoints[1]),paste0(cutpoints[length(cutpoints)],"<"))
    } else {
      TimeLabels <- c(paste0("0<:",cutpoints[1]),
                      sapply(cutpoints[2:length(cutpoints)],function(x,cutpoints){
                        paste0(cutpoints[(which(cutpoints==x)-1)],"<:",x)
                      },cutpoints),
                      paste0(cutpoints[length(cutpoints)],"<"))
    }
    #If time is truncated, drop timepoints
    if("0ToFirst" %in% truncatecutpoints){
      TimeLabels <- TimeLabels[-1]
    }
    if("lastToEnd" %in% truncatecutpoints){
      TimeLabels <- TimeLabels[-length(TimeLabels)]
    }
    
    #Identify any time-independent levels, and create their own time frame based on full interval modelled
    timeIndep<- is.na(coxres$Time)
    if(any(timeIndep)){
      coxres$Time[is.na(coxres$Time)] <- length(TimeLabels)+1
      
      if(grepl("<$",TimeLabels[length(TimeLabels)]) || length(TimeLabels)==1){
        indepLabel<- gsub(":.*","",TimeLabels[1])
      } else {
        indepLabel<- paste0(gsub(":.*","",TimeLabels[1]),":",gsub(".*:","",TimeLabels[length(TimeLabels)]))
      }
      TimeLabels <- c(TimeLabels,indepLabel)
      
    }
    coxres$Time <- factor(coxres$Time,levels=1:length(TimeLabels),labels = TimeLabels)
  }
  
  ####
  ## Generate summary plot
  ####    
  plotdataset<- coxres %>%
    mutate(`Pr(>|z|)` = if_else(grepl("<",`Pr(>|z|)`),gsub("<","< ",`Pr(>|z|)`),paste0("= ",`Pr(>|z|)`)), #Tidy p-value signs
           Variables = case_when(!is.na(var_levels) & !is.na(`exp(coef)`) ~ paste0(var_levels,"\n[p ",`Pr(>|z|)`,"; N = ",N,"]"),
                                 !is.na(var_levels) & is.na(`exp(coef)`) ~ paste0(var_levels,"\n[reference group; N = ", N,"]"),
                                 #is.na(var_levels)  ~ paste0(vars,"\n[reference group; N = ", N,"]"),
                                 TRUE ~ paste0(vars,"\n[p ",`Pr(>|z|)`,"; N = ", N,"]")),
           `exp(coef)` = case_when(!is.na(`exp(coef)`) ~ `exp(coef)`,
                                   TRUE ~ 1
           )
    )
  
  #If using time-dependence add the time variable to variable labels and arrange descending by time and second HR Otherwise, arrange by HR only.
  if(any(incTime)){ 
    plotdataset <- plotdataset %>%
      group_by(vars) %>%
      arrange(desc(Time),desc(`exp(coef)`), .by_group=TRUE) %>% 
      rowwise() %>%
      mutate(Variables=if_else(!is.na(Time),gsub("\n",paste0(" (Time ",Time,")\n"),Variables),
                               Variables
      )) %>%
      ungroup()
  } else {
    plotdataset <- plotdataset %>%
      group_by(vars) %>%
      arrange(desc(`exp(coef)`), .by_group=TRUE) %>%
      ungroup()
  }
  
  #Fix Variable labels as factor variable, and set alternating background colours
  plotdataset <- plotdataset %>%
    mutate(VarLabels=factor(row_number(),labels=Variables),
           background=if_else(as.integer(factor(vars)) %% 2 == 0,2,1),
           background=case_when(row_number()==1 ~ 1,
                                lag(background)!=background ~ background,
                                TRUE ~ NA_real_),
           panelstrip=if_else(!is.na(background),vars,NA_character_)
    )
  
  
  plot <- ggplot(plotdataset,aes(x=`exp(coef)`,y=VarLabels,fill = factor(background)))+
    facet_grid(rows=vars(factor(vars)),scales = "free", space="free")+
    geom_point(aes(size=N),shape="diamond")+#,position = position_dodge(width=0.4))+
    geom_errorbar(aes(xmin=`lower .95`,xmax=`upper .95`),width=0.2)+#,position=position_dodge(width=0.4))+
    scale_size(guide=NULL,range=c(2.5,5))+
    geom_vline(xintercept=1,lty=2)+
    theme_minimal()+
    geom_rect(xmin = -Inf,xmax = Inf, ### Colour the panels alternately
              ymin = -Inf,ymax = Inf,alpha = 0.4) +
    geom_text(aes(x = Inf, y = Inf, label = ifelse(is.na(panelstrip), "",panelstrip)), hjust="right", vjust="top",fontface = "bold")+
    scale_fill_manual(values=c("1"="azure3","2"="white"),na.value = NA,guide='none')+
    labs(x="Hazard ratio [95% CI]")+
    theme(axis.title.y = element_blank(),
          panel.spacing.y=unit(0, "lines"),
          strip.background = element_blank(),# element_rect(fill="azure3",colour = "azure3"),
          strip.text.y =  element_blank()# element_text(margin = margin(5,5,20,5, "mm"))
    )
  
  ####
  ## Generate tabular summary
  ####  
  startcol<- which(colnames(coxres)=="coef")
  roundcols<- startcol:ncol(coxres)
  roundcols <- roundcols[!roundcols==which(colnames(coxres)=="Pr(>|z|)")] 
  
  sumry <- coxres %>% #tidy to 3 signficiant figures, then combine Coef with SE and exp coef with CI. Rearrange column names accordingly
    mutate(across(all_of(roundcols), ~signif(.,3)),
           coef=if_else(!is.na(coef),paste0(coef, " [",`se(coef)`,"]"),NA_character_),
           `exp(coef)`=if_else(!is.na(coef),paste0(`exp(coef)`, " [",`lower .95`,", ", `upper .95`,"]")
                               ,NA_character_)
    ) %>%
    rename(Variable=vars,`Factor level`=var_levels,`coef [SE]`=coef,`Hazard Ratio [95% CI]`=`exp(coef)`,
           `Inverse hazard ratio`=`exp(-coef)`,`Z-score`=z,`P-value`=`Pr(>|z|)`) %>%
    select(-c(Variables,`se(coef)`,`lower .95`, `upper .95`,`Z-score`,`P-value`),`Z-score`,`P-value`)
  
  
  return(list(tabulation=sumry,figure=plot))
} 
