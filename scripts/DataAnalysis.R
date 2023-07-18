library(tidyverse)
library(latex2exp)
library(xtable)
library(RCurl)
library(ggplot2)
library(scales)
require(gridExtra)
library(facetscales)


setwd("../Results") 

#small instances
scenarios=data.frame(InstSize="small",
            m=c(1,4),
            n=c(10),
            l=c(10),
            k=c(10),
            num_instance=50,
            time_horizon=600)

#medium instances
scenarios=rbind(scenarios,data.frame(InstSize="medium",
                                     m=c(1,4),
                                     n=c(600),
                                     l=c(15),
                                     k=c(25),
                                     num_instance=50,
                                     time_horizon=1200))


 
#large instances
scenarios=rbind(scenarios,data.frame(InstSize="large",
                                     m=c(1,4),
                                     n=c(3600),
                                     l=c(16),
                                     k=c(30),
                                     num_instance=50,
                                     time_horizon=3600))

min_time=0.1

#collecting the data files
data_OCO=c()
data_SGSP=c()
for (row in seq(1,nrow(scenarios))){
  m=scenarios$m[row]
  n=scenarios$n[row]
  l=scenarios$l[row]
  K=scenarios$k[row]
  num_instance=scenarios$num_instance[row]
  for (i in seq(1,num_instance,1)){
    print(paste(m,n,l,K,i))
    for (method in c("PDOCOTF","PDOCOFF","PDOCOTT")){
      file_name=sprintf("m=%d_n=%d_l=%d_K=%d_i=%d_%s_Iteration.txt",m,n,l,K,i,method)
      print(paste0(file_name))
      temp=read.table(file_name,quote = "",sep="",header=FALSE)
      names(temp) <- c("outer_iter","inner_iter","time","avg_obj","avg_slack_const","last_obj","last_slack_const","lb_obj","ub_obj","cur_obj")
      temp<-temp %>% mutate(m=m,n=n,l=l,K=K,rep=i,method=method)
      if (length(data_OCO)==0){
        data_OCO<-temp
      }else{
            data_OCO=rbind(data_OCO,temp)
      }
    }
    for (method in c("SGSP")){
      file_name=sprintf("m=%d_n=%d_l=%d_K=%d_i=%d_%s_Iteration.txt",m,n,l,K,i,method)
      print(paste0(file_name))
      temp=read.table(file_name,quote = "",sep="",header=FALSE)
      names(temp) <- c("outer_iter","inner_iter","time","avg_obj","avg_slack_const","last_obj","last_slack_const")
      temp<-temp %>% mutate(m=m,n=n,l=l,K=K,rep=i,method=method)
      if (length(data_SGSP)==0){
        data_SGSP<-temp
      }else{
        data_SGSP=rbind(data_SGSP,temp)
      }
    }
  }
}

tol=0.001
BigM=1e+6
data_SGSP_new=data_SGSP%>%group_by(m,n,l,K,method,rep)%>%arrange(time, .by_group = TRUE)%>%mutate(best_obj=cummin(last_obj*(last_slack_const<=tol)+BigM*(last_slack_const>tol)),best_avg_obj=cummin(avg_obj*(avg_slack_const<=tol)+BigM*(avg_slack_const>tol)))
data_opt=rbind(data_OCO%>%group_by(m,n,l,K,rep)%>%summarize(opt_obj=max(lb_obj))%>%ungroup(),
               data_SGSP_new%>%group_by(m,n,l,K,rep)%>%arrange(time, .by_group = TRUE)%>%summarize(opt_obj=min(pmin(best_obj,best_avg_obj)))%>%ungroup())%>%
          group_by(m,n,l,K,rep)%>%summarize(opt_obj=min(opt_obj))%>%ungroup()
data_opt_ratio_line=data_opt%>%mutate(tol=tol/abs(opt_obj))%>%group_by(m,n,l,K)%>%summarize(med_tol=median(tol))
data_OCO_new=left_join(data_OCO,data_opt)
data_OCO_new=data_OCO_new%>%group_by(m,n,l,K,rep,method)%>%arrange(time, .by_group = TRUE)%>%
                            mutate(obj=last_obj,obj_min=pmin(obj,avg_obj),
                                   best_obj=cummin((outer_iter>1)*(last_slack_const>tol)*ub_obj+(last_slack_const<=tol)*pmin(last_obj,ub_obj)+(outer_iter<=1)*(last_slack_const>tol)*BigM),
                                   best_avg_obj=cummin((outer_iter>1)*(avg_slack_const>tol)*ub_obj+(avg_slack_const<=tol)*pmin(avg_obj,ub_obj)+(outer_iter<=1)*(avg_slack_const>tol)*BigM),
                                   best_min_obj=pmin(best_obj,best_avg_obj),
                                   obj_ratio=(obj-opt_obj)/abs(opt_obj),opt_ratio=(best_obj-opt_obj)/abs(opt_obj),feas_gap=last_slack_const,iter=inner_iter,
                                   avg_opt_ratio=(best_avg_obj-opt_obj)/abs(opt_obj),avg_feas_gap=avg_slack_const,
                                   min_opt_ratio=(best_min_obj-opt_obj)/abs(opt_obj),min_feas_gap=pmin(feas_gap,avg_feas_gap) )%>%ungroup()%>%
              select(m,n,l,K,rep,method,time,obj_ratio,opt_ratio,best_obj,feas_gap,avg_opt_ratio,best_avg_obj,avg_feas_gap,min_opt_ratio,best_min_obj,min_feas_gap,iter)
data_SGSP_new=left_join(data_SGSP,data_opt)
data_SGSP_new=data_SGSP_new%>%group_by(m,n,l,K,method,rep)%>%arrange(time, .by_group = TRUE)%>%
              mutate(obj=last_obj,best_obj=cummin(last_obj*(last_slack_const<=tol)+BigM*(last_slack_const>tol)),
                     avg_obj=avg_obj,best_avg_obj=cummin(avg_obj*(avg_slack_const<=tol)+BigM*(avg_slack_const>tol)),
                     obj_min=pmin(obj,avg_obj),best_min_obj=cummin(pmin(best_obj,best_avg_obj)))%>%ungroup()%>%
              mutate(obj_ratio=(obj-opt_obj)/abs(opt_obj),opt_ratio=(best_obj-opt_obj)/abs(opt_obj),feas_gap=last_slack_const,
              avg_opt_ratio=(best_avg_obj-opt_obj)/abs(opt_obj),avg_feas_gap=avg_slack_const,
              min_opt_ratio=(best_min_obj-opt_obj)/abs(opt_obj),min_feas_gap=pmin(feas_gap,avg_feas_gap),
              iter=inner_iter)%>%
              select(m,n,l,K,rep,method,time,obj_ratio,opt_ratio,best_obj,feas_gap,avg_opt_ratio,best_avg_obj,avg_feas_gap,min_opt_ratio,best_min_obj,min_feas_gap,iter)
all_data_new=rbind(data_OCO_new,data_SGSP_new)
if (is_empty(all_data_new%>%filter(method=="SGSP",time==0))){
  temp=all_data_new%>%filter(method=="PDOCOTT",time==0)%>%mutate(method="SGSP")
  all_data_new=rbind(all_data_new,temp)
}


#computes performance statistics for each point in time
time_horizon=max(scenarios$time_horizon);
time_values=c(seq(0,1,0.1),seq(2,3600,1))
for (t in time_values){
  opt_by_time=all_data_new%>%filter(time<=t,best_obj<BigM)%>%group_by(m,n,l,K,method,rep)%>%slice_max(time, n = 1)%>% slice_min(opt_ratio, n = 1)%>%slice_max(iter, n=1)%>%mutate(time_round=t)
  temp=opt_by_time%>%group_by(m,n,l,K,time_round,method)%>%summarize(avg_obj_ratio=median(opt_ratio),upper_obj_ratio=quantile(opt_ratio,.8),lower_obj_ratio=quantile(opt_ratio,.2),
                                                                                    num_feas=n())
  min_opt_by_time=all_data_new%>%filter(time<=t,best_min_obj<BigM)%>%group_by(m,n,l,K,method,rep)%>% slice_max(time, n = 1)%>% slice_min(min_opt_ratio, n = 1)%>%slice_max(iter, n=1)%>%mutate(time_round=t)
  tempa=min_opt_by_time%>%group_by(m,n,l,K,time_round,method)%>%summarize(avg_min_obj_ratio=median(min_opt_ratio),upper_min_obj_ratio=quantile(min_opt_ratio,.95),lower_min_obj_ratio=quantile(min_opt_ratio,.05),
                                                                     num_min_feas=n())
  
  temp=left_join(tempa,temp)
  
  feas_by_time=all_data_new%>%filter(time<=t)%>%group_by(m,n,l,K,method,rep)%>%slice_max(time, n=1)%>% slice_min(feas_gap, n = 1)%>%slice_max(iter, n=1)%>%mutate(time_round=t,feas_gap=pmax(feas_gap,1e-10))%>%select(m,n,l,K,rep,time,method,time_round,feas_gap)%>%ungroup()
  temp1=feas_by_time%>%group_by(m,n,l,K,time_round,method)%>%summarize(avg_feas_gap=median(feas_gap),upper_feas_gap=quantile(feas_gap,.8),lower_feas_gap=quantile(feas_gap,.2))
  
  feas_by_time2=all_data_new%>%filter(time<=t)%>%group_by(m,n,l,K,method,rep)%>%mutate(feas_gap=cummin(feas_gap))%>%slice_max(time, n=1)%>% slice_min(feas_gap, n = 1)%>%slice_max(iter, n=1)%>%filter(best_obj>=1e+6)%>%mutate(time_round=t,feas_gap=max(feas_gap,1e-10))%>%select(m,n,l,K,rep,method,time,time_round,feas_gap,obj_ratio)%>%ungroup()
  temp2=feas_by_time2%>%group_by(m,n,l,K,time_round,method)%>%summarize(avg_feas_gap_not_feas=median(feas_gap),upper_feas_gap_not_feas=quantile(feas_gap,1),"9_feas_gap_not_feas"=quantile(feas_gap,.9),"8_feas_gap_not_feas"=quantile(feas_gap,.8),"7_feas_gap_not_feas"=quantile(feas_gap,.7),"6_feas_gap_not_feas"=quantile(feas_gap,.6),
                                                                        "4_feas_gap_not_feas"=quantile(feas_gap,.4),"3_feas_gap_not_feas"=quantile(feas_gap,.3),"2_feas_gap_not_feas"=quantile(feas_gap,.2),"1_feas_gap_not_feas"=quantile(feas_gap,.1),lower_feas_gap_not_feas=quantile(feas_gap,0),
                                                                        avg_obj_gap_not_feas=median(obj_ratio),upper_obj_gap_not_feas=quantile(obj_ratio,1),"9_obj_gap_not_feas"=quantile(obj_ratio,.9),"8_obj_gap_not_feas"=quantile(obj_ratio,.8),"7_obj_gap_not_feas"=quantile(obj_ratio,.7),"6_obj_gap_not_feas"=quantile(obj_ratio,.6),
                                                                        "4_obj_gap_not_feas"=quantile(obj_ratio,.4),
                                                                        "3_obj_gap_not_feas"=quantile(obj_ratio,.3),"2_obj_gap_not_feas"=quantile(obj_ratio,.2),"1_obj_gap_not_feas"=quantile(obj_ratio,.1),lower_obj_gap_not_feas=quantile(obj_ratio,0))
  feas_by_time2a=all_data_new%>%filter(time<=t)%>%group_by(m,n,l,K,method,rep)%>%mutate(min_feas_gap=cummin(min_feas_gap))%>%slice_max(time, n=1)%>% slice_min(min_feas_gap, n = 1)%>%slice_max(iter, n=1)%>%filter(best_min_obj>=1e+6)%>%mutate(time_round=t,min_feas_gap=max(min_feas_gap,1e-10))%>%select(m,n,l,K,rep,method,time,time_round,min_feas_gap)%>%ungroup()
  temp2a=feas_by_time2a%>%group_by(m,n,l,K,time_round,method)%>%summarize(avg_min_feas_gap_not_feas=median(min_feas_gap),upper_min_feas_gap_not_feas=quantile(min_feas_gap,1),"9_min_feas_gap_not_feas"=quantile(min_feas_gap,.9),"8_min_feas_gap_not_feas"=quantile(min_feas_gap,.8),"7_min_feas_gap_not_feas"=quantile(min_feas_gap,.7),"6_min_feas_gap_not_feas"=quantile(min_feas_gap,.6),
                                                                        "4_min_feas_gap_not_feas"=quantile(min_feas_gap,.4),"3_min_feas_gap_not_feas"=quantile(min_feas_gap,.3),"2_min_feas_gap_not_feas"=quantile(min_feas_gap,.2),"1_min_feas_gap_not_feas"=quantile(min_feas_gap,.2),lower_min_feas_gap_not_feas=quantile(min_feas_gap,0))#,
                                                                       
  temp2=left_join(temp2a,temp2)
  
  
  temp3=left_join(left_join(temp1,temp),temp2)
  if (t==0){
    sum_epsilon_by_time=temp3
  }else{
    sum_epsilon_by_time=rbind(sum_epsilon_by_time,temp3)
  }  
}
sum_epsilon_by_time$num_feas[is.na(sum_epsilon_by_time$num_feas)]=0
sum_epsilon_by_time$num_min_feas[is.na(sum_epsilon_by_time$num_min_feas)]=0



methods_names=c("Cut. Planes","FO-Pess.","OFO","SGSP")
methods_colors=c("blue","red","cyan","green")

#########################
#show performance over time - Figure 1
for (row in seq(1,nrow(scenarios))){
  if (scenarios$m[row]==1){
    rel_data_temp=sum_epsilon_by_time%>%filter(m==1,n==scenarios$n[row],time_round<=scenarios$time_horizon[row])
    if (row==1){
      rel_data=rel_data_temp
    } else {
      rel_data=rbind(rel_data,rel_data_temp)
     }
       
  }
}
p1=ggplot(rel_data)+geom_line(aes(x=time_round,y=pmax(avg_obj_ratio,1e-4),color=method))+
  geom_ribbon(aes(x=time_round,ymin=pmax(lower_obj_ratio,1e-4),ymax=pmax(upper_obj_ratio,1e-4),fill=method),alpha=.2)+
  scale_y_log10()+
  scale_x_log10()+
  theme_bw(11)+ 
  theme(legend.title = element_blank(),
        legend.position = "bottom",
        strip.background = element_rect(
          color="black", fill="white", size=1.5, linetype="solid"
        ),
        strip.text = element_text(
          size = 12, color = "black", face = "bold"
        ))+
  scale_color_manual(values=methods_colors,labels=methods_names)+
  scale_fill_manual(values=methods_colors,labels=methods_names)+
  geom_hline(data=data_opt_ratio_line%>%filter(m==1),aes(yintercept=med_tol),linetype="dashed",color="black")+
  #facet_grid(.~n,scales="free_x")+
  facet_grid(.~n,scales="free_x",,labeller=label_both)+
  annotation_logticks(sides="lb") +
  ylab("Optimality gap ratio")+
  xlab("Time (Sec)")
file_name=sprintf("opt_gap_noconstraints.eps")
ggsave(p1, file=file_name, device="eps",width = 8, height = 4)
ggsave(p1, file=paste(file_name,".png",sep=""),width = 8, height = 4)  

#Figure 2-4
#Only for problems with constraints
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

for (sc in seq(1,nrow(scenarios))){
    if (scenarios$m[sc]>1){
        m1=scenarios$m[sc]
        n1=scenarios$n[sc]
        k1=scenarios$k[sc]
        l1=scenarios$l[sc]
        time_horizon=scenarios$time_horizon[sc]
        rel_data=sum_epsilon_by_time%>%filter(m==scenarios$m[sc],n==scenarios$n[sc],time_round<=scenarios$time_horizon[sc])
        
        p1=ggplot(rel_data)+geom_line(aes(x=time_round,y=pmax(avg_obj_ratio,1e-4),color=method))+
          geom_ribbon(aes(x=time_round,ymin=pmax(lower_obj_ratio,1e-4),ymax=pmax(upper_obj_ratio,1e-4),fill=method),alpha=.2)+
          scale_y_log10()+
          scale_x_log10(limit=c(min_time,time_horizon))+
          theme_bw(11)+ 
          theme(legend.title = element_blank(),
                strip.background = element_rect(
                  color="black", fill="white", size=1.5, linetype="solid"
                ),
                strip.text = element_text(
                  size = 12, color = "black", face = "bold"
                ))+
          scale_color_manual(values=methods_colors,labels=methods_names)+
          scale_fill_manual(values=methods_colors,labels=methods_names)+
          geom_hline(data=data_opt_ratio_line%>%filter(m==m1,n==n1),aes(yintercept=med_tol),linetype="dashed",color="black")+
          #facet_grid(K~n+l,scales="free_y")+
          #facet_grid(K~n,scales="free_y",,labeller=label_both)+
          annotation_logticks(sides="lb") +
          ylab("Optimality gap ratio")+
          xlab("Time (Sec)")
        file_name=sprintf("opt_gap_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
        ggsave(p1, file=file_name, device="eps",width = 6, height = 4)
        ggsave(p1, file=paste(file_name,".png",sep=""),width = 6, height = 4)  
        
        p3=ggplot(rel_data)+geom_line(aes(x=time_round,y=(num_feas/num_instance)*100,color=method))+
          theme_bw(11)+ 
          scale_x_log10(limit=c(min_time,time_horizon))+
          scale_y_continuous(limit=c(0,100))+
          #facet_grid(K~n+l,scales="free_y")+
          #facet_grid(K~n,scales="free_y",labeller=label_both)+
          scale_color_manual(values=methods_colors,labels=methods_names)+
          ylab("% Feasible Realizations")+
          xlab("Time (Sec)")+
          annotation_logticks(sides="b") +
          theme(legend.title = element_blank(),
                strip.background = element_rect(
                  color="black", fill="white", size=1.5, linetype="solid"
                ),
                strip.text = element_text(
                  size = 12, color = "black", face = "bold"
                ))
        file_name=sprintf("feas_realizations_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
        ggsave(p3, file=file_name, device="eps")
        ggsave(p3, file=paste(file_name,".png",sep=""))
        
        
        
        if(sum(is.na(rel_data$avg_feas_gap_not_feas))<length(rel_data$avg_feas_gap_not_feas)){
          p4=ggplot(rel_data)+geom_line(aes(x=time_round,y=avg_feas_gap_not_feas,color=method))+
            #geom_ribbon(aes(x=time_round,ymin=lower_feas_gap_not_feas,ymax=upper_feas_gap_not_feas,fill=method),alpha=.1)+
            geom_ribbon(aes(x=time_round,ymin=`1_feas_gap_not_feas`,ymax=`9_feas_gap_not_feas`,fill=method),alpha=.1)+
            #geom_ribbon(aes(x=time_round,ymin=`2_feas_gap_not_feas`,ymax=`8_feas_gap_not_feas`,fill=method),alpha=.1)+
            #geom_ribbon(aes(x=time_round,ymin=`3_feas_gap_not_feas`,ymax=`7_feas_gap_not_feas`,fill=method),alpha=.1)+
            #geom_ribbon(aes(x=time_round,ymin=`4_feas_gap_not_feas`,ymax=`6_feas_gap_not_feas`,fill=method),alpha=.1)+
            scale_y_log10()+
            scale_x_log10(limit=c(min_time,time_horizon))+
            #scale_color_manual(values=methods_colors,labels=methods_names)+
            #scale_color_manual(name = "Gaps",values=c("green","blue"),labels=c("Feas gap","Opt ratio"))+
            #scale_fill_manual(name = "Gaps",values=c("green","blue"),labels=c("Feas gap","Opt ratio"))+
            scale_color_manual(values=methods_colors,labels=methods_names)+
            scale_fill_manual(values=methods_colors,labels=methods_names)+
            theme_bw(11)+ 
            theme(legend.title = element_blank(),
                  strip.background = element_rect(
                    color="black", fill="white", size=1.5, linetype="solid"
                  ),
                  strip.text = element_text(
                    size = 12, color = "black", face = "bold"
                  ))+
            ylab("Feasibility gap")+
            xlab("Time (Sec)")+
            annotation_logticks(sides="lb") +
            geom_hline(data=data_opt_ratio_line%>%filter(m==m1,n==n1),aes(yintercept=0.001),linetype="dashed",color="black")
          #facet_grid(K~n+l,scales="free_y")+
          #facet_grid(K~n,labeller=label_both)#,scales="free_y")
          file_name=sprintf("infeasible_feas_gap_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
          ggsave(p4, file=file_name, device="eps")
          ggsave(p4, file=paste(file_name,".png",sep=""))
          print(file_name)
          
          legend <- g_legend(p1 + theme(legend.position="bottom"))
          file_name=sprintf("all_graphs_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
          p5=grid.arrange(p3+ theme(legend.position = "none") , p4+ theme(legend.position = "none"), p1 + theme(legend.position="none"), legend, layout_matrix=rbind(c(1,2,3), c(4,4,4)),heights=c(6,1))
          ggsave(p5, file=file_name, device="eps",width = 8, height = 3)
          ggsave(p5, file=paste(file_name,".png",sep=""),width = 8, height = 3)
        }  
        
        # p1a=ggplot(rel_data)+geom_line(aes(x=time_round,y=pmax(avg_min_obj_ratio,1e-4),color=method))+
        # geom_ribbon(aes(x=time_round,ymin=pmax(lower_min_obj_ratio,1e-4),ymax=pmax(upper_min_obj_ratio,1e-4),fill=method),alpha=.2)+
        # scale_y_log10()+
        # scale_x_log10(limit=c(min_time,time_horizon))+
        # theme_bw(11)+ 
        # theme(legend.title = element_blank(),          
        #   strip.background = element_rect(
        #     color="black", fill="white", size=1.5, linetype="solid"
        #   ),
        #   strip.text = element_text(
        #     size = 12, color = "black", face = "bold"
        #   ))+
        # scale_color_manual(values=methods_colors,labels=methods_names)+
        # scale_fill_manual(values=methods_colors,labels=methods_names)+
        # geom_hline(data=data_opt_ratio_line%>%filter(m==m1,n==n1),aes(yintercept=med_tol),linetype="dashed",color="black")+
        # annotation_logticks(sides="lb") +
        # ylab("Optimality gap")+
        # xlab("Time (Sec)")
        # file_name=sprintf("min_opt_gap_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
        # ggsave(p1a, file=file_name, device="eps",width = 6, height = 4)
        # ggsave(p1a, file=paste(file_name,".png",sep=""),width = 6, height = 4)
        # 
        # p3a=ggplot(rel_data)+geom_line(aes(x=time_round,y=(num_min_feas/num_instance)*100,color=method))+
        #   theme_bw(11)+ 
        #   scale_x_log10(limit=c(0.1,time_horizon))+
        #   scale_y_continuous(limit=c(0,100))+
        #   #facet_grid(K~n+l,scales="free_y")+
        #   #facet_grid(K~n,scales="free_y",labeller=label_both)+
        #   scale_color_manual(values=methods_colors,labels=methods_names)+
        #   ylab("% Feasible Realizations")+
        #   xlab("Time (Sec)")+
        #   annotation_logticks(sides="b") +
        #   theme(legend.title = element_blank(),
        #         strip.background = element_rect(
        #           color="black", fill="white", size=1.5, linetype="solid"
        #         ),
        #         strip.text = element_text(
        #           size = 12, color = "black", face = "bold"
        #         ))
        # file_name=sprintf("min_feas_realizations_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
        # ggsave(p3a, file=file_name, device="eps")
        # ggsave(p3a, file=paste(file_name,".png",sep=""))
        # 
        # if(sum(is.na(rel_data$avg_min_feas_gap_not_feas))<length(rel_data$avg_min_feas_gap_not_feas)){
        #     p4a=ggplot(rel_data)+geom_line(aes(x=time_round,y=avg_min_feas_gap_not_feas,color=method))+
        #     geom_ribbon(aes(x=time_round,ymin=lower_feas_gap_not_feas,ymax=upper_feas_gap_not_feas,fill=method),alpha=.1)+
        #     #geom_ribbon(aes(x=time_round,ymin=`1_min_feas_gap_not_feas`,ymax=`9_min_feas_gap_not_feas`,fill=method),alpha=.1)+
        #     #geom_ribbon(aes(x=time_round,ymin=`2_feas_gap_not_feas`,ymax=`8_feas_gap_not_feas`,fill=method),alpha=.1)+
        #     #geom_ribbon(aes(x=time_round,ymin=`3_feas_gap_not_feas`,ymax=`7_feas_gap_not_feas`,fill=method),alpha=.1)+
        #     #geom_ribbon(aes(x=time_round,ymin=`4_feas_gap_not_feas`,ymax=`6_feas_gap_not_feas`,fill=method),alpha=.1)+
        #     scale_y_log10()+
        #     scale_x_log10(limit=c(.1,time_horizon))+
        #     #scale_color_manual(values=methods_colors,labels=methods_names)+
        #     #scale_color_manual(name = "Gaps",values=c("green","blue"),labels=c("Feas gap","Opt ratio"))+
        #     #scale_fill_manual(name = "Gaps",values=c("green","blue"),labels=c("Feas gap","Opt ratio"))+
        #     scale_color_manual(values=methods_colors,labels=methods_names)+
        #     scale_fill_manual(values=methods_colors,labels=methods_names)+
        #     theme_bw(11)+ 
        #     theme(legend.title = element_blank(),
        #       strip.background = element_rect(
        #         color="black", fill="white", size=1.5, linetype="solid"
        #       ),
        #       strip.text = element_text(
        #         size = 12, color = "black", face = "bold"
        #       ))+
        #       ylab("Feasibility gap")+
        #       xlab("Time (Sec)")+
        #       annotation_logticks(sides="lb") 
        #     #facet_grid(K~n+l,scales="free_y")+
        #     #facet_grid(K~n,labeller=label_both)#,scales="free_y")
        #     file_name=sprintf("min_infeasible_feas_gap_n%d_m%d_k%d_l%d.eps",n1,m1,k1,l1)
        #     ggsave(p4a, file=file_name, device="eps")
        #     ggsave(p4a, file=paste(file_name,".png",sep=""))
        #     print(file_name)
  
            # legend <- g_legend(p1a + theme(legend.position="bottom"))
            # file_name=sprintf("min_all_graphs_n%d_m%d.eps",n1,m1)
            # p5a=grid.arrange(p3a+ theme(legend.position = "none") , p4a+ theme(legend.position = "none"), p1a + theme(legend.position="none"), legend, layout_matrix=rbind(c(1,2,3), c(4,4,4)),heights=c(6,1))
            # ggsave(p5a, file=file_name, device="eps",width = 8, height = 3)
            # ggsave(p5a, file=paste(file_name,".png",sep=""),width = 8, height = 3)
        }
    }
}
