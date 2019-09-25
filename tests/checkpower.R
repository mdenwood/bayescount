library('bayescount')

if(FALSE){

out <- bayescount:::fecrt_power_wrap(reduction=0.95)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)

out <- bayescount:::fecrt_power_wrap(reduction=0.89)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)

out <- bayescount:::fecrt_power_wrap(reduction=0.89, paired_data=TRUE, paired_analysis=TRUE)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)
out <- bayescount:::fecrt_power_wrap(reduction=0.89, paired_data=TRUE, paired_analysis=FALSE)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)

out <- bayescount:::fecrt_power_wrap(use_truek=TRUE, preN=10, paired_data=TRUE, paired_analysis=FALSE, approximation=1, reduction=0.951, animalk=1, efficacyk=1, prek=1, postk=0.7)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)
out <- bayescount:::fecrt_power_wrap(use_truek=TRUE, preN=10, paired_data=TRUE, paired_analysis=TRUE, approximation=1, reduction=0.90, animalk=1, efficacyk=Inf, prek=1, postk=0.7)
sum(out$p_1 < 0.025) / nrow(out)
sum(out$p_2 < 0.025) / nrow(out)


out <- bayescount:::fecrt_power_comparison_wrap(reduction=0.95)


#stop('At a minimum I need to compare fec.power to nb dist')

stop('Also need to change geometric inits')

}
