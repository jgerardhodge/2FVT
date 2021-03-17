#Note: This script serves as the final step in the 2FVT analysis workflow and largely serves as a tool for
#providing image outputs of the parameters fitted during analysis by 'plm_modeling_AIC_2FVT.R'.  As such the 
#assorted '..._2FVT_fitted_parameters.csv', '..._whole_plant_FVTs.csv', and 'All_2FVT_fitted_parameters.csv'
#serve as required input files prior to running. 

path='/Users/johnhodge/Desktop/AIC_outputs/2FVT_fits_w_Wald_CIs/'
lines=c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100');

line_cols=c(rgb(0,0.39,0,0.9), rgb(0.63, 0.13, 0.94,0.9), rgb(1, 0.65, 0, 0.9), rgb(0,0,1,0.9), rgb(0.8,0.6,0.11,0.9))
CI_cols=c(rgb(0,0.39,0,0.3), rgb(0.63, 0.13, 0.94,0.3), rgb(1, 0.65, 0, 0.3), rgb(0,0,1,0.3), rgb(0.8,0.6,0.11,0.3))

N<-cbind(c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100'), c(10, 9, 9, 8, 10))

method_p1<-'A10'
method_p2<-'B100'

x_fit=seq(4,60,by=0.1)


sig_lvl=c(1.96,2.56, 3.29)


sig_palette=cbind(c(rev(-sig_lvl),0,sig_lvl), c(colorRampPalette(c('deepskyblue','black'))(n=length(sig_lvl)),'black',colorRampPalette(c('black','gold'))(n=length(sig_lvl))))
sig_palette<-as.data.frame(sig_palette); sig_palette[,1]<-as.numeric(as.character(sig_palette[,1])); sig_palette[,2]<-as.character(sig_palette[,2]);


if(!is.null(method_p1) & !is.null(method_p2)){	
	method_title='Transgression Test'
}else if(!is.null(method_p1)){
	method_title='Standard Test'
}else{
	method_title='Outlier Test'
}

whole_plant_attr=c()

for (line in c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')){
	whole_plant_fits<-read.csv(file=paste('~/Desktop/AIC_outputs/2FVT_parameters/', line, '_whole_plant_FVTs.csv', sep=''), row.names=1)
	cur_plant_attr<-c(whole_plant_fits[1,],whole_plant_fits[2,])
	names(cur_plant_attr)<-c('flw_xo', 'flw_k', 'flw_L', 'lfno_xo', 'lfno_k', 'lfno_L')
	whole_plant_attr<-rbind(whole_plant_attr, cur_plant_attr)
}

rownames(whole_plant_attr)=c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')

meta_attr<-read.csv(file='~/Desktop/AIC_outputs/2FVT_parameters/All_2FVT_fitted_parameters.csv')

obs_culm_hgt<-read.csv(file='~/Desktop/AIC_outputs/QC_hgt_flw.csv')

#################################################################
#1A. Leaf number with confidence intervals by genotype
#################################################################

#Merge fitted GLM data into a single list to enable calling of measures between genotypes

fit_list<-list()

iter = 1

for (line in 1:length(lines)){
	L    =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6])
	k    =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],5])
	xo   =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],4])
	y_fit        =L/(1+exp(-k*(x_fit-xo)))

	y_p          =y_fit/max(y_fit);
	p            =(y_p+2)/(length(y_p)+4)
	CI_weight    =max(y_fit)*sqrt((p*(1-p))/length(y_fit)); 


	cur_data=cbind(x_fit, y_fit, CI_weight)

	fit_list[[iter]]<-cur_data	
	print(iter)
	iter=iter+1
}

if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
	L<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,6]))
	k<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,5]))
	xo<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,4]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))
}else if(!is.null(method_p1)){											#>>> Standard test performed
	L<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,6]))
	k<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,5]))
	xo<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,4]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))
}else{																	#>>> Outlier test performed
	L<-mean(as.numeric(whole_plant_attr[,6]))
	k<-mean(as.numeric(whole_plant_attr[,5]))
	xo<-mean(as.numeric(whole_plant_attr[,4]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))																	
}

method_fits<-cbind(x_fit, test_method)
fit_list[[length(lines)+1]]<-method_fits

lfno_list<-fit_list

##############################
# Segmented Heights (Gaussian)
##############################


fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/AIC_outputs/2FVT_parameters/', lines[line], '_seg_culm_hgt_2FVT_fitted_parameters.csv', sep=''))
	cur_n=as.numeric(N[line,2])
	
	geno_list=c()

	geno_ceiling=floor(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6]))
	geno_ceiling=geno_ceiling+1							#Add one unit to account for peduncle

	for (phyt in 2:geno_ceiling){
		Lg   =cur_pars[1,14]
		kg   =cur_pars[1,13]
		xg   =cur_pars[1,12]

		L    =Lg*exp(-(phyt-xg)^2/(2*kg^2))
		k    =as.numeric(cur_pars[1,5])
		xo   =(phyt-as.numeric(cur_pars[1,8]))/as.numeric(cur_pars[1,9])

		y_fit        =L/(1+exp(-k*(x_fit-xo)))
		y_p          =y_fit/max(y_fit);
		p            =(y_p+2)/(length(y_p)+4)
		CI_weight    =max(y_fit)*sqrt((p*(1-p))/length(y_fit)); 

		cur_fits=cbind(y_fit, CI_weight)
		colnames(cur_fits)=c(paste('y_fit',phyt, sep=''), paste('y_Wald',phyt, sep=''))

		geno_list=cbind(geno_list, cur_fits)

	}

	fit_list[[iter]]<-cbind(x_fit, geno_list)
	iter=iter+1

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
		if(lines[line]==method_p1 | lines[line]==method_p2){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else if(!is.null(method_p1)){											#>>> Standard test performed
		if(lines[line]==method_p1){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else{		
		method_pars=rbind(method_pars,cur_pars[1,])
	}

}

method_fits=c()

method_ceiling=floor(max(as.numeric(whole_plant_attr[,6])))
method_ceiling=method_ceiling+1							#Add one unit to account for peduncle

for (phyt in 2:method_ceiling){

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
		Lg   =mean(as.numeric(method_pars[,14]))
		kg   =mean(as.numeric(method_pars[,13]))
		xg   =mean(as.numeric(method_pars[,12]))

		L    =Lg*exp(-(phyt-xg)^2/(2*kg^2))
		k    =mean(as.numeric(method_pars[,5]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else if(!is.null(method_p1)){											#>>> Standard test performed
		Lg   =as.numeric(method_pars[,14])
		kg   =as.numeric(method_pars[,13])
		xg   =as.numeric(method_pars[,12])

		L    =Lg*exp(-(phyt-xg)^2/(2*kg^2))
		k    =as.numeric(method_pars[,5])
		xo   =(phyt-as.numeric(method_pars[,8]))/as.numeric(method_pars[,9])

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else{																	#>>> Outlier test performed
		Lg   =mean(as.numeric(method_pars[,14]))
		kg   =mean(as.numeric(method_pars[,13]))
		xg   =mean(as.numeric(method_pars[,12]))

		L    =Lg*exp(-(phyt-xg)^2/(2*kg^2))
		k    =mean(as.numeric(method_pars[,5]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))																
	}

	method_fits=cbind(method_fits, test_method)
	colnames(method_fits)[ncol(method_fits)]=paste('standard',phyt,sep='')

	for (iter in 1:length(fit_list)){

		y_fit=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_fit', phyt, sep='')]
		y_wald=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_Wald', phyt, sep='')]
		
	}

}

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

sghgt_list<-fit_list

#Culm Height Tests by Genotype against 2FVT Model

test_x=c(19, 24, 23, 30, 27)
hgt_test=c()

for (line in 1:length(lines)){
	
	culm_hgt=c()

	for (cur_x in 1:length(x_fit)){

		geno_ceiling=floor(lfno_list[[line]][cur_x,2])
																		#Add one unit to account for peduncle
		sum_y_fit<-0
		sum_CI_weight<-0

		for (phyt in 2:geno_ceiling){

			y_fit       =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_fit', phyt, sep='')]
			CI_weight   =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_Wald', phyt, sep='')]

			sum_y_fit=sum_y_fit+y_fit
			sum_CI_weight=sum_CI_weight+CI_weight

			
		}

		culm_hgt<-c(culm_hgt, sum_y_fit)

	}

	test_y<-culm_hgt[x_fit==test_x[line]]

	obs_y<-obs_culm_hgt[obs_culm_hgt[,1]==lines[line],5]
	mean_obs_y<-mean(obs_y)

	delta_y<-mean_obs_y-test_y
	pval_y<-t.test(obs_y, alternative='two.sided', mu=test_y)$p.value

	hgt_test<-rbind(hgt_test, cbind(lines[line], 'gaussian', test_y, mean_obs_y, delta_y, pval_y))
}

##############################
# Segmented Heights (Linear)
##############################


fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/AIC_outputs/2FVT_parameters/', lines[line], '_seg_culm_hgt_2FVT_fitted_parameters.csv', sep=''))
	cur_n=as.numeric(N[line,2])
	
	geno_list=c()

	geno_ceiling=floor(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6]))
	geno_ceiling=geno_ceiling+1							#Add one unit to account for peduncle

	for (phyt in 2:geno_ceiling){

		L    =(phyt-as.numeric(cur_pars[1,10]))/as.numeric(cur_pars[1,11])
		k    =(phyt-as.numeric(cur_pars[1,6]))/as.numeric(cur_pars[1,7])
		xo   =(phyt-as.numeric(cur_pars[1,8]))/as.numeric(cur_pars[1,9])

		y_fit        =L/(1+exp(-k*(x_fit-xo)))
		y_p          =y_fit/max(y_fit);
		p            =(y_p+2)/(length(y_p)+4)
		CI_weight    =max(y_fit)*sqrt((p*(1-p))/length(y_fit)); 

		cur_fits=cbind(y_fit, CI_weight)
		colnames(cur_fits)=c(paste('y_fit',phyt, sep=''), paste('y_Wald',phyt, sep=''))

		geno_list=cbind(geno_list, cur_fits)

	}

	fit_list[[iter]]<-cbind(x_fit, geno_list)
	iter=iter+1

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
		if(lines[line]==method_p1 | lines[line]==method_p2){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else if(!is.null(method_p1)){											#>>> Standard test performed
		if(lines[line]==method_p1){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else{		
		method_pars=rbind(method_pars,cur_pars[1,])
	}

}

method_fits=c()

method_ceiling=floor(max(as.numeric(whole_plant_attr[,6])))
method_ceiling=method_ceiling+1							#Add one unit to account for peduncle

for (phyt in 2:method_ceiling){

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed


		L    =(phyt-mean(as.numeric(method_pars[,10])))/mean(as.numeric(method_pars[,11]))
		k    =(phyt-mean(as.numeric(method_pars[,6])))/mean(as.numeric(method_pars[,7]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else if(!is.null(method_p1)){											#>>> Standard test performed

		L    =(phyt-as.numeric(method_pars[,10]))/as.numeric(method_pars[,11])
		k    =(phyt-as.numeric(method_pars[,6]))/as.numeric(method_pars[,7])
		xo   =(phyt-as.numeric(method_pars[,8]))/as.numeric(method_pars[,9])

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else{																	#>>> Outlier test performed

		L    =(phyt-mean(as.numeric(method_pars[,10])))/mean(as.numeric(method_pars[,11]))
		k    =(phyt-mean(as.numeric(method_pars[,6])))/mean(as.numeric(method_pars[,7]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))																
	}

	method_fits=cbind(method_fits, test_method)
	colnames(method_fits)[ncol(method_fits)]=paste('standard',phyt,sep='')

	for (iter in 1:length(fit_list)){

		y_fit=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_fit', phyt, sep='')]
		y_wald=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_Wald', phyt, sep='')]
		
	}

}

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

sghgt_list<-fit_list

#Culm Height Tests by Genotype against 2FVT Model

test_x=c(19, 24, 23, 30, 27)

for (line in 1:length(lines)){
	
	culm_hgt=c()

	for (cur_x in 1:length(x_fit)){

		geno_ceiling=floor(lfno_list[[line]][cur_x,2])
																		#Add one unit to account for peduncle
		sum_y_fit<-0
		sum_CI_weight<-0

		for (phyt in 2:geno_ceiling){

			y_fit       =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_fit', phyt, sep='')]
			CI_weight   =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_Wald', phyt, sep='')]

			sum_y_fit=sum_y_fit+y_fit
			sum_CI_weight=sum_CI_weight+CI_weight

			
		}

		culm_hgt<-c(culm_hgt, sum_y_fit)

	}

	test_y<-culm_hgt[x_fit==test_x[line]]

	obs_y<-obs_culm_hgt[obs_culm_hgt[,1]==lines[line],5]
	mean_obs_y<-mean(obs_y)

	delta_y<-mean_obs_y-test_y
	pval_y<-t.test(obs_y, alternative='two.sided', mu=test_y)$p.value

	hgt_test<-rbind(hgt_test, cbind(lines[line], 'linear', test_y, mean_obs_y, delta_y, pval_y))
}

##############################
# Segmented Heights (Simplified k)
##############################


fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/AIC_outputs/2FVT_parameters/', lines[line], '_seg_culm_hgt_2FVT_fitted_parameters.csv', sep=''))
	cur_n=as.numeric(N[line,2])
	
	geno_list=c()

	geno_ceiling=floor(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6]))
	geno_ceiling=geno_ceiling+1							#Add one unit to account for peduncle

	for (phyt in 2:geno_ceiling){

		L    =(phyt-as.numeric(cur_pars[1,10]))/as.numeric(cur_pars[1,11])
		k    =as.numeric(cur_pars[1,5])
		xo   =(phyt-as.numeric(cur_pars[1,8]))/as.numeric(cur_pars[1,9])

		y_fit        =L/(1+exp(-k*(x_fit-xo)))
		y_p          =y_fit/max(y_fit);
		p            =(y_p+2)/(length(y_p)+4)
		CI_weight    =max(y_fit)*sqrt((p*(1-p))/length(y_fit)); 

		cur_fits=cbind(y_fit, CI_weight)
		colnames(cur_fits)=c(paste('y_fit',phyt, sep=''), paste('y_Wald',phyt, sep=''))

		geno_list=cbind(geno_list, cur_fits)

	}

	fit_list[[iter]]<-cbind(x_fit, geno_list)
	iter=iter+1

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
		if(lines[line]==method_p1 | lines[line]==method_p2){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else if(!is.null(method_p1)){											#>>> Standard test performed
		if(lines[line]==method_p1){
			method_pars=rbind(method_pars,cur_pars[1,])
		}
	}else{		
		method_pars=rbind(method_pars,cur_pars[1,])
	}

}

method_fits=c()

method_ceiling=floor(max(as.numeric(whole_plant_attr[,6])))
method_ceiling=method_ceiling+1							#Add one unit to account for peduncle

for (phyt in 2:method_ceiling){

	if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed


		L    =(phyt-mean(as.numeric(method_pars[,10])))/mean(as.numeric(method_pars[,11]))
		k    =mean(as.numeric(method_pars[,5]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else if(!is.null(method_p1)){											#>>> Standard test performed

		L    =(phyt-as.numeric(method_pars[,10]))/as.numeric(method_pars[,11])
		k    =as.numeric(method_pars[,5])
		xo   =(phyt-as.numeric(method_pars[,8]))/as.numeric(method_pars[,9])

		test_method<-L/(1+exp(-k*(x_fit-xo)))
	}else{																	#>>> Outlier test performed

		L    =(phyt-mean(as.numeric(method_pars[,10])))/mean(as.numeric(method_pars[,11]))
		k    =mean(as.numeric(method_pars[,5]))
		xo   =(phyt-mean(as.numeric(method_pars[,8])))/mean(as.numeric(method_pars[,9]))

		test_method<-L/(1+exp(-k*(x_fit-xo)))																
	}

	method_fits=cbind(method_fits, test_method)
	colnames(method_fits)[ncol(method_fits)]=paste('standard',phyt,sep='')

	for (iter in 1:length(fit_list)){

		y_fit=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_fit', phyt, sep='')]
		y_wald=fit_list[[iter]][,colnames(fit_list[[iter]])==paste('y_Wald', phyt, sep='')]
		
	}

}

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

sghgt_list<-fit_list

#Culm Height Tests by Genotype against 2FVT Model

test_x=c(19, 24, 23, 30, 27)

for (line in 1:length(lines)){
	
	culm_hgt=c()

	for (cur_x in 1:length(x_fit)){

		geno_ceiling=floor(lfno_list[[line]][cur_x,2])
																		#Add one unit to account for peduncle
		sum_y_fit<-0
		sum_CI_weight<-0

		for (phyt in 2:geno_ceiling){

			y_fit       =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_fit', phyt, sep='')]
			CI_weight   =sghgt_list[[line]][cur_x,colnames(sghgt_list[[line]])==paste('y_Wald', phyt, sep='')]

			sum_y_fit=sum_y_fit+y_fit
			sum_CI_weight=sum_CI_weight+CI_weight

			
		}

		culm_hgt<-c(culm_hgt, sum_y_fit)

	}

	test_y<-culm_hgt[x_fit==test_x[line]]

	obs_y<-obs_culm_hgt[obs_culm_hgt[,1]==lines[line],5]
	mean_obs_y<-mean(obs_y)

	delta_y<-mean_obs_y-test_y
	pval_y<-t.test(obs_y, alternative='two.sided', mu=test_y)$p.value

	hgt_test<-rbind(hgt_test, cbind(lines[line], 'simple', test_y, mean_obs_y, delta_y, pval_y))
}

colnames(hgt_test)<-c('Genotype', 'Method', '2FVT_fitted_hgt', 'Avg_obs_hgt', 'Delta', 'T.test_pvalue')

hgt_test<-as.data.frame(hgt_test)
hgt_test

pdf(file='~/Desktop/AIC_outputs/2FVT_model_graphs/2FVT_culm_height_tests.pdf', height=5, width=6)

plot(c(-0.8, length(lines)), c(-25, max(as.numeric(as.character(hgt_test[,3])))), main='Genotype Height Treatments', xlab='', ylab='Culm Height (cm)', col='white', xaxt='n')
abline(v=0:6, col='gray70', lwd=3); abline(h=0, col='gray70', lwd=3); axis(1, at=1:length(lines)-0.5, labels=lines); 
legend('topleft', pch=c(19, 22, 24, 1), legend=c('Mean Observed Heights', '2FVT (Gaussian)', '2FVT (Complete Linear)', '2FVT (Simplified k)'), cex=0.7, bg='white')


for (line in 1:length(lines)){
	cur_group<-hgt_test[hgt_test[,1]==lines[line],]
	points(line-0.8, as.numeric(as.character(cur_group[1,4])), pch=19, col=line_cols[line])
	points(line-0.6, as.numeric(as.character(cur_group[1,3])), pch=22, col=line_cols[line])
	points(line-0.4, as.numeric(as.character(cur_group[2,3])), pch=24, col=line_cols[line])
	points(line-0.2, as.numeric(as.character(cur_group[3,3])), pch=1, col=line_cols[line])

	pvals<-as.numeric(as.character(cur_group[,6]))

	pval_stars<-c()

	for (p in pvals){
		if(p>0.05){
			pval_stars<-c(pval_stars, 'N.S.')
		}else if(p<0.05 & p>0.01){
			pval_stars<-c(pval_stars, '*')
		}else if(p<0.01 & p>0.001){
			pval_stars<-c(pval_stars, '**')
		}else if(p<0.001){
			pval_stars<-c(pval_stars, '***')
		}
	}

	arrows(line-0.7, -5, line-0.8, -5, length=0.025, angle=90); arrows(line-0.7, -5, line-0.6, -5, length=0.025, angle=90);
	text(line-0.7, -8, labels=pval_stars[1], cex=0.45)
	arrows(line-0.6, -15, line-0.8, -15, length=0.025, angle=90); arrows(line-0.6, -15, line-0.4, -15, length=0.025, angle=90);
	text(line-0.6, -18, labels=pval_stars[2], cex=0.45)
	arrows(line-0.5, -25, line-0.8, -25, length=0.025, angle=90); arrows(line-0.5, -25, line-0.2, -25, length=0.025, angle=90);
	text(line-0.5, -28, labels=pval_stars[3], cex=0.45)

	if(line==1){
		text(-0.5, -6.5, labels='Gaussian (Sig. Level)', cex=0.45)
		text(-0.5, -16.5, labels='Linear (Sig. Level)', cex=0.45)
		text(-0.5, -26.5, labels='Simplified k (Sig. Level)', cex=0.45)
	}
}

dev.off()

write.csv(file='~/Desktop/AIC_outputs/2FVT_model_graphs/2FVT_culm_height_tests.csv', hgt_test, row.names=FALSE)
