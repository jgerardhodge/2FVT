#Note: This script serves as the final step in the 2FVT analysis workflow and largely serves as a tool for
#providing image outputs of the parameters fitted during analysis by 'plm_modeling_AIC_2FVT.R'.  As such the 
#assorted '..._2FVT_fitted_parameters.csv', '..._whole_plant_FVTs.csv', and 'All_2FVT_fitted_parameters.csv'
#serve as required input files prior to running. 

path='/Users/johnhodge/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/'
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
	whole_plant_fits<-read.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_parameters/', line, '_whole_plant_FVTs.csv', sep=''), row.names=1)
	cur_plant_attr<-c(whole_plant_fits[1,],whole_plant_fits[2,])
	names(cur_plant_attr)<-c('flw_xo', 'flw_k', 'flw_L', 'lfno_xo', 'lfno_k', 'lfno_L')
	whole_plant_attr<-rbind(whole_plant_attr, cur_plant_attr)
}

rownames(whole_plant_attr)=c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')

meta_attr<-read.csv(file='~/Desktop/2FVT_outputs/2FVT_parameters/All_2FVT_fitted_parameters.csv')

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


pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_lfno_w_WaldCI.pdf', height=7, width=9)

plot(c(0,max(x_fit)), c(0, round(max(as.numeric(whole_plant_attr[,6])))), xlab='Days', ylab='Leaf no', main='Leaf Number', col='white')
legend('topleft', legend=c(lines, method_title), col=c(line_cols,'black'), lwd=c(rep(2,length(line_cols)), 3), cex=0.8)
legend('bottomright', pch=22, legend='Flowering Day', cex=0.8)

for (iter in 1:length(fit_list)){
	CI_x=c(fit_list[[iter]][,1], rev(fit_list[[iter]][,1]))
	CI_y=c(fit_list[[iter]][,2]-(sig_lvl[3]*fit_list[[iter]][,3]), rev(fit_list[[iter]][,2]+(sig_lvl[3]*fit_list[[iter]][,3]))); CI_y[CI_y<0]=0; #Prevent CI bounds from extending below zero
	CI_bounds<-cbind(CI_x, CI_y)
	polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[iter], border=NA)	
	lines(fit_list[[iter]][,1], fit_list[[iter]][,2], col=line_cols[iter])

	flw_day=round(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[iter],1]),0)
	flw_size=fit_list[[iter]][as.numeric(fit_list[[iter]][,1])==flw_day,2]
	#Flowering Day
	points(flw_day, flw_size, col=line_cols[iter], pch=22)
}


flw_day=round(mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)=='A10' | rownames(whole_plant_attr)=='B100',1])),0)
flw_size=test_method[x_fit==flw_day]

lines(x_fit, test_method, lty=3)
points(flw_day, flw_size, pch=22)

dev.off()

method_fits<-cbind(x_fit, test_method)
fit_list[[length(lines)+1]]<-method_fits

lfno_list<-fit_list

#################################################################
#1B. Flowering time with confidence intervals by genotype
#################################################################

#Merge fitted GLM data into a single list to enable calling of measures between genotypes

fit_list<-list()

iter = 1

for (line in 1:length(lines)){
	L    =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],3])
	k    =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],2])
	xo   =as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],1])
	y_fit       =L/(1+exp(-k*(x_fit-xo)))

	y_p          =y_fit/max(y_fit);
	p            =(y_p+2)/(length(y_p)+4)
	CI_weight    =sqrt((p*(1-p))/length(y_fit)); 

	cur_data=cbind(x_fit, y_fit, CI_weight)

	fit_list[[iter]]<-cur_data	
	iter=iter+1
}

if(!is.null(method_p1) & !is.null(method_p2)){							#>>> Transgression test performed
	L<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,3]))
	k<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,2]))
	xo<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1 | rownames(whole_plant_attr)==method_p2,1]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))
}else if(!is.null(method_p1)){											#>>> Standard test performed
	L<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,3]))
	k<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,2]))
	xo<-mean(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==method_p1,1]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))
}else{																	#>>> Outlier test performed
	L<-mean(as.numeric(whole_plant_attr[,3]))
	k<-mean(as.numeric(whole_plant_attr[,2]))
	xo<-mean(as.numeric(whole_plant_attr[,1]))
	test_method<-L/(1+exp(-k*(x_fit-xo)))																	
}

pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_flw_w_WaldCI.pdf', height=7, width=9)

plot(c(0,max(x_fit)), c(0, round(max(as.numeric(whole_plant_attr[,3])))), xlab='Days', ylab='Flowering', main=paste(method_title, ' Genotype Leaf Number',sep=''), col='white')
legend('topleft', legend=c(lines, method_title), col=c(line_cols,'black'), lwd=c(rep(2,length(line_cols)), 3), cex=0.7)

for (iter in 1:length(fit_list)){
	CI_x=c(fit_list[[iter]][,1], rev(fit_list[[iter]][,1]))
	CI_y=c(fit_list[[iter]][,2]-(sig_lvl[3]*fit_list[[iter]][,3]), rev(fit_list[[iter]][,2]+(sig_lvl[3]*fit_list[[iter]][,3]))); CI_y[CI_y<0]=0; CI_y[CI_y>1]=1;#Prevent CI bounds from extending beyond bounds
	CI_bounds<-cbind(CI_x, CI_y)
	polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[iter], border=NA)	
	lines(fit_list[[iter]][,1], fit_list[[iter]][,2], col=line_cols[iter])
}

lines(x_fit, test_method, lty=3)

dev.off()

method_fits<-cbind(x_fit, test_method)
fit_list[[length(lines)+1]]<-method_fits

flw_list<-fit_list

#################################################################
#2. Segmented height with confidence intervals by genotype
#################################################################

fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/2FVT_outputs/2FVT_parameters/', lines[line], '_seg_culm_hgt_2FVT_fitted_parameters.csv', sep=''))
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

pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_seg_hgt_w_WaldCI.pdf', height=7, width=9)

for (phyt in 2:method_ceiling){

	plot(c(0,60), c(0,3), col='white', xlab='Days', ylab='Height (cm)', main=paste('Segmented Height at Phytomer ',phyt, sep=''))
	legend('topleft', legend=c(lines, method_title), col=c(line_cols,'black'), lwd=c(rep(2,length(line_cols)), 3), cex=0.7)

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
		
		if(length(y_fit)>0){		
			CI_x=c(fit_list[[iter]][,1], rev(fit_list[[iter]][,1]))
			CI_y=c(y_fit-(sig_lvl[3]*y_wald),rev(y_fit+(sig_lvl[3]*y_wald))); CI_y[CI_y<0]<-0
			CI_bounds<-cbind(CI_x, CI_y)

			polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[iter], border=NA)	
			lines(x_fit, y_fit, col=line_cols[iter])
		}
	}

	lines(x_fit, test_method, lty=3)

}

dev.off()

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

sghgt_list<-fit_list

#################################################################
# CULM HEIGHT OUTPUT
#################################################################

pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_total_hgt_w_WaldCI_pub.pdf', height=6, width=5)

plot(c(0,30),c(0,30),col='white', xlab='Days', ylab='Length (cm)', main='Main Culm Height')
legend('topleft', legend=c(lines), col=c(line_cols), lwd=c(rep(2,length(line_cols))), cex=0.7)

for (line in 1:length(lines)){
	
	culm_hgt=c()
	culm_CI=c()

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
		culm_CI<-c(culm_CI, sum_CI_weight)

	}

	CI_x=c(sghgt_list[[line]][,1], rev(sghgt_list[[line]][,1]))
	CI_y=c(culm_hgt-(sig_lvl[3]*culm_CI),rev(culm_hgt+(sig_lvl[3]*culm_CI))); CI_y[CI_y<0]<-0
	CI_bounds<-cbind(CI_x, CI_y)

	polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[line], border=NA)	

	lines(x_fit, culm_hgt, col=line_cols[line])

}

abline(v=c(15,25),lty=2)

dev.off()

#################################################################
#3. Leaf Distances with confidence intervals by genotype
#################################################################

fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/2FVT_outputs/2FVT_parameters/', lines[line], '_leaf_dist_2FVT_fitted_parameters.csv', sep=''))
	cur_n=as.numeric(N[line,2])
	
	geno_list=c()

	geno_ceiling=floor(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6]))
	geno_ceiling=geno_ceiling							#Add one unit to account for peduncle

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

pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_leaf_dist_w_WaldCI.pdf', height=7, width=9)

for (phyt in 2:method_ceiling){

	plot(c(0,30), c(0,12), col='white', xlab='Days', ylab='Leaf Tip-Ligule Distances (cm)', main=paste('Leaf Distance at Phytomer ',phyt, sep=''))
	legend('topleft', legend=c(lines, method_title), col=c(line_cols,'black'), lwd=c(rep(2,length(line_cols)), 3), cex=0.7)

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
		
		if(length(y_fit)>0){		
			CI_x=c(fit_list[[iter]][,1], rev(fit_list[[iter]][,1]))
			CI_y=c(y_fit-(sig_lvl[3]*y_wald),rev(y_fit+(sig_lvl[3]*y_wald))); CI_y[CI_y<0]<-0
			CI_bounds<-cbind(CI_x, CI_y)

			polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[iter], border=NA)	
			lines(x_fit, y_fit, col=line_cols[iter])
		}
	}

	lines(x_fit, test_method, lty=3)

}

dev.off()

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

lfdist_list<-fit_list

#################################################################
#4. Leaf Angles with confidence intervals by genotype
#################################################################

fit_list<-list(); method_pars=c()

iter = 1

for (line in 1:length(lines)){

	cur_pars=read.csv(paste('~/Desktop/2FVT_outputs/2FVT_parameters/', lines[line], '_leaf_angles_2FVT_fitted_parameters.csv', sep=''))
	cur_n=as.numeric(N[line,2])
	
	geno_list=c()

	geno_ceiling=floor(as.numeric(whole_plant_attr[rownames(whole_plant_attr)==lines[line],6]))
	geno_ceiling=geno_ceiling							#Add one unit to account for peduncle

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

pdf(file='~/Desktop/2FVT_outputs/2FVT_model_graphs/2FVT_leaf_angle_w_WaldCI.pdf', height=7, width=9)

for (phyt in 2:method_ceiling){

	plot(c(0,30), c(0,110), col='white', xlab='Days', ylab='Angle', main=paste('Leaf Angles at Phytomer ',phyt, sep=''))
	legend('topleft', legend=c(lines, method_title), col=c(line_cols,'black'), lwd=c(rep(2,length(line_cols)), 3), cex=0.7)

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
		
		if(length(y_fit)>0){		
			CI_x=c(fit_list[[iter]][,1], rev(fit_list[[iter]][,1]))
			CI_y=c(y_fit-(sig_lvl[3]*y_wald),rev(y_fit+(sig_lvl[3]*y_wald))); CI_y[CI_y<0]<-0
			CI_bounds<-cbind(CI_x, CI_y)

			polygon(CI_bounds[,1], CI_bounds[,2], col=CI_cols[iter], border=NA)	
			lines(x_fit, y_fit, col=line_cols[iter])
		}
	}

	lines(x_fit, test_method, lty=3)

}

dev.off()

method_fits<-cbind(x_fit, method_fits)
fit_list[[length(lines)+1]]<-method_fits

lfangle_list<-fit_list

#################################################################
#5. Tiller growth with confidence intervals by genotype
#################################################################




#################################################################
# IN SILICO MORPH OUTPUT
#################################################################

pdf.fn=paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/Simulated_morphs.pdf')
#pdf.fn=c()
if(length(pdf.fn)>0){
	pdf(file=pdf.fn, height=5, width=14)
}

plotsize=40*length(lines)

for (day in 11:(length(x_fit)-1)){	

	if(length(pdf.fn)<=0){
		if((day-10)<10){
			fn=paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/frames/frame_00', day-10,'.png',sep='')
		}else if((day-10)<100){
			fn=paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/frames/frame_0', day-10,'.png',sep='')
		}else{
			fn=paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/frames/frame_', day-10,'.png',sep='')	
		}
	
		png(file=fn, height=500, width=1500)
	}

	plot(c(0,plotsize), c(0,35), main=paste(method_title, ' Simulated morphs ', x_fit[day], ' days', sep=''), xlab='', ylab='(cm)', col='white', xaxt='n')
	axis(1, at=(1:length(lines)*40)-20, labels=lines)

	legend_x=4; legend_y=30
	for(heat in 1:nrow(sig_palette)){
		rect(legend_x, legend_y, legend_x+2, legend_y+2, col=sig_palette[heat,2], border=NA)
		legend_x=legend_x+2

		if(heat==1){
			text(legend_x-1, legend_y-0.8, sig_palette[heat,1], cex=0.8)
		}
		if(heat==ceiling(nrow(sig_palette)/2)){
			text(legend_x-1, legend_y-0.8, sig_palette[heat,1], cex=0.8)
			text(legend_x-1, legend_y+2.8, 'z', cex=1.25)
		}
		if(heat==nrow(sig_palette)){
			text(legend_x-1, legend_y-0.8, sig_palette[heat,1], cex=0.8)
		}
	}

	for (line in 1:length(lines)){

		apex=0
		base=(line*40)-20

		flw=flw_list[[line]][day,]
		flw_std=flw_list[[length(lines)+1]][day,2]

		if(flw[2]>=0.95){
			heading=1
		}else{
			heading=0
		}

		leaves=floor(lfno_list[[line]][day,2])

		width=seq(1.2+(leaves-2)*0.05, 1.2, by=-0.05)

		culm_segments=sghgt_list[[line]][day,]
		culm_ligules=lfangle_list[[line]][day,]
		leaf_dists=lfdist_list[[line]][day,]

		hgt_standard=sghgt_list[[length(lines)+1]][day,]
		lig_standard=lfangle_list[[length(lines)+1]][day,]
		dis_standard=lfdist_list[[length(lines)+1]][day,]

		for (phyt in 2:(leaves-1)){

			hgt     =culm_segments[names(culm_segments)==paste('y_fit', phyt, sep='')]
			hgt_wald=culm_segments[names(culm_segments)==paste('y_Wald', phyt, sep='')]
			ang     =culm_ligules[names(culm_ligules)==paste('y_fit', phyt, sep='')]
			ang_wald=culm_ligules[names(culm_ligules)==paste('y_Wald', phyt, sep='')]
			dis     =leaf_dists[names(leaf_dists)==paste('y_fit', phyt, sep='')]
			dis_wald=leaf_dists[names(leaf_dists)==paste('y_Wald', phyt, sep='')]

			hgt_std=hgt_standard[names(hgt_standard)==paste('standard', phyt, sep='')]
			ang_std=lig_standard[names(lig_standard)==paste('standard', phyt, sep='')]
			dis_std=dis_standard[names(dis_standard)==paste('standard', phyt, sep='')]	

			hgt_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=hgt-sig_lvl[Z]*hgt_wald
				CI_hi=hgt+sig_lvl[Z]*hgt_wald
				if (CI_hi<hgt_std){hgt_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>hgt_std){hgt_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			ang_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=ang-sig_lvl[Z]*ang_wald
				CI_hi=ang+sig_lvl[Z]*ang_wald
				if (CI_hi<ang_std){ang_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>ang_std){ang_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			dis_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=dis-sig_lvl[Z]*dis_wald
				CI_hi=dis+sig_lvl[Z]*dis_wald
				if (CI_hi<dis_std){dis_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>dis_std){dis_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			rect(base-width[phyt-1]/2, apex, base+width[phyt-1]/2, apex+hgt, col=hgt_col)
			apex=apex+hgt
			rect(base-width[phyt-1]/2, apex, base+width[phyt-1]/2, apex+0.3, col=ang_col)
			apex=apex+0.3

			if(phyt %% 2){
				lig=base-width[phyt-1]/2
				theta<-ang*(pi/180)
				#hypotenuse<-sim_lf_dist[sim_row,colnames(sim_lf_dist) %in% colnames(sim_seg_culm_hgt)[s]]
				lf_tip_y_offset<-abs(cos(theta)*dis)	#Based on Sin(theta)=H/O Rule
				lf_tip_x_offset<-abs(sin(theta)*dis)	#Based on Cos(theta)=H/A Rule
				lines(c(lig, lig-lf_tip_x_offset), c(apex, apex+lf_tip_y_offset), col=dis_col)
			} else{
				lig=base+width[phyt-1]/2
				theta<-ang*(pi/180)
				#hypotenuse<-sim_lf_dist[sim_row,colnames(sim_lf_dist) %in% colnames(sim_seg_culm_hgt)[s]]
				lf_tip_y_offset<-abs(cos(theta)*dis)	#Based on Sin(theta)=H/O Rule
				lf_tip_x_offset<-abs(sin(theta)*dis)	#Based on Cos(theta)=H/A Rule		
				lines(c(lig, lig+lf_tip_x_offset), c(apex, apex+lf_tip_y_offset), col=dis_col)
			}
		}



		if (heading==0){
			dis=leaf_dists[names(leaf_dists)==paste('y_fit', leaves, sep='')]
			dis_wald=leaf_dists[names(leaf_dists)==paste('y_Wald', leaves, sep='')]
			dis_std=dis_standard[names(dis_standard)==paste('standard', leaves, sep='')]	

			dis_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=dis-sig_lvl[Z]*dis_wald
				CI_hi=dis+sig_lvl[Z]*dis_wald
				if (CI_hi<dis_std){dis_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>dis_std){dis_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			lines(c(base, base), c(apex, apex+dis), col=dis_col)
		} else {
			hgt     =culm_segments[names(culm_segments)==paste('y_fit', leaves, sep='')]
			hgt_wald=culm_segments[names(culm_segments)==paste('y_Wald', leaves, sep='')]
			ang     =culm_ligules[names(culm_ligules)==paste('y_fit', leaves, sep='')]
			ang_wald=culm_ligules[names(culm_ligules)==paste('y_Wald', leaves, sep='')]
			dis     =leaf_dists[names(leaf_dists)==paste('y_fit', leaves, sep='')]
			dis_wald=leaf_dists[names(leaf_dists)==paste('y_Wald', leaves, sep='')]

			dis_std=dis_standard[names(dis_standard)==paste('standard', leaves, sep='')]	
			ang_std=lig_standard[names(lig_standard)==paste('standard', leaves, sep='')]
			dis_std=dis_standard[names(dis_standard)==paste('standard', leaves, sep='')]	

			hgt_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=hgt-sig_lvl[Z]*hgt_wald
				CI_hi=hgt+sig_lvl[Z]*hgt_wald
				if (CI_hi<hgt_std){hgt_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>hgt_std){hgt_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			ang_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=ang-sig_lvl[Z]*ang_wald
				CI_hi=ang+sig_lvl[Z]*ang_wald
				if (CI_hi<ang_std){ang_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>ang_std){ang_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			dis_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=dis-sig_lvl[Z]*dis_wald
				CI_hi=dis+sig_lvl[Z]*dis_wald
				if (CI_hi<dis_std){dis_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]}
				if(CI_hi>dis_std){dis_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]}
			}

			#Final Vegetative phytomer
			rect(base-width[phyt-1]/2, apex, base+width[phyt-1]/2, apex+hgt, col=hgt_col)
			apex=apex+hgt
			rect(base-width[phyt-1]/2, apex, base+width[phyt-1]/2, apex+0.3, col=ang_col)
			apex=apex+0.3	

			#Peduncle and inflorescence
			ped=culm_segments[names(culm_segments)==paste('y_fit', leaves+1, sep='')]
			rect(base-0.2, apex, base+0.2, apex+ped)
			apex=apex

			if(leaves %% 2){
				lig=base-0.6
				theta<-ang*(pi/180)
				#hypotenuse<-sim_lf_dist[sim_row,colnames(sim_lf_dist) %in% colnames(sim_seg_culm_hgt)[s]]
				lf_tip_y_offset<-abs(cos(theta)*dis)	#Based on Sin(theta)=H/O Rule
				lf_tip_x_offset<-abs(sin(theta)*dis)	#Based on Cos(theta)=H/A Rule
				lines(c(lig, lig-lf_tip_x_offset), c(apex, apex+lf_tip_y_offset), col=dis_col)
			} else{
				lig=base+0.6
				theta<-ang*(pi/180)
				#hypotenuse<-sim_lf_dist[sim_row,colnames(sim_lf_dist) %in% colnames(sim_seg_culm_hgt)[s]]
				lf_tip_y_offset<-abs(cos(theta)*dis)	#Based on Sin(theta)=H/O Rule
				lf_tip_x_offset<-abs(sin(theta)*dis)	#Based on Cos(theta)=H/A Rule		
				lines(c(lig, lig+lf_tip_x_offset), c(apex, apex+lf_tip_y_offset), col=dis_col)
			}

			flw_col='black'

			for (Z in 1:length(sig_lvl)){	
				CI_lo=flw[2]-sig_lvl[Z]*flw[3]
				CI_hi=flw[2]+sig_lvl[Z]*flw[3]
				if(flw_std>CI_lo & flw_std<CI_hi){
					flw_col='black'
				}else if(CI_hi<flw_std){
					flw_col=sig_palette[sig_palette[,1]==-sig_lvl[Z],2]
				}else if(CI_lo>flw_std){
					flw_col=sig_palette[sig_palette[,1]==sig_lvl[Z],2]
				}
			}

			points(x=base, y=apex+ped+0.1, pch=19, cex=2)
			points(x=base, y=apex+ped+0.1, pch=19, cex=1.55, col=flw_col)	
		}
	}
	if(length(pdf.fn)<=0){
		dev.off()
	}
}

if(length(pdf.fn)>0){
	dev.off()
}


