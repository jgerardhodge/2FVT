#Note: This code is used to estimate the overall accuracy of the original indiviually fitted models for each organ unit found within the 
#'plm_culm_segmentation.R' script as well as the two competing 2FVT models and compare them against one another through the use of 
#jack-knifing each trait group so that the delta AIC's can be compared.  This script is largely for the use of validating the legitimacy of 
#the use of 2FVT traits rather than for any actual down stream uses.  Buyer beware!

#####
# Phase 0: Specify the current line to test AIC values for model fits on
#####

Sec.FVT.pars<-c()
Sec.FVT.AICs<-c()
Sec.FVT.pred.flw<-c()

verbose=TRUE

scaler=62.7 #pixel/cm conversion

for (line in c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')){

	print(paste('Beginning run on genotype: ', line))

	if (line=='A10'){reps=c(1:10)}
	if (line=='B100'){reps=c(1:10)}
	if (line=='RIL159'){reps=c(1:6,9:10)}
	if (line=='RIL110'){reps=c(1:5,7:10)}
	if (line=='RIL39'){reps=c(1:6,8:10)}


	anno_data<-read.csv(paste('~/Dropbox/Acute_manuscript/Data/homology_groups/',line,'_homology_groups_cady_annotated.csv',sep=''))

	#####
	# Phase 1: Extract Biologically meaningful measures to test from annotated plm homology lists for the RIL
	#####

	dup_sanity_check<-0

	#Sanity check for duplicates
	for (rep in reps){
		#Extract filenames for images corresponding to the current replicate
		img_filenames<-unique(anno_data[grepl(paste('fixed_mask_', line, '_rep', rep, '_',sep=''), anno_data[,3]),3])
		#Remove file suffix
		img_names<-strsplit(as.character(img_filenames), '.png')
		img_names<-unlist(img_names)
		#Split name at day substring in name
		days_list<-strsplit(img_names, '_d')
		#Extract days and convert to a numeric vector, remove first day '[-1]' to allow for ease of iteration later at ends of vector
		days<-as.numeric(unlist(days_list)[2*1:length(img_filenames)])[-1]

		for (day in days){
			slice<-anno_data[grepl(paste(line,'_rep', rep, '_d', day, sep=''), anno_data[,2]),]
			#Identify duplicates which aren't unannotated (i.e. '-' values)
			dups<-duplicated(slice[,2]) & slice[,2]!='-'
			if(sum(dups)>0){
				dup_sanity_check<-1
				print(as.factor(paste('Rep', rep, 'day', day, 'duplicates:')), max.levels=0)
				print(slice[dups,2], max.levels=0)
			}
		}
		#print(as.factor(''), max.levels=0)
		
	}

	if (dup_sanity_check==1){
		print('Duplicate sanity check failed, screen duplicates from dataset then try to rerun analysis...')
	} else {
		combined_leaf_dist<-c()
		combined_internode_lengths<-c()
		combined_ligule_angles<-c()
		combined_heading_counts<-c()
		combined_leaf_no<-c()

		for (rep in reps){
			#Extract filenames for images corresponding to the current replicate
			img_filenames<-unique(anno_data[grepl(paste('fixed_mask_', line, '_rep', rep, '_',sep=''), anno_data[,3]),3])
			#Remove file suffix
			img_names<-strsplit(as.character(img_filenames), '.png')
			img_names<-unlist(img_names)
			#Split name at day substring in name
			days_list<-strsplit(img_names, '_d')
			#Extract days and convert to a numeric vector, remove first day '[-1]' to allow for ease of iteration later at ends of vector
			days<-as.numeric(unlist(days_list)[2*1:length(img_filenames)])[-1]

			heading_counts<-c()

			for (day in days){
				#img<-load.image(paste('/Volumes/Image_Seq_Data/Oct_Nov_tillering_RILs/image_masks/',line,'/',line,'_rep',rep,'/',line,'_rep',rep,'_d',day,'.png',sep=''))

				slice<-anno_data[grepl(paste(line,'_rep', rep, '_d', day, sep=''), anno_data[,2]),]

				culm_lfs<-slice[grepl('leaf', slice[,1]) & !grepl('till', slice[,1]),1]
				culm_ligs<-slice[grepl('lig', slice[,1]) & !grepl('till', slice[,1]),1]

				max_leaf<-0

				for (cur_lf in culm_lfs){
					cur_max<-as.numeric(strsplit(as.character(cur_lf), 'leaf')[[1]][2])
					if(cur_max>max_leaf){
						max_leaf<-cur_max
					}
				}

				for (cur_lig in culm_ligs){
					cur_max<-as.numeric(strsplit(as.character(cur_lig), 'lig')[[1]][2])
					if(cur_max>max_leaf){
						max_leaf<-cur_max
					}
				}				

				combined_leaf_no<-rbind(combined_leaf_no, c(day, max_leaf))
				rownames(combined_leaf_no)[nrow(combined_leaf_no)]<-paste('rep',rep,'_d', day,sep='')

				#Check if replicate has headed at this time point 
				heading_counts<-rbind(heading_counts, c(day, sum(slice$group=='inflo')))
				rownames(heading_counts)[nrow(heading_counts)]<-paste('rep',rep,'_d', day,sep='')

				#Flip orientation of y-axis from pythonic display
				slice[,5]<-2592-slice[,5];  #NOTE THIS IS SPECIFIC TO THE RESOLUTION OF THE CAMERA USED!
				slice[,7]<-2592-slice[,7];	#NOTE THIS IS SPECIFIC TO THE RESOLUTION OF THE CAMERA USED! 
				slice[,9]<-2592-slice[,9];	#NOTE THIS IS SPECIFIC TO THE RESOLUTION OF THE CAMERA USED!

				culm_lig<-slice[grepl('lig', slice$group) & !grepl('till', slice$group),]
		        
				if(nrow(culm_lig)>=3){
					x<-culm_lig[,4]
					y<-culm_lig[,5]

					lig_lm<-lm(y~x)

					b=lig_lm$coef[1]
					m=lig_lm$coef[2]
					ort_m=-(1/m)
					ort_b=y-(ort_m*x)

					segments<-c()

					for(n in 1:length(ort_b)){
						#abline(c(ort_b[n], ort_m), col='red', lty=3)
						cm<-rbind(c(b, m), c(ort_b[n], ort_m))
						segments<-rbind(segments, c(-solve(cbind(cm[,2],-1)) %*% cm[,1]))
					}

					#points(segments, col='red', pch=20)

					rownames(segments)<-culm_lig$group
					segments<-segments[order(rownames(segments)),]

					internode_lengths<-c()

					for(p in 1:(nrow(segments)-1)){
						dist<-sqrt(((segments[p,1]-segments[p+1,1])^2)+((segments[p,2]-segments[p+1,2])^2))

						pos1<-as.numeric(strsplit(rownames(segments)[p], 'lig')[[1]][2])
						pos2<-as.numeric(strsplit(rownames(segments)[p+1], 'lig')[[1]][2])

						#Sanity check to prevent measures from across multiple phytomers
						if(abs(pos1-pos2)==1){ 
							internode_lengths<-rbind(internode_lengths, dist)
							rownames(internode_lengths)[nrow(internode_lengths)]<-paste('rep', rep, '_', rownames(segments)[p],'_', rownames(segments)[p+1], sep='')
						}
					}

					if (!is.null(internode_lengths)){

						time<-rep(day, nrow(internode_lengths))
						internode_lengths<-cbind(time, internode_lengths)
						internode_lengths[,1]<-internode_lengths[,1]
						internode_lengths[,2]<-internode_lengths[,2]/scaler
						combined_internode_lengths<-rbind(combined_internode_lengths, internode_lengths)

					}

					culm_apex=c(((segments[nrow(segments),2]+100)-b)/m,segments[nrow(segments),2]+100)
					
					for (l in 1:nrow(segments)){
						cur_leaf<-paste('leaf',strsplit(rownames(segments)[l], 'lig')[[1]][2],sep='')
						cur_lig<-rownames(segments)[l]
						if (sum(slice$group %in% cur_leaf)==1 & sum(slice$group %in% cur_lig)==1){
							leaf_coord<-as.vector(slice[slice$group %in% cur_leaf,4:5])
							lig_coord<-as.vector(slice[slice$group %in% cur_lig,4:5])
							P12 = sqrt(((leaf_coord[1]-lig_coord[1])^2)+((leaf_coord[2]-lig_coord[2])^2))
					        P13 = sqrt(((leaf_coord[1]-culm_apex[1])^2)+((leaf_coord[2]-culm_apex[2])^2))
					        P23 = sqrt(((lig_coord[1]-culm_apex[1])^2)+((lig_coord[2]-culm_apex[2])^2))
					        dot = (P12^2 + P13^2 - P23^2)/(2*P12*P13)
					        ang = acos(dot)*(180/pi)
					        combined_leaf_dist<-rbind(combined_leaf_dist, cbind(day, P12/scaler))
					        rownames(combined_leaf_dist)[nrow(combined_leaf_dist)]<-paste('rep',rep,'_d', day, '_', cur_leaf,sep='')
					        combined_ligule_angles<-rbind(combined_ligule_angles, cbind(day, ang))
					        rownames(combined_ligule_angles)[nrow(combined_ligule_angles)]<-paste('rep',rep,'_d', day, '_', cur_lig,sep='')
						#} else {
						#	print(paste('Rejecting rep', rep, 'day', day, rownames(segments)[l]))
						}
					}	
				}
			}

			if (sum(heading_counts[,2])==0 & line!='B100'){
				heading_counts[nrow(heading_counts),2]<-1
				combined_heading_counts<-rbind(combined_heading_counts, heading_counts)
			} else {
				combined_heading_counts<-rbind(combined_heading_counts, heading_counts)
			}
		}


		#####
		# Phase 2: Jack-knife AIC values for each of the individualized model traits then save to files
		#####

		#####
		# Whole Plant Trait Screening
		#####

		if (verbose==TRUE){print('Fitting for Whole Stem Attributes (Flowering and Leaf number)')}

		pdf(paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/', line,'_whole_plant_FVTs.pdf',sep=''), height=8, width=6)

		layout(matrix(c(1,2,3), 3, 1, byrow=TRUE))

		#################
		#Flowering Time #
		#################

		plot(c(0,30),c(0, 1), col='white', xlab='Days', ylab='Flowering', main=paste(line, 'Days to Heading'))

		if(line=='B100'){     #Special case given B100 was not grown to heading
			L_min=0; L_max=1;
		}else{
			L_min=1; L_max=1;
		}

		growth=list(dist = as.numeric(combined_heading_counts[,2]), age = as.numeric(combined_heading_counts[,1]))

		#1/lm(growth$hgt[growth$hgt!=0]~growth$age[growth$hgt!=0])$coef[2]
		start_params=list(xo = 10, k = 0.001, L = mean(growth$dist[growth$dist!=0]))
		lo_par_bound=list(xo = 2, k = 0, L = 1)
		hi_par_bound=list(xo = summary(growth$age[growth$dist!=0])[[5]]+4, k = 3, L = 1)

		if(line!='B100'){     #Special case given B100 was not grown to heading	
			nlm_fit<-nls(growth$dist ~ L/(1+exp(-k*(growth$age-xo))), data = growth, start=start_params, algorithm='port', 
				trace=F, lower=lo_par_bound, upper=hi_par_bound)
			flw_coefs = coef(nlm_fit)	
		}else{
			flw_coefs = c(xo=53, k=3, L=1)
		}
				
		x_fit = seq(0, 30, by = 0.1); y_fit = flw_coefs[3]/(1+exp(-flw_coefs[2]*(x_fit-flw_coefs[1])));
		y_p         =  y_fit/max(y_fit); 
		CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
		flw_model_fits  =  cbind(x_fit, y_fit, CI_weight)
		colnames(flw_model_fits)=c('days', 'flw', 'CI_weight')	
		lines(x_fit, y_fit, col=rgb(0.2,0.2,0.2,0.4), lwd=2)
		points(combined_heading_counts, col=rgb(0.7,0.7,0.7,0.2), pch=19, cex=0.4)	



		##############
		#Leaf Number #
		##############

		plot(c(0,30),c(0, max(combined_leaf_no[,2])*1.2), col='white', xlab='Days', ylab='Leaf no.', main=paste(line, 'Leaf Number'))

		growth=list(dist = as.numeric(combined_leaf_no[,2]), age = as.numeric(combined_leaf_no[,1]))

		#1/lm(growth$hgt[growth$hgt!=0]~growth$age[growth$hgt!=0])$coef[2]
		start_params=list(xo = 4, k = 1, L = mean(growth$dist[growth$dist!=0]))
		lo_par_bound=list(xo = 2, k = 0.0001, L = 0)
		hi_par_bound=list(xo = summary(growth$age[growth$dist!=0])[[5]]+4, k = 10000, L = max(combined_leaf_no[,2])*1.5)

		#print(paste(line, 'Leaf number'))
		#print(rbind(start_params, lo_par_bound, hi_par_bound))

		nlm_fit<-nls(growth$dist ~ L/(1+exp(-k*(growth$age-xo))), data = growth, start=start_params, algorithm='port', 
			trace=F, lower=lo_par_bound, upper=hi_par_bound)

		lfno_coefs = coef(nlm_fit)			
		x_fit = seq(0, 30, by = 0.1); y_fit = lfno_coefs[3]/(1+exp(-lfno_coefs[2]*(x_fit-lfno_coefs[1]))); flw_y_fit=y_fit;
		y_p         =  y_fit/max(y_fit); 
		CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
		lfno_model_fits  =  cbind(x_fit, y_fit, CI_weight)
		colnames(lfno_model_fits)=c('days', 'leaf_no', 'CI_weight')			
		lines(x_fit, y_fit, col=rgb(0.2,0.6,0.2,0.4), lwd=2)
		points(combined_leaf_no, col=rgb(0.3,0.7,0.3,0.2), pch=19, cex=0.4)	



		#########################
		#Prediction starts here!#
		#########################

		if (line!='B100'){
			plot(c(0,30),c(0, max(combined_leaf_no[,2])*1.5), col='white', xlab='Days', ylab='Leaf no.', main=paste(line, 'Predicted Days to Heading'))

			gaps=c(0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
			x_fit = seq(0, 30, by = 0.1); y_fit = lfno_coefs[3]/(1+exp(-lfno_coefs[2]*(x_fit-lfno_coefs[1]))); flw_y_fit=y_fit;
			fitted_flw_date=min(x_fit[flw_y_fit>=max(floor(flw_y_fit*1))])

		}else{
			plot(c(0,60),c(0, max(combined_leaf_no[,2])*1.5), col='white', xlab='Days', ylab='Leaf no.', main=paste(line, 'Predicted Days to Heading'))

			gaps=1
			x_fit = seq(0, 60, by = 0.1); y_fit = lfno_coefs[3]/(1+exp(-lfno_coefs[2]*(x_fit-lfno_coefs[1]))); flw_y_fit=y_fit;
			fitted_flw_date=min(x_fit[flw_y_fit>=max(floor(flw_y_fit*1))])
		}

		legend('topleft', lty=rep(1,length(gaps)), col=rgb(1-gaps,gaps-0.7,0.5,0.5), 
			legend=paste(gaps*100, '% Ontogeny Used'), cex=0.55)

		for (gap in gaps){

			pred_growth=list(dist=growth$dist[growth$dist<=max(floor(flw_y_fit*gap))], age=growth$age[growth$dist<=max(floor(flw_y_fit*gap))])

			start_params=list(xo = 8, k = 1, L = mean(pred_growth$dist[pred_growth$dist!=0]))
			lo_par_bound=list(xo = 3, k = 0.0001, L = 0)
			hi_par_bound=list(xo = summary(pred_growth$age[pred_growth$dist!=0])[[5]]+7, k = 10000, L = max(combined_leaf_no[,2])*1.5)

			#print(paste(line, 'Leaf number'))
			#print(rbind(start_params, lo_par_bound, hi_par_bound))

			nlm_fit<-nls(pred_growth$dist ~ L/(1+exp(-k*(pred_growth$age-xo))), data = pred_growth, start=start_params, algorithm='port', 
				trace=F, lower=lo_par_bound, upper=hi_par_bound)

			pred_coefs = coef(nlm_fit)	

			y_fit = pred_coefs[3]/(1+exp(-pred_coefs[2]*(x_fit-pred_coefs[1])));
			cur_pred=min(x_fit[round(y_fit)>=max(round(y_fit))])

			lines(x_fit, y_fit, col=rgb(1-gap,gap-0.7,0.5,0.5), lwd=2)
			lines(c(cur_pred,cur_pred), c(0,y_fit[which(cur_pred==x_fit)]), lty=3, col=rgb(1-gap,gap-0.7,0.5,0.5))

			if (gap==0.7 & line!='B100'){
				flw_predict<-c(fitted_flw_date, cur_pred)
			}else{
				flw_predict<-c(flw_predict, cur_pred)
			}
			if (line=='B100'){
				flw_predict<-c(fitted_flw_date, '-', '-', '-', '-', '-', '-', cur_pred)
			}

		}

		Sec.FVT.pred.flw<-rbind(Sec.FVT.pred.flw, flw_predict)

		dev.off()

		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_parameters/', line,'_whole_plant_FVTs.csv',sep=''), rbind(flw_coefs, lfno_coefs))
		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/', line, '_flw_fits_w_Wald_CI_weights.csv', sep=''), flw_model_fits, row.names=FALSE)
		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/', line, '_lfno_fits_w_Wald_CI_weights.csv', sep=''), lfno_model_fits, row.names=FALSE)


		#####
		# Subcomponent Block A: Segmented Culm Height Screening
		#####

		pdf(paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/seg_culm_hgt_', line,'_2FVT.pdf',sep=''), height=12, width=10)

		if (verbose==TRUE){print('Initializing block A')}

		layout(matrix(c(1,1,1,7,7,7,6,6,6,5,5,5,2,2,3,3,4,4), 3, 6, byrow=TRUE))

		plot(c(0,30),c(0, max(combined_internode_lengths[,2])), col='white', xlab='Days', ylab='Ligular distance (cm)', main=paste(line, 'Segmented Height (Independent Models)'))

		cols=colorRampPalette(c("darkgreen", "orange", "red"))(8)
		legend('topleft', pch=rep(19, 8), col=cols, legend=paste('Phytomer ', 2:9), cex=0.8)

		sghgt_coefs<-c();

		if(line=='A10'){xo_min=4; gap_days<-2; k_min=5; k_bounds=c(0.0001,10000); kg=35}
		if(line=='RIL39'){xo_min=3; gap_days<-3; k_min=1.8; k_bounds=c(0.0001,10000); kg=20}
		if(line=='RIL110'){xo_min=3; gap_days<-3; k_min=1; k_bounds=c(0.0001,10000); kg=35}
		if(line=='RIL159'){xo_min=3; gap_days<-2; k_min=1; k_bounds=c(0.0001,10000); kg=35}
		if(line=='B100'){xo_min=3; gap_days<-2; k_min=1; k_bounds=c(0.0001,10000); kg=10}

		for (u in 3:8){
			
			unit_measures<-c()
			unit<-paste('_lig', u, '_lig', u+1, sep='')
			time_courses<-combined_internode_lengths[grepl(unit, rownames(combined_internode_lengths)),]

			if (nrow(time_courses)>0){ #Ensure units don't advance beyond the limits of phytomers that exist within the genotype

				for (rep in reps){
					matches<-grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
					time_series<-time_courses[matches,]
					
					#print(paste('Unit',u,'Rep',rep))
					if(sum(matches)==0){
						#print(paste('Unit',u,'Rep',rep,': No hits'))
						next
					}else if(sum(matches)/length(days)<0.5){
						if (sum(matches)==1){
							#print(paste('Unit',u,'Rep',rep,': One hit'))
							if(4<(min(time_series[1])-3)){
								missing_days<-seq(4,min(time_series[1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							time_series<-cbind(c(missing_days, time_series[1]), c(missing_dist, time_series[2]))
							rownames(time_series)<-rep(paste('rep',rep,'_lig',u,'_lig',u+1, sep=''), nrow(time_series))
							unit_measures<-rbind(unit_measures, time_series)	
						}else{
							#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
							if(4<(min(time_series[,1])-3)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							time_series<-rbind(cbind(missing_days, missing_dist), time_series)
							rownames(time_series)<-rep(paste('rep',rep,'_lig',u,'_lig',u+1, sep=''), nrow(time_series))
							unit_measures<-rbind(unit_measures, time_series)
						}
					}else if(sum(matches)/length(days)>=0.5){	
						#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
						if(4<(min(time_series[,1])-3)){
							missing_days<-seq(4,min(time_series[,1])-gap_days)
							missing_dist<-c(rep(0, length(missing_days)))
						}else{
							missing_days<-4; missing_dist<-0;
						}
						time_series<-rbind(cbind(missing_days, missing_dist), time_series)
						rownames(time_series)<-rep(paste('rep',rep,'_lig',u,'_lig',u+1, sep=''), nrow(time_series))
						unit_measures<-rbind(unit_measures, time_series)	
					}
				}

				growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
				
				#1/lm(growth$hgt[growth$hgt!=0]~growth$age[growth$hgt!=0])$coef[2]
				start_params=list(xo = xo_min, k = k_min, L = mean(growth$dist[growth$dist!=0]))
				lo_par_bound=list(xo = floor(xo_min-1), k = k_bounds[1], L = 0)
				hi_par_bound=list(xo = summary(growth$age[growth$dist!=0])[[5]]+4, k = k_bounds[2], L = max(growth$dist[growth$dist!=0])*1.5)

				print(paste(line, 'segmented heights'))
				print(rbind(start_params, lo_par_bound, hi_par_bound))

				nlm_fit<-nls(growth$dist ~ L/(1+exp(-k*(growth$age-xo))), data = growth, start=start_params, algorithm='port', 
					trace=T, lower=lo_par_bound, upper=hi_par_bound)
				
				coefs = coef(nlm_fit)
				
				#Estimate a fitted curve based on coefficients of nlm
				x_fit = seq(0, 30, by = 0.1); y_fit = coefs[3]/(1+exp(-coefs[2]*(x_fit-coefs[1])));
				x_scale = x_fit[y_fit>=1][1]
				adj_x_fit = x_fit[y_fit>=1]-x_scale; adj_y_fit = y_fit[y_fit>=1];

				points(growth$age, growth$dist, cex=0.3, col=cols[u-1])
				lines(x_fit, y_fit, col=cols[u-1])
				coefs = coef(nlm_fit); sghgt_coefs<-rbind(sghgt_coefs,c(u , coefs));

				xo_min<-sghgt_coefs[nrow(sghgt_coefs),2]

			}
		}

		colnames(sghgt_coefs)<-c('Phytomer', 'xo (ind)', 'k (ind)', 'L (ind)')

		#Use linear regression to estimate base k parameter for a given phytomer
		lm_k_fit<-lm(sghgt_coefs[,3]~sghgt_coefs[,1])
		cooks_dist_index_k<-cooks.distance(lm_k_fit)>4*mean(cooks.distance(lm_k_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_k)>0){
			lm_k_fit<-lm(sghgt_coefs[!cooks_dist_index_k,3]~sghgt_coefs[!cooks_dist_index_k,1])
		}

		add_coefs<-cbind(rep(mean(sghgt_coefs[!cooks_dist_index_k,3]), nrow(sghgt_coefs)), rep(lm_k_fit$coef[1], nrow(sghgt_coefs)), rep(lm_k_fit$coef[2], nrow(sghgt_coefs)))
		colnames(add_coefs)<-c('Simple k (2FVT)', 'k-b (2FVT)', 'k-m (2FVT)')
		sghgt_coefs<-cbind(sghgt_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(sghgt_coefs[,c(1,3)], xlim=c(0,max(sghgt_coefs[,1])), ylim=c(0,max(sghgt_coefs[,3])), main='Max. Phytomer Growth Rate', xlab='Phytomer no.', ylab='k')
		points(sghgt_coefs[cooks_dist_index_k,1], sghgt_coefs[cooks_dist_index_k,3], pch=19, col='red', cex=0.55); abline(lm_k_fit); abline(h=mean(sghgt_coefs[!cooks_dist_index_k,3]), lty=3)
		text(summary(1:max(sghgt_coefs[,1]))[2]*0.5, max(sghgt_coefs[,3])*0.8, paste('adj. R2:', round(summary(lm_k_fit)$adj.r.squared,3)))
		text(summary(1:max(sghgt_coefs[,1]))[2]*0.5, max(sghgt_coefs[,3])*0.7, paste('slope:', round(coef(lm_k_fit)[2],3)))
		legend('topleft', lty=c(1,3,0), pch=c(-1,-1,19), col=c('black','black','red'), legend=c('Complete', 'Simplified', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_xo_fit<-lm(sghgt_coefs[,1]~sghgt_coefs[,2])
		cooks_dist_index_xo<-cooks.distance(lm_xo_fit)>4*mean(cooks.distance(lm_xo_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_xo)>0){
			lm_xo_fit<-lm(sghgt_coefs[!cooks_dist_index_xo,1]~sghgt_coefs[!cooks_dist_index_xo,2])
		}

		add_coefs<-cbind(rep(lm_xo_fit$coef[1], nrow(sghgt_coefs)), rep(lm_xo_fit$coef[2], nrow(sghgt_coefs)))
		colnames(add_coefs)<-c('xo-b (2FVT)', 'xo-m (2FVT)')
		sghgt_coefs<-cbind(sghgt_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(sghgt_coefs[,c(2,1)], xlim=c(0,max(sghgt_coefs[,2])), ylim=c(0,max(sghgt_coefs[,1])), main='Rate of Phytomer Deposition', xlab='xo (days)', ylab='Phytomer no.')
		points(sghgt_coefs[cooks_dist_index_xo,2], sghgt_coefs[cooks_dist_index_xo,1], pch=19, col='red', cex=0.55); abline(lm_xo_fit);
		text(summary(1:max(sghgt_coefs[,2]))[2]*0.5, summary(1:max(sghgt_coefs[,1]))[5], paste('adj. R2:', round(summary(lm_xo_fit)$adj.r.squared,3)))
		text(summary(1:max(sghgt_coefs[,2]))[2]*0.5, summary(1:max(sghgt_coefs[,1]))[5]*.9, paste('slope:', round(coef(lm_xo_fit)[2],3)))
		legend('topleft', lty=c(1,0), pch=c(-1,19), col=c('black','red'), legend=c('Complete', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_L_fit<-lm(sghgt_coefs[,4]~sghgt_coefs[,1])
		cooks_dist_index_L<-cooks.distance(lm_L_fit)>4*mean(cooks.distance(lm_L_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		complex_L=list(y = as.numeric(sghgt_coefs[,4]), x = as.numeric(sghgt_coefs[,1]))

		gauss_params=list(xg = sghgt_coefs[which.max(sghgt_coefs[,4]),1], kg = kg, Lg = mean(sghgt_coefs[,4]))
		gauss_lo_par=list(xg = sghgt_coefs[which.max(sghgt_coefs[,4]),1]-2, kg = -1000, Lg = 0)
		gauss_hi_par=list(xg = sghgt_coefs[which.max(sghgt_coefs[,4]),1]+2, kg = 1000, Lg = max(sghgt_coefs[,4])*1.5)

		complex_L_fit<-nls(complex_L$y ~ Lg*exp(-(complex_L$x-xg)^2/(2*kg^2)), data = complex_L, start=gauss_params, algorithm='port', 
					trace=T, lower=gauss_lo_par, upper=gauss_hi_par)

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_L)>0){
			lm_L_fit<-lm(sghgt_coefs[!cooks_dist_index_L,4]~sghgt_coefs[!cooks_dist_index_L,1])
		}

		add_coefs<-cbind(rep(lm_L_fit$coef[1], nrow(sghgt_coefs)), rep(lm_L_fit$coef[2], nrow(sghgt_coefs)), rep(coef(complex_L_fit)[1], nrow(sghgt_coefs)), rep(coef(complex_L_fit)[2], nrow(sghgt_coefs)), rep(coef(complex_L_fit)[3], nrow(sghgt_coefs)))
		colnames(add_coefs)<-c('L-b (2FVT)', 'L-m (2FVT)', 'L-xg (2FVT Gaussian)', 'L-kg (2FVT Gaussian)', 'L-g (2FVT Gaussian)')
		sghgt_coefs<-cbind(sghgt_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(sghgt_coefs[,c(1,4)], xlim=c(0,max(sghgt_coefs[,1])), ylim=c(0,max(sghgt_coefs[,4])*1.5), main='Max. Segmented', xlab='Phytomer no.', ylab='L')
		points(sghgt_coefs[cooks_dist_index_L,1], sghgt_coefs[cooks_dist_index_L,4], pch=19, col='red', cex=0.55); abline(lm_L_fit);  abline(h=mean(sghgt_coefs[!cooks_dist_index_L,4]), lty=3); 
		text(summary(1:max(sghgt_coefs[,1]))[2]*0.5, summary(1:max(sghgt_coefs[,4]))[5]*1.5, paste('Linear adj. R2:', round(summary(lm_L_fit)$adj.r.squared,3)))
		text(summary(1:max(sghgt_coefs[,1]))[2]*0.5, summary(1:max(sghgt_coefs[,4]))[5]*1.25, paste('slope:', round(coef(lm_L_fit)[2],3)))
		x_fit=seq(0, max(complex_L$x), by=0.1); y_fit=coef(complex_L_fit)[3]*exp(-(x_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)); lines(x_fit, y_fit, lwd=2);
		legend('topleft', lty=c(1,3,1,0), lwd=c(1,1,2,0), pch=c(-1,-1,-1,19), col=c('black','black','black','red'), legend=c('Complete', 'Simplified', 'Gaussian', 'Outliers'))

		AIC_vals<-c()
		x_fit = seq(0, 30, by = 0.1); 
		
		for (FVT_round in c('Complex', 'Complete', 'Simple')){

			grand_model_fits<-c()

			#Replot the Segmented Heights using the parameters generated from the linear regressions (i.e. second-order FVT estimation)
			plot(c(0,30),c(0, max(combined_internode_lengths[,2])), col='white', xlab='Days', ylab='Ligular distance (cm)', main=paste(line, 'Segmented Height (2FVT)'))
			
			for (phyt_fit in 2:max(sghgt_coefs[,1])){

				unit_measures<-c()
				unit<-paste('_lig', phyt_fit, '_lig', phyt_fit+1, sep='')
				time_courses<-combined_internode_lengths[grepl(unit, rownames(combined_internode_lengths)),]

				#Ensure units don't advance beyond the limits of phytomers that exist within the genotype
				#If statement provides an escape should the unit list exceed the number of phytomers that exist within the background...
				if(nrow(time_courses)>0){
					#Fill in missing data for reps by setting early timepoints to zero when organs being modeled don't exist yet...
					for (rep in reps){
						matches<-!grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
						time_series<-time_courses[matches,]
						
						#print(paste('Unit',u,'Rep',rep))
						if(sum(matches)==0){
							#print(paste('Unit',u,'Rep',rep,': No hits'))
							next
						}else if(sum(matches)/length(days)<0.5){
							if (sum(matches)==1){
								#print(paste('Unit',u,'Rep',rep,': One hit'))
								if(4<(min(time_series[1])-3)){
									missing_days<-4:(min(time_series[1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(time_series)<-rep(paste('rep',rep,'_lig',phyt_fit,'_lig',phyt_fit+1, sep=''), nrow(time_series))

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}else{
								#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
								if(4<(min(time_series[,1])-3)){
									missing_days<-4:(min(time_series[,1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(time_series)<-rep(paste('rep',rep,'_lig',phyt_fit,'_lig',phyt_fit+1, sep=''), nrow(time_series))

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}
						}else if(sum(matches)/length(days)>=0.5){
							#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
							if(4<(min(time_series[,1])-3)){
								missing_days<-4:(min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(time_series)<-rep(paste('rep',rep,'_lig',phyt_fit,'_lig',phyt_fit+1, sep=''), nrow(time_series))

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)	
						}	
					}

					growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
					points(growth$age, growth$dist, cex=0.3, col=cols[phyt_fit-1])
				}

				if (length(growth$age)!=0){ #Only model glm's if the testing reps have measures for the current phytomer
					#Fit parameters for given phytomer given the second order FVT linear regression models...
					fitted_k <-  as.numeric(lm_k_fit$coefficients[1]+phyt_fit*lm_k_fit$coefficients[2])
					fitted_xo <- as.numeric((phyt_fit-lm_xo_fit$coefficients[1])/lm_xo_fit$coefficients[2])
					fitted_L <-  as.numeric(lm_L_fit$coefficients[1]+phyt_fit*lm_L_fit$coefficients[2])
					avg_k <-     as.numeric(mean(sghgt_coefs[!cooks_dist_index_k,3]))
					avg_L <-     as.numeric(mean(sghgt_coefs[!cooks_dist_index_L,4]))
					gauss_L <-   as.numeric(coef(complex_L_fit)[3]*exp(-(phyt_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)))

					#Extract NLS fitted parameters for comparison...
					if (length(sghgt_coefs[sghgt_coefs[,1] %in% phyt_fit,2])>0){
						ind_k <- sghgt_coefs[sghgt_coefs[,1] %in% phyt_fit,3]
						ind_xo <- sghgt_coefs[sghgt_coefs[,1] %in% phyt_fit,2]
						ind_L <- sghgt_coefs[sghgt_coefs[,1] %in% phyt_fit,4]
					}

					if (FVT_round=='Complex'){
						#Test various Generalized Linear Models using Second order parameters (Complete and Simplified) as well as Individual parameters derived from NLS...
						Sec.FVT.glm.gaussian.L<-    glm(growth$dist ~ rep(gauss_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						Sec.FVT.glm<-               glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-rep(fitted_k, length(growth$age))*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.simplified<-    glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-avg_k*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.oversimplified<-glm(growth$dist ~ rep(avg_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						if (length(sghgt_coefs[sghgt_coefs[,1] %in% phyt_fit,2])>0){
							Ind.par.glm<-           glm(growth$dist ~ rep(ind_L, length(growth$age))/(1+exp(-rep(ind_k, length(growth$age))*(growth$age-rep(ind_xo, length(growth$age))))), growth, family=gaussian)
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), round(extractAIC(Ind.par.glm)[2],0))
							dAICs<-AICs-min(AICs)
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), as.numeric(Ind.par.glm[10]))
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], logLik(Ind.par.glm)[1], AIC_weights))
						}else{
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), NA)
							dAICs<-AICs-min(na.omit(AICs))
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), NA)
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], NA, AIC_weights))
						}
						
						y_fit_gauss =  gauss_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_p         =  y_fit_gauss/max(y_fit_gauss); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_gauss, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};
						
						lines(x_fit, y_fit_gauss, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, lwd=2, legend=c('L gaussian and k simplified'))
					}else if(FVT_round=='Complete'){
						y_fit = fitted_L/(1+exp(-fitted_k*(x_fit-fitted_xo)));
						y_p         =  y_fit/max(y_fit); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};
						
						lines(x_fit, y_fit, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('Complete linear'))
					}else if(FVT_round=='Simple'){
						y_fit_simple = fitted_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_fit_oversimple = avg_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_p         =  y_fit_simple/max(y_fit_simple); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_simple, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=c(1,2); if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};						  
						
						lines(x_fit, y_fit_simple, col=cols[phyt_fit-1], lty=line_types[1], lwd=1+wdth_wgt)
						lines(x_fit, y_fit_oversimple, col=cols[phyt_fit-1], lty=line_types[2], lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('k simplified', 'k and L simplified'))
					}
				}
			}

			grand_model_fits<-cbind(x_fit, grand_model_fits)
			write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/', FVT_round, '_', line, '_seg_culm_hgt_fits_w_Wald_CI_weights.csv', sep=''), grand_model_fits, row.names=FALSE)
		}

		dev.off()

		colnames(AIC_vals)<-c('Phytomer', 'Null Dev.', 
			'2FVT (Gaussian L) Res. Dev.', '2FVT (Complete) Res. Dev.', '2FVT (Simplified) Res. Dev.', '2FVT (Over-simplified) Res. Dev.', 'Ind. Par. Res. Dev.', 
			'2FVT (Gaussian L) AIC', '2FVT (Complete) AIC', '2FVT (Simplified) AIC', '2FVT (Over-simplified) AIC', 'Ind. Par. AIC',
			'2FVT (Gaussian L) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC',
			'2FVT (Gaussian L) log Lik.', '2FVT (Complete) log Lik.', '2FVT (Simplified) log Lik.', '2FVT (Over-simplified) log Lik.', 'Ind. Par. log Lik.', 
			'2FVT (Gaussian L) weights', '2FVT (Complete) weights', '2FVT (Simplified) weights', '2FVT (Over-simplified) weights', 'Ind. Par. weights')

		AIC_vals

		sigmaAICs<-cbind(sum(na.omit(AIC_vals[,8])), sum(na.omit(AIC_vals[,9])),  sum(na.omit(AIC_vals[,10])),  sum(na.omit(AIC_vals[,11])),  sum(na.omit(AIC_vals[,12])))
		sigmaAICs<-cbind(line, 'seg.hgt', sigmaAICs, sigmaAICs-min(sigmaAICs))
		colnames(sigmaAICs)<-c('genotype', 'attribute', '2FVT (Gaussian) Sigma AIC', '2FVT (Complete) Sigma AIC', '2FVT (Simplified) Sigma AIC', '2FVT (Over-simplified) Sigma AIC', 'Ind. Par. Sigma AIC', '2FVT (Gaussian) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC')		

		Sec.FVT.AICs<-rbind(Sec.FVT.AICs, sigmaAICs)

		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_parameters/', line,'_seg_culm_hgt_2FVT_fitted_parameters.csv',sep=''), sghgt_coefs, row.names=FALSE)
		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_model_fits/', line,'_seg_culm_hgt_AICs.csv',sep=''), AIC_vals)


		#####
		# Subcomponent Block B: Leaf Tip-Ligule Distance Screening
		#####


		pdf(paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/leaf_dist_', line,'_2FVT.pdf',sep=''), height=12, width=10)

		if (verbose==TRUE){print('Initializing block B')}

		layout(matrix(c(1,1,1,7,7,7,6,6,6,5,5,5,2,2,3,3,4,4), 3, 6, byrow=TRUE))

		plot(c(0,30),c(0, max(combined_leaf_dist[,2])), col='white', xlab='Days', ylab='Leaf Tip-Ligule distance (cm)', main=paste(line, 'Leaf Distances (Independent Models)'))

		cols=colorRampPalette(c("darkgreen", "orange", "red"))(8)
		legend('topleft', pch=rep(19, 8), col=cols, legend=paste('Phytomer ', 2:9), cex=0.8)

		if(line=='A10'){xo_min=4; gap_days<-2; k_min=15; k_bounds=c(0.001,10000); kg=1}
		if(line=='RIL39'){xo_min=4; gap_days<-3; k_min=1; k_bounds=c(0.001, 10000); kg=1}
		if(line=='RIL110'){xo_min=4; gap_days<-2; k_min=1; k_bounds=c(0.001,10000); kg=1}
		if(line=='RIL159'){xo_min=4; gap_days<-3; k_min=1; k_bounds=c(0.001,10000); kg=1}
		if(line=='B100'){xo_min=4; gap_days<-2; k_min=1; k_bounds=c(0.001,10000); kg=1}

		lfdist_coefs<-c();

		for (u in 3:8){
			
			unit_measures<-c()
			unit<-paste('leaf', u, sep='')
			time_courses<-combined_leaf_dist[grepl(unit, rownames(combined_leaf_dist)),]

			if (nrow(time_courses)>0){ #Ensure units don't advance beyond the limits of phytomers that exist within the genotype

				for (rep in reps){
					matches<-grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
					time_series<-time_courses[matches,]
					
					#print(paste('Unit',u,'Rep',rep))
					if(sum(matches)==0){
						#print(paste('Unit',u,'Rep',rep,': No hits'))
						next
					}else if(sum(matches)/length(days)<0.5){
						if (sum(matches)==1){
							#print(paste('Unit',u,'Rep',rep,': One hit'))
							if(4<(min(time_series[1])-gap_days)){
								missing_days<-seq(4,min(time_series[1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',u,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)
						}else{
							#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
							if(4<(min(time_series[,1])-gap_days)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',u,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)
						}
					}else if(sum(matches)/length(days)>=0.5){
						#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
						if(4<(min(time_series[,1])-gap_days)){
							missing_days<-seq(4,min(time_series[,1])-gap_days)
							missing_dist<-c(rep(0, length(missing_days)))
						}else{
							missing_days<-4; missing_dist<-0;
						}
						missing_series<-cbind(missing_days, missing_dist)
						colnames(missing_series)<-colnames(time_series)
						rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',u,sep='')

						time_series<-rbind(missing_series, time_series)
						unit_measures<-rbind(unit_measures, time_series)	
					}
				}

				growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
				points(growth$age, growth$dist, cex=0.3, col=cols[u-1])

				#1/lm(growth$hgt[growth$hgt!=0]~growth$age[growth$hgt!=0])$coef[2]
				start_params=list(xo = xo_min, k = k_min, L = mean(growth$dist[growth$dist!=0]))
				lo_par_bound=list(xo = floor(xo_min-1), k = k_bounds[1], L = 0)
				hi_par_bound=list(xo = summary(growth$age[growth$dist!=0])[[5]]+4, k = k_bounds[2], L = max(growth$dist[growth$dist!=0])*1.5)

				print(paste(line, 'leaf distances'))
				print(rbind(start_params, lo_par_bound, hi_par_bound))

				nlm_fit<-nls(growth$dist ~ L/(1+exp(-k*(growth$age-xo))), data = growth, start=start_params, algorithm='port', 
					trace=T, lower=lo_par_bound, upper=hi_par_bound, control=list(maxiter=100))
				
				coefs = coef(nlm_fit)
				
				#Estimate a fitted curve based on coefficients of nlm
				x_fit = seq(0, 30, by = 0.1); y_fit = coefs[3]/(1+exp(-coefs[2]*(x_fit-coefs[1])));
				x_scale = x_fit[y_fit>=1][1]
				adj_x_fit = x_fit[y_fit>=1]-x_scale; adj_y_fit = y_fit[y_fit>=1];

				lines(x_fit, y_fit, col=cols[u-1])
				coefs = coef(nlm_fit); lfdist_coefs<-rbind(lfdist_coefs,c(u , coefs));

				xo_min<-lfdist_coefs[nrow(lfdist_coefs),2]+2

			}
		}

		colnames(lfdist_coefs)<-c('Phytomer', 'xo (ind)', 'k (ind)', 'L (ind)')

		#Use linear regression to estimate base k parameter for a given phytomer
		lm_k_fit<-lm(lfdist_coefs[,3]~lfdist_coefs[,1])
		cooks_dist_index_k<-cooks.distance(lm_k_fit)>3*mean(cooks.distance(lm_k_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_k)>0){
			lm_k_fit<-lm(lfdist_coefs[!cooks_dist_index_k,3]~lfdist_coefs[!cooks_dist_index_k,1])
		}

		add_coefs<-cbind(rep(mean(lfdist_coefs[!cooks_dist_index_k,3]), nrow(lfdist_coefs)), rep(lm_k_fit$coef[1], nrow(lfdist_coefs)), rep(lm_k_fit$coef[2], nrow(lfdist_coefs)))
		colnames(add_coefs)<-c('Simple k (2FVT)', 'k-b (2FVT)', 'k-m (2FVT)')
		lfdist_coefs<-cbind(lfdist_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(lfdist_coefs[,c(1,3)], xlim=c(0,max(lfdist_coefs[,1])), ylim=c(0,max(lfdist_coefs[,3])), main='Max. Phytomer Growth Rate', xlab='Phytomer no.', ylab='k')
		points(lfdist_coefs[cooks_dist_index_k,1], lfdist_coefs[cooks_dist_index_k,3], pch=19, col='red', cex=0.55); abline(lm_k_fit); abline(h=mean(lfdist_coefs[!cooks_dist_index_k,3]), lty=3)
		text(summary(1:max(lfdist_coefs[,1]))[2]*.5, summary(1:max(lfdist_coefs[,3]))[5]*0.9, paste('adj. R2:', round(summary(lm_k_fit)$adj.r.squared,3)))
		text(summary(1:max(lfdist_coefs[,1]))[2]*.5, summary(1:max(lfdist_coefs[,3]))[5]*0.8, paste('slope:', round(coef(lm_k_fit)[2],3)))
		legend('topleft', lty=c(1,3,0), pch=c(-1,-1,19), col=c('black','black','red'), legend=c('Complete', 'Simplified', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_xo_fit<-lm(lfdist_coefs[,1]~lfdist_coefs[,2])
		cooks_dist_index_xo<-cooks.distance(lm_xo_fit)>3*mean(cooks.distance(lm_xo_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_xo)>0){
			lm_xo_fit<-lm(lfdist_coefs[!cooks_dist_index_xo,1]~lfdist_coefs[!cooks_dist_index_xo,2])
		}

		add_coefs<-cbind(rep(lm_xo_fit$coef[1], nrow(lfdist_coefs)), rep(lm_xo_fit$coef[2], nrow(lfdist_coefs)))
		colnames(add_coefs)<-c('xo-b (2FVT)', 'xo-m (2FVT)')
		lfdist_coefs<-cbind(lfdist_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(lfdist_coefs[,c(2,1)], xlim=c(0,max(lfdist_coefs[,2])), ylim=c(0,max(lfdist_coefs[,1])), main='Rate of Phytomer Deposition', xlab='xo (days)', ylab='Phytomer no.')
		points(lfdist_coefs[cooks_dist_index_xo,2], lfdist_coefs[cooks_dist_index_xo,1], pch=19, col='red', cex=0.55); abline(lm_xo_fit);
		text(summary(1:max(lfdist_coefs[,2]))[2]*0.5, summary(1:max(lfdist_coefs[,1]))[5], paste('adj. R2:', round(summary(lm_xo_fit)$adj.r.squared,3)))
		text(summary(1:max(lfdist_coefs[,2]))[2]*0.5, summary(1:max(lfdist_coefs[,1]))[5]*0.9, paste('slope:', round(coef(lm_xo_fit)[2],3)))
		legend('topleft', lty=c(1,0), pch=c(-1,19), col=c('black','red'), legend=c('Complete', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_L_fit<-lm(lfdist_coefs[,4]~lfdist_coefs[,1])
		cooks_dist_index_L<-cooks.distance(lm_L_fit)>3*mean(cooks.distance(lm_L_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		complex_L=list(y = as.numeric(lfdist_coefs[,4]), x = as.numeric(lfdist_coefs[,1]))

		gauss_params=list(xg = lfdist_coefs[which.max(lfdist_coefs[,4]),1], kg = kg, Lg = mean(lfdist_coefs[,4]))
		gauss_lo_par=list(xg = lfdist_coefs[which.max(lfdist_coefs[,4]),1]-2, kg = -1000, Lg = 0)
		gauss_hi_par=list(xg = lfdist_coefs[which.max(lfdist_coefs[,4]),1]+2, kg = 1000, Lg = max(lfdist_coefs[,4])*1.5)

		complex_L_fit<-nls(complex_L$y ~ Lg*exp(-(complex_L$x-xg)^2/(2*kg^2)), data = complex_L, start=gauss_params, algorithm='port', 
					trace=T, lower=gauss_lo_par, upper=gauss_hi_par)

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_L)>0){
			lm_L_fit<-lm(lfdist_coefs[!cooks_dist_index_L,4]~lfdist_coefs[!cooks_dist_index_L,1])
		}

		add_coefs<-cbind(rep(lm_L_fit$coef[1], nrow(lfdist_coefs)), rep(lm_L_fit$coef[2], nrow(lfdist_coefs)), rep(coef(complex_L_fit)[1], nrow(lfdist_coefs)), rep(coef(complex_L_fit)[2], nrow(lfdist_coefs)), rep(coef(complex_L_fit)[3], nrow(lfdist_coefs)))
		colnames(add_coefs)<-c('L-b (2FVT)', 'L-m (2FVT)', 'L-xg (2FVT Gaussian)', 'L-kg (2FVT Gaussian)', 'L-g (2FVT Gaussian)')
		lfdist_coefs<-cbind(lfdist_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(lfdist_coefs[,c(1,4)], xlim=c(0,max(lfdist_coefs[,1])), ylim=c(0,max(lfdist_coefs[,4])*1.5), main='Max. Leaf Distances', xlab='Phytomer no.', ylab='L')
		points(lfdist_coefs[cooks_dist_index_L,1], lfdist_coefs[cooks_dist_index_L,4], pch=19, col='red', cex=0.55); abline(lm_L_fit); abline(h=mean(lfdist_coefs[!cooks_dist_index_L,4]), lty=3)
		text(summary(1:max(lfdist_coefs[,1]))[2]*.5, max(lfdist_coefs[,4])*1.1, paste('Linear adj. R2:', round(summary(lm_L_fit)$adj.r.squared,3)))
		text(summary(1:max(lfdist_coefs[,1]))[2]*.5, max(lfdist_coefs[,4]), paste('slope:', round(coef(lm_L_fit)[2],3)))
		x_fit=seq(0, max(complex_L$x), by=0.1); y_fit=coef(complex_L_fit)[3]*exp(-(x_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)); lines(x_fit, y_fit, lwd=2);
		legend('topleft', lty=c(1,3,1,0), lwd=c(1,1,2,0), pch=c(-1,-1,-1,19), col=c('black','black','black','red'), legend=c('Complete', 'Simplified', 'Gaussian', 'Outliers'))

		AIC_vals<-c()
		x_fit = seq(0, 30, by = 0.1); 
		
		for (FVT_round in c('Complex', 'Complete', 'Simple')){

			grand_model_fits<-c()

			#Replot the Segmented Heights using the parameters generated from the linear regressions (i.e. second-order FVT estimation)
			plot(c(0,30),c(0, max(combined_leaf_dist[,2])), col='white', xlab='Days', ylab='Leaf Tip-Ligule distance (cm)', main=paste(line, 'Leaf Tip-Ligule Distances (2nd Order FVTs)'))

			for (phyt_fit in 2:max(lfdist_coefs[,1])){

				unit_measures<-c()
				unit<-paste('leaf', phyt_fit, sep='')
				time_courses<-combined_leaf_dist[grepl(unit, rownames(combined_leaf_dist)),]

				#Ensure units don't advance beyond the limits of phytomers that exist within the genotype
				#If statement provides an escape should the unit list exceed the number of phytomers that exist within the background...
				if(nrow(time_courses)>0){
					#Fill in missing data for reps by setting early timepoints to zero when organs being modeled don't exist yet...
					for (rep in reps){
						matches<-!grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
						time_series<-time_courses[matches,]
						
						#print(paste('Unit',u,'Rep',rep))
						if(sum(matches)==0){
							#print(paste('Unit',u,'Rep',rep,': No hits'))
							next
						}else if(sum(matches)/length(days)<0.5){
							if (sum(matches)==1){
								#print(paste('Unit',u,'Rep',rep,': One hit'))
								if(4<(min(time_series[1])-3)){
									missing_days<-seq(4,min(time_series[,1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',phyt_fit,sep='')

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}else{
								#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
								if(4<(min(time_series[,1])-3)){
									missing_days<-seq(4,min(time_series[,1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',phyt_fit,sep='')

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}
						}else if(sum(matches)/length(days)>=0.5){
							#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
							if(4<(min(time_series[,1])-3)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_leaf',phyt_fit,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)	
						}	
					}

					growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
					points(growth$age, growth$dist, cex=0.3, col=cols[phyt_fit-1])
				}



				if (length(growth$age)!=0){ #Only model glm's if the testing reps have measures for the current phytomer
					#Fit parameters for given phytomer given the second order FVT linear regression models...
					fitted_k <-  as.numeric(lm_k_fit$coefficients[1]+phyt_fit*lm_k_fit$coefficients[2])
					fitted_xo <- as.numeric((phyt_fit-lm_xo_fit$coefficients[1])/lm_xo_fit$coefficients[2])
					fitted_L <-  as.numeric(lm_L_fit$coefficients[1]+phyt_fit*lm_L_fit$coefficients[2])
					avg_k <-     as.numeric(mean(lfdist_coefs[!cooks_dist_index_k,3]))
					avg_L <-     as.numeric(mean(lfdist_coefs[!cooks_dist_index_L,4]))
					gauss_L <-   as.numeric(coef(complex_L_fit)[3]*exp(-(phyt_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)))

					#Extract NLS fitted parameters for comparison...
					if (length(lfdist_coefs[lfdist_coefs[,1] %in% phyt_fit,2])>0){
						ind_k <- lfdist_coefs[lfdist_coefs[,1] %in% phyt_fit,3]
						ind_xo <- lfdist_coefs[lfdist_coefs[,1] %in% phyt_fit,2]
						ind_L <- lfdist_coefs[lfdist_coefs[,1] %in% phyt_fit,4]
					}

					
					if (FVT_round=='Complex'){
						#Test various Generalized Linear Models using Second order parameters (Complete and Simplified) as well as Individual parameters derived from NLS...
						Sec.FVT.glm.gaussian.L<-    glm(growth$dist ~ rep(gauss_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						Sec.FVT.glm<-               glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-rep(fitted_k, length(growth$age))*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.simplified<-    glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-avg_k*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.oversimplified<-glm(growth$dist ~ rep(avg_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						if (length(lfdist_coefs[lfdist_coefs[,1] %in% phyt_fit,2])>0){
							Ind.par.glm<-           glm(growth$dist ~ rep(ind_L, length(growth$age))/(1+exp(-rep(ind_k, length(growth$age))*(growth$age-rep(ind_xo, length(growth$age))))), growth, family=gaussian)
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), round(extractAIC(Ind.par.glm)[2],0))
							dAICs<-AICs-min(AICs)
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), as.numeric(Ind.par.glm[10]))
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], logLik(Ind.par.glm)[1], AIC_weights))
						}else{
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), NA)
							dAICs<-AICs-min(na.omit(AICs))
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), NA)
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], NA, AIC_weights))
						}
						
						y_fit_gauss =  gauss_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_p         =  y_fit_gauss/max(y_fit_gauss); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_gauss, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};
						
						lines(x_fit, y_fit_gauss, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, lwd=2, legend=c('L gaussian and k simplified'))
					}else if(FVT_round=='Complete'){
						y_fit = fitted_L/(1+exp(-fitted_k*(x_fit-fitted_xo)));
						y_p         =  y_fit/max(y_fit); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};
						
						lines(x_fit, y_fit, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('Complete linear'))
					}else if(FVT_round=='Simple'){
						y_fit_simple = fitted_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_fit_oversimple = avg_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_p         =  y_fit_simple/max(y_fit_simple); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_simple, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)		
						
						line_types=c(1,2); if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};						  
						
						lines(x_fit, y_fit_simple, col=cols[phyt_fit-1], lty=line_types[1], lwd=1+wdth_wgt)
						lines(x_fit, y_fit_oversimple, col=cols[phyt_fit-1], lty=line_types[2], lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('k simplified', 'k and L simplified'))				  
					}
				}
			}

			grand_model_fits<-cbind(x_fit, grand_model_fits)
			write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/', FVT_round, '_', line, '_leaf_dist_fits_w_Wald_CI_weights.csv', sep=''), grand_model_fits, row.names=FALSE)

		}

		dev.off()

		colnames(AIC_vals)<-c('Phytomer', 'Null Dev.', 
			'2FVT (Gaussian L) Res. Dev.', '2FVT (Complete) Res. Dev.', '2FVT (Simplified) Res. Dev.', '2FVT (Over-simplified) Res. Dev.', 'Ind. Par. Res. Dev.', 
			'2FVT (Gaussian L) AIC', '2FVT (Complete) AIC', '2FVT (Simplified) AIC', '2FVT (Over-simplified) AIC', 'Ind. Par. AIC',
			'2FVT (Gaussian L) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC',
			'2FVT (Gaussian L) log Lik.', '2FVT (Complete) log Lik.', '2FVT (Simplified) log Lik.', '2FVT (Over-simplified) log Lik.', 'Ind. Par. log Lik.', 
			'2FVT (Gaussian L) weights', '2FVT (Complete) weights', '2FVT (Simplified) weights', '2FVT (Over-simplified) weights', 'Ind. Par. weights')

		AIC_vals	

		sigmaAICs<-cbind(sum(na.omit(AIC_vals[,8])), sum(na.omit(AIC_vals[,9])),  sum(na.omit(AIC_vals[,10])),  sum(na.omit(AIC_vals[,11])),  sum(na.omit(AIC_vals[,12])))
		sigmaAICs<-cbind(line, 'lf.dist', sigmaAICs, sigmaAICs-min(sigmaAICs))
		colnames(sigmaAICs)<-c('genotype', 'attribute', '2FVT (Gaussian) Sigma AIC', '2FVT (Complete) Sigma AIC', '2FVT (Simplified) Sigma AIC', '2FVT (Over-simplified) Sigma AIC', 'Ind. Par. Sigma AIC', '2FVT (Gaussian) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC')		

		Sec.FVT.AICs<-rbind(Sec.FVT.AICs, sigmaAICs)

		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_parameters/', line,'_leaf_dist_2FVT_fitted_parameters.csv',sep=''), lfdist_coefs, row.names=FALSE)
		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_model_fits/', line,'_leaf_dist_AICs.csv',sep=''), AIC_vals)

		#####
		# Subcomponent Block C: Leaf Tip-Ligule Angle Screening
		#####

		pdf(paste('~/Desktop/2FVT_outputs/2FVT_model_graphs/leaf_angle_', line,'_2FVT.pdf',sep=''), height=12, width=10)

		if (verbose==TRUE){print('Initializing block C')}

		layout(matrix(c(1,1,1,7,7,7,6,6,6,5,5,5,2,2,3,3,4,4), 3, 6, byrow=TRUE))

		plot(c(0,30),c(0, 300), col='white', xlab='Days', ylab='Leaf Angles (\u00B0)', main=paste(line, 'Leaf Tip-Ligule Angles (Independent Models)'))

		cols=colorRampPalette(c("darkgreen", "orange", "red"))(8)
		legend('topleft', pch=rep(19, 8), col=cols, legend=paste('Phytomer ', 2:9), cex=0.8)

		if(line=='A10'){xo_min=3; gap_days<-1; k_min=2.5; k_bounds=c(0.001,100)}
		if(line=='RIL39'){xo_min=3; gap_days<-2; k_min=0.5; k_bounds=c(0.01, 100)}
		if(line=='RIL110'){xo_min=4; gap_days<-3; k_min=1; k_bounds=c(0.01,50)}
		if(line=='RIL159'){xo_min=6; gap_days<-3; k_min=16; k_bounds=c(0.0001,50)}
		if(line=='B100'){xo_min=4; gap_days<-2; k_min=1; k_bounds=c(0.01,50)}

		angle_coefs<-c();

		for (u in 3:8){
			
			unit_measures<-c()
			unit<-paste('_lig', u, sep='')
			time_courses<-combined_ligule_angles[grepl(unit, rownames(combined_ligule_angles)),]

			if (nrow(time_courses)>0){ #Ensure units don't advance beyond the limits of phytomers that exist within the genotype

				for (rep in reps){
					matches<-grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
					time_series<-time_courses[matches,]
					
					#print(paste('Unit',u,'Rep',rep))
					if(sum(matches)==0){
						#print(paste('Unit',u,'Rep',rep,': No hits'))
						next
					}else if(sum(matches)/length(days)<0.5){
						if (sum(matches)==1){
							#print(paste('Unit',u,'Rep',rep,': One hit'))
							if(4<(min(time_series[1])-3)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',u,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)
						}else{
							#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
							if(4<(min(time_series[,1])-3)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',u,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)
						}
					}else if(sum(matches)/length(days)>=0.5){
						#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
						if(4<(min(time_series[,1])-3)){
							missing_days<-seq(4,min(time_series[,1])-gap_days)
							missing_dist<-c(rep(0, length(missing_days)))
						}else{
							missing_days<-4; missing_dist<-0;
						}
						missing_series<-cbind(missing_days, missing_dist)
						colnames(missing_series)<-colnames(time_series)
						rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',u,sep='')

						time_series<-rbind(missing_series, time_series)
						unit_measures<-rbind(unit_measures, time_series)	
					}	
				}

				growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
				
				#1/lm(growth$hgt[growth$hgt!=0]~growth$age[growth$hgt!=0])$coef[2]
				start_params=list(xo = xo_min, k = k_min, L = mean(growth$dist[growth$dist!=0])/2)
				lo_par_bound=list(xo = floor(xo_min-1), k = k_bounds[1], L = 5)
				hi_par_bound=list(xo = summary(growth$age[growth$dist!=0])[[5]]+7, k = k_bounds[2], L = max(growth$dist[growth$dist!=0]))

				print(paste(line, 'leaf angles'))
				print(rbind(start_params, lo_par_bound, hi_par_bound))

				nlm_fit<-nls(growth$dist ~ L/(1+exp(-k*(growth$age-xo))), data = growth, start=start_params, algorithm='port', 
					trace=T, lower=lo_par_bound, upper=hi_par_bound)
				
				coefs = coef(nlm_fit)
				
				#Estimate a fitted curve based on coefficients of nlm
				x_fit = seq(0, 30, by = 0.1); y_fit = coefs[3]/(1+exp(-coefs[2]*(x_fit-coefs[1])));
				x_scale = x_fit[y_fit>=1][1]
				adj_x_fit = x_fit[y_fit>=1]-x_scale; adj_y_fit = y_fit[y_fit>=1];

				points(growth$age, growth$dist, cex=0.3, col=cols[u-1])
				lines(x_fit, y_fit, col=cols[u-1])
				coefs = coef(nlm_fit); angle_coefs<-rbind(angle_coefs,c(u , coefs));

				xo_min<-angle_coefs[nrow(angle_coefs),2]+1

			}
		}

		colnames(angle_coefs)<-c('Phytomer', 'xo (ind)', 'k (ind)', 'L (ind)')

		#Use linear regression to estimate base k parameter for a given phytomer
		lm_k_fit<-lm(angle_coefs[,3]~angle_coefs[,1])
		cooks_dist_index_k<-cooks.distance(lm_k_fit)>4*mean(cooks.distance(lm_k_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_k)>0){
			lm_k_fit<-lm(angle_coefs[!cooks_dist_index_k,3]~angle_coefs[!cooks_dist_index_k,1])
		}

		add_coefs<-cbind(rep(mean(angle_coefs[!cooks_dist_index_k,3]), nrow(angle_coefs)), rep(lm_k_fit$coef[1], nrow(angle_coefs)), rep(lm_k_fit$coef[2], nrow(angle_coefs)))
		colnames(add_coefs)<-c('Simple k (2FVT)', 'k-b (2FVT)', 'k-m (2FVT)')
		angle_coefs<-cbind(angle_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(angle_coefs[,c(1,3)], xlim=c(0,max(angle_coefs[,1])), ylim=c(0,max(angle_coefs[,3])), main='Max. Leaf Angle Deflection Rate', xlab='Phytomer no.', ylab='k')
		points(angle_coefs[cooks_dist_index_k,1], angle_coefs[cooks_dist_index_k,3], pch=19, col='red', cex=0.55); abline(lm_k_fit); abline(h=mean(angle_coefs[!cooks_dist_index_k,3]), lty=3)
		text(summary(1:max(angle_coefs[,1]))[2]*0.5, max(angle_coefs[,3])*0.8, paste('adj. R2:', round(summary(lm_k_fit)$adj.r.squared,3)))
		text(summary(1:max(angle_coefs[,1]))[2]*0.5, max(angle_coefs[,3])*0.7, paste('slope:', round(coef(lm_k_fit)[2],3)))
		legend('topleft', lty=c(1,3,0), pch=c(-1,-1,19), col=c('black','black','red'), legend=c('Complete', 'Simplified', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_xo_fit<-lm(angle_coefs[,1]~angle_coefs[,2])
		cooks_dist_index_xo<-cooks.distance(lm_xo_fit)>4*mean(cooks.distance(lm_xo_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_xo)>0){
			lm_xo_fit<-lm(angle_coefs[!cooks_dist_index_xo,1]~angle_coefs[!cooks_dist_index_xo,2])
		}

		add_coefs<-cbind(rep(lm_xo_fit$coef[1], nrow(angle_coefs)), rep(lm_xo_fit$coef[2], nrow(angle_coefs)))
		colnames(add_coefs)<-c('xo-b (2FVT)', 'xo-m (2FVT)')
		angle_coefs<-cbind(angle_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(angle_coefs[,c(2,1)], xlim=c(0,max(angle_coefs[,2])), ylim=c(0,max(angle_coefs[,1])), main='Rate of Leaf Deflection', xlab='xo (days)', ylab='Phytomer no.')
		points(angle_coefs[cooks_dist_index_xo,2], angle_coefs[cooks_dist_index_xo,1], pch=19, col='red', cex=0.55); abline(lm_xo_fit);
		text(summary(1:max(angle_coefs[,2]))[2]*0.5, summary(1:max(angle_coefs[,1]))[5], paste('adj. R2:', round(summary(lm_xo_fit)$adj.r.squared,3)))
		text(summary(1:max(angle_coefs[,2]))[2]*0.5, summary(1:max(angle_coefs[,1]))[5]*.9, paste('slope:', round(coef(lm_xo_fit)[2],3)))
		legend('topleft', lty=c(1,0), pch=c(-1,19), col=c('black','red'), legend=c('Complete', 'Outliers'))

		#Use linear regression to estimate base xo parameter for a given phytomer
		lm_L_fit<-lm(angle_coefs[,4]~angle_coefs[,1])
		cooks_dist_index_L<-cooks.distance(lm_L_fit)>4*mean(cooks.distance(lm_L_fit)) #Test for outliers that exceed 4 times the mean Cook's Distance

		complex_L=list(y = as.numeric(angle_coefs[,4]), x = as.numeric(angle_coefs[,1]))

		gauss_params=list(xg = angle_coefs[which.max(angle_coefs[,4]),1], kg = 5, Lg = mean(angle_coefs[,4]))
		gauss_lo_par=list(xg = angle_coefs[which.max(angle_coefs[,4]),1]-2, kg = -1000, Lg = 0)
		gauss_hi_par=list(xg = angle_coefs[which.max(angle_coefs[,4]),1]+2, kg = 1000, Lg = max(angle_coefs[,4])*1.5)

		complex_L_fit<-nls(complex_L$y ~ Lg*exp(-(complex_L$x-xg)^2/(2*kg^2)), data = complex_L, start=gauss_params, algorithm='port', 
					trace=T, lower=gauss_lo_par, upper=gauss_hi_par)

		#If outliers detected based on Cook's distance, scrub from list and rerun regression...
		if(sum(cooks_dist_index_L)>0){
			lm_L_fit<-lm(angle_coefs[!cooks_dist_index_L,4]~angle_coefs[!cooks_dist_index_L,1])
		}

		add_coefs<-cbind(rep(lm_L_fit$coef[1], nrow(angle_coefs)), rep(lm_L_fit$coef[2], nrow(angle_coefs)), rep(coef(complex_L_fit)[1], nrow(angle_coefs)), rep(coef(complex_L_fit)[2], nrow(angle_coefs)), rep(coef(complex_L_fit)[3], nrow(angle_coefs)))
		colnames(add_coefs)<-c('L-b (2FVT)', 'L-m (2FVT)', 'L-xg (2FVT Gaussian)', 'L-kg (2FVT Gaussian)', 'L-g (2FVT Gaussian)')
		angle_coefs<-cbind(angle_coefs, add_coefs)

		#Plot linear regression of coefficients to visualize fit against the parameters used to estimate model
		plot(angle_coefs[,c(1,4)], xlim=c(0,max(angle_coefs[,1])), ylim=c(0,max(angle_coefs[,4])*1.5), main='Max. Leaf Deflection Angles', xlab='Phytomer no.', ylab='L')
		points(angle_coefs[cooks_dist_index_L,1], angle_coefs[cooks_dist_index_L,4], pch=19, col='red', cex=0.55); abline(lm_L_fit); abline(h=mean(angle_coefs[!cooks_dist_index_L,4]), lty=3)
		text(summary(1:max(angle_coefs[,1]))[5], max(angle_coefs[,4])*1.1, paste('Linear adj. R2:', round(summary(lm_L_fit)$adj.r.squared,3)))
		text(summary(1:max(angle_coefs[,1]))[5], max(angle_coefs[,4]), paste('slope:', round(coef(lm_L_fit)[2],3)))
		x_fit=seq(0, max(complex_L$x), by=0.1); y_fit=coef(complex_L_fit)[3]*exp(-(x_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)); lines(x_fit, y_fit, lwd=2);
		legend('topright', lty=c(1,2,1,0), lwd=c(1,1,2,0), pch=c(-1,-1,-1,19), col=c('black','black','black','red'), legend=c('Complete', 'Simplified', 'Gaussian', 'Outliers'))

		AIC_vals<-c()
		x_fit = seq(0, 30, by = 0.1); 
		
		for (FVT_round in c('Complex', 'Complete', 'Simple')){

			grand_model_fits<-c()

			#Replot the Segmented Heights using the parameters generated from the linear regressions (i.e. second-order FVT estimation)
			plot(c(0,30),c(0, 300), col='white', xlab='Days', ylab='Leaf Angles (\u00B0)', main=paste(line, 'Leaf Tip-Ligule Angles (2nd Order FVTs)'))

			for (phyt_fit in 2:max(angle_coefs[,1])){

				unit_measures<-c()
				unit<-paste('_lig', phyt_fit, sep='')
				time_courses<-combined_ligule_angles[grepl(unit, rownames(combined_ligule_angles)),]

				#Ensure units don't advance beyond the limits of phytomers that exist within the genotype
				#If statement provides an escape should the unit list exceed the number of phytomers that exist within the background...
				if(nrow(time_courses)>0){
					#Fill in missing data for reps by setting early timepoints to zero when organs being modeled don't exist yet...
					for (rep in reps){
						matches<-!grepl(paste('rep', rep, '_', sep=''), rownames(time_courses))
						time_series<-time_courses[matches,]
						
						#print(paste('Unit',u,'Rep',rep))
						if(sum(matches)==0){
							#print(paste('Unit',u,'Rep',rep,': No hits'))
							next
						}else if(sum(matches)/length(days)<0.5){
							if (sum(matches)==1){
								#print(paste('Unit',u,'Rep',rep,': One hit'))
								if(4<(min(time_series[1])-3)){
									missing_days<-seq(4,min(time_series[,1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',phyt_fit,sep='')

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}else{
								#print(paste('Unit',u,'Rep',rep,': Few (<50%) hits'))
								if(4<(min(time_series[,1])-3)){
									missing_days<-seq(4,min(time_series[,1])-gap_days)
									missing_dist<-c(rep(0, length(missing_days)))
								}else{
									missing_days<-4; missing_dist<-0;
								}
								missing_series<-cbind(missing_days, missing_dist)
								colnames(missing_series)<-colnames(time_series)
								rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',phyt_fit,sep='')

								time_series<-rbind(missing_series, time_series)
								unit_measures<-rbind(unit_measures, time_series)
							}
						}else if(sum(matches)/length(days)>=0.5){
							#print(paste('Unit',u,'Rep',rep,': Many (>=50%) hits'))
							if(4<(min(time_series[,1])-3)){
								missing_days<-seq(4,min(time_series[,1])-gap_days)
								missing_dist<-c(rep(0, length(missing_days)))
							}else{
								missing_days<-4; missing_dist<-0;
							}
							missing_series<-cbind(missing_days, missing_dist)
							colnames(missing_series)<-colnames(time_series)
							rownames(missing_series)<-paste('rep',rep,'_d',missing_days,'_lig',phyt_fit,sep='')

							time_series<-rbind(missing_series, time_series)
							unit_measures<-rbind(unit_measures, time_series)	
						}	
					}

					growth=list(dist = as.numeric(unit_measures[,2]), age = as.numeric(unit_measures[,1]))
					points(growth$age, growth$dist, cex=0.3, col=cols[phyt_fit-1])
				}


				if (length(growth$age)!=0){ #Only model glm's if the testing reps have measures for the current phytomer
					#Fit parameters for given phytomer given the second order FVT linear regression models...
					fitted_k <-  as.numeric(lm_k_fit$coefficients[1]+phyt_fit*lm_k_fit$coefficients[2])
					fitted_xo <- as.numeric((phyt_fit-lm_xo_fit$coefficients[1])/lm_xo_fit$coefficients[2])
					fitted_L <-  as.numeric(lm_L_fit$coefficients[1]+phyt_fit*lm_L_fit$coefficients[2])
					avg_k <-  	 as.numeric(mean(angle_coefs[!cooks_dist_index_k,3]))
					avg_L <-     as.numeric(mean(angle_coefs[!cooks_dist_index_L,4]))
					gauss_L <-   as.numeric(coef(complex_L_fit)[3]*exp(-(phyt_fit-coef(complex_L_fit)[1])^2/(2*coef(complex_L_fit)[2]^2)))

					#Extract NLS fitted parameters for comparison...
					if (length(angle_coefs[angle_coefs[,1] %in% 2,2])>0){
						ind_k <- angle_coefs[angle_coefs[,1] %in% phyt_fit,3]
						ind_xo <- angle_coefs[angle_coefs[,1] %in% phyt_fit,2]
						ind_L <- angle_coefs[angle_coefs[,1] %in% phyt_fit,4]
					}

					if (FVT_round=='Complex'){
						#Test various Generalized Linear Models using Second order parameters (Complete and Simplified) as well as Individual parameters derived from NLS...
						Sec.FVT.glm.gaussian.L<-    glm(growth$dist ~ rep(gauss_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						Sec.FVT.glm<-               glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-rep(fitted_k, length(growth$age))*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.simplified<-    glm(growth$dist ~ rep(fitted_L, length(growth$age))/(1+exp(-avg_k*(growth$age-rep(fitted_xo, length(growth$age))))), growth, family=gaussian)
						Sec.FVT.glm.oversimplified<-glm(growth$dist ~ rep(avg_L, length(growth$age))/(1+exp(-rep(avg_k*(growth$age-rep(fitted_xo, length(growth$age)))))), growth, family=gaussian)
						if (length(angle_coefs[angle_coefs[,1] %in% phyt_fit,2])>0){
							Ind.par.glm<-           glm(growth$dist ~ rep(ind_L, length(growth$age))/(1+exp(-rep(ind_k, length(growth$age))*(growth$age-rep(ind_xo, length(growth$age))))), growth, family=gaussian)
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), round(extractAIC(Ind.par.glm)[2],0))
							dAICs<-AICs-min(AICs)
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), as.numeric(Ind.par.glm[10]))
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], logLik(Ind.par.glm)[1], AIC_weights))
						}else{
							AICs<-c(round(extractAIC(Sec.FVT.glm.gaussian.L)[2],0), round(extractAIC(Sec.FVT.glm)[2],0), round(extractAIC(Sec.FVT.glm.simplified)[2],0), round(extractAIC(Sec.FVT.glm.oversimplified)[2],0), NA)
							dAICs<-AICs-min(na.omit(AICs))
							AIC_weights<-exp(- 0.5 * dAICs) / sum(exp(- 0.5 * dAICs))
							res.dev<-c(as.numeric(Sec.FVT.glm.gaussian.L[10]), as.numeric(Sec.FVT.glm[10]), as.numeric(Sec.FVT.glm.simplified[10]), as.numeric(Sec.FVT.glm.oversimplified[10]), NA)
							nul.dev<-as.numeric(Sec.FVT.glm[12])
							AIC_vals<-rbind(AIC_vals, c(phyt_fit, nul.dev, res.dev, AICs, dAICs,logLik(Sec.FVT.glm.gaussian.L)[1], logLik(Sec.FVT.glm)[1], logLik(Sec.FVT.glm.simplified)[1], logLik(Sec.FVT.glm.oversimplified)[1], NA, AIC_weights))
						}
						
						y_fit_gauss =  gauss_L/(1+exp(-avg_k*(x_fit-fitted_xo)));  
						y_p         =  y_fit_gauss/max(y_fit_gauss); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_gauss, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};
						
						lines(x_fit, y_fit_gauss, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, lwd=2, legend=c('L gaussian and k simplified'))
					}else if(FVT_round=='Complete'){
						y_fit = fitted_L/(1+exp(-fitted_k*(x_fit-fitted_xo)));
						y_p         =  y_fit/max(y_fit); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)
						
						line_types=1; if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};						
						
						lines(x_fit, y_fit, col=cols[phyt_fit-1], lty=line_types, lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('Complete linear'))
					}else if(FVT_round=='Simple'){
						y_fit_simple = fitted_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_fit_oversimple = avg_L/(1+exp(-avg_k*(x_fit-fitted_xo)));
						y_p         =  y_fit_simple/max(y_fit_simple); 
						CI_weight   =  sqrt((y_p*(1-y_p))/(length(growth$age)+4)); 
						model_fits  =  cbind(y_fit_simple, CI_weight)
						colnames(model_fits)=c(paste('y_phyt_', phyt_fit, sep=''), 'CI_weight')
						grand_model_fits<-cbind(grand_model_fits, model_fits)						  
						
						line_types=c(1,2); if(phyt_fit==2){wdth_wgt=2}else{wdth_wgt=0};						  
						
						lines(x_fit, y_fit_simple, col=cols[phyt_fit-1], lty=line_types[1], lwd=1+wdth_wgt)
						lines(x_fit, y_fit_oversimple, col=cols[phyt_fit-1], lty=line_types[2], lwd=1+wdth_wgt)
						legend('topleft', lty=line_types, legend=c('k simplified', 'k and L simplified'))
					}
				}
			}

			grand_model_fits<-cbind(x_fit, grand_model_fits)
			write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_fits_w_Wald_CIs/', FVT_round, '_', line, '_leaf_angle_fits_w_Wald_CI_weights.csv', sep=''), grand_model_fits, row.names=FALSE)

		}

		dev.off()

		colnames(AIC_vals)<-c('Phytomer', 'Null Dev.', 
			'2FVT (Gaussian L) Res. Dev.', '2FVT (Complete) Res. Dev.', '2FVT (Simplified) Res. Dev.', '2FVT (Over-simplified) Res. Dev.', 'Ind. Par. Res. Dev.', 
			'2FVT (Gaussian L) AIC', '2FVT (Complete) AIC', '2FVT (Simplified) AIC', '2FVT (Over-simplified) AIC', 'Ind. Par. AIC',
			'2FVT (Gaussian L) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC',
			'2FVT (Gaussian L) log Lik.', '2FVT (Complete) log Lik.', '2FVT (Simplified) log Lik.', '2FVT (Over-simplified) log Lik.', 'Ind. Par. log Lik.', 
			'2FVT (Gaussian L) weights', '2FVT (Complete) weights', '2FVT (Simplified) weights', '2FVT (Over-simplified) weights', 'Ind. Par. weights')

		AIC_vals

		sigmaAICs<-cbind(sum(na.omit(AIC_vals[,8])), sum(na.omit(AIC_vals[,9])),  sum(na.omit(AIC_vals[,10])),  sum(na.omit(AIC_vals[,11])),  sum(na.omit(AIC_vals[,12])))
		sigmaAICs<-cbind(line, 'lf.angle', sigmaAICs, sigmaAICs-min(sigmaAICs))
		colnames(sigmaAICs)<-c('genotype', 'attribute', '2FVT (Gaussian) Sigma AIC', '2FVT (Complete) Sigma AIC', '2FVT (Simplified) Sigma AIC', '2FVT (Over-simplified) Sigma AIC', 'Ind. Par. Sigma AIC', '2FVT (Gaussian) dAIC', '2FVT (Complete) dAIC', '2FVT (Simplified) dAIC', '2FVT (Over-simplified) dAIC', 'Ind. Par. dAIC')		

		Sec.FVT.AICs<-rbind(Sec.FVT.AICs, sigmaAICs)

		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_parameters/', line,'_leaf_angles_2FVT_fitted_parameters.csv',sep=''), angle_coefs, row.names=FALSE)
		write.csv(file=paste('~/Desktop/2FVT_outputs/2FVT_model_fits/', line,'_leaf_angle_model_AICs.csv',sep=''), AIC_vals)

		Sec.FVT.pars<-rbind(Sec.FVT.pars, c(line, 'seg.hgt', sghgt_coefs[1,5:ncol(sghgt_coefs)]))
		Sec.FVT.pars<-rbind(Sec.FVT.pars, c(line, 'lf.dist', lfdist_coefs[1,5:ncol(lfdist_coefs)]))
		Sec.FVT.pars<-rbind(Sec.FVT.pars, c(line, 'lf.angle', angle_coefs[1,5:ncol(angle_coefs)]))
	}
}

colnames(Sec.FVT.pars)[1:2]<-c('genotype', 'attribute')
Sec.FVT.pars<-as.table(Sec.FVT.pars)
write.csv(file='~/Desktop/2FVT_outputs/2FVT_parameters/All_2FVT_fitted_parameters.csv', Sec.FVT.pars, row.names=FALSE)

Sec.FVT.AICs<-as.table(Sec.FVT.AICs)
write.csv(file='~/Desktop/2FVT_outputs/2FVT_model_fits/All_2FVT_sigmaAICs.csv', Sec.FVT.AICs, row.names=FALSE)

colnames(Sec.FVT.pred.flw)<-c('Obs_Flw', '70%_pred', '75%_pred', '80%_pred', '85%_pred', '90%_pred', '95%_pred', '100%_pred')
rownames(Sec.FVT.pred.flw)<-c('A10', 'RIL39', 'RIL110', 'RIL159', 'B100')
write.csv(file='~/Desktop/2FVT_outputs/2FVT_model_fits/2FVT_predicted_flowering.csv', Sec.FVT.pred.flw, row.names=TRUE)


