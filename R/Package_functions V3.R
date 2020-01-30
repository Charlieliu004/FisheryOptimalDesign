#Function 1
CV_est = function(Landings_pct,Landings_mean,Landings_var){
  if( Landings_pct < 0 | Landings_pct > 1 ) stop('Landings_pct not between 0 and 1')
  if( Landings_mean < 0 ) stop('Landings_mean not greater than 0')
  if( Landings_var < 0  ) stop('Landings_var not greater than 0')
  var = Landings_pct*(1-Landings_pct)*Landings_mean^2+Landings_pct*Landings_var
  mean = Landings_pct*Landings_mean
  CV = sqrt(var/mean^2)
  results = list("Mean" = mean, "Var" = var, "CV" = CV)
  #colnames(results) = c("Mean","Var","CV")
  return(results)
}

#Function 2
#Comparision function
#compares three different estimators for certain reporting rate and desired PSE

InterceptSampleSize = function(CVy,Mean_dockside,target_p1,target_PSE,
                               p1_obs=NULL,n_obs=NULL,
                               Mean_report=NULL, CVy_report=NULL,
                               EstMean_report=NULL,EstCVy_report=NULL,
                               deff=NULL,R=NULL,type=NULL){

  MeanLandings_dockside = Mean_dockside
  MeanLandings_s_report = Mean_report
  CVy_s_obs = CVy_report
  MeanLandings_report = EstMean_report
  CV_1y_obs = EstCVy_report


  if( target_p1 < 0 | target_p1 > 1 ) stop('target_p1 not between 0 and 1')
  if( CVy < 0 ) stop('CVy not greater than')
  if( target_PSE < 0 | target_PSE > 1 ) stop('target_PSE not between 0 and 1')


  if(!is.null(p1_obs)) {
    if( p1_obs < 0 | p1_obs > 1) stop('p1_obs not between 0 and 1')
  }
  if(!is.null(MeanLandings_dockside)) {
    if( MeanLandings_dockside < 0 ) stop('MeanLandings_dockside not greater than 0')
  }
  if(!is.null( MeanLandings_s_report)) {
    if(  MeanLandings_s_report < 0 ) stop(' MeanLandings_s_report not greater than 0')
  }
  if(!is.null( MeanLandings_report)) {
    if(  MeanLandings_report < 0 ) stop(' MeanLandings_report not greater than 0')
  }
  if(!is.null(n_obs)) {
    if( n_obs < 0 ) stop('n_obs not greater than 0')
  }
  if(!is.null(CV_1y_obs)) {
    if( CV_1y_obs < 0 ) stop('CV_1y not greater than 0')
  }
  if(!is.null(CVy_s_obs)) {
    if( CVy_s_obs < 0 ) stop('CVy_s not greater than 0')
  }
  if(!is.null(deff)) {
    if( deff < 0 ) stop('deff not greater than 0')
  }
  if(!is.null(R)) {
    if( R < 0 | R > 1) stop('R not between 0 and 1')
  }


  if(is.null(R)) {
    R = 1
    warning("R is set to 1 by default")}
  #if(is.null(k_obs)){k_obs = 1}
  #if(is.null(p1_obs)) {p1_obs = target_p1}
  if(is.null(deff)){
    deff = 2.5
    warning("Design effect is set to 2.5 by default")}

  #calculate k k_s CV_1y CVy_s
  #this condition specifies when some info about self reports exist, the corresponding rest info should also exist
  if (!is.null(MeanLandings_report)|!is.null(MeanLandings_s_report)|!is.null(CV_1y_obs)|!is.null(CVy_s_obs)){
    if (is.null(p1_obs)){ #p1_obs is null
      stop('p1_obs is missing')
    } else { #else for p1_obs exist
      #calculate k and CV1y
      #guarantee both MeanLandings_report and CV_1y_obs exist/missing
      if (!is.null(MeanLandings_report)){
        if (is.null(CV_1y_obs)){
          stop('CV_1y_obs is missing')
        } else if (target_p1 <= p1_obs){
          k =  MeanLandings_report/MeanLandings_dockside
          CV_1y = CV_1y_obs
        }#end for target p_1 < p1_obs
        else { #p_1 > p1_obs
          diff = target_p1-p1_obs
          MeanLandings_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_report)/(1-p1_obs)
          target_MeanLandings_report = (p1_obs*MeanLandings_report+diff*MeanLandings_report_c)/target_p1

          k = target_MeanLandings_report/MeanLandings_dockside

          Vy = (CVy*MeanLandings_dockside)^2
          V1y = (CV_1y_obs*MeanLandings_report)^2
          V1y_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y+MeanLandings_report^2)
                   -(1-p1_obs)*MeanLandings_report_c^2)/(1-p1_obs)
          if ( V1y_c<0) stop("Impractical inputs of self reports")
          target_V1y = (p1_obs*V1y+diff*V1y_c)/target_p1
          +p1_obs*diff*(MeanLandings_report-MeanLandings_report_c)^2/target_p1^2
          if (target_V1y<0) stop("Impractical inputs of self reports")
          CV_1y = sqrt(target_V1y)/target_MeanLandings_report
        }#end for target p_1 > p1_obs

      } else if (!is.null(CV_1y_obs)){
        stop('MeanLandings_report is missing')
      }

      #calculate k_s and CVy_s
      #guarantee both MeanLandings_s_report and CVy_s_obs exist/missing
      if (!is.null(MeanLandings_s_report)){
        if (is.null(CVy_s_obs)){
          stop('CVy_s_obs is missing')
        } else if (target_p1 <= p1_obs){
          k_s =  MeanLandings_s_report/MeanLandings_dockside
          CVy_s = CVy_s_obs
        }#end for target p_1 < p1_obs
        else { #p_1 > p1_obs
          diff = target_p1-p1_obs
          MeanLandings_s_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_s_report)/(1-p1_obs)
          target_MeanLandings_s_report = (p1_obs*MeanLandings_s_report+diff*MeanLandings_s_report_c)/target_p1

          k_s = target_MeanLandings_s_report/MeanLandings_dockside

          Vy = (CVy*MeanLandings_dockside)^2
          V1y_s = (CVy_s_obs*MeanLandings_s_report)^2
          V1y_s_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y_s+MeanLandings_s_report^2)
                     -(1-p1_obs)*MeanLandings_s_report_c^2)/(1-p1_obs)
          if ( V1y_s_c<0) stop("Impractical inputs of self reports")
          target_V1y_s = (p1_obs*V1y_s+diff*V1y_s_c)/target_p1
          +p1_obs*diff*(MeanLandings_s_report-MeanLandings_s_report_c)^2/target_p1^2
          if (target_V1y_s<0) stop("Impractical inputs of self reports")
          CVy_s = sqrt(target_V1y_s)/target_MeanLandings_s_report
        }#end for target p_1 > p1_obs

      } else if (!is.null(CVy_s_obs)){
        stop('MeanLandings_s_report is missing')
      }

      #for the case when MeanLandings_report and CV_1y_obs is missing
      if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
        k = k_s
        warning("The self reports are assumed to be accurate in terms of mean")
        if(is.null(type)){
          CV_1y = CVy_s
          warning("The self reports are assumed to be accurate in terms of variance")
        } else if (type == "accurate"){
          CV_1y = CVy_s
        }
        else if(type == "CME"){
          alpha = (1/R)^2-1
          CV_1y = round(CVy_s/sqrt(1+alpha),2)
        } else if(type == "Berkson"){
          beta = (1/R)^2-1
          CV_1y = round(CVy_s*sqrt(1+beta),2)
        }

      }
      #for the case when MeanLandings_s_report and CVy_s_obs is missing
      if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){
        k_s = k
        warning("The self reports are assumed to be accurate in terms of mean")
        if(is.null(type)){
          CVy_s = CV_1y
          warning("The self reports are assumed to be accurate in terms of variance")
        } else if (type == "accurate"){
          CVy_s = CV_1y
        }
        else if(type == "CME"){
          alpha = (1/R)^2-1
          CVy_s = round(CV_1y*sqrt(1+alpha),2)
        } else if(type == "Berkson"){
          beta = (1/R)^2-1
          CVy_s = round(CV_1y/sqrt(1+beta),2)
        }

      }
    } #end for p1_obs exist
  } else if (!is.null(p1_obs)){  #only p1_obs exist
    stop('Information about self reports is missing')
  } else {#else for no information about self reports
    k = 1
    k_s = 1
    CV_1y = CVy
    warning("The self reports are assumed to be representative of the population in terms of mean and variance")
    warning("The self reports are assumed to be accurate in terms of mean")
    #warning("The self reports are assumed to have constant mean and variance for different reporting rates p1")
    if(is.null(type)){
      CVy_s = CV_1y
      warning("The self reports are assumed to be accurate in terms of variance")
    } else if (type == "accurate"){
      CVy_s = CV_1y
    }
    else if(type == "CME"){
      alpha = (1/R)^2-1
      CVy_s = round(CV_1y*sqrt(1+alpha),2)
    } else if(type == "Berkson"){
      beta = (1/R)^2-1
      CVy_s = round(CV_1y/sqrt(1+beta),2)
    }
  }

  k = round(k,2)
  k_s = round(k_s,2)
  CV_1y = round(CV_1y,2)
  CVy_s = round(CVy_s,2)

  grid = 800
  pse = matrix(0,grid,3)
  pse_comparison = c(0,0,0)
  n2_comparison = c(0,0,0)
  p1 = target_p1

  for (ii in c(1:grid)){
    n2 = c(1:grid)[ii]

    #t_yp_PSE
    if( (CVy^2 +(1+1/p1)-2*k) < 0 ) stop('PSE less than zero')
    pse[ii,1] = sqrt((CVy^2 +(1+1/p1)-2*k)/(n2/deff))
    #t_yc_PSE
    if( (CVy^2 +(1+1/p1)-2*k+CVy_s^2/p1-2*k*R*CV_1y*CVy_s) < 0 ) stop('PSE less than zero')
    pse[ii,2] = sqrt((CVy^2 +(1+1/p1)-2*k+CVy_s^2/p1-2*k*R*CV_1y*CVy_s)/(n2/deff))
    #t_y2_PSE
    if( (CVy^2 +(1+1/p1)-2*k+p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y)) < 0 ) stop('PSE less than zero')
    pse[ii,3] = sqrt((CVy^2 +(1+1/p1)-2*k+p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y))/(n2/deff))
  }#end for ii

  n2_comparison[1] = round((deff/(target_PSE^2))*(CVy^2 +(1+1/p1)-2*k),0)
  n2_comparison[2] = round((deff/(target_PSE^2))*(CVy^2 +(1+1/p1)-2*k+CVy_s^2/p1-2*k*R*CV_1y*CVy_s),0)
  n2_comparison[3] = round((deff/(target_PSE^2))*(CVy^2 +(1+1/p1)-2*k+p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y)),0)

  plot(c(1:grid ),pse[,1],xaxt = "n",ylim = c(0,0.8),ylab = "PSE",xlab = "Intercept Sample Size",type = "l",lty=1,col = 153,las =1)
  axis(1, at=c(0,200,400,600,800), labels=c(0,200,400,600,800),cex.axis = 1,pch = 0)
  mtext(n2_comparison[1],side =1,line =-1.1, at= n2_comparison[1]-25,cex=0.8,col = 153)
  mtext(n2_comparison[2],side =1,line =-1.8, at= n2_comparison[2]-25,cex=0.8, col = 34)
  mtext(n2_comparison[3],side =1,line =-2.5, at= n2_comparison[3]-25,cex=0.8,col = 139)
  title(main = "Effect of Intercept Sample Size on PSE")

  #axis(1, at=n2_comparison, labels=n2_comparison,cex.axis = 0.7)
  axis(2, at=target_PSE, labels=target_PSE,cex.axis = 1,las =1,cex=1)
  points(c(1:20)*40,pse[,1][c(1:20)*40],pch =0,col = 153,cex = 0.7)
  lines(c(1:grid ),pse[,2],type = "l",lty = 1,col=34)
  points(20+c(1:19)*40,pse[,2][20+c(1:19)*40],pch = 1,col = 34)
  lines(c(1:grid ),pse[,3],type = "l",lty = 1,col=139)
  points(60+c(1:19)*40,pse[,3][60+c(1:19)*40],pch = 2,cex = 0.7,col = 139)
  #lines(c(1:grid ),pse[,4],type = "l",lty = 1,col=28)

  legend("topright", inset = 0.02,
         legend=c(expression(hat(t)[yp]),expression(hat(t)[yc]),expression(hat(t)[y2]))
         ,cex=1,lty=c(1,1,1),pch = c(0,1,2),col = c(153,34,139))

  legend("top", inset = 0.02,
         legend=c(paste("Target PSE = ",target_PSE),as.expression(bquote("Target"~p["1"]~"="~.(target_p1)))
                  ,as.expression(bquote(CV["y"]~"="~.(CVy)))
                  ,as.expression(bquote(CV["1y"]~"="~.(CV_1y))),as.expression(bquote(CV["1y"^'*']~"="~.(CVy_s)))
                  ,paste("Design Effect = ",deff),paste("R = ",R)),cex=1)

  lines(c(-1000,n2_comparison),rep(target_PSE,4),lty=2,col= 153)
  lines(c(rep(n2_comparison[1],2)),c(-1,target_PSE),lty=2,col=153)
  lines(c(rep(n2_comparison[2],2)),c(-1,target_PSE),lty=2,col=34)
  lines(c(rep(n2_comparison[3],2)),c(-1,target_PSE),lty=2,col=139)

  if(!is.null(n_obs)){
    pse_comparison[1] = round(sqrt((CVy^2 +(1+1/p1)-2*k)/(n_obs/deff)),3)
    pse_comparison[2] = round(sqrt((CVy^2 +(1+1/p1)-2*k+CVy_s^2/p1-2*k*R*CV_1y*CVy_s)/(n_obs/deff)),3)
    pse_comparison[3] = round(sqrt((CVy^2 +(1+1/p1)-2*k+p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y))/(n_obs/deff)),3)
    mtext(c(pse_comparison[1]),side =2,line =-1.5, at= c(pse_comparison[1]+0.02),cex=0.7,las = 1,col = 153)
    mtext(c(pse_comparison[2]),side =2,line =-3, at= c(pse_comparison[2]+0.02),cex=0.7,las = 1,col = 34)
    mtext(c(pse_comparison[3]),side =2,line =-4.5, at= c(pse_comparison[3]+0.02),cex=0.7,las = 1,col = 139)
    mtext(n_obs,side =1,line =-1.1, at= n_obs-25,cex=0.8,col = 153)
    lines(rep(n_obs,4),c(-1,pse_comparison),lty=2,col=28)
    lines(c(-1000,n_obs),rep(pse_comparison[1],2),lty=2,col=153)
    lines(c(-1000,n_obs),rep(pse_comparison[2],2),lty=2,col=34)
    lines(c(-1000,n_obs),rep(pse_comparison[3],2),lty=2,col=139)
  }

}#end for PSE_comparison


#Function 3

#Comparision function
#compares three different estimators for certain n_obs and desired PSE
ReportingRate = function(CVy,Mean_dockside,target_n2,target_PSE,p1_obs=NULL,
                         Mean_report=NULL, CVy_report=NULL,
                         EstMean_report=NULL,EstCVy_report=NULL,
                         deff=NULL,R=NULL,type=NULL){

  MeanLandings_dockside = Mean_dockside
  MeanLandings_s_report = Mean_report
  CVy_s_obs = CVy_report
  MeanLandings_report = EstMean_report
  CV_1y_obs = EstCVy_report


  if( target_n2 < 0 ) stop('target_n2 not greater than 0')
  if( CVy < 0 ) stop('CVy not greater than')
  if( target_PSE < 0 | target_PSE > 1 ) stop('target_PSE not between 0 and 1')

  if(!is.null(p1_obs)) {
    if( p1_obs < 0 | p1_obs > 1) stop('p1_obs not between 0 and 1')
  }
  if(!is.null(MeanLandings_dockside)) {
    if( MeanLandings_dockside < 0 ) stop('MeanLandings_dockside not greater than 0')
  }
  if(!is.null( MeanLandings_s_report)) {
    if(  MeanLandings_s_report < 0 ) stop(' MeanLandings_s_report not greater than 0')
  }
  if(!is.null( MeanLandings_report)) {
    if(  MeanLandings_report < 0 ) stop(' MeanLandings_report not greater than 0')
  }

  if(!is.null(CV_1y_obs)) {
    if( CV_1y_obs < 0 ) stop('CV_1y not greater than 0')
  }
  if(!is.null(CVy_s_obs)) {
    if( CVy_s_obs < 0 ) stop('CVy_s not greater than 0')
  }
  if(!is.null(deff)) {
    if( deff < 0 ) stop('deff not greater than 0')
  }
  if(!is.null(R)) {
    if( R < 0 | R > 1) stop('R not between 0 and 1')
  }


  if(is.null(R)) {
    R = 1
    warning("R is set to 1 by default")}
  #if(is.null(k_obs)){k_obs = 1}
  if(is.null(deff)){
    deff = 2.5
    warning("Design effect is set to 2.5 by default")}

  if (is.null(MeanLandings_report)&is.null(CV_1y_obs)&is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

    warning("The self reports are assumed to be representative of the population in terms of mean and variance")
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){
      #CVy_s = CV_1y
      warning("The self reports are assumed to be accurate in terms of variance")
    }

  } else if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
    #k = k_s
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){

      warning("The self reports are assumed to be accurate in terms of variance")
    }} else if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

      warning("The self reports are assumed to be accurate in terms of mean")
      if(is.null(type)){

        warning("The self reports are assumed to be accurate in terms of variance")
      }}


  grid = 100
  pse = matrix(0,grid,3)
  pse_comparison = c(0,0,0)
  p1_comparison = c(0,0,0)
  n2 = target_n2

  for (ii in c(1:grid)){

    target_p1 = c(1:grid)[ii]/100

    #calculate k k_s CV_1y CVy_s
    #this condition specifies when some info about self reports exist, the corresponding rest info should also exist
    if (!is.null(MeanLandings_report)|!is.null(MeanLandings_s_report)|!is.null(CV_1y_obs)|!is.null(CVy_s_obs)){
      if (is.null(p1_obs)){ #p1_obs is null
        stop('p1_obs is missing')
      } else { #else for p1_obs exist
        #calculate k and CV1y
        #guarantee both MeanLandings_report and CV_1y_obs exist/missing
        if (!is.null(MeanLandings_report)){
          if (is.null(CV_1y_obs)){
            stop('CV_1y_obs is missing')
          } else if (target_p1 <= p1_obs){
            k =  MeanLandings_report/MeanLandings_dockside
            CV_1y = CV_1y_obs
          }#end for target p_1 < p1_obs
          else { #p_1 > p1_obs
            diff = target_p1-p1_obs
            MeanLandings_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_report)/(1-p1_obs)
            target_MeanLandings_report = (p1_obs*MeanLandings_report+diff*MeanLandings_report_c)/target_p1

            k = target_MeanLandings_report/MeanLandings_dockside

            Vy = (CVy*MeanLandings_dockside)^2
            V1y = (CV_1y_obs*MeanLandings_report)^2
            V1y_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y+MeanLandings_report^2)
                     -(1-p1_obs)*MeanLandings_report_c^2)/(1-p1_obs)
            if ( V1y_c<0) stop("Impractical inputs of self reports")
            target_V1y = (p1_obs*V1y+diff*V1y_c)/target_p1
            +p1_obs*diff*(MeanLandings_report-MeanLandings_report_c)^2/target_p1^2
            if (target_V1y<0) stop("Impractical inputs of self reports")
            CV_1y = sqrt(target_V1y)/target_MeanLandings_report
          }#end for target p_1 > p1_obs

        } else if (!is.null(CV_1y_obs)){
          stop('MeanLandings_report is missing')
        }

        #calculate k_s and CVy_s
        #guarantee both MeanLandings_s_report and CVy_s_obs exist/missing
        if (!is.null(MeanLandings_s_report)){
          if (is.null(CVy_s_obs)){
            stop('CVy_s_obs is missing')
          } else if (target_p1 <= p1_obs){
            k_s =  MeanLandings_s_report/MeanLandings_dockside
            CVy_s = CVy_s_obs
          }#end for target p_1 < p1_obs
          else { #p_1 > p1_obs
            diff = target_p1-p1_obs
            MeanLandings_s_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_s_report)/(1-p1_obs)
            target_MeanLandings_s_report = (p1_obs*MeanLandings_s_report+diff*MeanLandings_s_report_c)/target_p1

            k_s = target_MeanLandings_s_report/MeanLandings_dockside

            Vy = (CVy*MeanLandings_dockside)^2
            V1y_s = (CVy_s_obs*MeanLandings_s_report)^2
            V1y_s_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y_s+MeanLandings_s_report^2)
                       -(1-p1_obs)*MeanLandings_s_report_c^2)/(1-p1_obs)
            if (V1y_s_c<0) stop("Impractical inputs of self reports")
            target_V1y_s = (p1_obs*V1y_s+diff*V1y_s_c)/target_p1
            +p1_obs*diff*(MeanLandings_s_report-MeanLandings_s_report_c)^2/target_p1^2
            if (target_V1y_s<0) stop("Impractical inputs of self reports")
            CVy_s = sqrt(target_V1y_s)/target_MeanLandings_s_report
          }#end for target p_1 > p1_obs

        } else if (!is.null(CVy_s_obs)){
          stop('MeanLandings_s_report is missing')
        }

        #for the case when MeanLandings_report and CV_1y_obs is missing
        if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
          k = k_s
          #warning("The self reports are assumed to be accurate in terms of mean")
          if(is.null(type)){
            CV_1y = CVy_s
            # warning("The self reports are assumed to be accurate in terms of variance")
          } else if (type == "accurate"){
            CV_1y = CVy_s
          }
          else if(type == "CME"){
            alpha = (1/R)^2-1
            CV_1y = round(CVy_s/sqrt(1+alpha),2)
          } else if(type == "Berkson"){
            beta = (1/R)^2-1
            CV_1y = round(CVy_s*sqrt(1+beta),2)
          }

        }
        #for the case when MeanLandings_s_report and CVy_s_obs is missing
        if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){
          k_s = k
          # warning("The self reports are assumed to be accurate in terms of mean")
          if(is.null(type)){
            CVy_s = CV_1y
            #   warning("The self reports are assumed to be accurate in terms of variance")
          } else if (type == "accurate"){
            CVy_s = CV_1y
          }
          else if(type == "CME"){
            alpha = (1/R)^2-1
            CVy_s = round(CV_1y*sqrt(1+alpha),2)
          } else if(type == "Berkson"){
            beta = (1/R)^2-1
            CVy_s = round(CV_1y/sqrt(1+beta),2)
          }

        }
      } #end for p1_obs exist
    } else if (!is.null(p1_obs)){  #only p1_obs exist
      stop('Information about self reports is missing')
    } else {#else for no information about self reports
      k = 1
      k_s = 1
      CV_1y = CVy
      #warning("The self reports are assumed to be representative of the population in terms of mean and variance")
      # warning("The self reports are assumed to be accurate in terms of mean")
      #warning("The self reports are assumed to have constant mean and variance for different reporting rates p1")
      if(is.null(type)){
        CVy_s = CV_1y
        #warning("The self reports are assumed to be accurate in terms of variance")
      } else if (type == "accurate"){
        CVy_s = CV_1y
      }
      else if(type == "CME"){
        alpha = (1/R)^2-1
        CVy_s = round(CV_1y*sqrt(1+alpha),2)
      } else if(type == "Berkson"){
        beta = (1/R)^2-1
        CVy_s = round(CV_1y/sqrt(1+beta),2)
      }
    }

    k = round(k,2)
    k_s = round(k_s,2)
    CV_1y = round(CV_1y,2)
    CVy_s = round(CVy_s,2)

    #t_yp_PSE
    if( (CVy^2 +(1+1/target_p1)-2*k) < 0 ) stop('PSE less than zero')
    pse[ii,1] = sqrt((CVy^2 +(1+1/target_p1)-2*k)/(n2/deff))
    #t_yc_PSE
    if( (CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s) < 0 ) stop('PSE less than zero')
    pse[ii,2] = sqrt((CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s)/(n2/deff))
    #t_y2_PSE
    if( (CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y)) < 0 ) stop('PSE less than zero')
    pse[ii,3] = sqrt((CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y))/(n2/deff))
  }#end for ii

  approx = abs(pse - target_PSE)
  p1_comparison[1]=which(approx[,1]==min(approx[,1]))/100
  p1_comparison[2]=which(approx[,2]==min(approx[,2]))/100
  p1_comparison[3]=which(approx[,3]==min(approx[,3]))/100


  plot(c(1:grid)/100,pse[,1],xaxt = "n",ylim = c(0,0.8),ylab = "PSE",xlab = "Reporting Rate",type = "l",lty=1,col = 153,las =1)
  axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1),cex.axis = 1,pch = 0)
  mtext(p1_comparison[1],side =1,line =-1.1, at= p1_comparison[1]-25/800,cex=0.8,col = 153)
  mtext(p1_comparison[2],side =1,line =-1.8, at= p1_comparison[2]-25/800,cex=0.8, col = 34)
  mtext(p1_comparison[3],side =1,line =-2.5, at= p1_comparison[3]-25/800,cex=0.8,col = 139)
  title(main = "Effect of Reporting Rate on PSE")

  #axis(1, at=n2_comparison, labels=n2_comparison,cex.axis = 0.7)
  axis(2, at=target_PSE, labels=target_PSE,cex.axis = 1,las =1,cex=1)
  points(c(1:10)*10/100,pse[,1][c(1:10)*10],pch =0,col = 153,cex = 0.7)
  lines(c(1:grid)/100,pse[,2],type = "l",lty = 1,col=34)
  points((3+c(0:9)*10)/100,pse[,2][3+c(0:9)*10],pch = 1,col = 34)
  lines(c(1:grid)/100,pse[,3],type = "l",lty = 1,col=139)
  points((6+c(0:9)*10)/100,pse[,3][6+c(0:9)*10],pch = 2,cex = 0.7,col = 139)
  #lines(c(1:grid ),pse[,4],type = "l",lty = 1,col=28)

  legend("topright", inset = 0.02,
         legend=c(expression(hat(t)[yp]),expression(hat(t)[yc]),expression(hat(t)[y2]))
         ,cex=1,lty=c(1,1,1),pch = c(0,1,2),col = c(153,34,139))

  legend("top", inset = 0.02,
         legend=c(paste("Target PSE = ",target_PSE),as.expression(bquote("Target"~n["2"]~"="~.(target_n2)))
                  ,as.expression(bquote(CV["y"]~"="~.(CVy))),paste("Design Effect = ",deff),paste("R = ",R)),cex=1)

  lines(c(-1000,p1_comparison),rep(target_PSE,4),lty=2,col= 153)
  lines(c(rep(p1_comparison[1],2)),c(-1,target_PSE),lty=2,col=153)
  lines(c(rep(p1_comparison[2],2)),c(-1,target_PSE),lty=2,col=34)
  lines(c(rep(p1_comparison[3],2)),c(-1,target_PSE),lty=2,col=139)

}#end for PSE_comparison

#Function 4
#Comparision function
#compares three different estimators for certain n_obs and desired PSE

Tradeoff = function(CVy,Mean_dockside,target_PSE,target_n2=NULL,p1_obs=NULL,
                    Mean_report=NULL, CVy_report=NULL,
                    EstMean_report=NULL,EstCVy_report=NULL,
                    deff=NULL,R=NULL,type=NULL){

  MeanLandings_dockside = Mean_dockside
  MeanLandings_s_report = Mean_report
  CVy_s_obs = CVy_report
  MeanLandings_report = EstMean_report
  CV_1y_obs = EstCVy_report


  if( CVy < 0 ) stop('CVy not greater than')
  if( target_PSE < 0 | target_PSE > 1 ) stop('target_PSE not between 0 and 1')

  if(!is.null(target_n2)) {
    if( target_n2 < 0 ) stop('target_n2 not greater than 0')
  }

  if(!is.null(p1_obs)) {
    if( p1_obs < 0 | p1_obs > 1) stop('p1_obs not between 0 and 1')
  }
  if(!is.null(MeanLandings_dockside)) {
    if( MeanLandings_dockside < 0 ) stop('MeanLandings_dockside not greater than 0')
  }
  if(!is.null( MeanLandings_s_report)) {
    if(  MeanLandings_s_report < 0 ) stop(' MeanLandings_s_report not greater than 0')
  }
  if(!is.null( MeanLandings_report)) {
    if(  MeanLandings_report < 0 ) stop(' MeanLandings_report not greater than 0')
  }

  if(!is.null(CV_1y_obs)) {
    if( CV_1y_obs < 0 ) stop('CV_1y not greater than 0')
  }
  if(!is.null(CVy_s_obs)) {
    if( CVy_s_obs < 0 ) stop('CVy_s not greater than 0')
  }
  if(!is.null(deff)) {
    if( deff < 0 ) stop('deff not greater than 0')
  }
  if(!is.null(R)) {
    if( R < 0 | R > 1) stop('R not between 0 and 1')
  }

  if(is.null(R)) {
    R = 1
    warning("R is set to 1 by default")}
  #if(is.null(k_obs)){k_obs = 1}
  if(is.null(deff)){
    deff = 2.5
    warning("Design effect is set to 2.5 by default")}

  if (is.null(MeanLandings_report)&is.null(CV_1y_obs)&is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

    warning("The self reports are assumed to be representative of the population in terms of mean and variance")
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){
      #CVy_s = CV_1y
      warning("The self reports are assumed to be accurate in terms of variance")
    }

  } else if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
    #k = k_s
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){

      warning("The self reports are assumed to be accurate in terms of variance")
    }} else if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

      warning("The self reports are assumed to be accurate in terms of mean")
      if(is.null(type)){

        warning("The self reports are assumed to be accurate in terms of variance")
      }}

  rr_comparison  = c()
  grid = 100
  SampleSize = matrix(0,grid,3)
  n2 = target_n2

  for (ii in c(1:grid)){
    target_p1 = c(1:grid)[ii]/100

    #calculate k k_s CV_1y CVy_s
    #this condition specifies when some info about self reports exist, the corresponding rest info should also exist
    if (!is.null(MeanLandings_report)|!is.null(MeanLandings_s_report)|!is.null(CV_1y_obs)|!is.null(CVy_s_obs)){
      if (is.null(p1_obs)){ #p1_obs is null
        stop('p1_obs is missing')
      } else { #else for p1_obs exist
        #calculate k and CV1y
        #guarantee both MeanLandings_report and CV_1y_obs exist/missing
        if (!is.null(MeanLandings_report)){
          if (is.null(CV_1y_obs)){
            stop('CV_1y_obs is missing')
          } else if (target_p1 <= p1_obs){
            k =  MeanLandings_report/MeanLandings_dockside
            CV_1y = CV_1y_obs
          }#end for target p_1 < p1_obs
          else { #p_1 > p1_obs
            diff = target_p1-p1_obs
            MeanLandings_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_report)/(1-p1_obs)
            target_MeanLandings_report = (p1_obs*MeanLandings_report+diff*MeanLandings_report_c)/target_p1

            k = target_MeanLandings_report/MeanLandings_dockside

            Vy = (CVy*MeanLandings_dockside)^2
            V1y = (CV_1y_obs*MeanLandings_report)^2
            V1y_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y+MeanLandings_report^2)
                     -(1-p1_obs)*MeanLandings_report_c^2)/(1-p1_obs)
            if (V1y_c<0) stop("Impractical inputs of self reports")
            target_V1y = (p1_obs*V1y+diff*V1y_c)/target_p1
            +p1_obs*diff*(MeanLandings_report-MeanLandings_report_c)^2/target_p1^2
            if (target_V1y<0) stop("Impractical inputs of self reports")
            CV_1y = sqrt(target_V1y)/target_MeanLandings_report
          }#end for target p_1 > p1_obs

        } else if (!is.null(CV_1y_obs)){
          stop('MeanLandings_report is missing')
        }

        #calculate k_s and CVy_s
        #guarantee both MeanLandings_s_report and CVy_s_obs exist/missing
        if (!is.null(MeanLandings_s_report)){
          if (is.null(CVy_s_obs)){
            stop('CVy_s_obs is missing')
          } else if (target_p1 <= p1_obs){
            k_s =  MeanLandings_s_report/MeanLandings_dockside
            CVy_s = CVy_s_obs
          }#end for target p_1 < p1_obs
          else { #p_1 > p1_obs
            diff = target_p1-p1_obs
            MeanLandings_s_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_s_report)/(1-p1_obs)
            target_MeanLandings_s_report = (p1_obs*MeanLandings_s_report+diff*MeanLandings_s_report_c)/target_p1

            k_s = target_MeanLandings_s_report/MeanLandings_dockside

            Vy = (CVy*MeanLandings_dockside)^2
            V1y_s = (CVy_s_obs*MeanLandings_s_report)^2
            V1y_s_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y_s+MeanLandings_s_report^2)
                       -(1-p1_obs)*MeanLandings_s_report_c^2)/(1-p1_obs)
            if (V1y_s_c<0) stop("Impractical inputs of self reports")
            target_V1y_s = (p1_obs*V1y_s+diff*V1y_s_c)/target_p1
            +p1_obs*diff*(MeanLandings_s_report-MeanLandings_s_report_c)^2/target_p1^2
            if (target_V1y_s<0) stop("Impractical inputs of self reports")
            CVy_s = sqrt(target_V1y_s)/target_MeanLandings_s_report
          }#end for target p_1 > p1_obs

        } else if (!is.null(CVy_s_obs)){
          stop('MeanLandings_s_report is missing')
        }

        #for the case when MeanLandings_report and CV_1y_obs is missing
        if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
          k = k_s
          #warning("The self reports are assumed to be accurate in terms of mean")
          if(is.null(type)){
            CV_1y = CVy_s
            # warning("The self reports are assumed to be accurate in terms of variance")
          } else if (type == "accurate"){
            CV_1y = CVy_s
          }
          else if(type == "CME"){
            alpha = (1/R)^2-1
            CV_1y = round(CVy_s/sqrt(1+alpha),2)
          } else if(type == "Berkson"){
            beta = (1/R)^2-1
            CV_1y = round(CVy_s*sqrt(1+beta),2)
          }

        }
        #for the case when MeanLandings_s_report and CVy_s_obs is missing
        if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){
          k_s = k
          # warning("The self reports are assumed to be accurate in terms of mean")
          if(is.null(type)){
            CVy_s = CV_1y
            #   warning("The self reports are assumed to be accurate in terms of variance")
          } else if (type == "accurate"){
            CVy_s = CV_1y
          }
          else if(type == "CME"){
            alpha = (1/R)^2-1
            CVy_s = round(CV_1y*sqrt(1+alpha),2)
          } else if(type == "Berkson"){
            beta = (1/R)^2-1
            CVy_s = round(CV_1y/sqrt(1+beta),2)
          }

        }
      } #end for p1_obs exist
    } else if (!is.null(p1_obs)){  #only p1_obs exist
      stop('Information about self reports is missing')
    } else {#else for no information about self reports
      k = 1
      k_s = 1
      CV_1y = CVy
      #warning("The self reports are assumed to be representative of the population in terms of mean and variance")
      # warning("The self reports are assumed to be accurate in terms of mean")
      #warning("The self reports are assumed to have constant mean and variance for different reporting rates p1")
      if(is.null(type)){
        CVy_s = CV_1y
        #warning("The self reports are assumed to be accurate in terms of variance")
      } else if (type == "accurate"){
        CVy_s = CV_1y
      }
      else if(type == "CME"){
        alpha = (1/R)^2-1
        CVy_s = round(CV_1y*sqrt(1+alpha),2)
      } else if(type == "Berkson"){
        beta = (1/R)^2-1
        CVy_s = round(CV_1y/sqrt(1+beta),2)
      }
    }

    k = round(k,2)
    k_s = round(k_s,2)
    CV_1y = round(CV_1y,2)
    CVy_s = round(CVy_s,2)

    PSE = target_PSE
    #t_yp_PSE
    SampleSize[ii,1] = (deff/(PSE^2))*(CVy^2 +(1+1/target_p1)-2*k)
    #t_yc_PSE
    SampleSize[ii,2] = (deff/(PSE^2))*(CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s)
    #t_y2_PSE
    SampleSize[ii,3] = (deff/(PSE^2))*(CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CVy_s))

  }#end for ii

  plot(c(1:grid)/100,SampleSize[,1],ylim = c(0,800),xaxt = "n",ylab = "Intercept Sample Size",xlab = "Report Rate",type = "l",lty=1,col = 153,las =1)
  axis(1, at=c(0,0.25,0.5,0.75,1), labels=c(0,0.25,0.5,0.75,1),cex.axis = 1,pch = 0)
  points(c(1:10)*10/100,SampleSize[,1][c(1:10)*10],pch =0,col = 153,cex = 0.7)
  lines(c(1:grid)/100,SampleSize[,2],type = "l",lty = 1,col=34)
  points((3+c(0:9)*10)/100,SampleSize[,2][3+c(0:9)*10],pch = 1,col = 34)
  lines(c(1:grid)/100,SampleSize[,3],type = "l",lty = 1,col=139)
  points((6+c(0:9)*10)/100,SampleSize[,3][6+c(0:9)*10],pch = 2,cex = 0.7,col = 139)
  #lines(c(1:grid ),pse[,4],type = "l",lty = 1,col=28)
  title(main = "Intercept Sample Size - Reporting Rate Tradeoff",cex = 0.8)

  legend("topright", inset = 0.02,
         legend=c(expression(hat(t)[yp]),expression(hat(t)[yc]),expression(hat(t)[y2]))
         ,cex=1,lty=c(1,1,1),pch = c(0,1,2),col = c(153,34,139))

  if(!is.null(target_n2)){

    legend("top", inset = 0.02,
           legend=c(paste("Target PSE = ",target_PSE),as.expression(bquote("Target"~n["2"]~"="~.(target_n2)))
                    ,as.expression(bquote(CV["y"]~"="~.(CVy))),paste("Design Effect = ",deff),paste("R = ",R)),cex=1)
    approx = abs(SampleSize - target_n2)
    rr_comparison[1]=which(approx[,1]==min(approx[,1]))
    rr_comparison[2]=which(approx[,2]==min(approx[,2]))
    rr_comparison[3]=which(approx[,3]==min(approx[,3]))
    mtext(rr_comparison[1]/100,side =1,line =-1.1, at= rr_comparison[1]/100-25/800,cex=0.8,col = 153)
    mtext(rr_comparison[2]/100,side =1,line =-1.8, at= rr_comparison[2]/100-25/800,cex=0.8, col = 34)
    mtext(rr_comparison[3]/100,side =1,line =-2.5, at= rr_comparison[3]/100-25/800,cex=0.8,col = 139)
    axis(2, at=target_n2, labels=target_n2,cex.axis = 1,las =1,cex=1)

    lines(c(-1000,rr_comparison/100),rep(target_n2,4),lty=2,col= 153)
    lines(c(rep(rr_comparison[1]/100,2)),c(-100,target_n2),lty=2,col=153)
    lines(c(rep(rr_comparison[2]/100,2)),c(-100,target_n2),lty=2,col=34)
    lines(c(rep(rr_comparison[3]/100,2)),c(-100,target_n2),lty=2,col=139)

  } else {

    legend("top", inset = 0.02,
           legend=c(paste("Target PSE = ",target_PSE),as.expression(bquote(CV["y"]~"="~.(CVy)))
                    ,paste("Design Effect = ",deff),paste("R = ",R)),cex=1)

  }


}#end for TradeOff


#Function 5
OptimalDesign = function(CVy, Mean_dockside,cost_ratio,RelBudget,
                         p1_obs=NULL,n_obs = NULL,
                         Mean_report=NULL, CVy_report=NULL,
                         EstMean_report=NULL,EstCVy_report=NULL,
                         deff=NULL,R=NULL,type=NULL){

  MeanLandings_dockside = Mean_dockside
  MeanLandings_s_report = Mean_report
  CVy_s_obs = CVy_report
  MeanLandings_report = EstMean_report
  CV_1y_obs = EstCVy_report

  if( CVy < 0 ) stop('CVy not greater than')
  if( RelBudget < 0 ) stop('RelBudget not greater than 0')
  if( cost_ratio < 0 ) stop('cost_ratio not greater than 0')

  if(!is.null(n_obs)) {
    if( n_obs < 0 ) stop('target_n2 not greater than 0')
  }
  if(!is.null(p1_obs)) {
    if( p1_obs < 0 | p1_obs > 1) stop('p1_obs not between 0 and 1')
  }
  if(!is.null(MeanLandings_dockside)) {
    if( MeanLandings_dockside < 0 ) stop('MeanLandings_dockside not greater than 0')
  }
  if(!is.null( MeanLandings_s_report)) {
    if(  MeanLandings_s_report < 0 ) stop(' MeanLandings_s_report not greater than 0')
  }
  if(!is.null( MeanLandings_report)) {
    if(  MeanLandings_report < 0 ) stop(' MeanLandings_report not greater than 0')
  }

  if(!is.null(CV_1y_obs)) {
    if( CV_1y_obs < 0 ) stop('CV_1y not greater than 0')
  }
  if(!is.null(CVy_s_obs)) {
    if( CVy_s_obs < 0 ) stop('CVy_s not greater than 0')
  }
  if(!is.null(deff)) {
    if( deff < 0 ) stop('deff not greater than 0')
  }
  if(!is.null(R)) {
    if( R < 0 | R > 1) stop('R not between 0 and 1')
  }

  if(is.null(R)) {
    R = 1
    warning("R is set to 1 by default")}
  #if(is.null(k_obs)){k_obs = 1}
  if(is.null(deff)){
    deff = 2.5
    warning("Design effect is set to 2.5 by default")}

  if (is.null(MeanLandings_report)&is.null(CV_1y_obs)&is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

    warning("The self reports are assumed to be representative of the population in terms of mean and variance")
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){
      #CVy_s = CV_1y
      warning("The self reports are assumed to be accurate in terms of variance")
    }

  } else if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
    #k = k_s
    warning("The self reports are assumed to be accurate in terms of mean")
    if(is.null(type)){

      warning("The self reports are assumed to be accurate in terms of variance")
    }} else if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){

      warning("The self reports are assumed to be accurate in terms of mean")
      if(is.null(type)){

        warning("The self reports are assumed to be accurate in terms of variance")
      }}

  if (is.null(n_obs)){#pilot study


    p1_max = floor(RelBudget/cost_ratio)-1
    grid = p1_max
    pse = matrix(0,p1_max,3)
    pse_comparison = c(0,0,0)
    p1_comparison = c(0,0,0)
    rr_comparison  = c()
    Opt_p1 = c()
    Opt_n2 = c()

    for (ii in c(1:grid)){
      target_p1 = c(1:grid)[ii]/100
      n2 = floor(RelBudget - target_p1*100*cost_ratio)

      #calculate k k_s CV_1y CVy_s
      #this condition specifies when some info about self reports exist, the corresponding rest info should also exist
      if (!is.null(MeanLandings_report)|!is.null(MeanLandings_s_report)|!is.null(CV_1y_obs)|!is.null(CVy_s_obs)){
        if (is.null(p1_obs)){ #p1_obs is null
          stop('p1_obs is missing')
        } else { #else for p1_obs exist
          #calculate k and CV1y
          #guarantee both MeanLandings_report and CV_1y_obs exist/missing
          if (!is.null(MeanLandings_report)){
            if (is.null(CV_1y_obs)){
              stop('CV_1y_obs is missing')
            } else if (target_p1 <= p1_obs){
              k =  MeanLandings_report/MeanLandings_dockside
              CV_1y = CV_1y_obs
            }#end for target p_1 < p1_obs
            else { #p_1 > p1_obs
              diff = target_p1-p1_obs
              MeanLandings_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_report)/(1-p1_obs)
              target_MeanLandings_report = (p1_obs*MeanLandings_report+diff*MeanLandings_report_c)/target_p1

              k = target_MeanLandings_report/MeanLandings_dockside

              Vy = (CVy*MeanLandings_dockside)^2
              V1y = (CV_1y_obs*MeanLandings_report)^2
              V1y_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y+MeanLandings_report^2)
                       -(1-p1_obs)*MeanLandings_report_c^2)/(1-p1_obs)
              if (V1y_c<0) stop("Impractical inputs of self reports")
              target_V1y = (p1_obs*V1y+diff*V1y_c)/target_p1
              +p1_obs*diff*(MeanLandings_report-MeanLandings_report_c)^2/target_p1^2
              if (target_V1y<0) stop("Impractical inputs of self reports")
              CV_1y = sqrt(target_V1y)/target_MeanLandings_report
            }#end for target p_1 > p1_obs

          } else if (!is.null(CV_1y_obs)){
            stop('MeanLandings_report is missing')
          }

          #calculate k_s and CVy_s
          #guarantee both MeanLandings_s_report and CVy_s_obs exist/missing
          if (!is.null(MeanLandings_s_report)){
            if (is.null(CVy_s_obs)){
              stop('CVy_s_obs is missing')
            } else if (target_p1 <= p1_obs){
              k_s =  MeanLandings_s_report/MeanLandings_dockside
              CVy_s = CVy_s_obs
            }#end for target p_1 < p1_obs
            else { #p_1 > p1_obs
              diff = target_p1-p1_obs
              MeanLandings_s_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_s_report)/(1-p1_obs)
              target_MeanLandings_s_report = (p1_obs*MeanLandings_s_report+diff*MeanLandings_s_report_c)/target_p1

              k_s = target_MeanLandings_s_report/MeanLandings_dockside

              Vy = (CVy*MeanLandings_dockside)^2
              V1y_s = (CVy_s_obs*MeanLandings_s_report)^2
              V1y_s_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y_s+MeanLandings_s_report^2)
                         -(1-p1_obs)*MeanLandings_s_report_c^2)/(1-p1_obs)
              if (V1y_s_c<0) stop("Impractical inputs of self reports")
              target_V1y_s = (p1_obs*V1y_s+diff*V1y_s_c)/target_p1
              +p1_obs*diff*(MeanLandings_s_report-MeanLandings_s_report_c)^2/target_p1^2
              if (target_V1y_s<0) stop("Impractical inputs of self reports")
              CVy_s = sqrt(target_V1y_s)/target_MeanLandings_s_report
            }#end for target p_1 > p1_obs

          } else if (!is.null(CVy_s_obs)){
            stop('MeanLandings_s_report is missing')
          }

          #for the case when MeanLandings_report and CV_1y_obs is missing
          if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
            k = k_s
            #warning("The self reports are assumed to be accurate in terms of mean")
            if(is.null(type)){
              CV_1y = CVy_s
              # warning("The self reports are assumed to be accurate in terms of variance")
            } else if (type == "accurate"){
              CV_1y = CVy_s
            }
            else if(type == "CME"){
              alpha = (1/R)^2-1
              CV_1y = round(CVy_s/sqrt(1+alpha),2)
            } else if(type == "Berkson"){
              beta = (1/R)^2-1
              CV_1y = round(CVy_s*sqrt(1+beta),2)
            }

          }
          #for the case when MeanLandings_s_report and CVy_s_obs is missing
          if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){
            k_s = k
            # warning("The self reports are assumed to be accurate in terms of mean")
            if(is.null(type)){
              CVy_s = CV_1y
              #   warning("The self reports are assumed to be accurate in terms of variance")
            } else if (type == "accurate"){
              CVy_s = CV_1y
            }
            else if(type == "CME"){
              alpha = (1/R)^2-1
              CVy_s = round(CV_1y*sqrt(1+alpha),2)
            } else if(type == "Berkson"){
              beta = (1/R)^2-1
              CVy_s = round(CV_1y/sqrt(1+beta),2)
            }

          }
        } #end for p1_obs exist
      } else if (!is.null(p1_obs)){  #only p1_obs exist
        stop('Information about self reports is missing')
      } else {#else for no information about self reports
        k = 1
        k_s = 1
        CV_1y = CVy
        #warning("The self reports are assumed to be representative of the population in terms of mean and variance")
        # warning("The self reports are assumed to be accurate in terms of mean")
        #warning("The self reports are assumed to have constant mean and variance for different reporting rates p1")
        if(is.null(type)){
          CVy_s = CV_1y
          #warning("The self reports are assumed to be accurate in terms of variance")
        } else if (type == "accurate"){
          CVy_s = CV_1y
        }
        else if(type == "CME"){
          alpha = (1/R)^2-1
          CVy_s = round(CV_1y*sqrt(1+alpha),2)
        } else if(type == "Berkson"){
          beta = (1/R)^2-1
          CVy_s = round(CV_1y/sqrt(1+beta),2)
        }
      }

      k = round(k,2)
      k_s = round(k_s,2)
      CV_1y = round(CV_1y,2)
      CVy_s = round(CVy_s,2)

      #t_yp_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k) < 0 ) stop('PSE less than zero')
      pse[ii,1] = sqrt((CVy^2 +(1+1/target_p1)-2*k)/(n2/deff))
      #t_yc_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s) < 0 ) stop('PSE less than zero')
      pse[ii,2] = sqrt((CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s)/(n2/deff))
      #t_y2_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y)) < 0 ) stop('PSE less than zero')
      pse[ii,3] = sqrt((CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y))/(n2/deff))
    }#end for ii

    y_max = max(pse)+0.1
    Opt_p1[1] = which(pse[,1]==min(pse[,1]))
    Opt_p1[2] = which(pse[,2]==min(pse[,2]))
    Opt_p1[3] = which(pse[,3]==min(pse[,3]))

    Opt_n2[1] = floor(RelBudget - Opt_p1[1]*cost_ratio)
    Opt_n2[2] = floor(RelBudget - Opt_p1[2]*cost_ratio)
    Opt_n2[3] = floor(RelBudget - Opt_p1[3]*cost_ratio)

    plot(c(1:grid)/100,pse[,1],ylim = c(0,min(y_max,1)),ylab = "PSE",xlab = "Reporting Rate",type = "l",lty=1,col = 153,las =1)
    mtext(Opt_p1[1]/100,side =1,line =-1.1, at= Opt_p1[1]/100-0.01,cex=0.8,col = 153)
    mtext(Opt_p1[2]/100,side =1,line =-1.8, at= Opt_p1[2]/100-0.01,cex=0.8, col = 34)
    mtext(Opt_p1[3]/100,side =1,line =-2.5, at= Opt_p1[3]/100-0.01,cex=0.8,col = 139)
    title(main = "Optimal Sampling Design")

    points(c(1:round(grid/3))*3/100,pse[,1][c(1:round(grid/3))*3],pch =0,col = 153,cex = 0.7)
    lines(c(1:grid)/100,pse[,2],type = "l",lty = 1,col=34)
    points(2/100+c(1:round(grid/3))*3/100,pse[,2][2+c(1:round(grid/3))*3],pch = 1,col = 34)
    lines(c(1:grid)/100,pse[,3],type = "l",lty = 1,col=139)
    points(4/100+c(1:round(grid/3))*3/100,pse[,3][4+c(1:round(grid/3))*3],pch = 2,cex = 0.7,col = 139)

    pse_comparison[1] = min(pse[,1])
    pse_comparison[2] = min(pse[,2])
    pse_comparison[3] = min(pse[,3])
    mtext(c(round(pse_comparison[1],2)),side =2,line =-1.5, at= c(pse_comparison[1]+0.02),cex=0.7,las = 1,col = 153)
    mtext(c(round(pse_comparison[2],2)),side =2,line =-3, at= c(pse_comparison[2]+0.02),cex=0.7,las = 1,col = 34)
    mtext(c(round(pse_comparison[3],2)),side =2,line =-4.5, at= c(pse_comparison[3]+0.02),cex=0.7,las = 1,col = 139)


    lines(c(-1000,Opt_p1[1]/100),rep(pse_comparison[1],2),lty=2,col= 153)
    lines(c(-1000,Opt_p1[2]/100),rep(pse_comparison[2],2),lty=2,col= 34)
    lines(c(-1000,Opt_p1[3]/100),rep(pse_comparison[3],2),lty=2,col= 139)
    lines(c(rep(Opt_p1[1]/100,2)),c(-10000,pse_comparison[1]),lty=2,col=153)
    lines(c(rep(Opt_p1[2]/100,2)),c(-10000,pse_comparison[2]),lty=2,col=34)
    lines(c(rep(Opt_p1[3]/100,2)),c(-10000,pse_comparison[3]),lty=2,col=139)

    legend("top", inset = 0.02,
           legend=c(as.expression(bquote(hat(t)[yp]~" optimal p1 = " ~  .(Opt_p1[1]/100)~", n2 = " ~  .(Opt_n2[1]))),
                    as.expression(bquote(hat(t)[yc]~" optimal p1 = " ~  .(Opt_p1[2]/100)~", n2 = " ~  .(Opt_n2[2]))),
                    as.expression(bquote(hat(t)[y2]~" optimal p1 = " ~  .(Opt_p1[3]/100)~", n2 = " ~  .(Opt_n2[3]))))
           ,cex=0.8,lty=c(1,1,1),pch = c(0,1,2),col = c(153,34,139))

    legend("bottomright", inset = 0.02,
           legend=c(paste("RelBudget = ",RelBudget),paste("Cost Ratio = ",cost_ratio )
                    ,as.expression(bquote(CV["y"]~"="~.(CVy))),paste("Design Effect = ",deff),paste("R = ",R)),cex=0.8)

  }#end for pilot study
  else { #Optimal resources allocation

    p1_max = p1_obs*100+floor(RelBudget/cost_ratio)-1
    grid = floor(RelBudget/cost_ratio)
    pse = matrix(0,grid,3)
    pse_comparison = c(0,0,0)
    p1_comparison = c(0,0,0)
    rr_comparison  = c()
    Opt_p1 = c()
    Opt_n2 = c()


    for (ii in c(1:grid)){
      target_p1 = p1_obs+c(1:grid)[ii]/100-0.01
      n2 = n_obs+floor(RelBudget - (c(1:grid)[ii]-1)*cost_ratio)


      #calculate k k_s CV_1y CVy_s
      #this condition specifies when some info about self reports exist, the corresponding rest info should also exist
      if (!is.null(MeanLandings_report)|!is.null(MeanLandings_s_report)|!is.null(CV_1y_obs)|!is.null(CVy_s_obs)){
        if (is.null(p1_obs)){ #p1_obs is null
          stop('p1_obs is missing')
        } else { #else for p1_obs exist
          #calculate k and CV1y
          #guarantee both MeanLandings_report and CV_1y_obs exist/missing
          if (!is.null(MeanLandings_report)){
            if (is.null(CV_1y_obs)){
              stop('CV_1y_obs is missing')
            } else if (target_p1 <= p1_obs){
              k =  MeanLandings_report/MeanLandings_dockside
              CV_1y = CV_1y_obs
            }#end for target p_1 < p1_obs
            else { #p_1 > p1_obs
              diff = target_p1-p1_obs
              MeanLandings_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_report)/(1-p1_obs)
              target_MeanLandings_report = (p1_obs*MeanLandings_report+diff*MeanLandings_report_c)/target_p1

              k = target_MeanLandings_report/MeanLandings_dockside

              Vy = (CVy*MeanLandings_dockside)^2
              V1y = (CV_1y_obs*MeanLandings_report)^2
              V1y_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y+MeanLandings_report^2)
                       -(1-p1_obs)*MeanLandings_report_c^2)/(1-p1_obs)
              if (V1y_c<0) stop("Impractical inputs of self reports")
              target_V1y = (p1_obs*V1y+diff*V1y_c)/target_p1
              +p1_obs*diff*(MeanLandings_report-MeanLandings_report_c)^2/target_p1^2
              if (target_V1y<0) stop("Impractical inputs of self reports")
              CV_1y = sqrt(target_V1y)/target_MeanLandings_report
            }#end for target p_1 > p1_obs

          } else if (!is.null(CV_1y_obs)){
            stop('MeanLandings_report is missing')
          }

          #calculate k_s and CVy_s
          #guarantee both MeanLandings_s_report and CVy_s_obs exist/missing
          if (!is.null(MeanLandings_s_report)){
            if (is.null(CVy_s_obs)){
              stop('CVy_s_obs is missing')
            } else if (target_p1 <= p1_obs){
              k_s =  MeanLandings_s_report/MeanLandings_dockside
              CVy_s = CVy_s_obs
            }#end for target p_1 < p1_obs
            else { #p_1 > p1_obs
              diff = target_p1-p1_obs
              MeanLandings_s_report_c = (MeanLandings_dockside - p1_obs*MeanLandings_s_report)/(1-p1_obs)
              target_MeanLandings_s_report = (p1_obs*MeanLandings_s_report+diff*MeanLandings_s_report_c)/target_p1

              k_s = target_MeanLandings_s_report/MeanLandings_dockside

              Vy = (CVy*MeanLandings_dockside)^2
              V1y_s = (CVy_s_obs*MeanLandings_s_report)^2
              V1y_s_c = (Vy+MeanLandings_dockside^2-p1_obs*(V1y_s+MeanLandings_s_report^2)
                         -(1-p1_obs)*MeanLandings_s_report_c^2)/(1-p1_obs)
              if (V1y_s_c<0) stop("Impractical inputs of self reports")
              target_V1y_s = (p1_obs*V1y_s+diff*V1y_s_c)/target_p1
              +p1_obs*diff*(MeanLandings_s_report-MeanLandings_s_report_c)^2/target_p1^2
              if (target_V1y_s<0) stop("Impractical inputs of self reports")
              CVy_s = sqrt(target_V1y_s)/target_MeanLandings_s_report
            }#end for target p_1 > p1_obs

          } else if (!is.null(CVy_s_obs)){
            stop('MeanLandings_s_report is missing')
          }

          #for the case when MeanLandings_report and CV_1y_obs is missing
          if (is.null(MeanLandings_report)&is.null(CV_1y_obs)){
            k = k_s
            #warning("The self reports are assumed to be accurate in terms of mean")
            if(is.null(type)){
              CV_1y = CVy_s
              # warning("The self reports are assumed to be accurate in terms of variance")
            } else if (type == "accurate"){
              CV_1y = CVy_s
            }
            else if(type == "CME"){
              alpha = (1/R)^2-1
              CV_1y = round(CVy_s/sqrt(1+alpha),2)
            } else if(type == "Berkson"){
              beta = (1/R)^2-1
              CV_1y = round(CVy_s*sqrt(1+beta),2)
            }

          }
          #for the case when MeanLandings_s_report and CVy_s_obs is missing
          if (is.null(MeanLandings_s_report)&is.null(CVy_s_obs)){
            k_s = k
            # warning("The self reports are assumed to be accurate in terms of mean")
            if(is.null(type)){
              CVy_s = CV_1y
              #   warning("The self reports are assumed to be accurate in terms of variance")
            } else if (type == "accurate"){
              CVy_s = CV_1y
            }
            else if(type == "CME"){
              alpha = (1/R)^2-1
              CVy_s = round(CV_1y*sqrt(1+alpha),2)
            } else if(type == "Berkson"){
              beta = (1/R)^2-1
              CVy_s = round(CV_1y/sqrt(1+beta),2)
            }

          }
        } #end for p1_obs exist
      } else if (!is.null(p1_obs)){  #only p1_obs exist
        stop('Information about self reports is missing')
      } else {#else for no information about self reports
        k = 1
        k_s = 1
        CV_1y = CVy
        #warning("The self reports are assumed to be representative of the population in terms of mean and variance")
        # warning("The self reports are assumed to be accurate in terms of mean")
        #warning("The self reports are assumed to have constant mean and variance for different reporting rates p1")
        if(is.null(type)){
          CVy_s = CV_1y
          #warning("The self reports are assumed to be accurate in terms of variance")
        } else if (type == "accurate"){
          CVy_s = CV_1y
        }
        else if(type == "CME"){
          alpha = (1/R)^2-1
          CVy_s = round(CV_1y*sqrt(1+alpha),2)
        } else if(type == "Berkson"){
          beta = (1/R)^2-1
          CVy_s = round(CV_1y/sqrt(1+beta),2)
        }
      }

      k = round(k,2)
      k_s = round(k_s,2)
      CV_1y = round(CV_1y,2)
      CVy_s = round(CVy_s,2)


      #t_yp_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k) < 0 ) stop('PSE less than zero')
      pse[ii,1] = sqrt((CVy^2 +(1+1/target_p1)-2*k)/(n2/deff))
      #t_yc_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s) < 0 ) stop('PSE less than zero')
      pse[ii,2] = sqrt((CVy^2 +(1+1/target_p1)-2*k+CVy_s^2/target_p1-2*k*R*CV_1y*CVy_s)/(n2/deff))
      #t_y2_PSE
      if( (CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y)) < 0 ) stop('PSE less than zero')
      pse[ii,3] = sqrt((CVy^2 +(1+1/target_p1)-2*k+target_p1*k_s*CVy_s*(k_s*CVy_s-2*k*R*CV_1y))/(n2/deff))
    }#end for ii

    y_max = max(pse)+0.1
    Opt_p1[1] = p1_obs*100+which(pse[,1]==min(pse[,1]))-1
    Opt_p1[2] = p1_obs*100+which(pse[,2]==min(pse[,2]))-1
    Opt_p1[3] = p1_obs*100+ which(pse[,3]==min(pse[,3]))-1

    Opt_n2[1] = n_obs+floor(RelBudget - (Opt_p1[1]-p1_obs*100)*cost_ratio)
    Opt_n2[2] = n_obs+floor(RelBudget - (Opt_p1[2]-p1_obs*100)*cost_ratio)
    Opt_n2[3] = n_obs+floor(RelBudget - (Opt_p1[3]-p1_obs*100)*cost_ratio)

    plot(p1_obs+c(1:grid)/100-0.01,pse[,1],xlim = c(p1_obs,p1_obs+(grid)/100),ylim = c(0,min(y_max,1)),ylab = "PSE",xlab = "Reporting Rate",type = "l",lty=1,col = 153,las =1)
    mtext(Opt_p1[1]/100,side =1,line =-1.1, at= Opt_p1[1]/100-0.005,cex=0.8,col = 153)
    mtext(Opt_p1[2]/100,side =1,line =-1.8, at= Opt_p1[2]/100-0.005,cex=0.8, col = 34)
    mtext(Opt_p1[3]/100,side =1,line =-2.5, at= Opt_p1[3]/100-0.005,cex=0.8,col = 139)
    title(main = "Optimal Sampling Design")

    points(p1_obs+c(1:round(grid/3))*3/100,pse[,1][c(1:round(grid/3))*3],pch =0,col = 153,cex = 0.7)
    lines(p1_obs+c(1:grid)/100,pse[,2],type = "l",lty = 1,col=34)
    points(p1_obs+2/100+c(1:round(grid/3))*3/100,pse[,2][2+c(1:round(grid/3))*3],pch = 1,col = 34)
    lines(p1_obs+c(1:grid)/100,pse[,3],type = "l",lty = 1,col=139)
    points(p1_obs+4/100+c(1:round(grid/3))*3/100,pse[,3][4+c(1:round(grid/3))*3],pch = 2,cex = 0.7,col = 139)

    pse_comparison[1] = min(pse[,1])
    pse_comparison[2] = min(pse[,2])
    pse_comparison[3] = min(pse[,3])
    mtext(c(round(pse_comparison[1],2)),side =2,line =-1.5, at= c(pse_comparison[1]+0.02),cex=0.7,las = 1,col = 153)
    mtext(c(round(pse_comparison[2],2)),side =2,line =-3, at= c(pse_comparison[2]+0.02),cex=0.7,las = 1,col = 34)
    mtext(c(round(pse_comparison[3],2)),side =2,line =-4.5, at= c(pse_comparison[3]+0.02),cex=0.7,las = 1,col = 139)

    lines(c(-1000,Opt_p1[1]/100),rep(pse_comparison[1],2),lty=2,col= 153)
    lines(c(-1000,Opt_p1[2]/100),rep(pse_comparison[2],2),lty=2,col= 34)
    lines(c(-1000,Opt_p1[3]/100),rep(pse_comparison[3],2),lty=2,col= 139)
    lines(c(rep(Opt_p1[1]/100,2)),c(-10000,pse_comparison[1]),lty=2,col=153)
    lines(c(rep(Opt_p1[2]/100,2)),c(-10000,pse_comparison[2]),lty=2,col=34)
    lines(c(rep(Opt_p1[3]/100,2)),c(-10000,pse_comparison[3]),lty=2,col=139)

    legend("top", inset = 0.02,
           legend=c(as.expression(bquote(hat(t)[yp]~" optimal p1 = " ~  .(Opt_p1[1]/100)~", n2 = " ~  .(Opt_n2[1]))),
                    as.expression(bquote(hat(t)[yc]~" optimal p1 = " ~  .(Opt_p1[2]/100)~", n2 = " ~  .(Opt_n2[2]))),
                    as.expression(bquote(hat(t)[y2]~" optimal p1 = " ~  .(Opt_p1[3]/100)~", n2 = " ~  .(Opt_n2[3]))))
           ,cex=0.8,lty=c(1,1,1),pch = c(0,1,2),col = c(153,34,139))

    legend("bottomright", inset = 0.02,
           legend=c(paste("RelBudget = ",RelBudget),paste("Cost Ratio = ",cost_ratio )
                    ,as.expression(bquote(CV["y"]~"="~.(CVy))),paste("Current p1 = ",p1_obs),
                    paste("Current n2 = ",n_obs),paste("Design Effect = ",deff),paste("R = ",R)),cex=0.8)


  } #end for resources allocation
} #End for OptimalDesign
