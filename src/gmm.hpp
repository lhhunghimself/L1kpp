#include <limits>
#include <ctime>
#include <fstream>
#include <cmath>
#include <random>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/digamma.hpp>

#define LNINVSQRT2PI -0.918938533
#define INVSQRT2PI 0.3989422804
#define LN2 0.69314718056
#define PEAKWIDTH 0.325
#define MAX_EM_ITERATIONS 500
#define EM_CONVERGE 1e-3
//this is used for the minimum logValue as well - to distinguish it from an undef value which is 0
#define MIN_NOISE 1e-8
 //lookup tables
static double logx[32768]; //log base2
static double lnx[32768];  //natural logs
 
void init_logx();
template <class T> int fitLogNormEM(T *smooth, int nBins,vector<T> &peakBins,vector<T> &peakValues, vector<T> &peakMeans);
using namespace std;
using namespace boost::math;


template <class T> class MedianRecord{
	public:
	 vector<T> means0;
	 vector<T> means1;
	 T median0;
	 T median1;
		MedianRecord(){		
   blank();
	 }
  void blank(){
			const T null=std::numeric_limits<T>::quiet_NaN();
			means0.clear();
			means1.clear();
			median0=null;
		 median1=null;				
		}

	 void addValues(T mean0, T mean1){
			if(!mean0 || !mean1 || isnan(mean0) || isnan(mean1)) return;
			means0.push_back(mean0);
			means1.push_back(mean1);
		}
		void calculate_medians()	{
			if( means0.size() &&  means1.size() ){
    median0=median(means0.begin(),means0.end());
	   median1=median(means1.begin(),means1.end());
			}
		}
		private:	
  template <typename Iterator>  T median(Iterator begin, Iterator end) {
	 Iterator middle = begin + (end - begin)/2;
	 std::nth_element(begin, middle, end);
  if ((end - begin) % 2) return *middle;
		Iterator lower_middle = std::max_element(begin, middle);
		return (*middle + *lower_middle) / 2.0;
 }	   
};

template <class T> class DeconRecord{
	public:
	T means0 [501];
	T means1 [501];
	T stdevs0 [501];
	T stdevs1 [501];
	T a0[501];
	T a1[501];
	T noiseFractions[501];
	T meansEst0[501];
	T meansEst1[501];
	int dataSizes[501];
	T loglikelihoods[501];
	bool flipped[501];
	bool converged[501];
	DeconRecord(){
		blank();
	}
	void blank(){
		const T null=std::numeric_limits<T>::quiet_NaN();
		fill(means0,means0+501,null);
		fill(means1,means1+501,null);
		fill(stdevs0,stdevs0+501,null);
		fill(stdevs1,stdevs1+501,null);
		fill(a0,a0+501,null);
		fill(a1,a1+501,null);
		fill(loglikelihoods,loglikelihoods+501,null);
		fill(meansEst0,meansEst0+501,null);
		fill(meansEst1,meansEst1+501,null);
		fill(dataSizes,dataSizes+501,0);
		fill(noiseFractions,noiseFractions+501,0);
		fill(flipped,flipped+501,0);			
		fill(converged,converged+501,0);					
	}
	void transfer(int i,T *meanVars){
		meansEst0[i]=meanVars[0];
		meansEst1[i]=meanVars[1];
	}
	
	void print_binary(FILE *fp){
		//prints out a short binary file with just the fields
		for (int i=11;i<501;i++){
			fwrite(means0+i,sizeof(T),1,fp);
			fwrite(means1+i,sizeof(T),1,fp);
		}
	}		
	void print_binary(FILE *fp,int colorStart,int colorEnd){
		//prints out a short binary file with just the fields
		for (int i=colorStart;i<=colorEnd;i++){
			fwrite(means0+i,sizeof(T),1,fp);
			fwrite(means1+i,sizeof(T),1,fp);
		}
	}
	void print_text(FILE *fp,int verbosity){
		if(verbosity <0){
			print_binary(fp);
			return;
  }	
		string headers[14]={"Analyte","Mean1","Mean2","StDev1","StDev2","a1","a2","LogLike","DataSize","NoiseFraction","Converged","Flip","CMean1","CMean2"};
		int nfields=3;
		if(verbosity ==1) nfields=8;
		else if (verbosity >1) nfields=14;
		for(int i=0;i<nfields-1;i++){
			fprintf(fp,"%s\t",headers[i].c_str());
		}	
		fprintf(fp,"%s\n",headers[nfields-1].c_str());
		for (int i=1;i<501;i++){
			if(verbosity==0) fprintf(fp,"%3d\t%f\t%f\n",i,means0[i],means1[i]);	
			else if(verbosity==1)	fprintf(fp,"%3d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n",i,means0[i],means1[i],stdevs0[i],stdevs1[i],a0[i],a1[i],loglikelihoods);
   else{
				if((means0[i] > means1 [i] && 	meansEst0[i] < meansEst1[i]) ||  (means0[i] < means1 [i] && 	meansEst0[i] > meansEst1[i])){
					flipped[i]=1;
				}	
				fprintf(fp,"%3d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%d\t%f\t%s\t%s\t%f\t%f\n",i,means0[i],means1[i],stdevs0[i],stdevs1[i],a0[i],a1[i],loglikelihoods[i],dataSizes[i],noiseFractions[i],"false\0true"+6*(int)converged[i],"false\0true"+6*(int)flipped[i],meansEst0[i],meansEst1[i]);
			}		
		}	
	}
	void print_text(FILE *fp,int colorStart,int colorEnd,int verbosity){
		if(verbosity <0){
			print_binary(fp,colorStart,colorEnd);
			return;
  }	
		string headers[14]={"Analyte","Mean1","Mean2","StDev1","StDev2","a1","a2","LogLike","DataSize","NoiseFraction","Converged","Flip","CMean1","CMean2"};
		int nfields=3;
		if(verbosity ==1) nfields=8;
		else if (verbosity >1) nfields=14;
		for(int i=0;i<nfields-1;i++){
			fprintf(fp,"%s\t",headers[i].c_str());
		}	
		fprintf(fp,"%s\n",headers[nfields-1].c_str());
		for (int i=colorStart;i<=colorEnd;i++){
			if(verbosity==0) fprintf(fp,"%3d\t%f\t%f\n",i,means0[i],means1[i]);	
			else if(verbosity==1)	fprintf(fp,"%3d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\n",i,means0[i],means1[i],stdevs0[i],stdevs1[i],a0[i],a1[i],loglikelihoods);
   else{
				if((means0[i] > means1 [i] && 	meansEst0[i] < meansEst1[i]) ||  (means0[i] < means1 [i] && 	meansEst0[i] > meansEst1[i])){
					flipped[i]=1;
				}
				//something wrong with format string near space <-fix	
				fprintf(fp,"%3d\t%f\t%f\t%f\t%f\t%f\t%f\t%e\t%d\t%f\t%s\t%s \t%f\t%f\n",i,means0[i],means1[i],stdevs0[i],stdevs1[i],a0[i],a1[i],loglikelihoods[i],dataSizes[i],noiseFractions[i],"false\0true"+6*(int)converged[i],"false\0true"+6*(int)flipped[i],meansEst0[i],meansEst1[i]);
			}		
		}	
	}
};
template <class T> class EM_GMM{
	public:
	 //these are for lognormals usually expressed using base 2 - we change to natural logs so that normalization works
	 T var0,var1,var_est,var_lower,var_upper,mean0,mean1,meanEst0,meanEst1,mean_lower,mean_upper,a0,a1,p_est,p_upper,p_sigma;
	 T *data;
	 int nBins;
	 int nPoints;
	 int maxValue;
	 int minValue;
	 T loglikelihood;
	 bool converged;
	 bool normed; //normalized data?

	 EM_GMM(T *_data,bool _normed,int _nBins,int _nPoints, int _minValue,int _maxValue,T _mus[2],T _vars[2],T _a, T _p_est, T _p_sigma,T _p_upper, T _sigma_lower, T _sigma_est, T _sigma_upper){
			nPoints=_nPoints;
			mean0=_mus[0]*LN2;
			mean1=_mus[1]*LN2;
			meanEst0=mean0;
			meanEst1=mean1;
			var0=_vars[0]*LN2*LN2;
			var1=_vars[1]*LN2*LN2;;
			var_est=_sigma_est*_sigma_est*LN2*LN2;   
			var_lower=_sigma_lower*_sigma_lower*LN2*LN2;
			var_upper=_sigma_upper*_sigma_upper*LN2*LN2;
			a0=_a;
			a1=1.0-a0;
			data=_data;
			nBins=_nBins;
			minValue=_minValue;
			mean_lower=lnx[minValue];
			maxValue=_maxValue;
			mean_upper=lnx[maxValue];
			p_est=_p_est;	 
   p_upper=_p_upper;
   p_sigma=_p_sigma;
   loglikelihood=0;
   normed=_normed;
			//hard estimate for p distorts the distribution
			//with large number of points (i.e. composite of many wells) the spread of possible p values needs incresased as using a binomial correction is the same as a hard limit
			//we assume that the mixture value is not exactly _p since we have estimates from .625 - .700 and there is probably some mixture error in addition to sampling variation which would be binomial
			//for lower numbers (i.e. individual wells) we use the estimate of mixture error to determine the beta distribution and add binomial correction 
			//for higher numbers the data swamps the prior if we do this - so we explicitly calculate beta to limit the possible values 
			
			//when given limiting based on distribution
			//calculate alpha and beta of beta distribution
			//alpha=((1-_p)/(_pErr*_pErr) -1/_p)*_p*_p
			//beta=alpha*(1/_p-1)
			//update is now (alpha+k)/(alpha+beta+n) instead of k/n
			//use this formula to derive beta-binomial correction to loglikelihood 
			
			//will not use k/n as estimate of partition ratio - instead calculate an alpha and beta so that the
			//so that the maximum value if the entire sample partitioned 1 way the value equal to the limit
   //alpha =-(N*p_est*p_upper-N*p_est)/(p_upper-p_est)
   //beta=((N*p_est-N)*p_upper-N*p_est+N)/(p_upper-p_est)
   
   //for betabinomial loglikelihood for y successes in n trials 
   //lgamma(n+a+b)-lgamma(y+alpha) -lgamma(n-y+beta) +(alpha+y-1)log(p)+(b+n-y-1)*log(1-p)
   
   //NB check and flip a1 a0 if a1>a0 - check - define alarger asmaller  
	 }
	 void reset(T _mean0,T _mean1, T _var0,T _var1, T _a){
			mean0=_mean0*LN2;;
			mean1=_mean1*LN2;;
			var0=_var0*LN2*LN2;
			var1=_var1*LN2*LN2;
			a0=_a;
			a1=1.0-a0;
			loglikelihood=0;
		}	
	 int seed_mclust(T input_p){
   //use m-clust to provide initial points
   //very similar to optimize function but with different estimates
   const int maxIterations=MAX_EM_ITERATIONS;
   const double converge=EM_CONVERGE;
   bool done=0;
   int nIterations=0;
   const int last=(nBins-1<maxValue)? nBins-1 : maxValue;
   T df=3.0,kp=0.01,svar=var_est;
   
   
  	//T alpha1=((1.0-p_est)/(p_sigma*p_sigma) -1.0/p_est)*p_est*p_est;
			//T beta1=alpha1*(1.0/(p_est-1));

   //fit mean stdevs to normals norma a*invsqrtpi/sigma * exp((-(x-mean)**2)/(2*var*var))
   //however the  invsqrtpi/sigma needs to be in natural log scale to properly normalize
   
   //adjust ln
   vector <int> map;
   vector <T> mylnx; //mapped copy of lnx
   vector <T> mydata;
   for (int j=minValue;j<=last; j++){
			//check for underflow
			 if(data[j] && !isnan(data[j])){
					map.push_back(j);
					mylnx.push_back(lnx[j]);
					mydata.push_back(data[j]);
				}	
			} 
			
   while(!done){
    double term0=log(a0)+LNINVSQRT2PI-log(sqrt(var0));
    double term1=log(a1)+LNINVSQRT2PI-log(sqrt(var1));
    double inv2Var0=1.0/(2.0*var0);
    double inv2Var1=1.0/(2.0*var1); 

   
  //calculate weighted percentage that points in j come from distribution 1 or 2
  //then calculate stdevs and mean for expectation maximization
   double L=0,sum0=0,sumsq0=0,sum1=0,sumsq1=0,counts0=0,counts1=0;
   const int setSize=map.size();
   for (int m=0;m<setSize; m++){
				const int j=map[m];
			//check for underflow
			 	T logLike0=term0-inv2Var0*(mylnx[m]-mean0)*(mylnx[m]-mean0);
			 	T logLike1=term1-inv2Var1*(mylnx[m]-mean1)*(mylnx[m]-mean1);

			 	//check for underflows
			 	if((!logLike1 || std::isnan(logLike1)) && logLike0 && !std::isnan(logLike0)){
		     L+=mydata[m]*logLike0;
			 		 counts0+=mydata[m];
			 		 sum0+=mylnx[m]*mydata[m];
			 		 sumsq0+=mylnx[m]*mylnx[m]*mydata[m];
			 	}
			 	else if	((!logLike0 || std::isnan(logLike0)) && logLike1 && !std::isnan(logLike1)){
		     L+=mydata[m]*logLike1;
			 		 counts1+=mydata[m];
			 		 sum1+=mylnx[m]*mydata[m];
			 		 sumsq1+=mylnx[m]*mylnx[m]*mydata[m];
			 	}
			 	else{
			 		const T logRatio=logLike0-logLike1;
			 	 if(logRatio > 15){
		     L+=mydata[m]*logLike0;
			 		 counts0+=mydata[m];
			 		 sum0+=mylnx[m]*mydata[m];
			 		 sumsq0+=mylnx[m]*mylnx[m]*mydata[m];	
						}
						else if (logRatio < -15){
		     L+=mydata[m]*logLike1;
			 		 counts1+=mydata[m];
			 		 sum1+=mylnx[m]*mydata[m];
			 		 sumsq1+=mylnx[m]*mylnx[m]*mydata[m];
						}
						else{
			 		 const T ratio=exp(logRatio);
			 		 const T r0=ratio/(ratio+1.0);
			 		 const T r1=1.0-r0;
			 	
					  T nPoints0=mydata[m]*r0;
				 	 T nPoints1=mydata[m]*r1;
				 	 

				   counts0+=nPoints0;
				   counts1+=nPoints1;
				   L+=nPoints0*logLike0+nPoints1*logLike1;
				 	 sum0+=nPoints0*mylnx[m];
					  sumsq0+=nPoints0*mylnx[m]*mylnx[m];
				 	 sum1+=nPoints1*mylnx[m];
					  sumsq1+=nPoints1*mylnx[m]*mylnx[m];
					 }
					}		
			}
		 //calculate estimate of mean
		 T mean0MLE=(sum0+kp*meanEst0)/(kp+counts0); 
   T mean1MLE=(sum1+kp*meanEst1)/(kp+counts1); 

   mean0MLE= (mean0MLE > mean_upper)? mean_upper : mean0MLE;  
   mean0MLE= (mean0MLE < mean_lower)? mean_lower: mean0MLE;  
   mean1MLE= (mean1MLE > mean_upper)? mean_upper: mean1MLE; 
   mean1MLE= (mean1MLE < mean_lower)? mean_lower : mean1MLE; 

   //estimate a0 a1
   //check if estimated a1 > a0 - if so switch means and var
   T a0MLE,a1MLE;
   if(counts1) a0MLE = counts0/(counts0+counts1);
   else a0MLE=1;
  
   if(a0MLE < .5){
			//swap terms 
			 swap(mean0,mean1);
			 swap(var0,var1);
			 swap(mean0MLE,mean1MLE);
			 swap(counts1,counts0);
			 swap(a0,a1);
			 a0MLE=1.0-a0;
		 }
		 //find unormalized counts
		 T ucounts,ucounts0,ucounts1;
		 if(normed){
				T totalCounts=counts0+counts1;
				ucounts=(T) nPoints * totalCounts;
				ucounts0=(T) nPoints * counts0;	
				ucounts1=(T) nPoints * counts1;
			}	
		 else{
				ucounts=counts0+counts1;
				ucounts0=counts0;
				ucounts1=counts1;
			}	
   
   a0MLE=(a0MLE > MIXTURE_RATIO + 2*MIXTURE_STD)?MIXTURE_RATIO + 2*MIXTURE_STD :a0MLE;
   a0MLE=(a0MLE < MIXTURE_RATIO - 2*MIXTURE_STD)?MIXTURE_RATIO - 2*MIXTURE_STD :a0MLE;
   a1MLE=1.0-a0MLE;
   T rawVar0MLE= sumsq0-mean0MLE*mean0MLE*counts0/(counts0+6.0);
   T rawVar1MLE= sumsq1-counts1-mean1MLE*mean1MLE*counts1/(counts1+6.0);
   
   if(isnan(rawVar0MLE)){
  		rawVar0MLE=(var_est+(kp*counts0)*(meanEst0-mean0MLE)*(meanEst0-mean0MLE))/(counts0+kp);	
			}
			else{
				rawVar0MLE+=(var_est+(kp*counts0)*(meanEst0-mean0MLE)*(meanEst0-mean0MLE))/(counts0+kp);	
			}
   if(isnan(rawVar1MLE)){
  		rawVar1MLE=(var_est+(kp*counts1)*(meanEst1-mean1MLE)*(meanEst1-mean1MLE))/(counts1+kp);	
			}
			else{
				rawVar1MLE+=(var_est+(kp*counts1)*(meanEst1-mean1MLE)*(meanEst1-mean1MLE))/(counts1+kp);	
			}				

   T var0MLE=rawVar0MLE;
   T var1MLE=rawVar1MLE;

  //enforce limits
   var0MLE=(var0MLE > var_upper)? var_upper : var0MLE;
   var0MLE=(var0MLE < var_lower)? var_lower : var0MLE;
   var1MLE=(var1MLE > var_upper)? var_upper : var1MLE;
   var1MLE=(var1MLE < var_lower)? var_lower : var1MLE;

   if(fabs(mean0MLE-mean0) > converge || fabs(mean1MLE-mean1) > converge){
	   mean0=mean0MLE;mean1=mean1MLE;
	   var0=var0MLE;var1=var1MLE;
	   a0=a0MLE;a1=a1MLE;	   
		 }
		 else{	
				loglikelihood=L;
				converged=1;
    var0=rawVar0MLE;var1=rawVar1MLE;
			 var0=(var0 > var_upper*2.0)? var_upper*2.0 : var0;
    var0=(var0 < var_lower/2.0)? var_lower/2.0 : var0;
    var1=(var1 > var_upper*2.0)? var_upper*2.0 : var1;
    var1=(var1< var_lower/2.0)? var_lower/2.0 : var1;    
    a1=1.0-a0;
     mean0=mean0MLE;mean1=mean1MLE;
    return(1);
			}
				//adjustment for binomial - seems to give worse results - enough to soft constrain p 
				//loglikelihood=L+(ucounts0+alpha1-1)*log(p_est)+(ucounts1+beta1-1)*log(1.0-p_est);				return(1);

		 if(++nIterations > maxIterations){
				loglikelihood=L;
				var0=rawVar0MLE;var1=rawVar1MLE;   
				
			 var0=(var0 > var_upper*2.0)? var_upper*2.0 : var0;
    var0=(var0 < var_lower/2.0)? var_lower/2.0 : var0;
    var1=(var1 > var_upper*2.0)? var_upper*2.0 : var1;
    var1=(var1< var_lower/2.0)? var_lower/2.0 : var1;
    a1=1.0-a0;
				mean0=mean0MLE;mean1=mean1MLE;
				return(0);           
			}
		}
	}  
  int optimize(T input_p){
   //optimize mean and use common variance
    //return loglikelihood after convergence with betabinomial term   
   const int maxIterations=MAX_EM_ITERATIONS;
   const double converge=EM_CONVERGE;
   bool done=0;
   int nIterations=0;
   const int last=(nBins-1<maxValue)? nBins-1 : maxValue;
   
  	//T alpha1=((1.0-p_est)/(p_sigma*p_sigma) -1.0/p_est)*p_est*p_est;
			//T beta1=alpha1*(1.0/(p_est-1));

   //fit mean stdevs to normals norma a*invsqrtpi/sigma * exp((-(x-mean)**2)/(2*var*var))
   //however the  invsqrtpi/sigma needs to be in natural log scale to properly normalize
   
   //adjust ln
   vector <int> map;
   vector <T> mylnx; //mapped copy of lnx
   vector <T> mydata;
   for (int j=minValue;j<=last; j++){
			//check for underflow
			 if(data[j] && !isnan(data[j])){
					map.push_back(j);
					mylnx.push_back(lnx[j]);
					mydata.push_back(data[j]);
				}	
			} 
   while(!done){
    double term0=log(a0)+LNINVSQRT2PI-log(sqrt(var0));
    double term1=log(a1)+LNINVSQRT2PI-log(sqrt(var1));
    double inv2Var0=1.0/(2.0*var0);
    double inv2Var1=1.0/(2.0*var1); 

   
  //calculate weighted percentage that points in j come from distribution 1 or 2
  //then calculate stdevs and mean for expectation maximization
   double L=0,sum0=0,sumsq0=0,sum1=0,sumsq1=0,counts0=0,counts1=0;
   const int setSize=map.size();
   for (int m=0;m<setSize; m++){
			//check for underflow
			 	T logLike0=term0-inv2Var0*(mylnx[m]-mean0)*(mylnx[m]-mean0);
			 	T logLike1=term1-inv2Var1*(mylnx[m]-mean1)*(mylnx[m]-mean1);

			 	//check for underflows
			 	if((!logLike1 || std::isnan(logLike1)) && logLike0 && !std::isnan(logLike0)){
		     L+=mydata[m]*logLike0;
			 		 counts0+=mydata[m];
			 		 sum0+=mylnx[m]*mydata[m];
			 		 sumsq0+=mylnx[m]*mylnx[m]*mydata[m];
			 	}
			 	else if	((!logLike0 || std::isnan(logLike0)) && logLike1 && !std::isnan(logLike1)){
		     L+=mydata[m]*logLike1;
			 		 counts1+=mydata[m];
			 		 sum1+=mylnx[m]*mydata[m];
			 		 sumsq1+=mylnx[m]*mylnx[m]*mydata[m];
			 	}
			 	else{
			 		const T logRatio=logLike0-logLike1;
			 	 if(logRatio > 15){
		     L+=mydata[m]*logLike0;
			 		 counts0+=mydata[m];
			 		 sum0+=mylnx[m]*mydata[m];
			 		 sumsq0+=mylnx[m]*mylnx[m]*mydata[m];	
						}
						else if (logRatio < -15){
		     L+=mydata[m]*logLike1;
			 		 counts1+=mydata[m];
			 		 sum1+=mylnx[m]*mydata[m];
			 		 sumsq1+=mylnx[m]*mylnx[m]*mydata[m];
						}
						else{
			 		 const T ratio=exp(logRatio);
			 		 const T r0=ratio/(ratio+1.0);
			 		 const T r1=1.0-r0;
			 	
					  T nPoints0=mydata[m]*r0;
				 	 T nPoints1=mydata[m]*r1;
				   counts0+=nPoints0;
				   counts1+=nPoints1;
				   L+=nPoints0*logLike0+nPoints1*logLike1;
				 	 sum0+=nPoints0*mylnx[m];
					  sumsq0+=nPoints0*mylnx[m]*mylnx[m];
				 	 sum1+=nPoints1*mylnx[m];
					  sumsq1+=nPoints1*mylnx[m]*mylnx[m];
					 }
					}		
			}
		 //calculate estimate of mean
		 T mean0MLE=sum0/(counts0); 
   T mean1MLE=sum1/(counts1);
	
   mean0MLE= (mean0MLE > mean_upper)? mean_upper : mean0MLE;  
   mean0MLE= (mean0MLE < mean_lower)? mean_lower: mean0MLE;  
   mean1MLE= (mean1MLE > mean_upper)? mean_upper: mean1MLE; 
   mean1MLE= (mean1MLE < mean_lower)? mean_lower : mean1MLE; 
 
   //estimate a0 a1
   //check if estimated a1 > a0 - if so switch means and var
   T a0MLE,a1MLE;
   if(counts1) a0MLE = counts0/(counts0+counts1);
   else a0MLE=1;
  
   if(a0MLE < .5){
			//swap terms 
			 swap(mean0,mean1);
			 swap(var0,var1);
			 swap(mean0MLE,mean1MLE);
			 swap(counts1,counts0);
			 swap(a0,a1);
			 a0MLE=1.0-a0;
		 }
		 //find unormalized counts
		 T ucounts,ucounts0,ucounts1;
		 if(normed){
				T totalCounts=counts0+counts1;
				ucounts=(T) nPoints * totalCounts;
				ucounts0=(T) nPoints * counts0;	
				ucounts1=(T) nPoints * counts1;
			}	
		 else{
				ucounts=counts0+counts1;
				ucounts0=counts0;
				ucounts1=counts1;
			}	
   
   a0MLE=(a0MLE > MIXTURE_RATIO + 2*MIXTURE_STD)?MIXTURE_RATIO + 2*MIXTURE_STD :a0MLE;
   a0MLE=(a0MLE < MIXTURE_RATIO - 2*MIXTURE_STD)?MIXTURE_RATIO - 2*MIXTURE_STD :a0MLE;
   a1MLE=1.0-a0MLE;
   T rawVar0MLE=(sumsq0/(counts0)-mean0MLE*mean0MLE)*.6;
   T rawVar1MLE=(sumsq1/(counts1)-mean1MLE*mean1MLE)*.6;
   if(isnan( rawVar0MLE))rawVar0MLE=var_lower/2.0;  
   if(isnan( rawVar1MLE))rawVar1MLE=var_lower/2.0;
   T var0MLE=rawVar0MLE;
   T var1MLE=rawVar1MLE;

  //enforce limits
   var0MLE=(var0MLE > var_upper)? var_upper : var0MLE;
   var0MLE=(var0MLE < var_lower)? var_lower : var0MLE;
   var1MLE=(var1MLE > var_upper)? var_upper : var1MLE;
   var1MLE=(var1MLE < var_lower)? var_lower : var1MLE;

   if(fabs(mean0MLE-mean0) > converge || fabs(mean1MLE-mean1) > converge){
	   mean0=mean0MLE;mean1=mean1MLE;
	   var0=var0MLE;var1=var1MLE;
	   a0=a0MLE;a1=a1MLE;	   
		 }
		 else{	
				loglikelihood=L;
				converged=1;
    var0=rawVar0MLE;var1=rawVar1MLE;
			 var0=(var0 > var_upper*2.0)? var_upper*2.0 : var0;
    var0=(var0 < var_lower/2.0)? var_lower/2.0 : var0;
    var1=(var1 > var_upper*2.0)? var_upper*2.0 : var1;
    var1=(var1< var_lower/2.0)? var_lower/2.0 : var1;    
    a1=1.0-a0;
     mean0=mean0MLE;mean1=mean1MLE;
    //fprintf(stderr,"%f %f %f %f %f %f\n",mean0MLE*INVLOG2,mean1MLE*INVLOG2,var0MLE*INVLOG2*INVLOG2,var1MLE*INVLOG2*INVLOG2,a0MLE,a1MLE);

    return(1);
			}
				//adjustment for binomial - seems to give worse results - enough to soft constrain p 
				//loglikelihood=L+(ucounts0+alpha1-1)*log(p_est)+(ucounts1+beta1-1)*log(1.0-p_est);				return(1);

		 if(++nIterations > maxIterations){
				loglikelihood=L;
				var0=rawVar0MLE;var1=rawVar1MLE;   
				
			 var0=(var0 > var_upper*2.0)? var_upper*2.0 : var0;
    var0=(var0 < var_lower/2.0)? var_lower/2.0 : var0;
    var1=(var1 > var_upper*2.0)? var_upper*2.0 : var1;
    var1=(var1< var_lower/2.0)? var_lower/2.0 : var1;
    a1=1.0-a0;
				mean0=mean0MLE;mean1=mean1MLE;
			 //fprintf(stderr,"%f %f %f %f %f %f\n",mean0MLE*INVLOG2,mean1MLE*INVLOG2,var0MLE*INVLOG2*INVLOG2,var1MLE*INVLOG2*INVLOG2,a0MLE,a1MLE);

				return(0);           
			}
		}
	}
};	

void init_logx(){
	logx[0]=0;
	lnx[0]=0;
 for (int i=1;i<32768;i++){ 
		lnx[i]=log((double)i);  
  logx[i]=lnx[i]*INVLOG2;
	}
}
template <class T> void transferData(int color,DeconRecord<T> &deconRecord,EM_GMM<T> &EM){
	
	deconRecord.means0[color]=(EM.mean0*INVLOG2 < MIN_NOISE)? MIN_NOISE : EM.mean0*INVLOG2;
	deconRecord.means1[color]=(EM.mean1*INVLOG2 < MIN_NOISE)? MIN_NOISE : EM.mean1*INVLOG2;
	deconRecord.meansEst0[color]=EM.meanEst0*INVLOG2;
	deconRecord.meansEst1[color]=EM.meanEst1*INVLOG2;
	deconRecord.stdevs0[color]=sqrt(EM.var0)*INVLOG2;
	deconRecord.stdevs1[color]=sqrt(EM.var1)*INVLOG2;
	deconRecord.a0[color]=EM.a0;
	deconRecord.a1[color]=EM.a1;			
	deconRecord.loglikelihoods[color]=-EM.loglikelihood;
	deconRecord.dataSizes[color]=EM.nPoints;
	deconRecord.converged[color]=EM.converged;
	if((EM.meanEst0-EM.meanEst1)*(EM.mean0-EM.mean1) <0) deconRecord.flipped[color]=1;
	else	 deconRecord.flipped[color]=0;
}

template <typename T1, typename T2> class wellRecord{
	//qnorm is done by well
	//these are then combined
	public:
	 vector <T1> values[NCOLORS1OFF]; //original values in int16
	 vector <string> groupNames;
	 vector <T2> offsets[NCOLORS1OFF];
	 vector <int> nBeads;
	 vector <int> groupOffsets; //the well where the group starts
	 int nTotalWells;
	 
	 int minValue; //minimum value for peak 
	 int maxValue; //maximum value for peak
	 float peakGridSize; //gridSize for peakSearch
	 float smoothHalfWindowSize; //half window size for smoothing
	 wellRecord(string list,float minLogValue, float maxLogValue,float _smoothHalfWindowSize,float _peakGridSize){
	  int16_t header[2];
		 int16_t v;
		 nTotalWells=0;
		 minValue=pow(2.0,minLogValue); //round down bin
		 maxValue=pow(2.0,maxLogValue)+1;
		 if(maxValue > 32767){
				maxValue=32767;
			}
			smoothHalfWindowSize=_smoothHalfWindowSize;
			peakGridSize=_peakGridSize;	  
		 FILE *fp=fopen(list.c_str(),"r");
		 //can have very long lines
		 fseek(fp,0,SEEK_END);
		 int size=ftell(fp)+1;
		 fseek(fp,0,SEEK_SET);
		 char *line=new char[size];
		 //cerr << "allocated " << size << " bytes"<<endl;
		 if(!line){
				cerr << "unable to allocate " << size << " bytes" <<endl;
				exit(0);
			}	
		 int n=0;
		 for(int i=1;i<NCOLORS1OFF;i++){
		 	offsets[i].push_back(0);
		 }
		 groupOffsets.push_back(0);
   while(fgets(line,size,fp) != NULL){
				//format is name\tfileName1\tfileName2...fileNamen\n
				vector <string> fileNames;
				string groupName;
				const char delimit[]=" \t\r\n\v\f";
				char *pch = strtok (line,delimit);
			 groupName=(pch);
			 pch = strtok (NULL,delimit); 
			 //start extracting list of fileNames
			 
			 //check for case when there is only one element per line in which case the groupName is also the file
			 if(pch == NULL){
					fileNames.push_back(groupName);
				}
				else{	  
				 while (pch != NULL){
      fileNames.push_back(pch);
      pch = strtok (NULL,delimit);  										
				 }
				}
				//now extract the information from the filenames
				int nFilesRead=0;
				for (int k=0;k<fileNames.size();k++){
				//	fprintf(stderr,"opening %s\n",fileNames[k].c_str());
			  FILE *fp1=fopen(fileNames[k].c_str(),"r");
			  if(fp1 == NULL){
			   fprintf(stderr,"unable to open %s\n",fileNames[k].c_str());
			   continue;
			  }
			  while(fread(&header,2*sizeof(int16_t),1,fp1)){
			 	 const int id=header[0];
      for(int k=0;k<header[1];k++){
		     if(fread(&v,sizeof(int16_t),1,fp1)){
				  	 if(v>=0){values[id].push_back(v);}	
				 	 }
			   }
			  }
		   fclose(fp1);
		   //save sizes to offset
		   int totalRead=0;
		   for(int i=1;i<NCOLORS1OFF;i++){
			  	offsets[i].push_back(values[i].size());
			  	totalRead+=values[i].size()-offsets[i][nBeads.size()];
			  }	
		   nBeads.push_back(totalRead);			  
					nFilesRead++;	
			 }
			 if(nFilesRead){
				 groupOffsets.push_back(groupOffsets[groupOffsets.size()-1]+nFilesRead);
				 groupNames.push_back(groupName);
				 nTotalWells+=nFilesRead;	
				}
		 }
		 fprintf(stderr,"closing file\n");
		 fclose(fp);
		 delete[] line;
	 }	
  wellRecord(const wellRecord &A): values(A.values),groupNames(A.groupNames),groupOffsets(A.groupOffsets),offsets(A.offsets),maxValue(A.maxValue),minValue(A.minValue),smoothHalfWindowSize(A.smoothHalfWindowSize),peakGridSize(A.peakGridSize),nTotalWells(A.nTotalWells),nBeads(A.nBeads){}
  
	void dumpEnsemble(string outputFile){
		cerr << "dumping density to " << outputFile <<endl;
		FILE *fp=fopen(outputFile.c_str(),"w");
		int counts[32768];
		fill(counts,counts+32768,0);		
		//the 0th slice can be used to store errors in the old PERL script - keep it in case this is to be implemented later in C++
		fwrite(counts,sizeof(int),32768,fp);
		for(int id=1;id<NCOLORS1OFF;id++){
			fill(counts,counts+32768,0); 
			for(int k=0;k<values[id].size();k++){
				counts[values[id][k]]++;
			}
			fwrite(counts,sizeof(int),32768,fp);
		}
		fclose(fp);
	}
	void dumpWells(string outputDir){
		cerr << "dumping density to directory " << outputDir <<endl;
		for(int n=0;n<nTotalWells;n++){
			ostringstream nstream;
			nstream << "." << n;
			
			string outputFile=outputDir+"/"+groupNames[n]+nstream.str()+".den";	
	 	FILE *fp;
	 	if(fp=fopen(outputFile.c_str(),"w")){
	 	 int counts[32768];
	 	 fill(counts,counts+32768,0);
		  //the 0th slice can be used to store errors in the old PERL script - keep it in case this is to be implemented later in C++
    fwrite(counts,sizeof(int),32768,fp);		
		  for(int id=1;id<NCOLORS1OFF;id++){
					fill(counts,counts+32768,0); 
			  for(int k=offsets[id][n];k<offsets[id][n+1];k++){
			  	counts[values[id][k]]++;
			  }	
		  	fwrite(counts,sizeof(int),32768,fp);
		  }
		  fclose(fp);
	  }
		}
	}	
	void dumpSingleWell(string outputFile,int slice){
		cerr << "dumping density of well " << slice << " to file " << outputFile <<endl;
	 FILE *fp;
	 if(fp=fopen(outputFile.c_str(),"w")){
	  int counts[32768];
	  fill(counts,counts+32768,0);
		 //the 0th slice can be used to store errors in the old PERL script - keep it in case this is to be implemented later in C++
   fwrite(counts,sizeof(int),32768,fp);		
		 for(int id=1;id<NCOLORS1OFF;id++){ 
				fill(counts,counts+32768,0);
			 for(int k=offsets[id][slice];k<offsets[id][slice+1];k++){
			  counts[values[id][k]]++;
			 }	
		  fwrite(counts,sizeof(int),32768,fp);
		 }
		 fclose(fp);
	 }
	}
	void dumpDensity(string outputFile,float densitySmoothHalfWindowSize){
	 FILE *fp;
	 if(fp=fopen(outputFile.c_str(),"w")){
	 	int counts[32768];
	 	fill(counts,counts+32768,0);
	 	float smooth[32768];
	 	fill(smooth,smooth+32768,0);
   fwrite(smooth,1,sizeof(float)*32768,fp);//0 well - keep it for errors
		 int upper[32768];
		 int lower [32768];
		 memset (lower,0,32768*sizeof(int));
		 memset (upper,0,32768*sizeof(int));
		 const float upperConstant=pow(2,densitySmoothHalfWindowSize);
		 const float lowerConstant=1.0/upperConstant;
		 for (int i=0;i<32768;i++){
		 	lower[i]= (int) (((float)i)*lowerConstant+.5);
    upper[i]= (int) (((float)i)*upperConstant+.5);	
    if(upper[i] > 32767)upper[i]=32767;
	  }
		 for (int i=1; i<=NCOLORS; i++){
			 fill(smooth,smooth+32768,0);	
			 int nBins=logWindowSmooth(lower,upper,smooth,&(values[i][0]),values[i].size());
		 	fwrite(smooth,1,32768*sizeof(float),fp);
			}
		 fclose(fp);
		}
	}
template<class Ta, class Tb> int logWindowSmooth(int *lower,int *upper,Ta *smooth,Tb *colorValues,int nValues){
		int maxBin=0;
		int counts[32768];
		fill(counts,counts+32768,0);	
		double totalWeight=0;
		for(int i=0;i<nValues;i++){
		 counts[colorValues[i]]++;
	 }
		for(int i=1;i<32768;i++){
		 if(counts[i])maxBin=i;	
	 }		
		for(int i=1;i<=maxBin;i++){
   int sum=0;
   //const int limit=(upper[i]<maxBin)? maxBin:upper[i];
	 	for(int k=lower[i];k<=upper[i];k++){
				sum+=counts[k];
			}
			smooth[i]=(Ta)sum/(Ta)(upper[i]-lower[i]+1);
			//fprintf(stderr,"%d %d %d %e\n",i,lower[i],upper[i],smooth[i]);	
		 totalWeight+=smooth[i];
		}
		double invTotalWeight=1.0/totalWeight;
		for(int i=0;i<=maxBin;i++){
		 smooth[i]*=invTotalWeight;
	 }
		return(maxBin+1);
	}
	template <class T> double dampenLowNoise(T *smooth,int nBins){
		//version for a single file in wells
	 double lowNoiseFraction=0;
  for(int i=0;i<minValue;i++){
  	const double factor= exp(-0.3*(minValue-i));
  	lowNoiseFraction +=(1.0-factor)*smooth[i];
  	smooth[i]*=factor;
		}	
  //renormalize
  double normFactor=1.0/(1.0-lowNoiseFraction);
  for(int i=minValue;i<nBins;i++)
   smooth[i]*=normFactor; 
	 return (lowNoiseFraction);
	}
	template <class T> void findPeaks(T *smooth, int nBins, T gridSize,vector<int> &peakBins, vector<T> &peakValues, vector<T> &peakMeans){
		//finds peaks using a interval (grid) search
		//uses equal log intervals and logSmoothed data
		//an interval contains a peak if the previous interval and and the following interval have lower mean values
		//the bin with the maxCounts in that interval is the location of a peak
		//peaks are then sorted by mean values (max values is more vulnerable to spikes)
		//a good value for grid size is roughly the std deviation of the expected peaks
		
  vector<T> maxValues,meanValues;
  vector<int>maxBins;
  if(nBins <4){
			//edge case - return the middle bin
			peakBins.resize(1);
		 peakMeans.resize(1);
		 peakValues.resize(1);
			peakBins[0]=nBins/2;
			peakMeans[0]=smooth[nBins/2];
			peakValues[0]=smooth[nBins/2];
			return;
		}
	 
	 int minBin=pow(2,minValue);int maxBin=(pow(2,maxValue) < nBins)?pow(2,maxValue) : nBins;
	 T gridFactor=pow(2,gridSize);
	 

	 int i=minBin;
		while(i<nBins){
			int upperBin=(int) ((T)i*gridFactor + 0.5);
	  if (upperBin >nBins-1) upperBin=nBins-1;
		 T maxValue=smooth[i];
	  int maxBin=i;
   double sum=0;	 
   for (int k=i+1;k<=upperBin;k++){
		 	if(smooth[k] > maxValue){
				 maxValue=smooth[k];
				 maxBin=k;
			 }
 		 sum+=smooth[k];		 
		 }
		 maxValues.push_back(maxValue);
		 maxBins.push_back(maxBin);
		 //meanValues.push_back((T)sum/(T)(upperBin-i+1));	
		 meanValues.push_back((T)sum);			 		 
		 i=upperBin+1;
		}	 
		vector <T> myPeakValues,myPeakBins,myPeakMeans;
	 for  (int i=2;i<maxBins.size()-2;i++){
			if(maxBins[i-2] > minValue && maxBins[i+1] < maxValue){ //-2 necessary otherwise there is an edge effect
 		 if(meanValues[i] > meanValues[i-1] && meanValues[i] > meanValues[i+1]){	
		  	myPeakBins.push_back(maxBins[i]);
			  myPeakValues.push_back(maxValues[i]);	
			  myPeakMeans.push_back(meanValues[i]);
		  }
			}
		}
		//check for empty case - no peaks - seed with middle
		if(!myPeakBins.size()){
			for  (int i=maxBins.size()/3;i<maxBins.size();i+=maxBins.size()/2){
			 myPeakBins.push_back(maxBins[i]);
			 myPeakValues.push_back(maxValues[i]);	
			 myPeakMeans.push_back(meanValues[i]);
			}
		}	
		 
		vector<int> index;
		index.resize(myPeakMeans.size());
		peakBins.resize(1);
		peakMeans.resize(1);
		peakValues.resize(1);
		sort_by_scores(myPeakMeans.size(),myPeakMeans.data(),index.data(),0);
				 
		peakBins[0]=myPeakBins[index[0]];
		peakValues[0]=myPeakValues[index[0]];
		peakMeans[0]=myPeakMeans[index[0]];
	 for (int i=1;i<index.size();i++){
			const int s=index[i];
			if(myPeakMeans[s] < 0.2*myPeakMeans[0])break; 
		 peakBins.push_back(myPeakBins[s]);
		 peakValues.push_back(myPeakValues[s]);
		 peakMeans.push_back(myPeakMeans[s]);
		}	 
	}
	template <class T> void findPeaks(int *counts, int nBins, T gridSize,vector<int> &peakBins, vector<T> &peakValues, vector<T> &peakMeans){
		//finds peaks using a interval (grid) search
		//uses equal log intervals
		//an interval contains a peak if the previous interval and and the following interval have lower mean values
		//the bin with the maxCounts in that interval is the location of a peak
		//peaks are then sorted by mean values (max values is more vulnerable to spikes)
		//a good value for grid size is roughly the std deviation of the expected peaks
		
  vector<T> maxValues,meanValues;
  vector<int>maxBins;
  if(nBins <4){
			//edge case - return the middle bin
			peakBins.resize(1);
		 peakMeans.resize(1);
		 peakValues.resize(1);
			peakBins[0]=nBins/2;
			peakMeans[0]=counts[nBins/2];
			peakValues[0]=counts[nBins/2];
			return;
		}
	 
	 
	 int minBin=pow(2,minValue);int maxBin=(pow(2,maxValue) < nBins)?pow(2,maxValue) : nBins;
	 T gridFactor=pow(2,gridSize);
	 

	 int i=minBin;
		while(i<nBins){
			int upperBin=(int) ((T)i*gridFactor + 0.5);
	  if (upperBin >nBins-1) upperBin=nBins-1;
		 T maxValue=counts[i];
	  int maxBin=i;
   double sum=0;	 
   for (int k=i+1;k<=upperBin;k++){
		 	if(counts[k] > maxValue){
				 maxValue=counts[k];
				 maxBin=k;
			 }
 		 sum+=counts[k];		 
		 }
		 maxValues.push_back(maxValue);
		 maxBins.push_back(maxBin);
		 //meanValues.push_back((T)sum/(T)(upperBin-i+1));	
		 meanValues.push_back((T)sum);			 		 
		 i=upperBin+1;
		}	 
		vector <T> myPeakValues,myPeakBins,myPeakMeans;
	 for  (int i=2;i<maxBins.size()-2;i++){
			if(maxBins[i-2] > minValue && maxBins[i+1] < maxValue){ //-2 necessary otherwise there is an edge effect
 		 if(meanValues[i] > meanValues[i-1] && meanValues[i] > meanValues[i+1]){	
		  	myPeakBins.push_back(maxBins[i]);
			  myPeakValues.push_back(maxValues[i]);	
			  myPeakMeans.push_back(meanValues[i]);
		  }
			}
		}
		//check for empty case - no peaks - seed with middle
		if(!myPeakBins.size()){
			for  (int i=maxBins.size()/3;i<maxBins.size();i+=maxBins.size()/2){
			 myPeakBins.push_back(maxBins[i]);
			 myPeakValues.push_back(maxValues[i]);	
			 myPeakMeans.push_back(meanValues[i]);
			}
		}	
		 
		vector<int> index;
		index.resize(myPeakMeans.size());
		peakBins.resize(1);
		peakMeans.resize(1);
		peakValues.resize(1);
		sort_by_scores(myPeakMeans.size(),myPeakMeans.data(),index.data(),0);
				 
		peakBins[0]=myPeakBins[index[0]];
		peakValues[0]=myPeakValues[index[0]];
		peakMeans[0]=myPeakMeans[index[0]];
	 for (int i=1;i<index.size();i++){
			const int s=index[i];
			if(myPeakMeans[s] < 0.2*myPeakMeans[0])break; 
		 peakBins.push_back(myPeakBins[s]);
		 peakValues.push_back(myPeakValues[s]);
		 peakMeans.push_back(myPeakMeans[s]);
		}	 
	}
	template <class T> int deconvolute(DeconRecord<T> &deconRecord,FILE *fpRaw,FILE *fpSmooth){
	 if(!logx[2]){init_logx();}
		T smooth[32768];
		int counts[32768];
		fill(smooth,smooth+32768,0);
		fill(counts,counts+32768,0);
		if(fpSmooth) fwrite(smooth,1,sizeof(T)*32768,fpSmooth); 
	 if(fpRaw) fwrite(counts,1,sizeof(T)*32768,fpRaw);
		int upper[32768];
		int lower [32768];
		memset (lower,0,32768*sizeof(int));
		memset (upper,0,32768*sizeof(int));
		const T upperConstant=pow(2,smoothHalfWindowSize);
		const T lowerConstant=1.0/upperConstant;
		for (int i=0;i<32768;i++){
			lower[i]= (int) (((T)i)*lowerConstant+.5);
   upper[i]= (int) (((T)i)*upperConstant+.5);	
   if(upper[i] > 32767)upper[i]=32767;
	 }
		for (int i=1; i<=NCOLORS; i++){
			fill(smooth,smooth+32768,0);
			vector<T> peakValues,peakMeans;
			vector<int>peakBins;			
			int nBins=logWindowSmooth(lower,upper,smooth,&(values[i][0]),values[i].size());
			double lowNoiseFraction=dampenLowNoise(smooth,nBins);
			if(fpSmooth)fwrite(smooth,1,32768*sizeof(T),fpSmooth);
			fill(counts,counts+32768,0);
			for(int k=0;k<values[i].size();k++)
			 counts[values[i][k]]++;
			 
			if(fpRaw){
    fwrite(counts,1,32768*sizeof(T),fpRaw);
			}
   findPeaks<T>(smooth,nBins,0.2,peakBins,peakValues,peakMeans);
   if(i >10){
				T nPoints=(double)values[i].size()-lowNoiseFraction*(double)values[i].size();
			 for(int k=0;k<32768;k++){
	 			smooth[k]=counts[k];
 			} 
    fitLogNormEM(smooth,0,nBins,nPoints,peakBins,deconRecord,i,1,0);
    deconRecord.noiseFractions[i]=lowNoiseFraction;
			}
		}
	}
	template <class T> int deconvolute(int nFold,int n,int nGroups,vector<DeconRecord<T>> &deconRecords,vector<MedianRecord<T>>&medianRecords,FILE *fpRaw){
		return deconvolute_peakSeed(nFold,n,nGroups,deconRecords,medianRecords,fpRaw);
		//return deconvolute_peakSeed(n,nGroups,deconRecords[0],fpRaw);
	}	
 template <class T> int deconvolute_peakSeed(int n,int nGroups,DeconRecord<T> &deconRecord,FILE *fpRaw){
	 if(!logx[2]){init_logx();}
		int counts[32768];
		fill(counts,counts+32768,0);
	 if (fpRaw) fwrite(counts,1,sizeof(T)*32768,fpRaw);
	 int finalGroupOffset= (n+nGroups >= groupOffsets.size())? groupOffsets[groupOffsets.size()-1] : groupOffsets[n+nGroups];
	 //fprintf(stderr,"%d start %d final offset %d\n",n,groupOffsets[n],finalGroupOffset);
		for (int i=1; i<=NCOLORS; i++){
			const int start=offsets[i][groupOffsets[n]];
			const int finish= offsets[i][finalGroupOffset];
			int nBins=0;
			T nPoints=finish-start;
			if(nPoints >= MINBEADS){ 
			 vector<int>peakBins;
			 peakBins.resize(2);		
			 
			 fill(counts,counts+32768,0);
			 for(int k=start;k<finish;k++){
					if(values[i][k]>nBins) nBins=values[i][k]+1;
			  counts[values[i][k]]++;
				}
			 if(fpRaw) fwrite(counts,1,32768*sizeof(int),fpRaw);
    peakBins[0]=pow(2,8.1);
    peakBins[1]=pow(2,7.9);		 
    if(i >10){
					T tcounts[32768];
					fill(tcounts,tcounts+32768,0);
					int k=0;
					while (k<minValue){
						nPoints-=counts[k];
						//tcounts[k]=counts[k];
						k++;
					}
					while (k<32768){
	 			 tcounts[k]=counts[k];
	 			 k++;
 			 }
     fitLogNormEM(tcounts,0,nBins,nPoints,peakBins,deconRecord,i,1,1);
     deconRecord.noiseFractions[i]=1.0-(nPoints/(T)(finish-start));
				}
			}
			else{
			 deconRecord.dataSizes[i]=finish-start;
				if(fpRaw){
		   fill(counts,counts+32768,0);
		   for(int k=start;k<finish;k++)
			   counts[values[i][k]]++;
	    fwrite(counts,1,sizeof(int)*32768,fpRaw);
			 }				
			}	
		}
		return(1);
	}
	template <class T> int deconvolute_peakSeed(int nFold,int n,int nGroups,vector<DeconRecord<T>> &deconRecords,vector<MedianRecord<T>>&medianRecords,FILE *fpRaw){
  //this overload takes 1/10 slices and finds the median value
	 if(!logx[2]){init_logx();}
		int counts[32768];
	 if (fpRaw) fwrite(counts,1,sizeof(T)*32768,fpRaw);
	 int finalGroupOffset= (n+nGroups >= groupOffsets.size())? groupOffsets[groupOffsets.size()-1] : groupOffsets[n+nGroups];
	 //fprintf(stderr,"%d start %d final offset %d\n",n,groupOffsets[n],finalGroupOffset);
		for (int i=1; i<=NCOLORS; i++){
			const int start=offsets[i][groupOffsets[n]];
			const int finish= offsets[i][finalGroupOffset];			
			if(i >10){
			 int nBins=0;
			 int totalPoints=finish-start;
			 if(totalPoints >= MINBEADS){ 
			  vector<int>peakBins(2);
			  vector<unsigned int>foldIndices=getFolds(totalPoints,nFold);; 
			  for(int f=0;f<nFold;f++){
						int nPoints=0;
						fill(counts,counts+32768,0);
						unsigned int j=0;
		  	 for(int k=start;k<finish;k++){
							if(nFold >1 && foldIndices[j++] == f ) continue;
		  			if(values[i][k]>nBins) nBins=values[i][k]+1;
		 	   counts[values[i][k]]++;
		 	   nPoints++;
		 	 	}
		 	  if(fpRaw) fwrite(counts,1,32768*sizeof(int),fpRaw);
      peakBins[0]=pow(2,8.1);
      peakBins[1]=pow(2,7.9);		 

			 		T tcounts[32768];
			 		fill(tcounts,tcounts+32768,0);
			 		int k=0;
			 		while (k<minValue){
				 		nPoints-=counts[k];
			 			//tcounts[k]=counts[k];
				 		k++;
				 	}
				 	while (k<32768){
	 		 	 tcounts[k]=counts[k];
	 		 	 k++;
 			  }
      fitLogNormEM(tcounts,0,nBins,nPoints,peakBins,deconRecords[f],i,1,1);
      medianRecords[i].addValues(deconRecords[f].means0[i], deconRecords[f].means1[i]);
      deconRecords[f].noiseFractions[i]=1.0-(nPoints/(T)(totalPoints));
				 }
				}
			}
			else{
				//write a blank Deconrecord with nPoints
				for(int f=0;f<nFold;f++){
				 deconRecords[f].dataSizes[i]=finish-start;
				}
				if(fpRaw){
		   fill(counts,counts+32768,0);
		   for(int k=start;k<finish;k++)
			   counts[values[i][k]]++;
	    fwrite(counts,1,sizeof(int)*32768,fpRaw);
			 }				
			}
		}
		return(1);
	}
 template <class T> int fitLogNormEM(T *smooth,bool normed, int nBins,int nPoints,vector<int> &peakBins,DeconRecord<T> &deconRecord,int color,bool isMode,bool seed){
  const T s0=PEAKWIDTH;
  const T s1=PEAKWIDTH;
  T sigmasqs[2]={s0*s0,s1*s1};
  vector<T> mus;
  mus.resize(peakBins.size());
  if(isMode){
   for(int i=0;i<peakBins.size();i++)
    mus[i]=logx[peakBins[i]]-s0*s0; 
		}
		else{
	  for(int i=0;i<peakBins.size();i++)
    mus[i]=logx[peakBins[i]]; 		
		}	
  //starting for search is both lognormals at the estimated medians
  
  //keep record of initial mode estimate
		deconRecord.meansEst0[color]=logx[peakBins[0]];
		deconRecord.meansEst1[color]=logx[peakBins[1]];
		T input_mus[2]={mus[0]+ (T).05,mus[0]-(T).05};
		EM_GMM <float> EM (smooth,normed,nBins,nPoints,minValue,maxValue,input_mus,sigmasqs,.5, MIXTURE_RATIO, 0.15,0.7,.275,PEAKWIDTH,.355); 
		if (seed){
			EM.reset(mus[0],mus[1],sigmasqs[0],sigmasqs[1],MIXTURE_RATIO);       			
			EM.seed_mclust(0.5);
			mus.resize(2);
			if(isnan(EM.mean0)){
				mus[0]=EM.mean1*INVLOG2+.05;
			 mus[1]=EM.mean1*INVLOG2-.05;
			}
			else if(isnan(EM.mean1)){
			 mus[0]=EM.mean0*INVLOG2+.05;
			 mus[1]=EM.mean0*INVLOG2-.05;
			}
			else{
				mus[0]=EM.mean0*INVLOG2+.05;
			 mus[1]=EM.mean1*INVLOG2-.05;	
			}
			EM.reset(mus[0]+.05,mus[1]-.05,sigmasqs[0],sigmasqs[1],MIXTURE_RATIO); 
		}
		
  int converged=EM.optimize(MIXTURE_RATIO);
  transferData(color,deconRecord,EM);
  double best_likelihood=EM.loglikelihood;

  for(int j=0;j<mus.size()-1;j++){
			for (int k=j+1;k<mus.size();k++){ 
				EM.reset(mus[j],mus[k],sigmasqs[0],sigmasqs[1],MIXTURE_RATIO); 
				int converged=EM.optimize(MIXTURE_RATIO);
				double likelihood=EM.loglikelihood;   
				if(likelihood < best_likelihood){
					best_likelihood=likelihood;
				 transferData(color,deconRecord,EM);
				 deconRecord.meansEst0[color]=mus[j];deconRecord.meansEst1[color]=mus[k];
				}
				EM.reset(mus[k],mus[j],sigmasqs[1],sigmasqs[0],MIXTURE_RATIO); 
				converged=EM.optimize(MIXTURE_RATIO);
				likelihood=EM.loglikelihood;   
				if(likelihood < best_likelihood){
					best_likelihood=likelihood;
				 transferData(color,deconRecord,EM);
				 deconRecord.meansEst0[color]=mus[j];deconRecord.meansEst1[color]=mus[k];
				}					
			}	 
		}
 }	  
	vector<unsigned int> getFolds(unsigned int size, unsigned int nFold){
		//for assigning nfolds
		vector<unsigned int> shuffledIndices(size);
  std::random_device rd;
  std::mt19937 rng(rd());
		for(unsigned int i=0;i<size;i++){
			shuffledIndices[i]=i%nFold;
		}
		for	(unsigned int i=0;i<size;i++){
			std::uniform_int_distribution<int> uni(i,size-1);
			unsigned int j=uni(rng);
 	 swap(shuffledIndices[i],shuffledIndices[j]);
		}
		return shuffledIndices;
	}		
};

