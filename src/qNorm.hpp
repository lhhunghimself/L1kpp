#define INVLOG2 1.44269504089

using namespace std;

class density{
	public:
	float counts [32768]; //counts
	float ccounts [32768];//cumulative counts
	int maxBin;
	
	//calculate density from well record
	density(){
	 maxBin=0;
	}
	//copy constructor
	density(const density  &A){
	 maxBin=A.maxBin;
  if(maxBin){
			memmove(counts,A.counts,32768*sizeof(float));
			memmove(ccounts,A.ccounts,32768*sizeof(float));
		}  		 
	}
	density & operator = (const density  &rhs){
	 maxBin=rhs.maxBin;
  if(maxBin){
			memmove(counts,rhs.counts,32768*sizeof(float));
			memmove(ccounts,rhs.ccounts,32768*sizeof(float));
		}  		 
	}
	density(wellRecord<int16_t,int> wells){
		//generates an averaged density file
		fill(counts,counts+32768,0.0f);
		fill(ccounts,ccounts+32768,0.0f);
		
		int *wellcCounts=new int[32768*wells.nTotalWells];
		memset(wellcCounts,0,32768*wells.nTotalWells*sizeof(int));
		maxBin=0;
		
		//define quantiles using the maxCounts in the different wells
		//go through each of these values and map the counts using the weighted mean value of the value of the well 
		//basically using histogram as a sort
		int maxCounts=wells.nBeads[0];
		int maxCountIndex=0;
		vector <float> ratios;
		vector <float> currentCounts; //keep track of currentCounts for each of the samples - needed for averaging the distribution
		vector <int> currentBins; 
		ratios.resize(wells.nTotalWells);
		currentCounts.resize(wells.nTotalWells);
		currentBins.resize(wells.nTotalWells);
		for (int i=0;i<wells.nTotalWells;i++){
		 if(maxCounts < wells.nBeads[i]){
				maxCounts=wells.nBeads[i];
				maxCountIndex=i;
			}	
		}
		
		float totalRatio=0;//normalization factor
		for (int i=0;i<wells.nTotalWells;i++){
		 ratios[i]=(float)wells.nBeads[i]/(float)maxCounts;
	  totalRatio+=ratios[i];
		}
		//put all the wells in their own histogram
		for(int i=1;i<=NCOLORS;i++){
			for(int n=0;n<wells.nTotalWells;n++){
				int *slice=wellcCounts+n*32768;
		 	for (int j=wells.offsets[i][n];j<wells.offsets[i][n+1];j++){
				 slice[wells.values[i][j]]++;
				}
			}	
		}
		//convert to cumulative counts and initialize currentCounts and currentBins
	 for(int n=0;n<wells.nTotalWells;n++){
		 int *slice=wellcCounts+n*32768;
		 sumCounts(slice);
		 currentCounts[n]=0;
		 currentBins[n]=0;
		}
		//combine histograms 
		//start from 0..to maxcounts
		//find equivalent bin for each sample	
		for(int i=0;i<maxCounts;i++){
			float sum=0;
			for(int n=0;n<wells.nTotalWells;n++){
				int *rawcCounts=wellcCounts+n*32768;
				currentCounts[n]+=ratios[n];
				//round Counts
				int currentCountsInt=(int) currentCounts[n]+.5;
				if (currentCountsInt > wells.nBeads[n])currentCountsInt=wells.nBeads[n]; //check for roundoff error NB otherwise this can cause search to run to maxValue
				while(rawcCounts[currentBins[n]]< currentCountsInt ){
					currentBins[n]++;
				}	
	   sum+=ratios[n]*currentBins[n];
			}
			
			float fBin=sum/totalRatio;
		 int lowerBin=(int)fBin;
		 int upperBin =lowerBin+1;
		 if(lowerBin > 32767)lowerBin=32767;
		 if(upperBin > 32767)upperBin=32767;
		 //increment weighted average bin for the first bead; - should do weighted average...
		 if(fBin - lowerBin > 1e-7){
		  counts[upperBin]+=(fBin-lowerBin);
		  counts[lowerBin]+=(1.0-(fBin-lowerBin));
			}
			else{
				counts[lowerBin]++;
			}
		}	
		for(int i=0;i<32768;i++){
		 if(counts[i]){
				maxBin=i;
			}
	 } 	
		delete[] wellcCounts;
		sumCounts();
	}
	
	density(wellRecord<int16_t,int> wells,int slice){
		//converts a slice of the wellRecord into a density file without normalization
		//if slice parameter is set to zero -the record is converted
		fill(counts,counts+32768,0);
		fill(ccounts,ccounts+32768,0);		

		//put all the wells in their own histogram
		if(slice){
	 	for(int i=1;i<=NCOLORS;i++){
		  for (int j=wells.offsets[i][slice];j<wells.offsets[i][slice+1];j++){
		 		counts[wells.values[i][j]]++;
		 	}	
		 }
		}
		else{	 	
		 for(int i=1;i<=NCOLORS;i++){
		  for (int j=0;j<wells.values[i].size();j++){
		 		counts[wells.values[i][j]]++;
		 	}	
		 }
		}
 	for(int i=0;i<32768;i++){
		 if(counts[i]){
				maxBin=i;
			}
	 }	
		sumCounts(); 
	}
	//calculate density from binary expression file or list of expression files 
	density(std::string inputFile,int fileType){
		fill(counts,counts+32768,0);
		fill(ccounts,ccounts+32768,0);
		maxBin=0;	
	 if (fileType ==1){
		 //text matrix file - one value on each line
	  FILE *fp=fopen(inputFile.c_str(),"r");
		 char line[1024];
		 int i=0;
	  while(fgets(line,1024,fp)){
			 sscanf(line,"%f %f",counts+i,ccounts+i);
			 i++;
		 }	
	  fclose(fp);
	 }
	 //binary matrix file
 	else if (fileType ==2){
	  FILE *fp=fopen(inputFile.c_str(),"r");
	 	fread(counts,sizeof(float),32768,fp);
	 	fread(ccounts,sizeof(float),32768,fp);
	 	fclose(fp);
	 	for(int i=0;i<32768;i++)
	 	 if(counts[i])maxBin=i;
	 }
	}

	//read density directly from a density file

	void sumCounts(){
		ccounts[0]=counts[0];
		for(int i=1;i<32768;i++){
			ccounts[i]=ccounts[i-1]+counts[i];
		}	
	}		
	void sumCounts(int maxBin){
		ccounts[0]=counts[0];
		for(int i=1;i<=maxBin;i++){
			ccounts[i]=ccounts[i-1]+counts[i];
		}	
		for(int i=1;i>maxBin;i++){
			ccounts[i]=ccounts[maxBin];
		}	
	}
	void sumCounts(int *myCounts){
		for(int i=1;i<32768;i++){
			myCounts[i]=myCounts[i-1]+myCounts[i];
		}	
	}

	//write to a binary reference file (one slice)
	void write_refFile(string outputFile){
		FILE *fp=fopen(outputFile.c_str(),"w");
	 fwrite(counts,1,32768*sizeof(float),fp);
	 fwrite(ccounts,1,32768*sizeof(float),fp);
	 fclose(fp);
	}
	//write to a text file	
	void write_textRefFile(string outputFile){
		FILE *fp=fopen(outputFile.c_str(),"w");
		for(int i=0;i<32768;i++)
		 fprintf(fp,"%f %f\n",counts[i],ccounts[i]);
	 fclose(fp);
	}
	int qNorm(wellRecord <int16_t,int> &wells){
		//no sort constant time quantile normalization
		//a sorted method might be faster for very small samples and reference distributions	 
		for(int n=0;n<wells.nTotalWells;n++){
			int rawCounts[32768];
		 int qmap[32768];
		 int rawMaxBin=0;
		 int totalRawCounts=0;
		 fill(qmap,qmap+32768,0);
		 fill(rawCounts,rawCounts+32768,0);
   for(int id=1;id<NCOLORS1OFF;id++){
		  for(int k=wells.offsets[id][n];k<wells.offsets[id][n+1];k++){
					rawCounts[wells.values[id][k]]++;
					totalRawCounts++;
				}
			}
			for (int i=0;i<32768;i++){
				if(rawCounts[i]){
					rawMaxBin=i;
				}	
			}	
	  //calculate the mapping from the raw bins to the normalized bins (qmap)
	  int rawc=0;//cumulative raw counts
	  int refBin=0;//the current reference bin
	  if(!totalRawCounts || !ccounts[32767])continue;
	  double ratio=(double)ccounts[32767]/(double)totalRawCounts;//ratio * rawcounts tells us where to find the next corresponding reference bin
	  for (int rawBin=0;rawBin<=rawMaxBin;rawBin++){
	  	if(rawCounts[rawBin]){
	  		const double limit=((double)rawc+0.5*(double)rawCounts[rawBin])*ratio; //find number of corresponding cumuluative counts in reference distribution
	  		while(ccounts[refBin]<limit && refBin <32767){ //find corresponding bin
	     refBin++;
	  		}
	  		//fprintf(stdout,"%f %d %d %d %f %d\n",ratio,rawBin,refBin,rawCounts[rawBin],limit,ccounts[refBin]);
	 	 	qmap[rawBin]=refBin;
	 		 rawc+=rawCounts[rawBin];
	 	 }
	  }
	 //apply the map to the counts
	 //log2 transform the counts later
		 for(int id=1;id<NCOLORS1OFF;id++){  
			 for(int k=wells.offsets[id][n];k<wells.offsets[id][n+1];k++){
     wells.values[id][k]=qmap[wells.values[id][k]];
				}
			}
		}
	}
};	
