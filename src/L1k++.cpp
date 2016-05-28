#include<fstream>
#include<stdio.h>
#include<string>
#include<tclap/CmdLine.h>
#include<common.hpp>
#include<gmm.hpp>
#include<qNorm.hpp>
#include<omp.h>

void show_help ();
void process_lxb_file(std::string inputFile,vector<float> *values,int nBins,bool logTransform);
void process_binary_file(std::string inputFile,vector<float> *values,int nBins,bool logTransform);

int main (int argc, char *argv[]){
	using namespace std;
	string inputFile,outputFile,inputDir,list,refList,inqRefTextFile,inqRefBinFile,outDensityFile,outqRefBinFile,outqRefTextFile,outputDir,scoreFile;
 bool outBinary=0,inBinary=0,logTransform=0,outqRef=0,qNorm=0,deCon=0,pruneFlag=0,ensemble=0,clusterSeed=0,outputRawData=0;
 int minBeads=MINBEADS,nThreads=1,metaData=0,nGroups=1;
 float minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize,densityHalfWindowSize;
	try{  
  using namespace TCLAP;
	 CmdLine cmd("lxbdecon flags", ' ', "0.1");
	 ValueArg<string>	scoreFileArg ("s","scoreFile","score the wells using this a given density file",false,"","string"); 
	 ValueArg<string> inputDirArg ("d","dir","input lxb directory",false,"","string");
	 ValueArg<string> outputDirArg ("","outputDir","output directory",false,"","string");		
	 ValueArg<string> inputFileArg ("i","inFile","input lxb file",false,"","string");	 
	 ValueArg<string> outputFileArg ("o","outputFile","output file",false,"","string");
	 ValueArg<string> outDensityFileArg ("","outDensityFile","output density file",false,"","string");
	 ValueArg<string> outqRefBinFileArg ("","outqRefBinFile","make binary quantile reference density file from list of binary level 1.5 data",false,"","string");
	 ValueArg<string> outqRefTextFileArg ("","outqRefTextFile","make text quantile reference density file from list of binary level 1.5 data",false,"","string");
	 ValueArg<string> listArg ("","list","list of files to be processed",false,"","string");	 	 
	 ValueArg<string> refListArg ("","refList","list of files to be used as a reference for qNormalization",false,"","string");	 
	 ValueArg<string> inqRefTextFileArg ("","inqRefText","defines reference distribution for quantile normalization file in text format",false,"","string");
	 ValueArg<string> inqRefBinFileArg ("","inqRefBin","defines reference distribution for quantile normalization filmat",false,"","string");	 
	 ValueArg<int> metaDataArg ("m","metaData","indicates the amount of metadata to be output with deconvolution ",false,0,"int"); 
	 ValueArg<int> nThreadsArg ("n","nThreads","number of threads ",false,1,"int");	 
	 ValueArg<float> minLogValueArg ("","minLogValue","minimum expression value in log2 scale  ",false,3.0,"float");
	 ValueArg<float> maxLogValueArg ("","maxLogValue","maximum expression value in log2 scale ",false,14.5,"float");		 
	 ValueArg<float> smoothHalfWindowSizeArg ("","smoothHalfWindowSize","windowing size for smoothing data ",false,0.2,"float");
	 ValueArg<float> densityHalfWindowSizeArg ("","densityHalfWindowSize","windowing size for smoothing density ",false,0.1,"float");
	 ValueArg<float> peakGridSizeArg ("","peakGridSize","maximum expression value in log2 scale ",false,0.2,"float");
	 ValueArg<int>nGroupsArg("","nGroups","number of groups to group together - for testing purposes",false,1,"int");	

	 SwitchArg clusterSeedArg ("","clusterSeed","fit using cluster centers as start point",cmd,false);
	 SwitchArg outBinaryArg ("","outBinary","indicates that the output file will be in binary - default is text",cmd,false);
	 SwitchArg inBinaryArg ("","inBinary","indicates that the input file is a binary file - default is lxb",cmd,false);
	 SwitchArg logTransformArg ("l","logTransform","indicates that input file is to be logTransformed",cmd,false);	 	 
	 SwitchArg qNormArg ("q","qNorm","indicates that input file or list is to be quantile normalized",cmd,false);
	 SwitchArg ensembleArg ("e","ensemble","indicates that input file or list of files will be treated as one large experiment",cmd,false);	 
	 SwitchArg deConArg ("","deCon","deConvolute spectra",cmd,false);
  SwitchArg outputRawDataArg ("","outputRawData","the raw expression numbers are outputted for each experiment",cmd,false); 
  cmd.add(scoreFileArg);
	 cmd.add(inputFileArg);
	 cmd.add(inputDirArg);	 
	 cmd.add(outputDirArg);
	 cmd.add(outputFileArg);
	 cmd.add(metaDataArg);
	 cmd.add(nThreadsArg);	 	 
	 cmd.add(nGroupsArg);	 
	 cmd.add(listArg);
	 cmd.add(refListArg);
	 cmd.add(inqRefTextFileArg);
	 cmd.add(inqRefBinFileArg);
	 cmd.add(outqRefTextFileArg);
	 cmd.add(outqRefBinFileArg);	 
	 cmd.add(outDensityFileArg);	 	 
	 cmd.add(minLogValueArg);
	 cmd.add(maxLogValueArg);
	 cmd.add(smoothHalfWindowSizeArg);
	 cmd.add(densityHalfWindowSizeArg);
	 cmd.parse( argc, argv );

	 outBinary=outBinaryArg.getValue();
	 inBinary=inBinaryArg.getValue();
	 logTransform=logTransformArg.getValue();
  inputFile= inputFileArg.getValue();  
  outputFile= outputFileArg.getValue();
  inputDir=inputDirArg.getValue();  
  outputDir=outputDirArg.getValue();
  metaData=metaDataArg.getValue();
  nThreads=nThreadsArg.getValue();
  outDensityFile=outDensityFileArg.getValue();
  inqRefTextFile=inqRefTextFileArg.getValue();
  inqRefBinFile=inqRefBinFileArg.getValue();
  outqRefTextFile=outqRefTextFileArg.getValue();
  outqRefBinFile=outqRefBinFileArg.getValue();
  list=listArg.getValue();
  refList=refListArg.getValue();
  qNorm=qNormArg.getValue();
  deCon=deConArg.getValue();
  ensemble=ensembleArg.getValue();
  minLogValue=minLogValueArg.getValue();
  maxLogValue=maxLogValueArg.getValue();  
  peakGridSize=peakGridSizeArg.getValue();
  smoothHalfWindowSize=smoothHalfWindowSizeArg.getValue();
  densityHalfWindowSize=densityHalfWindowSizeArg.getValue();
  clusterSeed=clusterSeedArg.getValue();
  nGroups=nGroupsArg.getValue();
  outputRawData=outputRawDataArg.getValue();
	} 
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
 cerr << list <<endl;
 
	//routines to make qNormed baseline files 
 if(outqRefBinFile != "" && refList == "" && list != ""){ //make qNorm files from list of files and save
  fprintf (stderr,"refList not given - using list %s to generate quantile normalization reference file\n",list.c_str()); 
  wellRecord <int16_t,int> wells(list,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
		density ref(wells);
		ref.write_refFile(outqRefBinFile.c_str());
		return(1);
	}
	else if(outqRefBinFile != "" && refList != ""){ //make qNorm files from list of files and save
  wellRecord <int16_t,int> wells(refList,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
		density ref(wells);
		ref.write_refFile(outqRefBinFile.c_str());
		return(1);
	}
	else if(outqRefTextFile !="" && refList != ""){
  wellRecord <int16_t,int> wells(list,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
		density ref(wells);
		ref.write_textRefFile(outqRefTextFile.c_str());
	 return(1);
	}
	else if(outqRefTextFile !="" && refList == "" && list != ""){
		 fprintf (stderr,"refList not given - using list %s to generate quantile normalization reference file\n",list.c_str());
  wellRecord <int16_t,int> wells(list,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
		density ref(wells);
		ref.write_textRefFile(outqRefTextFile.c_str());
	 return(1);
	}
		//read the list of files in binary format
	wellRecord <int16_t,int> wells(list,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
 for(int i=0;i<wells.groupNames.size();i++){
		fprintf(stderr,"%d %d name %s offset %d\n",i,wells.nTotalWells,wells.groupNames[i].c_str(),wells.groupOffsets[i]);
	}	 

	density ref;
 //are we doing quantile normalization? 	
 		
 if(qNorm){
		//where are we getting the reference distribution from?  
		if(inqRefTextFile !=""){
			ref=density(inqRefTextFile,1);
		}
		else if	(inqRefBinFile != ""){
   ref=density(inqRefBinFile,2);
		}
		else if(refList !=""){  
			wellRecord <int16_t,int> qwells(refList,minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize);
		 ref=density(qwells);
		}
		else{
			ref=density(wells);
		}
		ref.qNorm(wells);
	}
	//do we deconvolute?
	if(deCon){
		if(outBinary){
			//check if we are outputing binary files
			metaData=-1;
		}	
		//do we treat this as an ensemble or as individual files 
  if(ensemble){
			DeconRecord <float> deconRecord; 
   if(outputFile != ""){
				string rawOutputFile=outputFile+".raw";
				string smoothOutputFile=outputFile+".smooth";
				string peaksFile=outputFile+".peaks";
				FILE *fpRaw=fopen(rawOutputFile.c_str(),"w");
				FILE *fpSmooth=fopen(smoothOutputFile.c_str(),"w");
				wells.deconvolute(deconRecord,fpRaw,fpSmooth);
		  fclose(fpRaw);
		  fclose(fpSmooth);
				FILE *fpPeaks=fopen(peaksFile.c_str(),"w");
				deconRecord.print_text(fpPeaks,metaData);
				fclose(fpPeaks);
			}   
			else{	
				wells.deconvolute(deconRecord,0,0);
				deconRecord.print_text(stdout,metaData);
			}	
	 }
	 else{
			//decon individual wells
			vector <DeconRecord <float>> deconRecords(nThreads);
			FILE *fpAllPeaks=0;
			if(outputFile != ""){					
				fpAllPeaks=fopen(outputFile.c_str(),"w");
			}				
			//deconvolute individual groups
			//make sure that there are no remainders - i.e we are really using 
			int nIterations=(wells.groupNames.size()/nGroups)*nGroups;
			if(outputDir != ""){
				 #pragma omp parallel for num_threads(nThreads)  
			 for(int i=0;i<nIterations;i+=nGroups){
					int t=0;
					if(nThreads >1){
						t=omp_get_thread_num();
					}	
				 string groupName=wells.groupNames[i];
					FILE *fpRaw=0;	
					if(outputRawData){				
					 string rawOutputFile=outputDir+"/"+groupName+".raw";
					 fpRaw=fopen(rawOutputFile.c_str(),"w");
					}
				 string peaksFile=outputDir+"/"+groupName+".peaks";
				 if(wells.deconvolute(i,nGroups,deconRecords[t],fpRaw)){
						FILE *fpPeaks=fopen(peaksFile.c_str(),"w");
						deconRecords[t].print_text(fpPeaks,metaData);		
						fclose(fpPeaks);
					}
					if(fpRaw)	fclose(fpRaw);
				}
			}
			else{
				for(int i=0;i<nIterations;i+=nGroups){
					int t=0;
				 if(wells.deconvolute(i,nGroups,deconRecords[t],0)){
						string groupName=wells.groupNames[i];
						if(fpAllPeaks){
						 fprintf(fpAllPeaks,"%s\n",groupName.c_str());
						 deconRecords[t].print_text(fpAllPeaks,metaData);
						}
						else{
						 fprintf(stdout,"%s\n",groupName.c_str());
						 deconRecords[t].print_text(stdout,metaData);
						}		
					}		
				}
			}		
			if(fpAllPeaks)fclose(fpAllPeaks);
		}
	}
 //are we dumping out a density file
 if(outDensityFile != ""){
		wells.dumpDensity(outDensityFile,smoothHalfWindowSize);
	}	
	return(1);
}
void process_binary_file(std::string inputFile,vector<float> *values,int nBins,bool logTransform){
 int16_t header[2];
 FILE *fp=fopen(inputFile.c_str(),"r");
 while(fread(&header,2*sizeof(int16_t),1,fp)){
		values[header[0]].resize(0);
		values[header[0]].reserve(header[1]);
		for(int k=0;k<values[header[0]].size();k++){
			int16_t v;
		 int nRead=fread(&v,sizeof(int16_t),1,fp);
		 if(nRead){
		  if(!logTransform)values[header[0]][k]=(float)v;
		  else if(v<=1)values[header[0]][k]=0;
		  else values[header[0]][k]=log((float)v)*INVLOG2;
			}
		}
	}
	fclose(fp);
}
