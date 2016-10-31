#include<fstream>
#include<stdio.h>
#include<string>
#include<FCSIO.hpp>
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
	string inputFile,outputFile,inputDir,list,refList,inqRefTextFile,inqRefBinFile,outDensityFile,outqRefBinFile,outqRefTextFile,outputDir;
 bool outBinary=0,inBinary=0,logTransform=0,convert=0,outqRef=0,qNorm=0,deCon=0,pruneFlag=0,ensemble=0,outputRawData=0;
 int minBeads=MINBEADS,nThreads=1,metaData=0,nGroups=1,nFold=1;
 float minLogValue,maxLogValue,smoothHalfWindowSize,peakGridSize,densityHalfWindowSize;
	try{  
  using namespace TCLAP;
	 CmdLine cmd("lxbdecon flags", ' ', "0.1");
	 ValueArg<string> inputDirArg ("d","dir","input lxb directory",false,"","string");
	 ValueArg<string> outputDirArg ("o","outputDir","output directory",false,"","string");		
	 ValueArg<string> inputFileArg ("i","inFile","input lxb file",false,"","string");	 
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
  ValueArg<int> nBootStrapsArg ("","nBootStraps","number of bootstrap runs ",false,0,"int");
  ValueArg<int> nFoldArg ("","nFold","number of repetitions for calculating median values (1/nFold samples are excluded from each run)",false,1,"int"); 
  ValueArg<float> bootstrapSizeArg ("","bootstrapSize","proportion of set to sample per bootstrap",false,0.8,"float");
	 SwitchArg convertArg ("c","convert","indicates that the input file will be converted to another format - default is lxb to binary",cmd,false);	 
	 SwitchArg outBinaryArg ("","outBinary","indicates that the output file will be in binary - default is text",cmd,false);
	 SwitchArg inBinaryArg ("","inBinary","indicates that the input file is a binary file - default is lxb",cmd,false);
	 SwitchArg logTransformArg ("l","logTransform","indicates that input file is to be logTransformed",cmd,false);	 	 
	 SwitchArg qNormArg ("q","qNorm","indicates that input file or list is to be quantile normalized",cmd,false);
	 SwitchArg ensembleArg ("e","ensemble","indicates that input file or list of files will be treated as one large experiment",cmd,false);	 
	 SwitchArg deConArg ("","deCon","deConvolute spectra",cmd,false);
  SwitchArg outputRawDataArg ("","outputRawData","the raw expression numbers are outputted for each experiment",cmd,false); 
	 cmd.add(inputFileArg);
	 cmd.add(inputDirArg);	 
	 cmd.add(outputDirArg);
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
	 cmd.add(nFoldArg);
	 cmd.parse( argc, argv );

	 outBinary=outBinaryArg.getValue();
	 inBinary=inBinaryArg.getValue();
	 logTransform=logTransformArg.getValue();
	 convert=convertArg.getValue();
  inputFile= inputFileArg.getValue();  
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
  nGroups=nGroupsArg.getValue();
  outputRawData=outputRawDataArg.getValue();
  nFold=nFoldArg.getValue();
	} 
	catch (TCLAP::ArgException &e)  // catch any exceptions
	{ cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }
 cerr << list <<endl;
 
 //conversion routines
 //these routines were used to convert lxb files to level1.5 binary files
 //should move these into a separate toolbox for file conversions
 
 if(convert){
		//convert a list of lxb files to a binary file
 	size_t counts[NCOLORS1OFF]; //0 analyte is present - indicates error?
	 memset(counts,0,(NCOLORS1OFF)*sizeof(size_t));
	 if(inputDir != ""){ 
		 namespace fs = boost::filesystem;
   if ( fs::exists(inputDir) && fs::is_directory(inputDir)){
			 cerr << "working on directory " << inputDir << endl;

    fs::directory_iterator end_itr;
    fs::path p(inputDir); 
    fs::create_directory(inputDir+"_b");

  			 //gather all the filenames - this is necessary for OpenMP parallelization 
    vector <string> files;
    for (fs::directory_iterator itr(p); itr != end_itr; ++itr){
				 if (is_regular_file(itr->path())) {
   		 if(inBinary){
					 	files.push_back(itr->path().stem().string()+itr->path().extension().string());
					 }	      
					 else if(itr->path().extension().string() == ".lxb"){
					 	files.push_back(itr->path().stem().string()+itr->path().extension().string());
					 }
					 else{
				 		continue;
			 		}	 
			 	}
			 }
		  #pragma omp parallel for num_threads(nThreads)  
    for (int f=0;f<files.size();f++){
					//read in each lxb file into values vector
				 vector<float> values[NCOLORS1OFF];
     string file = files[f];
     cerr << "checking file " << inputDir+"/"+file <<" of type binary" <<endl;
     if(file.substr(file.length()-4,4) == ".lxb"){
					 process_lxb_file(inputDir+"/"+file,values,NCOLORS1OFF,logTransform);
				 }   
				 string outputFile=inputDir+"_b/"+file+".bin";
				 cerr << "converting file to outputFile " << outputFile <<endl;
		 	 FILE *outfp=fopen(outputFile.c_str(),"w");
		 	 for (int i=1; i<NCOLORS1OFF; ++i){
		 	  if(values[i].size()){
		 	  	int16_t buffer[2];
		 	  	buffer[0]=i;
		 	  	buffer[1]=(int16_t)values[i].size();
		 	  	fwrite (buffer ,sizeof(int16_t),2,outfp);
	      for(int k=0;k<values[i].size();k++){
					 	 int16_t v =(int16_t) (values[i][k]+.5);
		 	  		fwrite (&(v),sizeof(int16_t),1,outfp);
					  }
		 	  }
		 	 }
		 	 fclose(outfp);
			 }
			}
	 }	
	 else if(inputFile!=""){
   cerr << "checking file " << inputFile <<endl;
   vector<float> values[NCOLORS1OFF];
   if(inputFile.substr(inputFile.length()-4) == ".lxb"){
			 process_lxb_file(inputFile,values,NCOLORS1OFF,logTransform);
		 }	
		 string outputFile=inputFile.substr(0,inputFile.length()-4) +".bin";
	  cerr << "converting file to outputFile " << outputFile <<endl;
		 FILE *outfp=fopen(outputFile.c_str(),"w");
		 for (int i=1; i<NCOLORS1OFF; ++i){
		 	if(values[i].size()){
		 	 int16_t buffer[2];
		 	 buffer[0]=i;
		 	 buffer[1]=(int16_t)values[i].size();
		 	 fwrite (buffer ,sizeof(int16_t),2,outfp);
	    for(int k=0;k<values[i].size();k++){
						int16_t v =(int16_t) values[i][k]+.5;
		 	 	fwrite (&(v),sizeof(int16_t),1,outfp);
					}
		 	}
		 }
		 fclose(outfp);
		}
	}
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
	//if we are not doing a conversion or making baseline files - we are deconvoluting
	
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
		//decon individual wells
		if(outputDir == "" ){
			char buffer [L_tmpnam];
   tmpnam (buffer);
   string bufferStr=buffer;
   outputDir=="output_"+bufferStr;
		}
		//deconvolute individual groups
		//make sure that there are no remainders
		 int nIterations=(wells.groupNames.size()/nGroups)*nGroups;
		 #pragma omp parallel num_threads(nThreads)
			{
				vector <DeconRecord <float>> deconRecords(nFold); 
				vector <MedianRecord<float>> medianRecords(501);   
		  for(int i=0;i<nIterations;i+=nGroups){
			 	int t=0;
			 	if(nThreads >1){
			 		t=omp_get_thread_num();
				 }
				 //initialize the deconRecord to nans and zeros 
				 for(int f=0;f<nFold;f++){
						deconRecords[f].blank();
					}	
					for(int k=0;k<501;k++){
					 medianRecords[k].blank();
					} 
			  string groupName=wells.groupNames[i];
				 FILE *fpRaw=0;	
				 if(outputRawData){				
				  string rawOutputFile=outputDir+"/"+groupName+".raw";
				  fpRaw=fopen(rawOutputFile.c_str(),"w");
				 }
			  string peaksFile=outputDir+"/"+groupName+".peaks";
			  string medianFile=outputDir+"/"+groupName+".medians";
			  if(wells.deconvolute(nFold,i,nGroups,deconRecords,medianRecords,fpRaw)){
				 	FILE *fpPeaks=fopen(peaksFile.c_str(),"w");
				 	for(int f=0;f<nFold;f++){
				 	 deconRecords[f].print_text(fpPeaks,metaData);
						}
						fclose(fpPeaks);
						FILE *fpMedian=fopen(medianFile.c_str(),"w");
						for	(int k=0;k<501;k++){
							medianRecords[k].calculate_medians();
							fprintf(fpMedian,"%d %f ",k,medianRecords[k].median0);
							for(int f=0;f<medianRecords[k].means0.size();f++){
							 fprintf(fpMedian,"%f:",medianRecords[k].means0[f]);
							}
							fprintf(fpMedian," %f ",medianRecords[k].median1);
							for(int f=0;f<medianRecords[k].means0.size();f++){
							 fprintf(fpMedian,"%f:",medianRecords[k].means1[f]);
							}
							fprintf(fpMedian,"\n");
						}		
				 	fclose(fpMedian);
				 }
				 if(fpRaw)	fclose(fpRaw);
				}
			}
		}
	
 //are we dumping out a density file
 if(outDensityFile != ""){
		wells.dumpDensity(outDensityFile,densityHalfWindowSize);
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
void process_lxb_file(std::string inputFile,vector<float> *values,int nBins,bool logTransform){
 using namespace FCSTools;
 std::fstream file (inputFile, std::ios::binary|std::ios::in);
 FCS<std::size_t> fcs= Reader<std::size_t> (file, 1);
 int nRecords=fcs.Data.size ();
 int colRID,colRP1;
 //find the RID and RP1 records 
 for (int i=0; i<fcs.Head.Parameter.size();++i){
		if(fcs.Head.Parameter[i].Name == "RID"){
			colRID=i;
		}
		else if(fcs.Head.Parameter[i].Name == "RP1"){
			colRP1=i;
		}
	}
	if(logTransform){
 	for (int i=0; i<fcs.Data.size(); ++i){
 	 if(fcs.Data[i][colRID] && fcs.Data[i][colRID]<nBins){
 			if(fcs.Data[i][colRP1] <= 2) values[fcs.Data[i][colRID]].push_back(1);
 			else values[fcs.Data[i][colRID]].push_back(log(fcs.Data[i][colRP1])*INVLOG2);
 		}
 	}
	}
	else{	
	 for (int i=0; i<fcs.Data.size(); ++i){
	  if(fcs.Data[i][colRID] && fcs.Data[i][colRID]<nBins){
	   values[fcs.Data[i][colRID]].push_back(fcs.Data[i][colRP1]);
		 }
		}
	}
}

