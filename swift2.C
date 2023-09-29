/****************************************************************
 swift2.C
 Copyright (C)2023 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/File.H"
#include "BOOM/Regex.H"
#include "BOOM/Exceptions.H"
#include "BOOM/GSL/BetaDistribution.H"
#include "BOOM/GSL/Random.H"
#include "BOOM/Array1DSorter.H"
#include "BOOM/File.H"
#include "SwiftSample.H"
#include "Swift.H"
#include "Experiment.H"
#include "Simulator.H"
#include "EmpiricalPvalues.H"
using namespace std;
using namespace BOOM;

const int PSEUDOCOUNT=1;

class Application {
  Experiment data;
  void readReps(const Vector<String> &fields,int start,Replicates &);
  void skipLines(int num,File &);
  float getMedian(const Array1D<SwiftSample> &samples);
  void getCI(const Array1D<SwiftSample> &samples,float percent,float &left,
	     float &right);
  void getP_reg(const Array1D<SwiftSample> &samples,float lambda,
		float &leftP,float &rightP,float &P);
  void addPseudocounts(Replicates &);
  void save(const Array1D<SwiftSample> &samples,File &) const;
public:
  Application();
  int main(int argc,char *argv[]);
  bool loadInputs(File &,String &variantID);
  void computePValues(Swift &,const float lambda,
		      Array1D<SwiftSample> &realSamples,
		      const Experiment &realData,const int numNulls,
		      float &medianP,float &areaP);
  void reportSummary(Array1D<SwiftSample> &samples,const String &id,
		     float lambda,const float medianP,const float areaP);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
{
  GSL::Random::randomize();
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"s:c:p:");
  if(cmd.numArgs()!=4)
    throw String("\n\
swift2 <input.txt> <lambda> <first-last> <#samples>\n\
   variant indices are 0-based and inclusive\n\
   1.25 is recommended for lambda (min effect size)\n\
   1000 is recommended for #samples\n\
   optional: -s mcmc.txt = save samples in mcmc.txt\n\
             -c conc = use shrinkage with the given concentration\n\
                       (100 is recommended for concentration, min is 2)\n\
             -p = compute empirical p-values\n\
");
  const String infile=cmd.arg(0);
  const float lambda=cmd.arg(1).asFloat();
  const String variantRange=cmd.arg(2);
  const int numSamples=cmd.arg(3).asInt();
  if(lambda<1.0) throw "lambda must be >= 1";
  String sampleFilename=cmd.optParm('s');
  const float concentration=cmd.option('c') ? cmd.optParm('c').asFloat() : 0;
  const bool wantPValues=cmd.option('p');
  const int numNulls=wantPValues ? cmd.optParm('p').asInt() : 0;

  // Get ready to run on input file
  Regex reg("(\\d+)-(\\d+)");
  if(!reg.match(variantRange)) throw "can't parse variant index range";
  const int firstVariant=reg[1].asInt();
  const int lastVariant=reg[2].asInt();
  File f(infile);
  skipLines(firstVariant,f);
  String id;
  File *sampleFile=sampleFilename.empty() ? NULL : 
    new File(sampleFilename,"w");

  Swift swift(concentration);
  Array1D<SwiftSample> samples;
  for(int i=firstVariant ; i<=lastVariant ; ++i) {
    data.DNA.clear(); data.RNA.clear();

    // Load inputs
    if(!loadInputs(f,id)) break;

    // Draw samples
    swift.run(data.DNA,data.RNA,numSamples,samples);

    // Compute empirical p-values if requested
    float medianP, areaP;
    if(wantPValues)
      computePValues(swift,lambda,samples,data,numNulls,medianP,areaP);
    
    // Report median and 95% CI
    reportSummary(samples,id,lambda,medianP,areaP);

    // Save samples if requested
    if(sampleFile) save(samples,*sampleFile);
  }

  delete sampleFile;
  return 0;
}



void Application::skipLines(int num,File &file)
{
  for(int i=0 ; i<num ; ++i) file.getline();
}



void Application::readReps(const Vector<String> &fields,int countField,
			   Replicates &reps)
{
  const int numReps=fields[countField].asInt();
  if(numReps<1) throw "Invalid number of replicates";
  const int begin=countField+1;
  const int end=begin+numReps*2;
  for(int i=begin ; i<end ; i+=2) {
    const int ref=fields[i].asInt(), alt=fields[i+1].asInt();
    reps.add(Replicate(ref,alt));
  }
}



bool Application::loadInputs(File &f,String &variantID)
{
  if(f.eof()) throw "End of file";
  String line=f.getline();
  line.trimWhitespace();
  if(line.isEmpty()) return false;
  Vector<String> fields;
  line.getFields(fields);
  if(fields.size()<7) throw line+" : Not enough fields";
  variantID=fields[0];
  readReps(fields,1,data.DNA);
  const int numDnaReps=fields[1].asInt();
  readReps(fields,2+numDnaReps*2,data.RNA);
  addPseudocounts(data.DNA);
  addPseudocounts(data.RNA);
  return true;
}



float Application::getMedian(const Array1D<SwiftSample> &samples)
{
  // PRECONDITION: samples have been sorted by theta

  int n=samples.size();
  if(n<2) throw "Too few samples to identify median";
  int mid=n/2;
  float median;
  if(n%2==0)
    median=(samples[mid-1].getTheta()+samples[mid].getTheta())/2.0;
  else
    median=samples[mid].getTheta();
  return median;
}



void Application::addPseudocounts(Replicates &reps)
{
  const int n=reps.size();
  for(int i=0 ; i<n ; ++i)
    reps[i].addPseudocount(PSEUDOCOUNT);
}



void Application::getCI(const Array1D<SwiftSample> &samples,
			float percent,float &left,float &right)
{
  // PRECONDITION: samples have been sorted by theta

  float halfAlpha=(1.0-percent)/2.0;
  const int n=samples.size();
  const int countIn=int(n*halfAlpha+5.0/9.0);
  left=samples[countIn].getTheta();
  right=samples[n-countIn].getTheta();
}



void Application::getP_reg(const Array1D<SwiftSample> &samples,
			   float lambda,float &leftP,float &rightP,float &P)
{
  const int n=samples.size();
  float invLambda=1.0/lambda;
  int numLess=0;
  int numGreater=0;
  for(int i=0 ; i<n ; ++i) {
    if(samples[i].getTheta()<invLambda) ++numLess;
    if(samples[i].getTheta()>lambda) ++numGreater;
    leftP=float(numLess)/float(n);
    rightP=float(numGreater)/float(n);
    P=leftP>rightP ? leftP : rightP;
  }
}



void Application::reportSummary(Array1D<SwiftSample> &samples,
				const String &id,float lambda,
				const float medianP,const float areaP)
{
  // Sort the samples
  SwiftSampleComparator cmp;
  Array1DSorter<SwiftSample> sorter(samples,cmp);
  sorter.sortAscendInPlace();
  
  // Get the median
  const float median=getMedian(samples);

  // Get the 95% CI
  float left, right;
  getCI(samples,0.95,left,right);

  // Get p-value-like statistics
  float leftP, rightP, P_reg;
  getP_reg(samples,lambda,leftP,rightP,P_reg);

  // Generate output
  cout<<id<<"\t"<<median<<"\t"<<left<<"\t"<<right<<"\t"<<P_reg
      <<"\t"<<medianP<<"\t"<<areaP<<endl;
}



void Application::save(const Array1D<SwiftSample> &samples,File &f) const
{
  const int n=samples.size();
  for(int i=0 ; i<n ; ++i) {
    const SwiftSample &sample=samples[i];
    float theta=int(sample.getTheta()*10000+5.0/9.0)/10000.0;
    f.print(String(theta));
    if(i+1<n) f.print("\t");
  }
  f.print("\n");
}



void Application::computePValues(Swift &swift,const float lambda,
				 Array1D<SwiftSample> &realSamples,
				 const Experiment &realData,
				 const int numNulls,
				 float &medianP,float &areaP)
{
  // Simulate some nulls
  Array1D<Experiment> nulls(numNulls);
  Simulator simulator;
  simulator.sim(realData,numNulls,nulls);

  // Perform posterior inference on nulls using Swift
  const int numSamples=realSamples.size();
  Array1D< Array1D<SwiftSample> > nullSamples(numNulls);
  for(int i=0 ; i<numNulls ; ++i)
    swift.run(nulls[i].DNA,nulls[i].RNA,numSamples,nullSamples[i]);
  
  // Compute empirical p-values
  EmpiricalPvalues emp;
  medianP=emp.median_p(realSamples,nullSamples);
  areaP=emp.area_p(realSamples,lambda,nullSamples);
}



