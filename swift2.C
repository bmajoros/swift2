/****************************************************************
 swift.C
 Copyright (C)2018 William H. Majoros (bmajoros@alumni.duke.edu)
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
#include "BOOM/VectorSorter.H"
#include "BOOM/File.H"
#include "Replicates.H"
#include "SwiftSample.H"
using namespace std;
using namespace BOOM;

const int PSEUDOCOUNT=1;

class Application {
  Replicates DNA, RNA;
  Vector<SwiftSample> samples;
  void readReps(const Vector<String> &fields,int start,Replicates &);
  void skipLines(int num,File &);
  float getMedian();
  void getCI(float percent,float &left,float &right);
  void getP_reg(float lambda,float &leftP,float &rightP,float &P);
  void addPseudocounts(Replicates &);
  void save(File &) const;
public:
  Application();
  int main(int argc,char *argv[]);
  bool loadInputs(File &,String &variantID);
  void performSampling(int numSamples,float concentration);
  void reportSummary(const String &id,float lambda);
  void chooseAlphaBeta(float c,const float p,float &alpha,float &beta);
  void advanceReps(int &dnaIndex,int &rnaIndex,const int numDnaReps,
		   const int numRnaReps);
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
  CommandLine cmd(argc,argv,"s:c:");
  if(cmd.numArgs()!=4)
    throw String("\n\
swift2 <input.txt> <lambda> <first-last> <#samples>\n\
   variant indices are 0-based and inclusive\n\
   1.25 is recommended for lambda (min effect size)\n\
   1000 is recommended for #samples\n\
   optional: -s mcmc.txt = save samples in mcmc.txt\n\
             -c conc = use shrinkage with the given concentration\n\
                       (100 is recommended for concentration, min is 2)\n\
");
  const String infile=cmd.arg(0);
  const float lambda=cmd.arg(1).asFloat();
  const String variantRange=cmd.arg(2);
  const int numSamples=cmd.arg(3).asInt();
  if(lambda<1.0) throw "lambda must be >= 1";
  String sampleFilename=cmd.optParm('s');
  const float concentration=cmd.option('c') ? cmd.optParm('c').asFloat() : 0;

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

  for(int i=firstVariant ; i<=lastVariant ; ++i) {
    DNA.clear(); RNA.clear(); samples.clear();

    // Load inputs
    if(!loadInputs(f,id)) break;

    // Draw samples
    performSampling(numSamples,concentration);
    //cout<<samples.size()<<" samples were drawn"<<endl;
    
    // Report median and 95% CI
    reportSummary(id,lambda);

    // Save samples if requested
    if(sampleFile) save(*sampleFile);
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
  readReps(fields,1,DNA);
  const int numDnaReps=fields[1].asInt();
  readReps(fields,2+numDnaReps*2,RNA);
  //DNA.collapse();
  //RNA.collapse();
  addPseudocounts(DNA);
  addPseudocounts(RNA);
  return true;
}



void Application::performSampling(int numSamples,float conc)
{
  //const int a=DNA[0].getAlt(), b=DNA[0].getRef();
  //const int k=RNA[0].getAlt(), m=RNA[0].getRef();
  const int numDnaReps=DNA.size(), numRnaReps=RNA.size();
  int nextDnaRep=0, nextRnaRep=0;
  for(int i=0 ; i<numSamples ; ++i) {
    // Get counts from current reps
    const int a=DNA[nextDnaRep].getAlt(), b=DNA[nextDnaRep].getRef();
    const int k=RNA[nextRnaRep].getAlt(), m=RNA[nextRnaRep].getRef();
    
    // Sample p from P(p|a,b)
    GSL::BetaDistribution beta1(a+1,b+1); // posterior with uniform prior
    float p;
    do { p=beta1.random(); } while(p==0.0 || p==1.0);
    //cout<<"a="<<a<<" b="<<b<<" p="<<p;

    // Sample q from P(q|p,k,m)
    float alpha, beta;
    chooseAlphaBeta(conc,p,alpha,beta);
    GSL::BetaDistribution beta2(k+alpha,m+beta); // posterior
    float q;
    do { q=beta2.random(); } while(q==0.0 || q==1.0);
    //cout<<" k="<<k<<" m="<<m<<" q="<<q<<endl;

    // Create a new sample via (q/(1-q))/(p/(1-p))
    samples.push_back(SwiftSample(p,q));

    // Advance to next pairing of DNA & RNA replicates
    advanceReps(nextDnaRep,nextRnaRep,numDnaReps,numRnaReps);
  }
}



void Application::chooseAlphaBeta(float c,const float p,
				  float &alpha,float &beta)
{
  // If c>2, this will use a beta prior with p as the mode.
  // For c=2, this will result in a uniform prior.
  if(c<2) c=2;

  // ### ADAPTIVE CONCENTRATION METHOD:
  //const float c=5.5/p;
  // ###

  // Compute parameters for a beta prior with mode and concentration
  alpha=p*(c-2)+1;
  beta=(1-p)*(c-2)+1;
}



void Application::advanceReps(int &dnaIndex,int &rnaIndex,const int numDnaReps,
			      const int numRnaReps)
{
  ++dnaIndex;
  if(dnaIndex>=numDnaReps) {
    dnaIndex=0;
    rnaIndex=(rnaIndex+1)%numRnaReps;
  }
}



float Application::getMedian()
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



void Application::getCI(float percent,float &left,float &right)
{
  // PRECONDITION: samples have been sorted by theta

  float halfAlpha=(1.0-percent)/2.0;
  const int n=samples.size();
  const int countIn=int(n*halfAlpha+5.0/9.0);
  left=samples[countIn].getTheta();
  right=samples[n-countIn].getTheta();
}



void Application::getP_reg(float lambda,float &leftP,float &rightP,float &P)
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



void Application::reportSummary(const String &id,float lambda)
{
  // Sort the samples
  SwiftSampleComparator cmp;
  VectorSorter<SwiftSample> sorter(samples,cmp);
  sorter.sortAscendInPlace();
  
  // Get the median
  const float median=getMedian();

  // Get the 95% CI
  float left, right;
  getCI(0.95,left,right);

  // Get p-value-like statistics
  float leftP, rightP, P_reg;
  getP_reg(lambda,leftP,rightP,P_reg);

  // Generate output
  cout<<id<<"\t"<<median<<"\t"<<left<<"\t"<<right<<"\t"<<P_reg<<endl;
}



void Application::save(File &f) const
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



