#include "ML_util.h"

using namespace bpp;
using namespace std;

scalar_type tree_LL(string tree,string aln_filename,bool optimize_bls,scalar_type tolerance)
{

	const Alphabet* alphabet = new ProteicAlphabet();
	OrderedSequenceContainer *alignment;
	VectorSiteContainer* sites;
	Fasta Reader;
	//Phylip * Reader=new Phylip(true,true,100,true,"\r");
	alignment = Reader.read(aln_filename, alphabet);
	sites = new VectorSiteContainer(*alignment);
	SiteContainerTools::changeGapsToUnknownCharacters(*sites);
	
	TreeTemplate<Node>* ttree1=TreeTemplateTools::parenthesisToTree(tree,false,"ID");

	//Newick newick1;
	//ttree1 = newick1.read(tree);

	DiscreteRatesAcrossSitesTreeLikelihood* tl1;
	SubstitutionModel*    model    = 0;
	DiscreteDistribution* rDist    = 0;	

	model = new LG08(&AlphabetTools::PROTEIN_ALPHABET, new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET), true);
	model->setFreqFromData(*sites);

	rDist = new GammaDiscreteDistribution(4, 1, 1);

	tl1 = new RHomogeneousTreeLikelihood(*ttree1, *sites, model, rDist, true, false, false);
	tl1->initialize();
		/*

	if (optimize_bls)
	  {
	    Optimizer* optimizer = new PseudoNewtonOptimizer(tl1);
	    //	  Optimizer* optimizer = new PseudoNewtonOptimizer(tl1);

	    ParameterList * parameters= new ParameterList();
	    parameters->addParameters( tl1->getBranchLengthsParameters());
	    parameters->addParameters( tl1->getRateDistributionParameters());
	    //Newton..
	    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
	    optimizer->setProfiler(0);
	    optimizer->setMessageHandler(0);
	    optimizer->setVerbose(0);
	    optimizer->getStopCondition()->setTolerance(0.01);
	    optimizer->init(*parameters);
	    //optimizer->init(tl1->getParameters());
	    optimizer->setMaximumNumberOfEvaluations(1000);
	    optimizer->optimize();
	    delete  parameters;
	    delete optimizer;       
	
	  }
		*/
	if (optimize_bls)
	  {
	    //Newton..
	    ParameterList * parameters= new ParameterList();
	    parameters->addParameters( tl1->getBranchLengthsParameters());
	    parameters->addParameters( tl1->getRateDistributionParameters());
	    unsigned int k = OptimizationTools::optimizeNumericalParameters(
									     dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>  (tl1),
									     //tl1->getParameters(),
									     *parameters,
									     0,
									     1,
									     tolerance,
									     1000,
									     0,
									     0,
									     false,
									     0,
									     OptimizationTools::OPTIMIZATION_NEWTON,
									     //OptimizationTools::OPTIMIZATION_BRENT);
									     OptimizationTools::OPTIMIZATION_BFGS);
	
	    delete parameters;
	      }
	scalar_type LL=- tl1->getValue(); //Here's your log likelihood value !

	delete sites;
	delete alphabet;
	delete model;
	delete rDist;
	delete tl1;
	return 	LL;
}
	
	

vector<string> all_NNIs(string Sstring,bool rooted)
{
  tree_type * S = TreeTemplateTools::parenthesisToTree(Sstring);
  vector<string> NNIs;
  vector < Node * > nodes = S->getNodes();
  NNIs.push_back(TreeTemplateTools::treeToParenthesis(*S));
  if (nodes.size()<4)
    return NNIs;
  Node * root=S->getRootNode();
  for (vector < Node *> :: iterator ni=nodes.begin();ni!=nodes.end(); ni++)
    if (!((*ni)==root) && !((*ni)->isLeaf()) && (!((*ni)->getFather()==root) || rooted))
      {
	Node * node0=(*ni)->getSon(0);
	Node * node1=(*ni)->getSon(1);
	Node * father=(*ni)->getFather();
	Node * node2;
	if (father->getSon(0)==(*ni))
	  node2=father->getSon(1);
	else
	  node2=father->getSon(0);
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	//NNI 1.
	(*ni)->addSon(node0);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node1);
       
	NNIs.push_back(TreeTemplateTools::treeToParenthesis(*S));
	//disolve
	(*ni)->removeSons();
	father->removeSons();
	
	//NNI 2.
	(*ni)->addSon(node1);
	(*ni)->addSon(node2);
	father->addSon((*ni));
	father->addSon(node0);

	NNIs.push_back(TreeTemplateTools::treeToParenthesis(*S));
	//disolve
	(*ni)->removeSons();
	father->removeSons();

	//restore
	(*ni)->addSon(node0);
	(*ni)->addSon(node1);
	father->addSon((*ni));
	father->addSon(node2);
      }
  
  return NNIs; 
}




pair<string,scalar_type> tree_LL_nucl(string tree,string aln_filename,bool optimize_bls,scalar_type tolerance)
{
  //const Alphabet* alphabet = new ProteicAlphabet();
  const Alphabet* alphabet = new RNA();
	OrderedSequenceContainer *alignment;
	VectorSiteContainer* sites;
	Fasta Reader;
	//NexusIOSequence Reader;
	//Phylip * Reader=new Phylip(true,true,100,true,"\r");
	alignment = Reader.read(aln_filename, alphabet);
	sites = new VectorSiteContainer(*alignment);
	SiteContainerTools::removeGapOnlySites(*sites);	
	SiteContainerTools::changeGapsToUnknownCharacters(*sites);	

	TreeTemplate<Node>* ttree1=TreeTemplateTools::parenthesisToTree(tree,false,"ID");
	DiscreteRatesAcrossSitesTreeLikelihood* tl1;
	SubstitutionModel*    model    = 0;
	DiscreteDistribution* rDist    = 0;	
	model = new GTR(&AlphabetTools::RNA_ALPHABET);
	model->setFreqFromData(*sites);
	rDist = new GammaDiscreteDistribution(8, 1, 1);
	tl1 = new RHomogeneousTreeLikelihood(*ttree1, *sites, model, rDist, true, false, false);
	tl1->initialize();
	if (optimize_bls)
	  {
	    //Newton..
	    ParameterList * parameters= new ParameterList();
	    parameters->addParameters( tl1->getBranchLengthsParameters());
	    parameters->addParameters( tl1->getRateDistributionParameters());
	    unsigned int k = OptimizationTools::optimizeNumericalParameters(
									     dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>  (tl1),
									     //tl1->getParameters(),
									     *parameters,
									     0,
									     1,
									     tolerance,
									     1000,
									     0,
									     0,
									     false,
									     0,
									     OptimizationTools::OPTIMIZATION_NEWTON,
									     //OptimizationTools::OPTIMIZATION_BRENT);
									     OptimizationTools::OPTIMIZATION_BFGS);
	
	    delete parameters;
	      }
	scalar_type LL=- tl1->getValue(); //Here's your log likelihood value !
	//tl1->getParameters().printParameters(cout);
	//cout << TreeTemplateTools::treeToParenthesis( tl1->getTree() ) <<endl;
	pair<string,scalar_type> return_pair;
	return_pair.first= TreeTemplateTools::treeToParenthesis( tl1->getTree() ) ;
	return_pair.second=LL;
	delete sites;
	delete alphabet;
	delete model;
	delete rDist;
	delete tl1;
	return 	return_pair;
}
