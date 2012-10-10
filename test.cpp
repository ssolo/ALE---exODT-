
#include "exODT.h"
#include "ALE_util.h"
#include "ML_util.h"

using namespace std;
using namespace bpp;

void dtl_estimate(  exODT_model* model,scalar_type N_ales_norm)
{
  scalar_type c=0,delta_avg=0,tau_avg=0,lambda_avg=0;
  scalar_type call=model->last_branch;
  for (int branch=0;branch<model->last_branch;branch++)	
    {
  scalar_type w=0.1;
  scalar_type t_branch=model->t_end[branch];
  scalar_type all_Os=0;
  scalar_type below_Os=0;

  for (int o_branch=0;o_branch<model->last_branch;o_branch++)	
    {
      scalar_type t_o_branch=model->t_begin[o_branch];
      scalar_type Os=model->branch_counts["Os"][o_branch];
      if (t_o_branch>t_branch)
	below_Os+=Os;
      all_Os+=Os;      
    }
  scalar_type fT=1;//below_Os/all_Os;
  
  //cout << "fT " << branch << " "  << fT << endl;
  scalar_type delta=0;
  scalar_type lambda=0;
  scalar_type tau=0;
  scalar_type t=model->t_begin[branch]-model->t_end[branch];
  //obs. prob. being lost
  scalar_type p10=model->branch_counts["Ls"][branch]/model->branch_counts["count"][branch];
  //obs. avg. copy
  scalar_type k=(model->branch_counts["copies"][branch]-model->branch_counts["Ts"][branch])/(model->branch_counts["count"][branch]);
  if (k==1 or k==0 or p10==1  or model->branch_counts["count"][branch]==0)
    {
      delta=1e-6;
    }
  else 
    {
      delta=(p10 + k - 1)* log(k)/(1 - p10)/(k - 1)/t;
    }
  
  if (k==0 or p10==1)
    {
      lambda=model->branch_counts["Ls"][branch]/t;
    }
  else if (k==1 or p10==0 or model->branch_counts["Ls"][branch]==0)
    {
      lambda=1e-6;
    }
  else
    {
      lambda=p10*k*log(k)/(1 - p10)/(k - 1)/t;
    }
  if (fT>0)
    tau=max(  (scalar_type)( model->branch_counts["Ts"][branch]/(N_ales_norm*fT) * lambda) / (1-exp(-lambda*t))  ,(scalar_type)1e-6);	
  else
    tau=1e-6;
  
  if (delta>1e-6 or lambda>1e-6 or tau>1e-6 or true )
    {
      cout << branch <<" " <<delta << " " << tau << " " << lambda << " " << model->branch_counts["copies"][branch] << " " << model->branch_counts["count"][branch] << " " << model->branch_counts["Ds"][branch] << " " << model->branch_counts["Ts"][branch] << " " << model->branch_counts["Ls"][branch]<<endl;
      delta_avg+=delta;
      tau_avg+=tau;
      lambda_avg+=lambda;
      c++;
    }
    }
  //cout <<"cavg "<< delta_avg/c << " " <<  tau_avg/c << " " << lambda_avg/c << endl;
  //cout <<"callavg "<< delta_avg/call << " " <<  tau_avg/call << " " << lambda_avg/call  << " " << model->MLRec_events["D"] << " " << model->MLRec_events["T"] << " " << model->MLRec_events["L"]<< " " << model->MLRec_events["S"]<< endl;

}




int main(int argc, char ** argv)
{
  string Sstring;
  ifstream file_stream_S (argv[1]);
  getline (file_stream_S,Sstring);

  exODT_model* model=new exODT_model();
  model->set_model_parameter("D",2);
  model->set_model_parameter("DD",20);

  model->set_model_parameter("event_node",0);
  model->construct(Sstring);
  //cout << model->string_parameter["S_with_ranks"]<<endl;
  scalar_type N=atof(argv[6]);
  model->set_model_parameter("N",N);
  model->set_model_parameter("Delta_bar",N*2.);
  model->set_model_parameter("Lambda_bar",N*2.);
  model->set_model_parameter("tau", 0.205257);
  model->set_model_parameter("delta",1.5493e-07);
  model->set_model_parameter("lambda",0.22224);
  model->calculate_EGb();
  model->set_model_parameter("leaf_events",1);

  string ale_file=argv[2];
  
  ifstream file_stream (argv[2]);
  vector<string> ale_names;
  if (file_stream.is_open())  //  ########## read trees ############
    {
      while (! file_stream.eof())
	    {
	      string line;
	      getline (file_stream,line);
	      if (line.find("(")!=line.npos || line.find("ale")!=line.npos)
		ale_names.push_back(line);			     
	    }
    }
  for (int i=0;i<(int)ale_names.size();i++)
    {
      //cout << ale_names[i] << endl;
      approx_posterior * ale=load_ALE_from_file(ale_names[i]);//obsorve_ALE_from_string(ale_names[i]);//load_ALE_from_file(ale_file);
      
      scalar_type delta=atof(argv[3]);
      scalar_type tau=atof(argv[4]);
      scalar_type lambda=atof(argv[5]);
      scalar_type diff=0.1;
      model->set_model_parameter("N",N);
      model->set_model_parameter("delta",delta);
      model->set_model_parameter("tau",tau);
      model->set_model_parameter("lambda",lambda);
      model->set_model_parameter("event_node",0);
      model->calculate_EGb();
      
      //model->set_model_parameter("DD",100);
      
      model->set_model_parameter("leaf_events",1);
      
      pair<string, scalar_type> mpp_res = ale->mpp_tree();
      pair<string, scalar_type> res = model->p_MLRec(ale);  
      
      scalar_type  max_new_delta=delta, max_new_tau=tau, max_new_lambda=lambda,max_new_ll=-1e20;
      string max_rtree;
      scalar_type old_ll;
      scalar_type old_p;//=model->p(ale);

      scalar_type j=0;
      /*
      //cout << old_p << " " << delta << " " << tau << " " << lambda << endl;
      cout << res.first << endl;
      //cout << res.second  << endl;
      //return 1;
      //model->MLRec_events.clear();
      
      */
      //cout << TreeTemplateTools::parenthesisToTree(ale_names[i])->getNumberOfLeaves();
      //cout << " " << model->MLRec_events["genes"]<< endl;
      
      //if (TreeTemplateTools::parenthesisToTree(ale_names[i])->getNumberOfLeaves()!=model->MLRec_events["genes"])
      cout << "#R " << res.first << endl;
      cout << "*";
      for (vector<string> :: iterator it=model->Ttokens.begin();it!=model->Ttokens.end();it++)
	cout << *it << "*" ;
      cout << endl;

      //model->show_counts("copies");
      //model->show_counts("count");
      //model->show_counts("Ls");
      //model->show_counts("Ts");
      //model->show_counts("Os");
      cout << "#E " << " " << model->MLRec_events["D"] << " " << model->MLRec_events["T"] << " " << model->MLRec_events["L"]<< " " << model->MLRec_events["S"] <<endl;
      dtl_estimate(model,1);
      model->MLRec_events.clear();
      for (map<string, vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)//del_loc
	for ( vector<scalar_type>::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
	  (*jt)=0;

    }
  

  // int ii=0;
  // while(ii<atoi(argv[6]))
  //   {      
  //     scalar_type new_delta,new_tau,new_lambda;
      
  //     scalar_type r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  //     scalar_type d;
  //     if (r<1./3.) d=RandomTools::randExponential(0.001)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
  //     else if (r<2./3.) d=RandomTools::randExponential(0.01)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
  //     else d=RandomTools::randExponential(0.1)*2*(0.5-RandomTools::giveRandomNumberBetweenZeroAndEntry(1));
  //     r=RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  //     new_delta=delta;
  //     new_tau=tau;
  //     new_lambda=lambda;

  //     if (r<1./3.) new_delta=delta+d;
  //     else if (r<2./3.)new_tau=tau+d;
  //     else new_lambda=lambda+d;

  //     model->set_model_parameter("delta",min(max(new_delta,(scalar_type)1e-6),(scalar_type)10.));
  //     model->set_model_parameter("tau",min(max(new_tau,(scalar_type)1e-6),(scalar_type)10.));
  //     model->set_model_parameter("lambda",min(max(new_lambda,(scalar_type)1e-6),(scalar_type)10.));
  //     model->calculate_EGb();
  //     scalar_type new_p=(model->p(ale));
  //     //scalar_type new_p=model->p_MLRec(ale).second;
  //     j++;
  //     //cout << "try" << new_p << " " << delta << " " << tau << " " << lambda <<  " " << i/j << endl;
  //     if (new_p>=old_p or new_p/old_p>RandomTools::giveRandomNumberBetweenZeroAndEntry(1))
  // 	{
  // 	  i++;
  // 	  delta=min(max(new_delta,(scalar_type)1e-6),(scalar_type)10.);
  // 	  tau=min(max(new_tau,(scalar_type)1e-6),(scalar_type)10.);
  // 	  lambda=min(max(new_lambda,(scalar_type)1e-6),(scalar_type)10.);	  
  // 	  if ((int)i%10==0) 
  // 	    {
  // 	      pair<string,scalar_type> new_res=model->p_MLRec(ale);
  // 	      boost::trim(new_res.first);
  // 	      cout << i << " " << new_p << " " << delta << " " << tau << " " << lambda <<  " " << " " << model->MLRec_events["D"] << " " << model->MLRec_events["T"] << " " << model->MLRec_events["L"] << " " << new_res.second << " " << new_res.first<< endl;
  // 	      model->MLRec_events.clear();
  // 	    }	   
  // 	  old_p=new_p;
  // 	}

  //   }

  // while (abs(old_ll-max_new_ll)>0.1 or diff>0.001)
  //   {

  //     if (old_ll>max_new_ll)
  // 	{
  // 	  diff/=2;
  // 	}
  //     else
  // 	{
  // 	  delta=max_new_delta;
  // 	  tau=max_new_tau;
  // 	  lambda=max_new_lambda;
  // 	}
  //     old_ll=max_new_ll;
  //     max_new_ll=-1e20;      
  //     for (int dd=-1;dd<=+1;dd+=1)
  // 	for (int dt=-1;dt<=+1;dt+=1)
  // 	  for (int dl=-1;dl<=+1;dl+=1)
  // 	    {
  // 	      model->set_model_parameter("delta",max(delta+dd*diff,(scalar_type)1e-6));
  // 	      model->set_model_parameter("tau",max(tau+dt*diff,(scalar_type)1e-6));
  // 	      model->set_model_parameter("lambda",max(lambda+dl*diff,(scalar_type)1e-6));
  // 	      model->calculate_EG();
  // 	      //model->MLRec_events["L"]
  // 	      pair<string,scalar_type> try_res=model->p_MLRec(ale);
  // 	      //scalar_type new_ll=log( try_res.second);
  // 	      scalar_type new_ll=log(model->p(ale));
  // 	      cout << old_ll-new_ll << " " << max(delta+dd*diff,(scalar_type)1e-6) << " " << max(tau+dt*diff,(scalar_type)1e-6) << " " << max(lambda+dl*diff,(scalar_type)1e-6) << endl;
  // 	      if (new_ll>max_new_ll)
  // 		{
  // 		  max_new_ll=new_ll;
  // 		  max_new_delta=max(delta+dd*diff,(scalar_type)1e-6);
  // 		  max_new_tau=max(tau+dt*diff,(scalar_type)1e-6);
  // 		  max_new_lambda=max(lambda+dl*diff,(scalar_type)1e-6);
  // 		  max_rtree=try_res.first;
  // 		}
  // 	    }
  //     cout <<"! "<< old_ll-max_new_ll << " " << max_new_delta << " " << max_new_tau << " " << max_new_lambda << endl;
  //     cout << max_rtree << endl;
  //   }
  // pair<string,scalar_type> final_res = model->p_MLRec(ale);
  // //cout << res.first;
  // cout <<"p="<< res.second << endl;
  // for (vector<string>::iterator it = model->Ttokens.begin();it!=model->Ttokens.end();it++ )
  //   cout <<"T: " <<(*it) << endl;
  
  // //pair<string,scalar_type> mpp_res=ale->mpp_tree();
  // string mpp_tree=mpp_res.first;
  // string jML_tree=res.first;
  // return 1;

  // vector<string > tokens;
  // boost::split(tokens, ale_file,boost::is_any_of("."),boost::token_compress_on);
  
  // string aln_name="/home/ssolo/sccy36/"+tokens[0]+".aln-gb.fasta";
  // cout << "ALN:" << aln_name << endl;
  // cout << mpp_tree << endl;  
  // cout << -tree_LL(mpp_tree,aln_name) <<endl;
  // cout << jML_tree << endl;  
  // cout << -tree_LL(jML_tree,aln_name) <<endl;
  
  // delete ale;
  // delete model;
  // return 1;
}

