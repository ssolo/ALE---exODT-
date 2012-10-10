#include "mpi_tree.h"
#include "ALE_util.h"

using namespace std;
using namespace bpp;
using namespace boost::mpi;


vector<scalar_type> mpi_tree::dtl_estimate(int branch, scalar_type N_ales_norm)
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
  scalar_type fT=below_Os/all_Os;
  
  //cout << "fT " << branch << " "  << fT << endl;
  scalar_type delta=0;
  scalar_type lambda=0;
  scalar_type tau=0;
  scalar_type t=model->t_begin[branch]-model->t_end[branch];
  //obs. prob. being lost
  scalar_type p10=model->branch_counts["Ls"][branch]/model->branch_counts["count"][branch];
  //obs. avg. copy
  scalar_type k=(model->branch_counts["copies"][branch]-model->branch_counts["Ts"][branch])/model->branch_counts["count"][branch];
  if (k==1 or k==0 or p10==1  )
    {
      delta=scalar_parameter["min_delta"];
    }
  else 
    {
      delta=(p10 + k - 1)* log(k)/(1 - p10)/(k - 1)/t;
    }
  
  if (k==0 or p10==1)
    {
      lambda=model->branch_counts["Ls"][branch]/t;
    }
  else if (k==1 or p10==0 )
    {
      lambda=scalar_parameter["min_lambda"];
    }
  else
    {
      lambda=p10*k*log(k)/(1 - p10)/(k - 1)/t;
    }
  if (fT>0)
    tau=max(  ( model->branch_counts["Ts"][branch]/(N_ales_norm*fT) * lambda) / (1-exp(-lambda*t))  ,(scalar_type)scalar_parameter["min_tau"]);	
  else
    tau=(scalar_type)scalar_parameter["min_tau"];
  vector<scalar_type> dtl;
  scalar_type old_h=model->vector_parameter["delta"][branch]+model->vector_parameter["tau"][branch]+model->vector_parameter["lambda"][branch];
  scalar_type new_h=delta+tau+lambda;

  dtl.push_back( (new_h*w+(1-w)*old_h)/new_h*delta );
  dtl.push_back( (new_h*w+(1-w)*old_h)/new_h*tau );
  dtl.push_back( (new_h*w+(1-w)*old_h)/new_h*lambda );
  return dtl;
}

scalar_type mpi_tree::estimate_rates(string mode)
{
  
  scalar_type ll=-1e20;
  scalar_type old_ll=-2e20;
  vector <scalar_type> delta_last;//del-loc
  vector <scalar_type> tau_last;//del-loc
  vector <scalar_type> lambda_last;//del-loc
  map<string,vector<scalar_type> > last_counts;//del-loc
  vector<pair<string,scalar_type> > last_MLRec_res;//del-loc

  for (int branch=0;branch<model->last_branch;branch++)	
    {
      delta_last.push_back(0);
      tau_last.push_back(0);
      lambda_last.push_back(0);
      for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
	{
	  string count_name=(*it).first;
	  last_counts[count_name].push_back(0);
	}
    }
  int iters=0;
  while (ll>old_ll)
    {
      iters++;
      //calculation
      model->calculate_EG();
      old_ll=ll;
      ll=calculate_MLRecs();
      if (rank==server) {cout << "calc. done" <<endl;}
      //calculation

      //estimate rates
      if (rank==server) {cout << "estimate" <<endl;}
      
      vector <scalar_type> delta_estimate;//del-loc
      vector <scalar_type> tau_estimate;//del-loc
      vector <scalar_type> lambda_estimate;//del-loc
      if (rank==server)
	{
	  for (int branch=0;branch<model->last_branch;branch++)	
	    {
	      vector <scalar_type> estimate=dtl_estimate(branch,N_ales);	      
	      delta_estimate.push_back(estimate[0]);
	      tau_estimate.push_back(estimate[1]);	
	      lambda_estimate.push_back(estimate[2]);

	      estimate.clear();
	    }    
	}
      if (rank==server) {cout << "done" <<endl;}

      broadcast(world,delta_estimate,server);
      broadcast(world,tau_estimate,server);
      broadcast(world,lambda_estimate,server);
      
      //store previous rates
      if (ll>old_ll or iters<3)
	  {
	    if (rank==server) cout << "stored" <<endl;
	    show_rates();	  

	    for (int branch=0;branch<model->last_branch;branch++)	
	      {	      
		delta_last[branch]=model->vector_parameter["delta"][branch];
		tau_last[branch]=model->vector_parameter["tau"][branch]*model->vector_parameter["N"][0];
		lambda_last[branch]=model->vector_parameter["lambda"][branch];
		for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
		  {
		    string count_name=(*it).first;
		    last_counts[count_name][branch]=model->branch_counts[count_name][branch];
		  }
		last_MLRec_res.clear();
		for (int i=0;i< MLRec_res.size(); i++)
		  {
		    last_MLRec_res.push_back(MLRec_res[i]);
		  }
	      }
	    
	    //update rates            
	    if (mode=="uniform")
	      {
		scalar_type avg_delta=0;
		scalar_type avg_lambda=0;
		scalar_type avg_tau=0;
		scalar_type c=0;
		
		for (int branch=0;branch<model->last_branch;branch++)	
		  {
		    c+=1;
		    avg_delta+=delta_estimate[branch];
		    avg_tau+=tau_estimate[branch];
		    avg_lambda+=lambda_estimate[branch];
		  }
		avg_delta/=c;
		avg_lambda/=c;
		avg_tau/=c;
		
		model->set_model_parameter("delta",avg_delta);
		model->set_model_parameter("tau",avg_tau);
		model->set_model_parameter("lambda",avg_lambda);
		
	      }
	    else if (mode=="full_bw")
	      {
		model->set_model_parameter("delta",delta_estimate);
		model->set_model_parameter("tau",tau_estimate);      
		model->set_model_parameter("lambda",lambda_estimate);
	      }
	    show_rates();	  
	    show_branch_counts();
	  }
      //del-locs
      delta_estimate.clear();
      tau_estimate.clear();
      lambda_estimate.clear();
      if (rank==server) cout << "LL of est.  " << ll << " vs. previous " << old_ll<<endl;
    } 
    
  for (int branch=0;branch<model->last_branch;branch++)	
    for (map<string,vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
      {
	string count_name=(*it).first;
	model->branch_counts[count_name][branch]=last_counts[count_name][branch];	
      }
  for (int i=0;i< MLRec_res.size(); i++)
    {
      MLRec_res[i]=last_MLRec_res[i];
    }
  
  model->set_model_parameter("delta",delta_last);
  model->set_model_parameter("tau",tau_last);      
  model->set_model_parameter("lambda",lambda_last);
  if (rank==server) cout << "#last"<<endl;
  show_rates();	  
  model->calculate_EG();
  ll=old_ll;

  //del-locs

  delta_last.clear();
  tau_last.clear();
  lambda_last.clear();
  for (map<string,vector<scalar_type> >::iterator it=last_counts.begin();it!=last_counts.end();it++)
    (*it).second.clear();
  last_MLRec_res.clear();
  last_counts.clear();
  broadcast(world,ll,server);
  return ll;

}

void mpi_tree::show_rates()
{
  if (rank==server) 
    {
      // model->show_rates("delta");
      //model->show_rates("tau");
      //model->show_rates("lambda");
      scalar_type avg_delta=0;
      scalar_type avg_lambda=0;
      scalar_type avg_tau=0;
      scalar_type c=0;
      
      for (int branch=0;branch<model->last_branch;branch++)	
	{
	  c+=1;
	  avg_delta+=model->vector_parameter["delta"][branch];
	  avg_tau+=model->vector_parameter["tau"][branch];
	  avg_lambda+=model->vector_parameter["lambda"][branch];
	}
      avg_delta/=c;
      avg_lambda/=c;
      avg_tau/=c;
      cout << "# avg. delta= " << avg_delta;
      cout << " avg. tau= " << avg_tau*model->vector_parameter["N"][0];
      cout << " avg. lambda= " << avg_lambda;
      cout <<endl;
    }
  int done=1;
  broadcast(world,done,server);
}

void mpi_tree::set_rates(bool branchwise)
{

  if (!branchwise)
    {		
      vector<scalar_type> gathered_delta_avg,gathered_tau_avg,gathered_lambda_avg;
      vector<scalar_type> gathered_delta_norm,gathered_tau_norm,gathered_lambda_norm;
      
      gather(world,delta_avg,gathered_delta_avg,server);
      gather(world,tau_avg,gathered_tau_avg,server);
      gather(world,lambda_avg,gathered_lambda_avg,server);
      gather(world,delta_norm,gathered_delta_norm,server);
      gather(world,tau_norm,gathered_tau_norm,server);
      gather(world,lambda_norm,gathered_lambda_norm,server);
      if (rank==server)	
	{
	  delta_avg=0;
	  tau_avg=0;
	  lambda_avg=0;
	  delta_norm=0;
	  tau_norm=0;
	  lambda_norm=0;      
	  for (int i=0;i<(int)gathered_delta_avg.size();i++)
	    {
	      delta_avg+=gathered_delta_avg[i];
	      tau_avg+=gathered_tau_avg[i];
	      lambda_avg+=gathered_lambda_avg[i];
	      delta_norm+=gathered_delta_norm[i];
	      tau_norm+=gathered_tau_norm[i];
	      lambda_norm+=gathered_lambda_norm[i];      
	    }
	}
      gathered_delta_avg.clear();
      gathered_tau_avg.clear();
      gathered_lambda_avg.clear();
      gathered_delta_norm.clear();
      gathered_tau_norm.clear();
      gathered_lambda_norm.clear();

      broadcast(world,delta_avg,server);
      broadcast(world,tau_avg,server);
      broadcast(world,lambda_avg,server);
      broadcast(world,delta_norm,server);
      broadcast(world,tau_norm,server);
      broadcast(world,lambda_norm,server);
      
      model->set_model_parameter("delta",delta_avg/delta_norm);
      model->set_model_parameter("tau",tau_avg/tau_norm);
      model->set_model_parameter("lambda",lambda_avg/lambda_norm);
      //cout << rank << delta_avg/delta_norm << " " << tau_avg/tau_norm << " " << lambda_avg/lambda_norm << endl;
    }	    
  else 
    {
      
      vector<vector<scalar_type> > gathered_delta_avg,gathered_tau_avg,gathered_lambda_avg;
      vector<vector<scalar_type> > gathered_delta_norm,gathered_tau_norm,gathered_lambda_norm;
      
      gather(world,delta_branch_avg,gathered_delta_avg,server);
      gather(world,tau_branch_avg,gathered_tau_avg,server);
      gather(world,lambda_branch_avg,gathered_lambda_avg,server);
      gather(world,delta_branch_norm,gathered_delta_norm,server);
      gather(world,tau_branch_norm,gathered_tau_norm,server);
      gather(world,lambda_branch_norm,gathered_lambda_norm,server);

      if (rank==server)	
	{
	  for (int branch=0;branch<model->last_branch;branch++)	
	    {
	      delta_branch_avg[branch]=0;
	      tau_branch_avg[branch]=0;
	      lambda_branch_avg[branch]=0;
	      delta_branch_norm[branch]=0;
	      tau_branch_norm[branch]=0;
	      lambda_branch_norm[branch]=0;      
	      for (int i=0;i<(int)gathered_delta_avg.size();i++)
		{
		  delta_branch_avg[branch]+=gathered_delta_avg[i][branch];
		  tau_branch_avg[branch]+=gathered_tau_avg[i][branch];
		  lambda_branch_avg[branch]+=gathered_lambda_avg[i][branch];
		  delta_branch_norm[branch]+=gathered_delta_norm[i][branch];
		  tau_branch_norm[branch]+=gathered_tau_norm[i][branch];
		  lambda_branch_norm[branch]+=gathered_lambda_norm[i][branch];      
		}
	    }
	}

      broadcast(world,delta_branch_avg,server);
      broadcast(world,tau_branch_avg,server);
      broadcast(world,lambda_branch_avg,server);
      broadcast(world,delta_branch_norm,server);
      broadcast(world,tau_branch_norm,server);
      broadcast(world,lambda_branch_norm,server);

      vector<scalar_type> delta_estimate,tau_estimate,lambda_estimate;
      for (int branch=0;branch<model->last_branch;branch++)	
	{
	  delta_estimate.push_back(delta_branch_avg[branch]/delta_branch_norm[branch]);
	  tau_estimate.push_back(tau_branch_avg[branch]/tau_branch_norm[branch]);
	  lambda_estimate.push_back(lambda_branch_avg[branch]/lambda_branch_norm[branch]);
	}
      
      model->set_model_parameter("delta",delta_estimate);
      model->set_model_parameter("tau",tau_estimate);      
      model->set_model_parameter("lambda",lambda_estimate);
    }
  model->calculate_EG();
  show_rates();	  

}
void mpi_tree::clear_rate_register()
{
  for (int branch=0;branch<model->last_branch;branch++)	
    {
      delta_branch_avg.push_back(0);
      delta_branch_norm.push_back(0);
      tau_branch_avg.push_back(0);
      tau_branch_norm.push_back(0);
      lambda_branch_avg.push_back(0);
      lambda_branch_norm.push_back(0);	        
    }
  delta_avg=0;
  delta_norm=0;
  tau_avg=0;
  tau_norm=0;
  lambda_avg=0;
  lambda_norm=0;	        
}
void mpi_tree::register_rates()
{

  for (int branch=0;branch<model->last_branch;branch++)	
    {
      scalar_type t_branch=model->t_end[branch];
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
	  delta=scalar_parameter["min_delta"];
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
	  lambda=scalar_parameter["min_lambda"];
	}
      else
	{
	  lambda=p10*k*log(k)/(1 - p10)/(k - 1)/t;
	}
      //tau=max(  (scalar_type)( model->branch_counts["Ts"][branch]/(1.) * lambda) / (1-exp(-lambda*t))  ,(scalar_type)scalar_parameter["min_tau"]);	
      tau=max(  (scalar_type)( model->branch_counts["Ts"][branch]/t)  ,(scalar_type)scalar_parameter["min_tau"]);	
      if(model->branch_counts["count"][branch]>0)
	{
	  delta_avg+=delta;
	  delta_norm+=1;
	  delta_branch_avg[branch]+=delta;
	  delta_branch_norm[branch]+=1;
	  
	  lambda_avg+=lambda;
	  lambda_norm+=1;
	  lambda_branch_avg[branch]+=lambda;
	  lambda_branch_norm[branch]+=1;
	}
      tau_avg+=tau;     
      tau_norm+=1;
      tau_branch_avg[branch]+=tau;      
      tau_branch_norm[branch]+=1;
      
    }
  for (map<string, vector<scalar_type> >::iterator it=model->branch_counts.begin();it!=model->branch_counts.end();it++)
    for ( vector<scalar_type>::iterator jt=(*it).second.begin();jt!=(*it).second.end();jt++)
      (*jt)=0;
  
     
}
