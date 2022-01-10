#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Lapack.h>




#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()


void score_interf(int *N, int *P, double *Y, int *Tscore, int *stdize, double *Thre, double *EPS, int *STEP_MAX, int *PRINT, int *ini_adj_matrix, int *B_ini_adj, int *blacklist, int *whitelist, int *N_boot, int *stdize_boot, int *RANDOM_forest, double *rand_step_length, int *Restart, int *Perturb, int *movement, int *adj_matrix, double *delta_min, int *bic_OK, int *final_step);
void copy_data(int p, int *cur_par, int *num_par, int *cur_par_t, int *num_par_t, int *cur_child, int *num_child, int *cur_child_t, int *num_child_t, int *acyc_all, int *acyc_all_t, int *ini_adj_matrix, int *adj_temp);
void whitelist_blacklist(int *P, int *STEP_max, int *adj_matrix, int *movement, int *whitelist, int *blacklist, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child, int *Step_rec);
void adj_ini_fn(int *P, int *adj_matrix, int *ini_adj_matrix, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child);
void random_restart(int p, int n, int step_max, int step_rec, double thre, double eps, int nrestart, int nperturb, double *Y, int *Tscore, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *acyc_all, double *bic_delta_all, double *cur_BIC, int *adj_matrix, int *movement, double *delta_min, int *bic_OK, int *best_movement, int *best_adj, int *final_step);

void score_step(int *N, int *P, double *Y, int *Tscore, double *Thre, double *EPS, int *STEP_MAX, int *PRINT, int *last_move, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *acyc_all, double *bic_delta_all, int *RANDOM_forest, double *rand_step_length, int *movement, int *adj_matrix, double *delta_min, double *cur_BIC, int *bic_OK, int *Step_rec, int *final_step);
void oper_adj(int p, int *opt_oper, int *acyc_all, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *adj_matrix);
void operation_list(int *N, int *P, double *Y, int *Tscore, double *EPS,  int *adj_matrix, int *last_oper, int *acyc_all, double *bic_delta_all, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, int *opt_oper, double *MIN_Delta, int *NUM_oper);
void oper_add(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int *opt_oper, double *MIN_Delta, int *NUM_oper);
void oper_del(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int *opt_oper, double *MIN_Delta, int *NUM_oper);
void oper_rev(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int *opt_oper, double *MIN_Delta, int *NUM_oper);

void score_type_bic(int n, int p, double *Y, int tscore, int recal_ind, int *cand_oper, double *bic_delta_all, double *cur_BIC, int *cur_par, int *num_par, int *opt_oper, double *oper_BIC, double *Min_delta, double eps);
void comp_delta(double *Delta_temp, double *Min_delta, double *EPS,int *temp_oper, int *opt_oper);
void find_path(int *P, int *Node, int *cur_child, int *num_child, int *res, int *res_length);

void bic_score(int *N, int *P, double *Y, int *Tscore, int *oper, double *cur_BIC,  int *cur_par, int *num_par, double *node_BIC, double * bic_delta, double *EPS, int *info);
void bic_add(int *N, int *P, int *Fnode, int *Tnode, double *Y, int *Tscore, int *cur_par, int *num_par, double *bic_res, int *info);
void bic_del(int *N, int *P, int *Fnode, int *Tnode, double *Y, int *Tscore, int *cur_par, int *num_par, double *bic_res, int *info);
//void bic_cal(int *x_nrow, int *x_ncol, double *X, double *y, int *Tscore, double *bic_res, int *ok);
void bic_cal(int *x_nrow, int *x_ncol, int *Npar, double *X,double *y, int *Tscore, double *bic_res, int *info);
void bic_par(int *N, int *P, double *Y, int *Tscore, int *cur_par, int *num_par, double *cur_BIC);

void check_acyclic(int *P, int *adj_matrix, int *cur_child_copy, int *num_child_copy, int *oper, int *res);
void del_parchild(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child);
void add_parchild(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child);
void acyclic_add(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child, int *acyclic);
void SampleNoReplace(int k, int n, int *y, int *x);
void sampling(int *NUM, int *index, int *n_sample, double *prob, int *res);
void randforest(int *P, int *NUM_oper, double *MIN_delta, double *RAND_step_length, int *acyc_all, double *bic_delta_all, double *EPS, int *cur_step, int *LL, double *MIN_valve, int *e_oper);
void mean(int *N, double *x, double *res);
void data_stdize(int *N, int *P, double *x, double *x_std);


/////////////////////////////////////////////////
///////////////////// hill climbing DAG learning
void score_interf(int *N, int *P, double *Y, int *Tscore, int *stdize, double *Thre, double *EPS, int *STEP_MAX, int *PRINT, int *ini_adj_matrix, int *B_ini_adj, int *blacklist, int *whitelist, int *N_boot, int *stdize_boot, int *RANDOM_forest, double *rand_step_length, int *Restart, int *Perturb, int *movement, int *adj_matrix, double *delta_min, int *bic_OK, int *final_step){


  int i,j,k,temp_row,s=0,p=*P, n=*N, step_max=*STEP_MAX, pr=*PRINT, n_boot=*N_boot, b_ini_adj=*B_ini_adj;
  double thre=*Thre, eps=*EPS, min_delta=1e+06, prob=-1.0;
  int rf=*RANDOM_forest, step_rec=0, nrestart=*Restart, nperturb=*Perturb, last_move[3], tscore=*Tscore;

  last_move[0]=-1;
  last_move[1]=-1;
  last_move[2]=-1;

  if(*stdize==1){
      data_stdize(&n,&p,Y,Y);
  }

	int *adj_t=(int *) malloc(p*p*sizeof(int));
	int *cur_par=(int *)malloc(p*p*sizeof(int));
	int *cur_child=(int *)malloc(p*p*sizeof(int));
	int *num_par=(int *)malloc(p*sizeof(int));
	int *num_child=(int *)malloc(p*sizeof(int));
	int *acyc_all=(int *)malloc(p*p*3*sizeof(int));
	double *cur_BIC=(double *)malloc(p*sizeof(double));

    double *bic_delta_all=(double *)malloc(p*p*3*sizeof(double));

	for(i=0;i<p;i++){
		num_child[i]=0;
		num_par[i]=0;
		for(j=0;j<p;j++){
			cur_par[j*p+i]=0;
			cur_child[j*p+i]=0;
			adj_t[j*p+i]=0;
		}
	}

	whitelist_blacklist(&p, &step_max, adj_t, movement, whitelist, blacklist, acyc_all, cur_par, num_par, cur_child, num_child, &step_rec);
	if(b_ini_adj==1){
		adj_ini_fn(&p, adj_t, ini_adj_matrix, acyc_all, cur_par, num_par, cur_child, num_child);
	}

  if(n_boot>0){
    double *y_boot=(double *) malloc(n*p*sizeof(double));
    int *boot_ind=(int *) malloc(n*sizeof(int));
    int *ori_ind=(int *) malloc(n*sizeof(int));
    double *y_std=(double *)malloc(n*p*sizeof(double));

	int *adj_temp=(int *) malloc(p*p*sizeof(int));
	int *cur_par_t=(int *)malloc(p*p*sizeof(int));
	int *cur_child_t=(int *)malloc(p*p*sizeof(int));
	int *num_par_t=(int *)malloc(p*sizeof(int));
	int *num_child_t=(int *)malloc(p*sizeof(int));
	int *acyc_all_t=(int *)malloc(p*p*3*sizeof(int));
	int *movement_t=(int *)malloc(step_max*3*sizeof(int));

    for(j=0;j<n;j++){
      ori_ind[j]=j;
    }

	  //start boostrapping
    for(i=0;i<n_boot;i++){

	  last_move[0]=-1;
	  last_move[1]=-1;
	  last_move[2]=-1;

      sampling(&n, ori_ind, &n, &prob, boot_ind);
      for(j=0;j<n;j++){
        temp_row=boot_ind[j];
        for(k=0;k<p;k++){
          y_boot[k*n+j]=Y[k*n+temp_row];
        }
      }//end bootstrapping Y

     if(*stdize_boot==1){
       data_stdize(&n,&p,y_boot,y_boot);
     }

	copy_data(p, cur_par, num_par, cur_par_t, num_par_t, cur_child, num_child, cur_child_t, num_child_t, acyc_all, acyc_all_t, adj_t, adj_temp);
	bic_par(&n, &p, y_boot, &tscore, cur_par_t, num_par_t, cur_BIC);

	score_step(&n, &p, y_boot, &tscore, &thre, &eps, &step_max, &pr, last_move, blacklist, cur_par_t, num_par_t, cur_child_t, num_child_t, acyc_all_t, bic_delta_all, &rf, rand_step_length, movement, adj_temp, delta_min, cur_BIC, bic_OK, &step_rec, &s);

	 for(j=0;j<p;j++){
	   for(k=0;k<p;k++){
		 adj_matrix[i*p*p+j*p+k]=adj_temp[j*p+k];
	   }
	 }

   }//end n_boot loop


   free(acyc_all_t);
   free(cur_child_t);
   free(cur_par_t);
   free(num_child_t);
   free(num_par_t);
   free(y_boot);
   free(y_std);
   free(adj_temp);
   free(boot_ind);
   free(ori_ind);
   free(movement_t);
  }//end if

  else{//no bootstrap
	  int m,b_feasible=1;

	  for(j=0;j<p;j++){
		  for(k=0;k<p;k++){
			  adj_matrix[j*p+k]=adj_t[j*p+k];
		  }//end for(k)
	  }//end for(j)

	bic_par(&n, &p, Y, &tscore, cur_par, num_par, cur_BIC);
    score_step(&n, &p, Y, &tscore, &thre, &eps, &step_max, &pr, last_move, blacklist, cur_par, num_par, cur_child, num_child, acyc_all, bic_delta_all, &rf, rand_step_length, movement, adj_matrix, delta_min, cur_BIC, bic_OK, &step_rec, &s);

	  //random start, 11-09-2012
	if(s>step_max||(s+nperturb)>step_max){
		b_feasible=0;
///6/23		printf("Not enough steps to do random restart! Larger step.max is needed\n");
	}
	if((nrestart>0)&&(b_feasible==1)){

		int *best_movement=(int *)malloc(step_max*3*sizeof(int));
		int *best_adj=(int *)malloc(p*p*sizeof(int));

		for(i=0;i<step_max;i++){
		  best_movement[i]=movement[i];
		  best_movement[step_max+i]=movement[step_max+i];
		  best_movement[2*step_max+i]=movement[2*step_max+i];
		}//end for i
		for(j=0;j<p;j++){
		  for(i=0;i<p;i++){
			best_adj[j*p+i]=adj_matrix[j*p+i];
		  }//end for i
		}//end for j
		random_restart(p, n, step_max, s, thre, eps, nrestart, nperturb, Y, &tscore, blacklist, cur_par, num_par, cur_child, num_child, acyc_all, bic_delta_all, cur_BIC, adj_matrix, movement, delta_min, bic_OK, best_movement, best_adj, &s);

		for(i=0;i<step_max;i++){
			movement[i]=best_movement[i];
			movement[step_max+i]=best_movement[step_max+i];
			movement[2*step_max+i]=best_movement[2*step_max+i];
		}//end for i

		for(i=0;i<p;i++){
			for(j=0;j<p;j++){
				adj_matrix[j*p+i]=best_adj[j*p+i];
			}//end for j
		}//end i

		free(best_movement);
		free(best_adj);
	}

	*final_step=s;
  }//end else

    free(cur_BIC);
	free(adj_t);
	free(acyc_all);
	free(cur_child);
	free(cur_par);
	free(num_child);
	free(num_par);
	free(bic_delta_all);

}//r_interf

void copy_data(int p, int *cur_par, int *num_par, int *cur_par_t, int *num_par_t, int *cur_child, int *num_child, int *cur_child_t, int *num_child_t, int *acyc_all, int *acyc_all_t, int *ini_adj_matrix, int *adj_matrix){
	int i,j,k;
	for(j=0;j<p;j++){
		num_par_t[j]=num_par[j];
		num_child_t[j]=num_child[j];

		for(i=0;i<p;i++){
			cur_par_t[j*p+i]=cur_par[j*p+i];
			cur_child_t[j*p+i]=cur_child[j*p+i];
			for(k=0;k<3;k++){
				acyc_all_t[k*p*p+j*p+i]=acyc_all[k*p*p+j*p+i];
			}
			adj_matrix[j*p+i]=ini_adj_matrix[j*p+i];
		}

	}

}//end fn copy_data


void random_restart(int p, int n, int step_max, int step_rec, double thre, double eps, int nrestart, int nperturb, double *Y, int *Tscore, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *acyc_all, double *bic_delta_all, double *cur_BIC, int *adj_matrix, int *movement, double *delta_min, int *bic_OK, int *best_movement, int *best_adj, int *final_step){
    int *adj_t=(int *)malloc(p*p*sizeof(int));
	int *cur_par_t=(int *)malloc(p*p*sizeof(int));
	int *cur_child_t=(int *)malloc(p*p*sizeof(int));
	int *num_par_t=(int *)malloc(p*sizeof(int));
	int *num_child_t=(int *)malloc(p*sizeof(int));
	int *acyc_all_t=(int *)malloc(p*p*3*sizeof(int));
	int *movement_t=(int *)malloc(step_max*3*sizeof(int));
	double *cur_BIC_t=(double *)malloc(p*sizeof(double));
	double *best_BIC=(double *)malloc(p*sizeof(double));
	int *temp_rand_oper=(int *)malloc(3*sizeof(int));
	double *bic_all_t=(double *)malloc(p*p*3*sizeof(double));
    int *index_p=(int *)malloc(p*sizeof(int));
	int m, i, j, k, s, best_s, nodev=2, nodes[2], count, temp, temp2, num_oper,pr=0,rf=0;
	int rand_oper, index_oper[3],operv=1,opertotal=3, opt_oper[3], tscore=*Tscore;

	double min_delta, mean_bic_1, mean_bic_2,*rand_step_length;

	for(j=0;j<p;j++){
	   best_BIC[j]=cur_BIC[j];
    }
	mean(&p, best_BIC, &mean_bic_2);

	///start
	for(m=0;m<nrestart;m++){
		copy_data(p, cur_par, num_par, cur_par_t, num_par_t, cur_child, num_child, cur_child_t, num_child_t, acyc_all, acyc_all_t, adj_matrix, adj_t);

		for(i=0;i<step_max;i++){
			movement_t[i]=movement[i];
			movement_t[step_max+i]=movement[step_max+i];
			movement_t[2*step_max+i]=movement[2*step_max+i];
		}

		for(j=0;j<p;j++){
		  for(i=0;i<p;i++){
			for(k=0;k<3;k++){
			  bic_all_t[k*p*p+j*p+i]=bic_delta_all[k*p*p+j*p+i];
			}
		  }
		  cur_BIC_t[j]=cur_BIC[j];
		}

		///choose nperturb operation to be perturbed
		temp2=step_rec;
		if(movement[temp2]==0){
			temp2=step_rec-1;
		}
		else{
			temp2=step_rec;
		}

        count=0;
		while(count<nperturb){
			SampleNoReplace(nodev, p, nodes,index_p);
			SampleNoReplace(operv, opertotal, &rand_oper, index_oper);
			temp=rand_oper*p*p+nodes[1]*p+nodes[0];
			if(acyc_all_t[temp]==0){
			  count=count+1;

			  temp_rand_oper[0]=nodes[0]+1;
			  temp_rand_oper[1]=nodes[1]+1;
			  temp_rand_oper[2]=rand_oper+1;

	          oper_adj(p, temp_rand_oper, acyc_all_t, blacklist, cur_par_t, num_par_t, cur_child_t, num_child_t, adj_t);

			  movement_t[temp2]=temp_rand_oper[0];
			  movement_t[step_max+temp2]=temp_rand_oper[1];
			  movement_t[2*step_max+temp2]=temp_rand_oper[2];

			  temp2=temp2+1;

				if(temp_rand_oper[2]==1){
					oper_add(&n, &p, Y, &tscore, adj_t, temp_rand_oper, cur_par_t, num_par_t, cur_child_t, num_child_t, cur_BIC_t, &eps, acyc_all_t, bic_all_t,opt_oper,&min_delta,&num_oper);
				}//end addition

				//deletion operations
				if(temp_rand_oper[2]==2){
					oper_del(&n, &p, Y, &tscore, adj_t, temp_rand_oper, cur_par_t, num_par_t, cur_child_t, num_child_t, cur_BIC_t, &eps, acyc_all_t, bic_all_t,opt_oper,&min_delta,&num_oper);
				}//end if(last_move[2]==2)

				//reversal operations
				if(temp_rand_oper[2]==3){
					oper_rev(&n, &p, Y, &tscore, adj_t, temp_rand_oper, cur_par_t, num_par_t, cur_child_t, num_child_t, cur_BIC_t, &eps, acyc_all_t, bic_all_t,opt_oper,&min_delta,&num_oper);
				}//end if(last_move[2]==3)

			}//end if(acyc_all_t[temp]==0)

		}//end while(count<nperturb)

		bic_par(&n, &p, Y, &tscore, cur_par_t, num_par_t, cur_BIC);
		score_step(&n, &p, Y, &tscore, &thre, &eps, &step_max, &pr, temp_rand_oper, blacklist, cur_par_t, num_par_t, cur_child_t, num_child_t, acyc_all_t, bic_all_t, &rf, rand_step_length, movement_t, adj_t, delta_min, cur_BIC_t, bic_OK, &temp2, &s);

		mean(&p, cur_BIC_t, &mean_bic_1);
		if(mean_bic_1<mean_bic_2){//if better, copy, adjacency matrix, movement
			for(i=0;i<p;i++){
				for(j=0;j<p;j++){
					best_adj[j*p+i]=adj_t[j*p+i];
				}
				best_BIC[i]=cur_BIC_t[i];
			}
			mean_bic_2=mean_bic_1;

			for(i=0;i<step_max;i++){
				best_movement[i]=movement_t[i];
				best_movement[step_max+i]=movement_t[step_max+i];
				best_movement[2*step_max+i]=movement_t[2*step_max+i];
			}
			best_s=s;
		}
		//
	}//end for loop

	free(temp_rand_oper);
	free(cur_BIC_t);
	free(movement_t);
	free(acyc_all_t);
	free(num_child_t);
	free(num_par_t);
	free(cur_child_t);
	free(cur_par_t);
	free(best_BIC);
	free(adj_t);
	free(bic_all_t);
	free(index_p);

	*final_step=best_s;

}//end function

//main function
void score_step(int *N, int *P, double *Y, int *Tscore, double *Thre, double *EPS, int *STEP_MAX, int *PRINT, int *last_move, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *acyc_all, double *bic_delta_all, int *RANDOM_forest, double *rand_step_length, int *movement, int *adj_matrix, double *delta_min, double *cur_BIC, int *bic_OK, int *Step_rec, int *final_step){
  //written on 02-08-2012
  //main function for DAG constructure using BIC score
  //bic_delta_all: a p by p by 3 array, recording the bic delta for addition, deletion and reversal
  //acyc_all: a p by p by 3 array, recording the acyclic information for each operation
  //cur_par(cur_child): a p by p matrix, recording the parents (children) indices for each node,
  //num_par(num_child): a p-length vector, recording the number of parents (children) for each node
  //movement: a STEP_MAX by 3 matrix
  //adj_matrix: adjacency matrix
  //adj_result: a p by p by STEP_MAX array
  //score_result: a STEP_MAX by p matrix, recording the BIC score vector for each step
  //delta_min: a STEP_MAX-length vector, recording the minimum delta for each step
  //bic_OK: a STEP_MAX-length vector, recording whether BIC is calcuated correctly

  int i,j=0,k=0, itemp, jtemp, s, step_rec=*Step_rec, ok=0, p=*P, n=*N, step_max=*STEP_MAX, pr=*PRINT, rf=*RANDOM_forest;
  int opt_oper[3], num_oper=0, LL=p*p*3, tscore=*Tscore;
  double thre=*Thre, eps=*EPS, min_delta=1e+06, tnew_BIC_temp[2], delta_temp=0.0, min_valve=0.0,delta_diff=0.0;

//  double *cur_BIC=(double *)malloc(p*sizeof(double));
//  bic_par(&n, &p, Y, cur_par, num_par, cur_BIC);
  //start DAG
  for(s=0;s<step_max;s++){
    if(pr==1){
 ////6/23     printf("step -- %d \n",s);
    }//end if

    opt_oper[0]=0;
    opt_oper[1]=0;
    opt_oper[2]=0;

    //acyc_all updated, bic_delta_all updated, return opt_oper, and min_delta
    operation_list(&n, &p, Y, &tscore, &eps, adj_matrix, last_move, acyc_all, bic_delta_all, cur_par, num_par, cur_child, num_child, cur_BIC, opt_oper, &min_delta, &num_oper);

    //random forest idea, added on 03-01-2012
    if((min_delta<thre)&&(rf==1)&&(num_oper>0)){
        randforest(&p, &num_oper, &min_delta, &rand_step_length[s], acyc_all, bic_delta_all, &eps, &s, &LL, &min_valve, opt_oper);
		//min_delta may be changed because of random forest, here min_delta is set to be the corresponding bic delta value for opt_oper selected.
    }//end if(random_forest==1)

    if(min_delta>=thre){
      opt_oper[0]=0;
      opt_oper[1]=0;
      opt_oper[2]=0;
      delta_min[s]=0.0;
      bic_OK[s]=0;
    }

    //because we only store the BIC delta, not BIC scores, after choosing optimal operation, BIC needs to be calcuated to update cur_BIC vector
    if(opt_oper[2]>0){// not "no movement"
      bic_score(&n, &p, Y, &tscore, opt_oper, cur_BIC, cur_par, num_par, tnew_BIC_temp, &delta_temp, &eps, &ok);
      cur_BIC[opt_oper[0]-1]=tnew_BIC_temp[0];
      cur_BIC[opt_oper[1]-1]=tnew_BIC_temp[1];
      delta_min[s]=delta_temp;
      bic_OK[s]=ok;
    }

// update, cur_par, num_par, cur_child, num_child
// update movment, adj_result and score_result
    oper_adj(p, opt_oper, acyc_all, blacklist, cur_par, num_par, cur_child, num_child, adj_matrix);

    //record optimal movement
    movement[step_rec]=opt_oper[0];
    movement[step_max+step_rec]=opt_oper[1];
    movement[2*step_max+step_rec]=opt_oper[2];
    step_rec=step_rec+1;
    //update last movement
    last_move[0]=opt_oper[0];
    last_move[1]=opt_oper[1];
    last_move[2]=opt_oper[2];
    //////////////////////////////////////////////////////////////

    if((opt_oper[2]==0)){
 ////6/23     printf("DAG established, no more improvement");
      break;
    }

  }//end for(step)

  *final_step=step_rec;


}///end function dag_step()

void oper_adj(int p, int *opt_oper, int *acyc_all, int *blacklist, int *cur_par, int *num_par, int *cur_child, int *num_child, int *adj_matrix){

	if(opt_oper[2]==1){//current optimal is addition
		add_parchild(&p, &opt_oper[0], &opt_oper[1], cur_child, num_child);
		add_parchild(&p, &opt_oper[1], &opt_oper[0], cur_par, num_par);
		adj_matrix[(opt_oper[1]-1)*p+opt_oper[0]-1]=1;

		acyc_all[(opt_oper[1]-1)*p+(opt_oper[0]-1)]=99;//cannot add the same edge again
		acyc_all[p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=0;//deletion is valid now

		if(blacklist[(opt_oper[0]-1)*p+(opt_oper[1]-1)]==0){//if the reverse edge is not in blacklist
			acyc_all[2*p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=0;//reversal is possible valid now, but has not been considered before this operation
		}
    }

    if(opt_oper[2]==2){//deletion
		del_parchild(&p, &opt_oper[0], &opt_oper[1], cur_child, num_child);
		del_parchild(&p, &opt_oper[1], &opt_oper[0], cur_par, num_par);
		adj_matrix[(opt_oper[1]-1)*p+opt_oper[0]-1]=0;

		acyc_all[(opt_oper[1]-1)*p+(opt_oper[0]-1)]=0;// add this edge is possible valid
		acyc_all[p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=-99;// edge does not exist

		if(blacklist[(opt_oper[0]-1)*p+(opt_oper[1]-1)]==0){//if the reverse edge is not in blacklist, otherwise stick to -1e+06 value
			acyc_all[2*p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=-99;//meanless
		}
    }

    if(opt_oper[2]==3){//reversal
		del_parchild(&p, &opt_oper[0], &opt_oper[1], cur_child, num_child);
		del_parchild(&p, &opt_oper[1], &opt_oper[0], cur_par, num_par);
		add_parchild(&p, &opt_oper[1], &opt_oper[0], cur_child, num_child);
		add_parchild(&p, &opt_oper[0], &opt_oper[1], cur_par, num_par);

		adj_matrix[(opt_oper[1]-1)*p+opt_oper[0]-1]=0;
		adj_matrix[(opt_oper[0]-1)*p+opt_oper[1]-1]=1;

		acyc_all[(opt_oper[1]-1)*p+(opt_oper[0]-1)]=99;//acyclic
		acyc_all[p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=-99;//meanless
		acyc_all[2*p*p+(opt_oper[1]-1)*p+(opt_oper[0]-1)]=-99;//meanless

		acyc_all[(opt_oper[0]-1)*p+(opt_oper[1]-1)]=-99;//cannot add the same edge again
		acyc_all[p*p+(opt_oper[0]-1)*p+(opt_oper[1]-1)]=0;//valid
		acyc_all[2*p*p+(opt_oper[0]-1)*p+(opt_oper[1]-1)]=0; //possible valid
    }

}

void operation_list(int *N, int *P, double *Y, int *Tscore, double *EPS,  int *adj_matrix, int *last_oper, int *acyc_all, double *bic_delta_all, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, int *opt_oper, double *MIN_Delta, int *NUM_oper){
  //modified on 02-08-2012
  //added on 01-25-2012
  //last_oper: lastest operation, -1:no previous step
  //acyc_all: p by p by 3 array recording acyaclic information, where 3 is addition, deletion and reversal
  //bic_delta_all: p by p by 3 array recording updating delta of BIC, 3 is addition, deletion and reversal
  //cur_par, cur_child are p by p matrices
  //num_par, num_child are p-length vector, recording how many parents and children for each node
  //cur_BIC: p-length vector, recording the BIC scores for the p nodes
  //opt_oper: the best movement in the current context, 3-length vector
  //MIN_Delta: the minimum delta of all possible operations
  //return: opt_oper and MIN_Delta
  int n=*N, p=*P, num_oper=0, tscore=*Tscore;
  double eps=*EPS, min_delta=1e+06;

  if(last_oper[2]==-1){//first operation, do all additions
    int i,j,k,ok;
    int temp_oper[3];
    double tnew_BIC_temp[2];
    double delta_temp;

    temp_oper[2]=1;
    for(j=0;j<p;j++){
      temp_oper[1]=j+1;
      for(i=0;i<p;i++){
       if(j!=i){
        temp_oper[0]=i+1;//(i+1)-->(j+1)
        delta_temp=1e+06;

        bic_score(&n, &p, Y, &tscore, temp_oper, cur_BIC, cur_par, num_par, tnew_BIC_temp, &delta_temp, &eps, &ok);
        bic_delta_all[j*p+i]=delta_temp;//bic_delta_all initiation
        if(acyc_all[j*p+i]==0){
          comp_delta(&delta_temp, &min_delta, &eps, temp_oper, opt_oper);
          num_oper=num_oper+1;
        }
       }//end if(j!=i)
      }//end for
    }//end for

  }//end if(last_oper==-1)

  else{
    //addition operations
    if(last_oper[2]==1){
      oper_add(&n, &p, Y, &tscore, adj_matrix, last_oper, cur_par, num_par, cur_child, num_child, cur_BIC, &eps, acyc_all, bic_delta_all,opt_oper,&min_delta,&num_oper);
    }//end addition

  //deletion operations
    if(last_oper[2]==2){
      oper_del(&n, &p, Y, &tscore, adj_matrix, last_oper, cur_par, num_par, cur_child, num_child, cur_BIC, &eps, acyc_all, bic_delta_all,opt_oper,&min_delta,&num_oper);
    }//end if(last_move[2]==2)

  //reversal operations
    if(last_oper[2]==3){
      oper_rev(&n, &p, Y, &tscore, adj_matrix, last_oper, cur_par, num_par, cur_child, num_child, cur_BIC, &eps, acyc_all, bic_delta_all,opt_oper,&min_delta,&num_oper);
    }//end if(last_move[2]==3)

   }//end else--last_oper!=-1

   *MIN_Delta=min_delta;
   *NUM_oper=num_oper;
}//end function opeartion_list


void oper_add(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int *opt_oper, double *MIN_Delta, int *NUM_oper){
  //written on 02-12-2012
  //if last operation is addition, update acyclic information, bic_delta, choosing optimal operation, and minimum delta
  //when last_operation is addition, then if an operation creates a cycle before the last operation, it will create a cycle still.
  int n=*N, p=*P, temp_oper[3], tscore=*Tscore;
  temp_oper[0]=0;
  temp_oper[1]=0;
  temp_oper[2]=0;

  double oper_BIC[2], delta_temp=1e+06, min_delta=1e+06, eps=*EPS;
  int i,j,k,ok, recal_ind, cur_acyc, acyc_res, num_oper=0, last_from=last_oper[0], last_to=last_oper[1];

  int *from_am=(int *) malloc(p*sizeof(int));
  int *from_dm=(int *) malloc(p*sizeof(int));
  int *to_am=(int *) malloc(p*sizeof(int));
  int *to_dm=(int *) malloc(p*sizeof(int));

  int from_ancestor=0, from_descendant=0, to_ancestor=0, to_descendant=0;

  find_path(&p, &last_from, cur_par, num_par, from_am, &from_ancestor);// find last_from's ancestor
  find_path(&p, &last_from, cur_child, num_child, from_dm, &from_descendant);// find last_from's descendant
  find_path(&p, &last_to, cur_par, num_par, to_am, &to_ancestor);//find last_to's ancestor
  find_path(&p, &last_to, cur_child, num_child, to_dm, &to_descendant);//find last_to's descendant

  for(k=0;k<3;k++){//k loop
    temp_oper[2]=k+1;//operations: 1,2,3
	for(j=0;j<p;j++){//j loop,to node
	  temp_oper[1]=j+1;
	  for(i=0;i<p;i++){//i loop, from node
	    if(i!=j){
	      temp_oper[0]=i+1;
          recal_ind=1;//0--re-calculate, 1--OK and no need to calculate, 1e+06--cycle and no need to calculate
          cur_acyc=acyc_all[k*p*p+j*p+i];

	      if(cur_acyc==0){//last time no cycle
	      //k: current operation under consideration, k==0 means addition
	        if((k==0)&&(from_am[j]==1)&&(to_dm[i]==1)){//create cycle
			  acyc_all[k*p*p+j*p+i]=99;
			  recal_ind=1e+06;
			}
	        else if((k==0)&&(last_to==j+1)){//does not create cycle, but the neighborhood is changed
		      recal_ind=0;
			}

			if((k==1)&&(last_to==j+1)){//deletion and change neighborhood
			  recal_ind=0;
			}

	        if((k==2)&&(from_am[i]==1)&&(to_dm[j]==1)){//reveral, create cycles
			  if((last_from==i+1)&&(last_to==j+1)){//last_from->last_to reversal
			    check_acyclic(&p, adj_matrix, cur_child, num_child, temp_oper, &acyc_res);
				if(acyc_res==0){
				  recal_ind=0;
				}
				else{
				  acyc_all[k*p*p+j*p+i]=99;
				  recal_ind=1e+06;
				}
			  }
			  else{
				acyc_all[k*p*p+j*p+i]=99;
				recal_ind=1e+06;
			  }
			}
			else if((k==2)&&((last_to==i+1)||(last_to==j+1))){
			  recal_ind=0;
			}

			if(acyc_all[k*p*p+j*p+i]==0){
				score_type_bic(n, p, Y, tscore, recal_ind, temp_oper, bic_delta_all, cur_BIC, cur_par, num_par, opt_oper, oper_BIC, &min_delta, eps);
				num_oper=num_oper+1;
			}
	      }//end if(cur_acyc==0)

		}//end if(j!=i)
	  } //end for(i)
	}//end for(j)
  }//end for(k)


   *MIN_Delta=min_delta;
   *NUM_oper=num_oper;

   free(from_am);
   free(from_dm);
   free(to_am);
   free(to_dm);

 }//end oper_add()


void score_type_bic(int n, int p, double *Y, int tscore, int recal_ind, int *cand_oper, double *bic_delta_all, double *cur_BIC, int *cur_par, int *num_par, int *opt_oper, double *oper_BIC, double *Min_delta, double eps){

	double delta_temp, min_delta=*Min_delta;
	int k=cand_oper[2]-1, j=cand_oper[1]-1,i=cand_oper[0]-1, ok;

	if(recal_ind==0){//re-calculate bic_delta
		bic_score(&n, &p, Y, &tscore, cand_oper, cur_BIC, cur_par, num_par, oper_BIC, &delta_temp, &eps, &ok);
		bic_delta_all[k*p*p+j*p+i]=delta_temp;
		comp_delta(&delta_temp, &min_delta, &eps, cand_oper, opt_oper);
	}//end if(last_oper[1]==j)
	else if(recal_ind==1){//good operation but no need to re-calculate bic_delta, using recorded bic_delta
		delta_temp=bic_delta_all[k*p*p+j*p+i];
		comp_delta(&delta_temp, &min_delta, &eps, cand_oper, opt_oper);
	}

	*Min_delta=min_delta;
}

void oper_del(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int *opt_oper, double *MIN_Delta, int *NUM_oper){
  //if last operation is deletion, update acyclic information, bic_delta, choosing optimal operations, and min_delta
  //if last operatin is deletion, if an operation does not create cycle before last operation, then it will not create a cycle still;
  //if an opertion creates a cycle before last operation, then it still creates a cycle after last operation, or it stops creating a cycle if the cycle is breaked by last the deletion
  //  printf("Last Operation is Deletion \n");

  int n=*N, p=*P, temp_oper[3], tscore=*Tscore;
  temp_oper[0]=0;
  temp_oper[1]=0;
  temp_oper[2]=0;
  double oper_BIC[2], delta_temp=1e+06, min_delta=1e+06, eps=*EPS;

  int i,j,k,ok;

  int last_from=last_oper[0], last_to=last_oper[1];
  int cur_acyc, recal_ind, acyc_res,num_oper=0, temp1, temp2;

  int *from_am=(int *) malloc(p*sizeof(int));
  int *from_dm=(int *) malloc(p*sizeof(int));
  int *to_am=(int *) malloc(p*sizeof(int));
  int *to_dm=(int *) malloc(p*sizeof(int));

  int from_ancestor=0, from_descendant=0, to_ancestor=0, to_descendant=0;

  find_path(&p, &last_from, cur_par, num_par, from_am, &from_ancestor);
  find_path(&p, &last_from, cur_child, num_child, from_dm, &from_descendant);
  find_path(&p, &last_to, cur_par, num_par, to_am, &to_ancestor);
  find_path(&p, &last_to, cur_child, num_child, to_dm, &to_descendant);

  for(k=0;k<3;k++){//k loop
    temp_oper[2]=k+1;
	for(j=0;j<p;j++){//j loop
	  temp_oper[1]=j+1;
	  for(i=0;i<p;i++){//i loop

	    if(j!=i){
	      temp_oper[0]=i+1;
		  recal_ind=1;
          cur_acyc=acyc_all[k*p*p+j*p+i];

          if(cur_acyc==0){//last operation does not create a cycle, no need to update acyc_all
            if((k!=2)&&(last_to==j+1)){//for addition and deletion, which change the neighborhood, need to re-calculate
	          recal_ind=0;
            }
            if((k==2)&&((last_to==i+1)||(last_to==j+1))){//reversal which change the neighborhood
              recal_ind=0;
            }

			  score_type_bic(n, p, Y, tscore, recal_ind, temp_oper, bic_delta_all, cur_BIC, cur_par, num_par, opt_oper, oper_BIC, &min_delta, eps);
			  num_oper=num_oper+1;
          }//end if(acyc_all==0)

          else if(cur_acyc!=-1e+06){//operations that are not in blacklist and the operations that create cycles before last operation, or meanless operations
            temp1=((k==0)&&(to_dm[i]==1)&&(from_am[j]==1));// consider add i->j, there is a path from i to j before the last deletion
            temp2=((k==2)&&(from_am[i]==1)&&(to_dm[j]==1));// reversal i->j i.e. j->i

            if(temp1||temp2){
              check_acyclic(&p, adj_matrix, cur_child, num_child, temp_oper, &acyc_res);
              if(acyc_res==0){//if no cycle, re-calculate bic_delta and update bic_delta
                acyc_all[k*p*p+j*p+i]=0;
				  recal_ind=0;
				  score_type_bic(n, p, Y, tscore, recal_ind, temp_oper, bic_delta_all, cur_BIC, cur_par, num_par, opt_oper, oper_BIC, &min_delta, eps);
				  num_oper=num_oper+1;
              }
              else{
                acyc_all[k*p*p+j*p+i]=99;
              }
            }//end if(temp1||temp2)
          }//end else

		}//end if(j!=i)

	  }//for(i)
	}//for(j)
  }//for(k)

   *MIN_Delta=min_delta;
   *NUM_oper=num_oper;

   free(from_am);
   free(from_dm);
   free(to_am);
   free(to_dm);
 }//end oper_del()



void oper_rev(int *N, int *P, double *Y, int *Tscore, int *adj_matrix, int *last_oper, int *cur_par, int *num_par, int *cur_child, int *num_child, double *cur_BIC, double *EPS, int *acyc_all, double *bic_delta_all, int*opt_oper, double *MIN_Delta, int *NUM_oper){
   // if last operation is reversal, then udpate acyclic information, bic_delta, choosing optimal operation and min_delta
   // if an operation does not create a cycle, then it will only create a cycle after the last operation
   // if an operation creates a cycle before last operation, then
   //  printf("Last Operation is Reversal \n");

  int n=*N, p=*P, temp_oper[3], tscore=*Tscore;
  temp_oper[0]=0;
  temp_oper[1]=0;
  temp_oper[2]=0;
  double oper_BIC[2], delta_temp=0.0, min_delta=1e+06, eps=*EPS;

  int i,j,k,ok;

  int last_from=last_oper[0], last_to=last_oper[1];
  int cur_acyc, acyc_res, recal_ind, num_oper=0;

  int *from_am=(int *) malloc(p*sizeof(int));
  int *from_dm=(int *) malloc(p*sizeof(int));
  int *to_am=(int *) malloc(p*sizeof(int));
  int *to_dm=(int *) malloc(p*sizeof(int));

  int from_ancestor=0, from_descendant=0, to_ancestor=0, to_descendant=0;

  find_path(&p, &last_from, cur_par, num_par, from_am, &from_ancestor);
  find_path(&p, &last_from, cur_child, num_child, from_dm, &from_descendant);
  find_path(&p, &last_to, cur_par, num_par, to_am, &to_ancestor);
  find_path(&p, &last_to, cur_child, num_child, to_dm, &to_descendant);

  for(k=0;k<3;k++){//k loop
    temp_oper[2]=k+1;
	for(j=0;j<p;j++){//j loop
	  temp_oper[1]=j+1;
	  for(i=0;i<p;i++){//i loop
	    if(j!=i){
	      temp_oper[0]=i+1;
          cur_acyc=acyc_all[k*p*p+j*p+i];
	      recal_ind=1;

          if(cur_acyc==0){
	        if((k==0)&&(from_dm[i]==1)&&(to_am[j]==1)){//addition i->j create cycles because of the last reversal
			  acyc_all[k*p*p+j*p+i]=99;
			  recal_ind=1e+06;
			}
	        else if((k==0)&&((last_from==j+1)||(last_to==j+1))){
		      recal_ind=0;
			}

			if((k==1)&&((last_from==j+1)||(last_to==j+1))){
		      recal_ind=0;
			}
			if((k==2)&&((to_am[i]==1)&&(from_dm[j]==1))){
			  if((last_to==i+1)&&(last_from==j+1)){
				check_acyclic(&p, adj_matrix, cur_child, num_child, temp_oper, &acyc_res);
				if(acyc_res==0){
				  recal_ind=0;
				}
				else{
				  acyc_all[k*p*p+j*p+i]=99;
				  recal_ind=1e+06;
				}
			  }
			  else{
				acyc_all[k*p*p+j*p+i]=99;
				recal_ind=1e+06;
			  }
			}
			else if((k==2)&&((last_from==i+1)||(last_from==j+1)||(last_to==i+1)||(last_to==j+1))){
		      recal_ind=0;
			}

			  if(acyc_all[k*p*p+j*p+i]==0){
				  score_type_bic(n, p, Y, tscore, recal_ind, temp_oper, bic_delta_all, cur_BIC, cur_par, num_par, opt_oper, oper_BIC, &min_delta, eps);
				  num_oper=num_oper+1;
			  }
          }//if(cur_acyc==0)
          else if((cur_acyc!=0)&&(cur_acyc!=-1e+06)&&(((from_am[i]==1)&&(to_dm[j]==1))||((to_dm[i]==1)&&(from_am[j]==1)))){
      // edges that are not in blacklist and create a cycle because i and j are on the path of last_from, last_to,re-check acyclic, otherwise no need to check because, it still creates cycles
			check_acyclic(&p, adj_matrix, cur_child, num_child, temp_oper, &acyc_res);
			if(acyc_res==0){
			  acyc_all[k*p*p+j*p+i]=0;
				recal_ind=0;
				score_type_bic(n, p, Y, tscore, recal_ind, temp_oper, bic_delta_all, cur_BIC, cur_par, num_par, opt_oper, oper_BIC, &min_delta, eps);
				num_oper=num_oper+1;
			}
		    else{
			  acyc_all[k*p*p+j*p+i]=99;
			}
		  }//end else--cur_acyc!=0
	    }//end if(j!=i)
	  }//end for(i)
	}//end for(j)
  }//end for(k)

   *MIN_Delta=min_delta;
   *NUM_oper=num_oper;

   free(from_am);
   free(from_dm);
   free(to_am);
   free(to_dm);
 }//end oper_rev()


void comp_delta(double *Delta_temp, double *Min_delta, double *EPS,int *temp_oper, int *opt_oper){

    double temp_min=*Delta_temp-*Min_delta;
    //    double temp=0.0;

    if(fabs(temp_min)<*EPS){
       temp_min=0.0;
       //       temp=(rand()+0.0)/(RAND_MAX+1.0);// for equivalent bic_delta, flip a coin to choose a operation
    }//end if

    if(temp_min<0.0){

       opt_oper[0]=temp_oper[0];
       opt_oper[1]=temp_oper[1];
       opt_oper[2]=temp_oper[2];   // zero movement

       *Min_delta=*Delta_temp;
     }//end if
}

void find_path(int *P, int *Node, int *cur_child, int *num_child, int *res, int *res_length){
  // starting from the node,find ancestor(cur_par) or descentdant(cur_child) of the node
  // return a p-length vector, 0-1, indicating the node is in the set(1--yes, 0--no)
  int p=*P;
  int node=*Node;
  int i;
  int stop=0,cur=0,locus=0;
  int temp, temp_t, cur_node;

  int *path=(int *) malloc(p*sizeof(int));

  for(i=0;i<p;i++){
    path[i]=-1; //path: node indices
    res[i]=0; // have been visited or not
  }

  path[locus]=node;//node itself is its own parent and child
  res[node-1]=1;

  while(stop==0){
    if(cur>locus){// exhaust the whole path
      stop=1;
    }
    else{
     cur_node=path[cur];
     temp=num_child[cur_node-1];
     if(temp>0){
       for(i=0;i<temp;i++){
         temp_t=cur_child[p*i+cur_node-1];
         if(res[temp_t-1]==0){// not in the path yet
           locus=locus+1;
           path[locus]=temp_t;
           res[temp_t-1]=1;
         }
       }//end for
     }//end if(temp>0)
     cur=cur+1;
    }//else (i.e. cur<=locus)
  }//end while

  *res_length=locus+1;

  free(path);

}//end function find_path()



/// main function to calculate bic score for a given operation *////
 void bic_score(int *N, int *P, double *Y, int *Tscore, int *oper, double *cur_BIC,  int *cur_par, int *num_par, double *node_BIC, double *bic_delta, double *EPS, int *info){
  //written on 12-22-2011
  //sample_size: number of samples, i.e., n; num_pred: number of predictors, i.e. p.
  // Y: n by p data matrix
  //oper: 1 by 3 integer vector, (node1, node2, oper): oper=1  ("add"), 2("del"), 3 ("reverse" )
  //cur_BIC: length p vector of current BIC  score for each node
  //cur_par, num_par: current parenets information for each node
  //node_BIC: length 2 vector:  updated  BIC of the from-node and to-node of the cur_oper
  // return: node_BIC and bic_delta: updated BIC scores for the current from and to nodes and the change of overall BIC score (scalar)

  int n=*N, p=*P, from_node=oper[0], to_node=oper[1], cur_oper=oper[2], tscore=*Tscore;
  int i,j,k,ok,ok_two;
  double eps=*EPS, bic_p=1e+6;

  node_BIC[0]=cur_BIC[from_node-1];
  node_BIC[1]=cur_BIC[to_node-1];

  /***addition***/
  if(cur_oper==1){
    bic_add(&n, &p, &from_node, &to_node, Y, &tscore, cur_par,num_par, &bic_p, &ok);
    node_BIC[1]=bic_p; //update to_node BIC
    *bic_delta=node_BIC[1]-cur_BIC[to_node-1];
    *info=ok;
  }/*end if(cur_oper==1)*/


  /***deletion***/
  else if(cur_oper==2){
    bic_del(&n,&p,&from_node, &to_node, Y, &tscore, cur_par,num_par, &bic_p, &ok);
    node_BIC[1]=bic_p; //update to_node BIC
    *bic_delta=node_BIC[1]-cur_BIC[to_node-1];
    *info=ok;
  }/*end deletion*/


  /***reversal== deletion + addition; best way: write addition and deletion as functions and call here***/
  else if(cur_oper==3){
    /*to_node: deletion*/
    bic_del(&n,&p, &from_node, &to_node, Y,&tscore, cur_par,num_par, &bic_p, &ok);
    node_BIC[1]=bic_p; //update to_node BIC

    /*from_node: addition the reseved node*/
    bic_add(&n, &p, &to_node, &from_node, Y,&tscore, cur_par,num_par, &bic_p, &ok_two);
    node_BIC[0]=bic_p; //update from_node BIC
    *info=abs(ok)+abs(ok_two);
///
    *bic_delta=node_BIC[0]+node_BIC[1]-(cur_BIC[to_node-1]+cur_BIC[from_node-1]);
  }/*end reversal*/



  if(fabs(*bic_delta)<eps){
    *bic_delta=0.0;
  }//end if

}/*end function bic_score*/



//sub-function to calculate BIC score for addition operation
//By Ru 12-26-2011
//revised on 11-19-2012, which adds intercept
void bic_add(int *N, int *P, int *Fnode, int *Tnode, double *Y, int *Tscore, int *cur_par, int *num_par, double *bic_res, int *info){
  //if Fnode=-1, then calculate BIC score for current neighborhood, instead for possible addition oepration
  int n=*N, p=*P, from_node=*Fnode, to_node=*Tnode, tscore=*Tscore,termloc=0;
  int i,j,ok,temp,ncol,n_parents, npar_node=num_par[to_node-1];
  double res=0.0;

  double *to_y=(double *) malloc(n*sizeof(double));
  for(i=0; i<n; i++){
    to_y[i]=Y[n*(to_node-1)+i];
  }

  if(from_node>0){
    n_parents=npar_node+1;
    ncol=n_parents+1;
  }
  else{
	n_parents=npar_node;
	ncol=n_parents+1;
  }

  double *x_multi=(double *) malloc(n*ncol*sizeof(double));  //get data corresponding to new_par_index
  for(i=0;i<n;i++){
    x_multi[i]=1;
  }
  termloc=termloc+1;

  if(npar_node>0){
	for(i=0;i<npar_node;i++){
		temp=cur_par[i*p+to_node-1]; ///parent index
		for(j=0;j<n;j++){
			x_multi[termloc*n+j]=Y[(temp-1)*n+j];
		}
		termloc=termloc+1;
	}
  }// end if (npar_node>0)

  if(from_node>0){
	for(j=0;j<n;j++){
	  x_multi[termloc*n+j]=Y[(from_node-1)*n+j];
	}
  }

  bic_cal(&n, &ncol, &n_parents, x_multi, to_y, &tscore, &res, &ok);
  //bic_cal(&n, &n_parents, x_multi, to_y, &tscore, &res, &ok);  // calcualte current BIC for to_y
  *bic_res=res;
  *info=ok;

  free(x_multi);
  free(to_y);

}/*end bic_add function*/

//sub-function to calculate BIC score for deletion operation
//on 12-26-2011
//revised on 11-19-2012, which adds intercept
void bic_del(int *N, int *P, int *Fnode, int *Tnode, double *Y, int *Tscore, int *cur_par, int *num_par, double *bic_res, int *info){

  int n=*N, p=*P, from_node=*Fnode, to_node=*Tnode, tscore=*Tscore;
  int i,j,ok,temp,termloc=0, npar_node=num_par[to_node-1],n_parents=npar_node-1,ncol=n_parents+1;
  double res=0.0;

  double *to_y=(double *) malloc(n*sizeof(double));
  for(i=0; i<n; i++){
    to_y[i]=Y[i+n*(to_node-1)];
  }

  double *x_multi=(double *) malloc(n*(ncol)*sizeof(double));

  for(i=0;i<n;i++){
	x_multi[i]=1;
  }
  termloc=termloc+1;

  if(npar_node>1){
    for(i=0; i<npar_node;i++){
      temp=cur_par[i*p+to_node-1];
      if(temp!=from_node){
        for(j=0;j<n;j++){
          x_multi[termloc*n+j]=Y[(temp-1)*n+j];
        }//end for(j)
        termloc=termloc+1;
      }//end if
    }//end for(i)
  }//end if

  bic_cal(&n, &ncol, &n_parents, x_multi, to_y, &tscore, &res, &ok);
  //bic_cal(&n, &n_parents, x_multi, to_y, &tscore, &res, &ok);
  *bic_res=res;
  *info=ok;

  free(to_y);
  free(x_multi);

}/*end bic_add function*/

/// for a given data matrix Y (n by p) and a given set of indices, get the subset matrix corresponding to the indices


void bic_cal(int *x_nrow, int *x_ncol, int *Npar, double *X,double *y, int *Tscore, double *bic_res, int *info){
	/// variables
	/// x_nrow: integer, number of rows of the design matrix: i.e., number of observations
	/// x_ncol: integer, number of columns of the design matrix: i.e., number of predictors
	/// X: double, design matrix itself
	/// y: double, response vector
	/// return: bic_res: bic score
    //revised on 11-19-2012, which adds intercept
	int i, j, k, ok, nrow=*x_nrow, ncol=*x_ncol, tscore=*Tscore,npar=*Npar;
	double rss, resi;

	double *AAT=(double *) malloc(ncol*ncol*sizeof(double));
	/////////////////////////////////////////////////////
	/// X^T X: consider using existing fucntion in C
	////////////////////////////////////////////////////
	for(i = 0; i < ncol; i++){
		for( j = i; j<ncol; j++){
		    AAT[i*ncol+j]=0.0;   //initialization
			for( k = 0; k < nrow; k++){
				AAT[i*ncol+j] = AAT[i*ncol+j]+ X[k+nrow*i] * X[k+nrow*j];
			}
			AAT[j*ncol+i]=AAT[i*ncol+j];
		}
	}
	/// X'y
	double *Ay=(double *) malloc(ncol*sizeof(double));

	for(j=0; j<ncol; j++){
		Ay[j]=0.0;   ///initialization
		for(i=0; i<nrow;i++){
			Ay[j]=Ay[j]+X[i+nrow*j]*y[i];
		}
	}

	int c2=1;       //number of column of Ay/


	F77_CALL(dposv)("U", &ncol, &c2, AAT, &ncol, Ay, &ncol, &ok);  /// Ay is replace by (X'X)^{-1}X'y, and &1 is wrong
	*info=ok;

	if(ok==0){
		double fitted[nrow];

		rss=0.0;  //inialization
		for(i=0; i<nrow; i++){
			fitted[i]=0.0;
			for(j=0; j<ncol; j++){
				fitted[i]=fitted[i]+X[i+nrow*j]*Ay[j];
			}
			resi=y[i]-fitted[i];
			rss=rss+resi*resi;
		}

		if(rss>0.0){
			if(tscore==1){
				*bic_res=((double) nrow)*log(rss/((double) nrow))+((double) npar)*log((double) nrow);
			}
			else if(tscore==2){
				*bic_res=((double) nrow)*log(rss/((double) nrow));
			}
		}
		else{// if rss=0
			*bic_res=-(1e+6);
		}
	}//end if(ok==0)
	else{// if error
		*bic_res=1e+6;
	}

	free(AAT);
	free(Ay);
}/*end bic_score function*/

void bic_par(int *N, int *P, double *Y, int *Tscore, int *cur_par, int *num_par, double *cur_BIC){
	int p=*P, n=*N, i, tscore=*Tscore, ok=1,fnode,tnode;
	double bic_res;

	for(i=0;i<p;i++){
		tnode=i+1;
 		fnode=-1;
        bic_add(&n, &p, &fnode, &tnode, Y, &tscore, cur_par, num_par, &bic_res, &ok);
		cur_BIC[i]=bic_res;
    }
}


void check_acyclic(int *P, int *adj_matrix, int *cur_child_copy, int *num_child_copy, int *oper, int *res){
//written on 01-08-2012
//P: number of predictor
//adj_matrix: is p*p length vector, recording the adjacency matrix
//oper: 3-length vector, from_node, to_node, operation
//cur_child, num_child
//return res: indicate whether oper is acyclic or not. (0--acyclic, 1--cyclic)

  int p=*P;
  int from_node=oper[0];
  int to_node=oper[1];
  int cur_oper=oper[2];
  int acyclic=0, err=0, i, temp, locus;

  int edge_exist=adj_matrix[(to_node-1)*p+from_node-1];
  int tedge_exist=adj_matrix[(from_node-1)*p+to_node-1];

  if((cur_oper==1)&&((edge_exist==1)||(tedge_exist==1))){
    //    printf("%s \n","The edge is already in the graph, can not do addition."); //The edge is already in the graph, can not do addition
    err=1;
    acyclic=1;
  }

  if((cur_oper==3)&&((edge_exist==0)||(tedge_exist==1))){
    //    printf("%s \n", "The edge is not in the graph, or the reversed edge is already in the graph, can not do reversal");
    err=1;
    acyclic=1;
  }

  if((cur_oper==2)&&(edge_exist==0)){
    err=1;
    acyclic=1;
  }

  if(err==0){
    if(cur_oper==1){
      acyclic_add(&p, &from_node, &to_node, cur_child_copy, num_child_copy, &acyclic);
    }

    else if(cur_oper==3){

      for(i=0;i<num_child_copy[from_node-1];i++){
        temp=cur_child_copy[i*p+from_node-1];
        if(temp==to_node){
          cur_child_copy[i*p+from_node-1]=-1;
          locus=i;
          break;
        }
      }


      acyclic_add(&p, &to_node, &from_node,cur_child_copy, num_child_copy, &acyclic);
      cur_child_copy[locus*p+from_node-1]=to_node;
      //  del_parchild(&p,&from_node,&to_node,cur_child_copy, num_child_copy);
      //add_parchild(&p,&from_node,&to_node,cur_child_copy, num_child_copy);
    }//end ifelse

  }//end if(err==0)

  *res=acyclic;

}//end function check_acyclic



void del_parchild(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child){
  //written on 01-17-2012
  //to update cur_child and num_child when we do deletion opeation
  // P: p
  // Fnode: from node; Tnode: to node

  int p=*P;
  int from_node=*Fnode;
  int to_node=*Tnode;
  int i, temp, temp_t=0;
  int child_num=num_child[from_node-1]; //current number of children of from_node

  for(i=0;i<child_num;i++){//update cur_child
    temp=cur_child[i*p+from_node-1];

    if(temp==to_node){
      cur_child[i*p+from_node-1]=cur_child[(child_num-1)*p+from_node-1];
      cur_child[(child_num-1)*p+from_node-1]=0;
      break;
    }
  }//end for

  num_child[from_node-1]=child_num-1;
}



void add_parchild(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child){
  //written on 01-17-2012
  //to update cur_child and num_child when we do addition opeation
  // P: p
  // Fnode: from node; Tnode: to node

  int p=*P;
  int from_node=*Fnode;
  int to_node=*Tnode;
  int child_num=num_child[from_node-1]; //current number of children of from_node

  num_child[from_node-1]=child_num+1;
  cur_child[child_num*p+from_node-1]=to_node;
}

void acyclic_add(int *P, int *Fnode, int *Tnode, int *cur_child, int *num_child, int *acyclic){
//updated on 01-16-2012
//written on 01-08-2012
// Fnode: from_node, Tnode: to_node
//cur_child:
//num_child
//return acyclic, 0--acyclic, 1--cyclic

  int from_node=*Fnode;
  int to_node=*Tnode;
  int res=0,i=0;
  int p=*P;

  int *path=(int *) malloc(p*sizeof(int));
  int *exist_v=(int *) malloc(p*sizeof(int));

  for(i=0;i<p;i++){
    path[i]=-1; //path
    exist_v[i]=0; // have been visited or not
  }

  path[0]=to_node;
  exist_v[to_node-1]=1;

  int stop=0, temp, temp_t, cur_node;

  int cur=1; //current position in the path
  int locus=1; //last position in the path

  while(stop==0){
    if(cur>locus){// exhaust the whole path
      stop=1;
    }
    else{
       cur_node=path[cur-1];
       if(cur_node==from_node){//find a cycle
         stop=1;
         res=1;
       }
       else{
         temp=num_child[cur_node-1];
         if(temp>0){
           for(i=0;i<temp;i++){
             temp_t=cur_child[p*i+cur_node-1];
             if(temp_t==-1){
               continue;
             }
             if(exist_v[temp_t-1]==0){// not in the path yet
               path[locus]=temp_t;
               locus=locus+1;
               exist_v[temp_t-1]=1;
             }
           }//end for
         }//end if
        cur=cur+1;
      }//else (i.e. cur_node!=from_node)
    }//end else(i.e. cur<=locus)
  }//end while

  *acyclic=res;

  free(path);
  free(exist_v);
}//end acyclic_add()

//      sampling(&num_oper,r_index,&r_temp,e_delta,&r_res);
void sampling(int *NUM, int *index, int *n_sample, double *prob, int *res){
  //written on 02-29-2012
  //correponding to sample() in R
  //replacement is not set as a parameter, then viewed as "TRUE"
  //NUM: the length of index,
  //index: index sampled from, for example, 1:n
  //n_sample: the number of samples

  int n=*NUM, n_s=*n_sample;
  double temp;
  int i,j,flag;
  double *cum_prob=(double *) malloc((n+1)*sizeof(double));
  cum_prob[0]=0.0;


  if(prob[0]<0.0){//uniform sampling
    double p=(double ) 1/ (double) n;
    for(j=1;j<=n;j++){
      cum_prob[j]=j*p;
    }
  }

  else{
    for(j=1;j<=n;j++){
      cum_prob[j]=cum_prob[j-1]+prob[j-1];
    }
  }

  RANDIN;
  for(i=0;i<n_s;i++){
    temp=UNIF;
    j=0;
    flag=0;

    while(flag==0){

      if((temp<=cum_prob[j+1])&&(temp>=cum_prob[j])){
        res[i]=index[j];
        flag=1;
      }//end if
      j=j+1;
    }//end while
  }//end for
  RANDOUT;

  free(cum_prob);
}

void SampleNoReplace(int k, int n, int *y, int *x)
{
    int i, j;
	RANDIN;
    for (i = 0; i < n; i++)
		x[i] = i;
    for (i = 0; i < k; i++) {
		j = n * unif_rand();
		y[i] = x[j];
		x[j] = x[--n];
    }
	RANDOUT;
}

void randforest(int *P, int *NUM_oper, double *MIN_delta, double *RAND_step_length, int *acyc_all, double *bic_delta_all, double *EPS, int *cur_step, int *LL, double *MIN_valve, int *e_oper){
  //written on 03-02-2012 to implement random forest idea
  //RAND_step_length:temperature
  //LL:p*(p-1)*3*step_max
  //return e_oper: selected operation

//  if(*MIN_delta>*MIN_valve){
//    e_oper[0]=0;
//    e_oper[1]=0;
//    e_oper[2]=0;
//  }

//  else{
    int p=*P, num_oper=*NUM_oper+1,r_i,r_j,r_k, s=*cur_step, temp_LL=*LL;
    double min_delta=*MIN_delta, rand_step_length=*RAND_step_length, eps=*EPS;

    int *rand_oper=(int *) malloc(num_oper*3*sizeof(int));//+1 for no-movement (0,0,0)
    double *e_delta=(double *) malloc(num_oper*sizeof(double));
    int *r_index=(int *) malloc(num_oper*sizeof(int));
    int r_temp=1, r_res, max_index=0;
    double e_sum=0.0, diff=0.0;

    rand_oper[0]=0;
    rand_oper[num_oper]=0;
    rand_oper[2*num_oper]=0;

    e_delta[0]=exp(min_delta*rand_step_length);//selection probility for no-movement
    e_sum=e_sum+e_delta[0];

    int r_count=1;

      for(r_k=0;r_k<3;r_k++){//k loop
		for(r_j=0;r_j<p;r_j++){//j loop,to node
         for(r_i=0;r_i<p;r_i++){//i loop, from node

	       if(r_i!=r_j){
            if(acyc_all[r_k*p*p+r_j*p+r_i]==0){
              rand_oper[r_count]=r_i+1;
              rand_oper[num_oper+r_count]=r_j+1;
              rand_oper[num_oper*2+r_count]=r_k+1;

              diff=bic_delta_all[r_k*p*p+r_j*p+r_i]-min_delta;
              if(fabs(diff)<=eps){
                diff=0.0;
              }//end if(fabs(diff)<=eps)
	          e_delta[r_count]=exp(0.0-diff*rand_step_length);

/*			  if(diff==0.0){
				printf("which one is smallest %d\n",r_count);
			  }*/

              e_sum=e_sum+e_delta[r_count];
              r_count=r_count+1;

            }//end if(acyc_all[r_k*p*p+r_j*p+r_i]==0)
          }//end if(r_i!=r_j)

         }//end for r_i
       }//end for r_j
      }//end for r_k

      for(r_i=0;r_i<num_oper;r_i++){
        e_delta[r_i]=e_delta[r_i]/e_sum;
        r_index[r_i]=r_i+1;
      }

      sampling(&num_oper,r_index,&r_temp,e_delta,&r_res);

	  e_oper[0]=rand_oper[r_res-1];
      e_oper[1]=rand_oper[num_oper+r_res-1];
      e_oper[2]=rand_oper[num_oper*2+r_res-1];
	  *MIN_delta=bic_delta_all[(e_oper[2]-1)*p*p+(e_oper[1]-1)*p+e_oper[0]-1];

/*	  printf("which number is samples, supported to be smallest%d\n",r_res);
	  printf("e_delta value is %f\n",e_delta[r_res-1]);
	  printf("BIC delta value is %f--",bic_delta_all[(e_oper[2]-1)*p*p+(e_oper[1]-1)*p+e_oper[0]-1]);
	  printf("BIC delta value is %f\n",min_delta);
      printf("------------------\n");*/

      free(rand_oper);
      free(e_delta);
      free(r_index);

//  }//end else

}


void whitelist_blacklist(int *P, int *STEP_max, int *adj_matrix, int *movement, int *whitelist, int *blacklist, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child, int *Step_rec){
		// acyc_all: 0:OK,99:acyclic,-99:meanless,-1e+06:blacklist
		// check for allowable operations according to whitelist and blacklist

		int p=*P, step_max=*STEP_max, i,j,step_rec=0,itemp,jtemp,nwhite=0;

		for(j=0;j<p;j++){
			for(i=0;i<p;i++){
				acyc_all[j*p+i]=0;  ///allowable -- 0; not allowable: -(1e+06)
				acyc_all[p*p+j*p+i]=-99;
				acyc_all[2*p*p+j*p+i]=-99;

				if(whitelist[j*p+i]==1){
					itemp=i+1;
					jtemp=j+1;
					add_parchild(&p, &itemp, &jtemp, cur_child, num_child);
					add_parchild(&p, &jtemp, &itemp, cur_par, num_par);
					adj_matrix[j*p+i]=1;

					acyc_all[j*p+i]=-(1e+06);
					acyc_all[p*p+j*p+i]=-(1e+06);
					acyc_all[2*p*p+j*p+i]=-(1e+06);

					acyc_all[i*p+j]=-(1e+06);
					acyc_all[p*p+i*p+j]=-(1e+06);
					acyc_all[2*p*p+i*p+j]=-(1e+06);

					movement[step_rec]=i+1;
					movement[step_max+step_rec]=j+1;
					movement[2*step_max+step_rec]=1;
					step_rec=step_rec+1;
					nwhite=nwhite+1;
				}

				if(whitelist[i*p+j]==1){
					acyc_all[j*p+i]=-(1e+06);
					acyc_all[p*p+j*p+i]=-(1e+06);
					acyc_all[2*p*p+j*p+i]=-(1e+06);
				}

				if(blacklist[j*p+i]==1){
					acyc_all[j*p+i]=-(1e+06);
					acyc_all[p*p+j*p+i]=-(1e+06);
					acyc_all[2*p*p+j*p+i]=-(1e+06);
				}

				if(blacklist[i*p+j]==1){
					acyc_all[2*p*p+j*p+i]=-(1e+06);//reverse of the reverse edge is forbidden
				}

			}//end for
		}//end for

		int cand_oper[3], acyc_temp=0;

		if(nwhite>0){//check for additiona nonallowable additions  after whilelist and blacklist
			for(j=0;j<p;j++){
				for(i=0;i<p;i++){
					if((j!=i)&&(acyc_all[j*p+i]==0)){//check acylc after adding edges from whitelist
						cand_oper[0]=i+1;
						cand_oper[1]=j+1;
						cand_oper[2]=1;
						check_acyclic(&p, adj_matrix, cur_child, num_child, cand_oper, &acyc_temp);
						if(acyc_temp==1){// if acyclic, then set the operation nonallowable
							acyc_all[j*p+i]=-(1e+06);
							acyc_all[p*p+j*p+i]=-(1e+06);
							acyc_all[2*p*p+j*p+i]=-(1e+06);
						}//end if
					}//end if
				}//end for
			}//end for
		}//end if

		*Step_rec=step_rec;
}//end function whitelist_blacklist

void adj_ini_fn(int *P, int *adj_matrix, int *ini_adj_matrix, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child){
		int p=*P, i, j, itemp, jtemp, acyc_temp=0, cand_oper[3];

		for(j=0;j<p;j++){
			for(i=0;i<p;i++){
				adj_matrix[j*p+i]=ini_adj_matrix[j*p+i];
				if(ini_adj_matrix[j*p+i]==1){
					itemp=i+1;
					jtemp=j+1;
					add_parchild(&p, &itemp, &jtemp, cur_child, num_child);
					add_parchild(&p, &jtemp, &itemp, cur_par, num_par);
					acyc_all[j*p+i]=99;
					acyc_all[p*p+j*p+i]=0;
				}//end if
			}
		}

		for(j=0;j<p;j++){
			for(i=0;i<p;i++){
				if((j!=i)&&(acyc_all[j*p+i]==0)){
					cand_oper[0]=i+1;
					cand_oper[1]=j+1;
					cand_oper[2]=1;
					check_acyclic(&p, adj_matrix, cur_child, num_child, cand_oper, &acyc_temp);
					if(acyc_temp==1){// if acyclic, then set the operation nonallowable
						acyc_all[j*p+i]=99;
					}//end if
				}//end if

				if(ini_adj_matrix[j*p+i]==1){
					cand_oper[0]=i+1;
					cand_oper[1]=j+1;
					cand_oper[2]=3;
					check_acyclic(&p, adj_matrix, cur_child, num_child, cand_oper, &acyc_temp);
					if(acyc_temp==1){// if acyclic, then set the operation nonallowable
						acyc_all[2*p*p+j*p+i]=99;
					}//end if
				}
			}//end for
		}//end for

}//end function adj_ini_fn


/// calculate mean of a vector
void mean(int *N, double *x, double *res){
  // length: length of the vector
  // x: vector
  //return: res: mean(x)
  int n=*N;
  int i;
  double sum=0.0;

  for(i=0;i<n;i++){
    sum=sum+x[i];
  }

  *res=sum/(n+0.0);

}

void data_stdize(int *N, int *P, double *x, double *x_std){

  int n=*N,p=*P,i,j;
  double x_sum, x_mean, var, sd;


    for(j=0;j<p;j++){
      x_sum=0.0,var=0.0;

      for(i=0;i<n;i++){
        x_sum=x_sum+x[j*n+i];
      }
      x_mean=x_sum/(double) n;

      for(i=0;i<n;i++){
	var=var+pow((x[j*n+i]-x_mean),2);
      }
      sd=sqrt(var/((double) (n-1)));

      for(i=0;i<n;i++){
        x_std[j*n+i]=(x[j*n+i]-x_mean)/sd;
      }//end for(i)
    }//end for(j)


}//end data_stdize

/* Unequal probability sampling; without-replacement case */

static void ProbSampleNoReplace(int n, double *p, int *perm,
								int nans, int *ans)
{
    double rT, mass, totalmass;
    int i, j, k, n1;

    /* Record element identities */
    for (i = 0; i < n; i++)
		perm[i] = i + 1;

    /* Sort probabilities into descending order */
    /* Order element identities in parallel */
    revsort(p, perm, n);

    /* Compute the sample */
    totalmass = 1;
    for (i = 0, n1 = n-1; i < nans; i++, n1--) {
		rT = totalmass * unif_rand();
		mass = 0;
		for (j = 0; j < n1; j++) {
			mass += p[j];
			if (rT <= mass)
				break;
		}
		ans[i] = perm[j];
		totalmass -= p[j];
		for(k = j; k < n1; k++) {
			p[k] = p[k + 1];
			perm[k] = perm[k + 1];
		}
    }
}


///////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////
//// DAGBagging: gshd

struct gshdm
{
	int fnode;
	int tnode;
	int do_oper;
	double freq;
};


//// declare functions for aggregation
void score_gshd(int *P, int *boot_adj, double *Alpha, int *NB, int *STEP_MAX, int *PRINT, int *blacklist, int *whitelist, double *Threshold, int *movement, int *adj_matrix, int *final_step, int *c_size);
int compare_double(double a, double b);
int compare_oper(const struct gshdm *a, const struct gshdm *b);
//int compare_oper(const void *a, const void *b);

void GSHDoper_last(int *P, int *last_oper, int *cand_oper, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child, int *res, int *c_size);



void score_gshd(int *P, int *boot_adj, double *Alpha, int *NB, int *STEP_MAX, int *PRINT, int *blacklist, int *whitelist, double *Threshold, int *movement, int *adj_matrix, int *final_step, int *c_size){

	int i, j=0, k=0, nb=*NB, s, p=*P, cp=p*(p-1), step_max=*STEP_MAX, pr=*PRINT, temp_gshd1=0, temp_gshd2=0, csize=*c_size;
	double threshold=*Threshold, alpha=*Alpha;

	//initiation of cur_par, cur_child, num_par, num_child, bic_delta_all and acyc_all
	int *cur_par=(int *)malloc(p*p*sizeof(int));
	int *cur_child=(int *)malloc(p*p*sizeof(int));
	int *num_par=(int *)malloc(p*sizeof(int));
	int *num_child=(int *)malloc(p*sizeof(int));
	int *acyc_all=(int *)malloc(p*p*sizeof(int));
	struct gshdm *gshd_all=(struct gshdm *)malloc(cp*sizeof(struct gshdm));

	int count=0, step_rec=0 ,nwhite=0, itemp, jtemp;
	//initiation of cur_par, cur_child, acyc_all
    for(j=0;j<p;j++){
		num_par[j]=0;
		num_child[j]=0;
		for(i=0;i<p;i++){
			cur_par[j*p+i]=0;
			cur_child[j*p+i]=0;
			acyc_all[j*p+i]=0;
		}
    }

    for(j=0;j<p;j++){
		for(i=0;i<p;i++){
			if(whitelist[j*p+i]==1){//edges in whitelist
				itemp=i+1;
				jtemp=j+1;
				add_parchild(&p, &itemp, &jtemp, cur_child, num_child);
				add_parchild(&p, &jtemp, &itemp, cur_par, num_par);
				adj_matrix[j*p+i]=1;

				acyc_all[j*p+i]=-(1e+06);
				acyc_all[i*p+j]=-(1e+06);
				movement[step_rec]=i+1;
				movement[step_max+step_rec]=j+1;
				movement[2*step_max+step_rec]=1;
				step_rec=step_rec+1;
				nwhite=nwhite+1;
			}//end if
			if(whitelist[i*p+j]==1){//edges in whitelist
				acyc_all[j*p+i]=-(1e+06);
			}
			if(blacklist[j*p+i]==1){
				acyc_all[j*p+i]=-(1e+06);
			}

		}//end for i
    }//end for j


    int acyc_temp=0, cand_oper[3];
    for(j=0;j<p;j++){
		for(i=0;i<p;i++){

			if(j!=i){
				if((nwhite>0)&&(acyc_all[j*p+i]==0)){//check acylc after adding edges from whitelist
					cand_oper[0]=i+1;
					cand_oper[1]=j+1;
					cand_oper[2]=1;
					check_acyclic(&p, adj_matrix, cur_child, num_child, cand_oper, &acyc_temp);
					if(acyc_temp==1){
						acyc_all[j*p+i]=-(1e+06);
					}
				}

				gshd_all[count].fnode=i+1;
				gshd_all[count].tnode=j+1;
				gshd_all[count].do_oper=1;

				temp_gshd1=0;
				temp_gshd2=0;
				for(k=0;k<nb;k++){
					temp_gshd1=temp_gshd1+boot_adj[k*p*p+j*p+i];
					temp_gshd2=temp_gshd2+boot_adj[k*p*p+i*p+j];
				}
				gshd_all[count].freq=((double) temp_gshd1)/((double) nb)+(1.0-alpha/2.0)*((double) temp_gshd2)/((double) nb);
				count=count+1;
			}//end if
		}//end for i
    }//end for j


    qsort(gshd_all,cp,sizeof(struct gshdm),compare_oper);

//	for(i=0;i<cp;i++){
//	  printf("%d--",gshd_all[i].fnode);
//	  printf("%d--",gshd_all[i].tnode);
//	  printf("%f\n",gshd_all[i].freq);
//	}

	//start DAG
    //to prevent error

    int cur_fnode, cur_tnode, last_oper[3], opt_oper=0;
	double cur_freq;
    last_oper[2]=-1;

	for(s=0;s<step_max;s++){
		if(pr==1){
///6/23			printf("step -- %d \n",s);
		}//end if

		cur_fnode=gshd_all[s].fnode;
		cur_tnode=gshd_all[s].tnode;

		cand_oper[0]=cur_fnode;
		cand_oper[1]=cur_tnode;
		cand_oper[2]=1;

		cur_freq=gshd_all[s].freq;

		acyc_temp=acyc_all[(cand_oper[1]-1)*p+(cand_oper[0]-1)];

		if(acyc_temp==0){
			if(last_oper[2]>0){
				GSHDoper_last(&p, last_oper, cand_oper, acyc_all, cur_par, num_par, cur_child, num_child, &acyc_temp, &csize);
				//        acyc_all[(cand_oper[1]-1)*p+cand_oper[0]-1]=acyc_temp; 09-19-2012
			}

			if((acyc_temp==0)&&(cur_freq>threshold)){
				opt_oper=1;
			}

			if(opt_oper==1){//current optimal is addition
				add_parchild(&p, &cur_fnode, &cur_tnode, cur_child, num_child);
				add_parchild(&p, &cur_tnode, &cur_fnode, cur_par, num_par);
				adj_matrix[(cur_tnode-1)*p+cur_fnode-1]=1;

				acyc_all[(cur_tnode-1)*p+(cur_fnode-1)]=-99;//cannot add the same edge again
				acyc_all[(cur_fnode-1)*p+(cur_tnode-1)]=-99;//cannot add the reverse edge again

				//record optimal movement
				movement[step_rec]=cur_fnode;
				movement[step_max+step_rec]=cur_tnode;
				movement[2*step_max+step_rec]=1;
				//update last movement

				last_oper[0]=cur_fnode;
				last_oper[1]=cur_tnode;
				last_oper[2]=1;
				step_rec=step_rec+1;
				opt_oper=0;
			}

		}//end if(acyc_temp==0)
		//////////////////////////////////////////////////////////////
		if((cur_freq<=threshold)){
///6/23			printf("DAG established, no more improvement");
			break;
		}

	}//end for loop
	*final_step=step_rec;

	free(gshd_all);
	free(acyc_all);
	free(cur_child);
	free(cur_par);
	free(num_child);
	free(num_par);

}


void GSHDoper_last(int *P, int *last_oper, int *cand_oper, int *acyc_all, int *cur_par, int *num_par, int *cur_child, int *num_child, int *res, int *c_size){
	//written on 09-19-2012
	//if last operation is addition, update acyclic information, bic_delta, choosing optimal operation, and minimum delta
	//when last_operation is addition, then if an operation creates a cycle before the last operation, it will create a cycle still.
	int p=*P,cur_acyc,fnode=cand_oper[0]-1, tnode=cand_oper[1]-1, csize=*c_size;
	int i,j;
	int last_from=last_oper[0],last_to=last_oper[1];

	int *from_am=(int *) malloc(p*sizeof(int));
	int *from_dm=(int *) malloc(p*sizeof(int));
	int *to_am=(int *) malloc(p*sizeof(int));
	int *to_dm=(int *) malloc(p*sizeof(int));

	int from_ancestor=0, from_descendant=0, to_ancestor=0, to_descendant=0;

	find_path(&p, &last_from, cur_par, num_par, from_am, &from_ancestor);// find last_from's ancestor
	find_path(&p, &last_from, cur_child, num_child, from_dm, &from_descendant);// find last_from's descendant
	find_path(&p, &last_to, cur_par, num_par, to_am, &to_ancestor);//find last_to's ancestor
	find_path(&p, &last_to, cur_child, num_child, to_dm, &to_descendant);//find last_to's descendant

	for(j=0;j<p;j++){//j loop,to node
		for(i=0;i<p;i++){//i loop, from node
			if(i!=j){
				cur_acyc=acyc_all[j*p+i];

				if(cur_acyc==0){//last time no cycle
					if((from_am[j]==1)&&(to_dm[i]==1)){//create cycle
						acyc_all[j*p+i]=99;
						csize=csize+1;
					}
				}//
			}
		}
	}
	*res=acyc_all[tnode*p+fnode];
	*c_size=csize;
	free(from_am);
	free(from_dm);
	free(to_am);
	free(to_dm);

}//end oper_last()

int compare_double(double a, double b)
{
	double temp = a - b;
	if (temp>=0.0){
		return -1;
	}
	else{
		return 1;
	}
}

int compare_oper(const struct gshdm *a, const struct gshdm *b){
	return compare_double(a->freq,b->freq);
}



//int compare_oper(const void *a, const void *b){
//	return compare_double(*(struct gshdm*)a->freq,*(struct gshdm*)b->freq);
//}
