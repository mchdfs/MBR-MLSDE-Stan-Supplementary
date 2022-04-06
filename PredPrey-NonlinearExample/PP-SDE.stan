functions {
real[] P_mat_to_array(matrix Pmat){
	real Pvec[3];
	Pvec[1] = Pmat[1,1];
	Pvec[2] = Pmat[1,2];
	Pvec[3] = Pmat[2,2];
	return Pvec;
}


matrix P_array_to_mat(real[] Pvec){
	matrix[2,2] Pmat;
	Pmat[1,1] = Pvec[1];
	Pmat[1,2] = Pvec[2];
	Pmat[2,1] = Pvec[2];
	Pmat[2,2] = Pvec[3];
	return Pmat;
}
real[] PPSDE(real t,
real[] yp,
real[] params,
real[] x,
int[] x_int) {
  real dyPdt[5];
  real a;
  real b;
  matrix[2,2] Jacob;
  

  matrix[2,2] sigmaP;
  real agrowth;
  real aconsumed;
  real bgrowth;
  real bdeath;
  real Pvec[3];

  agrowth=params[1];
  aconsumed=params[2];
  bgrowth=params[3];
  bdeath=params[4];
  Pvec[1]=params[5];
  Pvec[2]=params[6];
  Pvec[3]=params[7];

  sigmaP = P_array_to_mat(Pvec);
  a = yp[1];
  b = yp[2];
Jacob[1,1]=agrowth - aconsumed * b;Jacob[1,2]=-(a * aconsumed);Jacob[2,1]=b * bgrowth;Jacob[2,2]=a * bgrowth - bdeath;

  dyPdt[1] = agrowth * a - aconsumed * a * b;
  dyPdt[2] = bgrowth * a * b - bdeath * b;
  dyPdt[3:5] = P_mat_to_array(Jacob*P_array_to_mat(yp[3:5])+P_array_to_mat(yp[3:5])*Jacob'+sigmaP*sigmaP');
  return dyPdt;
}
real ACDEKF(real[,] y, matrix Lambda,
     		   real[] y0,real t0,real[,] ts,
		   real[] params, matrix sigmaY, 
		   real[] x,int[] x_int, int T){
  real y_pred[T,1,5];
  real y_up[T,5];
  matrix[2,1] v[T];
  matrix[2,2] V[T];
  real ll=0;
  matrix[2,2] K;

  y_pred[1,,] = integrate_ode_rk45(PPSDE, y0, t0, ts[1,], params, x, x_int); //first time point prediction step
  
  v[1] = to_matrix(y[1,],2,1)-to_matrix(y_pred[1,,1:2])'; // prediction error
  V[1] = Lambda*P_array_to_mat(y_pred[1,1,3:5])*Lambda'+sigmaY;
  
  K = P_array_to_mat(y_pred[1,1,3:5])*Lambda'*inverse(V[1]);
  y_up[1,1:2] = to_array_1d((to_matrix(y_pred[1,,1:2])'+K*v[1])');
  y_up[1,3:5] = P_mat_to_array(P_array_to_mat(y_pred[1,1,3:5])-K*Lambda*P_array_to_mat(y_pred[1,1,3:5]));
  
  ll+= -(log(determinant(V[1]))+v[1]'*inverse(V[1])*v[1])[1,1];

  for(t in 2:T){
    y_pred[t,,] = integrate_ode_rk45(PPSDE, y_up[t-1,], ts[t-1,1], ts[t,], params, x, x_int);
	
	//stan requires that the first time index is of type real and the second of type real[] (intended for an entire sequence)
	//take in starting point y as real[], returns real[ , ] (again, intended for the entire sequence of y's)
    
    v[t] = to_matrix(y[t,],2,1)-to_matrix(y_pred[t,,1:2])'; // prediction error
    V[t] = Lambda*P_array_to_mat(y_pred[t,1,3:5])*Lambda'+sigmaY;
  
    K = P_array_to_mat(y_pred[t,1,3:5])*Lambda'*inverse(V[t]);
    y_up[t,1:2] = to_array_1d((to_matrix(y_pred[t,,1:2])'+K*v[t])');
    y_up[t,3:5] = P_mat_to_array(P_array_to_mat(y_pred[t,1,3:5])-K*Lambda*P_array_to_mat(y_pred[t,1,3:5]));
  
  ll+= -(log(determinant(V[t]))+v[t]'*inverse(V[t])*v[t])[1,1];
  }
  return ll;
 }

}


data {
  int<lower=1> T;
  real y[T,2];
  real t0;
  real ts[T,1];
  matrix[2,2] Lambda;
  real pstart[3];
  real ystart[2];
}

transformed data {
  real x[0];
  int x_int[0];
}
parameters { 
  real<lower=0> agrowth;
  real<lower=0> aconsumed;
  real<lower=0> bgrowth;
  real<lower=0> bdeath;

  real log_sigmaP[2];
  real rho;
  real log_sigmaY[2];
}
transformed parameters{
  real<lower=0> sigmaP[2]=exp(log_sigmaP);
  real<lower=0> sigmaY[2]=exp(log_sigmaY);
}
model {
  real ll;
  real params[7];
  real y0[5];
  matrix[2,2] SigmaY;
  
//  sigmaY ~ cauchy(0,2.5);
//  sigmaP ~ cauchy(0,2.5);
  
  params[1]= agrowth;
  params[2]= aconsumed;
  params[3]= bgrowth;
  params[4]= bdeath;
  
  params[5]=sigmaP[1];
  params[6]=rho;
  params[7]= sigmaP[2];
 
  SigmaY[1,1]=sigmaY[1]^2;
  SigmaY[2,2]=sigmaY[2]^2;
  SigmaY[1,2]=0;
  SigmaY[2,1]=0;
  
  y0[1:2]=ystart;
  y0[3:5]=pstart;
  ll=ACDEKF(y, Lambda, y0, t0, ts,params, SigmaY, x, x_int, T);
  target+=ll;
  target+= log_sigmaP[1] + log_sigmaP[2] + log_sigmaY[1] + log_sigmaY[2]; 
}





