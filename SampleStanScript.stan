// Start of Function Block
functions {

// P_mat_to_array() puts unique elements in a symmetric matrix into a real-number array
real[] P_mat_to_array(matrix Pmat){
	real Pvec[3];
	Pvec[1] = Pmat[1,1];
	Pvec[2] = Pmat[1,2];
	Pvec[3] = Pmat[2,2];
	return Pvec;
}

// P_mat_to_array() puts unique elements in a symmetric matrix into a real-number array
matrix P_array_to_mat(real[] Pvec){
	matrix[2,2] Pmat;
	Pmat[1,1] = Pvec[1];
	Pmat[1,2] = Pvec[2];
	Pmat[2,1] = Pvec[2];
	Pmat[2,2] = Pvec[3];
	return Pmat;
}

real[] DyadicInteraction(real t, // time
real[] eta, // states (this includes both the latent states of the system and the ones for the prediction error covariance matrices)
real[] params, // parameters
real[] x, // real number covariates, if any
int[] x_int) { // integer covariates, if any
  real dFdt[5];
  matrix[2,2] Jacob;
  matrix[2,2] sigma_p;

  real a1;
  real a2;
  real b1;
  real b2;
  real mu1;
  real mu2;

  a1=params[1];
  a2=params[2];
  b1=params[3];
  b2=params[4];
  mu1=params[5];
  mu2=params[6];

  sigma_p[1,1] = params[7];
  sigma_p[2,2] = params[9];
  sigma_p[2,1] = params[8];
  sigma_p[1,2] = 0;

  Jacob[1,1]=-(a1 + a2);Jacob[1,2]=a2;Jacob[2,1]=b2;Jacob[2,2]=-(b1 + b2);

  dFdt[1] = a1 * (mu1 - eta[1]) + a2 * (eta[2] - eta[1]);
  dFdt[2] = b1 * (mu2 - eta[2]) + b2 * (eta[1] - eta[2]);
  dFdt[3:5] = P_mat_to_array(Jacob*P_array_to_mat(eta[3:5])+
               P_array_to_mat(eta[3:5])*Jacob'+sigma_p*sigma_p');
  return dFdt;
}

real CDEKF(real[,] y, //individual time-series y_i
                   matrix Lambda, // measurement loading matrix
                   matrix sigma_y, // measurement error variance matrix
                   real[] y0, // initial values for the system
                   real t0, // initial time index
                   real[] params, // other parameters in the model
                   real[,] ts, // array of observed time indices 
                   real[] x, // real number covariates, if any
                   int[] x_int, // integer covariates, if any
                   int T, // an interger for total number of observations
     		       real rel_tol, real abs_tol, real max_step // numerical solver arguments
                    ){
    real y_pred[T,1,5];
    real y_up[T,5];
    matrix[2,1] v[T];
    matrix[2,2] V[T];
    real ll=0;
    matrix[2,2] K;
    matrix[2,2] sigma_p;
    
    sigma_p[1,1] = params[7];
    sigma_p[2,2] = params[9];
    sigma_p[2,1] = params[8];
    sigma_p[1,2] = 0;
 
// first time point prediction step
// This function requires that the first time index is of type real and the second of type real[] (as it is intended for solving an entire sequence).
// It takes in starting point y_up[t-1,] as real[], returns real[ , ] (again, intended for the entire sequence of y's)
    y_pred[1,,] = integrate_ode_rk45(DyadicInteraction, y0, t0, ts[1,], params, x, x_int,rel_tol,abs_tol,max_step); 

// prediction error
    v[1] = to_matrix(y[1,],2,1)-to_matrix(y_pred[1,,1:2])'; 
    V[1] = Lambda*P_array_to_mat(y_pred[1,1,3:5])*Lambda'+sigma_y;

    K = P_array_to_mat(y_pred[1,1,3:5])*Lambda'*inverse(V[1]);

// first time point update step
    y_up[1,1:2] = to_array_1d((to_matrix(y_pred[1,,1:2])'+K*v[1])');
    y_up[1,3:5] = P_mat_to_array(P_array_to_mat(y_pred[1,1,3:5])-K*Lambda*P_array_to_mat(y_pred[1,1,3:5]));
// loglikelihood for the first time point
    ll+= -(log(determinant(V[1]))+v[1]'*inverse(V[1])*v[1])[1,1];

    for(t in 2:T){
    // time point t prediction step
        y_pred[t,,] = integrate_ode_rk45(DyadicInteraction, y_up[t-1,], ts[t-1,1], ts[t,], params, x, x_int,rel_tol,abs_tol, max_step);
    
    // prediction error at time point t
        v[t] = to_matrix(y[t,],2,1)-to_matrix(y_pred[t,,1:2])'; 
        V[t] = Lambda*P_array_to_mat(y_pred[t,1,3:5])*Lambda'+sigma_y;
    
        K = P_array_to_mat(y_pred[t,1,3:5])*Lambda'*inverse(V[t]);
    
    // time point t update step
        y_up[t,1:2] = to_array_1d((to_matrix(y_pred[t,,1:2])'+K*v[t])');
        y_up[t,3:5] = P_mat_to_array(P_array_to_mat(y_pred[t,1,3:5])- K*Lambda*P_array_to_mat(y_pred[t,1,3:5]));
    
    // adds to the loglikelihood object
        ll+= -(log(determinant(V[t]))+v[t]'*inverse(V[t])*v[t])[1,1];
    }
    return ll;
} // end of the CDEKF function
} // end of the function chunk

data {
int<lower=1> N; // total number of units
int<lower=1> Tmax; // the largest number of observations in all units

real y[Tmax,2,N]; // the observed data in a 3-d real-numbered array

real t0; // the initial time index of observation (we have recoded the data such that everyone started at t0=0)

int<lower=1> T[N]; // person specific max time count

real ts[Tmax,N,1]; // person specific time indices

real mu1[N]; // person specific x* and y* (defined as mean of baseline)
real mu2[N];

matrix[2,2] Lambda; // measurement matrix
real pstart[3]; // initial values of the prediction error matrix components for the CDEKF to start at (we set them all to zero in our implementation)
real ystart[2,N]; // initial values of the system latent states for the CDEKF to start at (we used the individual-specific first observed value)

// The three objects below are used as control parameters in the ODE solver.
// For more information please refer to https://mc-stan.org/docs/2_26/functions-reference/functions-ode-solver.html
real rel_tol; 
real abs_tol;
real max_step;
}

// We didn't actually need transformed data in our script. Since we didn't have covariates in the SDE system but Stan required those two arguments in the solver function, we defined them as two empty arrays

transformed data {
real x[0];
int x_int[0];
}
parameters { 
// Hyperparameters - fixed effects
	real mu_a1;
	real mu_a2;
	real mu_b1;
	real mu_b2;
// Hyperparameters - random effect standard deviations
	real<lower=0> sigma_a1;
	real<lower=0> sigma_a2;
	real<lower=0> sigma_b1;
	real<lower=0> sigma_b2;
// Individual parameters (fix + random effects)
	real a1[N];
	real a2[N];
	real b1[N];
	real b2[N];
// Process noise parameters
	real<lower=0> sigma_p1;
	real<lower=0> sigma_p2;
	real rho;
// Measurement error parameters
	real<lower=0> sigma_y1;
	real<lower=0> sigma_y2;
}
model {
// Priors and transformed parameters (that are not tracked for sampling)
int Tn;
real params[9];
real y0[5,N];
matrix[2,2] sigma_y;

sigma_y1 ~ cauchy(0,2.5);
sigma_y2 ~ cauchy(0,2.5);
sigma_p1 ~ cauchy(0,2.5);
sigma_p2 ~ cauchy(0,2.5);
sigma_a1 ~ cauchy(0,2.5);
sigma_a2 ~ cauchy(0,2.5);
sigma_b1 ~ cauchy(0,2.5);
sigma_b2 ~ cauchy(0,2.5);

a1 ~ normal(mu_a1,sigma_a1);
a2 ~ normal(mu_a2,sigma_a2);
b1 ~ normal(mu_b1,sigma_b1);
b2 ~ normal(mu_b2,sigma_b2);

sigma_y[1,1]=sigma_y1^2;
sigma_y[2,2]=sigma_y2^2;
sigma_y[1,2]=0;
sigma_y[2,1]=0;

params[7]= sigma_p1; 
params[8]= rho; 
params[9]= sigma_p2; 

// Calculating the individual-specific loglikelhood and adding it to the target density
for(n in 1:N){
  real lln;
	params[1]= a1[n];
	params[2]= a2[n]; 
	params[3]= b1[n]; 
	params[4]= b2[n]; 
	
	params[5]= mu1[n]; 
	params[6]= mu2[n]; 

    Tn=T[n];

  y0[1:2,n]=ystart[,n];
  y0[3:5,n]=pstart;

  lln=CDEKF(y[1:Tn,,n], Lambda,sigma_y, y0[,n], t0, ts[1:Tn,n,], params, x, x_int, Tn, rel_tol, abs_tol, max_step);
  target+=lln;
  }
}


