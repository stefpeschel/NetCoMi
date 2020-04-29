#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

int modulo(int x, int y){
  return x - floor(x / y) * y;
}


// Draw n random samples from a dirichlet distribution
// 'a' is a vector with shape parameters
// code is taken from R-function 'rdirichlet' from 'gtools' package
NumericMatrix C_rdirichlet(int n, NumericVector a) {
  int l = a.size();
  NumericVector sm(n);
  NumericMatrix mat(n, l);
  NumericMatrix out(n, l);

  for(int i = 0; i < n; i++){
    for(int j = 0; j < l; j++){
      mat(i, j) = rgamma(1, a[j], 1)[0];
    }
    sm[i] = sum(mat(i, _));
  }

  for(int j = 0; j < l; j++){
    out(_, j) = mat(_, j) / sm;
  }
  return out;
}



// Calculate basis variances
List C_basis_var(NumericMatrix fracs, arma::mat V, arma::mat M){
  int ncol = fracs.ncol();
  arma::vec Vvec(ncol);
  arma::vec Vbase(ncol);
  arma::mat Minv(ncol, ncol);

  for(int i = 0; i < ncol; i++){
    Vvec[i] = sum(V.row(i));
  }

  Minv = inv(M);
  Vbase = Minv * Vvec;
  for(int i = 0; i < ncol; i++){
    if(Vbase[i] < 0) Vbase[i] = 1e-4;
  }

  return List::create(Named("Vbase") = Vbase, Named("M") = M);
}


// Calculate Correlations and Covariances
List C_cor_from_basis(arma::mat V, arma::vec Vbase, int p){
  int ncomb = p*(p-1)/2;
  int l = 0;
  mat Cormat(p, p);
  mat Covmat(p, p);
  uvec idx1(ncomb);
  uvec idx2(ncomb);

  for(int i = 0; i < (p-1); i++){
    for(int j = i+1; j < p; j++){
      idx1[l] = i;
      idx2[l] = j;
      l += 1;
    }
  }

  for(int i = 0; i < (p-1); i++){
    uvec select = idx2.elem(find(idx1 == i));
    vec Vtmp = V.col(i);

    vec covtmp = 0.5 * (Vbase[i] + Vbase.elem(select) - Vtmp.elem(select));
    vec denom = sqrt(Vbase[i]) * sqrt(Vbase.elem(select));
    vec cortmp = covtmp / denom;
    vec abscor = abs(cortmp);

    if(any(abscor > 1)){
      uvec idxthr = find(abscor > 1);
      int isize = idxthr.size();
      cortmp.elem(idxthr) = sign(cortmp.elem(idxthr));
      for(int k = 0; k < isize; k ++){
        int indx = idxthr[k];
        covtmp[indx] = cortmp[indx] * denom[indx];
      }
    }

    vec z = zeros<vec>(i+1);
    vec cortmpcompl = join_cols(z, cortmp);
    vec covtmpcompl = join_cols(z, covtmp);
    Cormat.col(i) = cortmpcompl;
    Covmat.col(i) = covtmpcompl;
  }

  for(int i = 0; i < p; i++){
    for(int j = 0; j < p; j++){
      if(i == j){
        Cormat(i, j) = 1;
        Covmat(i, j) = Vbase[i];
      } else{
        Cormat(i, j) = Cormat(j, i);
        Covmat(i, j) = Covmat(j, i);
      }
    }
  }

  return List::create(Named("Cormat") = Cormat, Named("Covmat") = Covmat);
}



// Find pairs with highest correlation and exclude them
List C_exclude_pairs(mat Cormat, mat M, double thresh, uvec excluded) {
  bool belowthr = false;
  int p = M.n_cols;
  mat Cortmp = abs(Cormat);

  // Remove autocorrelation
  Cortmp.diag() = Cortmp.diag() - Cormat.diag();


  // Set rows belonging to excluded OTU pairs in the correlation matrix to 0
  if(!excluded.is_empty()) {
    int exsize = excluded.size();

    for(int i = 0; i < exsize; i++){
      int rowtoex = modulo(excluded[i], p);
      Cortmp.row(rowtoex).zeros();
    }
  }

  // find OUT pair with highest correlation
  double mm = Cortmp.max();
  uvec idxtorm = find(Cortmp == mm);
  int rmsize = idxtorm.size();

  if(mm > thresh) {
    // Substract 1 to elements in M that belong to highly correlatet pairs
    M.elem(idxtorm) = M.elem(idxtorm) - 1;

    // Subtract one to the diagonal of M
    uvec rowtoex(2 * rmsize);
    for(int i = 0; i < rmsize; i++) {
      rowtoex[2*(i+1)-2] = floor(idxtorm[i] / p);
      rowtoex[2*(i+1)-1] = modulo(idxtorm[i], p);
    }
    uvec rowtoexuni = find_unique(rowtoex);
    uvec diagidx = rowtoex.elem(rowtoexuni);
    vec diagnew = M.diag();
    diagnew.elem(diagidx) = diagnew.elem(diagidx) - 1;
    M.diag() = diagnew;

    excluded = idxtorm;
  } else{
    // last iteration because correlation is below threshold
    belowthr = true;
  }

  return List::create(Named("M") = M, Named("excluded") = excluded,
                      Named("belowthr") = belowthr);

}



// Calculate correlations by iteratively excluding pairs with highest correlation
List C_compute_corr(NumericMatrix x, int exiter, double thresh) {
  int nrow = x.nrow();
  int ncol = x.ncol();
  int ncomb = ncol*(ncol-1)/2;
  int l = 0;
  double vartmp = 0;

  NumericMatrix fracs(nrow, ncol);
  NumericMatrix logfracs(nrow, ncol);
  NumericMatrix ttmp(nrow, ncomb);
  arma::mat V(ncol, ncol);
  arma::mat M(ncol, ncol);
  List ll3(3);


  // Compute matrix with relative fractions from dirichlet distribution
  for(int i = 0; i < nrow; i++) {
    fracs(i, _) = C_rdirichlet(1, (x(i, _)+1));
  }

  // Calculate logarithm of fractions
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      logfracs(i, j) = log(fracs(i, j));
    }
  }

  // create matrix Ti,j
  for(int r = 0; r < nrow; r++){
    l = 0;
    for(int i = 0; i < (ncol-1); i++){
      for(int j = i+1; j < ncol; j++){
        ttmp(r, l) = logfracs(r, i) - logfracs(r, j);
        l = l+1;
      }
    }
  }

  // Create variance matrix
  l = 0;
  for(int i = 0; i < ncol; i++){
    for(int j = i; j < ncol; j++){
      if(i == j){
        V(i, j) = 1;
      } else{
        vartmp = var(ttmp(_, l));
        V(i, j) = vartmp;
        V(j, i) = vartmp;
        l = l+1;
      }
    }
  }

  // initialise M matrix
  for(int i = 0; i < ncol; i++){
    for(int j = 0; j < ncol; j++){
      if(i == j){
        M(i, j) = ncol - 1;
      } else{
        M(i, j) = 1;
      }
    }
  }


  List ll1 = C_basis_var(fracs, V, M);
  List ll2 = C_cor_from_basis(V, ll1["Vbase"], ncol);
  uvec excluded(1);
  excluded.reset();

  for(int i = 0; i < exiter; i++) {
    ll3 = C_exclude_pairs(ll2["Cormat"], ll1["M"], thresh, excluded);
    uvec excl = ll3["excluded"];
    if(excluded.is_empty()) {
      excluded = excl;
    } else {
      excluded = join_cols(excluded, excl);
    }


    bool belowthr = ll3["belowthr"];
    if(belowthr == FALSE) {
      ll1 = C_basis_var(fracs, V, ll3["M"]);
      ll2 = C_cor_from_basis(V, ll1["Vbase"], ncol);
    }
  }

  return List::create(Named("Vbase") = ll1["Vbase"],
                      Named("Cormat") = ll2["Cormat"],
                      Named("Covmat") = ll2["Covmat"],
                      Named("fracs") = fracs);
}


//' Estimate correlations between taxa via SparCC algorithm
//'
//' @param x matrix with read counts
//' @param exiter maximum number of iterations, where in each step the
//'  strongest correlation is excluded and the basis variances and correlations
//'  are re-estimated based on the reduced data
//' @param thresh threshold for excluding component pairs in the
//'  iteration steps (pairs are excluded if the magnitude of their correlation
//'  value is upon the threshold)
//' @param repEst estimation process is repeated 'repEst' times to
//'  account for sampling noise; new fraction samples are drawn from Dirichlet
//'  distribution in each iteration step
//'
//'  @references{\insertRef{friedman2012inferring}{NetCoMi}}
// [[Rcpp::export]]
List C_sparcc(NumericMatrix x, int exiter, double thresh, int repEst) {
  int p = x.ncol();
  mat Vmat(p, repEst);
  List res(3);
  List Corlist(repEst);
  List Covlist(repEst);
  List fraclist(repEst);

  for(int e = 0; e < repEst; e++){
    res = C_compute_corr(x, exiter, thresh);
    vec Vbase = res["Vbase"];
    Vmat.col(e) = Vbase;
    Corlist[e] = res["Cormat"];
    Covlist[e] = res["Covmat"];
    fraclist[e] = res["fracs"];
  }

  return List::create(Named("Vmat") = Vmat,
                      Named("Corlist") = Corlist,
                      Named("Covlist") = Covlist,
                      Named("fraclist") = fraclist);
}


