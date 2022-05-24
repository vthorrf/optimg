// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
#include <RcppArmadilloExtensions/sample.h>

using namespace Rcpp;
using namespace arma;

// ====================================Supporting functions====================================
// [[Rcpp::export]]
Rcpp::NumericVector grad(Function Model, NumericVector par, double h=1e-6) {
    NumericMatrix mat(par.length(), par.length());
    NumericVector out(par.length());
    double        f_x = as<double>(Model(par));

    for (int i = 0; i < mat.ncol(); i++) {
        mat(i, _) = par;
        mat(i, i) = mat(i, i) + h;
        out[i] = (as<double>(Model(mat(i,_))) - f_x) / h;
    }

    return out;
}

int binCoef(int n, int k) {
  if (k == 0 || k == n)
    return 1;
  return binCoef(n - 1, k - 1) + binCoef(n - 1, k);
}

// [[Rcpp::export]]
Rcpp::NumericVector gradN(Function f, NumericVector par, double h=1e-6, int order=1) {
  int m = par.length();
  NumericMatrix df(m,order);
  IntegerVector posit(m);
  for(int i = 0; i < m; ++i) {
    posit[i] = i;
  }
  for (int i = 0; i < m; ++i) {
    for(int j = 0; j < order; ++j) {
      NumericVector temp(m);
      for(int k = 0; k < m; ++k) {
        double t0;
        if (posit[k] == i) {
          t0 = 1.0;
        } else {
          t0 = 0.0;
        }
        temp[k] = par[k] + (t0 * double(j+1) * h);
      }
      df(i,j) = pow(-1.0, double(order-j+1)) * double(binCoef(order,j+1)) * as<double>(f(temp));
      if(!arma::is_finite(df(i,j))) {
        df(i,j) = 0;
      }
    }
  }
  
  NumericMatrix Final(m,order+1);
  double fix = as<double>(f(par)) * pow(-1, order % 2);
  for(int i = 0; i < m; ++i) {
    Final(i,0) = fix;
  }
  for(int j = 1; j < (order+1); ++j) {
    Final(_,j) = df(_,(j-1));
  }
  NumericVector out(m);
  for(int i = 0; i < m; ++i) {
    out[i] = sum(Final(i,_))/h;
  }
  
  return out;
}

double reltol(double cost, double tol=1e-8) {
    return tol * (abs(cost) + tol);
}

IntegerVector order_(NumericVector x) {
    // Picking up order() function from base package
    Function f("order");
    return as<IntegerVector>(f(x))-1;
}

double mean_v(NumericVector vec, int ele) {
    int vs = vec.size();
    double res = 0;
    for (int i = 0; i < vs; i++) {
        if (i == ele) {
            continue;
        }
        res += vec[i];
    }
    return res / (vs-1);
}

double crit_v(NumericMatrix Obj) {
    int ms = Obj.cols();
    NumericVector x1(ms);
    NumericMatrix matx(ms, ms);
    for (int iter = 0; iter < ms; iter++) {
        for (int j = 0; j < ms; j++) {
            matx(iter, j) = abs(Obj(iter + 1, j) - Obj(iter, j));
        }
        x1[iter] = max(matx(iter, _));
    }
    return max(x1);
}

double steepSGD(NumericVector par, List Data) {
  Rcpp::NumericVector der = Data["Gradients"];
  Rcpp::NumericVector init = Data["Estimates"];
  Rcpp::Function FUN = Data["FUN"];
  
  return as<double>(FUN(init - (par * der)));
}

double steepSTGD(NumericVector par, List Data) {
  Rcpp::NumericVector der = Data["Gradients"];
  Rcpp::NumericVector init = Data["Estimates"];
  Rcpp::Function FUN = Data["FUN"];
  Rcpp::NumericVector New(init.length());
  double SDD = sd(der);
  for (int iter = 0; iter < New.length(); iter++) {
    New[iter] = abs(der[iter]) > SDD ? init[iter] - (par[0] * der[iter]) : init[iter] - (par[1] * der[iter]);
  }
  
  return as<double>(FUN(New));
}

double steepLMM(NumericVector par, List Data) {
  arma::vec der = as<arma::vec>(wrap(Data["Gradients"]));
  Rcpp::NumericVector init = Data["Estimates"];
  Rcpp::Function FUN = Data["FUN"];
  Rcpp::NumericMatrix I = NumericMatrix::diag(init.length(), exp(par[1]));
  arma::mat SOG = (der * der.t()) + as<arma::mat>(wrap(I));
  Rcpp::NumericVector LM = as<Rcpp::NumericVector>(wrap(arma::inv(SOG) * der));
  
  return as<double>(FUN(init - (par[0] * LM)));
}

double fibonacci(double lower, double upper, List Data, int n=100) {
  double epsilon = .01;
  double a = lower;
  double b = upper;
  double phi = (1.0 + sqrt(5.0)) / 2.0;
  double S = (1.0 - sqrt(5.0)) / phi;
  double rho = 1.0 / (phi * (1 - pow(S, n + 1)) / (1 - pow(S, n)));
  double d = (rho * b) + ((1 - rho) * a);
  double yd = steepSGD(d, Data);
  double c, yc;
  for (int i = 0; i < n-1; i++) {
    if(a == b) {
      break;
    }
    if (i == (n - 1)) {
      c = (epsilon * a) + ((1 - epsilon) * d);
    }
    else {
      c = (rho * a) + ((1 - rho) * b);
    }
    yc = steepSGD(c, Data);
    if (yc < yd) {
      b = d, d = c, yd = yc;
    }
    else {
      a = b, b = c;
    }
    rho = 1.0 / (phi * (1 - pow(S, n - i + 1)) / (1 - pow(S, n - i)));;
  }
  return (a + b) / 2.0;
}

typedef double (*funcPtr)(NumericVector par, List Data);

XPtr<funcPtr> putFunPtrInXPtr(std::string fstr) {
  if (fstr == "steepSGD")
    return XPtr<funcPtr>(new funcPtr(&steepSGD));
  else if (fstr == "steepSTGD")
    return(XPtr<funcPtr>(new funcPtr(&steepSTGD)));
  else if (fstr == "steepLMM")
    return(XPtr<funcPtr>(new funcPtr(&steepLMM)));
  else
    return XPtr<funcPtr>(R_NilValue); // runtime error as NULL no XPtr
}

Rcpp::List nelder_mead(std::string funname, NumericVector init, List Data) {
    XPtr<funcPtr> xpfun = putFunPtrInXPtr(funname);
    funcPtr fn = *xpfun;
    double epsilon = 1e-10;
    int d = init.length();
    int d1 = d + 1;
    double alpha, beta, gamma, sigma;
    int maxfeval = 50 * pow(d, 2);
    if (d > 1) {
        alpha = 1.0, beta = 1.0 + 2.0 / d, gamma = .75 - (.5 / d), sigma = 1.0 - (1.0/d);
    }
    else {
        alpha = 1.0, beta = 2.0, gamma = .50, sigma=.50;
    }
    double xi = min(NumericVector::create(max(NumericVector::create(max(abs(init)), 1.0)), 10.0));

    // Initial settings
    NumericMatrix S(d1, d);
    for (int iter = 0; iter < d1; iter++) {
        for (int j = 0; j < d; j++) {
            S(iter, j) = iter == j ? 1.0 : iter == d ? (1.0 - sqrt(d1))/d : 0.0;
        }
    }
    NumericMatrix X(d1, d);
    NumericVector y_arr(d1);
    for (int iter = 0; iter < d1; iter++) {
        X(iter,_) = init + xi * S(iter,_);
        y_arr[iter] = fn(X(iter,_), Data);
    }
    // Ordering
    IntegerVector p = order_(y_arr);
    NumericMatrix SI(d1, d);
    NumericVector y_I(d1);
    for (int iter = 0; iter < d1; iter++) {
        SI(iter,_) = X(p[iter],_);
        y_I[iter]  = y_arr[p[iter]];
    }
    X = SI;
    y_arr = y_I;
    int ct = d1;
    double delta = crit_v(X);

    while (delta >= (xi * epsilon)) {
        if (ct > maxfeval) {
            break;
        }
        // Simplex
        NumericVector xl = X(0,_);
        double yl = y_arr[0];
        NumericVector xh = X(d,_);
        double yh = y_arr[d];
        NumericVector xs = X(d-1, _);
        double ys = y_arr[d-1];
        // Centroid
        NumericVector xm(d);
        for (int iter = 0; iter < d; iter++) {
            xm[iter] = mean_v(X(_,iter), d);
        }
        // Reflection
        NumericVector xr = xm + alpha * (xm - xh);
        double yr = fn(xr, Data);
        ct += 1;

        // Expansion
        if (yr < yl) {
            NumericVector xe = xm + beta * (xr - xm);
            double ye = fn(xe, Data);
            X(d, _) = ye < yr ? xe : xr;
            y_arr[d] = ye < yr ? ye : yr;
            ct += 1;
        }
        else if (yr >= ys) {
            if (!(yr >= yh)) {
                xh = xr;
                yh = yr;
                X(d, _) = xr;
                y_arr[d] = yr;
            }
            // Contraction
            NumericVector xc = xm + gamma * (xh - xm);
            double yc = fn(xc, Data);
            ct += 1;
            // Shrinkage
            if (yc > yh) {
                for (int iter = 1; iter < d1; iter++) {
                    X(iter, _) = xl + sigma * (X(iter, _) - xl);
                    y_arr[iter] = fn(X(iter, _), Data);
                }
                ct += d;
            }
            else {
                X(d,_) = xc;
                y_arr[d] = yc;
            }
        }
        else {
            X(d, _) = xr;
            y_arr[d] = yr;
        }
        p = order_(y_arr);
        NumericMatrix SI(d1, d);
        NumericVector y_I(d1);
        for (int iter = 0; iter < d1; iter++) {
            SI(iter, _) = X(p[iter], _);
            y_I[iter] = y_arr[p[iter]];
        }
        X = SI;
        y_arr = y_I;
        delta = crit_v(X);
    }

    return wrap(Rcpp::List::create(Rcpp::Named("par") = as<NumericVector>(wrap(X(0,_))),
                Rcpp::Named("value") = y_arr[0],
                Rcpp::Named("nfeval") = ct,
                Rcpp::Named("convergence") = delta > epsilon));
}

// ====================================Gradient Descent Algorithms====================================
// [[Rcpp::export]]
SEXP SGD(Function fn, NumericVector startvalue, Function gr, double h=1e-6,
         int maxit=10, double tol=1e-8, bool verbose=false) {
    // Opening message
    if (verbose == true) {
        Rcpp::Rcout << "Steepest Gradient Descent will run for " << maxit << " iterations at most." << std::endl;
    }

    // Initial settings
    if(h > tol) {
      tol = h;
    }
    NumericMatrix par(maxit, startvalue.length());
    NumericMatrix df(maxit, startvalue.length());
    NumericVector f(maxit);
    NumericVector step(maxit);
    bool convergence = false;
    int final = 0;
    
    // First step
    df(0, _) = as<NumericVector>(gr(startvalue));
    List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = startvalue,
                                      Rcpp::Named("Gradients") = df(0, _),
                                      Rcpp::Named("FUN") = fn);
    step[0] = fibonacci(-1.0, 1.0, DataSGD, 1000);
    NumericVector sss(1);
    sss[0] = step[0];
    par(0,_) = startvalue - (step[0] * df(0,_));
    f[0] = as<double>(fn(par(0, _)));
    if (verbose == true) {
      Rcpp::Rcout << "Iteration: " << 1 << ". Cost: " << f[0] << std::endl;
    }

    // Next steps
    for (int run = 1; run < maxit; run++) {
        // Gradients
        df(run, _) = as<NumericVector>(gr(par(run-1, _)));
        // Steepest step
        List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = par(run-1,_),
                                          Rcpp::Named("Gradients") = df(run, _),
                                          Rcpp::Named("FUN") = fn);
        //step[run] = as<double>(wrap(nelder_mead("steepSGD", sss, DataSGD)["par"]));
        step[run] = fibonacci(-1.0, 1.0, DataSGD, 1000);
        // New parameter
        par(run, _) = par(run-1,_) - (step[run] * df(run, _));
        f[run] = as<double>(fn(par(run, _)));
        if (verbose == true) {
          Rcpp::Rcout << "Iteration: " << run + 1 << ". Cost: " << f[run] << std::endl;
        }
        final += 1;
        // Check convergence
        if ( reltol(f[run], tol) > abs(f[run] - f[run-1]) ) {
          convergence = true;
          break;
        }
    }
    if (verbose == true) {
      if (convergence == true) {
        Rcpp::Rcout << "Convergence achieved!" << std::endl;
      }
      else {
        Rcpp::Rcout << "Convergence may not have been achieved!" << std::endl;
      }
    }
    f    = f[Range(0,final)];
    par  = par(Range(0,final),_);
    df   = df(Range(0,final),_);
    step = step[Range(0,final)];
    // Final Result
    return wrap(Rcpp::List::create(Rcpp::Named("Cost") = f,
                Rcpp::Named("Estimates") = par,
                Rcpp::Named("Gradients") = df,
                Rcpp::Named("steps") = step,
                Rcpp::Named("MaxInt") = final + 1,
                Rcpp::Named("convergence") = convergence));
}

// [[Rcpp::export]]
SEXP STGD(Function fn, NumericVector startvalue, Function gr, double h=1e-6,
          int maxit=10, double tol=1e-8, bool verbose=false) {
    // Opening message
    if (verbose == true) {
        Rcpp::Rcout << "Steepest 2-group Gradient Descent will run for " << maxit <<
            " iterations at most." << std::endl;
    }

    // Initial settings
    if(h > tol) {
      tol = h;
    }
    NumericMatrix par(maxit, startvalue.length());
    NumericMatrix df(maxit, startvalue.length());
    NumericVector f(maxit);
    NumericMatrix step(maxit, 2);
    NumericVector sss = { .001,.001 };
    bool convergence = false;
    int final = 0;

    // First step
    df(0, _) = as<NumericVector>(gr(startvalue));
    // Steepest step
    List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = startvalue,
                                      Rcpp::Named("Gradients") = df(0, _),
                                      Rcpp::Named("FUN") = fn);
    step(0, _) = as<NumericVector>(wrap(nelder_mead("steepSTGD", sss, DataSGD)["par"]));
    for (int iter = 0; iter < startvalue.length(); iter++) {
        par(0, iter) = abs(df(0, iter)) > sd(df(0, _)) ? startvalue[iter] - (step(0, 0) * df(0, iter)) : startvalue[iter] - (step(0, 1) * df(0, iter));
    }
    f[0] = as<double>(fn(par(0, _)));
    if (verbose == true) {
      Rcpp::Rcout << "Iteration: " << 1 << ". Cost: " << f[0] << std::endl;
    }

    // Next steps
    for (int run = 1; run < maxit; run++) {
        // Gradients
        df(run, _) = as<NumericVector>(gr(par(run - 1, _)));
        // Steepest steps
        List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = par(run-1,_),
                                          Rcpp::Named("Gradients") = df(run, _),
                                          Rcpp::Named("FUN") = fn);
        step(run, _) = as<NumericVector>(wrap(nelder_mead("steepSTGD", sss, DataSGD)["par"]));
        // New parameter
        for (int iter = 0; iter < startvalue.length(); iter++) {
            par(run, iter) = abs(df(run, iter)) > sd(df(run, _)) ? par(run - 1, iter) - (step(run, 0) * df(run, iter)) : par(run - 1, iter) - (step(run, 1) * df(run, iter));
        }
        f[run] = as<double>(fn(par(run, _)));
        if (verbose == true) {
          Rcpp::Rcout << "Iteration: " << run + 1 << ". Cost: " << f[run] << std::endl;
        }
        final += 1;
        // Check convergence
        if ( reltol(f[run], tol) > abs(f[run] - f[run-1]) ) {
          convergence = true;
          break;
        }
    }
    if (verbose == true) {
      if (convergence == true) {
        Rcpp::Rcout << "Convergence achieved!" << std::endl;
      }
      else {
        Rcpp::Rcout << "Convergence may not have been achieved!" << std::endl;
      }
    }
    f    = f[Range(0,final)];
    par  = par(Range(0,final),_);
    df   = df(Range(0,final),_);
    step = step(Range(0,final),_);
    // Final Result
    return wrap(Rcpp::List::create(Rcpp::Named("Cost") = f,
                Rcpp::Named("Estimates") = par,
                Rcpp::Named("Gradients") = df,
                Rcpp::Named("steps") = step,
                Rcpp::Named("MaxInt") = final+1,
                Rcpp::Named("convergence") = convergence));
}

// [[Rcpp::export]]
SEXP LMM(Function fn, NumericVector startvalue, Function gr, double h = 1e-6,
         int maxit = 10, double tol = 1e-8, bool verbose = false) {
    // Opening message
    if (verbose == true) {
        Rcpp::Rcout << "Levenberg-Marquardt method will run for " << maxit << " iterations at most." << std::endl;
    }

    // Initial settings
    if (h > tol) {
        tol = h;
    }
    NumericMatrix par(maxit, startvalue.length());
    NumericMatrix df(maxit, startvalue.length());
    NumericVector f(maxit);
    NumericMatrix step(maxit, 2);
    NumericVector sss = { .001,.001 };
    bool convergence = false;
    int final = 0;

    // First step
    df(0, _) = as<NumericVector>(gr(startvalue));
    List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = startvalue,
                                      Rcpp::Named("Gradients") = df(0, _),
                                      Rcpp::Named("FUN") = fn);
    //sss = fibonacci(-1.0, 1.0, DataSGD, 1000);
    //par(0, _) = startvalue - (sss * df(0, _));
    step(0, _) = as<NumericVector>(wrap(nelder_mead("steepLMM", sss, DataSGD)["par"]));
    Rcpp::NumericMatrix I = NumericMatrix::diag(startvalue.length(), exp(step(0, 1)));
    arma::mat SOG = (as<arma::vec>(wrap(df(0, _))) * as<arma::vec>(wrap(df(0, _))).t()) + as<arma::mat>(wrap(I));
    Rcpp::NumericVector LM = as<Rcpp::NumericVector>(wrap(arma::inv(SOG) * as<arma::vec>(wrap(df(0, _)))));
    par(0, _) = startvalue - (step(0, 0) * LM);
    f[0] = as<double>(fn(par(0, _)));
    if (verbose == true) {
        Rcpp::Rcout << "Iteration: " << 1 << ". Cost: " << f[0] << std::endl;
    }

    // Next steps
    for (int run = 1; run < maxit; run++) {
        // Gradients
        df(run, _) = as<NumericVector>(gr(par(run - 1, _)));
        // Steepest step
        List DataSGD = Rcpp::List::create(Rcpp::Named("Estimates") = par(run-1,_),
                                          Rcpp::Named("Gradients") = df(run, _),
                                          Rcpp::Named("FUN") = fn);
        step(run, _) = as<NumericVector>(wrap(nelder_mead("steepLMM", sss, DataSGD)["par"]));
        // New parameter
        Rcpp::NumericMatrix I = NumericMatrix::diag(startvalue.length(), exp(step(run,1)));
        arma::mat SOG = (as<arma::vec>(wrap(df(run, _))) * as<arma::vec>(wrap(df(run, _))).t()) + as<arma::mat>(wrap(I));
        Rcpp::NumericVector LM = as<Rcpp::NumericVector>(wrap(arma::inv(SOG) * as<arma::vec>(wrap(df(run, _)))));
        par(run, _) = par(run-1, _) - (step(run,0) * LM);
        f[run] = as<double>(fn(par(run, _)));
        if (verbose == true) {
            Rcpp::Rcout << "Iteration: " << run + 1 << ". Cost: " << f[run] << std::endl;
        }
        final += 1;
        // Check convergence
        if (reltol(f[run], tol) > abs(f[run] - f[run - 1])) {
            convergence = true;
            break;
        }
    }
    if (verbose == true) {
        if (convergence == true) {
            Rcpp::Rcout << "Convergence achieved!" << std::endl;
        }
        else {
            Rcpp::Rcout << "Convergence may not have been achieved!" << std::endl;
        }
    }
    f = f[Range(0, final)];
    par = par(Range(0, final), _);
    df = df(Range(0, final), _);
    step = step(Range(0, final), _);
    // Final Result
    return wrap(Rcpp::List::create(Rcpp::Named("Cost") = f,
                Rcpp::Named("Estimates") = par,
                Rcpp::Named("Gradients") = df,
                Rcpp::Named("steps") = step,
                Rcpp::Named("MaxInt") = final+1,
                Rcpp::Named("convergence") = convergence));
}

// ====================================THE END====================================