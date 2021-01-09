#include <bits/stdc++.h>
using namespace std;
void minusdtlz1(vector<double> &x, vector<double> &f)
{
  int nobj =f.size(), nvar = x.size();
  int k = nvar - nobj + 1;
  double g = 0.0 ;
  for (int i = nvar - k; i < nvar; i++)
    g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * M_PI * (x[i] - 0.5));
  g = 100 * (k + g);
  for (int i = 0; i < nobj; i++)
    f[i] = (1.0 + g) * 0.5;

  for (int i = 0; i < nobj; i++)
  {
    for (int j = 0; j < nobj - (i + 1); j++)  f[i] *= x[j];
      if (i != 0)
      {
        int aux = nobj - (i + 1);
        f[i] *= 1.0 - x[aux];
      } 
    f[i] *=-1.0;
  }
}
void minusdtlz2(vector<double> &x, vector<double> &f)
{
  int nvar = x.size();
  int nobj = x.size();
  int k = nvar - nobj + 1;
   double g = 0.0;
  for (int i = nvar- k; i < nvar; i++)
    g += (x[i] - 0.5)*(x[i] - 0.5);

  for (int i = 0; i < nobj; i++)
    f[i] = 1.0 + g;

  for (int i = 0; i < nobj; i++){
    for (int j = 0; j < nobj - (i + 1); j++)
      f[i] *= cos(x[j]*0.5*M_PI);
      if (i != 0){
        int aux = nobj - (i + 1);
        f[i] *= sin(x[aux]*0.5*M_PI);
      } //if
    f[i] *=-1.0;
  } // for
}
void minusdtlz3(vector<double> &x, vector<double> &f)
{
  int nvar = x.size();
  int nobj = f.size();

  int k = nvar - nobj + 1;


  double g = 0.0;
  for (int i = nvar - k; i < nvar; i++)
    g += (x[i] - 0.5)*(x[i] - 0.5) - cos(20.0 * M_PI * (x[i] - 0.5));

  g = 100.0 * (k + g);
  for (int i = 0; i < nobj; i++)
    f[i] = 1.0 + g;

  for (int i = 0; i < nobj; i++){
    for (int j = 0; j < nobj - (i + 1); j++)
      f[i] *= cos(x[j]*0.5*M_PI);
      if (i != 0){
        int aux = nobj - (i + 1);
        f[i] *= sin(x[aux]*0.5*M_PI);
      } // if
    f[i] *=-1.0;
  } //for
}
void minusdtlz4(vector<double> &x, vector<double> &f)
{
  int nvar = x.size();
  int nobj = f.size();
  int k = nvar - nobj + 1;
  double alpha = 100.0;


  double g = 0.0;
  for (int i = nvar - k; i < nvar; i++)
    g += (x[i] - 0.5)*(x[i] - 0.5);

  for (int i = 0; i < nobj; i++)
    f[i] = 1.0 + g;

  for (int i = 0; i < nobj; i++) {
    for (int j = 0; j < nobj - (i + 1); j++)
      f[i] *= cos(pow(x[j],alpha)*(M_PI/2.0));
      if (i != 0){
        int aux = nobj - (i + 1);
        f[i] *= sin(pow(x[aux],alpha)*(M_PI/2.0));
      } //if
    f[i] *=-1.0;
  } // for
}
void minusdtlz5(vector<double> &x, vector<double> &f)
{
 double g = 0.0;
  int nvar = x.size();
  int nobj = f.size();
  std::vector<double> theta_(nobj-1,0);
  int k = nvar - nobj + 1;
  double alpha = 100.0;


  for (int i = nvar - k; i < nvar; i++)
    g += (x[i] - 0.5)*(x[i] - 0.5);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = x[0] * M_PI / 2.0;
  for (int i = 1; i < (nobj-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * x[i]);

  for (int i = 0; i < nobj; i++)
    f[i] = 1.0 + g;

  for (int i = 0; i < nobj; i++){
    for (int j = 0; j < nobj - (i + 1); j++)
      f[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = nobj - (i + 1);
        f[i] *= sin(theta_[aux]);
      } // if
    f[i] *=-1.0;
  } //for
}
void minusdtlz6(vector<double> &x, vector<double> &f)
{
  double g = 0.0;
  int nvar = x.size();
  int nobj = f.size();

  std::vector<double> theta_(nobj-1,0);
  int k = nvar - nobj + 1;
  double alpha = 100.0;

  for (int i = nvar - k; i < nvar; i++)
    g += pow(x[i],0.1);

  double t = M_PI / (4.0 * (1.0 + g));

  theta_[0] = x[0] * M_PI / 2.0;
  for (int i = 1; i < (nobj-1); i++)
    theta_[i] = t * (1.0 + 2.0 * g * x[i]);

  for (int i = 0; i < nobj; i++)
    f[i] = 1.0 + g;

  for (int i = 0; i < nobj; i++){
    for (int j = 0; j < nobj - (i + 1); j++)
      f[i] *= cos(theta_[j]);
      if (i != 0){
        int aux = nobj - (i + 1);
        f[i] *= sin(theta_[aux]);
      } // if
    f[i] *=-1.0;
  } //for
}
void minusdtlz7(vector<double> &x, vector<double> &f)
{
  double g = 0.0;
  int nvar = x.size();
  int nobj = f.size();
  int k = nvar - nobj + 1;
  double alpha = 100.0;

  for (int i = nvar - k; i < nvar; i++)
    g += x[i] ;

  g = 1 + (9.0 * g)/k ;


  for (int i = 0; i < nobj - 1; i++)
    f[i] = x[i] ;

  double h = 0.0 ;
  for (int i = 0; i < nobj - 1; i++){
    h+=(f[i]/(1.0+g))*(1 + sin(3.0*M_PI*f[i])) ;
  } //for
  h = nobj - h ;
  f[nobj - 1] = (1+g)*h ;
  for (int i = 0; i < nobj; i++) f[i] *=-1.0;
}
int main()
{
  vector<double> f(2), x(10);
  for(int i = 0; i < 100000; i++)
  {
   for(auto &a:x) a = rand()/(double)RAND_MAX;
   minusdtlz1(x, f);
   for(auto &m:f) cout << -m<<" ";
   cout <<endl;
  }
  return 0;
}
