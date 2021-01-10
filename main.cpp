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
    f[i] =-f[i] + 0.5 + 111.0*k;
  }
}
void minusdtlz2(vector<double> &x, vector<double> &f)
{
  int nvar = x.size();
  int nobj = f.size();
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
    f[i] = -f[i] + 1.0 + k/4.0;;
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
    f[i] = -f[i] + 222*k;
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
    f[i] = -f[i] +1.0 + k/4.0 ;
  } // for
}

int main()
{
  vector<double> f(2), x(2);
  for(int i = 0; i < 1000000; i++)
  {
   for(auto &a:x) a = rand()/(double)RAND_MAX;
   minusdtlz4(x, f);
   for(auto &m:f) cout << m<<" ";
   cout <<endl;
  }
  return 0;
}
