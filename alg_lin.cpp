#include <iostream>
#include <vector>
#include <fstream>
// #include <bits/stdc++.h>

using namespace std;

class Matrix                      
{
  public:                     
    Matrix(vector<vector<float> > input);       
      void escalona_gauss();
      void trata_matriz(int i);
      float get_element(int i, int j){return mat[i][j];};
      int det_laplace(vector<vector<float>> matriz);
      vector<float> eliminacao_gauss();
      vector<vector<float>> fatoracao_lu();


      vector<float> sol_iter_jacobi(float tol = 0.000001, int numero_iteracoes = 0, bool usar_numero_iteracoes = false);
    
 private:                      // begin private section
    int n;                // member variable
    vector<vector<float> > mat;
    int dim;
    int det;
};


Matrix::Matrix(vector<vector<float> > input)
{
   mat = input;
   n = mat.size()-1; //nós fazemos -1 porque no lado direito da matrix temos o vetor de coeficientes b
}

void Matrix::trata_matriz(int i){

   for(int j = i + 1; j<n; j++){
      if(mat[j][i]!=0){ //achamos uma linha com pivo não nulo
         swap(mat[i], mat[j]);
         return; 
      }
   }
   //se não achar printar erro porque a matriz é LD e temos infinitas solucoes?
   cout << "A matiz é LD " << "\n";
   return;
}

int Matrix::det_laplace(vector<vector<float>> matriz){

   if(matriz.size()==1){
      return matriz[0][0];
   }
   if(matriz.size()==2){
      return matriz[0][0]*matriz[1][1] - matriz[0][1]*matriz[1][0];
   }
   else{
      vector<float> linha = matriz[0];
      vector<vector<float>> matriz_aux(n, vector<float>(n,0));
   }
   return 0;
}

void Matrix::escalona_gauss(){

   // primeiro checamos se é possível escalonar a matrix
   // det_laplace();
   // if (det_laplace()<0){
   //    cout << "Como o detemrinante é negativo não conseguimos encontrar solução" << "\n";
   //    return;
   // }
   for(int d = 0; d < n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
      for(int i = d+1; i <= n; i++){ //seleciona a linha 
         if (mat[d][d]==0){ //achamos uma linha com pivo zero
            trata_matriz(d); 
         }
         float multiplicador = mat[i][d]/mat[d][d];
         for(int j = 0; j <= n+1; j++){ //anda na linha 
            mat[i][j] = mat[i][j] - multiplicador * mat[d][j]; 
            }
      }
   }
   return;
}

vector<float> Matrix::eliminacao_gauss(){

   escalona_gauss();
   vector<float> x(n+1,0);
   // sabemos o último valor de x, o x_n
   x[n] = mat[n][n+1]/mat[n][n];

   for(int i = n-1; i>=0; i--){
      float sum = 0;
      for(int j = i+1; j<= n; j++){
         sum += mat[i][j]*x[j];
      }
      x[i] = (mat[i][n+1] - sum)/mat[i][i];
   }
   return x;
}

vector<float> Matrix::sol_iter_jacobi(float tol, int numero_iteracoes , bool usar_numero_iteracoes){
   
   // adicionar check se pode aplicar porque todos os elementos da diagonal tem que ser maior do que os elementos da linha e da coluna correpsondente 
   vector<float> x_old(n+1, 1.0); //iniciamos o vetor solucao inicial com todos os elementos como 1
   float r = 1;
   int iteracoes = 1;

   ofstream resultsFile("resultados_sol_iterativa.txt");

   if(resultsFile.is_open()){
      cout<< "arquivo criado com sucesso";
      cout <<"\n";
   }
   else{cerr << "Erro criando o arquivo " << endl;}

   if(usar_numero_iteracoes == false){
      while(r>tol){
         vector<float> x_new(n+1,1);
         for(int i = 0; i<= n; i++){
            float sum = 0;
            for(int j = 0; j <= n; j++){
               if(i!=j){
                  sum+= mat[i][j]*x_old[j];
               }
            }
            x_new[i] = (mat[i][n+1] - sum)/mat[i][i];
         }

         // calculamos o r para comparara com a tolerância
         float up_sum_square = 0;
         float down_sum_square = 0;
         for(int i = 0; i <= n-1; i++){
            up_sum_square += (x_old[i]-x_new[i])*(x_old[i]-x_new[i]);
            down_sum_square += x_new[i]*x_new[i];
         }
         r = (sqrt(up_sum_square)/sqrt(down_sum_square));
         resultsFile << "Numero de iterações: " << iteracoes << "  Norma do resíduo: " << r << endl;
         x_old = x_new; // agora usamos o mais novo como sendo o velho para a próxima iteração
         iteracoes++;
      }
   }
   else{
      while(iteracoes<= numero_iteracoes){
         vector<float> x_new(n+1,1);
         for(int i = 0; i<= n; i++){
            float sum = 0;
            for(int j = 0; j <= n; j++){
               if(i!=j){
                  sum+= mat[i][j]*x_old[j];
               }
            }
            x_new[i] = (mat[i][n+1] - sum)/mat[i][i];
         }

         // calculamos o r para comparara com a tolerância
         float up_sum_square = 0;
         float down_sum_square = 0;
         for(int i = 0; i <= n-1; i++){
            up_sum_square += (x_old[i]-x_new[i])*(x_old[i]-x_new[i]);
            down_sum_square += x_new[i]*x_new[i];
         }
         r = (sqrt(up_sum_square)/sqrt(down_sum_square));
         cout << "iteacoes : " << iteracoes << " norma do resíduo: " << r<< " " << endl;
         x_old = x_new; // agora usamos o mais novo como sendo o velho para a próxima iteração
         iteracoes++;
      }
   }
   resultsFile.close();

   // ná última iteração do algoritmo o x_old é o mais recente vetor de x, é a solução do sistema

   return x_old;
}

vector<vector<float>>Matrix::fatoracao_lu(){
   //
   vector<vector<float>> matriz_L(n+1,vector<float>(n+1));
   vector<vector<float>> matriz_U(n+1,vector<float>(n+1));
   for(int i = 0; i<=n;i++){
      // for(int j = 0; j <= n; j++){
      //    matriz_U[i][j] = mat[i][j];
      // }
      matriz_U[i]=mat[i];
   }
   
   for(int d = 0; d <= n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
      matriz_L[d][d] = 1;
      for(int i = d+1; i <= n; i++){ //seleciona a linha 


         float multiplicador = matriz_U[i][d]/matriz_U[d][d];
         matriz_L[i][d] = multiplicador;
         for(int j = 0; j <= n; j++){ //anda na linha 
            matriz_U[i][j] = matriz_U[i][j] - multiplicador * matriz_U[d][j]; 
            }
      }
   }
   return matriz_U;

}

int main()
{   
   vector<vector<float>> fds({ {3, -1, -1, 1}, 
                                {-1, 3, -1, 2}, 
                                {-1, -1, 3, 1} });
   cout << fds.size() << "\n";
   vector<vector<float>> aa({ {1,2,2,3},
                              {4,4,2,6},
                              {4,6,4,10}
   });
   Matrix m (fds);
   Matrix mmmm(aa);
   vector<float> sol;
   sol =  m.sol_iter_jacobi();
   for(int i = 0; i<sol.size();i++){
      cout << sol[i] << " ";
   }
   cout<< "\n";
   cout<< "\n";
   cout<< "\n";

   vector<float> resolveai;
   mmmm.escalona_gauss();

   for(int i = 0; i <=2;i++){
      for(int j = 0 ; j<=3; j++){
      cout << mmmm.get_element(i,j) << "  ";
      }
   cout << "\n"; 
   }

   resolveai = mmmm.eliminacao_gauss();

   for(int i = 0; i<resolveai.size();i++){
      cout << resolveai[i] << " ";
   }
   cout<< "\n";

   cout<< "\n";

   cout<< "\n";

   cout<< "\n";
 
   vector<vector<float>> jeson({
                              {1,2,2,3},
                              {4,4,2,6},
                              {4,6,4,10}
   });

   // Matrix jerson(jeson);
   // vector<vector<float>> aaaaaa;
   // aaaaaa = jerson.fatoracao_lu();

   // for(int i = 0; i<=2;i++){
   //    for(int j = 0; j<=2; j++){
   //       cout << aaaaaa[i][j]<< "  ";
   //    }
   //    cout<< "\n";
   // }
   return 0;
}  