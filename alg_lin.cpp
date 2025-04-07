#include <iostream>
#include <vector>
#include <fstream>
// #include <bits/stdc++.h>

using namespace std;

class Matrix                      
{
  public:                     
    Matrix(vector<vector<float> > input);       
      void trata_matriz(int i); // Tenta trocar a linha atual com uma abaixo que tenha pivô diferente de zero. Necessário para evitar divisões por zero durante a eliminação de Gauss
      float get_element(int i, int j){return mat[i][j];};
      vector<vector<float>> get_matriz_A(const vector<vector<float>>& mat_completa); //Obtem a matriz sem o lado direito para aplicar laplace. Retorna apenas a parte A da matriz aumentada [A | b], ignorando o vetor b
      int det_laplace(const vector<vector<float>>& submat); // Calcula o determinante de uma matriz de forma recursiva usando a expansão de Laplace. Método ineficiente para matrizes grandes, mas serve para fins didáticos

      void escalona_gauss(); // Aplica a eliminação de Gauss para triangularizar a matriz. Checa primeiro se o determinante é negativo como condição simplificada de existência de solução
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

vector<vector<float>> Matrix::get_matriz_A(const vector<vector<float>>& mat_completa) {
    int linhas = mat_completa.size();
    int colunas = mat_completa[0].size() - 1; // Ignora a última coluna
    vector<vector<float>> A(linhas, vector<float>(colunas));

    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            A[i][j] = mat_completa[i][j];
        }
    }

    return A;
}

int Matrix::det_laplace(const vector<vector<float>>& submat) {
    int tamanho = submat.size();

    if (tamanho == 1) {
        return submat[0][0];
    }

    if (tamanho == 2) {
        return submat[0][0] * submat[1][1] - submat[0][1] * submat[1][0];
    }

    float determinante = 0.0;

    for (int k = 0; k < tamanho; ++k) {
        vector<vector<float>> minor;

        for (int i = 1; i < tamanho; ++i) {
            vector<float> linha;
            for (int j = 0; j < tamanho; ++j) {
                if (j != k) {
                    linha.push_back(submat[i][j]);
                }
            }
            minor.push_back(linha);
        }

        float cofator = pow(-1, k) * submat[0][k] * det_laplace(minor);
        determinante += cofator;
    }

    return determinante;
}



void Matrix::escalona_gauss(){

   // primeiro checamos se é possível escalonar a matrix
   int determinante = det_laplace(get_matriz_A(mat));
   if (determinante<0){
      cout << "Como o detemrinante é negativo não conseguimos encontrar solução" << "\n";
      return;
   }
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
   
   // antes de iniciar o algoritmo devemos checar se cada elemento da diagonal é o maior elemento da sua linha e coluna correspondente
   for(int i = 0; i<n; i++){
      float dig_element = fabsf(mat[i][i]);
      float sum_col = 0;
      float sum_row = 0;
      for(int j = 0; j<n; j++){
         if(j != i){
            sum_col+= fabsf(mat[j][i]);
            sum_row+= fabs(mat[i][j]);
         }
      }
      if(!(dig_element>sum_col) || !(dig_element>sum_row)){
         throw runtime_error("ERRO: O ELEMENTO DÁ DIAGONAL DEVE SER MAIOR DO QUE A SOMA DOS ELEMENTOS DA COLUNA E DA LINHA");
      }
   }

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
      matriz_U[i]=mat[i];
   }
   
   for(int d = 0; d <= n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
      matriz_L[d][d] = 1;
      for(int i = d+1; i <= n; i++){ //seleciona a linha 


         float multiplicador = matriz_U[i][d]/matriz_U[d][d];
         matriz_L[i][d] = multiplicador;
         // Como L é a inversa do produto das matrizes, e sendo as matrizes M diagonais inferiores,
        // como o produto da inversa é a inversa do produto na ordem inversa, vale ir descobrindo M, e ir "colocando sua inversa" na matriz L
         for(int j = 0; j <= n; j++){ //anda na linha 
            matriz_U[i][j] = matriz_U[i][j] - multiplicador * matriz_U[d][j]; 
            }
      }
   }
   return matriz_U;

}

int main()
{   
   // vector<vector<float>> fds({ {3, 30, -1, 1}, 
   //                              {-1, 30, -1, 2}, 
   //                              {-1, 30, 3, 1} });
   // cout << fds.size() << "\n";
   // vector<vector<float>> aa({ {1,2,2,3},
   //                            {4,4,2,6},
   //                            {4,6,4,10}


   // });
   
   // cout << aa.size()-1;
   // Matrix m (fds);
   // Matrix mmmm(aa);
   // vector<float> sol;
   // sol =  m.sol_iter_jacobi();
   // for(int i = 0; i<sol.size();i++){
   //    cout << sol[i] << " ";
   // }
   // cout<< "\n";
   // cout<< "\n";
   // cout<< "\n";

   // vector<float> resolveai;
   // mmmm.escalona_gauss();

   // for(int i = 0; i <=2;i++){
   //    for(int j = 0 ; j<=3; j++){
   //    cout << mmmm.get_element(i,j) << "  ";
   //    }
   // cout << "\n"; 
   // }

   // resolveai = mmmm.eliminacao_gauss();

   // for(int i = 0; i<resolveai.size();i++){
   //    cout << resolveai[i] << " ";
   // }
   // cout<< "\n";

   // cout<< "\n";

   // cout<< "\n";

   // cout<< "\n";
 
   // vector<vector<float>> jeson({
   //                            {1,2,2,3},
   //                            {4,4,2,6},
   //                            {4,6,4,10}
   // });

   vector<vector<float>> gigi({
                              {-1000,2,3},
                              {1,1,1},
   });
   Matrix lk(gigi);
   lk.escalona_gauss();


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