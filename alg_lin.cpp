#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

// #include <bits/stdc++.h>

using namespace std;

class Matriz                      
{
  public:                     
    Matriz(vector<vector<float> > input);       
      void trata_matriz(int i); // Tenta trocar a linha atual com uma abaixo que tenha pivô diferente de zero. Necessário para evitar divisões por zero durante a eliminação de Gauss
      float get_element(int i, int j){return mat[i][j];};
      vector<vector<float>> get_matriz_A(const vector<vector<float>>& mat_completa); //Obtem a matriz sem o lado direito para aplicar laplace. Retorna apenas a parte A da matriz aumentada [A | b], ignorando o vetor b
      int det_laplace(const vector<vector<float>>& submat); // Calcula o determinante de uma matriz de forma recursiva usando a expansão de Laplace. Método ineficiente para matrizes grandes, mas serve para fins didáticos

      void escalona_gauss(); // Aplica a eliminação de Gauss para triangularizar a matriz. Checa primeiro se o determinante é negativo como condição simplificada de existência de solução
      vector<float> eliminacao_gauss(); // Resolve o sistema linear Ax = b usando eliminação de Gauss. Supõe que a matriz já está em forma escalonada triangular superior

      vector<vector<float>> fatoracao_lu();

      vector<float> sol_iter_jacobi(float tol = 0.000001, int numero_iteracoes = 1000);
    
 private:                      
    int n;                // n é o índice da útlima linha
    vector<vector<float> > mat;
    int dim;
    int det;
};


Matriz::Matriz(vector<vector<float> > input)
{
   mat = input;
   n = mat.size()-1; //nós fazemos -1 porque no lado direito da Matriz temos o vetor de coeficientes b
}

void Matriz::trata_matriz(int i){

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

vector<vector<float>> Matriz::get_matriz_A(const vector<vector<float>>& mat_completa) {
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

int Matriz::det_laplace(const vector<vector<float>>& submat) {
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



void Matriz::escalona_gauss(){

   // primeiro checamos se é possível escalonar a Matriz
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

vector<float> Matriz::eliminacao_gauss(){

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

vector<float> Matriz::sol_iter_jacobi(float tol, int numero_iteracoes){
   
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


   string filename = "resultados/resultados_sol_iterativa_dim_" + to_string(n+1) + ".txt";
   ofstream resultsFile(filename);
   if(resultsFile.is_open()){
      cout<< "arquivo criado com sucesso";
      cout <<"\n";
   }
   else{cerr << "Erro criando o arquivo " << endl;}

   while((iteracoes<= numero_iteracoes)&&(r>tol)){
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
      resultsFile << "Numero de iterações: " << iteracoes << "  Norma do resíduo: " << r << endl;
      x_old = x_new; // agora usamos o mais novo como sendo o velho para a próxima iteração
      iteracoes++;
   }
   // }
   resultsFile.close();

   // ná última iteração do algoritmo o x_old é o mais recente vetor de x, é a solução do sistema

   return x_old;
}

vector<vector<float>>Matriz::fatoracao_lu(){
   //
   // vector<vector<float>> matriz_L(n+1,vector<float>(n+1));
   // vector<vector<float>> matriz_U(n+1,vector<float>(n+1));
   // vector<float> b(n+1);
   // b = mat[n];

   // for(int i = 0; i<=n;i++){
   //    matriz_U[i]=mat[i];
   // }
   
   
   // for(int d = 0; d <= n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
   //    matriz_L[d][d] = 1;
   //    for(int i = d+1; i <= n; i++){ //seleciona a linha 


   //       float multiplicador = matriz_U[i][d]/matriz_U[d][d];
   //       matriz_L[i][d] = multiplicador;
   //       // Como L é a inversa do produto das matrizes, e sendo as matrizes M diagonais inferiores,
   //      // como o produto da inversa é a inversa do produto na ordem inversa, vale ir descobrindo M, e ir "colocando sua inversa" na matriz L
   //       for(int j = 0; j <= n; j++){ //anda na linha 
   //          matriz_U[i][j] = matriz_U[i][j] - multiplicador * matriz_U[d][j]; 
   //          }
   //    }
   // }

   vector<vector<float>> matriz_L(n+1,vector<float>(n+1));
   vector<vector<float>> matriz_U(n+1,vector<float>(n+1));
   vector<float> b(n+1);
   for(int i = 0; i<=n;i++){
    for (int j = 0;j<=n;j++){
        matriz_U[i][j]=mat[i][j];
      }
   }
   for (int i = 0; i<=n;i++){
      b[i] = mat[i][n+1];
   }

   for(int d = 0; d < n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
      matriz_L[d][d] = 1;
      for(int i = d+1; i <= n; i++){ //seleciona a linha 


         float multiplicador = matriz_U[i][d]/matriz_U[d][d];
         matriz_L[i][d] = multiplicador;
         // Como L é a inversa do produto das matrizes, e sendo as matrizes M diagonais inferiores,
        // como o produto da inversa é a inversa do produto na ordem inversa, vale ir descobrindo M, e ir "colocando sua inversa" na matriz L
         for(int j = 0; j < n; j++){ //anda na linha 
            matriz_U[i][j] = matriz_U[i][j] - multiplicador * matriz_U[d][j]; 
            }
      }
   }

   vector<float> x(n+1,0);
   vector<float> y(n+1,0);
   // sabemos o último valor de x, o x_n
   y[0] = b[0]/matriz_L[0][0];
   for(int i = 1; i<=n; i++){
      float sum = 0;
      for(int j = i+1; j<= n; j++){
         sum += matriz_L[i][j]*y[j];
      }
      y[i] = (matriz_L[i][n+1] - sum)/matriz_L[i][i];
   }

   x[n] = y[n]/matriz_U[n][n];
   for(int i = n-1; i>=0; i--){
      float sum = 0;
      for(int j = i+1; j<= n; j++){
         sum += matriz_U[i][j]*x[j];
      }
      x[i] = (matriz_U[i][n+1] - sum)/matriz_U[i][i];
   }

   for(int i = 0; i<=n;i++){
      for(int j = 0; j <= n; j++){
         cout << matriz_U[i][j] << " ";
      }
      cout << endl;
   }
   return matriz_U;

}

vector<vector<float>> gerar_matriz_jacobi_compat(int n) {
    srand(time(0));
    vector<vector<float>> A(n, vector<float>(n, 0.0));
    vector<float> b(n);
   
   if(n==1){
      A = {{4}};
   }
   else if(n==2){
      A = {{4,-1},
         {-1,4}};
   }
   else{
   for(int i = 0 ; i < n; i++){
      A[i][i] = 4;
   }
   for(int i = 1 ; i < n; i++){
      if(i%2 != 0){
         for(int j = 0; j < n-i; j++){
            A[j][j+i] = -1;
            A[j+i][j] = -1;
         }
      }
      else{
         for(int j = 0; j < n-i; j++){
            A[j][j+i] = 0;
            A[j+i][j] = 0;
         }
      }
   }
   }

    // Gera vetor b aleatório (valores de 1 a 10)
    for (int i = 0; i < n; ++i) {
        b[i] = rand() % 10 + 1;
    }

    // Cria matriz aumentada [A | b]
    vector<vector<float>> matriz_completa(n, vector<float>(n + 1));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            matriz_completa[i][j] = A[i][j];
        }
        matriz_completa[i][n] = b[i];
    }

   for(int i = 0; i<n;i++){
      for(int j = 0; j<= n; j++){
         cout << matriz_completa[i][j] << " ";
      }
      cout<< endl;
   }
   cout << endl;
    return matriz_completa;
}



int main()
{   
   // geramos as matrizes aleatórias
   vector<vector<float>> matriz3 = gerar_matriz_jacobi_compat(3);
   vector<vector<float>> matriz5 = gerar_matriz_jacobi_compat(5);
   vector<vector<float>> matriz7 = gerar_matriz_jacobi_compat(7);

   Matriz m3(matriz3);
   Matriz m5(matriz5);
   Matriz m7(matriz7);
   Matriz m({{},{},{}});

   
   // m3.sol_iter_jacobi();
   // m5.sol_iter_jacobi();
   // m7.sol_iter_jacobi();
   m3.fatoracao_lu();

   return 0;
}  