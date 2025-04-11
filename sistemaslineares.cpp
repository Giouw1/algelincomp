#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <cmath>
#include <bits/stdc++.h>
using namespace std;
typedef vector<vector<double>> matriz;

matriz trata_matriz(matriz m, int i){
    int n = m.size()-1; // índice do final da matriz

    for(int j = i + 1; j<n; j++){
       if(m[j][i]!=0){ //achamos uma linha com pivo não nulo
          swap(m[i], m[j]);
          return m; 
       }
    }
    //se não achar printar erro porque a matriz é LD e temos infinitas solucoes?
    cout << "A matiz é LD" << "\n";
    return m;
 }
 
 matriz get_matriz_A(const vector<vector<double>>& mat_completa) {
    int linhas = mat_completa.size();
    int colunas = mat_completa[0].size() - 1; // Ignora a última coluna(b)
    vector<vector<double>> A(linhas, vector<double>(colunas));

    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            A[i][j] = mat_completa[i][j];
        }
    }

    return A;
}
 int det_laplace(const vector<vector<double>>& submat) {
    int tamanho = submat.size();

    if (tamanho == 1) {
        return submat[0][0];
    }

    if (tamanho == 2) {
        return submat[0][0] * submat[1][1] - submat[0][1] * submat[1][0];
    }

    double determinante = 0.0;

    for (int k = 0; k < tamanho; ++k) {
        vector<vector<double>> minor;

        for (int i = 1; i < tamanho; ++i) {
            vector<double> linha;
            for (int j = 0; j < tamanho; ++j) {
                if (j != k) {
                    linha.push_back(submat[i][j]);
                }
            }
            minor.push_back(linha);
        }

        double cofator = pow(-1, k) * submat[0][k] * det_laplace(minor);
        determinante += cofator;
    }

    return determinante;
}
vector<double> eliminacao_baixo_cima(matriz m){
   int n = m.size() - 1;
   vector<double> x(n+1,0);
   // sabemos o último valor de x, o x_n
   x[n] = m[n][n+1]/m[n][n];

   for(int i = n-1; i>=0; i--){
      double sum = 0;
      for(int j = i+1; j<= n; j++){
         sum += m[i][j]*x[j];
      }
      x[i] = (m[i][n+1] - sum)/m[i][i];
   }
   return x;
}
vector<double> escalona_gauss(matriz m){
    int n = m.size()-1; // índice do final da matriz

    // primeiro checamos se é possível escalonar a Matriz
    int determinante = det_laplace(get_matriz_A(m));
    if (determinante<0){
       cout << "Como o detemrinante é negativo não conseguimos encontrar solução" << "\n";
       return {};
    }
    auto start = high_resolution_clock::now(); // marca o início
    for(int d = 0; d < n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
       for(int i = d+1; i <= n; i++){ //seleciona a linha 
          if (m[d][d]==0){ //achamos uma linha com pivo zero
             trata_matriz(m,d); 
          }
          double multiplicador = m[i][d]/m[d][d];
          for(int j = 0; j <= n+1; j++){ //anda na linha 
             m[i][j] = m[i][j] - multiplicador * m[d][j]; 
             }
       }
       auto end = high_resolution_clock::now(); // marca o fim
       auto duration = duration_cast<seconds>(end - start);
       string filename = "resultados/resultados_escalona_gauss" + to_string(m.size()) + ".txt";

       ofstream tempoescalona(filename);
       if(tempoescalona.is_open()){
         cout<< "arquivo criado com sucesso";
         cout <<"\n";
      }
      else{cerr << "Erro criando o arquivo " << endl;}
      tempoescalona<< "Tempo:"<< duration.count() <<"s"<< " "<< m.size()<<endl;
    }
    return eliminacao_baixo_cima(m);
}

 vector<double> sol_iter_jacobi(matriz m, double tol = 0.000001, int numero_iteracoes =1000){
    int n = m.size()-1; 
    // antes de iniciar o algoritmo devemos checar se cada elemento da diagonal é o maior elemento da sua linha e coluna correspondente
    for(int i = 0; i<n; i++){
       double dig_element = fabsf(m[i][i]);
       double sum_col = 0;
       double sum_row = 0;
       for(int j = 0; j<n; j++){
          if(j != i){
             sum_col+= fabsf(m[j][i]);
             sum_row+= fabs(m[i][j]);
          }
       }
       if(!(dig_element>sum_col) || !(dig_element>sum_row)){
          throw runtime_error("ERRO: O ELEMENTO DÁ DIAGONAL DEVE SER MAIOR DO QUE A SOMA DOS ELEMENTOS DA COLUNA E DA LINHA");
       }
    }
 
    vector<double> x_old(n+1, 1.0); //iniciamos o vetor solucao inicial com todos os elementos como 1
    double r = 1;
    int iteracoes = 1;
 
    string filename = "resultados/resultados_sol_iterativa_dim_" + to_string(n+1) + ".txt";
    ofstream resultsFile(filename);
    if(resultsFile.is_open()){
       cout<< "arquivo criado com sucesso";
       cout <<"\n";
    }
    else{cerr << "Erro criando o arquivo " << endl;}
    //loop para descobrir a solução de fato
    while((iteracoes<= numero_iteracoes)&&(r>tol)){
       vector<double> x_new(n+1,1);
       for(int i = 0; i<= n; i++){
          double sum = 0;
          for(int j = 0; j <= n; j++){
             if(i!=j){
                sum+= m[i][j]*x_old[j];
             }
          }
          x_new[i] = (m[i][n+1] - sum)/m[i][i];
       }
 
       // calculamos o r para comparar com a tolerância
       double up_sum_square = 0;
       double down_sum_square = 0;
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
 vector<double> eliminacao_cima_baixo_LU(matriz m, vector<double> b){
   int n = m.size() - 1;
   vector<double> x(n+1,0);
   // sabemos o primeiro valor
   x[0] = b[0]/m[0][0];
   for(int i = 1; i<=n; i++){
      double sum = 0.0;
      for(int j = 0; j < i; j++){
         sum += m[i][j]*x[j];
      }
      x[i] = (b[i] - sum)/m[i][i];
   }
   return x;
}
vector<double> eliminacao_baixo_cima_LU(matriz m,vector<double> y){
   int n = m.size() - 1;
   vector<double> x(n+1,0);
   // sabemos o último valor de x, o x_n
   x[n] = y[n]/m[n][n];

   for(int i = n-1; i>=0; i--){
      double sum = 0.0;
      for(int j = i+1; j<= n; j++){
         sum += m[i][j]*x[j];
      }
      x[i] = (y[i] - sum)/m[i][i];
   }
   return x;
}
 vector<double> fatoracao_lu(matriz m){
    int n = m.size()-1; //índice do final da matriz
    vector<vector<double>> matriz_L(n+1,vector<double>(n+1));
    vector<vector<double>> matriz_U(n+1,vector<double>(n+1));
    vector<double> b(n+1);
    vector<double> y(n+1);
    vector<double> x(n+1);
    for(int i = 0; i<=n;i++){
        for (int j = 0;j<=n;j++){
            matriz_U[i][j]=m[i][j];
            }
        }
    for (int i = 0; i<=n;i++){
        b[i] = m[i][n+1];
    }
    for(int d = 0; d < n; d++){ //quantidade de diagonais, ele não precisa escalonar embaixo do último pivo porque não tem ninguém para baixo
        matriz_L[d][d] = 1;
        for(int i = d+1; i <= n; i++){ //seleciona a linha 
            double multiplicador = matriz_U[i][d]/matriz_U[d][d];
            matriz_L[i][d] = multiplicador;
            // Como L é a inversa do produto das matrizes, e sendo as matrizes M diagonais inferiores,
            // como o produto da inversa é a inversa do produto na ordem inversa, vale ir descobrindo M, e ir "colocando sua inversa" na matriz L
            for(int j = 0; j <= n; j++){ //anda na linha 
                matriz_U[i][j] = matriz_U[i][j] - multiplicador * matriz_U[d][j]; 
            }
        }
    }
    matriz_L[n][n] = 1;

    y = eliminacao_cima_baixo_LU(matriz_L, b); // conseguir o vetor y, e fazer UX = Y
    return x = eliminacao_baixo_cima_LU(matriz_U, y);
    
}
 vector<vector<double>> gerar_matriz_jacobi_compat(int n) {
    srand(time(0));
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    vector<double> b(n);
   
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
    vector<vector<double>> matriz_completa(n, vector<double>(n + 1));
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
int main(){
    vector<vector<double>> mat = {{1,2,2,3}, {4,4,2,6},{4,6,4,10}};
    mat = gerar_matriz_jacobi_compat(6);
    vector<double> lu = fatoracao_lu(mat);

    for (int i =0;i<=lu.size()-1;i++){
        cout << lu[i]<<endl;
    }
    cout <<endl;
    vector<double> escalonaga = escalona_gauss(mat);
    for (int i =0;i<=escalonaga.size()-1;i++){
        cout << escalonaga[i]<<endl;
    }
    cout <<endl;
    vector<double> jaco =sol_iter_jacobi(mat);
    for (int i =0;i<=jaco.size()-1;i++){
        cout << jaco[i]<<endl;
    }
    cout <<endl;    
}


 

 