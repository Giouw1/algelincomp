#include <iostream>
#include <bits/stdc++.h>
using namespace std;
class matrix{
    private:
        vector<vector<int>> matriz; // Matriz de adjacência
        int n = 0;
    
    void construtor(){

    }
    public:
    matrix(vector<vector<int>> mat){
        matriz = mat;
        n = mat.size()-1; //índice do "lado direito" da matriz. a matriz tem índices n-1 por n-1
    }
    void trata_matriz(int i) {
        for (int j = i + 1; j < n; j++) {
            if (matriz[j][i] != 0) { // encontra uma linha com um pivô diferente de 0
                swap(matriz[i], matriz[j]); // troca as linhas
                return;
            }
        }
    }
    void escalona_gauss(){
        int multiplicador;
        for (int i =0; i<=n-2;i++){//não preciso escalonar a última coluna.
            for (int j = i+1;j<=n-1;j++){
                if (matriz[i][i] == 0){
                    trata_matriz(i);
                }
                multiplicador = matriz[j][i]/matriz[i][i];
                
                
                
                for (int k = 0;k<=n;k++){
                    matriz[j][k]-=matriz[i][k]*multiplicador;
                }

            }
        }
    }
    void solucaoiterativa(){
        int d = 1;
        for (int i = 0; i<= n-1;i++){
            d*= matriz[i][i];
        }
        if (d == 0){
            throw("Erro, diagonal nula")
        }
        vector<int> solucao(n-1, 0);
        vector<int> b = matriz[n];
    }
    void multmatriz(bool d /*se tiver 1, tem a diagonal, se não, ela é 0*/){

    }

}
int main(){
    return 1;
}
