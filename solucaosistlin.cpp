#include <iostream>
#include <bits/stdc++.h>
using namespace std;
class matrix{
    private:
        vector<vector<float>> matriz; // Matriz de adjacência
        int n = 0;
    
    void construtor(){

    }
    public:
    matrix(vector<vector<float>> mat){
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
        float multiplicador;
        for (int i =0; i<=n-2;i++){//coluna a coluna, eliminar os de baixo da diagonal
            for (int j = i+1;j<=n-1;j++){//para baixo da diagonal, eliminar
                if (matriz[i][i] == 0){
                    trata_matriz(i);
                }
                multiplicador = matriz[j][i]/matriz[i][i];//elemento menos "pivo" "vezes elemento" divido por pivo
                //matriz[j][i] é o elemento abaixo da diagonal/pivo, que é o matriz[i][i]
                //O loop com k faz a subtração para toda a linha em questão
                //Funciona pois, como a subtração é feita a partir do pivô, onde 
                //"já é resolvido" fica uma subtração de zero por zero:
                //não atrapalha o que já foi resolvido
                for (int k = 0;k<=n;k++){
                    matriz[j][k]-=matriz[i][k]*multiplicador;
                }

            }
        }
    }


    vector<float> solucaoiterativa(){
        int d = 1;
        for (int i = 0; i<= n-1;i++){
            d*= matriz[i][i];//Constante d
        }
        if (d == 0){
            throw("Erro, diagonal nula");//Divisão por zero: Usar o formato 2
        }

        vector<float> solucao(n-1, 0);//Solução anterior
        vector<float> solucaof(n-1,0);//Próxima solução
        vector<float> b = matriz[n];//Lado direito
        while (true){
            vector<float> vetoraux = multmatrizvetor(1,solucao); //pega o vetor pós multiplicação por L+U
            for (int j = 0;j<=n-1;j++){ //Faz a parte de cima da equação
                solucaof[j] = (b[j] - vetoraux[j])/d;
            }
            //Verifica se o resto já bate
            if(modulo(solucaof) - modulo(solucao) <= 0.01){
                return solucaof;
            }
            //A nova solução antiga é a antiga solução nova
            //E a solução nova vai ser redefenida
            solucao = solucaof;
        }
    }
    vector<float> multmatrizvetor/*para matriz matriz é só fazer n vezes*/
                                (bool d, /*se tiver 1, tem a diagonal, se não, ela é 0*/
                                vector<float> vetor){
        vector<float>vetorf(n-1,0);
        for (int i = 0;i<=n-1;i++){//para cada coluna da matriz
            for (int j = 0;j< matriz[i].size();j++){//para cada elemento da coluna
                if (d == 1 and i == j){//se for a diagonal, é soma por zero, não faz diferença
                    continue;
                }
                vetorf[j] += matriz[i][j]*vetor[i];//se não, multiplica o elemento, e soma à sua posição de linha correta no vetor final
            }
        }
        return vetorf;
    }
    float modulo(vector<float> vetor){
        float mod;
        for (int i = 0; i<= vetor.size();i++){
            mod = vetor[i]*vetor[i];
        }
        return sqrt(mod);
    }

}
void main(){
}
