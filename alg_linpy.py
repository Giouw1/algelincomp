#Elimnação de Gauss
import copy

def laplace(matriz):
    if len(matriz)==1:
        return matriz[0][0]
    elif len(matriz)==2: 
        return matriz[0][0]*matriz[1][1]-matriz[0][1]*matriz[1][0]
    else:
        linha = matriz[0]
        matriz_aux = []
        cofatores = [] 
        determinante = 0
        for i in range(len(matriz)):
            matriz_aux = matriz[1:]
            for j in range(len(matriz_aux)):
                matriz_aux[j] = matriz_aux[j][:i] + matriz_aux[j][i+1:]

            cofator = ((-1)**i)*laplace(matriz_aux) 
            cofatores.append(cofator)

        for k in range(len(cofatores)):
            determinante += linha[k]*cofatores[k]

        return determinante
    

def gauss(matriz, solucao): 
    m1 = copy.deepcopy(matriz)
    solucao1= copy.deepcopy(solucao)
    resultados=[]
    for r in range(len(matriz)):
        resultados.append(None)
    for a in range(len(m1)): #Linha
        for b in range(len(m1)): #Coluna
            if a>b:
                valor = m1[a][b]/m1[b][b]
                solucao1[a] = solucao1[a] - valor*solucao1[b]

                for c in range(len(m1)): #Percorrer a Linha
                    m1[a][c]= m1[a][c] - valor*m1[b][c]

    for d in range(len(m1)-1, -1, -1): #Linha
        equacao=0
        for e in range(len(m1)-1, -1, -1): #Coluna
            if e>=d:
                if resultados[e]!=None:
                    equacao+=m1[d][e]*resultados[e]
                else:
                    resultados[e] = (solucao1[e] - equacao)/m1[d][e]

    return m1, resultados

def LU(matriz, solucao):
    L = []
    U = []
    X = [] 
    Y = []
    B = copy.deepcopy(solucao)

    for w in range(len(matriz)):
        L.append([0]*len(matriz))
        U.append([0]*len(matriz))
        L[w][w]= 1
        X.append(None)
        Y.append(None)


    for t in range(len(matriz)): #Linha
        for r in range(len(matriz)): #Coluna
            soma = 0
            for s in range(t):
                soma += L[t][s]*U[s][r]
            if t<=r:
                U[t][r]= matriz[t][r] - soma
            elif t>r:
                L[t][r] = (matriz[t][r] - soma)/U[r][r]


    for x in range(len(matriz)): #Linha
        equacao1=0
        for z in range(len(matriz)): #Coluna
            if z<=x:
                if Y[z]!=None:
                    equacao1+=L[x][z]*Y[z]
                else:
                    Y[z] = (B[z] - equacao1)/L[x][z]


    for d in range(len(matriz)-1, -1, -1): #Linha
        equacao2=0
        for e in range(len(matriz)-1, -1, -1): #Coluna
            if e>=d:
                if X[e]!=None:
                    equacao2+=U[d][e]*X[e]
                else:
                    X[e] = (Y[e] - equacao2)/U[d][e]


    print(L, U, Y, X)
    return matriz, solucao

def Jacobi(matriz, solucao):
    return matriz, solucao

    

matriz=[]
solucao=[]
m = int(input("Número de linhas:"))


for x in range(m):
    matriz.append(list(map(int, input(f"Informe a linha {x+1}:").split())))


solucao = list(map(int, input("Informe a solução:").split()))

n = len(matriz[0])


if m!=n:
    print("Não é uma matriz quadrada.")

else:
    determinante=laplace(matriz)

    if determinante == 0:
        print("Determinante igual à zero")
    
    else:
        #print(gauss(matriz, solucao))
        LU(matriz, solucao)

        