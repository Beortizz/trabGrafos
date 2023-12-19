#include <iostream>
#include <vector>
#include <fstream>
#include <queue>
#include <stack>
#include <unordered_set>
#include <list>
#include <cstring>
#include <algorithm>
#include <limits>

#define INF numeric_limits<int>::max()

using namespace std;

void outputAdjacenciaToFile(const vector<vector<int>>& matrizAdjacencias, const string& filename) {
    ofstream outFile(filename);
    if (outFile.is_open()) {
        for (const auto& row : matrizAdjacencias) {
            for (int val : row) {
                outFile << val << " ";
            }
            outFile << "\n";
        }
        outFile.close();
    } else {
        cout << "Não foi possível abrir o arquivo.";
    }
}

void outputIncidenciaToFile(const vector<vector<int>>& matrizIncidencia, const string& filename) {
    ofstream outFile(filename);
    if (outFile.is_open()) {
        for (const auto& row : matrizIncidencia) {
            for (int val : row) {
                outFile << val << " ";
            }
            outFile << "\n";
        }
        outFile.close();
    } else {
        cout << "Não foi possível abrir o arquivo.";
    }
}

void outputTabelaIncidenciaToFile(const vector<vector<int>>& tabelaIncidencia, const string& filename) {
    ofstream outFile(filename);
    if (outFile.is_open()) {
        for (int i = 0; i < tabelaIncidencia.size(); ++i) {
            outFile << i << ": ";
            for (int j = 0; j < tabelaIncidencia[i].size(); ++j) {
                outFile << tabelaIncidencia[i][j] << " ";
            }
            outFile << "\n";
        }
        outFile.close();
    } else {
        cout << "Não foi possível abrir o arquivo.";
    }
}

vector<vector<int>> converterAdjacenciaParaTabelaIncidencia(const vector<vector<int>>& matrizAdjacencias) {
    int numVertices = matrizAdjacencias.size();
    int numArestas = 0;
    vector<vector<int>> tabelaIncidencia(numVertices);

    for (int i = 0; i < numVertices; ++i) {
        for (int j = i + 1; j < numVertices; ++j) {
            if (matrizAdjacencias[i][j] == 1) {
                tabelaIncidencia[i].push_back(j);
                tabelaIncidencia[j].push_back(i);
                numArestas++;
            }
        }
    }

    return tabelaIncidencia;
}

vector<vector<int>> converterMatrizIncidenciaParaTabelaIncidencia(const vector<vector<int>>& matrizIncidencia) {
    int numVertices = matrizIncidencia.size();
    int numArestas = matrizIncidencia[0].size();
    vector<vector<int>> tabelaIncidencia(numVertices);

    for (int j = 0; j < numArestas; ++j) {
        int firstVertex = -1, secondVertex = -1;
        for (int i = 0; i < numVertices; ++i) {
            if (matrizIncidencia[i][j] == 1) {
                if (firstVertex == -1) {
                    firstVertex = i;
                } else {
                    secondVertex = i;
                }
            }
        }
        if (firstVertex != -1 && secondVertex != -1) {
            tabelaIncidencia[firstVertex].push_back(secondVertex);
            tabelaIncidencia[secondVertex].push_back(firstVertex);
        }
    }
    return tabelaIncidencia;
}


vector<vector<int>> converterTabelaIncidenciaParaMatrizIncidencia(const vector<vector<int>>& tabelaIncidencia) {
    int numVertices = 0;
    for (const auto& verticesAdjacentes : tabelaIncidencia) {
        for (int v : verticesAdjacentes) {
            numVertices = max(numVertices, v + 1);
        }
    }

    int numArestas = tabelaIncidencia.size();
    vector<vector<int>> matrizIncidencia(numVertices, vector<int>(numArestas, 0));

    for (int j = 0; j < numArestas; ++j) {
        for (int i : tabelaIncidencia[j]) {
            matrizIncidencia[i][j] = 1;
        }
    }

    return matrizIncidencia;
}



vector<vector<int>> converterTabelaIncidenciaParaAdjacencia(const vector<vector<int>>& tabelaIncidencia) {
    int numVertices = tabelaIncidencia.size();
    vector<vector<int>> matrizAdjacencias(numVertices, vector<int>(numVertices, 0));

    for (int i = 0; i < numVertices; ++i) {
        for (int j : tabelaIncidencia[i]) {
            matrizAdjacencias[i][j] = 1;
        }
    }

    return matrizAdjacencias;
}



vector<vector<int>> converterAdjacenciaParaIncidencia(const vector<vector<int>>& matrizAdjacencias) {
    int numVertices = matrizAdjacencias.size();
    int numArestas = 0;

    // Conta o número de arestas no grafo
    for (int i = 0; i < numVertices; ++i) {
        for (int j = i + 1; j < numVertices; ++j) {
            if (matrizAdjacencias[i][j] == 1) {
                numArestas++;
            }
        }
    }

    vector<vector<int>> matrizIncidencia(numVertices, vector<int>(numArestas, 0));

    int coluna = 0;
    for (int j = 0; j < numVertices; ++j) {
        for (int i = j + 1; i < numVertices; ++i) {
            if (matrizAdjacencias[i][j] == 1) {
                matrizIncidencia[i][coluna] = 1;
                matrizIncidencia[j][coluna] = 1;
                coluna++;
            }
        }
    }

    return matrizIncidencia;
}

vector<vector<int>> converterIncidenciaParaAdjacencia(const vector<vector<int>>& matrizIncidencia) {
    int numVertices = matrizIncidencia.size();
    int numArestas = matrizIncidencia[0].size();

    vector<vector<int>> matrizAdjacencias(numVertices, vector<int>(numVertices, 0));

    for (int coluna = 0; coluna < numArestas; ++coluna) {
        int firstVertex = -1;
        int secondVertex = -1;

        for (int linha = 0; linha < numVertices; ++linha) {
            if (matrizIncidencia[linha][coluna] == 1) {
                if (firstVertex == -1)
                    firstVertex = linha;
                else
                    secondVertex = linha;
            }
        }

        if (firstVertex != -1 && secondVertex != -1) {
            matrizAdjacencias[firstVertex][secondVertex] = 1;
            matrizAdjacencias[secondVertex][firstVertex] = 1;
        }
    }

    return matrizAdjacencias;
}

void imprimirMatriz(const vector<vector<int>>& matriz) {
    for (int i = 0; i < matriz.size(); ++i) {
        for (int j = 0; j < matriz[i].size(); ++j) {
            cout << matriz[i][j] << " ";
        }
        cout << "\n";
    }
}

void buscaEmProfundidadeVisita(const vector<vector<int>>& grafo, int vertice, unordered_set<int>& visitados) {
    visitados.insert(vertice);
    cout << vertice << " ";

    int numVertices = grafo.size();
    for (int vizinho = 0; vizinho < numVertices; ++vizinho) {
        if (grafo[vertice][vizinho] != 0 && visitados.find(vizinho) == visitados.end()) {
            buscaEmProfundidadeVisita(grafo, vizinho, visitados);
        }
    }
}


void buscaEmProfundidade(const vector<vector<int>>& grafo, int verticeInicial) {
    int numVertices = grafo.size();

    if (verticeInicial >= numVertices || verticeInicial < 0) {
        cout << "Vértice inicial inválido." << endl;
        return;
    }

    unordered_set<int> visitados;
    buscaEmProfundidadeVisita(grafo, verticeInicial, visitados);
}

void buscaEmLargura(const vector<vector<int>>& grafo, int verticeInicial) {
    int numVertices = grafo.size();

    if (verticeInicial >= numVertices || verticeInicial < 0) {
        cout << "Vértice inicial inválido." << endl;
        return;
    }

    queue<int> fila;
    unordered_set<int> visitados;

    fila.push(verticeInicial);
    visitados.insert(verticeInicial);

    while (!fila.empty()) {
        int verticeAtual = fila.front();
        fila.pop();

        cout << verticeAtual << " ";

        for (int vizinho = 0; vizinho < numVertices; ++vizinho) {
            if (grafo[verticeAtual][vizinho] != 0 && visitados.find(vizinho) == visitados.end()) {
                visitados.insert(vizinho);
                fila.push(vizinho);
            }
        }
    }
}



int minKey(vector<int>& key, vector<bool>& mstSet) {
    int min = INT_MAX, min_index;

    for (int v = 0; v < key.size(); v++) {
        if (mstSet[v] == false && key[v] < min) {
            min = key[v];
            min_index = v;
        }
    }
    return min_index;
}

void printMST(vector<int>& parent, vector<vector<int>>& graph) {
    cout << "\n\nEdge \tWeight\n";
    for (int i = 1; i < graph.size(); i++)
        cout << parent[i] << " - " << i << "\t" << graph[i][parent[i]] << "\n";
}

void primMST(vector<vector<int>>& graph) {
    int V = graph.size();
    vector<int> parent(V);
    vector<int> key(V, INT_MAX);
    vector<bool> mstSet(V, false);

    key[0] = 0;
    parent[0] = -1;

    for (int count = 0; count < V - 1; count++) {
        int u = minKey(key, mstSet);
        mstSet[u] = true;

        for (int v = 0; v < V; v++) {
            if (graph[u][v] && mstSet[v] == false && graph[u][v] < key[v]) {
                parent[v] = u;
                key[v] = graph[u][v];
            }
        }
    }

    printMST(parent, graph);
}

int minDistance(vector<int>& dist, vector<bool>& visited) {
    int min = INT_MAX, min_index;

    for (int v = 0; v < dist.size(); v++) {
        if (visited[v] == false && dist[v] <= min) {
            min = dist[v];
            min_index = v;
        }
    }
    return min_index;
}

void dijkstra(const vector<vector<int>>& graph, int src) {
    int V = graph.size();
    vector<int> dist(V, INT_MAX);
    vector<bool> visited(V, false);

    dist[src] = 0;

    for (int count = 0; count < V - 1; ++count) {
        int u = minDistance(dist, visited); // Encontra o vértice de menor distância não visitado
        visited[u] = true;

        for (int v = 0; v < V; ++v) {
            if (!visited[v] && graph[u][v] != 0 &&
                dist[u] != INT_MAX && dist[u] + graph[u][v] < dist[v]) {
                dist[v] = dist[u] + graph[u][v];
            }
        }
    }

    // Mostra as distâncias mínimas a partir do vértice de origem
    cout << "Distancias minimas a partir do vertice " << src << ":\n";
    for (int i = 0; i < V; i++)
        cout << i << " \t" << dist[i] << "\n";
}



void dfs(int v, vector<vector<int>>& graph, vector<bool>& visited, stack<int>& stack) {
    visited[v] = true;

    for (int i = 0; i < graph.size(); ++i) {
        if (graph[v][i] == 1 && !visited[i]) {
            dfs(i, graph, visited, stack);
        }
    }

    stack.push(v);
}

vector<int> topologicalSort(vector<vector<int>>& graph) {
    int numVertices = graph.size();
    vector<bool> visited(numVertices, false);
    stack<int> stack;

    for (int v = 0; v < numVertices; ++v) {
        if (!visited[v]) {
            dfs(v, graph, visited, stack);
        }
    }

    vector<int> result;
    while (!stack.empty()) {
        result.push_back(stack.top());
        stack.pop();
    }

    return result;
}

void findEulerianCycleUtil(int u, vector<vector<int>>& graph, vector<int>& circuit) {
    for (int v = 0; v < graph.size(); ++v) {
        if (graph[u][v] > 0) {
            graph[u][v]--;
            graph[v][u]--;
            findEulerianCycleUtil(v, graph, circuit);
        }
    }
    circuit.push_back(u);
}

void findEulerianCycle(vector<vector<int>>& graph) {
    int vertices = graph.size();
    vector<int> degree(vertices, 0);
    vector<int> circuit;

    // Calculate degrees of vertices
    for (int i = 0; i < vertices; ++i) {
        for (int j = 0; j < vertices; ++j) {
            degree[i] += graph[i][j];
        }
    }

    // Find a vertex with odd degree as starting point
    int startVertex = 0;
    for (int i = 0; i < vertices; ++i) {
        if (degree[i] % 2 != 0) {
            startVertex = i;
            break;
        }
    }

    // Find Eulerian Cycle
    findEulerianCycleUtil(startVertex, graph, circuit);

    // Reverse the circuit to get the actual Eulerian cycle
    reverse(circuit.begin(), circuit.end());

    // Print Eulerian Cycle
    cout << "Eulerian Cycle: ";
    for (int i = 0; i < circuit.size(); ++i) {
        cout << circuit[i];
        if (i != circuit.size() - 1)
            cout << " -> ";
    }
    cout << endl;
}



int main() {
    // Exemplo de uma matriz de adjacências simples
    vector<vector<int>> matrizAdjacencias = {
            {0, 1, 0, 1, 0},
            {1, 0, 1, 0, 0},
            {0, 1, 0, 1, 1},
            {1, 0, 1, 0, 0},
            {0, 0, 1, 0, 0}
    };

    vector<vector<int>> grafo = {
            {0, 12, 3, 0},
            {12, 0, 1, 0},
            {3, 1, 0, 5},
            {0, 0, 5, 0}
    };

    vector<vector<int>> grafo2 = {
            {0, 1, 0, 0, 0, 0},
            {1, 0, 1, 0, 0, 0},
            {0, 0, 0, 1, 1, 0},
            {0, 0, 0, 0, 0, 1},
            {0, 0, 0, 0, 0, 1},
            {0, 0, 0, 0, 0, 0}
    };
    vector<vector<int>> grafo3 = {
            {0, 1, 1, 0, 0, 1, 0, 1, 0},
            {1, 0, 0, 0, 0, 0, 0, 1, 0},
            {1, 0, 0, 1, 0, 0, 1, 0, 1},
            {0, 0, 1, 0, 1, 1, 1, 0, 0},
            {0, 0, 0, 1, 0, 1, 0, 0, 0},
            {1, 0, 0, 1, 1, 0, 1, 0, 0},
            {0, 0, 1, 1, 0, 1, 0, 0, 1},
            {1, 1, 0, 0, 0, 0, 0, 0, 0},
            {0, 0, 1, 0, 0, 0, 1, 0, 0}

    };


    int verticeInicial;
    cout << "Digite o vértice inicial (de 0 a " << grafo.size() - 1 << "): ";
    cin >> verticeInicial;

    // Imprimir a matriz de adjacências original
    cout << "Matriz de Adjacências original:" << endl;
    imprimirMatriz(matrizAdjacencias);

    // Converter a matriz de adjacências para tabela de incidência
    vector<vector<int>> tabelaIncidenciaFromAdj = converterAdjacenciaParaTabelaIncidencia(matrizAdjacencias);
    // Imprimir a tabela de incidência obtida a partir da conversão de adjacências para tabela de incidência
    cout << "\nTabela de Incidência convertida de Adjacências:" << endl;
    for (int i = 0; i < tabelaIncidenciaFromAdj.size(); ++i) {
        cout << i << ": ";
        for (int j = 0; j < tabelaIncidenciaFromAdj[i].size(); ++j) {
            cout << tabelaIncidenciaFromAdj[i][j] << " ";
        }
        cout << endl;
    }

    // Converter a tabela de incidência para matriz de adjacências novamente
    vector<vector<int>> novaMatrizAdjFromTable = converterTabelaIncidenciaParaAdjacencia(tabelaIncidenciaFromAdj);

    // Imprimir a nova matriz de adjacências obtida a partir da conversão de tabela de incidência para adjacências
    cout << "\nNova Matriz de Adjacências convertida de Tabela de Incidência:" << endl;
    imprimirMatriz(novaMatrizAdjFromTable);

    // Converter a matriz de adjacências para matriz de incidência
    vector<vector<int>> matrizIncidenciaFromAdj = converterAdjacenciaParaIncidencia(matrizAdjacencias);

    // Imprimir a matriz de incidência obtida a partir da conversão de adjacências para matriz de incidência
    cout << "\nMatriz de Incidência convertida de Adjacências:" << endl;
    imprimirMatriz(matrizIncidenciaFromAdj);

    // Converter a matriz de incidência para tabela de incidência
    vector<vector<int>> tabelaIncidenciaFromInc = converterMatrizIncidenciaParaTabelaIncidencia(matrizIncidenciaFromAdj);

    // Imprimir a tabela de incidência obtida a partir da conversão de matriz de incidência para tabela de incidência
    cout << "\nTabela de Incidência convertida de Matriz de Incidência:" << endl;
    for (int i = 0; i < tabelaIncidenciaFromInc.size(); ++i) {
        cout << i << ": ";
        for (int j = 0; j < tabelaIncidenciaFromInc[i].size(); ++j) {
            cout << tabelaIncidenciaFromInc[i][j] << " ";
        }
        cout << endl;
    }

    // Converter a tabela de incidência para matriz de incidência
    vector<vector<int>> novaMatrizIncidenciaFromTable = converterTabelaIncidenciaParaMatrizIncidencia(tabelaIncidenciaFromInc);

    // Imprimir a nova matriz de incidência obtida a partir da conversão de tabela de incidência para matriz de incidência
    cout << "\nNova Matriz de Incidência convertida de Tabela de Incidência:" << endl;
    imprimirMatriz(novaMatrizIncidenciaFromTable);

    // Converter a matriz de incidência para matriz de adjacências
    vector<vector<int>> novaMatrizAdjFromInc = converterIncidenciaParaAdjacencia(matrizIncidenciaFromAdj);

    // Imprimir a nova matriz de adjacências obtida a partir da conversão de matriz de incidência para adjacências
    cout << "\nNova Matriz de Adjacências convertida de Matriz de Incidência:" << endl;
    imprimirMatriz(novaMatrizAdjFromInc);

    cout << "\n Nova Matriz de Incidencia convertida de Matriz de Adjacencia:" << endl;
    imprimirMatriz(converterAdjacenciaParaIncidencia(matrizAdjacencias));

    outputAdjacenciaToFile(matrizAdjacencias, "matriz_adjacencia.txt");
    outputIncidenciaToFile(matrizIncidenciaFromAdj, "matriz_incidencia.txt");
    outputTabelaIncidenciaToFile(tabelaIncidenciaFromAdj, "tabela_incidencia.txt");

    cout << "\nBusca em largura a partir do vértice " << verticeInicial << ": ";
    buscaEmLargura(grafo, verticeInicial);

    cout << "\nBusca em Profundidade a partir do vértice " << verticeInicial << ": ";
    buscaEmProfundidade(grafo, verticeInicial);

    primMST(grafo);
    cout << "Agora Vertice inicial para dijkstra: ";
    cin >> verticeInicial;
    dijkstra(grafo, verticeInicial);
    vector<int> result = topologicalSort(grafo2);

    cout << "Ordenação topológica: ";
    for (int vertex : result) {
        cout << vertex << " ";
    }
    cout << endl;

    findEulerianCycle(grafo3);

    return 0;
}


