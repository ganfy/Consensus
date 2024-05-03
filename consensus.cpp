#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <fstream>
#include <queue>
#include <stack> 
#include <bitset>

using namespace std;

string ReverseComplement(string s) {
    map<char, char> complements ={
        {'A', 'T'},
        {'T', 'A'},
        {'C', 'G'},
        {'G', 'C'}
    };
    for (char& nucleotide : s) {
        nucleotide = complements[nucleotide];
    }
    reverse(s.begin(), s.end());
    return s;
}

int FindOverlap(const string& seq1, const string& seq2) {
    int overlap = 0;
    int len = min(seq1.size(), seq2.size());
    for (int i = 0; i < len; ++i) {
        if (seq1.substr(seq1.size() - len + i) == seq2.substr(0, len - i)) {
            overlap = len - i;
            break;
        }
    }
    return overlap;
}

vector<pair<int, string>> GreedyHamiltonianPath(vector<string>& sequences) {

    vector<bool> visited(sequences.size(), false);

    vector<pair<int, string>> hamiltonianPath;

    vector<vector<int>> overlapMatrix(sequences.size(), vector<int>(sequences.size(), 0));
    for (int i = 0; i < sequences.size(); ++i) {
        for (int j = 0; j < sequences.size(); ++j) {
            if (i != j) {
                overlapMatrix[i][j] = FindOverlap(sequences[i], sequences[j]);
            }
        }
    }

    int maxOverlap = 0;
    int startVertex = -1;
    for (int i = 0; i < sequences.size(); ++i) {
        int currentOverlap = 0;
        for (int j = 0; j < sequences.size(); ++j) {
            if (i != j && overlapMatrix[i][j] > currentOverlap) {
                currentOverlap = overlapMatrix[i][j];
            }
        }
        if (currentOverlap > maxOverlap) {
            maxOverlap = currentOverlap;
            startVertex = i;
        }
    }

    if (startVertex == -1) {
        startVertex = 0;
    }

    hamiltonianPath.push_back({0, sequences[startVertex]});
    visited[startVertex] = true;

    for (int i = 0; i < sequences.size() - 1; ++i) {
        maxOverlap = 0;
        int nextVertex = -1;

        for (int j = 0; j < sequences.size(); ++j) {
            if (!visited[j]) {
                string currentSequence = sequences[j];
                string revComplement = ReverseComplement(currentSequence);
                int overlap = max(FindOverlap(hamiltonianPath.back().second, currentSequence), FindOverlap(hamiltonianPath.back().second, revComplement));
                if (overlap > maxOverlap) {
                    maxOverlap = overlap;
                    nextVertex = j;
                }
            }
        }
      
        if (nextVertex != -1) {
            string nextSequence = sequences[nextVertex];
            string revComplement = ReverseComplement(nextSequence);
            if (FindOverlap(hamiltonianPath.back().second, nextSequence) >= FindOverlap(hamiltonianPath.back().second, revComplement)) {
                hamiltonianPath.push_back({maxOverlap, nextSequence});
            } else {
                hamiltonianPath.push_back({maxOverlap, revComplement});
            }
            visited[nextVertex] = true;
        }
    }

    return hamiltonianPath;
}

vector<string> ReadSequences(const string& filename) {
  ifstream file(filename);
  if (!file.is_open()) {
    cerr << "Error: Unable to open file " << filename << endl;
    exit(1);
  }

  vector<string> sequences;
  string line, sequence, sequence2;

  while (getline(file, line)) {
    if (line.empty()) continue;
    sequences.push_back(line);
  }

  file.close();
  return sequences;
}

void VisualizeGraph_(const vector<vector<int>>& overlap_graph, string name) {
    ofstream dotFile(name + ".dot");
    if (!dotFile) {
        cerr << "Error: Unable to create DOT file." << endl;
        return;
    }

    dotFile << "digraph WeightedGraph {" << endl;
    int n = overlap_graph.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (overlap_graph[i][j] > 0) {
                dotFile << "    " << i << " -> " << j << " [label=\"" << overlap_graph[i][j] << "\"];" << endl;
            }
        }
    }
    dotFile << "}" << endl;

    dotFile.close();
}

string Hamiltonian_(const vector<string>& seqs_og, int t=0, bool consider_reverse = false) {
    vector<string> seqs = seqs_og;
    int n = seqs.size();
    vector<vector<int>> overlap_graph(n, vector<int>(n, 0));
    string consensus; 
    int overlap;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                overlap = FindOverlap(seqs[i], seqs[j]);
                if (overlap >= t)
                overlap_graph[i][j] = overlap;
            }
        }
    }

    VisualizeGraph_(overlap_graph, "acyclic");
    vector<int> order(n);
    for (int i = 0; i < n; ++i) {
        order[i] = i;
    }

    do {
        bool valid = true;
        string alignment;
      int pad = 0;
        for (int i = 0; i < n; ++i) {
            if (i > 0) {
                int overlap = overlap_graph[order[i - 1]][order[i]];
                if (overlap < t) {
                    valid = false;
                    break;
                }
                pad += seqs[order[i-1]].size() - overlap;
                alignment += string(pad, '-') + seqs[order[i]] + "\n";
                consensus += seqs[order[i]].substr(overlap);
            } else {
                alignment += seqs[order[i]] + "\n";
                consensus += seqs[order[i]];
            }
        }
        if (valid) {
            cout << alignment;
            return consensus;
        }
    } while (next_permutation(order.begin(), order.end()));

    return "No se encontró un camino hamiltoniano válido.";
}


int main() {
   vector<string> sequences = ReadSequences("B3.txt");
    vector<pair<int, string>> hamiltonianPath = GreedyHamiltonianPath(sequences);

    string alignment;
    int s_ant = 0;
  int pad = 0;
    for (const auto& pair : hamiltonianPath) {
      pad += s_ant - pair.first;
      alignment += string(pad, '-') + pair.second + '\n';
      s_ant = pair.second.size();
      cout << pad << " ";
        
    }

    cout << "Camino Hamiltoniano:" << endl;
    for (const auto& pair : hamiltonianPath) {
        cout << "Overlap: " << pair.first << ", Secuencia visitada: " << pair.second << endl;
    }

    cout << "Alineamiento:" << endl;
    cout << alignment << endl;

    return 0;
}
