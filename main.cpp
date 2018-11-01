#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <tuple>
#include <fstream>
#include <chrono>

using namespace std;

typedef pair<string, float> bet_pair;
typedef tuple<vector<string>, map<string, vector<string>>, map<string, float>> param_tuple;

vector<string> extract_keys(map<string, vector<string>> input_map) {
    vector<string> keys;
    keys.reserve(input_map.size());
    for (auto const& elem : input_map) {
        keys.push_back(elem.first);
    }

    return keys;
}

void print_top_five(map<string, float> bet_map) {
    vector<bet_pair> bet_vec;
    copy(bet_map.begin(), bet_map.end(), back_inserter<vector<bet_pair>>(bet_vec));

    sort(bet_vec.begin(), bet_vec.end(), [](const bet_pair& l, const bet_pair& r) {
        if (l.second != r.second) {
            return l.second < r.second;
        }
        return l.first < r.first;
    });

    reverse(bet_vec.begin(), bet_vec.end());

    for (auto it = bet_vec.begin(); it != bet_vec.begin() + 5; ++it) {
        cout << (*it).first << ": " << (*it).second << endl;
    }

    cout << endl;
}

tuple<vector<string>, map<string, vector<string>>, map<string, float>>
single_source_shortest_path_basic(vector<string> nodes, const string &s, map<string, vector<string>> connections) {
    vector<string> S;
    map<string, vector<string>> P;
    map<string, float> D;
    map<string, float> sigma;
    vector<string> Q;

    // create a path and distance for each node
    for (const auto &v : nodes) {
        P[v] = vector<string>();
        sigma[v] = 0.0;
    }

    // The distance of the source is 0
    D[s] = 0;
    sigma[s] = 1.0;

    Q.push_back(s);

    string v;
    float Dv;
    float sigmav;

    while (!Q.empty()) {
        v = Q[0];
        Q.erase(Q.begin());

        S.push_back(v);

        Dv = D[v];
        sigmav = sigma[v];

        for (const auto &w : connections[v]) {
            if (D.find(w) == D.end()) {
                Q.push_back(w);
                D[w] = Dv + 1;
            }
            if (D[w] == Dv + 1) {
                sigma[w] += sigmav;
                P[w].push_back(v);
            }
        }
    }

    tuple<vector<string>, map<string, vector<string>>, map<string, float>> return_tuple;
    return_tuple = make_tuple(S, P, sigma);
    return return_tuple;

};

map<string, float> accumulate_basic(map<string, float> betweenness, vector<string> S, map<string, vector<string>> P, map<string, float> sigma, string s) {

    map<string, float> delta;
    for (auto const &node : S) {
        delta[node] = 0.0;
    }

    string w;
    float coeff;

    while (!S.empty()) {
        w = S.back();
        S.pop_back();
        coeff = (1 + delta[w]) / sigma[w];
        for (auto const &v : P[w]) {
            delta[v] += sigma[v] * coeff;
        }
        if (w != s) {
            betweenness[w] += delta[w];
        }
    }



    return betweenness;

}

map<string, float> rescale(map<string, float> betweenness, int n) {
    double scale;
    if (n <= 2) {
        scale = 1.0;
    } else {
        cout << "scale = 1.0 / ((" << n << " - 1) * (" << n << " - 2))" << endl;
        scale = 1.0 / ((n - 1) * (n - 2));
    }
    cout << "scale = " << scale << endl;
    for (auto& v : betweenness) {
        v.second *= scale;
    }
    return betweenness;
}

int main() {
    map<string, float> betweenness;
    map<string, vector<string>> connections;
    param_tuple return_tuple;
    vector<string> nodes;

    string file_name = R"(/mnt/c/Users/Evan/PycharmProjects/fractal_proteins/connections/BIOGRID-ORGANISM-Homo_sapiens-3.5.165.tab2.connections)";

    // open file of connections
    ifstream inf(file_name);

    // check file could be opened
    if (!inf) {
        cout << "couldn't open file" << endl;
        return(1);
    }

    // read until end of file
    while (inf) {
        string strInput;
        string token;
        size_t pos = 0;
        string delimiter = "\t";
        vector<string> line;
        string first;

        getline(inf, strInput);

        // split line by tab characters
        while ((pos = strInput.find(delimiter)) != string::npos) {
            // get the token before the first tab
            token = strInput.substr(0, pos);
            // delete the token and the tab
            strInput.erase(0, pos + delimiter.length());
            // skip empty lines
            if (token.empty()) {
                continue;
            }
            // store the name of the parent protein
            if (first.empty()) {
                first = token;
            } else {
                // if not the parent, add to list of children
                line.push_back(token);
            }
        }
        // map children onto parents
        if (!first.empty()) {
            connections[first] = line;
        }
    }

    // get a list of the nodes present in the graph
    nodes = extract_keys(connections);
    int num_nodes = nodes.size();

    // set the betweenness of all nodes equal to 0.0
    for (const auto &node : nodes) {
        betweenness[node] = 0.0;
    }

    print_top_five(betweenness);

    int counter = 0;
    int report_iterations = 1;

    // record how long it takes to perform the process on a single node
    chrono::milliseconds elapsed{};
    auto start = chrono::high_resolution_clock::now();
    double remaining;


    for (const auto &s : nodes) {
        return_tuple = single_source_shortest_path_basic(nodes, s, connections);
        vector<string> S = get<0>(return_tuple);
        map<string, vector<string>> P = get<1>(return_tuple);
        map<string, float> sigma = get<2>(return_tuple);

        if (counter % report_iterations == 0) {
            cout << "node " << counter << " of " << num_nodes << endl;
        }

        ++counter;

        if (counter % report_iterations == 0) {
            elapsed = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
            cout << "elapsed = " << elapsed.count() << endl;
            remaining = (elapsed.count() / (1000.0 * counter)) * (num_nodes - counter);
            cout << "estimated remaining time: " << (int) remaining << endl;
        }

        betweenness = accumulate_basic(betweenness, S, P, sigma, s);
    }

    cout << "num_nodes = " << num_nodes << endl;
    betweenness = rescale(betweenness, num_nodes);

    cout << "\n\n";
    print_top_five(betweenness);

}