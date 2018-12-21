#include "easy_cui.h"
#include "box_cover.h"
using namespace std;
using namespace agl::box_cover_internal;

DEFINE_int32(rad_min, 0, "minimum radius");
DEFINE_int32(rad_max, 100000000, "maximum radius");
DEFINE_string(method, "sketch", "using method");
DEFINE_double(least_coverage, 1.0, "coverage");
DEFINE_int32(multipass, 1000000000, "Number of multi-pass");
DEFINE_int32(sketch_k, 128, "sketch k");
DEFINE_bool(rad_analytical, false, "Using analytical diameters for rads");
DEFINE_double(alpha, 1.0, "index size limit to use MEMB (alpha*n*k)");

vector<unweighted_edge_list>  extract_maximal_connected(const G& g_pre) {
  unweighted_edge_list es;
  V num_v = g_pre.num_vertices();
  size_t max_num_v = 0;
  int max_g = 0;
  vector<bool> vis(num_v, false);
  vector<vector<V>> connected_sets;
  for (V v = 0; v < num_v; ++v) {
    if (vis[v]) continue;
    vector<V> connected_v;
    queue<V> que;
    que.push(v);
    vis[v] = true;
    connected_v.push_back(v);
    while (!que.empty()) {
      V u = que.front();
      que.pop();
      for (V x : g_pre.neighbors(u)) {
        if (!vis[x]) {
          vis[x] = true;
          que.push(x);
          connected_v.push_back(x);
        }
      }
    }
    if (max_num_v < connected_v.size()) {
      max_num_v = connected_v.size();
      max_g = connected_sets.size();
      // sort(connected_v.begin(), connected_v.end());
      // connected_sets.push_back(connected_v);
    }
    sort(connected_v.begin(), connected_v.end());
    connected_sets.push_back(connected_v);
  }

  vector<V>& max_c = connected_sets[max_g];
  vector<V> inv(num_v, -1);
  for (V v : max_c) {
    inv[v] = (lower_bound(max_c.begin(), max_c.end(), v) - max_c.begin());
  }

  for (pair<V, V>& e : g_pre.edge_list()) {
    V from = e.first;
    V to = e.second;
    if (inv[from] == -1) continue;
    es.emplace_back(inv[from], inv[to]);
  }

  vector<unweighted_edge_list> connected_edges_list;
  connected_edges_list.push_back(es);
  for(unsigned int i = 0; i < connected_sets.size(); i++)
  {
      if((int)i != max_g)
      {
          es.clear();
          vector<V>& max_c = connected_sets[i];
          vector<V> inv(num_v, -1);

          for (V v : max_c) {
              inv[v] = (lower_bound(max_c.begin(), max_c.end(), v) - max_c.begin());
          }

          for (pair<V, V>& e : g_pre.edge_list()) {
              V from = e.first;
              V to = e.second;
              if (inv[from] == -1) continue;
              es.emplace_back(inv[from], inv[to]);
          }
          if(es.size() > 0)
          {
            connected_edges_list.push_back(es);
          }
      }
  }

  return connected_edges_list;
}


int main(int argc, char** argv) {
  auto str_replace = [](string& str, const string& from, const string& to) {
    string::size_type pos = 0;
    while (pos = str.find(from, pos), pos != string::npos) {
      str.replace(pos, from.length(), to);
      pos += to.length();
    }
  };

  // Extract maximal connected subgraph & Compress coordinates
  /*unweighted_edge_list es;
  {
    G g_pre = easy_cui_init(argc, argv);
    CHECK_MSG(FLAGS_force_undirected, "undirected only!!!");
    es = extract_maximal_connected(g_pre);
  }*/
  G g_pre = easy_cui_init(argc, argv);
  vector<unweighted_edge_list> es = extract_maximal_connected(g_pre);
  G g(es[0]);
  const int number_of_nodes = g_pre.num_vertices();

  // Output information of graph

  int total_nodes_in_connected_components = 0;
  for(unweighted_edge_list edge_list: es)
  {
    G g(edge_list);
    total_nodes_in_connected_components = total_nodes_in_connected_components + g.num_vertices();
  }
  const int total_number_of_disconnected_components = number_of_nodes - total_nodes_in_connected_components;
  const int total_number_of_components = total_number_of_disconnected_components + es.size();

  pretty_print(g);
  JLOG_ADD_OPEN("graph_info") {
    JLOG_PUT("vertices", number_of_nodes);
    JLOG_PUT("edges", g_pre.num_edges());
    JLOG_PUT("connected_components", es.size());
    JLOG_PUT("disconnected_components", total_number_of_disconnected_components);
    string gstr = FLAGS_graph;
    str_replace(gstr, " ", "-");
    JLOG_PUT("graph", gstr);
  }

  // Extract Graph Name
  const char* slash = strrchr(FLAGS_graph.c_str(), '/');
  string graph_name = slash ? slash + 1 : FLAGS_graph;
  str_replace(graph_name, " ", "-");

  // Prepare radius to calculate
  vector<W> rads;
  if (FLAGS_rad_analytical) {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    V required, u, v;
    if (!(iss >> required)) required = 44;
    if (!(iss >> u)) u = 2;
    if (!(iss >> v)) v = 2;
    auto as = find_analytical_solution(family, u, v, g);
    for (auto p : as) {
      W diameter = p.first;
      W rad = (diameter + 1) / 2;
      if (rad >= FLAGS_rad_max) break;
      rads.push_back(rad);
    }
  } else {
    for (W rad = FLAGS_rad_min; rad <= FLAGS_rad_max; ++rad)
      rads.push_back(rad);
  }

  string experiment_name = FLAGS_method;
  if (FLAGS_method == "sketch") {
    if (FLAGS_multipass == 1) FLAGS_least_coverage = 1.0;
    experiment_name.append("-k." + to_string(FLAGS_sketch_k));

    JLOG_PUT("name", experiment_name);
    string sk = to_string(FLAGS_sketch_k);
    while (sk.size() < 4) sk = "0" + sk;
    JLOG_PUT("k", to_string(FLAGS_sketch_k));
    JLOG_PUT("pass", to_string(FLAGS_multipass));
    JLOG_PUT("alpha", to_string(FLAGS_alpha));

    for (W rad : rads) {
      vector<V> res;
      coverage_manager cm(g, rad, FLAGS_least_coverage);
      JLOG_ADD_BENCHMARK("time")
      res = box_cover_sketch(g, rad, FLAGS_sketch_k, FLAGS_multipass, cm,
                             FLAGS_alpha);
      JLOG_ADD("size", res.size());
      JLOG_ADD("radius", rad);
      JLOG_ADD("coverage", cm.get_current_coverage());
      JLOG_ADD_OPEN("centers") {
        for (const auto& b : res) JLOG_ADD(to_string(rad).data(), b);
      }
      if (res.size() == 1) break;
    }
  } else if (FLAGS_method == "memb") {
    JLOG_PUT("name", "MEMB");
    for (W rad : rads) {
      vector<V> res;
      int total_boxes = 0;
      if(rad == 0)
      {
        total_boxes = number_of_nodes;
        JLOG_ADD_BENCHMARK("time");
        JLOG_ADD("size", total_boxes);
        JLOG_ADD("radius", rad);
        coverage_manager cm(g, rad, FLAGS_least_coverage);
        JLOG_ADD("coverage", 1);
      }
      else
      {
        JLOG_ADD_BENCHMARK("time") for(unweighted_edge_list edge_list: es)
        {
          G g(edge_list);
          res = box_cover_memb(g, rad);
          total_boxes += res.size();
        }
        total_boxes += total_number_of_disconnected_components;
        JLOG_ADD("size", total_boxes);
        JLOG_ADD("radius", rad);
        coverage_manager cm(g, rad, FLAGS_least_coverage);
        for (auto& v : res) cm.add(g, v);
        JLOG_ADD("coverage", cm.get_current_coverage());
      }
      if ((unsigned int)total_boxes == total_number_of_components) break;
    }
  } else if (FLAGS_method == "cbb") {
    JLOG_PUT("name", "CBB");
    for (W rad : rads) {
      vector<V> res;
      JLOG_ADD_BENCHMARK("time") res = box_cover_cbb(g, rad * 2);
      JLOG_ADD("size", res.size());
      JLOG_ADD("radius", rad);
      coverage_manager cm(g, rad, FLAGS_least_coverage);
      for (auto& v : res) cm.add(g, v);
      JLOG_ADD("coverage", cm.get_current_coverage());
      if (res.size() == 1) break;
    }
  } else if (FLAGS_method == "coloring") {
    JLOG_PUT("name", "Coloring");
    W rad_max = rads[rads.size() - 1];
    vector<pair<W, size_t>> res;

    JLOG_ADD_BENCHMARK("time") res = box_cover_coloring(g, rad_max * 2);
    if (FLAGS_rad_analytical) {
      for (auto p : res)
        if (p.first % 2 == 0 &&
            find(rads.begin(), rads.end(), p.first / 2) != rads.end()) {
          JLOG_ADD("size", p.second);
          JLOG_ADD("diameter", p.first);
        }
    } else
      for (auto p : res) {
        JLOG_ADD("size", p.second);
        JLOG_ADD("diameter", p.first);
      }
  } else if (FLAGS_method == "burning") {
    JLOG_PUT("name", "Burning");
    for (W rad : rads) {
      vector<V> res;
      JLOG_ADD_BENCHMARK("time") res = box_cover_burning(g, rad);
      JLOG_ADD("size", res.size());
      JLOG_ADD("radius", rad);
      coverage_manager cm(g, rad, FLAGS_least_coverage);
      for (auto& v : res) cm.add(g, v);
      JLOG_ADD("coverage", cm.get_current_coverage());
      if (res.size() == 1) break;
    }
  } else if (FLAGS_method == "analytical") {
    istringstream iss(FLAGS_graph);
    string family;
    iss >> family;
    JLOG_PUT("name", "analytical_solution_" + family);
    V required, u, v;
    if (!(iss >> required)) required = 44;
    if (!(iss >> u)) u = 2;
    if (!(iss >> v)) v = 2;
    auto as = find_analytical_solution(family, u, v, g);

    for (auto p : as) {
      JLOG_ADD("size", p.second);
      JLOG_ADD("diameter", p.first);
      JLOG_ADD("radius", p.first / 2 + p.first % 2);
    }
  } else {
    cerr << "Unknown method: " << FLAGS_method << endl;
    return 0;
  }
}
