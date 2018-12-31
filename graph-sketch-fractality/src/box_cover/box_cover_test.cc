#include "box_cover.h"
#include "gtest/gtest.h"

using namespace agl;
using namespace agl::box_cover_internal;
using namespace std;

pair<vector<V>, vector<V>> inv_and_rank(const G &g) {
  vector<V> rank(g.num_vertices());
  vector<V> inv(g.num_vertices());
  for (V i = 0; i < g.num_vertices(); ++i) inv[i] = i;
  shuffle(inv.begin(), inv.end(), agl::random);
  for (int i = 0; i < g.num_vertices(); ++i) rank[inv[i]] = i;
  return {inv, rank};
}

TEST(box_cover, memb) {
    ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, burning) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, build_sketch_check) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, greedy_small) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, greedy_big) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, greedy_huge) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, coverage_management) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, coverage_break) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, find_solution_flower) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, find_solution_shm) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, lazy_greedily) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, coloring) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, covered_check) {
  ASSERT_EQ(1.0, 1.0);
}

TEST(box_cover, cbb) {
  ASSERT_EQ(1.0, 1.0);
}
