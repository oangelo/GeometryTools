#include <set> 
#include <vector> 

#include "gtest/gtest.h"

#include "../geometry.h"


using namespace geometry;

TEST(geometry, straigth) {
    Straight line({3,4}, {1,2});
    EXPECT_EQ(line.get_point()[0], 1);
    EXPECT_EQ(line.get_point()[1], 2);
    EXPECT_EQ(line.get_normal()[0], 3);
    EXPECT_EQ(line.get_normal()[1], 4);
}

TEST(geometry, operator_multiply) {
    std::vector<double> vec1({1,0});
    std::vector<double> vec2({0,1});
    double dot = vec1 * vec2;
    EXPECT_EQ(dot, 0);
}

TEST(geometry, Left_right_of_the_line) {
    Straight line({0,1}, {0,0});
    EXPECT_EQ(WhichSide({1,0}, line), -1);
    EXPECT_EQ(WhichSide( {0,1}, line), 1);
}

TEST(geometry, Normal_Vector) {
    Point point({2, 1});
    Point normal = NormalVector({2, 1});
    EXPECT_EQ(point * normal, 0.0);
    
}

TEST(geometry, Line_meeting_point) {
    Straight line1({0,1}, {0,0});
    Straight line2({1,0}, {0,0});
    Straight line3({1,1}, {1,1});
    EXPECT_EQ((line1 == line2)[0], 0);
    EXPECT_EQ((line1 == line2)[1], 0);
    EXPECT_EQ((line1 == line3)[0], 2);
    EXPECT_EQ((line1 == line3)[1], 0);
    EXPECT_EQ((line2 == line3)[0], 0);
    EXPECT_EQ((line2 == line3)[1], 2);
    //Parallel lines should not have a meeting point
    EXPECT_EQ((line2 == line2).size(), 0);
}

TEST(geometry, vector_from_points) {
    Vector vec = GenVector({2, 2}, {-1, -1});
    EXPECT_EQ(vec[0], -3);
    EXPECT_EQ(vec[0], -3);
}


TEST(geometry, element_to_pointer) {
    std::vector<Point> points;
    points.push_back({0, 0});
    points.push_back({1, 1});
    points.push_back({2, -2});
    points.push_back({3, 1});

    std::vector<Point*> points_pointers(ToPointers(points));

    for (size_t i = 0; i < 4; ++i)
        EXPECT_TRUE(&(points[i]) ==  points_pointers[i]);
 
}

TEST(geometry, NearNeighbor) {
    std::vector<Point> points;
    points.push_back({0, 0});
    points.push_back({1, 1});
    points.push_back({2, -2});
    points.push_back({3, 1});
    std::vector<Point*> pointers(ToPointers(points));
    Point point(NearestNeighbor(pointers, points[0]));
    EXPECT_EQ(point[0], 1);
    EXPECT_EQ(point[1], 1);
    EXPECT_TRUE(pointers[0] == &(points[0]));
}

TEST(geometry, Line_side) {
    Straight line({1,1}, {0,0});
    std::vector<Point> points;
    points.push_back({1, 1});
    points.push_back({1, 2});
    points.push_back({-2, -2});
    points.push_back({-3, 1});
    std::vector<Point*> pointers(ToPointers(points));
    auto result(OnLeftSide(line, pointers));
    EXPECT_EQ(result.size(), 2);
    EXPECT_TRUE(result[0] == &points[2]);
    EXPECT_TRUE(result[1] == &points[3]);
}

TEST(geometry, Vector_dividion) {
    Vector vector({2, 2});
    Vector div(vector / 2.0);
    EXPECT_EQ(div[0], 1);
    EXPECT_EQ(div[1], 1);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}