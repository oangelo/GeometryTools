#ifndef geometry_H
#define geometry_H 

#include <math.h>
#include <float.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <ctime>

namespace geometry{

    typedef std::vector<double> Point;
    typedef std::vector<double> Vector;

    template <class Type>
        double EuclideanDistance(const Type & ac1,const Type & ac2);
    template <class Type>
        std::vector<double> Versor(const Type & ac1,const Type & ac2);
    template <class Type>
        bool InSight(const Type & point1, const Type & point2, std::vector<double> &direction, double cos_angle);
    template <class Type>
        const double DotProduct( const Type & vet1, const Type & vet2);
        Vector operator /(Vector& vector, double value);


    double operator*(Point lhs, Point rhs);
    Point NormalVector(Point point); //2D only

    class Straight {
        public:
            Straight(Vector normal_vector, Point point);

            const Point& get_point() const;
            const Vector& get_normal() const;
            //returns the meeting point of the lines
            Point operator==(const Straight & rhs) const;

        private:
                Vector normal_vector;
                Point point;
    };

    //returns the position where the dot is, 
    //on the left or right side of the line
    //right side is the direction of the normal 
    //vector. 1 = right and -1 = left
    int WhichSide(Point dot, Straight line);
    //generate a vector form two points
    Vector GenVector(Point point1, Point point2);
    //find the neares neighbor from a list of points
    Point& NearestNeighbor(std::vector<Point*> list, const Point& point);
    //Return the points that are not on the same side of the normal vector
    std::vector<Point*> OnLeftSide(Straight& line, std::vector<Point*> list);
     //return a vector of pointer
    std::vector<Point*> ToPointers(std::vector<Point>& list);
}

//#include "geometry.cpp"

#endif /* geometry_H */
