#ifndef geometry_H
#define geometry_H 

#include <math.h>
#include <float.h>

#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <ctime>
#include <iterator>

namespace geometry{

    const double epsilon = 0.00000001;

    typedef std::vector<double> Point;
    typedef std::vector<double> Vector;

    class Straight {
        public:
            Straight(Vector normal_vector, Point point);

            const Point& get_point() const;
            const Vector& get_normal() const;
            //returns the meeting point of the lines
            Point operator==(const Straight & rhs) const;
            //double operator==(const Point& rhs) const;

        private:
            Vector normal_vector;
            Point point;
    };

    template <class T>
        std::ostream& operator<<(std::ostream& os, const std::vector<T>& v); 
    template <class Type>
        double EuclideanDistance(const Type & ac1,const Type & ac2);
    template <class Type>
        Vector Rotate(const Type & ac1,double angle);
    template <class Type>
        Vector Versor(const Type & ac1,const Type & ac2);
    template <class Type>
        Vector Versor(const Type & ac1,const Type & ac2, double euclidean_distance);

    template <class Type>
        bool InSight(const Type & point1, const Type & point2, const Type& direction, double cos_angle);
    template <class Type>
        const double DotProduct( const Type & vet1, const Type & vet2);
    Vector operator /(Vector& vector, double value);
    Vector operator *(double value, Vector& vector);
    Vector operator +(Vector& vector1, Vector& vector2);
    Vector operator -(Vector& vector1, Vector& vector2);
    bool operator == (Vector& vector1, Vector& vector2);
    template <class T>
        std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);


    double operator*(Point lhs, Point rhs);
    Point NormalVector(Point point); //2D only
    //Return a nomal vector, but in the inverted direction
    Point NormalVectorInverted(Point point); //2D only


    //returns the position where the dot is, 
    //on the left or right side of the line
    //right side is the direction of the normal 
    //vector. 1 = right and -1 = left
    int WhichSide(Point dot, Straight line);
    bool NormalSide(const Point& dot,const Straight& line);
    //generate a vector form two points
    Vector GenVector(Point point1, Point point2);
    //find the neares neighbor from a list of points
    Point& NearestNeighbor(std::vector<Point*> list, const Point& point);
    //Return the points that are not on the same side of the normal vector
    std::vector<Point*> OnLeftSide(Straight& line, std::vector<Point*>& list);
    //Return the points that are on the same side of the normal vector
    std::vector<Point*> OnNormalSide(Straight& line, std::vector<Point*>& list);
    //divide a set of points in to sets by a line
    void DivideDots(Straight& line, std::vector<Point*>& list, std::vector<Point*>& result1, std::vector<Point*>& result2);
    //return a vector of pointer
    std::vector<Point*> ToPointers(std::vector<Point>& list);


    template <class Type>
        const double DotProduct(const Type& vet1, const Type& vet2){
            double vet_aux=0;
            for(unsigned i=0; i < vet1.size(); i++){
                vet_aux += vet1[i] * vet2[i];
            }
            return(vet_aux);
        }

    template <class Type>
        bool InSight(const Type & point1, const Type & point2, const Type& direction, double cos_angle){
            double dot = DotProduct(Versor(point1, point2), direction);
            if(dot > (cos_angle)){
                return(true);
            }else{
                return(false);
            }
        }

    template <class Type>
        double EuclideanDistance(const Type & v1,const Type & v2){
            std::vector<double> distance(v1.size());   
            for (size_t i = 0; i < v1.size(); ++i) {
                distance[i] = v1[i] - v2[i];
            }
            return sqrt((std::inner_product(distance.begin(),distance.end(),distance.begin(),0.0))); 
        }

    template <class Type>
        Vector Versor(const Type & v1,const Type & v2){
            double e_distance(EuclideanDistance(v1, v2));
            std::vector<double> distance(v1.size());   
            std::transform(v2.begin(),v2.end(),v1.begin(),distance.begin(),std::minus<double>());
            std::for_each(distance.begin(),distance.end(), [e_distance] (double &element){
                    element = (element / e_distance);
                    });
            return distance;
        }

    template <class Type>
        Vector Versor(const Type & v1, const Type & v2, double euclidean_distance){
            std::vector<double> distance(v1.size());   
            std::transform(v2.begin(),v2.end(),v1.begin(),distance.begin(),std::minus<double>());
            std::for_each(distance.begin(),distance.end(), [euclidean_distance] (double &element){
                    element = (element / euclidean_distance);
                    });
            return distance;
        }

    template <class Type>
        Vector Rotate(const Type& ac1, double angle){
            Vector new_vec(2); 
            angle = angle * M_PI / 180;
            new_vec[0] = (ac1[0] * cos(angle) - ac1[1] * sin(angle));
            new_vec[1] = (ac1[0] * sin(angle) + ac1[1] * cos(angle));
            return new_vec;
        }

    template <class T>
        std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) 
        { 
            os << "(";
            std::copy(v.begin(), v.end() - 1, std::ostream_iterator<T>(os, ", "));
            std::copy(v.end() - 1, v.end(), std::ostream_iterator<T>(os));
            return os<<")";
        }

}
#endif /* geometry_H */
