#include "geometry.h"

namespace geometry{



Vector operator /(Vector& vector, double value) { 
    Vector result;
    for(auto iten(vector.begin()); iten < vector.end(); ++iten)
           result.push_back(*iten / value); 
    return result;
}

Vector operator *(double value, Vector& vector) { 
    Vector result;
    for(auto iten(vector.begin()); iten < vector.end(); ++iten)
           result.push_back(value * (*iten)); 
    return result;
}

Vector operator +(Vector& vector1, Vector& vector2) { 
    Vector result;
    for (size_t i = 0; i < vector1.size(); ++i) {
        result.push_back(vector1[i] + vector2[i]);
    }
    return result;
}

Vector operator -(Vector& vector1, Vector& vector2) { 
    Vector result;
    for (size_t i = 0; i < vector1.size(); ++i) {
        result.push_back(vector1[i] - vector2[i]);
    }
    return result;
}


bool operator ==(Vector& vector1, Vector& vector2) { 
    for (size_t i = 0; i < vector1.size(); ++i) {
        if(fabs(vector1[i] - vector2[i]) > 0.00001);
            return false;
    }
    return true;
}

double operator*(Point lhs, Point rhs) {
   return DotProduct(lhs, rhs); 
}

//for 2D only, counterclock
Point NormalVector(Point point) {
    if(point.size() == 2) {
        double y = point[0];
        double x = -point[1];
        return {x, y};
    }
    return Point();
}

//for 2D only, counterclock
Point NormalVectorInverted(Point point) {
    if(point.size() == 2) {
        double y = -point[0];
        double x = point[1];
        return {x, y};
    }
    return Point();
}



Straight::Straight(Vector normal_vector, Point point):normal_vector(normal_vector), point(point) {};
const Point& Straight::get_point() const {return point;};
const Vector& Straight::get_normal() const {return normal_vector;};
//returns the meeting point of the lines
Point Straight::operator==(const Straight & rhs) const {
    Point director1 =  NormalVector(rhs.normal_vector);
    Point director2 =  NormalVector(this->normal_vector);

//    std::cout << "ddot "  << director1 * director2 << std::endl;
    //std::cerr << "dot line meet "<<  director2 * rhs.get_normal() << std::endl;
    //Case the lines are parallel
    if(fabs(director2 * rhs.get_normal()) < epsilon){
//        std::cerr << "linhas paralelas" << std::endl;
        return Point();
    }

    double d0 = director1[0]; double d1 = director1[1];
    double p0 = rhs.point[0]; double p1 = rhs.point[1];
    double _d0 = director2[0]; double _d1 = director2[1];
    double _p0 = this->point[0]; double _p1 = this->point[1];
    
//*
    if(fabs(d0) < epsilon and fabs(_d0) > epsilon) {
        double _t = (p0 - _p0) / _d0;
        double x = _d0 * _t + _p0;
        double y = _d1 * _t + _p1;
        return {x, y};
 
    }
    if(fabs(_d0) < epsilon and fabs(d0) > epsilon) {
        double t = (_p0 - p0) / d0;
        double x = d0 * t + p0;
        double y = d1 * t + p1;
        return {x, y};
    }
    if(fabs(d1) < epsilon and fabs(_d1) > epsilon) {
        double _t = (p1 - _p1) / _d1;
        double x = _d0 * _t + _p0;
        double y = _d1 * _t + _p1;
        return {x, y};
 
    }
    if(fabs(_d1) < epsilon and fabs(d1) > epsilon) {
        double t = (_p1 - p1) / d1;
        double x = d0 * t + p0;
        double y = d1 * t + p1;
        return {x, y};
    }

        double t = (_d0 * (p1 - _p1) - _d1 * (p0 - _p0))  / (_d1 * d0 - d1 * _d0);
        double x = d0 * t + p0;
        double y = d1 * t + p1;
        return {x, y};
//*/
}


int WhichSide(Point dot, Straight line) {
    double product = GenVector(line.get_point(), dot) * line.get_normal(); 
    if(product > 0)
        return 1; 
    else
        return -1;
}

bool NormalSide(const Point& dot, const Straight& line) {
    double product = GenVector(line.get_point(), dot) * line.get_normal(); 
    if(product > 0)
        return true; 
    else
        return false;
}

Vector GenVector(Point point1, Point point2) {
    Point aux;
    for (size_t i = 0; i < point1.size(); ++i)
    {
        aux.push_back(point2[i] - point1[i]);
    }
    return aux;
}

Point& NearestNeighbor(std::vector<Point*> list, const Point& point) {
    double max_dist = DBL_MAX;
    Point * neighbor = NULL;
    for(auto iten(list.begin()); iten < list.end(); iten++) {
        double distance = EuclideanDistance(*(*iten), point);
        if(*iten != &point && distance < max_dist) {
            max_dist = distance;
            neighbor = *iten;
        }
    }
    return(*neighbor);
}

std::vector<Point*> ToPointers(std::vector<Point>& list) {
    std::vector<Point*> pointers;
    for(auto iten(list.begin()); iten < list.end(); iten++){
        pointers.push_back(&(*iten));
    }
   return(pointers); 
}

std::vector<Point*> OnLeftSide(Straight& line, std::vector<Point*>& list) {
    std::vector<Point*> points;
    for(auto iten(list.begin()); iten < list.end(); ++iten)
        if(WhichSide(**iten, line) == -1 )
            points.push_back(*iten);
    return points;
}

std::vector<Point*> OnNormalSide(Straight& line, std::vector<Point*>& list) {
    std::vector<Point*> points;
    for(auto iten(list.begin()); iten < list.end(); ++iten)
        if(NormalSide(*(*iten), line))
            points.push_back(*iten);
    return points;
}

void DivideDots(Straight& line, std::vector<Point*>& list, std::vector<Point*>& result1, std::vector<Point*>& result2) {
    for(auto iten(list.begin()); iten < list.end(); ++iten)
        if(NormalSide(*(*iten), line))
            result1.push_back(*iten);
        else
            result2.push_back(*iten);
}

} //namespace geometry
