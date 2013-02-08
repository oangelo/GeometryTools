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


