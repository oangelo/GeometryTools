template <class Type>
const double DotProduct( const Type & vet1, const Type & vet2){
    double vet_aux=0;
    for(unsigned i=0; i < vet1.size(); i++){
        vet_aux += vet1[i] * vet2[i];
    }
    return(vet_aux);
}

template <class Type>
bool InSight(const Type & point1, const Type & point2, std::vector<double> & direction, double cos_angle){
//This function should receive cos(angle_of_vision / 2.0)
    std::vector<double> vet_ag(Versor(point1,point2));
    double dot = DotProduct(vet_ag, direction);
    if(dot > (cos_angle)){
        return(1);
    }else{
        return(0);
    }
}

template <class Type>
double EuclideanDistance(const Type & v1,const Type & v2){
    std::vector<double> distance(v1.size());   
    for (size_t i = 0; i < v1.size(); ++i)
    {
        distance[i] = v1[i] - v2[i];
    }
    double dot = (std::inner_product(distance.begin(),distance.end(),distance.begin(),0.0)); 
    return sqrt(dot); 
}

template <class Type>
std::vector<double> Versor(const Type & ac1,const Type & ac2){
    double e_distance = EuclideanDistance(ac1,ac2);
    std::vector<double> v1(ac1.get_position()), v2(ac2.get_position());
    std::vector<double> distance(v1.size());   
    std::transform(v2.begin(),v2.end(),v1.begin(),distance.begin(),std::minus<double>());
    std::for_each(distance.begin(),distance.end(), [e_distance] (double &element){
                element = (element / e_distance);
            });
    return distance;
}


