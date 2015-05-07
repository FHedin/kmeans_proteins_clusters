#ifndef CLUSTER_HPP
#define CLUSTER_HPP

/*
 * This class represents an atomic cluster
 */
class cluster{
    
private:
    double x;
    double y;
    double z;
    //unsigned int nStates;
    
public:
    cluster(double _x, double _y, double _z){ x=_x ; y=_y ; z=_z; /*nStates=1;*/}
    
    //void incrementnStates(){nStates++;}
    
    double getX(){return x;}
    double getY(){return y;}
    double getZ(){return z;}
    
    //unsigned int getNStates(){return nStates;}
    
    void setX(double _x){x=_x;}
    void setY(double _y){y=_y;}
    void setZ(double _z){z=_z;}
    
    void mergeWith(cluster& toMerge)
    {
        x += toMerge.getX();
        y += toMerge.getY();
        z += toMerge.getZ();
        x /= 2.0;
        y /= 2.0;
        z /= 2.0;
    }
    
};

#endif  /* CLUSTER_HPP */
