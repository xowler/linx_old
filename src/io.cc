#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <deque>
#include <Eigen/Dense>

using namespace std; 


vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) elems.push_back(item); 
    return elems;
}

vector<string> split(const string &s, char delim=' ') {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

vector<float> toFloat(vector<string> v){
    vector<float> res;
    for(vector<string>::iterator i=v.begin(); i!=v.end(); ++i){
        res.push_back(atof(i->c_str()));
    } 
    return res;
}



Eigen::MatrixXf loadtxt(string fn){
    Eigen::MatrixXf d1;
    ifstream f(fn.c_str());
    if(f.fail()) {
        cout << "Error: File could not be open";
        return d1;
    }

    
    deque<vector<float> > d0;
    string line;
    vector<float> fline;

    int dim = -1;
    while(getline(f,line)){ // ToDo: Move to get deque
        fline = toFloat(split(line));
        if(dim==-1) dim=fline.size();
        else{
            if(fline.size()!=dim){
                cout << "Dimension mismatch " << dim << " != " 
                     << fline.size() << "\n   ->" << line  << "\n";
            }
        }
        d0.push_back(fline); 
    }
    
    d1.resize(d0.size(), dim); // ToDo: Move to iter2DToEigen
    int rn=0, cn=0;
    for(deque<vector<float> >::iterator r=d0.begin();
        r!=d0.end(); ++r, ++rn){
        cn=0;
        for(vector<float>::iterator i=r->begin(); 
            i!=r->end(); ++i, ++cn){
            d1(rn,cn)= (*i);
        }
    }
    return d1;
}


template <typename T, 
    template<typename ELEM, typename ALLOC=std::allocator<ELEM> > 
        class Container >
ostream& operator<< (ostream& o, const Container<T>& container) {
        typename Container<T>::const_iterator beg = container.begin();
        int i = 0;
        o << "[\n" ; 
        while(beg != container.end()) {
            o << " " << *beg++ << "\n"; 
            if(++i>30){
                o << " ... \n";
                break;
            }
        }
        o << "]" << endl; 
        return o;
}


#ifdef _IO_TEST
int main(void){
    loadtxt("test/test2D.txt");
    
}
#endif
