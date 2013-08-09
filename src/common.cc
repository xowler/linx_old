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
    deque<vector<float> > d0;
    int dim = -1;

    { // Read File
        string line;
        vector<float> fline;

        ifstream f(fn.c_str()); // ToDo: Move to some file opener
        if(f.fail()) {
            cout << "Error: File could not be open";
            return d1;
        }

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
        f.close();
    }
    
    { // Construct Matrix
        d1.resize(d0.size(), dim); // ToDo: Move to iter2DToEigen
        int rn=0, cn=0;
        while( !d0.empty() ){
            cn=0;
            for(vector<float>::iterator i=(d0.front()).begin(); 
                i!=(d0.front()).end(); ++i, ++cn){
                d1(rn,cn)= (*i);
            }
            d0.pop_front();
            rn++;
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

#define check(a,b,c) { \
        float aa = a;\
        float bb = b;\
        if ( fabs(aa-bb)<0.001 ){\
                cerr << "PASSED: " << c << endl; \
        } else { \
                cerr << "FAILED: " << c << "("<<aa<<" "<<bb<<")"<<endl;\
        }\
}


#ifdef _COMMON_TEST
int main(void){
    Eigen::MatrixXf d = loadtxt("test/test2D.txt");
    
    check(1.697860215697601305e+01, d(0,0), "First Line");
    check(3.019403097946556613e+02, d(9999,1), "Last Line"); 
    
    
    
}
#endif
