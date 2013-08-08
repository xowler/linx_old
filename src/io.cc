#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>

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

vector<float> tofloat(vector<string> v){
    vector<float> res;
    for(vector<string>::iterator i=v.begin(); i!=v.end(); ++i){
        res.push_back(atof(i->c_str()));
    } 
    return res;
}






vector<vector<float> > loadtxt(string fn){
    vector<vector<float> > data;
    ifstream f(fn.c_str());
    if(f.fail()) {
        cout << "Error: File could not be open";
        return data;
    }

    string line;
    vector<float> fline;

    getline(f,line);
    cout << split(line);
    



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
