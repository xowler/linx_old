#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <deque>
#include <map>
#include <Eigen/Dense>

using namespace std; 

// iters and string stuff
string trim(string s){
    int n=s.size();
    int i=0, j=n-1; 
    for(;s[i]==' ' && i<n; ++i);
    for(;s[j]==' ' && j>0; --j);
    return (i<=j)? s.substr(i,j-i+1) : string("");
}

vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        if(delim!=' ' || item.size()>0)
            elems.push_back(item); 
    }
    return elems;
}
vector<string> split(const string &s, char delim=' ') {
    vector<string> elems;
    split(s, delim, elems);
    return elems;
}

vector<float> &to_float(vector<string> v,vector<float> &res){
    res.clear();
    for(vector<string>::iterator i=v.begin(); i!=v.end(); ++i) 
        res.push_back(atof(i->c_str()));
    return res;
    
}
vector<float> to_float(vector<string> v){
    vector<float> res;
    to_float(v,res);
    return res;
}

// IO
void fopen(string fn, ifstream &f){
    f.open(fn.c_str());
    if(f.fail()){
        cout << "ERROR: File could not be open";
    } 
}
void read_lines(string fn, deque<string> &lines){
    deque<string> res;
    string line;
    ifstream f;
    fopen(fn,f);
    while(getline(f,line)) res.push_back(line);
    res.swap(lines);
}
deque<string> read_lines(string fn){
    deque<string> d;
    read_lines(fn, d);
    return d;
}
Eigen::MatrixXf load_txt(string fn){
    Eigen::MatrixXf d1;
    deque<vector<float> > d0;
    int dim = -1;

    { // Read File
        string line;
        vector<float> fline;

        ifstream f;
        fopen(fn,f);

        while(getline(f,line)){ // ToDo: Move to get deque
            fline = to_float(split(line));
            if(dim==-1) dim=fline.size();
            else{
                if(fline.size()!=dim){ // ToDo: throw error
                    cout << "ERROR: Dimension mismatch " << dim << " != " 
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

// Models
struct atom {
    int id;
    int rid;
    string name;
    string rname;
    string chain;
    float mass;
    string b;
};
struct model {
    Eigen::Matrix<float,Eigen::Dynamic,3> x;
    vector<atom> p;
    int n;
    string extra;
};
vector<float> inter_residue_distance(){}

model read_gro(string fn){ // Not biggie copying it, rather small
    model res;
    deque<string> f = read_lines(fn);

    res.extra = f[0] + "\n" + f[f.size()-1];
    f.pop_front();
    f.pop_front();
    f.pop_back();
    res.x.resize(f.size(),3);
    res.p.reserve(f.size());
    res.n = f.size();

    string *l;
    atom a;
    int r=0;
    vector<float> x;
    while(!f.empty()){
        l = &(f.front());

        a.rid = atoi((l->substr(0,5)).c_str());
        a.rname = trim(l->substr(5,5));
        a.name = trim(l->substr(10,5));
        a.id = atoi((l->substr(15,6)).c_str());
        res.p.push_back(a);
        //cout << a.rid << "/" << a.rname << "/" << a.id << "/" << a.name << "/" << endl;

        to_float(split(l->substr(21,24)),x);
        res.x(r,0) = x[0];
        res.x(r,1) = x[1];
        res.x(r,2) = x[2];
        ++r;
        //cout << x[0] << '/' << x[1] << '/' <<x[2] << '/' << endl;
        f.pop_front();
    }
    
    return res;
}

vector<vector<int> > residuize(model &m){
    vector<vector<int> > res;
    if(m.n==0) cout << "ERROR: Something is wrong with this model. Zero lenght" << endl;
    string lr = "fake";
    for(int i=0; i<m.n; i++){
        if(lr!=m.p[i].rname) res.push_back(vector<int>());
        res[res.size()-1].push_back(i);
        lr = m.p[i].rname;
    }
    return res;
}

// Testing
#define ftest(a,b,c) { \
        float aa = a;\
        float bb = b;\
        if ( fabs(aa-bb)<0.001 ){\
                cerr << "PASSED: " << c << endl; \
        } else { \
                cerr << "!!!!!!!!!!FAILED: " << c << "("<<aa<<"|"<<bb<<")"<<endl;\
        }\
}
#define stest(a,b,c) { \
        string aa = a;\
        string bb = b;\
        if ( aa==bb ){\
                cerr << "PASSED: " << c << endl; \
        } else { \
                cerr << "!!!!!!!!!!FAILED: " << c << "("<<aa<<"|"<<bb<<")"<<endl;\
        }\
}

#ifdef _COMMON_TEST
int main(void){
    {
        cout << endl << "-- GRO " << endl;
        // readGRO
        model gro = read_gro("test/bhp.gro");
        ftest(gro.x(207,1),0.375,"Coordinates from file");
        stest(gro.p[184].name,"N","Attributes from file");

        // Model
        ftest(residuize(gro)[8][5],146,"Identify residues");
    }
    
    // Iter and strings
    {
        cout << endl << "-- split " << endl;
        ftest(7,to_float(split("0 1 2 3 4 5 6  7"))[7], "Split and toFloat ' '");
        ftest(7,to_float(split("0/1/2/3/4/5/6//7",'/'))[8], "Split and toFloat '/'");
    }

    {
        cout << endl << "-- trim " << endl;
        stest("",trim("       "),"Trim an empty string");
        stest("1 2",trim("      1 2  "),"Trim a string");
        stest("1 2",trim("      1 2"),"Trim a string left");
        stest("1 2",trim("1 2  "),"Trim a string right");
        stest("1",trim("1   "),"Trim a string right");
        stest("2",trim("    2"),"Trim a string left");
    }

    // Loadtxt
    {
        cout << endl << "-- loadtxt " << endl;
        Eigen::MatrixXf d = load_txt("test/test2D.txt");
        ftest(1.697860215697601305e+01, d(0,0), "First Line");
        ftest(3.019403097946556613e+02, d(9999,1), "Last Line"); 
    }
}
#endif
