#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <deque>
#include <map>
#include <cmath>
#include <cstdlib>


using namespace std; 

// iters and string stuff
string _(char *s){ return string(s); }
string _(int s){ 
    stringstream ss;
    ss << s;
    return ss.str(); 
}
string trim(string s){
    int n=s.size();
    int i=0, j=n-1; 
    for(;s[i]==' ' && i<n; ++i);
    for(;s[j]==' ' && j>0; --j);
    return (i<=j)? s.substr(i,j-i+1) : string("");
}
template <class T>
string join(string d, T v){
    if(v.size()==0) return string();
    stringstream ss;
    ss << v[0];
    for(int i=1;i<v.size();i++) ss << d << v[i];
    return ss.str();
}
template <class T>
void print_vector(T v){ cout << "[" << join(",",v) << "]" << endl; }
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
bool fopen(string fn, ifstream &f){
    f.open(fn.c_str());
    if(f.fail()){
        cout << "ERROR: File could not be open" << endl;
        return 0;
    } 
    return 1;
}
bool fopen(string fn){
    ifstream f;
    return fopen(fn, f);
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
void load_txt(string fn, vector<float> &d, int &s ){
    s = -1;
    string line;
    vector<float> fline;

    ifstream f;
    fopen(fn,f);

    while(getline(f,line)){
        fline = to_float(split(line));
        if(s==-1) s=fline.size();
        else{
            if(fline.size()!=s){ // ToDo: throw error
                cout << "ERROR: Dimension mismatch " << s << " != " 
                     << fline.size() << "\n   ->" << line  << "\n";
            }
        }
        for(int i=0; i<s; ++i) d.push_back(fline[i]);
    }
    f.close();
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
    vector<float> x;
    vector<atom> p;
    int n;
    string extra;
};
model read_gro(string fn){ // Not biggie copying it, rather small
    model res;
    deque<string> f = read_lines(fn);

    res.extra = f[0] + "\n" + f[f.size()-1];
    f.pop_front();
    f.pop_front();
    f.pop_back();
    res.x.reserve(f.size()*3);
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
        res.x.push_back(x[0]);
        res.x.push_back(x[1]);
        res.x.push_back(x[2]);
        ++r;
        //cout << x[0] << '/' << x[1] << '/' <<x[2] << '/' << endl;
        f.pop_front();
    }
    
    return res;
}
vector<vector<int> > residuize(model &m, vector<string> &rv, vector<bool> mask){
    vector<vector<int> > res;
    if(m.n==0) cout << "ERROR: Something is wrong with this model. Zero lenght" << endl;
    int lr=-1;
    for(int i=0; i<m.n; i++){
        if(i<mask.size()){
            if(!mask[i]) continue;
        }
        if(lr!=m.p[i].rid) {
            res.push_back(vector<int>());
            rv.push_back(_(m.p[i].rid) + m.p[i].rname);
        }
        res[res.size()-1].push_back(i);
        lr = m.p[i].rid;

        if(m.p[i].rname=="SOL"){
            cout << "ERROR: Solvent present, results not guaranteed" << endl;
            break;    
        }
    }
    return res;
}
vector<vector<int> > residuize(model &m, vector<string> &rv){
    return residuize(m,rv,vector<bool>());
}
vector<bool> heavy(model &m){
    vector<bool> res;
    for(int i=0; i<m.n; i++){
        if(m.p[i].name[0]=='H') res.push_back(0);
        else res.push_back(1);
    }
    return res;
}

// LinAlg
float norm2(float *x){
    return x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
}
float norm(float *x){
    return sqrt(norm2(x));
}
float dist2(float *x, float *y){
    float z[3];
    z[0] = x[0]-y[0];
    z[1] = x[1]-y[1];
    z[2] = x[2]-y[2];
    return norm2(z);
}
float dist(float *x, float *y){
    return sqrt(dist2(x,y));
}
float min(float x, float y){
    return x>y?y:x;
}
float min_dist(float *m, vector<int> g1, vector<int> g2){
    float res = 999999;
    int i,j,I=g1.size(),J=g2.size();
    for(i=0;i<I;i++){
        for(j=0;j<J;j++){
            res = min(res,dist(m+3*g1[i],m+3*g2[j]));
        }
    }
    return res;
}
vector<float> inter_group_distances(vector<vector<int> > r, float *m){
    vector<float> res;

    int N = r.size();
    for(int i=0; i<N; i++){
        for(int j=i+2; j<N;j++){
            res.push_back(min_dist(m,r[i],r[j]));  
        }
    }

    return res;
}


// Testing
#define ftest(a,b,c) { \
        float aa = a;\
        float bb = b;\
        if ( fabs(aa-bb)<0.001 ){\
                cout << "PASSED: " << c << endl; \
        } else { \
                cout << "!!!!!!!!!!FAILED: " << c << "("<<aa<<"|"<<bb<<")"<<endl;\
        }\
}
#define stest(a,b,c) { \
        string aa = a;\
        string bb = b;\
        if ( aa==bb ){\
                cout << "PASSED: " << c << endl; \
        } else { \
                cout << "!!!!!!!!!!FAILED: " << c << "("<<aa<<"|"<<bb<<")"<<endl;\
        }\
}

#ifdef _COMMON_TEST
int main(void){
    {
        cout << endl << "-- Linear Algebra " << endl;
        float x[3] = {1,0,0}, y[3] = {0,0,1};
        ftest(dist(x,y),1.41421,"Metric test");
    }

    { // Model
        cout << endl << "-- GRO " << endl;
        // readGRO
        model gro = read_gro("test/bhp.gro");
        ftest(gro.x[207*3+1],0.375,"Coordinates from file");
        stest(gro.p[184].name,"N","Attributes from file");

        vector<string> res;
        ftest(residuize(gro,res,heavy(gro))[8][5],132,"Identify residues");
        stest(res[12],"13GLY","Residue names");
        ftest(inter_group_distances(residuize(gro,res,heavy(gro)),gro.x.data())[1],
            6.674750931682766897e-01, "Contacts");
    }
    
    // Iter and strings
    {
        cout << endl << "-- split " << endl;
        ftest(4,to_float(split("0 1 2 3  4"))[4], "Split and toFloat ' '");
        ftest(4,to_float(split("0/1/2/3//4/",'/'))[5], "Split and toFloat '/'");
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
        
        int stride;
        vector<float> d;
        load_txt("test/test2D.txt",d,stride);
        ftest(1.697860215697601305e+01, d[0], "First Line");
        ftest(3.019403097946556613e+02, d[9999*stride + 1], "Last Line"); 
    }
}
#else
#include "ccxtc.h"
#include "gzstream.h"
void calc_pca(int argc, char *argv[]){
    bool usage=1;
    if(usage){
        cout << 
            "Usage: \n" <<
            "   " << argv[0] << " calc_pca -f file1.txt [file2.txt ...] -o pca.txt [-s skip] \n" <<
            "\n" <<
            "Description:\n" <<
            "   -f files    - input-    Input files to use as input for the PCA\n" <<
            "   -o pca.txt  -output-    Output file to store the PCA results\n" <<
            "   -s skip     integer     Read every [skip] lines, default 1 \n" <<
            "\n" << endl;
    }
}
void cts(string gro, string xtc, string gz, string txt){
    model m = read_gro(gro);
    
    vector<string> resn;
    vector<vector<int> > res = residuize(m,resn,heavy(m));
    cout << "Identifying residues: " << endl;
    print_vector(resn);
    int ncts = inter_group_distances(res,m.x.data()).size();
    cout << "Number of contacts: " << ncts << endl;

    ccxtc::xtc f(xtc.c_str());  
    ogzstream o(gz.c_str());
    while(f.next()){
        o << join(" ", inter_group_distances(res,*f.x)) << "\n"; 
    }
    o.close();

    
}
void cts(int argc, char *argv[]){
    bool usage = 1;
    
    string gro;    
    string xtc;    
    string gz;    
    string txt;    

    if( argc>2 ){
        int i=2;
        while(i<(argc-1)){
            if(_(argv[i])=="-s")  gro = _(argv[++i]);
            else if(_(argv[i])=="-f" ) xtc = _(argv[++i]);
            else if(_(argv[i])=="-o" ) gz = _(argv[++i]);
            else if(_(argv[i])=="-i" ) txt = _(argv[++i]);
            i++;
        }
        if(gro.size()==0 || !fopen(gro)) { cout << "ERROR: No GRO file" << endl; }
        else if(xtc.size()==0 || !fopen(xtc)) {cout << "ERROR: No XTC file" << endl; }
        else if(gz.size()==0) {cout << "ERROR: No output file" << endl; }
        else{
            cts(gro,xtc,gz,txt);
            usage = 0;
        }
    }
    
    if(usage){
        cout <<
        "Usage: \n" <<
        "   \n" <<
        "   " << argv[0] << " cts -s file.gro -f file.xtc -o file.gz [-i file.txt]\n" <<
        "\n" <<
        "Description:\n" <<
        "\n" <<
        "   -s file.gro  - input-   GRO file to identify Residues \n" <<
        "   -f file.xtc  - input-   XTC file to use for calculation \n" <<
        "   -o file.gz   -output-   output file (GZIP compressed) \n" <<
        "   -i file.txt  -output-   Index file with the description of the columns \n" <<
        "\n" << endl;
    }
    
}
int main(int argc, char *argv[]){
    {// Parse Options
        bool usage = 1;
        
        if( argc>1){
            if(_(argv[1])=="cts"){
                usage = 0;
                cts(argc,argv);
            }
            if(_(argv[1])=="calc_pca"){
                usage = 0;
                calc_pca(argc,argv);
            }
        }
        
        if(usage){
            cout <<
            "Usage: \n" <<
            "   \n"<<
            "   -> " << argv[0] <<" [cts|pca_calc|pca_do]\n"<<
            "   -> " << argv[0] <<" cts help\n"<<
            "   -> " << argv[0] <<" calc_pca help\n"<<
            "   -> " << argv[0] <<" use_pca help\n"<<
            ""<< endl;
            return 0;
        }
    }
    

    return 0;
}
#endif
