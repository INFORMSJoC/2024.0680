
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <tuple>
#include <queue>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>


using namespace std;

double objective(int n, int N,int Delta,const std::vector<double>& c1,const std::vector<int>& X1,const std::vector<double>& bestsol){
    double objval = 0;
    for (int i = 0;i<N;i++) {
        objval = objval + c1[i]*(bestsol[i]-X1[i]);
    }
    for (int i = 0; i<n; i++) {
        for (int j=0; j<n; j++) {
            int num = i*n+j;
            if (j<n-1) {
                objval=objval + abs(bestsol[num+1]-bestsol[num]);
            }
            if (num>=n) {
                objval = objval + abs(bestsol[num]-bestsol[num-n]);
            }
        }
    }
    return objval;

}
vector<double> rotate_double(const std::vector<double>& c,int n,int N){
    vector<double> rotated(N,0);
    for (int i = 0; i<n;i++){
        for (int j=0;j<n;j++) {
            rotated[j*n+(n-1-i)] = c[i*n+j];
        }
    }
    return rotated;
}
vector<int> rotate(const std::vector<int>& X,int n,int N) {
    vector<int> rotated(N,0);
    for (int i = 0; i<n;i++){
        for (int j=0;j<n;j++) {
            rotated[j*n+(n-1-i)] = X[i*n+j];
        }
    }
    return rotated;
}


vector<double> improvepath(
                     int n, int N, int Delta,
                     int m,const std::vector<double>& c1,
                     const std::vector<int>& X1, const std::vector<double>& X2,
                     const std::vector<int>& V, vector<double>& costvec1, vector<int>& prevvec1, vector<int>& valuevec1) {

    int d;
    int usedcap;
    int index;
    int predindex;
    int dabs;
    double bestcost;
    int Deltaplus1 = Delta+1;
    double finalcost =  std::numeric_limits<double>::max();
    int bestvert=0;
    double nodecost;
    double newcost;
    double precost;
  
    
    std::fill(costvec1.begin(),costvec1.end(), std::numeric_limits<double>::max());

    std::fill(prevvec1.begin(),prevvec1.end(), -1);

    std::fill(valuevec1.begin(),valuevec1.end(), 0);
    //vector<int> valuevec(nmalm*(Deltaplus1), 0);

    
    
    for (int o = 0; o < m; o++) {
        d = V[o] - X1[0];
        usedcap = abs(d);
        if (usedcap <= Delta) {
            index = o * Deltaplus1 + usedcap;
            costvec1[index] = c1[0] * d + abs(X2[n] - V[o]);
            valuevec1[index] = o;
        }
    }
    for (int k = 1; k < N; k++){
        int prevk = k-1;
        if (k%n==0) {
            prevk = prevk-n;
        }
        for (int o = 0; o < m; o++) {
            for (int p = 0; p < m; p++) {
                d = V[p] - X1[k];
                precost = c1[k]*d;
                if (k-prevk==1) {
                    precost = precost + abs(V[p]-V[o]);
                }
              
                if (k>=n) {
                    precost = precost + abs(V[p]- X2[k-n]);
         
                }
                if (k+n<N) {
                    precost = precost + abs(V[p]-X2[k+n]);
    
                }
                dabs = abs(V[p]-X1[k]);
                index = (k*m+p)*Deltaplus1+dabs;
                
                predindex = (prevk*m+o)*Deltaplus1;
                bestcost = std::numeric_limits<double>::max();
                for (int r = Delta-dabs; r >= 0; r--) {
                    nodecost = costvec1[predindex];
                    if (nodecost<bestcost) {
                        //bestcost = nodecost;
                        newcost = nodecost + precost;
                        if ( newcost<costvec1[index]) {
                            costvec1[index] = newcost;
                            prevvec1[index] = predindex;
                            valuevec1[index] = p;
                            if (k==N-1-n) {
                                if (newcost<finalcost) {     
                                    finalcost = newcost;
                                    bestvert = index;
                                }
                            }
                        }
                    }
                    predindex++;
                    index++;
                }
            }
            
        }
        if (k%n==n-1) { // ensures that only every other row is considered
            k=k+n;
        }
    }
    std::vector<double> result;
    result.reserve(N);
    int p;

    for (int j=n-1; j>=0; j--) {
        for (int i=n; i>0;i--) {
            int num = i+j*n;
            if (j%2==0) {
                
                //cout<<"bestvert"<<bestvert<<endl;
                p = valuevec1[bestvert];

                result.push_back(V[p]);
                
                bestvert = prevvec1[bestvert];
            }
            else {

                result.push_back(X2[num-1]);
            }
            
        }
    }
	std::reverse(result.begin(), result.end());
    for (int i=1;i<n; i=i+2) {
        for (int j=0;j<n;j++) {
            finalcost = finalcost + c1[i*n+j]*(result[i*n+j]-X1[i*n+j]);
            if (j>0){
                finalcost = finalcost+abs(result[i*n+j]-result[i*n+j-1]);
            }
        }
    }
    //cout<<"fc"<<finalcost<<endl;

    return result;
}

vector<double> top_improve(
                     int n, int N, int Delta,
                     int m, std::vector<double>& c1,
                     std::vector<int>& X1, std::vector<double>& bestsol,
                     const std::vector<int>& V) {
    bool improvement = true;
    int lastimprove =0;
    double bestval = objective(n,N,Delta,c1,X1,bestsol);
    double newobj;
    double remcap;
    int myloopcount=0;
    vector<double> costvec1(N*m*(Delta+1), std::numeric_limits<double>::max());
    
    vector<int> prevvec1(N*m*(Delta+1), -1);
    
    vector<int> valuevec1(N*m*(Delta+1), 0);
    while (improvement and myloopcount<25) {
        //myloopcount++;
        improvement=false;
        if (lastimprove!=1){
            remcap = Delta;
            for (int i=1; i<n; i=i+2) {
                for (int j=0;j<n; j++) {
                    remcap=remcap-abs(X1[i*n+j] -bestsol[i*n+j]);
                }
            }
            remcap = std::floor(remcap);
            bestsol = improvepath(n,N,remcap,m,c1,X1,bestsol,V,costvec1,prevvec1,valuevec1);
            newobj = objective(n,N,Delta,c1,X1,bestsol);
            if (newobj< bestval) {
                //cout<<"improved"<<bestval<<","<<newobj<<endl;
                bestval=newobj;
                improvement =true;
                lastimprove = 1;
                
            }
        } else {
            break;
        }
        
        X1 = rotate(X1,n,N);
        c1 = rotate_double(c1,n,N);
        bestsol = rotate_double(bestsol,n,N);
        if (lastimprove!=2) {
            remcap = Delta;
            for (int i=1; i<n; i=i+2) {
                for (int j=0;j<n; j++) {
                    remcap=remcap-abs(X1[i*n+j] -bestsol[i*n+j]);
                }
            }
            remcap = std::floor(remcap);
            bestsol = improvepath(n,N,remcap,m,c1,X1,bestsol,V,costvec1,prevvec1,valuevec1);
            newobj = objective(n,N,Delta,c1,X1,bestsol);
            if (newobj< bestval) {
                //cout<<"improved"<<bestval<<","<<newobj<<endl;
                bestval=newobj;
                improvement =true;
                lastimprove = 2;
            }
        } else {
            break;
        }
        
        X1 = rotate(X1,n,N);
        c1 = rotate_double(c1,n,N);
        bestsol = rotate_double(bestsol,n,N);
        if (lastimprove!=3) {
            remcap = Delta;
            for (int i=1; i<n; i=i+2) {
                for (int j=0;j<n; j++) {
                    remcap=remcap-abs(X1[i*n+j] -bestsol[i*n+j]);
                }
            }
            remcap = std::floor(remcap);
            bestsol = improvepath(n,N,remcap,m,c1,X1,bestsol,V,costvec1,prevvec1,valuevec1);
            newobj = objective(n,N,Delta,c1,X1,bestsol);
            if (newobj< bestval) {
                //cout<<"improved"<<bestval<<","<<newobj<<endl;
                bestval=newobj;
                improvement =true;
                lastimprove = 3;
                
            }
        } else {
            break;
        }
        X1 = rotate(X1,n,N);
        c1 = rotate_double(c1,n,N);
        bestsol = rotate_double(bestsol,n,N);
        if (lastimprove!= 4) {
        
            remcap = Delta;
            for (int i=1; i<n; i=i+2) {
                for (int j=0;j<n; j++) {
                    remcap=remcap-abs(X1[i*n+j] -bestsol[i*n+j]);
                }
            }
            remcap = std::floor(remcap);
            bestsol = improvepath(n,N,remcap,m,c1,X1,bestsol,V,costvec1,prevvec1,valuevec1);
            newobj = objective(n,N,Delta,c1,X1,bestsol);
            if (newobj< bestval) {
                //cout<<"improved"<<bestval<<","<<newobj<<endl;
                bestval=newobj;
                improvement =true;
                lastimprove = 4;
                
            }
        } else {
            break;
        }
        X1 = rotate(X1,n,N);
        c1 = rotate_double(c1,n,N);
        bestsol = rotate_double(bestsol,n,N);

    }
    return bestsol;
}



PYBIND11_MODULE(top_improve, m) {
    m.doc() = "pybind11 example plugin"; 
    m.def("top_improve", &top_improve, "top_improve");
}
