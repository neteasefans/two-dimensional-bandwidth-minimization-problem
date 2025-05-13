#define ILOUSESTL
using namespace std;
#ifndef TYPEDEF_H
#define TYPEDEF_H
//#define _CRT_SECURE_NO_DEPRECATE
#include <iostream>
#include <cstring>
#include <fstream>
#include <cassert>
#include <cmath>
#include <vector>
#include <list>
#include <ctime>
#include <sstream>

//using namespace std;

#define MAXVALUE 99999999
#define EMPTY (-1)

typedef struct MoveCandidate{
    explicit MoveCandidate(int smallVertex=EMPTY, int bigVertex=EMPTY);
    int vertex1;
    int vertex2;
    int objValue;
    int maxEdge;
    double moveValue;
}MoveCandidate;
//输入参数
typedef struct Input{
    explicit Input(int seed=0,int baseT = 0,int randomT=0,double timeLimit=1.0,int improveCut=1000);
    void printInformation(string& out) const;
    //parameter
    int seed{};
    int baseTenure{};
    int randomTenure{};
    double timeLimit{};
    int improveCut{};
    double reward{};
    double penalization{};
    double compensation{};
    double smoothC{};
    double smoothT{};
    int slack{};
    clock_t begin{};
    int stop{};
    int maxIter{};
}Input;

//记录迭代过程中的一些信息
typedef struct Output{
    Output();
    void printInformation(string& out) const;
    void copy(Output a);
    //result
    int iter{};
    int stagnate{};
    int obj{};
    int tie{};
    double moveValue{};
    //run time
    double totalTime{};     // 禁忌搜索总时间
    double findBestTime{};  // 找到最优解所需要的时间
    double searchTime{};
    double probabilityTime{};
    double regionTime{};
    double evaluateTime{};
    double tabuTime{};
    //frequency
    int searchNum{};
    int selectNum{};
    int emptyNum{};
    int tabuNum{};
}Output;

typedef struct PLOutput{
    PLOutput();
    void printInformation(string& out) const;
    void add(Output output);
    int iter{};
    int stagnate{};
    int obj{};
    int maxEdge{};
    int allIter{};
    double moveValue{};
    double totalTime{};
    double findBestTime{};
    double tabuTime{};
}PLOutput;

class InstanceData{
public:
    void readFile(string file);
    void initialization();
    void printGraph(string& out) const;
    ~InstanceData();

    int **graph;
    int **edge;
    int **vertexEdge;
    int **vertexVertex;
    int *degree;
    string insName;
    int vertexNum;
    int edgeNum;
    int hostGraphSize;
    int hostIndexNum;
    int accuracy;
    string graphInformation;
};

class Solution{
public:
    explicit Solution(const InstanceData& ins);
    ~Solution();
    void getLength(const InstanceData& ins);
    void initialingSolSequential(const InstanceData& ins);
    void initialingSolProbability(double ** probability,const InstanceData& ins);
    void initialingSolRandom(const InstanceData& ins);
    void initialingSolFromFile(const string& file,const InstanceData& ins);
    void initialingSolGreedyly(const InstanceData& ins);
    void initialingSolGreedyly_improve(const InstanceData& ins);
    void coverSol(const Solution& sol);
    void execute(MoveCandidate move,const InstanceData& ins);
    void checkSol(const InstanceData& ins) const;
    void printHostGraph(bool judge= false,bool judge2 = false) const;
    void printHostGraph(double ** pro,bool judge= false,bool judge2 = false) const;
    void outputSol(const string& file) const;
    void vogel(double ** probability,vector<int> &sequence) const;
    void roulette(double ** probability,const vector<int>& sequence) const;
    int *injection; //对应关系
    int *hostGraph;
    int **edgeLength;
    int objectValue;
    int maxEdge;
    double moveValue;
    int hostGraphSize;
    int hostIndexNum;
    int edgeNum;
    int vertexNum;
    int accuracy;
};

class TabuTable{
public:
    TabuTable(const InstanceData& ins,int bTenure,int rTenure);
    ~TabuTable();
    void tabu(MoveCandidate move,int iter,const Solution& sol) const;
    [[nodiscard]] bool judge(int vertex1,int vertex2,int iter,const Solution& sol) const;
    void printTabuTable() const;

    int tabuTableSize;
    int baseTenure;
    int randomTenure;
    int** tabuTable;
};

class RoI{
public:
    RoI(int size,int hostGraphSize,int obj,int slack=1);
    ~RoI();
    void updateRegions(int hostGraphSize,int obj,int slack=1) const;
    void intersection(const int* list) const;
    void vertexRegion(int vertex, const InstanceData& ins, const Solution& sol, double &runtime);
    void getBandwidth(const Solution& sol, double &runtime);
    int ** regions;
    int * region;
    int * bandwidth;
    vector<int> regionVertex;
    vector<int> maxEdgeVertex;
    int regionSize;
};

class Length{
public:
    explicit Length(int hostGraphSize,int n);
    ~Length();
    void initialize(const Solution& sol, double &runtime) const;
    void cutLength(const Solution& sol, const InstanceData& ins, int vertex, double &runtime) const;
    void addLength(const Solution& sol,const InstanceData& ins,MoveCandidate move,double &runtime) const;
    void getMoveValue(int accuracy,MoveCandidate &move,double &runtime) const;
    int *lengthNum;
    int *cutLengthNum ;
    int *addLengthNum ;
    int **distance;
    int size;
    int hostIndexNum;
};

class Func{
public:
    static int getDistance(int loc1,int loc2,int loc3);
    static float getDistance_o(int loc1,int loc2,int loc3);
    static int getRouteNum(int ** dis,int v1,int v2,const InstanceData& ins);
    static void updateMoves(vector<MoveCandidate> &moves, MoveCandidate operation);
    static void printMatrix(int** matrix,int row,int col);
    static void printMatrix(double** matrix,int row,int col);
    static void printMatrix(int** matrix,int ** matrix1,int row,int col);
    static void printMove(MoveCandidate move);
    static void printMoves(const vector<MoveCandidate>& moves);
    static void printList(int* list,int length);
    static void printList(const int *list1,const int *list2,int size);
};

class Probability{
public:
    explicit Probability(int hostIndexNum,Input input);
    explicit Probability(InstanceData &ins,int hostIndexNum,Input input);
    ~Probability();
//    void update(const Solution& sol,const MoveCandidate& move,double &runtime) const;
    void update(const Solution& sol0,const Solution& sol1,double &runtime) const;
    void smooth(int vertex,double pro) const;
    void initial();

    double ** probability;
    int size;
    double r{};
    double p{};
    double c{};
    double s{};
    double t{};
};

void searchNeighbor(Solution &sol, InstanceData &ins, const TabuTable& tabuTable,
    vector<MoveCandidate> &swaps, Output & output, RoI& roI, Length& length);
void tabuSearch(Solution &sol, InstanceData &ins, Input input,Output &output, Probability &pro,double baseTime);
void searchNeighbor_NRS(Solution &sol, InstanceData &ins, const TabuTable& tabuTable,
                    vector<MoveCandidate> &swaps, Output & output, RoI& roI, Length& length);
void tabuSearch_NRS(Solution &sol, InstanceData &ins, Input input,Output &output, Probability &pro,double baseTime);
void tabuSearch_nolearning(Solution &sol, InstanceData &ins, Input input,Output &output, double baseTime);
void PLTS(Solution &sol, InstanceData &ins, Input input,Output& inOutput,PLOutput &output);
void NRS(Solution &sol, InstanceData &ins, Input input,Output& inOutput,PLOutput &output);
void ITS(Solution &sol, InstanceData &ins, Input input,Output& inOutput,PLOutput &output);
void GITS(Solution &sol, InstanceData &ins, Input input,Output& inOutput,PLOutput &output);
void program2(const string input_file,const string out_file,const string type,int num);
void program(const string input_file,const string type,int num);
#endif