#include "TypeDef.h"
#define IFPRINT 1


int main(int argc, char* argv[]) {
    string input_file = argv[argc - 4];
	string out_file = argv[argc - 1];
    int num = atoi(argv[argc - 2]);
    int repeat = atoi(argv[argc - 3]);;
    string type = "RLTS";
    for (int i = 0;i < repeat;++i) {
		program2(input_file,out_file,type, (num-1)*1 + i +1);
    }
    return 0;
}


int main2(int argc, char* argv[]) {
    /*测试贪婪初始化*/ 
    // string input_file = "../ins/grids/mesh33_33.txt "; 
    // string input_file = "../ins/harwell-boeing/gr_30_30 ";
    // string input_file = "../ins/gr_3_3.txt "; 
    string input_file = "../ins/new_ins/ibm32.txt ";
    InstanceData ins{};
    ins.readFile(input_file);
    Solution sol(ins);
    sol.initialingSolGreedyly_improve(ins);
    sol.printHostGraph();
    return 0;
}


int main1(int argc, char* argv[]) {
    // string input_file = "../ins/gr_3_3.txt "; 
    // string input_file = "../ins/grids/mesh33_33.txt "; 
    string input_file = "../ins/harwell-boeing/gr_30_30 ";
    // string input_file = "../ins/regular/k10.mtx ";
    // string input_file = "../ins/new_ins/ibm32.txt ";
    int num = 66;
    int repeat = 1;
    string type = "RLTS";
    for (int i = 0;i < repeat;++i) {
        program(input_file,type, num * 100 + i);
    }
    return 0;
}

int main0(int argc, char* argv[]) {
    string input_file = argv[argc - 2];
    int num = atoi(argv[argc - 1]);
    int repeat;
    if (num == -1){
        repeat = 2;
    }else{
        repeat = 80;
    }
    // else if (num<9){
    //     repeat = 80;
    // }else{
    //     repeat = 40;
    // }
    string type = "ITS";
    for (int i = 0;i < repeat;++i) {
        int tt = i % 4;
        int ttt = i / 4;
        if (tt==0){
            type = "NRS";
        }else if (tt==1){
            type = "ITS";
        }else if (tt==2){
            type = "GITS";
        }else if (tt==3){
            type = "RLTS";
        }
	//if (ttt!=4){
	program(input_file,type, (num-1)*1 + ttt+1);
	//}else{
	//	program(input_file,type, 20);
	//}
    }
    return 0;
}

void program2(const string input_file,const string out_file,const string type,int num) {
    //读取算例，生成初始解
    clock_t begin = clock();
    
    InstanceData ins{};
    ins.readFile(input_file);
    Solution sol(ins);
    Input input;
    PLOutput output;
    Output inOutput;
    //输入参数-终止条件
    // sol.initialingSolGreedyly_improve(ins);
    input.seed = num; //1100;
    srand(input.seed);
    if (ins.vertexNum>=500){
        input.timeLimit = 500.0;//1.0*ins.vertexNum/2;
    }else if (ins.vertexNum>=0){
        input.timeLimit = 100.0;//1.0*ins.vertexNum/2;
    }else{
        input.timeLimit = 20.0;//1.0*ins.vertexNum/2;
    }
    input.improveCut = 2000;//10000 - 500*((num-1)%20);
    input.stop = 100;
    input.maxIter = 5000;
    input.begin = begin;
    //输入参数
    int neighborSize = ins.hostIndexNum;
    double meanTenure = 0.75;
    double devTenure = 0.5;
    input.baseTenure = 1 + int(meanTenure * (1 - devTenure) * neighborSize);
    input.randomTenure = 1 + int(2 * meanTenure * devTenure * neighborSize);
    input.reward = 0.4;
    input.penalization = 0.2;
    input.compensation = 0.3;
    input.smoothC = 0.7;
    input.smoothT = 0.9;
    string outInfo;
    //进行迭代
    if (type=="ITS"){
        ITS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   ITS\n";
    }else if (type=="NRS"){
        NRS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   NRS\n";
    }else if (type=="GITS"){
        GITS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   GITS\n";
    }else if (type=="RLTS"){
        PLTS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   PLTS\n";
    }else{
        printf("No define type!");
    }

    ins.printGraph(outInfo);
    input.printInformation(outInfo);
    output.printInformation(outInfo);
    inOutput.printInformation(outInfo);
    // printf("The improveCut :%d\n",input.improveCut);

    //检查并输出解
    sol.checkSol(ins);
    time_t now = time(NULL);
    tm* time_tm = localtime(&now);
    char strTime[25] = { 0 };
    sprintf(strTime, "%d%02d%02d%02d%02d%02d", time_tm->tm_year + 1900,
        time_tm->tm_mon + 1, time_tm->tm_mday, time_tm->tm_hour,
        time_tm->tm_min, time_tm->tm_sec);
    outInfo += "nowTime  :";
    outInfo += strTime;
    outInfo += "\n\n";
    for (int i = 0;i < sol.vertexNum;++i) {
        outInfo += to_string(sol.injection[i]) + " ";
    }
    outInfo += "\n---------------------\n";
    cout << outInfo << endl;
    string p = out_file;//"./" + ins.insName + "_" + type + "_" + to_string(num) + ".out";
    ofstream dataFile;
    dataFile.open(p);// 追加 ofstream::app
    if (!dataFile.is_open()) {
        fprintf(stderr, "Can not open file %s\n", p.c_str());
        exit(10);
    }
    dataFile << outInfo;
    dataFile.close();
}

void program(const string input_file,const string type,int num) {
    //读取算例，生成初始解
    clock_t begin = clock();
    
    InstanceData ins{};
    ins.readFile(input_file);
    Solution sol(ins);
    Input input;
    PLOutput output;
    Output inOutput;
    //输入参数-终止条件
    // sol.initialingSolGreedyly_improve(ins);
    input.seed = num; //1100;
    srand(input.seed);
    if (ins.vertexNum>=500){
        input.timeLimit = 500.0;//1.0*ins.vertexNum/2;
    }else if (ins.vertexNum>=0){
        input.timeLimit = 100.0;//1.0*ins.vertexNum/2;
    }else{
        input.timeLimit = 20.0;//1.0*ins.vertexNum/2;
    }
    input.improveCut = 2000;//10000 - 500*((num-1)%20);
    input.stop = 100;
    input.maxIter = 5000;
    input.begin = begin;
    //输入参数
    int neighborSize = ins.hostIndexNum;
    double meanTenure = 0.75;
    double devTenure = 0.5;
    input.baseTenure = 1 + int(meanTenure * (1 - devTenure) * neighborSize);
    input.randomTenure = 1 + int(2 * meanTenure * devTenure * neighborSize);
    input.reward = 0.4;
    input.penalization = 0.2;
    input.compensation = 0.3;
    input.smoothC = 0.7;
    input.smoothT = 0.9;
    string outInfo;
    //进行迭代
    if (type=="ITS"){
        ITS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   ITS\n";
    }else if (type=="NRS"){
        NRS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   NRS\n";
    }else if (type=="GITS"){
        GITS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   GITS\n";
    }else if (type=="RLTS"){
        PLTS(sol, ins, input, inOutput, output);
        outInfo = "exprType  :   PLTS\n";
    }else{
        printf("No define type!");
    }

    ins.printGraph(outInfo);
    input.printInformation(outInfo);
    output.printInformation(outInfo);
    inOutput.printInformation(outInfo);
    // printf("The improveCut :%d\n",input.improveCut);

    //检查并输出解
    sol.checkSol(ins);
    time_t now = time(NULL);
    tm* time_tm = localtime(&now);
    char strTime[25] = { 0 };
    sprintf(strTime, "%d%02d%02d%02d%02d%02d", time_tm->tm_year + 1900,
        time_tm->tm_mon + 1, time_tm->tm_mday, time_tm->tm_hour,
        time_tm->tm_min, time_tm->tm_sec);
    outInfo += "nowTime  :";
    outInfo += strTime;
    outInfo += "\n\n";
    for (int i = 0;i < sol.vertexNum;++i) {
        outInfo += to_string(sol.injection[i]) + " ";
    }
    outInfo += "\n---------------------\n";
    cout << outInfo << endl;
    string p = "./" + ins.insName + "_" + type + "_" + to_string(num) + ".out";
    ofstream dataFile;
    dataFile.open(p);// 追加 ofstream::app
    if (!dataFile.is_open()) {
        fprintf(stderr, "Can not open file %s\n", p.c_str());
        exit(10);
    }
    dataFile << outInfo;
    dataFile.close();
}


template <class T>
void change(T &a,T &b){
    T ch;
    ch = a;
    a = b;
    b = ch;
}

//构造函数
MoveCandidate::MoveCandidate(int smallVertex, int bigVertex) {
    vertex1 = smallVertex;
    vertex2 = bigVertex;
    objValue = MAXVALUE;
    maxEdge = MAXVALUE;
    moveValue = MAXVALUE;
}

//构造函数
Output::Output() {
    obj = MAXVALUE;
    tie = MAXVALUE;
}
//输出属性
void Output::printInformation(string& out) const {
    //times
//    printf("The region time      : %.6f / %.4f%%\n",regionTime, 100.0 * regionTime / totalTime);
//    printf("The evaluate time    : %.6f / %.4f%%\n",evaluateTime,100.0*evaluateTime/totalTime);
//    printf("The probability time : %.6f / %.4f%%\n",probabilityTime,100.0*probabilityTime/totalTime);
//    printf("Search Neighbor time : %.6f / %.4f%%\n",searchTime,100.0*searchTime/totalTime);

    //frequency
//    printf("The num of search solution :%d\n",searchNum);
//    printf("The num of select solution :%d\n",selectNum);
//    printf("The emptyNum        :%d\n",emptyNum);
//    printf("The num of tabu   solution :%d\n",tabuNum);
//    printf("The num of region solution :%d\n",regionSize);

    //results
    //printf("The inTotalTime       : %.6f \n", totalTime);
    //printf("The inFindBestTime    : %.6f \n", findBestTime);
    //printf("The inIterations      : %d\n", iter);
    //printf("The inStagnate        : %d\n", stagnate);

    out += "The inTotalTime       : " + to_string(totalTime) + "\n";
    out += "The inFindBestTime    : " + to_string(findBestTime) + "\n";
    out += "The inIterations      : " + to_string(iter) + "\n";
    out += "The inStagnate        : " + to_string(stagnate) + "\n";
//    printf("The solution values  : %d\n",obj);
//    printf("The maxEdge value        : %d\n",maxEdge);
//    printf("The move value       : %.12f\n",moveValue);

}

void Output::copy(Output a) {
    iter = a.iter;
    stagnate = a.stagnate;
    obj = a.obj;
    tie = a.tie;
    moveValue = a.moveValue;
    totalTime = a.totalTime;
    findBestTime = a.findBestTime;
    searchTime = a.searchTime;
    probabilityTime = a.probabilityTime;
    regionTime = a.regionTime;
    evaluateTime = a.evaluateTime;
    tabuTime = a.tabuTime;
    searchNum = a.searchNum;
    selectNum = a.selectNum;
    emptyNum = a.emptyNum;
    tabuNum = a.tabuNum;
}

//读取文件
void InstanceData::readFile(string input_file) {
    unsigned int loc = input_file.find_last_of('/') + 1;
    while(input_file[loc]!='\0'&&input_file[loc]!=' '){
        insName.push_back(input_file[loc]);
        loc++;
    }
    const char *file = input_file.c_str();
    ifstream inputFile(file);
    if (!inputFile) {
        fprintf(stderr, "Can not open file %s\n", file);
        exit(10);
    }
    getline(inputFile,graphInformation);
    inputFile >> vertexNum >> vertexNum >> edgeNum;
//    int ttt = (int)sqrt(vertexNum);
//    edgeNum = 2*ttt*(ttt-1);//grids
    edge = new int * [edgeNum];
    int edgeIndex = 0;
    string tempStr;
    int vertex1,vertex2;
    while (getline(inputFile,tempStr)){
        if (tempStr.size()>1){
            sscanf(tempStr.c_str(), "%d %d", &vertex1,&vertex2);
            assert(vertex1>0&&vertex1 <= vertexNum&&vertex2 > 0&&vertex2 <= vertexNum&&vertex1 != vertex2);
            if (vertex1>vertex2)
                edge[edgeIndex] = new int[2]{vertex1-1,vertex2-1};
            else
                edge[edgeIndex] = new int[2]{vertex2-1,vertex1-1};
            edgeIndex++;
        }
    }
    inputFile.close();
    assert(edgeIndex==edgeNum);
    initialization();
}
void InstanceData::initialization(){
    hostGraphSize = (int )sqrt(vertexNum-1) + 1;
    hostIndexNum = hostGraphSize * hostGraphSize;
    degree = new int [hostIndexNum];
    memset(degree,0,sizeof(int)*hostIndexNum);
    graph = new int * [hostIndexNum];
    vertexEdge = new int * [hostIndexNum];
    vertexVertex = new int * [hostIndexNum];
    for (int i = 0; i < hostIndexNum; ++i) {
        graph[i] = new int [hostIndexNum];
        memset(graph[i],EMPTY,sizeof(int)*hostIndexNum);
        vertexEdge[i] = new int [vertexNum];
        memset(vertexEdge[i],EMPTY,sizeof(int)*vertexNum);
        vertexVertex[i] = new int [vertexNum];
        memset(vertexVertex[i],EMPTY,sizeof(int)*vertexNum);
    }
    for (int i = 0; i < edgeNum; ++i) {
        int vertex1 = edge[i][0];
        int vertex2 = edge[i][1];
        graph[vertex1][vertex2] = graph[vertex2][vertex1] = 1;
        vertexVertex[vertex1][degree[vertex1]] = vertex2;
        vertexVertex[vertex2][degree[vertex2]] = vertex1;
        vertexEdge[vertex1][degree[vertex1]] = vertexEdge[vertex2][degree[vertex2]] = i;
        degree[vertex1]++;
        degree[vertex2]++;
    }
    accuracy = edgeNum;
}
//析构函数
InstanceData::~InstanceData() {
    for (int i = 0; i < hostIndexNum; i++) {
        delete[] graph[i];
        delete[] vertexEdge[i];
        delete[] vertexVertex[i];
    }
    for (int i=0;i<edgeNum;i++)
        delete[] edge[i];
    delete[] degree;
    delete[] graph;
    delete[] vertexEdge;
    delete[] vertexVertex;
    delete[] edge;
//    cout<<"The instance has been delete!"<<endl;
}
//输出算例信息
void InstanceData::printGraph(string& out) const {
    //cout << "The name            : " << insName << endl;
    //cout << "The vertexNum       : " << vertexNum << endl;
    //cout << "The EdgeNum         : " << edgeNum << endl;
    out += "The name            : " + insName + "\n";
    out += "The vertexNum       : " + to_string(vertexNum) + "\n";
    out += "The EdgeNum         : " + to_string(edgeNum) + "\n";
}


//构造函数 对solution进行初始化
Solution::Solution(const InstanceData& ins) {
    hostGraphSize = ins.hostGraphSize;
    hostIndexNum = ins.hostIndexNum;
    edgeNum = ins.edgeNum;
    vertexNum = ins.vertexNum;
    accuracy = ins.accuracy;
    objectValue = MAXVALUE;
    maxEdge = MAXVALUE;
    moveValue = MAXVALUE;
    injection = new int [hostIndexNum];
    hostGraph = new int [hostIndexNum];
    edgeLength = new int * [edgeNum];
    memset(injection,EMPTY, sizeof(int) * hostIndexNum);
    memset(hostGraph,EMPTY, sizeof(int) * hostIndexNum);
    for (int i=0;i<edgeNum;i++){
        edgeLength[i] = new int [3]{0};
        edgeLength[i][0] = ins.edge[i][0];
        edgeLength[i][1] = ins.edge[i][1];
    }
}
//析构函数
Solution::~Solution() {
    for (int i = 0;i < edgeNum;i++)
        delete[] edgeLength[i];
    delete[] injection;
    delete[] hostGraph;
    delete[] edgeLength;
//    cout<<"The solution has been delete!"<<endl;
}
//getLength() 可以对solution的边长整体更新，并且更新objectValue
void Solution::getLength(const InstanceData& ins) {
    objectValue = 0;
    maxEdge = 0;
    for (int i=0;i<edgeNum;i++){
        edgeLength[i][2] = Func::getDistance(injection[ins.edge[i][0]], injection[ins.edge[i][1]], hostGraphSize);
        if (edgeLength[i][2]>objectValue){
            objectValue = edgeLength[i][2];
            maxEdge = 1;
        }else if (edgeLength[i][2]==objectValue){
            maxEdge++;
        }
    }
    moveValue = 1.0*objectValue - 1.0*maxEdge/accuracy;
}
//顺序生成初始解
void Solution::initialingSolSequential(const InstanceData& ins) {
    for (int i = 0; i < hostIndexNum; i++)
        injection[i] = hostGraph[i] = i;
    getLength(ins);
}
//随机生成初始解
void Solution::initialingSolRandom(const InstanceData& ins){
    int *list;
    list = new int [hostIndexNum];
    for (int i = 0; i < hostIndexNum; i++)
        list[i] = i;
    int index = hostIndexNum;
    int randomIndex;
    for (int i = 0; i < hostIndexNum; i++) {
        randomIndex = rand()%index;
        hostGraph[list[randomIndex]] = i;
        injection[i] = list[randomIndex];
        list[randomIndex] = list[index-1];
        index--;
    }
    delete[] list;
    getLength(ins);
}
//从文件中生成初始解
void Solution::initialingSolFromFile(const string& file,const InstanceData& ins) {
    ifstream inputFile(file);
    if (!inputFile) {
        fprintf(stderr, "Can not open file %s\n", file.c_str());
        exit(10);
    }
    int fileVertex,fileEdge,fileObj;
    string temp = "0";
    while(temp[temp.size()-1]!=')'){
        inputFile >> temp;
    }
    inputFile >> fileVertex >> fileEdge >> fileObj;
    assert(fileVertex==vertexNum && fileEdge==edgeNum);
    string tempStr;
    while (getline(inputFile,tempStr)) {
        if (tempStr.size()>1) {
            int vertex, index;
            sscanf(tempStr.c_str(), "%d : %d", &vertex, &index);
            assert(vertex >= 0 && vertex < vertexNum);
            assert(index >= 0 && index < hostIndexNum);
            injection[vertex] = index;
            hostGraph[index] = vertex;
        }
    }
    getLength(ins);
    assert(fileObj==objectValue);
}

// Greedy initialing
void Solution::initialingSolGreedyly(const InstanceData& ins){
    memset(injection,EMPTY, sizeof(int) * hostIndexNum);
    memset(hostGraph,EMPTY, sizeof(int) * hostIndexNum);
    int max_vertex = rand() % hostIndexNum;
    int min_index  = rand() % hostIndexNum;
    injection[max_vertex] = min_index;
    hostGraph[min_index]  = max_vertex;
    int stay_vertex_num = hostIndexNum - 1;
    while (stay_vertex_num > 0){
        // g1:Find vertex
        float w1 = 0.5;
        float w2 = 0.5;
        float max_vertex_value = -1.0*MAXVALUE;
        max_vertex = EMPTY;
        for (int vertex=0;vertex<hostIndexNum;++vertex){
            if (injection[vertex]!=EMPTY){
                continue;
            }
            int A_num = 0;
            int B_num = 0;
            for (int vertex_v_i=0;vertex_v_i<ins.degree[vertex];++vertex_v_i){
                int ad_index = injection[ins.vertexVertex[vertex][vertex_v_i]];
                if (ad_index==EMPTY){
                    B_num ++;
                }else{
                    A_num ++;
                }
            }
            float vertex_value = w1*A_num - w2*B_num;
            if (vertex_value>max_vertex_value){
                max_vertex_value = vertex_value;
                max_vertex = vertex;
            // }else if (vertex_value==max_vertex_value){
            //     int ra = rand()%10;
            //     if (ra<5){
            //         max_vertex = vertex;
            //     }
            }
        }
        // g2:Find index
        min_index = EMPTY;
        int min_index_value = hostIndexNum;
        for (int index=0;index<hostIndexNum;++index){
            if (hostGraph[index]!=EMPTY){
                continue;
            }
            int index_value = -1;
            for (int vertex_v_i=0;vertex_v_i<ins.degree[max_vertex];++vertex_v_i){
                int ad_index = injection[ins.vertexVertex[max_vertex][vertex_v_i]];
                if (ad_index==EMPTY){
                    continue;
                }
                int dis = Func::getDistance(index,ad_index,hostGraphSize);
                if (dis > index_value){
                    index_value = dis;
                }
            }
            if (index_value < min_index_value){
                min_index_value = index_value;
                min_index = index;
            // }else if (index_value < min_index_value){
            //     int ra = rand()%10;
            //     if (ra<5){
            //         min_index = index;
            //     }
            }
        }
        injection[max_vertex] = min_index;
        hostGraph[min_index]  = max_vertex;
        stay_vertex_num --;
    }
    getLength(ins);
}


void Solution::initialingSolGreedyly_improve(const InstanceData& ins){
    memset(injection,EMPTY, sizeof(int) * hostIndexNum);
    memset(hostGraph,EMPTY, sizeof(int) * hostIndexNum);
    // step1: 计算每个点与点之间的距离
    int ** distance = new int * [hostIndexNum] ;
    for (int vertex=0;vertex<hostIndexNum;++vertex){
        distance[vertex] = new int [hostIndexNum];
        memset(distance[vertex],EMPTY,sizeof(int)*hostIndexNum);
        distance[vertex][vertex] = 0;
    }
    for (int edge_index=0;edge_index<ins.edgeNum;++edge_index){
        distance[ins.edge[edge_index][0]][ins.edge[edge_index][1]] = 1;
        distance[ins.edge[edge_index][1]][ins.edge[edge_index][0]] = 1;
    }
    for (int vertex1=0;vertex1<hostIndexNum;++vertex1){
        for (int vertex2=vertex1+1;vertex2<hostIndexNum;++vertex2){
            if (distance[vertex1][vertex2]<1) continue;
            for (int vertex3=0;vertex3<hostIndexNum;++vertex3){
                if (distance[vertex2][vertex3]>0 && distance[vertex1][vertex3]<0){
                    distance[vertex1][vertex3] = distance[vertex1][vertex2] + distance[vertex2][vertex3];
                    distance[vertex3][vertex1] = distance[vertex1][vertex2] + distance[vertex2][vertex3];
                }
            }
        }
    }
    // step2: 找到最长距离且度最小的点,距离相同选择度和最小的点,度和相同选择到达方案最少的点
    int * max_dis_vertex = new int [ins.vertexNum*2];
    int max_num = 0;
    int max_dis = EMPTY;
    int min_degree = 2*ins.vertexNum + 1;
    int min_routeNum = ins.vertexNum;
    for (int vertex1=0;vertex1<hostIndexNum;++vertex1){
        for (int vertex2=vertex1;vertex2<hostIndexNum;++vertex2){
            if (distance[vertex1][vertex2]>max_dis){
                max_dis = distance[vertex1][vertex2];
                max_num = 0;
                max_dis_vertex[max_num] = vertex1;
                max_dis_vertex[max_num+1] = vertex2;
                max_num += 2;
            }else if (distance[vertex1][vertex2]==max_dis){
                int sum_degree = ins.degree[vertex1] + ins.degree[vertex2];
                if (sum_degree==min_degree){
                    int rn = Func::getRouteNum(distance,vertex1,vertex2,ins);
                    if (rn == min_routeNum){
                        max_dis_vertex[max_num] = vertex1;
                        max_dis_vertex[max_num+1] = vertex2;
                        max_num += 2;
                    }else if (rn < min_routeNum){
                        min_routeNum = rn;
                        max_num = 0;
                        max_dis_vertex[max_num] = vertex1;
                        max_dis_vertex[max_num+1] = vertex2;
                        max_num += 2;
                    }
                }else if (sum_degree<min_degree){
                    max_num = 0;
                    min_degree = sum_degree;
                    max_dis_vertex[max_num] = vertex1;
                    max_dis_vertex[max_num+1] = vertex2;
                    max_num += 2;
                }

            }
        }
    }
    // for (int vertex=0;vertex<ins.vertexNum;++vertex){
    //     delete[] distance[vertex];
    //     distance[vertex] = nullptr;
    // }
    // delete[] distance;
    // distance = nullptr;
    // step3：为最前面的四个点分配host的四个角
    if (max_num>=4){
        if (max_dis_vertex[2]!=max_dis_vertex[0] && max_dis_vertex[2]!=max_dis_vertex[1] & max_dis_vertex[3]!=max_dis_vertex[0] && max_dis_vertex[3]!=max_dis_vertex[1]){
            max_num = 4;
        }else{
            max_num = 2;
        }
    }
    for (int line=0;line<max_num;++line){
        int be_index;
        if (line==0){
            be_index = 0;
        } else if (line == 1){
            be_index = ins.hostIndexNum - 1;
        } else if (line == 2){
            be_index = ins.hostGraphSize - 1;
        } else if (line == 3){
            be_index = ins.hostIndexNum - hostGraphSize;
        }
        printf(" %5d -> %5d \n",max_dis_vertex[line],be_index);
        injection[max_dis_vertex[line]] = be_index;
        hostGraph[be_index] = max_dis_vertex[line];
    }
    delete[] max_dis_vertex;
    max_dis_vertex = nullptr;
    int max_vertex;
    int min_index;
    int stay_vertex_num = hostIndexNum - max_num;
    while (stay_vertex_num > 0){
        // g1:Find vertex
        float w1 = 0.8;
        float w2 = 0.2;
        float max_vertex_value = -1.0*MAXVALUE;
        max_vertex = EMPTY;
        for (int vertex=0;vertex<hostIndexNum;++vertex){
            if (injection[vertex]!=EMPTY){
                continue;
            }
            int A_num = 0;
            int B_num = 0;
            for (int vertex_v_i=0;vertex_v_i<ins.degree[vertex];++vertex_v_i){
                int ad_index = injection[ins.vertexVertex[vertex][vertex_v_i]];
                if (ad_index==EMPTY){
                    B_num ++;
                }else{
                    A_num ++;
                }
            }
            float vertex_value = w1*A_num - w2*B_num;
            if (vertex_value>max_vertex_value){
                max_vertex_value = vertex_value;
                max_vertex = vertex;
            }
        }
        // g2:Find index
        // min_index = EMPTY;
        // int min_index_value = hostIndexNum;
        // for (int index=0;index<hostIndexNum;++index){
        //     if (hostGraph[index]!=EMPTY){
        //         continue;
        //     }
        //     int index_value = -1;
        //     for (int vertex_v_i=0;vertex_v_i<ins.degree[max_vertex];++vertex_v_i){
        //         int ad_index = injection[ins.vertexVertex[max_vertex][vertex_v_i]];
        //         if (ad_index==EMPTY){
        //             continue;
        //         }
        //         int dis = Func::getDistance(index,ad_index,hostGraphSize);
        //         if (dis > index_value){
        //             index_value = dis;
        //         }
        //     }
        //     if (index_value < min_index_value){
        //         min_index_value = index_value;
        //         min_index = index;
        //     }
        // }
        min_index = EMPTY;
        float min_index_value = MAXVALUE;
        for (int index=0;index<hostIndexNum;++index){
            if (hostGraph[index]!=EMPTY){
                continue;
            }
            float index_value = 0;
            int length = 1;
            int mi = 1;
            while (length < hostGraphSize){

                int size = 0;
                for (int vertex_v_i=0;vertex_v_i<hostIndexNum;++vertex_v_i){
                    int ad_index = injection[vertex_v_i];
                    if (ad_index==EMPTY || distance[max_vertex][vertex_v_i]!=length || vertex_v_i==max_vertex){
                        continue;
                    }
                    size += 1;
                    float dis = Func::getDistance(index,ad_index,hostGraphSize);
                    for (int line=0;line<mi-1;++line){
                        dis = 1.0*dis / (hostGraphSize-1) / 2;
                    }
                    index_value += dis ;
                }
                if (size>0){
                    mi += 1;
                }
                length += 1;

                if (index_value > min_index_value) break;
            }

            if (index_value < min_index_value){
                min_index_value = index_value;
                min_index = index;
            }
        }
        // printf(" %5d -> %5d (%3.5f)\n",max_vertex,min_index,min_index_value);
        injection[max_vertex] = min_index;
        hostGraph[min_index]  = max_vertex;
        stay_vertex_num --;
    }
    for (int vertex=0;vertex<ins.vertexNum;++vertex){
        delete[] distance[vertex];
        distance[vertex] = nullptr;
    }
    delete[] distance;
    distance = nullptr;
    getLength(ins);
}


void Solution::vogel(double **probability, vector<int> &sequence) const{
    //伏尔格法确定顺序
    //计算最大元素与此大元素的差值
    vector<double> deltaP(hostIndexNum);
    for (int vertex = 0; vertex < hostIndexNum; ++vertex) {
        double first = 0;
        double second = 0;
        for (int index = 0; index < hostIndexNum; ++index) {
            double temp = probability[vertex][index];
            if (temp>=first){
                second = first;
                first = temp;
            }else if (temp>second){
                second = temp;
            }
        }
        deltaP[vertex] = first - second;
    }
    for (int i = 0; i < deltaP.size(); ++i)
        sequence.push_back(i);
    //进行冒泡排序
    for (int i=1; i < deltaP.size(); ++i){
        for (int j = 0; j < deltaP.size() - i; ++j) {
            if (deltaP[j]<deltaP[j+1]){
                change(deltaP[j],deltaP[j+1]);
                change(sequence[j],sequence[j+1]);
            }
        }
    }
}

void Solution::roulette(double **probability, const vector<int>& sequence) const {
    vector<bool> indexes(hostIndexNum);
    for (int vertex:sequence) {
        vector<int> candidateIndex;
        vector<double> candidatePro;
        double candidateProSum = 0.0;
        int candidateNum = 0;
        for (int index = 0; index < hostIndexNum; ++index)
            if (!indexes[index]){
                candidateIndex.push_back(index);
                candidatePro.push_back(probability[vertex][index]);
                candidateProSum += probability[vertex][index];
                candidateNum ++;
            }
        double random = candidateProSum * (rand() % 1000) / 1000.0;
        double wheel = 0.0;
        int rIndex = candidateIndex[0];
        for (int i=0;i<candidateNum;++i){
            wheel += candidatePro[i];
            if (random<wheel){
                rIndex = candidateIndex[i];
                break;
            }
        }
        // printf("vertex %3d -> %3d (%.6f): random/sum %.6f/%.6f\n",vertex,rIndex,probability[vertex][rIndex],random,candidateProSum);
        assert(!indexes[rIndex]);
        indexes[rIndex] = true;
        injection[vertex] = rIndex;
        hostGraph[rIndex] = vertex;
    }
}

void Solution::initialingSolProbability(double **probability,const InstanceData& ins) {
    vector<int> sequence;
    vogel(probability,sequence);
    roulette(probability,sequence);
    getLength(ins);
}


//覆盖solution
void Solution::coverSol(const Solution& sol) {
    objectValue = sol.objectValue;
    maxEdge = sol.maxEdge;
    moveValue = sol.moveValue;
    memcpy(injection,sol.injection,sizeof(int)*hostIndexNum);
    memcpy(hostGraph,sol.hostGraph, sizeof(int) * hostIndexNum);
    for (int i = 0; i < edgeNum; i++)
        memcpy(edgeLength[i],sol.edgeLength[i],sizeof(int)*3);
}
//执行移动
void Solution::execute(MoveCandidate move, const InstanceData& ins) {
    int newIndex = injection[move.vertex2];
    int oldIndex = injection[move.vertex1];
    hostGraph[newIndex] = move.vertex1;
    hostGraph[oldIndex] = move.vertex2;
    injection[move.vertex1] = newIndex;
    injection[move.vertex2] = oldIndex;
    getLength(ins);
//    assert(maxEdge==move.maxEdge&&objectValue==move.objValue);//
    moveValue = move.moveValue;
}
//检查solution
void Solution::checkSol(const InstanceData& ins) const {
    bool judge = true;
    //检查对应关系
    for (int i=0; i < vertexNum; i++){
        int hostIndex = injection[i];
        if (hostIndex<0||hostIndex>=hostIndexNum){
            judge = false;
            cout << "The injection has overrun!!" << i << endl;
            break;
        }
        if (i != hostGraph[hostIndex]){
            cout << "The injection has mistake!!" << i << endl;
            judge = false;
        }
    }
    //检查最优值
    int obj = EMPTY;
    int tie = 1;
    for (int i = 0; i < hostIndexNum; i++)
        for (int j = 0; j < i; j++)
            if (ins.graph[i][j]==1){
                int dis = Func::getDistance(injection[i], injection[j], hostGraphSize);
                if (dis>obj){
                    obj = dis;
                    tie = 1;
                }else if (dis == obj)
                    tie++;
            }
    if (obj!=objectValue){
        cout<<"The objectValue has mistake!!"<<objectValue<<"/"<<obj<<endl;
        judge = false;
    }else if (tie != maxEdge){
        printf("The maxEdge has mistake!! %3d/%3d \n", maxEdge, tie);
    }
    if (!judge)
        cout<<endl<<"Check has mistakes!!"<<endl;
//    return judge;
}
//输出solution信息
void Solution::printHostGraph(bool judge,bool judge2) const {
    printf("The graphSize        :%d\n",hostGraphSize);
    printf("The objectValue      :%d\n",objectValue);
    printf("The maxEdgeNum       :%d\n",maxEdge);
    if (judge){
        cout<<"The host Graph:"<<endl;
        for (int i = 0; i < hostGraphSize; i++) {
            printf("%3d : ",i);
            for (int j = 0; j < hostGraphSize; j++){
                int index = i*hostGraphSize+j;
                if (hostGraph[index]>=vertexNum){
                    printf("(%3d) ",hostGraph[index]);
                }else{
                    printf(" %3d  ",hostGraph[index]);
                }
            }
            cout << endl;
        }
        if (judge2){
            cout<<endl<<"The edge Length:"<<endl;
            Func::printMatrix(edgeLength,edgeNum,3);
        }
    }
//    cout<<endl;
}

void Solution::printHostGraph(double **pro, bool judge, bool judge2) const {
    if (judge){
        cout<<"The host Graph:"<<endl;
        for (int i = 0; i < hostGraphSize; i++) {
            printf("%3d : ",i);
            for (int j = 0; j < hostGraphSize; j++){
                int index = i*hostGraphSize+j;
                if (hostGraph[index]>=vertexNum){
                    printf("%3d%.2f)",hostGraph[index],pro[hostGraph[index]][index]);
                }else{
                    printf("%3d%.2f ",hostGraph[index],pro[hostGraph[index]][index]);
                }
            }
            cout << endl;
        }
    }
}


void Solution::outputSol(const string& file) const {
    ofstream dataFile;
    dataFile.open(file);// 追加 ofstream::app
    if (!dataFile.is_open()) {
        fprintf(stderr, "Can not open file %s\n", file.c_str());
        exit(10);
    }
    dataFile << "vertexNumber "<<"edgeNumber "<<"objectValue "<<"injection(candidateGraphVertex : hostGraphIndex)"<<endl;
    dataFile << vertexNum<<" "<<edgeNum<<" "<<objectValue<<" "<<endl;
    for (int i = 0; i < vertexNum; ++i) {
        dataFile << i << " " << injection[i] <<endl;
    }
    dataFile.close();
}




TabuTable::TabuTable(const InstanceData& ins, int bTenure, int rTenure) {
    baseTenure = bTenure;
    randomTenure = rTenure;
    tabuTableSize = ins.hostIndexNum;
    tabuTable = new int * [tabuTableSize];
    for (int i = 0; i < tabuTableSize; i++) {
        tabuTable[i] = new int [tabuTableSize];
        memset(tabuTable[i],EMPTY, sizeof(int)*tabuTableSize);
    }
}
TabuTable::~TabuTable() {
    for (int i = 0; i < tabuTableSize; i++)
        delete[] tabuTable[i];
    delete[] tabuTable;
//    printf("The tabuTable has been delete!\n");
}
void TabuTable::tabu(MoveCandidate move,int iter,const Solution& sol) const{
    int hostIndex1 = sol.injection[move.vertex1];
    int hostIndex2 = sol.injection[move.vertex2];
    tabuTable[move.vertex1][hostIndex1] = iter + baseTenure + rand() % randomTenure;
    tabuTable[move.vertex2][hostIndex2] = iter + baseTenure + rand() % randomTenure;
}
void TabuTable::printTabuTable() const {
    printf("The tabu table size is %4d * %4d !\n",tabuTableSize,tabuTableSize);
    Func::printMatrix(tabuTable,tabuTableSize,tabuTableSize);
    cout<<endl;
}

bool TabuTable::judge(int vertex1, int vertex2, int iter, const Solution &sol) const{
    bool legal = true;
    int index1 = sol.injection[vertex1];
    int index2 = sol.injection[vertex2];
    //禁忌条件
    if (tabuTable[vertex1][index2] > iter || tabuTable[vertex2][index1] > iter)
        legal = false;
    return legal;
}


RoI::RoI(int size,int hostGraphSize,int obj,int slack) {
    regionSize = size;
    regions = new int * [regionSize];
    region = new int [regionSize];
    bandwidth = new int [regionSize];
    for (int i = 0; i < regionSize; ++i) {
        regions[i] = new int [regionSize];
    }
    updateRegions(hostGraphSize,obj,slack);
}

RoI::~RoI() {
    for (int i = 0; i < regionSize; ++i) {
        delete[] regions[i];
    }
    delete[] regions;
    delete[] region;
    delete[] bandwidth;
}

void RoI::updateRegions(int hostGraphSize,int obj,int slack) const {
    int dis;
    int size;
    int level = obj;
    if (obj<hostGraphSize - slack)
        level += slack;
    for (int i = 0; i < regionSize; ++i) {
        size = 0;
        for (int j = 0; j < regionSize; ++j) {
            regions[i][j] = EMPTY;
            dis = Func::getDistance(i,j,hostGraphSize);
            if (dis<=level)
                regions[i][size++] = j;
        }
    }
}
//对两个列表取交集，复杂度 n
void RoI::intersection(const int *list) const {
    int index = 0;
    while (list[index]!=EMPTY&&index<regionSize){
        region[list[index++]] ++;
    }
}
//对移动顶点index相连顶点的可接受区域regions取交集，复杂度 2m
void RoI::vertexRegion(int vertex, const InstanceData &ins, const Solution &sol, double &runtime) {
    clock_t begin,end;
    begin = clock();
    regionVertex.clear();
    //对顶点index的可接受区域进行初始化，复杂度 n
    for (int i = 0; i < regionSize; i++)
        region[i] = 0;
    //对顶点index的邻点的可接受区域取覆盖频率，复杂度 2m/n * n
    for (int i = 0; i < ins.degree[vertex]; i++)
        intersection(regions[sol.injection[ins.vertexVertex[vertex][i]]]);
    //将覆盖频率转化为交集，复杂度 n
    for (int i = 0;i<sol.hostIndexNum;i++)
        if (region[i]==ins.degree[vertex])
            regionVertex.push_back(sol.hostGraph[i]);
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}

void RoI::getBandwidth(const Solution &sol, double &runtime) {
    clock_t begin,end;
    begin = clock();
    for (int i = 0; i < regionSize; ++i)
        bandwidth[i] = 0;
    for (int edge=0;edge<sol.edgeNum;edge++)
        if (sol.objectValue == sol.edgeLength[edge][2]){
            bandwidth[sol.edgeLength[edge][0]] ++;
            bandwidth[sol.edgeLength[edge][1]] ++;
        }
    maxEdgeVertex.clear();
    for (int i = 0; i < regionSize; ++i) {
        if (bandwidth[i]>0)
            maxEdgeVertex.push_back(i);
    }
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}


Length::Length(int hostGraphSize,int n) {
    size = hostGraphSize*2 + 1;
    hostIndexNum = n;
    lengthNum = new int [size];
    cutLengthNum = new int [size];
    addLengthNum = new int [size];
    distance = new int * [hostIndexNum];
    for (int i = 0; i < hostIndexNum; ++i)
        distance[i] = new int [hostIndexNum];
    for (int i = 0; i < hostIndexNum; ++i) {
        for (int j = 0; j <= i; ++j)
            distance[i][j] = distance[j][i] = Func::getDistance(i,j,hostGraphSize);
    }
}

Length::~Length() {
    delete[] lengthNum;
    delete[] cutLengthNum;
    delete[] addLengthNum;
    for (int i = 0; i < hostIndexNum; ++i)
        delete[] distance[i];
    delete[] distance;
}

void Length::initialize(const Solution &sol, double &runtime) const {
    clock_t begin,end;
    begin = clock();
    for (int i = 0; i < size; ++i) {
        lengthNum[i] = 0;
    }
    for (int i = 0; i < sol.edgeNum; i++) {
        lengthNum[sol.edgeLength[i][2]]++;
    }
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}

void Length::cutLength(const Solution &sol, const InstanceData &ins, int vertex, double &runtime) const {
    clock_t begin,end;
    begin = clock();
    for (int i = 0; i < size; ++i)
        cutLengthNum[i] = lengthNum[i];
    for (int i=0;i<ins.degree[vertex]; i++)
        cutLengthNum[ sol.edgeLength[ ins.vertexEdge[vertex][i] ][2] ]--;
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}

void Length::addLength(const Solution &sol, const InstanceData &ins, MoveCandidate move, double &runtime) const {
    clock_t begin,end;
    begin = clock();
    for (int i = 0; i < size; ++i)
        addLengthNum[i] = cutLengthNum[i];
    for (int i=0;i<ins.degree[move.vertex2]; i++)
        addLengthNum[ sol.edgeLength[ ins.vertexEdge[move.vertex2][i] ][2] ]--;
    int dis;
    for (int i = 0; i < ins.degree[move.vertex1]; i++) {
        int otherVertex = ins.vertexVertex[move.vertex1][i];
//        dis = getDistance(sol.injection[otherVertex], sol.injection[move.vertex2], sol.hostGraphSize);
        dis = distance[sol.injection[otherVertex]][sol.injection[move.vertex2]];
        addLengthNum[dis]++;
    }
    for (int i = 0; i < ins.degree[move.vertex2]; i++) {
        int otherVertex = ins.vertexVertex[move.vertex2][i];
//        dis = getDistance(sol.injection[otherVertex], sol.injection[move.vertex1], sol.hostGraphSize);
        dis = distance[sol.injection[otherVertex]][sol.injection[move.vertex1]];
        addLengthNum[dis]++;
    }
    if (addLengthNum[0] != 0) {
//        dis = getDistance(sol.injection[move.vertex1], sol.injection[move.vertex2], sol.hostGraphSize);
        dis = distance[sol.injection[move.vertex1]][sol.injection[move.vertex2]];
        addLengthNum[dis] += addLengthNum[0];
        addLengthNum[0] = 0;
    }
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}

void Length::getMoveValue(int accuracy, MoveCandidate &move,double &runtime) const {
    clock_t begin,end;
    begin = clock();
    double moveValue = 0.0;
    int objValue = size - 1;
    int chu = 1;
    while(chu < 4 && objValue > 0){
        if (moveValue<1){
            if (addLengthNum[objValue]==0){
                objValue --;
            }else{
                moveValue += 1.0 * objValue;
                move.objValue = objValue;
                move.maxEdge = addLengthNum[objValue];
            }
        }else{
            double tieValue = 1.0 * addLengthNum[objValue];
            for (int i = 0; i < chu; ++i) {
                tieValue = tieValue / accuracy;
            }
            moveValue += tieValue;
            chu++;
            objValue--;
        }
    }
    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
    move.moveValue = moveValue;
}


int Func::getDistance(int loc1, int loc2, int loc3) {
    return abs(loc1/loc3 - loc2/loc3) + abs(loc1%loc3 - loc2%loc3);
}

int Func::getRouteNum(int ** dis,int v1,int v2,const InstanceData& ins){
    int line_dis = dis[v1][v2];
    int num = 0;
    for (int vertex=0;vertex<ins.vertexNum;++vertex){
        if (dis[v2][vertex] + dis[v1][vertex]==line_dis){
            num += 1;
        }
    }
    return num;
}

void Func::updateMoves(vector<MoveCandidate> &moves, MoveCandidate operation) {
    if (moves.empty() || operation.moveValue < moves[0].moveValue){
        moves.clear();
        moves.push_back(operation);
    }else if (operation.moveValue == moves[0].moveValue){
        moves.push_back(operation);
    }
}

void Func::printMatrix(int **matrix, int row, int col) {
    for (int i = 0; i < row; ++i) {
        printf("%4d :",i);
        for (int j = 0; j < col; ++j) {
            printf("%7d",matrix[i][j]);
        }
        cout<<endl;
    }
    printf("index:");
    for (int i = 0; i < col; ++i) {
        printf("%7d",i);
    }
    cout<<endl;
}

void Func::printMatrix(double **matrix, int row, int col) {
    for (int i = 0; i < row; ++i) {
        printf("%4d :",i);
        for (int j = 0; j < col; ++j) {
            printf("%.3f ",matrix[i][j]);
        }
        cout<<endl;
    }
    printf("index:");
    for (int i = 0; i < col; ++i) {
        printf("%5d ",i);
    }
    cout<<endl;
}

void Func::printMatrix(int **matrix, int **matrix1, int row, int col) {
    for (int i = 0; i < row; ++i) {
        printf("%4d :",i);
        for (int j = 0; j < col; ++j) {
            if (matrix[i][j]==matrix1[i][j])
                printf("%7d",matrix[i][j]);
            else
//                printf("%6d*",matrix[i][j]);
                printf("%4d/%4d",matrix[i][j],matrix1[i][j]);
        }
        cout<<endl;
    }
    printf("index:");
    for (int i = 0; i < col; ++i) {
        printf("%7d",i);
    }
    cout<<endl;
}

void Func::printMove(MoveCandidate move) {
    printf(" vertex %3d swap %3d    ,obj/maxEdge  : %3d / %3d (%.12f)!\n",
           move.vertex1, move.vertex2,move.objValue,move.maxEdge,move.moveValue);
}

void Func::printMoves(const vector<MoveCandidate> &moves) {
    for (MoveCandidate move:moves)
        printMove(move);
}

void Func::printList(int *list, int length) {
    for (int i = 0; i < length; ++i) {
        printf("%4d ",list[i]);
    }
    cout << endl ;
}

void Func::printList(const int *list1, const int *list2, int size) {
    bool judge = true;
    for (int i = 0; i < size; ++i) {
        if (list1[i]!=list2[i]){
            printf("%3d(%3d,%3d) ,",i,list1[i],list2[i]);
            judge = false;
        }
    }
    if (!judge)
        cout << "the two list are not same!" << endl;
}

Probability::Probability(int hostIndexNum,Input input) {
    size = hostIndexNum;
    r = input.reward;
    p = input.penalization;
    c = input.compensation;
    s = input.smoothC;
    t = input.smoothT;
    double initialP = 1.0/size;
    probability = new double * [size];
    for (int i = 0; i < size; ++i) {
        probability[i] = new double [size];
        for (int j = 0; j < size; ++j) {
            probability[i][j] = initialP;
        }
    }
}

Probability::Probability(InstanceData &ins,int hostIndexNum,Input input){
    size = hostIndexNum;
    r = input.reward;
    p = input.penalization;
    c = input.compensation;
    s = input.smoothC;
    t = input.smoothT;
    probability = new double * [size];
    for (int i = 0; i < size; ++i) {
        probability[i] = new double [size];
        for (int j = 0; j < size; ++j) {
            probability[i][j] = 0.0;
        }
    }
    // step1: 计算每个点与点之间的距离
    int ** distance = new int * [ins.hostIndexNum] ;
    for (int vertex=0;vertex<ins.hostIndexNum;++vertex){
        distance[vertex] = new int [ins.hostIndexNum];
        memset(distance[vertex],EMPTY,sizeof(int)*ins.hostIndexNum);
        distance[vertex][vertex] = 0;
    }
    for (int edge_index=0;edge_index<ins.edgeNum;++edge_index){
        distance[ins.edge[edge_index][0]][ins.edge[edge_index][1]] = 1;
        distance[ins.edge[edge_index][1]][ins.edge[edge_index][0]] = 1;
    }
    for (int vertex1=0;vertex1<ins.hostIndexNum;++vertex1){
        for (int vertex2=vertex1+1;vertex2<ins.hostIndexNum;++vertex2){
            if (distance[vertex1][vertex2]<1) continue;
            for (int vertex3=0;vertex3<ins.hostIndexNum;++vertex3){
                if (distance[vertex2][vertex3]>0 && distance[vertex1][vertex3]<0){
                    distance[vertex1][vertex3] = distance[vertex1][vertex2] + distance[vertex2][vertex3];
                    distance[vertex3][vertex1] = distance[vertex1][vertex2] + distance[vertex2][vertex3];
                }
            }
        }
    }
    // step2: 找到最长距离且度最小的点,距离相同选择度和最小的点,度和相同选择到达方案最少的点
    int * max_dis_vertex = new int [ins.vertexNum*4];
    int max_num = 0;
    int max_dis = EMPTY;
    int min_degree = 2*ins.vertexNum + 1;
    int min_routeNum = ins.vertexNum;
    for (int vertex1=0;vertex1<ins.hostIndexNum;++vertex1){
        for (int vertex2=vertex1;vertex2<ins.hostIndexNum;++vertex2){
            if (distance[vertex1][vertex2]>max_dis){
                max_dis = distance[vertex1][vertex2];
                max_num = 0;
                max_dis_vertex[max_num] = vertex1;
                max_dis_vertex[max_num+1] = vertex2;
                max_num += 2;
            }else if (distance[vertex1][vertex2]==max_dis){
                int sum_degree = ins.degree[vertex1] + ins.degree[vertex2];
                if (sum_degree==min_degree){
                    int rn = Func::getRouteNum(distance,vertex1,vertex2,ins);
                    if (rn == min_routeNum){
                        if (max_num > 4) continue;
                        max_dis_vertex[max_num] = vertex1;
                        max_dis_vertex[max_num+1] = vertex2;
                        max_num += 2;
                    }else if (rn < min_routeNum){
                        min_routeNum = rn;
                        max_num = 0;
                        max_dis_vertex[max_num] = vertex1;
                        max_dis_vertex[max_num+1] = vertex2;
                        max_num += 2;
                    }
                }else if (sum_degree<min_degree){
                    max_num = 0;
                    min_degree = sum_degree;
                    max_dis_vertex[max_num] = vertex1;
                    max_dis_vertex[max_num+1] = vertex2;
                    max_num += 2;
                }

            }
        }
    }
    // step3：为最前面的四个点分配host的四个角
    int * injection = new int [hostIndexNum];
    int * hostGraph = new int [hostIndexNum];
    memset(injection,EMPTY, sizeof(int) * hostIndexNum);
    memset(hostGraph,EMPTY, sizeof(int) * hostIndexNum);
    if (max_num>=4){
        if (max_dis_vertex[2]!=max_dis_vertex[0] && max_dis_vertex[2]!=max_dis_vertex[1] & max_dis_vertex[3]!=max_dis_vertex[0] && max_dis_vertex[3]!=max_dis_vertex[1]){
            max_num = 4;
        }else{
            max_num = 2;
        }
    }
    for (int line=0;line<max_num;++line){
        int be_index;
        if (line==0){
            be_index = 0;
        } else if (line == 1){
            be_index = ins.hostIndexNum - 1;
        } else if (line == 2){
            be_index = ins.hostGraphSize - 1;
        } else if (line == 3){
            be_index = ins.hostIndexNum - ins.hostGraphSize;
        }
        injection[max_dis_vertex[line]] = be_index;
        hostGraph[be_index] = max_dis_vertex[line];
    }
    delete[] max_dis_vertex;
    max_dis_vertex = nullptr;
    int max_vertex;
    int min_index;
    int stay_vertex_num = hostIndexNum - max_num;
    while (stay_vertex_num > 0){
        // g1:Find vertex
        float w1 = 0.8;
        float w2 = 0.2;
        float max_vertex_value = -1.0*MAXVALUE;
        max_vertex = EMPTY;
        for (int vertex=0;vertex<hostIndexNum;++vertex){
            if (injection[vertex]!=EMPTY){
                continue;
            }
            int A_num = 0;
            int B_num = 0;
            for (int vertex_v_i=0;vertex_v_i<ins.degree[vertex];++vertex_v_i){
                int ad_index = injection[ins.vertexVertex[vertex][vertex_v_i]];
                if (ad_index==EMPTY){
                    B_num ++;
                }else{
                    A_num ++;
                }
            }
            float vertex_value = w1*A_num - w2*B_num;
            if (vertex_value>max_vertex_value){
                max_vertex_value = vertex_value;
                max_vertex = vertex;
            }
        }
        // g2:Find index
        min_index = EMPTY;
        float min_index_value = MAXVALUE;
        for (int index=0;index<hostIndexNum;++index){
            if (hostGraph[index]!=EMPTY){
                continue;
            }
            float index_value = 0;
            int length = 1;
            int mi = 1;
            while (length < ins.hostGraphSize){
                int size = 0;
                for (int vertex_v_i=0;vertex_v_i<hostIndexNum;++vertex_v_i){
                    int ad_index = injection[vertex_v_i];
                    if (ad_index==EMPTY || distance[max_vertex][vertex_v_i]!=length || vertex_v_i==max_vertex){
                        continue;
                    }
                    size += 1;
                    float dis = Func::getDistance(index,ad_index,ins.hostGraphSize);
                    for (int line=0;line<mi-1;++line){
                        dis = 1.0*dis / (ins.hostGraphSize-1) / 2;
                    }
                    index_value += dis ;
                }
                if (size>0){
                    mi += 1;
                }
                length += 1;
                if (index_value > min_index_value){
                    break;
                }
            }
            if (index_value < min_index_value){
                min_index_value = index_value;
                min_index = index;
            }
        }
        injection[max_vertex] = min_index;
        hostGraph[min_index]  = max_vertex;
        stay_vertex_num --;
    }
    for (int vertex=0;vertex<hostIndexNum;++vertex){
        probability[vertex][injection[vertex]] = 1.0;
    }
    for (int vertex=0;vertex<hostIndexNum;++vertex){
        int s = 0;
        for (int v2=0;v2<hostIndexNum;++v2){
            s += probability[v2][vertex];
        }
    }
    // for (int v1=0;v1<hostIndexNum;++v1){
    //     for (int v2=0;v2<hostIndexNum;++v2){
    //         if (probability[v1][v2]==1){
    //             printf("%d %d : %f\n",v1,v2,probability[v1][v2]);
    //         }
    //     }
    // }
    for (int vertex=0;vertex<ins.vertexNum;++vertex){
        delete[] distance[vertex];
        distance[vertex] = nullptr;
    }
    delete[] distance;
    distance = nullptr;
    delete[] injection;
    delete[] hostGraph;
    injection = nullptr;
    hostGraph = nullptr;
}



Probability::~Probability() {
    for (int i = 0; i < size; ++i) {
        delete[] probability[i];
    }
    delete[] probability;
}

void Probability::update(const Solution &sol0, const Solution &sol1, double &runtime) const {
    //sol0 is new
    clock_t begin,end;
    begin = clock();
    double judgeSmooth;
    for (int vertex = 0; vertex < size; ++vertex) {
        int newIndex = sol0.injection[vertex];
        int oldIndex = sol1.injection[vertex];
        judgeSmooth = 0.0;
        for (int index = 0; index < size; ++index){
            if (newIndex==oldIndex){
                if (index==newIndex)
                    probability[vertex][index] = r + (1.0 - r)*probability[vertex][index];
                else
                    probability[vertex][index] = (1.0 - r)*probability[vertex][index];
            }else {
                if (index==newIndex)
                    probability[vertex][index] = c + (1.0 - c)*p/(size-1) + (1.0 - c) * (1.0 - p) * probability[vertex][index];
                else if (index==oldIndex)
                    probability[vertex][index] = (1.0 - c) * (1.0 - p) * probability[vertex][index];
                else
                    probability[vertex][index] = (1.0 - c)*p/(size-1) + (1.0 - c) * (1.0 - p) * probability[vertex][index];
            }
            if (probability[vertex][index]>t){
                probability[vertex][index] = probability[vertex][index] * s;
                judgeSmooth = probability[vertex][index];
            }
        }
        if (judgeSmooth>0.0)
            smooth(vertex,judgeSmooth);
    }

    end = clock();
    runtime += (double)(end - begin)/CLOCKS_PER_SEC;
}

void Probability::smooth(int vertex, double pro) const {
    double sumP = 1.0 - ( 1.0 - s )*pro;
    for (int index = 0; index < size; ++index)
        probability[vertex][index] = probability[vertex][index] / sumP ;
}

void Probability::initial() {
    double initialP = 1.0/size;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            probability[i][j] = initialP;
        }
    }
}


void Input::printInformation(string& out) const {
    //printf("The randomSeed      : %d\n", seed);
    //printf("The baseTenure      : %d\n", baseTenure);
    //printf("The randomTenure    : %d\n", randomTenure);
    //printf("The timeLimit       : %f\n", timeLimit);
    //printf("The improveCutOff   : %d\n", improveCut);
    //printf("The reward          : %f\n", reward);
    //printf("The penalization    : %f\n", penalization);
    //printf("The compensation    : %f\n", compensation);
    //printf("The smoothC         : %f\n", smoothC);
    //printf("The smoothT         : %f\n", smoothT);

    out += "The randomSeed      : " + to_string(seed) + "\n";
    out += "The baseTenure      : " + to_string(baseTenure) + "\n";
    out += "The randomTenure    : " + to_string(randomTenure) + "\n";
    out += "The timeLimit       : " + to_string(timeLimit) + "\n";
    out += "The improveCutOff   : " + to_string(improveCut) + "\n";
    out += "The reward          : " + to_string(reward) + "\n";
    out += "The penalization    : " + to_string(penalization) + "\n";
    out += "The compensation    : " + to_string(compensation) + "\n";
    out += "The smoothC         : " + to_string(smoothC) + "\n";
    out += "The smoothT         : " + to_string(smoothT) + "\n";
}

Input::Input(int seed, int baseT, int randomT, double timeLimit, int improveCut) {
    this->seed = seed;
    this->baseTenure = baseT;
    this->randomTenure = randomT;
    this->timeLimit = timeLimit;
    this->improveCut = improveCut;
    reward = 0.1;
    penalization = 0.3;
    compensation = 0.2;
    smoothC = 0.5;
    smoothT = 0.7;
    slack = 0;
}

void PLOutput::printInformation(string& out) const {
    //results
    //printf("The totalTime       : %.6f \n", totalTime);
    //printf("The tabuTime        : %.6f \n", tabuTime);
    //printf("The findBestTime    : %.6f \n", findBestTime);
    //printf("The iterations      : %d\n", iter);
    //printf("The allIterations   : %d\n", allIter);
    //printf("The stagnate        : %d\n", stagnate);
    //printf("The objectValues    : %d\n", obj);
    //printf("The maxEdge         : %d\n", maxEdge);
    //printf("The moveValue       : %.12f\n", moveValue);

    out += "The totalTime       : " + to_string(totalTime) + "\n";
    out += "The tabuTime        : " + to_string(tabuTime) + "\n";
    out += "The findBestTime    : " + to_string(findBestTime) + "\n";
    out += "The iterations      : " + to_string(iter) + "\n";
    out += "The allIterations   : " + to_string(allIter) + "\n";
    out += "The stagnate        : " + to_string(stagnate) + "\n";
    out += "The objectValues    : " + to_string(obj) + "\n";
    out += "The maxEdge         : " + to_string(maxEdge) + "\n";
    out += "The moveValue       : " + to_string(moveValue) + "\n";
}

void PLOutput::add(Output output) {
    allIter += output.iter;
    tabuTime += output.totalTime;
}

PLOutput::PLOutput() {
    obj = MAXVALUE;
    maxEdge = MAXVALUE;
}


void searchNeighbor_NRS(Solution &sol, InstanceData &ins, const TabuTable& tabuTable,
                    vector<MoveCandidate> &swaps, Output& output, RoI& roI, Length& length){
    clock_t begin,end;
    begin = clock();
    length.initialize(sol, output.evaluateTime);
    for (int vertex=0;vertex<ins.hostIndexNum;++vertex){
        length.cutLength(sol, ins, vertex, output.evaluateTime);
        roI.vertexRegion(vertex, ins, sol, output.regionTime);
        for (int otherVertex:roI.regionVertex){
            //定义移动操作
            MoveCandidate swap(vertex, otherVertex);
            //更新边长数量列表 n
            length.addLength(sol, ins, swap, output.evaluateTime);
            //计算目标值和平局评价值 n^0.5
            length.getMoveValue(sol.accuracy, swap, output.evaluateTime);
            if (tabuTable.judge(vertex, otherVertex, output.iter, sol) || swap.objValue<sol.objectValue)
                Func::updateMoves(swaps, swap);
            if (!swaps.empty() && int(swaps[0].moveValue * sol.edgeNum) < int(sol.moveValue * sol.edgeNum) )
                break;
        }
        if (!swaps.empty() && int(swaps[0].moveValue * sol.edgeNum) < int(sol.moveValue * sol.edgeNum) )
            break;
    }
    end = clock();
    output.searchTime += (double)(end - begin) / CLOCKS_PER_SEC;
}

void searchNeighbor(Solution &sol, InstanceData &ins, const TabuTable& tabuTable,
                    vector<MoveCandidate> &swaps, Output& output, RoI& roI, Length& length){
    clock_t begin,end;
    begin = clock();
    length.initialize(sol, output.evaluateTime);
    roI.getBandwidth(sol, output.regionTime);
    for (int vertex:roI.maxEdgeVertex){
        length.cutLength(sol, ins, vertex, output.evaluateTime);
        roI.vertexRegion(vertex, ins, sol, output.regionTime);
        for (int otherVertex:roI.regionVertex){
            //定义移动操作
            MoveCandidate swap(vertex, otherVertex);
            //更新边长数量列表 n
            length.addLength(sol, ins, swap, output.evaluateTime);
            //计算目标值和平局评价值 n^0.5
            length.getMoveValue(sol.accuracy, swap, output.evaluateTime);
            if (tabuTable.judge(vertex, otherVertex, output.iter, sol) || swap.objValue<sol.objectValue)
                Func::updateMoves(swaps, swap);
            if (!swaps.empty() && int(swaps[0].moveValue * sol.edgeNum) < int(sol.moveValue * sol.edgeNum) )
                break;
        }
        if (!swaps.empty() && int(swaps[0].moveValue * sol.edgeNum) < int(sol.moveValue * sol.edgeNum) )
            break;
    }
    end = clock();
    output.searchTime += (double)(end - begin) / CLOCKS_PER_SEC;
}

void tabuSearch(Solution &sol, InstanceData &ins, Input input,Output &output, Probability &pro,double baseTime){
    clock_t begin,end;
    begin = clock();
    Solution tempSol(ins);
    tempSol.coverSol(sol);
    TabuTable tabuTable(ins, input.baseTenure, input.randomTenure);
    RoI roI(sol.hostIndexNum,sol.hostGraphSize,sol.objectValue,input.slack);
    vector<MoveCandidate> moves;
    Length length(sol.hostGraphSize,sol.hostIndexNum);
    while (output.stagnate < input.improveCut &&  baseTime+output.totalTime < input.timeLimit){
        output.iter++;
        output.stagnate++;
        //搜索邻域
        searchNeighbor(tempSol, ins, tabuTable, moves, output, roI, length);
        //判断邻域是否为空
        if (!moves.empty()){
            output.selectNum += int(moves.size());
            int randomMove = int(rand()%moves.size());
            tabuTable.tabu(moves[randomMove], output.iter, tempSol);
            tempSol.execute(moves[randomMove], ins);
            moves.clear();
        }else
            output.emptyNum++;
        //更新解
        if (tempSol.objectValue<sol.objectValue){
           pro.update(tempSol, sol, output.probabilityTime);
            output.stagnate = 0;
            sol.coverSol(tempSol);
            end = clock();
            output.findBestTime = (double)(end - begin) / CLOCKS_PER_SEC;
            roI.updateRegions(sol.hostGraphSize,sol.objectValue,input.slack);
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
           pro.update(tempSol, sol, output.probabilityTime);
            sol.coverSol(tempSol);
        }
        end = clock();
        output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.tie = sol.maxEdge;
    output.moveValue = sol.moveValue;
    end = clock();
    output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
}

void tabuSearch_NRS(Solution &sol, InstanceData &ins, Input input,Output &output, Probability &pro,double baseTime){
    clock_t begin,end;
    begin = clock();
    Solution tempSol(ins);
    tempSol.coverSol(sol);
    TabuTable tabuTable(ins, input.baseTenure, input.randomTenure);
    RoI roI(sol.hostIndexNum,sol.hostGraphSize,sol.objectValue,input.slack);
    vector<MoveCandidate> moves;
    Length length(sol.hostGraphSize,sol.hostIndexNum);
    while (output.stagnate < input.improveCut &&  baseTime+output.totalTime < input.timeLimit){
        output.iter++;
        output.stagnate++;
        //搜索邻域
        searchNeighbor_NRS(tempSol, ins, tabuTable, moves, output, roI, length);
        //判断邻域是否为空
        if (!moves.empty()){
            output.selectNum += int(moves.size());
            int randomMove = int(rand()%moves.size());
            tabuTable.tabu(moves[randomMove], output.iter, tempSol);
            tempSol.execute(moves[randomMove], ins);
            moves.clear();
        }else
            output.emptyNum++;
        //更新解
        if (tempSol.objectValue<sol.objectValue){
           pro.update(tempSol, sol, output.probabilityTime);
            output.stagnate = 0;
            sol.coverSol(tempSol);
            end = clock();
            output.findBestTime = (double)(end - begin) / CLOCKS_PER_SEC;
            roI.updateRegions(sol.hostGraphSize,sol.objectValue,input.slack);
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
           pro.update(tempSol, sol, output.probabilityTime);
            sol.coverSol(tempSol);
        }
        end = clock();
        output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.tie = sol.maxEdge;
    output.moveValue = sol.moveValue;
    end = clock();
    output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
}

void tabuSearch_nolearning(Solution &sol, InstanceData &ins, Input input,Output &output,double baseTime){
    clock_t begin,end;
    begin = clock();
    Solution tempSol(ins);
    tempSol.coverSol(sol);
    TabuTable tabuTable(ins, input.baseTenure, input.randomTenure);
    RoI roI(sol.hostIndexNum,sol.hostGraphSize,sol.objectValue,input.slack);
    vector<MoveCandidate> moves;
    Length length(sol.hostGraphSize,sol.hostIndexNum);
    while (output.stagnate < input.improveCut &&  baseTime+output.totalTime < input.timeLimit){
        output.iter++;
        output.stagnate++;
        //搜索邻域
        searchNeighbor(tempSol, ins, tabuTable, moves, output, roI, length);
        //判断邻域是否为空
        if (!moves.empty()){
            output.selectNum += int(moves.size());
            int randomMove = int(rand()%moves.size());
            tabuTable.tabu(moves[randomMove], output.iter, tempSol);
            tempSol.execute(moves[randomMove], ins);
            moves.clear();
        }else
            output.emptyNum++;
        //更新解
        if (tempSol.objectValue<sol.objectValue){
            output.stagnate = 0;
            sol.coverSol(tempSol);
            end = clock();
            output.findBestTime = (double)(end - begin) / CLOCKS_PER_SEC;
            roI.updateRegions(sol.hostGraphSize,sol.objectValue,input.slack);
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
            sol.coverSol(tempSol);
        }
        end = clock();
        output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.tie = sol.maxEdge;
    output.moveValue = sol.moveValue;
    end = clock();
    output.totalTime = (double)(end - begin) / CLOCKS_PER_SEC;
}

void PLTS(Solution &sol, InstanceData &ins, Input input,Output &inOutput,PLOutput &output){
    Probability pro(ins,ins.hostIndexNum,input);
    Solution tempSol(ins);
    while(output.totalTime<input.timeLimit && output.stagnate<input.stop && output.iter<input.maxIter){
        output.iter ++;
        output.stagnate ++;
        tempSol.initialingSolProbability(pro.probability,ins);
        if (IFPRINT) printf("Iter.%4d , initialSol : %4d ( %4d) ->  ",output.iter,tempSol.objectValue,tempSol.maxEdge);
        Output tempOutput;
        tabuSearch(tempSol,ins,input,tempOutput,pro,output.totalTime);
        output.add(tempOutput);
        if (tempSol.objectValue<sol.objectValue){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
            output.findBestTime = output.totalTime + tempOutput.findBestTime;
            output.stagnate = 0;
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
        }
        if (IFPRINT) printf(" %4d ( %4d) \n",sol.objectValue,sol.maxEdge);
        clock_t end = clock();
        output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.maxEdge = sol.maxEdge;
    output.moveValue = sol.moveValue;
    clock_t end = clock();
    output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
}

void ITS(Solution &sol, InstanceData &ins, Input input,Output &inOutput,PLOutput &output){
    Solution tempSol(ins);
    while(output.totalTime<input.timeLimit && output.stagnate<input.stop && output.iter<input.maxIter){
        output.iter ++;
        output.stagnate ++;
        tempSol.initialingSolRandom(ins);
        if (IFPRINT) printf("Iter.%4d , initialSol : %4d ( %4d) ->  ",output.iter,tempSol.objectValue,tempSol.maxEdge);
        Output tempOutput;
        tabuSearch_nolearning(tempSol,ins,input,tempOutput,output.totalTime);
        output.add(tempOutput);
        if (tempSol.objectValue<sol.objectValue){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
            output.findBestTime = output.totalTime + tempOutput.findBestTime;
            output.stagnate = 0;
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
        }
        if (IFPRINT) printf(" %4d ( %4d) \n",sol.objectValue,sol.maxEdge);
        clock_t end = clock();
        output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.maxEdge = sol.maxEdge;
    output.moveValue = sol.moveValue;
    clock_t end = clock();
    output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
}

void GITS(Solution &sol, InstanceData &ins, Input input,Output &inOutput,PLOutput &output){
    Solution tempSol(ins);
    while(output.totalTime<input.timeLimit && output.stagnate<input.stop && output.iter<input.maxIter){
        output.iter ++;
        output.stagnate ++;
        tempSol.initialingSolGreedyly(ins);
        if (IFPRINT) printf("Iter.%4d , initialSol : %4d ( %4d) ->  ",output.iter,tempSol.objectValue,tempSol.maxEdge);
        Output tempOutput;
        tabuSearch_nolearning(tempSol,ins,input,tempOutput,output.totalTime);
        output.add(tempOutput);
        if (tempSol.objectValue<sol.objectValue){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
            output.findBestTime = output.totalTime + tempOutput.findBestTime;
            output.stagnate = 0;
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
        }
        if (IFPRINT) printf(" %4d ( %4d) \n",sol.objectValue,sol.maxEdge);
        clock_t end = clock();
        output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.maxEdge = sol.maxEdge;
    output.moveValue = sol.moveValue;
    clock_t end = clock();
    output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
}


void NRS(Solution &sol, InstanceData &ins, Input input,Output &inOutput,PLOutput &output){
    Probability pro(ins.hostIndexNum,input);
    Solution tempSol(ins);
    while(output.totalTime<input.timeLimit && output.stagnate<input.stop && output.iter<input.maxIter){
        output.iter ++;
        output.stagnate ++;
        tempSol.initialingSolProbability(pro.probability,ins);
        if (IFPRINT) printf("Iter.%4d , initialSol : %4d ( %4d) ->  ",output.iter,tempSol.objectValue,tempSol.maxEdge);
        Output tempOutput;
        tabuSearch_NRS(tempSol,ins,input,tempOutput,pro,output.totalTime);
        output.add(tempOutput);
        if (tempSol.objectValue<sol.objectValue){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
            output.findBestTime = output.totalTime + tempOutput.findBestTime;
            output.stagnate = 0;
        }else if (tempSol.objectValue == sol.objectValue && tempSol.maxEdge<sol.maxEdge){
            sol.coverSol(tempSol);
            inOutput.copy(tempOutput);
        }
        if (IFPRINT) printf(" %4d ( %4d) \n",sol.objectValue,sol.maxEdge);
        clock_t end = clock();
        output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
    }
    output.obj = sol.objectValue;
    output.maxEdge = sol.maxEdge;
    output.moveValue = sol.moveValue;
    clock_t end = clock();
    output.totalTime = (double)(end - input.begin) / CLOCKS_PER_SEC;
}
