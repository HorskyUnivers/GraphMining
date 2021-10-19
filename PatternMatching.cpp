#include "PatternMatching.h"
PatternMatching::~PatternMatching() {
    if (R_visited) {
        free(R_visited);
        R_visited = nullptr;
    }
    if (degree_R) {
        free(degree_R);
        degree_R = nullptr;
    }

    /*if (PMR) {
        free(PMR);
        PMR = nullptr;
    }*/

    if (R_adj) {
        free(R_adj);
        R_adj = nullptr;
    }
    if (R_reverse_adj) {
        free(R_reverse_adj);
        R_reverse_adj = nullptr;
    }
    if (R_adjIndex) {
        free(R_adjIndex);
        R_adjIndex = nullptr;
    }
    if (R_reverseAdjIndex) {
        free(R_reverseAdjIndex);
        R_reverseAdjIndex = nullptr;
    }
}

void PatternMatching::init() {
    R_visited = (unsigned*)malloc(vertexNum_R * sizeof(unsigned));
    memset(R_visited, 0, vertexNum_R * sizeof(unsigned));
}

bool PatternMatching::build_degree_R(std::string inputfile, unsigned vertexNum)//创建degree_R，测试成功
{
    maxID = 0;     //图中最大id
    edgeNum_R = 0; //图的边数
    vertexNum_R = vertexNum;
    degree_R = (Degree*)calloc(vertexNum_R, sizeof(Degree));
    FILE* inf = fopen(inputfile.c_str(), "r");
    if (inf == NULL)
    {
        std::cerr << "Could not load :" << inputfile << " error: " << strerror(errno)
            << std::endl;
    }
    assert(inf != NULL);
    std::cout << "Reading in adjacency list format!" << std::endl;
    int maxlen = 100000000;
    char* s = (char*)malloc(maxlen); //暂时存放一次读入的数据
    size_t bytesread = 0;             //读到的总字节数
    char delims[] = " \t";            //字符串的分隔符（tab键）
    size_t linenum = 0;               //行数
    size_t lastlog = 0;               //最后一个字节
    unsigned from = 0, to = 0;        //起点id  //终点id
    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                  //将换行符变成'\0'
        char* t = strtok(s, delims); //返回的是tab之前的那一段字符串的首地址
        from = atoi(t);              //提取起点id
        unsigned num = 0;            //出度个数
        while ((t = strtok(NULL, delims)) != NULL)
        {
            to = atoi(t);
            if (from != to)
            {
                maxID = max(to, maxID); //找出整张图中的最大id
                degree_R[from].outdeg++;       //出度
                degree_R[to].indeg++;          //入度
                num++;
            }
        }
        edgeNum_R += num; //统计所有from的总的出度个数 == 总的边数
    }
    std::cout << "finish first read R_adj... maxID:" << maxID << std::endl;
    free(s);
    fclose(inf);
    return true;
}
void PatternMatching::FIXLINE(char* s)
{
    int len = (int)strlen(s) - 1;
    if (s[len] == '\n')
        s[len] = 0;
}
bool PatternMatching::build_R_adj(std::string inputfile)//创建邻接表R_adj和逆邻接表R_reverse_adj，测试成功
{
    R_adj = new unsigned[edgeNum_R];
    R_reverse_adj = new unsigned[edgeNum_R];
    R_adjIndex = new unsigned[vertexNum_R];
    R_reverseAdjIndex = new unsigned[vertexNum_R];
    unsigned* R_reverseAdjIndex_tail = new unsigned[vertexNum_R];
    memset(R_reverseAdjIndex_tail, 0, vertexNum_R);
    R_adjIndex[0] = 0;
    R_reverseAdjIndex[0] = 0;
    R_reverseAdjIndex_tail[0] = 0;
    unsigned tail = 0;
    FILE* inf = fopen(inputfile.c_str(), "r");
    if (inf == NULL)
    {
        std::cerr << "Could not load :" << inputfile << " error: " << strerror(errno)
            << std::endl;
    }
    assert(inf != NULL);
    std::cout << "Reading in adjacency list format!" << std::endl;
    int maxlen = 100000000;
    char* s = (char*)malloc(maxlen); //暂时存放一次读入的数据
    char delims[] = " \t";            //字符串的分隔符（tab键）
    size_t linenum = 0;               //行数
    unsigned from = 0, to = 0;        //起点id  //终点id
    for (unsigned i = 1; i < vertexNum_R; i++) {
        R_adjIndex[i] = R_adjIndex[i - 1] + degree_R[i - 1].outdeg;
        R_reverseAdjIndex[i] = R_reverseAdjIndex[i - 1] + degree_R[i - 1].indeg;
        R_reverseAdjIndex_tail[i] = R_reverseAdjIndex[i];
    }

    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                                      //将换行符变成'\0'
        char* t = strtok(s, delims);                     //返回的是tab之前的那一段字符串的首地址
        from = atoi(t);                                  //提取起点id
        while ((t = strtok(NULL, delims)) != NULL)
        {
            to = atoi(t);
            if (from != to)
            {
                R_adj[tail++] = to;
                R_reverse_adj[R_reverseAdjIndex_tail[to]++] = from;
            }
        }
    }
    std::cout << "finish second read R_adj..." << std::endl;
    free(s);
    free(R_reverseAdjIndex_tail);
    fclose(inf);
    return true;
}
bool PatternMatching::build_P_adj(std::string inputfile, unsigned vertexNum) {
    vertexNum_P = vertexNum;
    degree_P = (Degree*)calloc(vertexNum_P, sizeof(Degree));
    std::vector<std::vector<P_ID>> tmp(vertexNum_P, std::vector<P_ID>(vertexNum_P));//临时存储，方便初始化
    FILE* inf = fopen(inputfile.c_str(), "r");
    if (inf == NULL)
    {
        std::cerr << "Could not load :" << inputfile << " error: " << strerror(errno)
            << std::endl;
    }
    assert(inf != NULL);
    std::cout << "Reading in adjacency list format!" << std::endl;
    int maxlen = 100000000;
    char* s = (char*)malloc(maxlen); //暂时存放一次读入的数据
    char delims[] = " \t";            //字符串的分隔符（tab键）
    size_t linenum = 0;               //行数
    unsigned from = 0, to = 0;        //起点id  //终点id
    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                                      //将换行符变成'\0'
        char* t = strtok(s, delims);                     //返回的是tab之前的那一段字符串的首地址
        from = atoi(t);                                  //提取起点id
        while ((t = strtok(NULL, delims)) != NULL)
        {
            to = atoi(t);
            if (from != to)
            {
                degree_P[from].outdeg++;       //出度
                degree_P[to].indeg++;          //入度
                tmp[from][to] = 1;
            }
        }
    }
    P_adj = tmp;
    std::cout << "finish first read P_adj..." << std::endl;
    free(s);
    fclose(inf);
    return true;
}

bool PatternMatching::matchPR_expand() {
    total_num = 0;
    sel.resize(vertexNum_P);//初始化选择度
    minMatchID = 0; //具有最小匹配数的模式图节点ID
    minMatchNum = UINT_MAX;  //最小匹配数
    unsigned current_size = 0;
    unsigned total_size = 0;//存储M(Vp)所需的空间大小

    for (unsigned i = 0; i < vertexNum_P; i++)//计算与每个模式图点匹配的实际图中节点数量
    {
        for (R_ID j = 0; j < vertexNum_R; j++) //j要修改成R中区域内的节点id区间
        {
            if (degree_R[j].indeg >= degree_P[i].indeg && degree_R[j].outdeg >= degree_P[i].outdeg)
            {
                current_size++;
            }
        }
        sel[i] = current_size;
        PMR_index.push_back(current_size);
        total_size += current_size;

        if (current_size < minMatchNum) {
            minMatchNum = current_size;
            minMatchID = i;
        }
        current_size = 0;
    }
    total_num = total_size;
    std::vector<int> tmp;
    tmp.push_back(0);
    for (unsigned i = 0; i < vertexNum_P - 1; i++) {
        tmp.push_back(PMR_index[i] + tmp[tmp.size() - 1]);
    }
    PMR_index = tmp;
    PMR.resize(vertexNum_P);
    int k = 0;

    //只更新最小匹配模式图结点的PMR集合，其他保持空的状态
    /*for (R_ID j = 0; j < vertexNum_R; j++) 
    {
        if (degree_R[j].indeg >= degree_P[minMatchID].indeg && degree_R[j].outdeg >= degree_P[minMatchID].outdeg)
        {
            PMR[minMatchID].emplace_back(j);
        }
    }*/

    std::cout << "finish compute minMatchID ..." << std::endl;
    return true;
}

//PatternMatching::matchPR()

//bool PatternMatching::matchPR() //创建初始匹配集合P_M_R,测试成功
//{
//    total_num = 0;
//    sel.resize(vertexNum_P);//初始化选择度
//    minMatchID = 0; //具有最小匹配数的模式图节点ID
//    minMatchNum = UINT_MAX;  //最小匹配数
//    unsigned current_size = 0;
//    unsigned total_size = 0;//存储M(Vp)所需的空间大小
//
//    for (unsigned i = 0; i < vertexNum_P; i++)//计算与每个模式图点匹配的实际图中节点数量
//    {
//        for (R_ID j = 0; j < vertexNum_R; j++) //j要修改成R中区域内的节点id区间
//        {
//            if (degree_R[j].indeg >= degree_P[i].indeg && degree_R[j].outdeg >= degree_P[i].outdeg)
//            {
//                current_size++;
//            }
//        }
//        sel[i] = current_size;
//        PMR_index.push_back(current_size);
//        total_size += current_size;
//
//        if (current_size < minMatchNum) {
//            minMatchNum = current_size;
//            minMatchID = i;
//        }
//        /*
//        minMatchNum = min(current_size, minMatchNum);
//        if (minMatchNum == current_size)
//            minMatchID = i;
//        */
//        current_size = 0;
//    }
//    total_num = total_size;
//    std::vector<int> tmp;
//    tmp.push_back(0);
//    for (unsigned i = 0; i < vertexNum_P - 1; i++) {
//        tmp.push_back(PMR_index[i] + tmp[tmp.size() - 1]);
//    }
//    PMR_index = tmp;
//    PMR = (unsigned*)calloc(total_size, sizeof(unsigned));
//    //FIXME,sky,PMR's value is wrong.
//    /*
//    int k = 0;
//    for (unsigned i = 0; i < vertexNum_P; i++) {
//        for (unsigned j = 0; j < vertexNum_R; j++) {
//            PMR[k++] = j;
//        }
//    }
//    */
//    int k = 0;
//    for (unsigned i = 0; i < vertexNum_P; i++)//计算与每个模式图点匹配的实际图中节点数量
//    {
//        for (R_ID j = 0; j < vertexNum_R; j++) //j要修改成R中区域内的节点id区间
//        {
//            if (degree_R[j].indeg >= degree_P[i].indeg && degree_R[j].outdeg >= degree_P[i].outdeg)
//            {
//                PMR[k++] = j;
//            }
//        }
//    }
//
//    std::cout << "finish compute minMatchID ..." << std::endl;
//    return true;
//}

bool PatternMatching::isNextEPatternEmpty()
{
    for (unsigned i = 0; i < P_adj.size(); i++)
    {
        for (unsigned j = 0; j < P_adj[i].size(); j++)
        {
            if (P_adj[i][j] == 1)
                return false;
        }
    }
    return true;
}

bool PatternMatching::isfinish(P_ID id) {
    for (P_ID i = 0; i < vertexNum_P; i++) {
        if (P_adj[id][i] == 1)
            return false;
    }
    for (P_ID i = 0; i < vertexNum_P; i++) {
        if (P_adj[i][id] == 1)
            return false;
    }
    return true;
}

//P_ID PatternMatching::getMaxSel() {
//    int maxSelID = 0;
//    P_ID curid = 0;
//    bool flag = true;
//    for (; curid < vertexNum_P; curid++) {//找到id最小的模式边没有被全部访问的点
//        P_ID j = 0;
//        for (; j < vertexNum_P; j++) {
//            if (P_adj[curid][j] == 1) {
//                flag = false;
//                break;
//            }
//        }
//        if (!flag)
//            break;
//    }
//    maxSelID = curid;
//    for (P_ID i = 0; i < sel.size(); i++) {
//        if (sel[i] < sel[maxSelID] && !isfinish(i)) {  //匹配数越小的模式图节点选择度越大
//            maxSelID = i;
//        }
//    }
//    return maxSelID;
//}


//这里是得到当前选择度表下匹配数最少的模式图点编号并返回
P_ID PatternMatching::getMaxSel_cur() {
    int maxSelID = 0;
    int maxSel = MAXINT;
    P_ID curid = minMatchID;

    for (int i = 0; i < vertexNum_P; ++i) {
        if (P_adj[curid][i] == 1 || P_adj[i][curid]) {
            if (sel[i] < maxSel) {
                maxSel = sel[i];
                maxSelID = i;
            }
        }
    }

    //bool flag = true;
    //for (; curid < vertexNum_P; curid++) {//找到id最小的模式边没有被全部访问的点
    //    P_ID j = 0;
    //    for (; j < vertexNum_P; j++) {
    //        if (P_adj[curid][j] == 1) {
    //            flag = false;
    //            break;
    //        }
    //    }
    //    if (!flag)
    //        break;
    //}
    //maxSelID = curid;
    //for (P_ID i = 0; i < sel.size(); i++) {
    //    if (sel[i] < sel[maxSelID] && !isfinish(i)) {  //匹配数越小的模式图节点选择度越大
    //        maxSelID = i;
    //    }
    //}

    return maxSelID;
}

void PatternMatching::searchAllPR() {
    std::cout << "Start graph mining..." << std::endl;

    std::vector<int> minMatchID_PMR;
    for (R_ID j = 0; j < vertexNum_R; j++)
    {
        if (degree_R[j].indeg >= degree_P[minMatchID].indeg && degree_R[j].outdeg >= degree_P[minMatchID].outdeg)
        {
            minMatchID_PMR.emplace_back(j);
        }
    }

    for (int i = 0; i < minMatchID_PMR.size(); ++i) {
        PMR[minMatchID].emplace_back(minMatchID_PMR[i]);
        R_visited[minMatchID_PMR[i]] = 1;  //标记vr,s为访问中心点

        //FIXME,这里或许可以优化，sel在此前的初始化计算是没有必要的
        //Initialize模式图每顶点的选择度Sel;除开vp,s之外，其它模式图顶点选择度为无穷大
        /*for (int i = 0; i < vertexNum_P; ++i) {
            if (i != minMatchID) {
                sel[i] = MAXINT;
            }
        }*/

        //由于需要复用PMR和sel，这里更改为需要参数复制传递的searchPG();
        //searchPR();        
        searchPG(PMR, sel, minMatchID_PMR[i]);
    }
    /*for (unsigned i = PMR_index[minMatchID],j=0; j<sel[minMatchID]; j++,i++) {
        searchPR();
        R_visited[i] = 1;
    }*/
}

void PatternMatching::searchPG(std::vector<std::vector<unsigned>>PMR_copy, std::vector<int>sel_copy, unsigned CenterPoint) {
    //Initialize模式图每顶点的选择度Sel;除开vp,s之外，其它模式图顶点选择度为无穷大
    for (int i = 0; i < vertexNum_P; ++i) {
        if (i != minMatchID) {
            sel_copy[i] = MAXINT;
        }
    }
    while (!isNextEPatternEmpty()) {
        P_ID maxSelId = getMaxSel_cur();
        //P_ID neighborID = 0; //由这两个点构成最小匹配的边模式
        unsigned minMatchNum = INT_MAX;//neighborID的最小匹配数
        bool isReverse = true;
        /*for (P_ID i = 0; i < vertexNum_P; i++) {
            if (P_adj[i][maxSelId] == 1 && sel[i] < minMatchNum) {
                minMatchNum = sel[i];
                neighborID = i;
            }
        }
        for (P_ID i = 0; i < vertexNum_P; i++) {
            if (P_adj[maxSelId][i] == 1 && sel[i] < minMatchNum) {
                minMatchNum = sel[i];
                neighborID = i;
                isReverse = false;
            }
        }*/
        if (P_adj[minMatchID][maxSelId] == 1) {
            isReverse = false;
        }
        else if (P_adj[maxSelId][minMatchID] == 1) {
            isReverse = true;
        }
        else {
            std::cout << "Error: unexpected edge in searchpg." << std::endl;
        }

        if (!isReverse) {
            extendEdgePattern(minMatchID,maxSelId);
        }
        else {
            reverse_extendEdgePattern(maxSelId,minMatchID);
        }
    }
}

//原版的searchPR函数，目前已不使用，现在请使用searchPG函数
void PatternMatching::searchPR() {
    while (!isNextEPatternEmpty()) {
        P_ID maxSelId = getMaxSel();
        P_ID neighborID = 0; //由这两个点构成最小匹配的边模式
        unsigned minMatchNum = INT_MAX;//neighborID的最小匹配数
        bool isReverse = true;
        for (P_ID i = 0; i < vertexNum_P; i++) {
            if (P_adj[i][maxSelId] == 1 && sel[i] < minMatchNum) {
                minMatchNum = sel[i];
                neighborID = i;
            }
        }
        for (P_ID i = 0; i < vertexNum_P; i++) {
            if (P_adj[maxSelId][i] == 1 && sel[i] < minMatchNum) {
                minMatchNum = sel[i];
                neighborID = i;
                isReverse = false;
            }
        }
        if (!isReverse) {
            extendEdgePattern(maxSelId, neighborID);
        }
        else {
            reverse_extendEdgePattern(neighborID, maxSelId);
        }
    }
}

//存储首值索引版本，得到M(v_p)的索引起止范围,筛除尾部NULL值，成功则返回true以及start与end值，未找到返回false
bool PatternMatching::get_Mvp_Index(P_ID v_p, R_ID& start, R_ID& end) {
    if (v_p == PMR_index.size() - 1) {
        start = PMR_index[v_p];
        end = total_num - 1;
    }
    else {
        start = PMR_index[v_p];
        end = PMR_index[v_p + 1] - 1;
    }
    /*while (end >= start) {
        if (PMR[end] != NULL) return true;
        else end--;
    }*/
    if (start <= end) {
        return true;
    }
    return false;
}

//存储尾值索引版本，得到M(v_p)的索引起止范围,找到返回true以及start与end值，未找到返回false
/*
bool PatternMatching::get_Mvp_Index(P_ID v_p, int& start, int& end) {
    if (v_p == 0) {
        start = 0;
        end = PMR_index[v_pt] - 1;
    }
    else {
        start = PMR_index[v_pt - 1];
        end = PMR_index[v_pt] - 1;
    }
    while (end >= start) {
        if (PMR[end] != NULL) return true;
        else end--;
    }
    return false;
}
*/

//得到图G上以vrs为源，与边模式vps, vpt相匹配的端点集合，即以vrs为源点与vpt相匹配的顶点集
std::vector<R_ID> PatternMatching::getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs) {
    int length = degree_R[v_rs].outdeg;
    std::vector<R_ID> ans;

    //从PMR_index中获取M(v_pt)的索引起止范围
    R_ID M_v_pt_start, M_v_pt_end;
    if (!get_Mvp_Index(v_pt, M_v_pt_start, M_v_pt_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
    }

    //获取R_adj中v_rs点一行的所有出度结点，遍历并判断是否属于M(v_pt)
    for (int j = 0; j < length; j++) {
        R_ID curr = R_adj[j + R_adjIndex[v_rs]];
        //如果报错，可以改成for循环遍历判断是否相等
        if (std::find(PMR + M_v_pt_start, PMR + M_v_pt_end + 1, curr) != PMR + M_v_pt_end + 1) {
            ans.emplace_back(curr);
        }
    }
    return ans;
}

//得到图G上以vrt为终点，与边模式vps, vpt相匹配的端点集合，即以vrt为终点点与vps相匹配的顶点集
std::vector<R_ID> PatternMatching::reverse_getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rt) {
    int length = degree_R[v_rt].outdeg;
    std::vector<R_ID> ans;

    //从PMR_index中获取M(v_ps)的索引起止范围
    R_ID M_v_ps_start, M_v_ps_end;
    if (!get_Mvp_Index(v_ps, M_v_ps_start, M_v_ps_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_ps << std::endl;
    }

    //获取R_adj中v_rt点一行的所有出度结点，遍历并判断是否属于M(v_ps)
    for (int j = 0; j < length; j++) {
        R_ID curr = R_adj[j + R_adjIndex[v_rt]];
        //如果报错，可以改成for循环遍历判断是否相等
        if (std::find(PMR + M_v_ps_start, PMR + M_v_ps_end + 1, curr) != PMR + M_v_ps_end + 1) {
            ans.emplace_back(curr);
        }
    }
    return ans;
}

//求Mtemp和Mp的交集，已针对大规模数据进行优化，返回值为交集是否为空，true交集为空，false交集非空
//FIXME,sky,此处不应该修改PMR，应该再单独存储结果
/*
bool PatternMatching::intersection(std::vector<R_ID>& Mtemp, R_ID& start, R_ID& end) {
    std::unordered_set<R_ID> temp(Mtemp.begin(), Mtemp.end());
    bool is_empty = true;
    for (unsigned i = start; i <= end; i++) {
        auto p = temp.find(PMR[i]);
        //TODO,sky,这里可以修改后将交集结果存储起来
        if (p != temp.end()) {
            temp.erase(PMR[i]);
            is_empty = false;
        }
        else {
            PMR[i] = NULL;
        }
    }
    return is_empty;
}
*/
bool PatternMatching::intersection(std::vector<R_ID>& Mtemp, R_ID& start, R_ID& end) {
    std::unordered_set<R_ID> temp(Mtemp.begin(), Mtemp.end());
    bool is_empty = true;

    std::cout << "Intersection set is: ";

    for (unsigned i = start; i <= end; ++i) {
        auto p = temp.find(PMR[i]);
        //TODO,sky,这里可以修改后将交集结果存储起来
        if (p != temp.end()) {
            std::cout << PMR[i] << " ";
            temp.erase(PMR[i]);
            is_empty = false;
        }
        //这里修改的话会对后续求解产生问题
        /*
        else {
            PMR[i] = NULL;
        }
        */
    }
    if (is_empty == true) {
        std::cout << "NULL";
    }
    std::cout << std::endl;

    return is_empty;
}

//M(v_p)的更新函数，输入v_p1和v_p2，用v_p1去更新v_p2
bool PatternMatching::updataPMR(P_ID& v_p1, P_ID& v_p2) {
    //先得到两个点的M(vp)索引范围
    unsigned vp1_left, vp1_right, vp2_left, vp2_right;
    if (!get_Mvp_Index(v_p1, vp1_left, vp1_right) || !get_Mvp_Index(v_p2, vp2_left, vp2_right)) {
        //std::cout << "Error: PatternMatching::intersection(P_ID& v_p1, P_ID& v_p2) can't get_Mvp_Index of v_p:" << v_p1 << ' ' << v_p2 << std::endl;
        return false;
    }
    for (unsigned i = vp2_left; i <= vp2_right; i++) {
        if (PMR[i] == NULL) {
            continue;
        }
        //std::vector<unsigned> temp(R_adj[PMR[i]], R_adj[PMR[i]] + degree_R[i].outdeg);//注意这里检查是否会报错，此处为利用vector构造函数使用指针范围初始化
        std::vector<unsigned> temp;
        for (int j = 0; j < degree_R[PMR[i]].outdeg; ++j) {
            temp.emplace_back(R_adj[R_adjIndex[PMR[i] + j]]);
        }
        //std::cout << "Using Vp" << v_p1 << " to updata Vp" << v_p2 << " ";
        //如果PMR[i]的邻接表与M(v_p1)交集为空，则删去PMR[i]点
        if (intersection(temp, vp1_left, vp1_right)) {
            //PMR[i] = NULL;
            if (sel[v_p2] > 0) sel[v_p2]--;//sky,20210906修改，存疑，之后检检查
        }
    }
    return true;
}

//给定源点vrs，计算与模式边vps, vpt端点vpt匹配的真实图顶点，一次扩展真实图中一顶点，返回是否更新成功
bool PatternMatching::extendEPatternwithsource(P_ID v_ps, P_ID v_pt, R_ID v_rs) {
    if (R_visited[v_rs] == 0) {
        std::vector<R_ID> Mtemp = getMV(v_ps, v_pt, v_rs);
        //从PMR_index中获取M(v_pt)的索引起止范围
        R_ID M_v_pt_start = 0, M_v_pt_end;
        if (!get_Mvp_Index(v_pt, M_v_pt_start, M_v_pt_end)) {
            //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
        }
        /*
        if (v_pt == 0) {
            M_v_pt_end = PMR_index[v_pt] - 1;
        }
        else {
            M_v_pt_start = PMR_index[v_pt - 1];
            M_v_pt_end = PMR_index[v_pt] - 1;
        }
        */
        //如果M(v_pt)不为空，即v_pt已访问,这里判断条件包含等号，因为相等时包含1个元素
        if (M_v_pt_end >= M_v_pt_start) {
            //std::cout << "When V_ps,V_pt,V_rs is: " << v_ps << ", " << v_pt << ", " << v_rs << std::endl;
            //这里对Mtemp和M(v_pt)求交集并缩减M(v_pt)，缩减过程中删去的部分赋值NULL
            if (intersection(Mtemp, M_v_pt_start, M_v_pt_end)) {
                return false;//交集得到的M(v_pt)为空 终止本轮搜索
            }
            else {
                //沿着已访问边逆向传播，这里考虑申请一个队列来完成，同时申请一个set来存储已插入过队列的节点ID
                std::queue<P_ID> que;
                std::unordered_set<P_ID> visited_set;
                que.push(v_ps);
                visited_set.insert(v_ps);
                P_ID currNode = v_pt;

                while (!que.empty()) {
                    for (P_ID i = 0; i < vertexNum_P; i++) {
                        //每一个结点的逆向向结点更新,这里加一个判断避免相同的点插入队列中
                        if (P_adj[i][currNode] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            //逆向传播缩减中间结果集过程
                            updataPMR(currNode, i);//逆向使用当前结点更新所有已访问的逆向边结点的真实图映射集合
                            //P_adj[i][currNode] = 3;//逆向更新过程，可能遇到一个环，会无限循环执行
                                //需要设置标记3代表两点有边且正向边，当前while循环中使用，while循环结束后变回2
                        }
                        //每一个结点的正向结点更新,这里加一个判断避免相同的点插入队列中
                        if (P_adj[currNode][i] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            updataPMR(currNode, i);
                            //updataPMR(i, currNode);
                            //P_adj[currNode][i] = 3;
                        }
                    }
                    currNode = que.front();
                    que.pop();
                }
                /*
                //逆向扩展的while循环结束，P_adj中标记为3的边全部变回2
                for (P_ID i = 0; i < vertexNum_P; i++) {
                    for (P_ID j = 0; j < vertexNum_P; j++) {
                        if (P_adj[i][j] == 3) {
                            P_adj[i][j] = 2;
                        }
                    }
                }
                */
            }
        }
        //v_pt未被访问时，直接更新为Mtemp
        else {
            //TODO 目前的程序思路不会出现这种情况
        }
    }
    return true;
}

//给定终点vrt，计算与模式边vps, vpt端点vps匹配的真实图顶点
bool PatternMatching::extendEPatternwithtarget(P_ID v_ps, P_ID v_pt, R_ID v_rt) {
    //TODO R_visited[]需要添加，已添加一个数组
    if (R_visited[v_rt] == 0) {
        std::vector<R_ID> Mtemp = reverse_getMV(v_ps, v_pt, v_rt);
        //从PMR_index中获取M(v_ps)的索引起止范围
        R_ID M_v_ps_start, M_v_ps_end;
        if (!get_Mvp_Index(v_pt, M_v_ps_start, M_v_ps_end)) {
            //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
        }

        //如果M(v_ps)不为空，即v_ps已访问，这里判断条件包含等号，因为相等时包含1个元素
        if (M_v_ps_end >= M_v_ps_start) {
            //std::cout << "When V_ps,V_pt,V_rt is: " << v_ps << ", " << v_pt << ", " << v_rt << std::endl;
            //这里对Mtemp和M(v_ps)求交集并缩减M(v_ps)，缩减过程中删去的部分赋值NULL
            if (intersection(Mtemp, M_v_ps_start, M_v_ps_end)) {
                return false;//交集得到的M(v_ps)为空 终止本轮搜索
            }
            else {
                //沿着已访问边逆向传播，这里的逆向并不是指正向边逆向边的逆向，考虑申请一个队列来完成
                std::queue<P_ID> que;
                std::unordered_set<P_ID> visited_set;
                que.push(v_pt);
                visited_set.insert(v_pt);
                P_ID currNode = v_ps;

                while (!que.empty()) {
                    for (P_ID i = 0; i < vertexNum_P; i++) {
                        //每一个结点的已访问正向结点更新这里加一个判断避免相同的点插入队列中
                        if (P_adj[currNode][i] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            updataPMR(currNode, i);
                            //P_adj[currNode][i] = 3;//正向更新过程，可能遇到一个环，会无限循环执行
                                    //需要设置标记3代表两点有边且正向边，当前while循环中使用，while循环结束后变回2
                        }
                        //每一个结点的已访问逆向结点更新,这里加一个判断避免相同的点插入队列中
                        if (P_adj[i][currNode] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            //逆向传播缩减中间结果集过程
                            updataPMR(currNode, i);
                            //P_adj[i][currNode] = 3;//逆向更新过程，可能遇到一个环，会无限循环执行
                                //需要设置标记3代表两点有边且正向边，当前while循环中使用，while循环结束后变回2
                        }
                    }
                    currNode = que.front();
                    que.pop();
                }
                /*
                //逆向扩展的while循环结束，P_adj中标记为3的边全部变回2
                for (P_ID i = 0; i < vertexNum_P; i++) {
                    for (P_ID j = 0; j < vertexNum_P; j++) {
                        if (P_adj[i][j] == 3) {
                            P_adj[i][j] = 2;
                        }
                    }
                }
                */
            }
        }
        //v_pt未被访问时，直接更新为Mtemp
        else {
            //TODO 目前的程序思路不会出现这种情况
        }
    }
    return true;
}

//计算与模式边vps, vpt端点vpt匹配的真实图顶点
void PatternMatching::extendEdgePattern(P_ID v_ps, P_ID v_pt) {
    //判断两点间是否存在正向边，有则标记并进行正向扩展
    if (P_adj[v_ps][v_pt] == 0) {
        //std::cout << "Error: can't reverse_extendEdgePattern, point " << v_pt << " and " << v_ps << " don't have edge." << std::endl;
        return;
    }
    else {
        P_adj[v_ps][v_pt] = 2; //2代表两点间有边且已访问
    }
    //下面开始从PMR中获取M(v_ps)
    R_ID M_v_ps_start, M_v_ps_end;
    if (!get_Mvp_Index(v_ps, M_v_ps_start, M_v_ps_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
    }
    /*
    if (v_ps == 0) {
        M_v_ps_end = PMR_index[v_ps] - 1;
    }
    else {
        M_v_ps_start = PMR_index[v_ps - 1];
        M_v_ps_end = PMR_index[v_ps] - 1;
    }
    */
    for (unsigned i = M_v_ps_start; i <= M_v_ps_end; i++) {
        R_ID v_rs = PMR[i];//按照PMR_index给出的索引逐一获取v_rs
        //匹配并逐步缩小中间结果集
        if (!extendEPatternwithsource(v_ps, v_pt, v_rs)) {
            //std::cout << "Error: can't extend E pattern with source (v_ps, v_pt, v_rs): " << v_ps << " " << v_pt << " " << v_rs << std::endl;
        }
    }
    return;
}

//逆向扩展，需在逆向存储中搜索
void PatternMatching::reverse_extendEdgePattern(P_ID v_pt, P_ID v_ps) {
    //判断两点间是否存在逆向边，有则标记并进行逆向扩展，注意这里P_adj[v_pt][v_ps]两个参数的顺序
    if (P_adj[v_pt][v_ps] == 0) {
        //std::cout << "Error: can't reverse_extendEdgePattern, point " << v_pt << " and " << v_ps << " don't have reverse edge." << std::endl;
        return;
    }
    else {
        P_adj[v_pt][v_ps] = 2; //2代表两点间有边且已访问
    }
    //下面开始从PMR中获取M(v_pt)
    R_ID M_v_pt_start = 0, M_v_pt_end;
    if (!get_Mvp_Index(v_pt, M_v_pt_start, M_v_pt_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
    }

    for (unsigned i = M_v_pt_start; i <= M_v_pt_end; i++) {
        R_ID V_rt = PMR[i];//按照PMR_index给出的索引逐一获取v_rt
        //匹配并逐步缩小中间结果集
        if (!extendEPatternwithtarget(v_ps, v_pt, V_rt)) {
            //std::cout << "Error: can't extend E pattern with target (v_ps, v_pt, V_rt): " << v_ps << " " << v_pt << " " << V_rt << std::endl;
        }
    }
    return;
}

//输出函数，将筛选出的结果集进行组合后以文本格式输出，输出文件名：out.txt
/*
void PatternMatching::output(string Dir) {
    ofstream ofile(Dir + "out.txt", ios::out);
}
*/

//打印PMR，用于测试
void PatternMatching::print_PMR() {
    std::cout << "Print current PMR collection." << std::endl;
    int row = 0;
    for (auto i : PMR) {
        std::cout << "V_p num: " << row << "    PMR nums: ";
        for (auto j : i) {
            std::cout << j << " ";
        }
        row++;
        std::cout << std::endl;
    }
    std::cout << "Finish print PMR collection" << std::endl;
}
