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

bool PatternMatching::build_degree_R(std::string inputfile, unsigned vertexNum)//����degree_R�����Գɹ�
{
    maxID = 0;     //ͼ�����id
    edgeNum_R = 0; //ͼ�ı���
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
    char* s = (char*)malloc(maxlen); //��ʱ���һ�ζ��������
    size_t bytesread = 0;             //���������ֽ���
    char delims[] = " \t";            //�ַ����ķָ�����tab����
    size_t linenum = 0;               //����
    size_t lastlog = 0;               //���һ���ֽ�
    unsigned from = 0, to = 0;        //���id  //�յ�id
    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                  //�����з����'\0'
        char* t = strtok(s, delims); //���ص���tab֮ǰ����һ���ַ������׵�ַ
        from = atoi(t);              //��ȡ���id
        unsigned num = 0;            //���ȸ���
        while ((t = strtok(NULL, delims)) != NULL)
        {
            to = atoi(t);
            if (from != to)
            {
                maxID = max(to, maxID); //�ҳ�����ͼ�е����id
                degree_R[from].outdeg++;       //����
                degree_R[to].indeg++;          //���
                num++;
            }
        }
        edgeNum_R += num; //ͳ������from���ܵĳ��ȸ��� == �ܵı���
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
bool PatternMatching::build_R_adj(std::string inputfile)//�����ڽӱ�R_adj�����ڽӱ�R_reverse_adj�����Գɹ�
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
    char* s = (char*)malloc(maxlen); //��ʱ���һ�ζ��������
    char delims[] = " \t";            //�ַ����ķָ�����tab����
    size_t linenum = 0;               //����
    unsigned from = 0, to = 0;        //���id  //�յ�id
    for (unsigned i = 1; i < vertexNum_R; i++) {
        R_adjIndex[i] = R_adjIndex[i - 1] + degree_R[i - 1].outdeg;
        R_reverseAdjIndex[i] = R_reverseAdjIndex[i - 1] + degree_R[i - 1].indeg;
        R_reverseAdjIndex_tail[i] = R_reverseAdjIndex[i];
    }

    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                                      //�����з����'\0'
        char* t = strtok(s, delims);                     //���ص���tab֮ǰ����һ���ַ������׵�ַ
        from = atoi(t);                                  //��ȡ���id
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
    std::vector<std::vector<P_ID>> tmp(vertexNum_P, std::vector<P_ID>(vertexNum_P));//��ʱ�洢�������ʼ��
    FILE* inf = fopen(inputfile.c_str(), "r");
    if (inf == NULL)
    {
        std::cerr << "Could not load :" << inputfile << " error: " << strerror(errno)
            << std::endl;
    }
    assert(inf != NULL);
    std::cout << "Reading in adjacency list format!" << std::endl;
    int maxlen = 100000000;
    char* s = (char*)malloc(maxlen); //��ʱ���һ�ζ��������
    char delims[] = " \t";            //�ַ����ķָ�����tab����
    size_t linenum = 0;               //����
    unsigned from = 0, to = 0;        //���id  //�յ�id
    while (fgets(s, maxlen, inf) != NULL)
    {
        linenum++;
        FIXLINE(s);                                      //�����з����'\0'
        char* t = strtok(s, delims);                     //���ص���tab֮ǰ����һ���ַ������׵�ַ
        from = atoi(t);                                  //��ȡ���id
        while ((t = strtok(NULL, delims)) != NULL)
        {
            to = atoi(t);
            if (from != to)
            {
                degree_P[from].outdeg++;       //����
                degree_P[to].indeg++;          //���
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
    sel.resize(vertexNum_P);//��ʼ��ѡ���
    minMatchID = 0; //������Сƥ������ģʽͼ�ڵ�ID
    minMatchNum = UINT_MAX;  //��Сƥ����
    unsigned current_size = 0;
    unsigned total_size = 0;//�洢M(Vp)����Ŀռ��С

    for (unsigned i = 0; i < vertexNum_P; i++)//������ÿ��ģʽͼ��ƥ���ʵ��ͼ�нڵ�����
    {
        for (R_ID j = 0; j < vertexNum_R; j++) //jҪ�޸ĳ�R�������ڵĽڵ�id����
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

    //ֻ������Сƥ��ģʽͼ����PMR���ϣ��������ֿյ�״̬
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

//bool PatternMatching::matchPR() //������ʼƥ�伯��P_M_R,���Գɹ�
//{
//    total_num = 0;
//    sel.resize(vertexNum_P);//��ʼ��ѡ���
//    minMatchID = 0; //������Сƥ������ģʽͼ�ڵ�ID
//    minMatchNum = UINT_MAX;  //��Сƥ����
//    unsigned current_size = 0;
//    unsigned total_size = 0;//�洢M(Vp)����Ŀռ��С
//
//    for (unsigned i = 0; i < vertexNum_P; i++)//������ÿ��ģʽͼ��ƥ���ʵ��ͼ�нڵ�����
//    {
//        for (R_ID j = 0; j < vertexNum_R; j++) //jҪ�޸ĳ�R�������ڵĽڵ�id����
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
//    for (unsigned i = 0; i < vertexNum_P; i++)//������ÿ��ģʽͼ��ƥ���ʵ��ͼ�нڵ�����
//    {
//        for (R_ID j = 0; j < vertexNum_R; j++) //jҪ�޸ĳ�R�������ڵĽڵ�id����
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
//    for (; curid < vertexNum_P; curid++) {//�ҵ�id��С��ģʽ��û�б�ȫ�����ʵĵ�
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
//        if (sel[i] < sel[maxSelID] && !isfinish(i)) {  //ƥ����ԽС��ģʽͼ�ڵ�ѡ���Խ��
//            maxSelID = i;
//        }
//    }
//    return maxSelID;
//}


//�����ǵõ���ǰѡ��ȱ���ƥ�������ٵ�ģʽͼ���Ų�����
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
    //for (; curid < vertexNum_P; curid++) {//�ҵ�id��С��ģʽ��û�б�ȫ�����ʵĵ�
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
    //    if (sel[i] < sel[maxSelID] && !isfinish(i)) {  //ƥ����ԽС��ģʽͼ�ڵ�ѡ���Խ��
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
        R_visited[minMatchID_PMR[i]] = 1;  //���vr,sΪ�������ĵ�

        //FIXME,�����������Ż���sel�ڴ�ǰ�ĳ�ʼ��������û�б�Ҫ��
        //Initializeģʽͼÿ�����ѡ���Sel;����vp,s֮�⣬����ģʽͼ����ѡ���Ϊ�����
        /*for (int i = 0; i < vertexNum_P; ++i) {
            if (i != minMatchID) {
                sel[i] = MAXINT;
            }
        }*/

        //������Ҫ����PMR��sel���������Ϊ��Ҫ�������ƴ��ݵ�searchPG();
        //searchPR();        
        searchPG(PMR, sel, minMatchID_PMR[i]);
    }
    /*for (unsigned i = PMR_index[minMatchID],j=0; j<sel[minMatchID]; j++,i++) {
        searchPR();
        R_visited[i] = 1;
    }*/
}

void PatternMatching::searchPG(std::vector<std::vector<unsigned>>PMR_copy, std::vector<int>sel_copy, unsigned CenterPoint) {
    //Initializeģʽͼÿ�����ѡ���Sel;����vp,s֮�⣬����ģʽͼ����ѡ���Ϊ�����
    for (int i = 0; i < vertexNum_P; ++i) {
        if (i != minMatchID) {
            sel_copy[i] = MAXINT;
        }
    }
    while (!isNextEPatternEmpty()) {
        P_ID maxSelId = getMaxSel_cur();
        //P_ID neighborID = 0; //���������㹹����Сƥ��ı�ģʽ
        unsigned minMatchNum = INT_MAX;//neighborID����Сƥ����
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

//ԭ���searchPR������Ŀǰ�Ѳ�ʹ�ã�������ʹ��searchPG����
void PatternMatching::searchPR() {
    while (!isNextEPatternEmpty()) {
        P_ID maxSelId = getMaxSel();
        P_ID neighborID = 0; //���������㹹����Сƥ��ı�ģʽ
        unsigned minMatchNum = INT_MAX;//neighborID����Сƥ����
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

//�洢��ֵ�����汾���õ�M(v_p)��������ֹ��Χ,ɸ��β��NULLֵ���ɹ��򷵻�true�Լ�start��endֵ��δ�ҵ�����false
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

//�洢βֵ�����汾���õ�M(v_p)��������ֹ��Χ,�ҵ�����true�Լ�start��endֵ��δ�ҵ�����false
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

//�õ�ͼG����vrsΪԴ�����ģʽvps, vpt��ƥ��Ķ˵㼯�ϣ�����vrsΪԴ����vpt��ƥ��Ķ��㼯
std::vector<R_ID> PatternMatching::getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs) {
    int length = degree_R[v_rs].outdeg;
    std::vector<R_ID> ans;

    //��PMR_index�л�ȡM(v_pt)��������ֹ��Χ
    R_ID M_v_pt_start, M_v_pt_end;
    if (!get_Mvp_Index(v_pt, M_v_pt_start, M_v_pt_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
    }

    //��ȡR_adj��v_rs��һ�е����г��Ƚ�㣬�������ж��Ƿ�����M(v_pt)
    for (int j = 0; j < length; j++) {
        R_ID curr = R_adj[j + R_adjIndex[v_rs]];
        //����������Ըĳ�forѭ�������ж��Ƿ����
        if (std::find(PMR + M_v_pt_start, PMR + M_v_pt_end + 1, curr) != PMR + M_v_pt_end + 1) {
            ans.emplace_back(curr);
        }
    }
    return ans;
}

//�õ�ͼG����vrtΪ�յ㣬���ģʽvps, vpt��ƥ��Ķ˵㼯�ϣ�����vrtΪ�յ����vps��ƥ��Ķ��㼯
std::vector<R_ID> PatternMatching::reverse_getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rt) {
    int length = degree_R[v_rt].outdeg;
    std::vector<R_ID> ans;

    //��PMR_index�л�ȡM(v_ps)��������ֹ��Χ
    R_ID M_v_ps_start, M_v_ps_end;
    if (!get_Mvp_Index(v_ps, M_v_ps_start, M_v_ps_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_ps << std::endl;
    }

    //��ȡR_adj��v_rt��һ�е����г��Ƚ�㣬�������ж��Ƿ�����M(v_ps)
    for (int j = 0; j < length; j++) {
        R_ID curr = R_adj[j + R_adjIndex[v_rt]];
        //����������Ըĳ�forѭ�������ж��Ƿ����
        if (std::find(PMR + M_v_ps_start, PMR + M_v_ps_end + 1, curr) != PMR + M_v_ps_end + 1) {
            ans.emplace_back(curr);
        }
    }
    return ans;
}

//��Mtemp��Mp�Ľ���������Դ��ģ���ݽ����Ż�������ֵΪ�����Ƿ�Ϊ�գ�true����Ϊ�գ�false�����ǿ�
//FIXME,sky,�˴���Ӧ���޸�PMR��Ӧ���ٵ����洢���
/*
bool PatternMatching::intersection(std::vector<R_ID>& Mtemp, R_ID& start, R_ID& end) {
    std::unordered_set<R_ID> temp(Mtemp.begin(), Mtemp.end());
    bool is_empty = true;
    for (unsigned i = start; i <= end; i++) {
        auto p = temp.find(PMR[i]);
        //TODO,sky,��������޸ĺ󽫽�������洢����
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
        //TODO,sky,��������޸ĺ󽫽�������洢����
        if (p != temp.end()) {
            std::cout << PMR[i] << " ";
            temp.erase(PMR[i]);
            is_empty = false;
        }
        //�����޸ĵĻ���Ժ�������������
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

//M(v_p)�ĸ��º���������v_p1��v_p2����v_p1ȥ����v_p2
bool PatternMatching::updataPMR(P_ID& v_p1, P_ID& v_p2) {
    //�ȵõ��������M(vp)������Χ
    unsigned vp1_left, vp1_right, vp2_left, vp2_right;
    if (!get_Mvp_Index(v_p1, vp1_left, vp1_right) || !get_Mvp_Index(v_p2, vp2_left, vp2_right)) {
        //std::cout << "Error: PatternMatching::intersection(P_ID& v_p1, P_ID& v_p2) can't get_Mvp_Index of v_p:" << v_p1 << ' ' << v_p2 << std::endl;
        return false;
    }
    for (unsigned i = vp2_left; i <= vp2_right; i++) {
        if (PMR[i] == NULL) {
            continue;
        }
        //std::vector<unsigned> temp(R_adj[PMR[i]], R_adj[PMR[i]] + degree_R[i].outdeg);//ע���������Ƿ�ᱨ���˴�Ϊ����vector���캯��ʹ��ָ�뷶Χ��ʼ��
        std::vector<unsigned> temp;
        for (int j = 0; j < degree_R[PMR[i]].outdeg; ++j) {
            temp.emplace_back(R_adj[R_adjIndex[PMR[i] + j]]);
        }
        //std::cout << "Using Vp" << v_p1 << " to updata Vp" << v_p2 << " ";
        //���PMR[i]���ڽӱ���M(v_p1)����Ϊ�գ���ɾȥPMR[i]��
        if (intersection(temp, vp1_left, vp1_right)) {
            //PMR[i] = NULL;
            if (sel[v_p2] > 0) sel[v_p2]--;//sky,20210906�޸ģ����ɣ�֮�����
        }
    }
    return true;
}

//����Դ��vrs��������ģʽ��vps, vpt�˵�vptƥ�����ʵͼ���㣬һ����չ��ʵͼ��һ���㣬�����Ƿ���³ɹ�
bool PatternMatching::extendEPatternwithsource(P_ID v_ps, P_ID v_pt, R_ID v_rs) {
    if (R_visited[v_rs] == 0) {
        std::vector<R_ID> Mtemp = getMV(v_ps, v_pt, v_rs);
        //��PMR_index�л�ȡM(v_pt)��������ֹ��Χ
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
        //���M(v_pt)��Ϊ�գ���v_pt�ѷ���,�����ж����������Ⱥţ���Ϊ���ʱ����1��Ԫ��
        if (M_v_pt_end >= M_v_pt_start) {
            //std::cout << "When V_ps,V_pt,V_rs is: " << v_ps << ", " << v_pt << ", " << v_rs << std::endl;
            //�����Mtemp��M(v_pt)�󽻼�������M(v_pt)������������ɾȥ�Ĳ��ָ�ֵNULL
            if (intersection(Mtemp, M_v_pt_start, M_v_pt_end)) {
                return false;//�����õ���M(v_pt)Ϊ�� ��ֹ��������
            }
            else {
                //�����ѷ��ʱ����򴫲������￼������һ����������ɣ�ͬʱ����һ��set���洢�Ѳ�������еĽڵ�ID
                std::queue<P_ID> que;
                std::unordered_set<P_ID> visited_set;
                que.push(v_ps);
                visited_set.insert(v_ps);
                P_ID currNode = v_pt;

                while (!que.empty()) {
                    for (P_ID i = 0; i < vertexNum_P; i++) {
                        //ÿһ�����������������,�����һ���жϱ�����ͬ�ĵ���������
                        if (P_adj[i][currNode] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            //���򴫲������м���������
                            updataPMR(currNode, i);//����ʹ�õ�ǰ�����������ѷ��ʵ�����߽�����ʵͼӳ�伯��
                            //P_adj[i][currNode] = 3;//������¹��̣���������һ������������ѭ��ִ��
                                //��Ҫ���ñ��3���������б�������ߣ���ǰwhileѭ����ʹ�ã�whileѭ����������2
                        }
                        //ÿһ���������������,�����һ���жϱ�����ͬ�ĵ���������
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
                //������չ��whileѭ��������P_adj�б��Ϊ3�ı�ȫ�����2
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
        //v_ptδ������ʱ��ֱ�Ӹ���ΪMtemp
        else {
            //TODO Ŀǰ�ĳ���˼·��������������
        }
    }
    return true;
}

//�����յ�vrt��������ģʽ��vps, vpt�˵�vpsƥ�����ʵͼ����
bool PatternMatching::extendEPatternwithtarget(P_ID v_ps, P_ID v_pt, R_ID v_rt) {
    //TODO R_visited[]��Ҫ��ӣ������һ������
    if (R_visited[v_rt] == 0) {
        std::vector<R_ID> Mtemp = reverse_getMV(v_ps, v_pt, v_rt);
        //��PMR_index�л�ȡM(v_ps)��������ֹ��Χ
        R_ID M_v_ps_start, M_v_ps_end;
        if (!get_Mvp_Index(v_pt, M_v_ps_start, M_v_ps_end)) {
            //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
        }

        //���M(v_ps)��Ϊ�գ���v_ps�ѷ��ʣ������ж����������Ⱥţ���Ϊ���ʱ����1��Ԫ��
        if (M_v_ps_end >= M_v_ps_start) {
            //std::cout << "When V_ps,V_pt,V_rt is: " << v_ps << ", " << v_pt << ", " << v_rt << std::endl;
            //�����Mtemp��M(v_ps)�󽻼�������M(v_ps)������������ɾȥ�Ĳ��ָ�ֵNULL
            if (intersection(Mtemp, M_v_ps_start, M_v_ps_end)) {
                return false;//�����õ���M(v_ps)Ϊ�� ��ֹ��������
            }
            else {
                //�����ѷ��ʱ����򴫲�����������򲢲���ָ���������ߵ����򣬿�������һ�����������
                std::queue<P_ID> que;
                std::unordered_set<P_ID> visited_set;
                que.push(v_pt);
                visited_set.insert(v_pt);
                P_ID currNode = v_ps;

                while (!que.empty()) {
                    for (P_ID i = 0; i < vertexNum_P; i++) {
                        //ÿһ�������ѷ�����������������һ���жϱ�����ͬ�ĵ���������
                        if (P_adj[currNode][i] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            updataPMR(currNode, i);
                            //P_adj[currNode][i] = 3;//������¹��̣���������һ������������ѭ��ִ��
                                    //��Ҫ���ñ��3���������б�������ߣ���ǰwhileѭ����ʹ�ã�whileѭ����������2
                        }
                        //ÿһ�������ѷ������������,�����һ���жϱ�����ͬ�ĵ���������
                        if (P_adj[i][currNode] == 2 && visited_set.count(i) == 0) {
                            que.push(i);
                            visited_set.insert(i);
                            //���򴫲������м���������
                            updataPMR(currNode, i);
                            //P_adj[i][currNode] = 3;//������¹��̣���������һ������������ѭ��ִ��
                                //��Ҫ���ñ��3���������б�������ߣ���ǰwhileѭ����ʹ�ã�whileѭ����������2
                        }
                    }
                    currNode = que.front();
                    que.pop();
                }
                /*
                //������չ��whileѭ��������P_adj�б��Ϊ3�ı�ȫ�����2
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
        //v_ptδ������ʱ��ֱ�Ӹ���ΪMtemp
        else {
            //TODO Ŀǰ�ĳ���˼·��������������
        }
    }
    return true;
}

//������ģʽ��vps, vpt�˵�vptƥ�����ʵͼ����
void PatternMatching::extendEdgePattern(P_ID v_ps, P_ID v_pt) {
    //�ж�������Ƿ��������ߣ������ǲ�����������չ
    if (P_adj[v_ps][v_pt] == 0) {
        //std::cout << "Error: can't reverse_extendEdgePattern, point " << v_pt << " and " << v_ps << " don't have edge." << std::endl;
        return;
    }
    else {
        P_adj[v_ps][v_pt] = 2; //2����������б����ѷ���
    }
    //���濪ʼ��PMR�л�ȡM(v_ps)
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
        R_ID v_rs = PMR[i];//����PMR_index������������һ��ȡv_rs
        //ƥ�䲢����С�м�����
        if (!extendEPatternwithsource(v_ps, v_pt, v_rs)) {
            //std::cout << "Error: can't extend E pattern with source (v_ps, v_pt, v_rs): " << v_ps << " " << v_pt << " " << v_rs << std::endl;
        }
    }
    return;
}

//������չ����������洢������
void PatternMatching::reverse_extendEdgePattern(P_ID v_pt, P_ID v_ps) {
    //�ж�������Ƿ��������ߣ������ǲ�����������չ��ע������P_adj[v_pt][v_ps]����������˳��
    if (P_adj[v_pt][v_ps] == 0) {
        //std::cout << "Error: can't reverse_extendEdgePattern, point " << v_pt << " and " << v_ps << " don't have reverse edge." << std::endl;
        return;
    }
    else {
        P_adj[v_pt][v_ps] = 2; //2����������б����ѷ���
    }
    //���濪ʼ��PMR�л�ȡM(v_pt)
    R_ID M_v_pt_start = 0, M_v_pt_end;
    if (!get_Mvp_Index(v_pt, M_v_pt_start, M_v_pt_end)) {
        //std::cout << "Error: Can't get_Mvp_Index of " << v_pt << std::endl;
    }

    for (unsigned i = M_v_pt_start; i <= M_v_pt_end; i++) {
        R_ID V_rt = PMR[i];//����PMR_index������������һ��ȡv_rt
        //ƥ�䲢����С�м�����
        if (!extendEPatternwithtarget(v_ps, v_pt, V_rt)) {
            //std::cout << "Error: can't extend E pattern with target (v_ps, v_pt, V_rt): " << v_ps << " " << v_pt << " " << V_rt << std::endl;
        }
    }
    return;
}

//�����������ɸѡ���Ľ����������Ϻ����ı���ʽ���������ļ�����out.txt
/*
void PatternMatching::output(string Dir) {
    ofstream ofile(Dir + "out.txt", ios::out);
}
*/

//��ӡPMR�����ڲ���
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
