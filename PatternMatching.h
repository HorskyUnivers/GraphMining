#ifndef PATTERNMATCHING_H_
#define PATTERNMATCHING_H_

#include <iostream>
//#include "HashID.h"
//#include "ThreadPool.h"
#include <vector>
#include <string>
#include <limits.h>
#include <windows.h>
#include <unordered_map>
#include <map>
#include <list>
//#include "BuildAdjGraph.h"
//#include "Bitmap.h"
#include <assert.h>
#include <fstream>
#include <algorithm>
#include <limits.H>
#include <unordered_set>
#include <queue>
//#include <queue>
#include <fstream>

typedef unsigned int P_ID;
typedef unsigned int R_ID;
struct Degree
{
    unsigned indeg;
    unsigned outdeg;
};

struct P_edge
{
    int to;
    int flag;
};

class PatternMatching
{
public:
    PatternMatching() {};                       //待补充，构造函数
    void init();
    virtual ~PatternMatching();              //待补充，析构函数
    bool build_degree_R(std::string inputfile, unsigned vertexNum); //得到degree_R
    bool build_R_adj(std::string inputfile);   //得到用一维数组存储的邻接表和逆邻接表
    bool build_P_adj(std::string inputfile, unsigned vertexNum);//用邻接矩阵存储模式图
    //bool matchPR();                            //计算P中各节点在G中的匹配集并返回
                                                 //211018sky需要修改为逐步拓展式
    bool matchPR_expand();
    void searchPR();                           //以模式图p中的s点为起点，Vr,s为中心点开始搜索
    void searchAllPR();
    //void reverse_extendEdgePattern();
    P_ID getMaxSel();                             //获取选择度最大的模式图点
    bool get_Mvp_Index(P_ID v_p, R_ID& start, R_ID& end);              //得到M(v_p)的索引起止范围
    std::vector<R_ID> reverse_getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rt);    //得到图G上以vrt为终点，与边模式vps, vpt相匹配的端点集合，即以vrt为终点点与vps相匹配的顶点集
    void extendEdgePattern(P_ID v_ps, P_ID v_pt);             //逐步扩展边进行匹配
    bool extendEPatternwithsource(P_ID v_ps, P_ID v_pt, R_ID v_rs);
    void reverse_extendEdgePattern(P_ID v_pt, P_ID v_ps);
    bool extendEPatternwithtarget(P_ID v_ps, P_ID v_pt, R_ID V_rt);
    //P_ID getNextEPattern(P_ID& v_ps);                            //返回未被访问的匹配顶点数最少的分支边终点ID

    std::vector<R_ID> getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs); //得到图G上以v_rs为源，与边模式(v_ps, v_pt)相匹配的顶点集
    bool intersection(std::vector<R_ID>& Mtemp, R_ID& start, R_ID& end);    //求Mtemp和Mpt的交集，已针对大规模数据进行优化，返回值为交集是否为空
    bool updataPMR(P_ID& v_p1, P_ID& v_p2);                              //M(v_p)的更新函数，输入v_p1和v_p2，用v_p1去更新v_p2，重载了函数PatternMatching::intersection
    // void searchSplits(P_ID &v_ps);                                        //遍历模式图中结点v_ps的后继，即该点的各分支
    // void MergeSplit(P_ID &v_ps);                                          //合并中间结果，缩小中间结果集
    bool isNextEPatternEmpty();                                  //判断剩下边模式是否为空
    //P_ID getNextEPattern(P_ID& v_ps);                            //返回未被访问的匹配顶点数最少的分支边终点ID
    //std::vector<R_ID> getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs); //得到图G上以v_rs为源，与边模式(v_ps, v_pt)相匹配的顶点集
    void FIXLINE(char* s);
    //void output(std::string Dir); //输出函数，将筛选出的结果集进行组合后以文本格式输出

    bool isfinish(P_ID id);
    void print_PMR();

    void searchPG(std::vector<std::vector<unsigned>>PMR_copy, std::vector<int>sel_copy, unsigned CenterPoint);
    P_ID getMaxSel_cur();

private:
    // std::vector<std::vector<int>> neighbor_p;//存储模式图中各节点在图上的邻居集合
    unsigned* R_visited; //标记真实图中结点是否已访问，用于排除中心点
    unsigned minMatchNum; //最小匹配的节点数
    P_ID minMatchID;//初始最小匹配节点
    R_ID maxID;
    Degree* degree_P;
    Degree* degree_R;
    P_ID maxSelId; //在getNextEPattern()中用到，记录每轮选择度最大的顶点id
    std::vector<int> sel;//每个模式图点的选择度，记录了每个模式图点的在实际图中的匹配数量
    unsigned vertexNum_R;
    unsigned vertexNum_P;
    unsigned edgeNum_R;
    unsigned edgeNum_P;
    // std::vector<std::vector<P_edge>> P_adj; //模式图的邻接表
    std::vector<std::vector<P_ID>> P_adj;//模式图邻接矩阵，0代表两点之间没有边，1代表两点之间有边且未被访问过，2代表两点之间有边且已被访问过
    unsigned R_adjSize;                     //大小为edgeNum_R*2+vertexNum_R*2
    std::vector<int> PMR_index; //每个模式图点在PMR中起始匹配位置

    //unsigned* PMR;//用一维数组存储M(Vp)
    //使用从小到大的PMR拓展过程20211008
    std::vector<std::vector<unsigned>> PMR;

    //std::vector<vector<unsigned>> PMR; //二维指针数组存储PMR

    unsigned* R_adj;//邻接表
    unsigned* R_reverse_adj;//逆邻接表
    unsigned* R_adjIndex;//邻接表from节点在邻接表中的起始位置
    unsigned* R_reverseAdjIndex;
    unsigned total_num;//存储M(Vp)的总的个数，在get_Mvp_Index中会用到
};
#endif // !PATTERNMATCHING_H_
