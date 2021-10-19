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
    PatternMatching() {};                       //�����䣬���캯��
    void init();
    virtual ~PatternMatching();              //�����䣬��������
    bool build_degree_R(std::string inputfile, unsigned vertexNum); //�õ�degree_R
    bool build_R_adj(std::string inputfile);   //�õ���һά����洢���ڽӱ�����ڽӱ�
    bool build_P_adj(std::string inputfile, unsigned vertexNum);//���ڽӾ���洢ģʽͼ
    //bool matchPR();                            //����P�и��ڵ���G�е�ƥ�伯������
                                                 //211018sky��Ҫ�޸�Ϊ����չʽ
    bool matchPR_expand();
    void searchPR();                           //��ģʽͼp�е�s��Ϊ��㣬Vr,sΪ���ĵ㿪ʼ����
    void searchAllPR();
    //void reverse_extendEdgePattern();
    P_ID getMaxSel();                             //��ȡѡ�������ģʽͼ��
    bool get_Mvp_Index(P_ID v_p, R_ID& start, R_ID& end);              //�õ�M(v_p)��������ֹ��Χ
    std::vector<R_ID> reverse_getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rt);    //�õ�ͼG����vrtΪ�յ㣬���ģʽvps, vpt��ƥ��Ķ˵㼯�ϣ�����vrtΪ�յ����vps��ƥ��Ķ��㼯
    void extendEdgePattern(P_ID v_ps, P_ID v_pt);             //����չ�߽���ƥ��
    bool extendEPatternwithsource(P_ID v_ps, P_ID v_pt, R_ID v_rs);
    void reverse_extendEdgePattern(P_ID v_pt, P_ID v_ps);
    bool extendEPatternwithtarget(P_ID v_ps, P_ID v_pt, R_ID V_rt);
    //P_ID getNextEPattern(P_ID& v_ps);                            //����δ�����ʵ�ƥ�䶥�������ٵķ�֧���յ�ID

    std::vector<R_ID> getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs); //�õ�ͼG����v_rsΪԴ�����ģʽ(v_ps, v_pt)��ƥ��Ķ��㼯
    bool intersection(std::vector<R_ID>& Mtemp, R_ID& start, R_ID& end);    //��Mtemp��Mpt�Ľ���������Դ��ģ���ݽ����Ż�������ֵΪ�����Ƿ�Ϊ��
    bool updataPMR(P_ID& v_p1, P_ID& v_p2);                              //M(v_p)�ĸ��º���������v_p1��v_p2����v_p1ȥ����v_p2�������˺���PatternMatching::intersection
    // void searchSplits(P_ID &v_ps);                                        //����ģʽͼ�н��v_ps�ĺ�̣����õ�ĸ���֧
    // void MergeSplit(P_ID &v_ps);                                          //�ϲ��м�������С�м�����
    bool isNextEPatternEmpty();                                  //�ж�ʣ�±�ģʽ�Ƿ�Ϊ��
    //P_ID getNextEPattern(P_ID& v_ps);                            //����δ�����ʵ�ƥ�䶥�������ٵķ�֧���յ�ID
    //std::vector<R_ID> getMV(P_ID& v_ps, P_ID& v_pt, R_ID& v_rs); //�õ�ͼG����v_rsΪԴ�����ģʽ(v_ps, v_pt)��ƥ��Ķ��㼯
    void FIXLINE(char* s);
    //void output(std::string Dir); //�����������ɸѡ���Ľ����������Ϻ����ı���ʽ���

    bool isfinish(P_ID id);
    void print_PMR();

    void searchPG(std::vector<std::vector<unsigned>>PMR_copy, std::vector<int>sel_copy, unsigned CenterPoint);
    P_ID getMaxSel_cur();

private:
    // std::vector<std::vector<int>> neighbor_p;//�洢ģʽͼ�и��ڵ���ͼ�ϵ��ھӼ���
    unsigned* R_visited; //�����ʵͼ�н���Ƿ��ѷ��ʣ������ų����ĵ�
    unsigned minMatchNum; //��Сƥ��Ľڵ���
    P_ID minMatchID;//��ʼ��Сƥ��ڵ�
    R_ID maxID;
    Degree* degree_P;
    Degree* degree_R;
    P_ID maxSelId; //��getNextEPattern()���õ�����¼ÿ��ѡ������Ķ���id
    std::vector<int> sel;//ÿ��ģʽͼ���ѡ��ȣ���¼��ÿ��ģʽͼ�����ʵ��ͼ�е�ƥ������
    unsigned vertexNum_R;
    unsigned vertexNum_P;
    unsigned edgeNum_R;
    unsigned edgeNum_P;
    // std::vector<std::vector<P_edge>> P_adj; //ģʽͼ���ڽӱ�
    std::vector<std::vector<P_ID>> P_adj;//ģʽͼ�ڽӾ���0��������֮��û�бߣ�1��������֮���б���δ�����ʹ���2��������֮���б����ѱ����ʹ�
    unsigned R_adjSize;                     //��СΪedgeNum_R*2+vertexNum_R*2
    std::vector<int> PMR_index; //ÿ��ģʽͼ����PMR����ʼƥ��λ��

    //unsigned* PMR;//��һά����洢M(Vp)
    //ʹ�ô�С�����PMR��չ����20211008
    std::vector<std::vector<unsigned>> PMR;

    //std::vector<vector<unsigned>> PMR; //��άָ������洢PMR

    unsigned* R_adj;//�ڽӱ�
    unsigned* R_reverse_adj;//���ڽӱ�
    unsigned* R_adjIndex;//�ڽӱ�from�ڵ����ڽӱ��е���ʼλ��
    unsigned* R_reverseAdjIndex;
    unsigned total_num;//�洢M(Vp)���ܵĸ�������get_Mvp_Index�л��õ�
};
#endif // !PATTERNMATCHING_H_
