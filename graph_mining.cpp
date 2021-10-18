#include <ctime>
#include "PatternMatching.h"

//#include "OSFile.h"
/*
void handleSegmentFault(int sig) {
    void* array[1024];
    size_t size;
    size = backtrace(array, sizeof(array) / sizeof(void*));
    backtrace_symbols_fd(array, size, fileno(stderr));
    abort();
    return;
}
*/
int main() {
    /*
    if (argc != 6) {
        fprintf(stderr, "Usage: %s <RDataset> <PDataset> <answer Directory> <RVertex Number> <PVertex Number>\n", argv[0]);//L:./bin/lrelease/graph_builder
                                                                                               //   Dataset/amazon-2008.txt
                                                                                               //   amazon/  //�����������Ŀ¼����������Ŀ¼
                                                                                               //   735322
        return -1;
    }
    */
    std::string R_input_file = "D:\\TestMining\\20210728TestMining\\R_input.txt";
    std::string P_input_file = "D:\\TestMining\\20210728TestMining\\P_input.txt";
    unsigned int R_vnums = 5;
    unsigned int P_vnums = 3;
    //signal(SIGSEGV, handleSegmentFault);

    //std::string Dir = argv[3];//datasbase directory�������������Ŀ¼��  //L:amazon/
    /*
    if (OSFile::DirectoryExists(argv[3]) == false) {//���dataset directory�����ڣ����Լ��������Ŀ¼
        OSFile::MkDir(argv[3]);
    }
    */

    //��¼������ʼʱ��
    //struct timeval start_time, end_time, tmp_start, tmp_end;
    //gettimeofday(&start_time, NULL);//��ʼʱ��
	PatternMatching* pm = new PatternMatching();
	pm->build_degree_R(R_input_file, R_vnums);
	pm->build_R_adj(R_input_file);
	pm->build_P_adj(P_input_file, P_vnums);
	pm->init();
	//pm->matchPR(); //211018sky��Ҫ�޸ĸú�������Ϊ����չʽ
    pm->matchPR_expand();
	pm->print_PMR();
    clock_t start, end;
    start = clock();
	pm->searchAllPR();
    end = clock();
    std::cout << "PatternMatching run time is: " << (end - start) / CLK_TCK * 1000 << "ms" << std::endl;
	//pm->print_PMR();
	//gettimeofday(&tmp_end, NULL);
	//cout << "graph mining time elapse:" << ((tmp_end.tv_sec - start_time.tv_sec) * 1000000 + (tmp_end.tv_usec - start_time.tv_usec)) / 1000000.0 << " s" << endl;
	return 0;
}
