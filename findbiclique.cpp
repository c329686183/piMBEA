/*
默认处理的二分图，左边节点集合为U，右边节点集合为V
*/

#include <stdio.h>
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
#include <iostream>
#include <fstream>  
#include <list>
#include <set>
#include <map>
#include <vector>
#include <algorithm>
#include <cstring>
#include <iterator>
#include <stack>
#include <omp.h>
#include <time.h>
#include <sstream>
#include <windows.h>
#include <sched.h> 
#include <time.h>

using namespace std;
/*线程参数*/
#define THREAD_NUM 1
#define TASK_THRESHOLD_VALUE 2//线程任务量阀值，比它大的才能负载均衡调度
#define PC_number 2//为了以后多PC 的scalable，对二分图分割

long R_NUM = 1;//为了计算p-q-biclique
long L_NUM = 1;
bool flag_p_q_biclique = false;
bool flag_loadbalance = true;//负载均衡标志位，默认关闭
bool flag_is_sparse = true;//稀疏图标志
bool flag_buffer_interval_value = true;//缓存中间结果标志位
bool flag_stack_monitor = true;//stack状态监控标志
bool flag_level_monitor = false;//层级状态监控标志
bool flag_read_monitor = false;//读数据状态监控标志
bool flag_VP_monitor = true;//读数据状态监控标志
bool flag_success_demo_1 = false;
bool flag_success_demo_2 = false;
bool flag_success_demo_3 = false;
bool flag_success_demo_4 = false;
bool flag_success_demo_5 = false;
long read_line_number = 0;
pthread_mutex_t mutex_thread[THREAD_NUM];//互斥变量
pthread_t stack_monitor_thread, level_monitor_thread, read_monitor_thread, VP_monitor_thread;
pthread_cond_t stack_monitor_cond, level_monitor_cond, read_monitor_cond, VP_monitor_cond;
pthread_mutex_t stack_monitor_mutex, level_monitor_mutex, read_monitor_mutex, VP_monitor_mutex;

/*文件读写变量*/

ofstream ofile;//用于输出结果文件
ofstream ofile_stack;//用于输出栈监控文件
ofstream ofile_level;//用于输出层级监控文件
ofstream ofile_read;//用于输出读数据监控文件
ofstream ofile_VP;//用于输出读数据监控文件
ifstream ifile;//用于读取原始图数据文件
ifstream ifile_info;//用于读取图数据的节点数量和边数量文件
string str_dataName;
string filePath_data = "e:\\bigraph-data\\dataName.txt";//图数据
string filePath_data_info = "e:\\bigraph-data\\dataName-info.txt";//图点和边数信息
string output_file = "e:\\bigraph-data\\dataName-result-parall-thread-"+std::to_string(THREAD_NUM)+".txt";//文件输出
string output_file_stack = "e:\\bigraph-data\\dataName-result-stack-thread-" + std::to_string(THREAD_NUM) + ".txt";//文件输出
string output_file_level = "e:\\bigraph-data\\dataName-result-level-thread-" + std::to_string(THREAD_NUM) + ".txt";//文件输出
string output_file_read = "e:\\bigraph-data\\dataName-result-read-thread-" + std::to_string(THREAD_NUM) + ".txt";//文件输出
string output_file_VP = "e:\\bigraph-data\\dataName-result-VP-thread-" + std::to_string(THREAD_NUM) + ".txt";//文件输出

/*图数据结构*/

typedef set<long> Type_v_sets;//普通节点集合
typedef vector<Type_v_sets> Type_Graph;//图数据类型

/*图全局变量*/

Type_Graph graph;//原始图数据，邻接表表示
Type_Graph p_q_graph;//原始图数据，邻接表表示
map<long, long> p_q_graph_map;//原始图数据，邻接表表示
map<long,long> p_q_graph_map_reverse;//原始图数据，邻接表表示
vector<long> degree;
long node_num;//总节点数
long p_q_node_num;//总节点数
long edge_num;//总边数

/*二分图数据结构*/
typedef vector<long> Type_ordered_sets;//有序集合
typedef vector<pair<long, bool> > Type_P;//P的类型有点特殊，为了保证有序遍历C，bool为false表示静止访问
typedef pair<int, pair<Type_v_sets, Type_v_sets> > Type_Biclique;//极大二分clique的类型,包含线程信息
Type_v_sets U;
long* OV;//V的有序集合
long* OVP;//记录节点经过排序后对应的原先的位置
Type_P initial_P;//初始化有序P，即OV加上flag，默认为1
Type_v_sets new_L_neibors[THREAD_NUM];//得到N[L']
Type_v_sets SQ[THREAD_NUM];//根据稀疏图特性选择的精简Q集合
Type_ordered_sets SP[THREAD_NUM];//根据稀疏图特性选择的精简P集合
Type_ordered_sets SPP[THREAD_NUM];//small_info中的元素对应在P中中的下标position
long U_num;//初始化总的U节点数
long p_q_U_num;//初始化经过q_p_core剪后总的p_q_U节点数
long V_num;//初始化总的V节点数
long p_q_V_num;//初始化经过q_p_core剪总的p_q_V节点数

struct TASK
{
	long x_index;//P(0)中元素的下标
	int level;//层级
};

/*线程全局数组*/
vector<Type_P> P[THREAD_NUM];
vector<vector<Type_v_sets> > PN[THREAD_NUM];//记录生成P'的过程中产生的邻居，用于下一层生成L'
vector<Type_ordered_sets> Q[THREAD_NUM];
vector<Type_v_sets> L[THREAD_NUM];
vector<Type_v_sets> R[THREAD_NUM];
stack<TASK> task_stack[THREAD_NUM];
Type_ordered_sets task_init[THREAD_NUM];//每个线程经过分配之后的初始化任务节点
int thread_ids[THREAD_NUM];
TASK nt[THREAD_NUM];//尝试获取栈中下一个元素
int nl[THREAD_NUM];//下一个元素层级
long ni[THREAD_NUM];//下一个元素层级
bool flag_maximal[THREAD_NUM];

set<Type_Biclique> solution[THREAD_NUM];

/*计时变量*/

time_t start_time, end_time;//全局计时开始、结束时间
time_t thread_stop_time[THREAD_NUM];//thread结束时间

/*函数声明*/
inline void read_monitor_create();//创建readData监控线程
inline void add_solution(int cl, int tid);//加入solution
inline void recall_memory(int tid, int nl);
inline void initUV();//普通图转换成二分图并且初始化L,R,P,Q
inline void print_item(Type_Biclique &item);//打印一个极大二分团
inline void print_thread_end_time();//打印线程结束时间
inline void push_stack(int tid, int cl);
inline void distribute_sorted_work(int tid);//P有序版本，根据tid分配任务
inline int getCmdLineArgumentInt(const int argc, const char **argv, const char *string_ref);//从输入中获得整数
inline bool getCmdLineArgumentString(const int argc, const char **argv, const char *string_ref, char **string_retval);//从输入中获取字符串
inline bool checkCmdLineFlag(const int argc, const char **argv, const char *string_ref);////从输入中获得布尔值
inline int stringRemoveDelimiter(char delimiter, const char *string);//删除界定符，例如字符串"-flag"删除"-"
inline void print_help_info(); //打印帮助信息
inline void print_parameter(); //打印当前参数信息

vector<string> split(string &str, char ch);//字符串分割，读入文件转邻接表用
void int2str(const long &int_temp, string &string_temp);//整数转成字符串
void str2int(long &int_temp, const string &string_temp);//字符串转成整数
void checkArugument(int argc, char** argv);
void initial();
void readData();
void end();
void create_threads();// 创建线程
void init_stack(int tid);//初始化栈
void do_task(int tid);//执行任务
void cut_P(long &ci, int &cl, int tid);
void cut_p_q_biclique();//提前p-q-core剪枝，模拟k-core分解，循环左右遍历，查找度小于p和q的点，删除点和相关的边
bool find_U_R_Num();
bool find_V_L_Num();
void get_task(int src_thread_id, int &des_thread_id);//线程执行完成，开始从其他线程获取任务
void sort_new_P(vector<Type_v_sets> &PN, Type_P &new_P);
void get_small_PQ(int tid, int ci, int cl);
void operate_small_Q(bool &flag_maximal, int tid, int cl);
void operate_big_Q(bool &flag_maximal, int tid, int cl, int ci);
void operate_small_P(int tid, int cx, int cl, Type_v_sets &complement_new_L);
void operate_big_P(int tid, int cl, int ci, Type_v_sets &complement_new_L);
void * thread_compute(void * args);//线程执行的函数
void * thread_level_monitor(void * arg);//监视线程tid的level大小
void * thread_stack_monitor(void * arg);//监视线程tid的stack大小
void * thread_VP_monitor(void * arg);//监视线程tid的VP内存消耗情况
void * thread_read_monitor(void * arg);//监视读入数据的大小情况
bool compare_ov(long ov_a, long ov_b);
bool loadBalance(int src_thread_id);//负载均衡


void test() {
	vector<int> a;
	vector<int> b;
	a.push_back(1);
	a.push_back(2);
	a.push_back(3);
	b.push_back(4);
	b.push_back(5);
	b = a;
}

inline void print_parameter() {
	printf("---------------------------\n");
	printf("thread num:%d\n", THREAD_NUM);
	printf("数据名称:%s\n", str_dataName);
	printf("flag_p_q_biclique:%d\n", flag_p_q_biclique);
	printf("flag_is_sparse:%d\n", flag_is_sparse);
	printf("flag_buffer_interval_value:%d\n", flag_buffer_interval_value);
	printf("flag_stack_monitor:%d\n", flag_stack_monitor);
	printf("flag_level_monitor:%d\n", flag_level_monitor);
	printf("flag_read_monitor:%d\n", flag_read_monitor);
	printf("flag_VP_monitor:%d\n", flag_VP_monitor);
	printf("flag_loadbalance:%d\n", flag_loadbalance);
	printf("R_NUM:%ld\n", R_NUM);
	printf("L_NUM:%ld\n", L_NUM);
	printf("---------------------------\n");
}
inline void print_help_info() {
	printf("---------------------------\n");
	printf("参数'-dataName xxx'表示处理xxx文件\n");
	printf("参数'-sparse 1'表示稀疏二分图，否则不是，默认为1\n");
	printf("参数'-pqbiclique 1'表示要枚举极大p-q-biclique，否则不要，默认为0\n");
	printf("参数'-RNum 和 -LNum'表示p-q-biclique参数，数量默认为1\n");
	printf("参数'-buffer 1'表示要缓存中间结果，否则不要，默认为1\n");
	printf("---------监控标志位--------\n");
	printf("参数'-stackMonitor 1'表示要监控stack，否则不要，默认为1\n");
	printf("参数'-levelMonitor 1'表示要监控level，否则不要，默认为0\n");
	printf("参数'-readMonitor 1'表示要监控read，否则不要，默认为0\n");
	printf("参数'-VPMonitor 1'表示要监控VP，否则不要，默认为1\n");
	printf("---------注意注意--------\n");
	printf("输入-info输出当前参数信息\n");
	printf("暂时，thread num要在cpp文件中修改！！\n");
	printf("---------------------------\n");
}
void checkArugument(int argc, char** argv) {
	//替换文件名字符串中的指定字符串
	string str_repalce = "dataName";//即将被替换的string
	string default_dataName = "data00";
	int pos;
	if (checkCmdLineFlag(argc, (const char**)argv, "dataName")) {
		char* dataName;
		getCmdLineArgumentString(argc, (const char**)argv, "dataName", &dataName);
		str_dataName = dataName;
		if ((pos = filePath_data.find(str_repalce)) != string::npos) {
			filePath_data.replace(pos, str_repalce.size(), str_dataName);
			filePath_data_info.replace(pos, str_repalce.size(), str_dataName);
			output_file.replace(pos, str_repalce.size(), str_dataName);
			output_file_stack.replace(pos, str_repalce.size(), str_dataName);
			output_file_level.replace(pos, str_repalce.size(), str_dataName);
			output_file_VP.replace(pos, str_repalce.size(), str_dataName);
		}
	}
	else {
		//默认都是计算default_dataName
		printf("-dataName 没有设置，默认处理data00.txt\n");
		if ((pos = filePath_data.find(str_repalce)) != string::npos) {
			filePath_data.replace(pos, str_repalce.size(), default_dataName);
			filePath_data_info.replace(pos, str_repalce.size(), default_dataName);
			output_file.replace(pos, str_repalce.size(), default_dataName);
			output_file_stack.replace(pos, str_repalce.size(), default_dataName);
			output_file_level.replace(pos, str_repalce.size(), default_dataName);
			output_file_read.replace(pos, str_repalce.size(), default_dataName);
			output_file_VP.replace(pos, str_repalce.size(), default_dataName);
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "sparse")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "sparse");
		if (flag==1) {
			flag_is_sparse = true;
		}
		else {
			flag_is_sparse = false;
		}
	}

	if (checkCmdLineFlag(argc, (const char**)argv, "LB")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "LB");
		if (flag == 1) {
			flag_loadbalance = true;
		}
		else {
			flag_loadbalance = false;
		}
	}

	if (checkCmdLineFlag(argc, (const char**)argv, "pqbiclique")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "pqbiclique");
		if (flag == 1) {
			flag_p_q_biclique = true;
		}
		else {
			flag_p_q_biclique = false;
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "stackMonitor")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "stackMonitor");
		if (flag == 1) {
			flag_stack_monitor = true;
		}
		else {
			flag_stack_monitor = false;
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "readMonitor")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "readMonitor");
		if (flag == 1) {
			flag_read_monitor = true;
		}
		else {
			flag_read_monitor = false;
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "VPMonitor")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "VPMonitor");
		if (flag == 1) {
			flag_VP_monitor = true;
		}
		else {
			flag_VP_monitor = false;
		}
	}

	if (checkCmdLineFlag(argc, (const char**)argv, "levelMonitor")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "levelMonitor");
		if (flag == 1) {
			flag_level_monitor = true;
		}
		else {
			flag_level_monitor = false;
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "buffer")) {
		int flag = getCmdLineArgumentInt(argc, (const char**)argv, "buffer");
		if (flag == 1) {
			flag_buffer_interval_value = true;
		}
		else {
			flag_buffer_interval_value = false;
		}
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "RNum")) {
		int R_NUM = getCmdLineArgumentInt(argc, (const char**)argv, "RNum");
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "LNum")) {
		int L_NUM = getCmdLineArgumentInt(argc, (const char**)argv, "LNum");
	}
}
void generateRandomData() {
	//打开输出文件
	ofstream ofile_demo_data;//用于输出文件
	ofstream ofile_demo_data_info;//用于读取原始图数据文件
	const char * output_file = "E:/bigraph-data/new_data00.txt";//图数据
	const char * output_file_info = "E:/bigraph-data/new_data00-info.txt";//图点和边数信息
	ofile_demo_data.open(output_file, ios::out);
	ofile_demo_data_info.open(output_file_info, ios::out);
	long uNum = 7;
	long vNum = 4;
	vector<long> U_global;
	vector<long> V_global;
	int edge_num;
	int total_edge = 0;
	for (long i = 0; i < vNum + uNum; ++i) {
		if (i<uNum) {
			U_global.push_back(i);
		}
		else {
			V_global.push_back(i);
		}
	}
	srand(time(NULL));
	vector<long> sequence;
	for (int i = 1; i <= vNum; ++i) {
		sequence.push_back(i);
	}
	for (int u = 0; u < uNum; ++u) {
		random_shuffle(sequence.begin(), sequence.end());
		edge_num = sequence[0];
		total_edge += edge_num;
		random_shuffle(V_global.begin(), V_global.end());
		for (int i = 0; i < edge_num; i++)
		{
			ofile_demo_data << u << " " << V_global[i] << endl;
		}
	}
	ofile_demo_data_info << vNum + uNum << " " << total_edge << endl;
	ofile_demo_data_info << uNum << endl;

	ofile_demo_data.close();// 关闭输出文件
	ofile_demo_data_info.close();// 关闭输出文件
}

int main(int argc, char** argv) {
	//while (!flag_success_demo_1|| !flag_success_demo_2|| !flag_success_demo_3||!flag_success_demo_4||!flag_success_demo_5) {
	//generateRandomData();

	if (checkCmdLineFlag(argc, (const char**)argv, "help")) {
		print_help_info();
		return 0;
	}
	if (checkCmdLineFlag(argc, (const char**)argv, "info")) {
		checkArugument(argc, argv);
		print_parameter();
		return 0;
	}

	checkArugument(argc, argv);
	initial();
	create_threads();
	end();
//}
	return 0;
}

void cut_p_q_biclique() {
	bool flag ;
	while (true) {
		flag = false;
		if (find_U_R_Num()) {
			flag = true;
		}
		if (find_V_L_Num()) {
			flag = true;
		}
		if (flag == false) {
			break;
		}
	}
	//处理degree
	long index_new_graph = 0;
	p_q_graph = Type_Graph(p_q_node_num);

	
	for (long u = 0; u< node_num; ++u) {
		if (degree[u] != 0) {
			//建立映射
			if (p_q_U_num>p_q_V_num) {
				p_q_graph_map[u] = index_new_graph;
				p_q_graph_map_reverse[index_new_graph] = u;
			}
			else {
				p_q_graph_map[u] = p_q_node_num - 1 -index_new_graph;
				p_q_graph_map_reverse[p_q_node_num - 1 - index_new_graph] = u;
			}
			
			set<long>::iterator iter = graph[u].begin();//遍历它的邻居
			while (iter != graph[u].end()) {
				long neibor = *iter;
				if (degree[neibor] != 0) {
					if (p_q_U_num>p_q_V_num) {
						p_q_graph[index_new_graph].insert(neibor);
					}
					else {
						p_q_graph[p_q_node_num - 1 - index_new_graph].insert(neibor);
					}
					
				}
				++iter;
			}
			index_new_graph++;
		}
	}

	//处理p_q_graph的节点序号，注意生成新旧之间的映射
	for (long u = 0; u < p_q_node_num; ++u) {
		set<long> temp;
		Type_v_sets::iterator iter = p_q_graph[u].begin();
		while (iter!= p_q_graph[u].end()) {
			temp.insert(p_q_graph_map[*iter]);
			++iter;
		}
		p_q_graph[u].clear();
		p_q_graph[u] = temp;
	}
	
}
bool find_U_R_Num() {
	//遍历U的点
	long v = 0;
	int offset = 1;
	bool flag = false;
	while (v<U_num) {
		if (degree[v]<R_NUM && degree[v]>0) {
			flag = true;
			degree[v] = 0;
			p_q_node_num--;
			p_q_U_num--;
			//遍历它的邻居
			set<long>::iterator iter = graph[v].begin();
			while (iter != graph[v].end()) {
				long neibor = *iter;
				if (degree[neibor] > 0) {
					degree[neibor]--;//邻居的度减一
				}
				++iter;
			}
		}
		v += offset;
	}
	return flag;
}
bool find_V_L_Num() {
	//遍历右边的点
	long v = U_num;
	int offset = 1;
	bool flag = false;
	while (v<node_num) {
		if (degree[v]<L_NUM && degree[v]>0) {
			flag = true;
			degree[v] = 0;
			p_q_node_num--;
			p_q_V_num--;
			//遍历它的邻居
			set<long>::iterator iter = graph[v].begin();
			while (iter != graph[v].end()) {
				long neibor = *iter;
				if (degree[neibor] > 0) {
					//邻居的度减一
					degree[neibor]--;
				}
				++iter;
			}
		}
		v += offset;
	}
	return flag;
}
void * thread_level_monitor(void * arg) {
	timespec outtime;
	pthread_mutex_lock(&level_monitor_mutex);
	while (flag_level_monitor) {
		outtime.tv_sec = time(NULL) + 3;
		pthread_cond_timedwait(&level_monitor_cond, &level_monitor_mutex, &outtime);
		for (int i = 0; i < THREAD_NUM; ++i) {
			ofile_level << P[i].size() << ":" << (P[i].back()).size() << "\t\t";
		}
		ofile_level << endl;
	}
	pthread_mutex_unlock(&level_monitor_mutex);
	return NULL;
}

void * thread_read_monitor(void * arg) {
	timespec outtime;
	pthread_mutex_lock(&read_monitor_mutex);
	while (flag_read_monitor) {
		outtime.tv_sec = time(NULL) + 3;
		pthread_cond_timedwait(&read_monitor_cond, &read_monitor_mutex, &outtime);
		ofile_read << read_line_number << endl;
	}
	pthread_mutex_unlock(&read_monitor_mutex);
	pthread_mutex_destroy(&read_monitor_mutex);
	pthread_cond_destroy(&read_monitor_cond);
	return NULL;
}

void * thread_stack_monitor(void * arg) {
	timespec outtime;
	pthread_mutex_lock(&stack_monitor_mutex);
	while (flag_stack_monitor) {
		outtime.tv_sec = time(NULL) + 3;
		pthread_cond_timedwait(&stack_monitor_cond, &stack_monitor_mutex, &outtime);
		for (int i = 0; i < THREAD_NUM; ++i) {
			ofile_stack << task_stack[i].size() << "\t\t";
		}
		ofile_stack << endl;
	}
	pthread_mutex_unlock(&stack_monitor_mutex);
	return NULL;
}

void * thread_VP_monitor(void * arg) {
	timespec outtime;
	pthread_mutex_lock(&VP_monitor_mutex);
	while (flag_stack_monitor) {
		outtime.tv_sec = time(NULL) + 5;
		pthread_cond_timedwait(&VP_monitor_cond, &VP_monitor_mutex, &outtime);
		//输出VP所有数组的size，即表示内存占用
		for (int i = 0; i < THREAD_NUM; ++i) {
			int total_size = 0;
			total_size += P[i].size()*sizeof(pair<long, bool>);
			for (int j = 0; j < PN[i].size(); ++j) {
				total_size += PN[i][j].size() * sizeof(long);
			}
			total_size += L[i].size() * sizeof(long);
			total_size += R[i].size() * sizeof(long);
			total_size += Q[i].size() * sizeof(long);
			ofile_VP << total_size << "\t\t";
		}
		ofile_VP << endl;
	}
	pthread_mutex_unlock(&VP_monitor_mutex);
	return NULL;
}

bool loadBalance(int src_thread_id) {
	bool task_got = false;
	int des_thread_id;
	get_task(src_thread_id, des_thread_id);
	if (des_thread_id != -1) {//表示其他线程的任务都大于阀值,开始get任务
		task_got = true;
		pthread_mutex_lock(&mutex_thread[des_thread_id]);//给线程加锁

														 //非常重要，考虑到同层的C的先后发现问题，如果同层的栈数据（为探测节点）被移到了other thread中，可能出现本地thread的节点在C中，没办法被发现，会重复探索出相同结果
														 //取出des_thread_id一半任务,连同*_info[des_thread_id]
														 //比如栈中有第1,2,4层，取完任务之后只有第1层，就应该pop 3次
		long size_task = task_stack[des_thread_id].size();
		stack<TASK> stack_temp;
		//这里不一定有栈顶元素,因为防止在判断des_thread_id之后到上锁的时间内可能栈为空
		if (!task_stack[des_thread_id].empty()) {//栈为空，表示选中的已经足够小了，可以继续执行负载均衡了。
			int cl = task_stack[des_thread_id].top().level;
			for (long i = 0; i < size_task / 2; ++i) {
				TASK task = task_stack[des_thread_id].top();
				task_stack[des_thread_id].pop();
				stack_temp.push(task);
			}
			//考虑到顺序问题，以后可以改成双端队列
			while (!stack_temp.empty()) {
				TASK task = stack_temp.top();
				stack_temp.pop();
				task_stack[src_thread_id].push(task);
			}

			P[src_thread_id] = P[des_thread_id];
			if (flag_buffer_interval_value) {
				PN[src_thread_id] = PN[des_thread_id];
			}
			
			Q[src_thread_id] = Q[des_thread_id];
			R[src_thread_id] = R[des_thread_id];
			L[src_thread_id] = L[des_thread_id];
			//RE[src_thread_id] = RE[des_thread_id];

			//修改响应端端的全局数组信息,被取走的层级的信息全部删掉

			TASK nt = task_stack[des_thread_id].top();//查看取完之后当前top层级
			recall_memory(des_thread_id, nt.level);
		}
		pthread_mutex_unlock(&mutex_thread[des_thread_id]);		//给线程释放锁
	}
	return task_got;
}
void get_task(int src_thread_id, int &des_thread_id) {
	//遍历其他所有线程，从最多任务的线程get task
	des_thread_id = -1;//默认目的线程id
	int max = TASK_THRESHOLD_VALUE;//设置一个阀值
	for (int i = 0; i < THREAD_NUM; ++i) {
		if (src_thread_id != i) {
			int task_num = task_stack[i].size();
			if (max < task_num) {
				max = task_num;
				des_thread_id = i;
			}
		}
	}
}
void InsertionSort_ordered_sets(Type_ordered_sets &v, Type_ordered_sets &CNSize)
{
	for (int j = 1; j<v.size(); j++)
	{
		int key = CNSize[j];
		int key_v = v[j];
		int i = j - 1;
		while (i >= 0 && CNSize[i]>key)
		{
			CNSize[i + 1] = CNSize[i];
			v[i + 1] = v[i];
			i--;
		}
		CNSize[i + 1] = key;
		v[i + 1] = key_v;
	}
}

inline void distribute_sorted_work(int tid) {
	long v = 0;//遍历所有节点
	while (v != V_num) {
		if (v%THREAD_NUM == tid) {//分配节点规则
			task_init[tid].push_back(v);//只存下标
		}
		++v;
	}
}

void init_stack(int tid) {
	P[tid].push_back(initial_P);
	L[tid].push_back(U);
	if (flag_buffer_interval_value) {
		PN[tid].push_back(vector<Type_v_sets>());
	}
	Q[tid].push_back(Type_ordered_sets());
	R[tid].push_back(Type_v_sets());

	distribute_sorted_work(tid);//分配工作

	long i = task_init[tid].size() - 1;
	while (i >= 0) {
		TASK task;
		task.level = 0;
		task.x_index = task_init[tid][i];//节点号即序号
		task_stack[tid].push(task);
		--i;
	}

}

inline void add_solution(int cl, int tid) {
	Type_v_sets RAndRE;
	for (int index = 0; index <= cl; ++index) {//获取R[tid]、RE[tid]每一行的所有元素
		Type_v_sets::iterator iter = R[tid][index].begin();
		while (iter != R[tid][index].end()) {
			RAndRE.insert(*iter);
			++iter;
		}
	}
	if (flag_p_q_biclique) {
		if (RAndRE.size() < R_NUM || L[tid][cl + 1].size() < L_NUM) {
			return;
		}
		else {
			solution[tid].insert(Type_Biclique(tid, pair<Type_v_sets, Type_v_sets>(L[tid][cl + 1], RAndRE)));
		}
	}
	else {
		solution[tid].insert(Type_Biclique(tid, pair<Type_v_sets, Type_v_sets>(L[tid][cl + 1], RAndRE)));
	}
}

inline void recall_memory(int tid,int nl) {
	while (P[tid].size()>nl + 2) {
		P[tid].pop_back();
		if (flag_buffer_interval_value) {
			PN[tid].pop_back();
		}
		L[tid].pop_back();
		Q[tid].pop_back();
		R[tid].pop_back();
	}
	//R[tid][nl] = Type_v_sets();
	if (flag_buffer_interval_value) {
		//PN[tid][nl] = vector<Type_v_sets>();
	}
}
inline void push_stack(int tid, int cl) {

	long i = P[tid][cl + 1].size() - 1;
	while (i >= 0) {
		TASK task;
		task.level = cl + 1;
		task.x_index = i;
		task_stack[tid].push(task);
		--i;
	}
}

void get_small_PQ(int tid, int ci, int cl) {
	new_L_neibors[tid] = Type_v_sets();
	SQ[tid] = Type_v_sets();
	SP[tid] = Type_ordered_sets();
	SPP[tid] = Type_ordered_sets();
	Type_v_sets::iterator iter_new_L = L[tid][cl + 1].begin();
	while (iter_new_L != L[tid][cl + 1].end()) {
		set_union(new_L_neibors[tid].begin(), new_L_neibors[tid].end(), graph[*iter_new_L].begin(), graph[*iter_new_L].end(), inserter(new_L_neibors[tid], new_L_neibors[tid].begin()));
		++iter_new_L;
	}
	Type_v_sets::iterator iter_new_L_neibors = new_L_neibors[tid].begin();
	while (iter_new_L_neibors != new_L_neibors[tid].end()) {
		long v = (*iter_new_L_neibors);
		long index = OVP[v - U_num];
		if (index < ci) {
			SQ[tid].insert(v);
		}
		else {
			SP[tid].push_back(v);
			SPP[tid].push_back(index);

		}
		++iter_new_L_neibors;
	}
	//根据index的大小排序，小的在前面
	InsertionSort_ordered_sets(SP[tid], SPP[tid]);
}
void operate_small_P(int tid, int cx, int cl, Type_v_sets &complement_new_L) {
	Type_ordered_sets::iterator iter_P1 = SP[tid].begin();
	int index = 0;
	while (iter_P1 != SP[tid].end()) {
		int v = (*iter_P1);
		if (v != cx) {
			Type_v_sets Nv;// Get the neighbors of v in L
			set_intersection(L[tid][cl + 1].begin(), L[tid][cl + 1].end(), graph[v].begin(), graph[v].end(), inserter(Nv, Nv.begin()));
			if (flag_p_q_biclique && Nv.size()<L_NUM) {
				++iter_P1;
				++index;
				continue;
			}
			if (Nv.size() == L[tid][cl + 1].size()) {
				R[tid][cl].insert(v);
				Type_v_sets S;
				set_intersection(complement_new_L.begin(), complement_new_L.end(), graph[v].begin(), graph[v].end(), inserter(S, S.begin()));
				if (S.size() == 0) {
					P[tid][cl][SPP[tid][index]].second = false;
					Q[tid][cl].push_back(v);
				}
			}
			else if (Nv.size()>0) {
				P[tid][cl + 1].push_back(pair<int, bool>(v, true));
				if (flag_buffer_interval_value) {
					PN[tid][cl].push_back(Type_v_sets(Nv));
				}
			}
		}
		++iter_P1;
		++index;
	}
}
void operate_big_P(int tid, int cl, int ci, Type_v_sets &complement_new_L) {
	for (int i = ci + 1; i < P[tid][cl].size(); ++i) {
		if (P[tid][cl][i].second) {
			int v = P[tid][cl][i].first;
			Type_v_sets Nv;// Get the neighbors of v in L
			set_intersection(L[tid][cl + 1].begin(), L[tid][cl + 1].end(), graph[v].begin(), graph[v].end(), inserter(Nv, Nv.begin()));
			if (flag_p_q_biclique && Nv.size()<L_NUM) {
				continue;
			}
			if (Nv.size() == L[tid][cl + 1].size()) {
				R[tid][cl].insert(v);
				flag_success_demo_2 = true;
				Type_v_sets S;
				set_intersection(complement_new_L.begin(), complement_new_L.end(), graph[v].begin(), graph[v].end(), inserter(S, S.begin()));
				if (S.size() == 0) {
					flag_success_demo_3 = true;
					P[tid][cl][i].second = false;
					Q[tid][cl].push_back(v);//当前层Q增加v
				}
			}
			else if (Nv.size()>0) {
				P[tid][cl + 1].push_back(pair<int, bool>(v, true));
				if (flag_buffer_interval_value) {
					PN[tid][cl].push_back(Type_v_sets(Nv));
				}
			}
		}
		else {
			int kk = 1;
		}
	}
}
void cut_P(long &ci, int &cl, int tid) {

	pair<long, bool> p = P[tid][cl][ci];
	//在这里控制分割
	while (!p.second) {
	//while (!p.second) {
		if (!task_stack[tid].empty()) {
			flag_success_demo_4 = true;
			TASK top_task = task_stack[tid].top();
			task_stack[tid].pop();
			nl[tid] = top_task.level;
			ni[tid] = top_task.x_index;
			recall_memory(tid, nl[tid]);//回收内存
			p = P[tid][nl[tid]][ni[tid]];
			cl = nl[tid];
			ci = ni[tid];
			
		}
		else {
			break;
		}
	}
}
void operate_small_Q(bool &flag_maximal, int tid, int cl) {
	Type_v_sets::iterator iter_Q = SQ[tid].begin();
	while (iter_Q != SQ[tid].end()) {
		int v = *iter_Q;
		Type_v_sets Nv;
		set_intersection(L[tid][cl + 1].begin(), L[tid][cl + 1].end(), graph[v].begin(), graph[v].end(), inserter(Nv, Nv.begin()));
		if (Nv.size() == L[tid][cl + 1].size()) {
			flag_maximal = false;
			if (!task_stack[tid].empty()) {
				nt[tid] = task_stack[tid].top();
				nl[tid] = nt[tid].level;
				recall_memory(tid, nl[tid]);
			}
			break;
		}
		else if (Nv.size()>0) {
			Q[tid][cl + 1].push_back(v);
		}
		++iter_Q;
	}
}
void operate_big_Q(bool &flag_maximal, int tid, int cl, int ci) {
	//先遍历Q[tid][cl]
	Type_ordered_sets::iterator iter_Q = Q[tid][cl].begin();
	while (iter_Q != Q[tid][cl].end()) {
		long v = *iter_Q;
		Type_v_sets Nv;
		set_intersection(L[tid][cl + 1].begin(), L[tid][cl + 1].end(), graph[v].begin(), graph[v].end(), inserter(Nv, Nv.begin()));
		if (Nv.size() == L[tid][cl + 1].size()) {
			flag_maximal = false;
			flag_success_demo_5 = true;
			if (!task_stack[tid].empty()) {
				nl[tid] = task_stack[tid].top().level;
				recall_memory(tid, nl[tid]);
			}
			return;
		}
		else if (Nv.size()>0) {
			Q[tid][cl + 1].push_back(v);
			flag_success_demo_1 = true;
		}
		++iter_Q;
	}
	//再遍历P[tid][cl][0],...,P[tid][cl][ci-1]
	for (int i = 0; i < ci; ++i) {
		long v = P[tid][cl][i].first;
		Type_v_sets Nv;
		set_intersection(L[tid][cl + 1].begin(), L[tid][cl + 1].end(), graph[v].begin(), graph[v].end(), inserter(Nv, Nv.begin()));
		if (Nv.size() == L[tid][cl + 1].size()) {
			flag_maximal = false;
			if (!task_stack[tid].empty()) {
				nl[tid] = task_stack[tid].top().level;
				recall_memory(tid, nl[tid]);
			}
			return;
		}
		else if (Nv.size()>0) {
			Q[tid][cl + 1].push_back(v);
			flag_success_demo_1 = true;
		}
	}
}

void do_task(int tid) {
	while (!task_stack[tid].empty()) {
		if (flag_loadbalance) {
			pthread_mutex_lock(&mutex_thread[tid]);//给线程加锁
		}

		Type_v_sets complement_new_L;
		flag_maximal[tid] = true;

		TASK ct = task_stack[tid].top();
		task_stack[tid].pop();
		int cl = ct.level;
		long ci = ct.x_index;
		cut_P(ci, cl, tid);//找出cf为false的节点剪枝
		
		if (cl + 1 == P[tid].size()) {
			P[tid].push_back(Type_P());
			if (flag_buffer_interval_value) {
				PN[tid].push_back(vector<Type_v_sets>());
			}
			L[tid].push_back(Type_v_sets());
			Q[tid].push_back(Type_ordered_sets());
			R[tid].push_back(Type_v_sets());
		}

		/*else {
			P[tid][cl + 1] = Type_P();
			if (flag_buffer_interval_value) {
				PN[tid][cl + 1] = vector<Type_v_sets>();
			}
			L[tid][cl + 1] = Type_v_sets();
			Q[tid][cl + 1] = Type_ordered_sets();
			R[tid][cl + 1] = Type_v_sets();
		}*/

		P[tid][cl + 1] = Type_P();
		L[tid][cl + 1] = Type_v_sets();
		Q[tid][cl + 1] = Type_ordered_sets();

		long cx = P[tid][cl][ci].first;

		//设置多PC，解除注释即可
		/*int thresh = (2 * node_num - V_num) / PC_number;
		if (cl==0&&cx < thresh) {
			if (flag_loadbalance) {
				pthread_mutex_unlock(&mutex_thread[tid]);
			}
			continue;
		}*/
		PN[tid][cl] = vector<Type_v_sets>();
		R[tid][cl] = set<long>{cx};//当前节点加入该层R
		//R[tid][cl].insert(cx);
		if (flag_p_q_biclique) {
			if (graph[cx].size()<L_NUM) {//不用向下探索了,剪枝
				if (!task_stack[tid].empty()) {
					nl[tid] = task_stack[tid].top().level;
					recall_memory(tid, nl[tid]);
					if (flag_loadbalance) {
						pthread_mutex_unlock(&mutex_thread[tid]);
					}
					continue;
				}
				else {
					if (flag_loadbalance) {
						pthread_mutex_unlock(&mutex_thread[tid]);
					}
					break;
				}
			}
		}
		if (cl == 0) {/*通过上一层的PN获得new_L，第一层则直接是全局邻居，减少集合交运算*/
			L[tid][cl + 1] = graph[cx];
		}
		else {
			if (flag_buffer_interval_value) {
				L[tid][cl + 1] = PN[tid][cl - 1][ci];//上一层
			}
			else {
				//自己计算L[tid][cl + 1]
				set_intersection(L[tid][cl].begin(), L[tid][cl].end(), graph[cx].begin(), graph[cx].end(), inserter(L[tid][cl + 1], L[tid][cl + 1].begin()));
			}
		}
		set_difference(L[tid][cl].begin(), L[tid][cl].end(), L[tid][cl + 1].begin(), L[tid][cl + 1].end(), inserter(complement_new_L, complement_new_L.begin()));

		if (flag_is_sparse&&cl == 0) {
			get_small_PQ(tid, ci, cl);
			operate_small_Q(flag_maximal[tid], tid, cl);
		}
		else {
			operate_big_Q(flag_maximal[tid], tid, cl, ci);
		}


		if (flag_maximal[tid]) {
			if (flag_is_sparse&&cl == 0) {//第一层针对稀疏图剪枝,生成new_P
				operate_small_P(tid, cx, cl, complement_new_L);
			}
			else {
				operate_big_P(tid, cl, ci, complement_new_L);
			}

			add_solution(cl, tid);
			if (P[tid][cl + 1].empty()) {

				if (task_stack[tid].empty()) {
					if (flag_loadbalance) {
						pthread_mutex_unlock(&mutex_thread[tid]);
					}
					break;//如果是最后一个元素，且没有孩子节点，退出整个栈循环
				}
				else {
					nl[tid] = task_stack[tid].top().level;
					recall_memory(tid, nl[tid]);
				}
			}
			else {
				if (flag_buffer_interval_value) {
					sort_new_P(PN[tid][cl], P[tid][cl + 1]);
				}
				push_stack(tid, cl);
			}

		}

		if (flag_loadbalance) {
			pthread_mutex_unlock(&mutex_thread[tid]);
		}

	}//while

	 //线程执行完，获取其他线程的任务，如果获得了任务，则继续执行do_task(tid)
	if (flag_loadbalance && loadBalance(tid)) {
		do_task(tid);
	}
	//记录结束时间
	thread_stop_time[tid] = time(NULL);
}
void * thread_compute(void * args) {
	int tid = *(int*)args;//获得参数信息，用于分配节点
	init_stack(tid);
	do_task(tid);
	return NULL;
}
inline void read_monitor_create() {
	pthread_mutex_init(&read_monitor_mutex, NULL);
	pthread_cond_init(&read_monitor_cond, NULL);
	pthread_create(&read_monitor_thread, NULL, thread_read_monitor, NULL);
}

void create_threads() {
	pthread_t threads[THREAD_NUM];
	int i = 0;
	pthread_attr_t attr;
	/*初始化互斥锁*/
	for (int j = 0; j < THREAD_NUM; ++j) {
		if (flag_loadbalance) {
			pthread_mutex_init(&mutex_thread[j], NULL);
		}
	}
	pthread_mutex_init(&stack_monitor_mutex, NULL);
	pthread_mutex_init(&level_monitor_mutex, NULL);
	pthread_mutex_init(&VP_monitor_mutex, NULL);

	pthread_cond_init(&stack_monitor_cond, NULL);
	pthread_cond_init(&level_monitor_cond, NULL);
	pthread_cond_init(&VP_monitor_cond, NULL);

	pthread_create(&stack_monitor_thread, NULL, thread_stack_monitor, NULL);
	pthread_create(&VP_monitor_thread, NULL, thread_VP_monitor, NULL);
	pthread_create(&level_monitor_thread, NULL, thread_level_monitor, NULL);


	while (i<THREAD_NUM) {
		pthread_attr_init(&attr);//初始化属性
		pthread_create(&threads[i], &attr, thread_compute, (void*)&thread_ids[i]);//创建线程
		++i;
	}

	for (int j = 0; j < THREAD_NUM; j++) {
		pthread_join(threads[j], NULL);
	}

	pthread_mutex_lock(&stack_monitor_mutex);
	flag_stack_monitor = false;
	pthread_cond_signal(&stack_monitor_cond);
	pthread_mutex_unlock(&stack_monitor_mutex);

	pthread_mutex_lock(&VP_monitor_mutex);
	flag_VP_monitor = false;
	pthread_cond_signal(&VP_monitor_cond);
	pthread_mutex_unlock(&VP_monitor_mutex);

	pthread_mutex_lock(&level_monitor_mutex);
	flag_level_monitor = false;
	pthread_cond_signal(&level_monitor_cond);
	pthread_mutex_unlock(&level_monitor_mutex);

	pthread_join(stack_monitor_thread, NULL);
	pthread_join(level_monitor_thread, NULL);
	pthread_join(VP_monitor_thread, NULL);

	pthread_attr_destroy(&attr);//注销线程属性
	if (flag_loadbalance) {
		for (int j = 0; j < THREAD_NUM; ++j) {
			pthread_mutex_destroy(&mutex_thread[j]);
		}
	}

	pthread_mutex_destroy(&stack_monitor_mutex);
	pthread_mutex_destroy(&level_monitor_mutex);
	pthread_mutex_destroy(&VP_monitor_mutex);

	pthread_cond_destroy(&stack_monitor_cond);
	pthread_cond_destroy(&level_monitor_cond);
	pthread_cond_destroy(&VP_monitor_cond);
}
vector<string> split(string &str, char ch)
{
	vector<string> res;//存储分割以后的结果  
	str += ch;//加入一个分割字符，方便操作  
	int start = 0;//分割的起始位置  
	int last = str.find(ch);//分割的终止位置  
	while (last < str.size())//找到最后一个分隔符，终止条件  
	{
		if (start != last) {//如果字符串不为空，则添加到结果中  
			res.push_back(str.substr(start, last - start));
		}
		start = last + 1;//起始位置更新  
		last = str.find(ch, start);//终止位置更新  
	}
	return res;
}
void int2str(const long &int_temp, string &string_temp)
{
	stringstream stream;
	stream << int_temp;
	string_temp = stream.str();
}
void str2int(long &int_temp, const string &string_temp)
{
	stringstream stream(string_temp);
	stream >> int_temp;
}
void readData() {
	ofile << "thread number: " << THREAD_NUM << endl;
	time_t start_read = time(NULL);
	ofile << "start read time: " << start_read << endl;
	ifile.open(filePath_data.c_str(), ios::in);
	ifile_info.open(filePath_data_info.c_str(), ios::in);

	if (ifile.is_open() && ifile_info.is_open()) {

		string data_info;//图数据的节点数量和边数量的一行，例如123 123
		getline(ifile_info, data_info);
		vector<string> str_split;
		str_split = split(data_info, ' ');// 拆分数据，例如1 2（中间空格）
		str2int(node_num, str_split[0]);//节点数
		str2int(edge_num, str_split[1]);//边数
		if (edge_num > (node_num / 2)) { //判断是否为稀疏图
			//flag_is_sparse = false;
		}
		getline(ifile_info, data_info);
		str_split = split(data_info, ' ');
		str2int(U_num, str_split[0]);
		V_num = node_num - U_num;
		p_q_U_num = U_num;
		p_q_V_num = V_num;
		p_q_node_num = node_num;
		degree = vector<long>(node_num);
		graph = Type_Graph(node_num);
		string data;//原始图数据的一行
		while (getline(ifile, data))
		{
			++read_line_number;
			if (data.empty()) {
				continue;
			}
			str_split = split(data, ' ');
			long v, neibor;
			str2int(v, str_split[0]);
			str2int(neibor, str_split[1]);
			graph[v].insert(neibor);
			degree[v]++;
			graph[neibor].insert(v);
			degree[neibor]++;
		}

		ifile.close();
		ifile_info.close();
		flag_read_monitor = false;
	}
	time_t end_read = time(NULL);
	ofile << "total read time: " << end_read - start_read << endl;
	start_time = time(NULL);//读完数据开始计时
	ofile << "start time: " << start_time << endl;
}
void initial() {
	if (flag_read_monitor) {
		ofile_read.open(output_file_read.c_str(), ios::out);//打开read监控输出文件
		ofile_read << "data-read" << endl;
		read_monitor_create();
	}
	if (flag_stack_monitor) {
		ofile_stack.open(output_file_stack.c_str(), ios::out);//打开stack监控输出文件
	}
	if (flag_VP_monitor) {
		ofile_VP.open(output_file_VP.c_str(), ios::out);//打开stack监控输出文件
	}
	if (flag_level_monitor) {
		ofile_level.open(output_file_level.c_str(), ios::out);//打开level监控输出文件
	}

	ofile.open(output_file.c_str(), ios::out);//打开输出文件

	readData();//读数据
	
	initUV();//初始化UV
	
	for (int i = 0; i < THREAD_NUM; ++i) {
		thread_ids[i] = i;//初始化threads_id信息
		if (flag_stack_monitor) {
			ofile_stack << "T" << i << "\t\t";
		}
		if (flag_level_monitor) {
			ofile_level << "T" << i << "\t\t";
		}
	}
	if (flag_stack_monitor) {
		ofile_stack << endl;
	}
	if (flag_level_monitor) {
		ofile_level << endl;
	}
	
}
inline void print_item(Type_Biclique &item) {
	int tid = item.first;
	Type_v_sets::iterator iter_L = item.second.first.begin();
	Type_v_sets::iterator iter_R = item.second.second.begin();
	ofile << "{< ";

	if (flag_p_q_biclique) {
		if (p_q_U_num>p_q_V_num) {
			for (iter_R; iter_R != item.second.second.end(); ++iter_R) {
				ofile << p_q_graph_map_reverse[*iter_R] << " ";
			}
		}
		else {
			for (iter_L; iter_L != item.second.first.end(); ++iter_L) {
				ofile << p_q_graph_map_reverse[*iter_L] << " ";
			}
		}
		
	}
	else {
		for (iter_R; iter_R != item.second.second.end(); ++iter_R) {
			ofile << *iter_R << " ";
		}
	}
	ofile << ">;< ";

	if (flag_p_q_biclique) {
		if (p_q_U_num>p_q_V_num) {
			for (iter_L; iter_L != item.second.first.end(); ++iter_L) {
				ofile << p_q_graph_map_reverse[*iter_L] << " ";
			}
		}
		else {
			for (iter_R; iter_R != item.second.second.end(); ++iter_R) {
				ofile << p_q_graph_map_reverse[*iter_R] << " ";
			}
		}
	}
	else {
		for (iter_L; iter_L != item.second.first.end(); ++iter_L) {
			ofile << *iter_L << " ";
		}
	}
	

	ofile << ">}" << " tid:" << tid << endl;
}
inline void print_thread_end_time() {
	for (int i = 0; i < THREAD_NUM; ++i) {
		ofile << " thread " << i << " total time: " << thread_stop_time[i] - start_time << endl;
	}
}
void end() {
	if (flag_level_monitor) {
		ofile_level.close();
	}
	if (flag_stack_monitor) {
		ofile_stack.close();
	}
	if (flag_read_monitor) {
		ofile_read.close();
	}
	end_time = time(NULL);
	ofile << "total time: " << end_time - start_time << endl;
	ofile << "total time: " << end_time - start_time << endl;
	for (int tid = 0; tid < THREAD_NUM; ++tid) {
		set<Type_Biclique>::iterator solution_iter = solution[tid].begin();
		for (solution_iter; solution_iter != solution[tid].end(); ++solution_iter) {
			Type_Biclique item = (*solution_iter);
			print_item(item);
		}
	}
	print_thread_end_time();
	ofile.close();// 关闭输出文件
}

bool compare_ov(long ov_a, long ov_b) {
	long a_L_neibor_num = 0;
	long b_L_neibor_num = 0;
	if (graph[ov_a].size() == graph[ov_b].size()) {
		Type_v_sets::iterator iter_neibor_a = graph[ov_a].begin();
		while (iter_neibor_a != graph[ov_a].end()) {
			a_L_neibor_num += graph[*iter_neibor_a].size();
			++iter_neibor_a;
		}
		Type_v_sets::iterator iter_neibor_b = graph[ov_b].begin();
		while (iter_neibor_b != graph[ov_b].end()) {
			b_L_neibor_num += graph[*iter_neibor_b].size();
			++iter_neibor_b;
		}
		return a_L_neibor_num < b_L_neibor_num;
	}
	return graph[ov_a].size() < graph[ov_b].size();
}

inline void initUV() {
	if (flag_p_q_biclique) {
		cut_p_q_biclique();
		//执行结束之后判断谁大谁小，小的给V
		if (p_q_U_num>p_q_V_num) {
			V_num = p_q_V_num;
			U_num = p_q_U_num;
		}
		else {
			V_num = p_q_U_num;
			U_num = p_q_V_num;
			//交换下R_NUM和L_Num条件
			int temp = R_NUM;
			R_NUM = L_NUM;
			L_NUM = temp;
		}
		node_num = p_q_node_num;
		graph.clear();
		graph = p_q_graph;
	}
	OV = new long[V_num];
	OVP = new long[V_num];
	for (long i = 0; i < U_num; ++i) {
		U.insert(i);
	}
	for (long j = 0; j < V_num; ++j) {
		OV[j] = j + U_num;
		if (flag_is_sparse) {
			OVP[j] = j + U_num;
		}
	}

	time_t  start_sort = time(NULL);
	ofile << "sort OV start time: " << start_sort << endl;
	sort(OV, OV + V_num, compare_ov);
	ofile << "sort OV total time: " << time(NULL) - start_sort << endl;
	
	initial_P = Type_P(V_num);
	
	for (int i = 0; i < V_num; ++i) {
		initial_P[i] = pair<long, bool>(OV[i], true);
		if (flag_is_sparse) {
			OVP[OV[i] - U_num] = i;
		}
	}
	delete[] OV;
}

void sort_new_P(vector<Type_v_sets> &PN, Type_P &new_P){
	for (long j = 1; j<PN.size(); j++)
	{
		Type_v_sets key = PN[j];
		pair<long, bool> key_v = new_P[j];
		long i = j - 1;
		while (i >= 0 && PN[i].size()>key.size())
		{
			new_P[i + 1] = new_P[i];
			PN[i + 1] = PN[i];
			i--;
		}
		PN[i + 1] = key;
		new_P[i + 1] = key_v;
	}
}
inline int stringRemoveDelimiter(char delimiter, const char *string)
{
	int string_start = 0;

	while (string[string_start] == delimiter)
	{
		string_start++;
	}

	if (string_start >= (int)strlen(string) - 1)
	{
		return 0;
	}

	return string_start;
}

inline bool checkCmdLineFlag(const int argc, const char **argv, const char *string_ref)
{
	bool bFound = false;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];

			const char *equal_pos = strchr(string_argv, '=');
			int argv_length = (int)(equal_pos == 0 ? strlen(string_argv) : equal_pos - string_argv);

			int length = (int)strlen(string_ref);

			if (length == argv_length && !strncmp(string_argv, string_ref, length))
			{
				bFound = true;
				continue;
			}
		}
	}

	return bFound;
}

inline bool getCmdLineArgumentString(const int argc, const char **argv,
	const char *string_ref, char **string_retval)
{
	bool bFound = false;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			char *string_argv = (char *)&argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!strncmp(string_argv, string_ref, length))
			{
				*string_retval = &string_argv[length + 1];
				bFound = true;
				continue;
			}
		}
	}

	if (!bFound)
	{
		*string_retval = NULL;
	}

	return bFound;
}

inline int getCmdLineArgumentInt(const int argc, const char **argv, const char *string_ref)
{
	bool bFound = false;
	int value = -1;

	if (argc >= 1)
	{
		for (int i = 1; i < argc; i++)
		{
			int string_start = stringRemoveDelimiter('-', argv[i]);
			const char *string_argv = &argv[i][string_start];
			int length = (int)strlen(string_ref);

			if (!strncmp(string_argv, string_ref, length))
			{
				if (length + 1 <= (int)strlen(string_argv))
				{
					int auto_inc = (string_argv[length] == '=') ? 1 : 0;
					value = atoi(&string_argv[length + auto_inc]);
				}
				else
				{
					value = 0;
				}

				bFound = true;
				continue;
			}
		}
	}

	if (bFound)
	{
		return value;
	}
	else
	{
		return 0;
	}
}

inline void string_replace(string&s1, const string&s2, const string&s3)
{
	string::size_type pos = 0;
	string::size_type a = s2.size();
	string::size_type b = s3.size();
	while ((pos = s1.find(s2, pos)) != string::npos)
	{
		s1.replace(pos, a, s3);
		pos += b;
	}
}