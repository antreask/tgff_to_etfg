#include <iostream>
#include <fstream>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/range/algorithm/replace_if.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <random>
#include <regex>
#include <filesystem>
#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <queue>

#ifdef __unix__
    #define IS_POSIX 1
    #include <unistd.h>
#else
    #define IS_POSIX 0
#endif

using namespace std;
namespace fs = std::filesystem;

const string red_text("\033[0;31m");
const string reset("\033[0m");
const std::string WHITESPACE = " \n\r\t\f\v";
const std::string WHITESPACE_EXCLUDING_SPACES = "\n\r\t\f\v";
const bool do_topological_sorting=1;
const bool unconstrained=1;


enum Executed { e=1,g=2,egc=3};//fixed allocation at the edge only (1), hub only (2), or no constraints (3) can run at the edge, hub (ground), or cloud
struct Vertex
{
    size_t id=0,tfg_id=0;
    string type = "Null";

    vector <int> ancestors{-1}; //means no ancestor
    vector <int> predecessors{-1}; //means no predecessor
    vector <int> etfg_ids{-1};// the corresponding ids in the ETFG, e.g 1 ->1e,1g,1c => 1 -> 1,2,3
    vector <int> etfg_ancestors{-1};// the corresponding ids in the ETFG, e.g 1 ->1e,1g,1c => 1 -> 1,2,3

    size_t where = egc;

    double cloud_latency = 0.0;
    double mi_latency = 0.0;
    double t490_latency = 0.0;
    double jetson_latency = 0.0;
    double odroid_latency = 0.0;
    double rpi_latency = 0.0;

    double cloud_power = 0.0; 
    double mi_power = 0.0;
    double t490_power = 0.0;
    double jetson_power = 0.0;
    double odroid_power = 0.0;
    double rpi_power = 0.0;

    double cloud_energy = 0.0; 
    double mi_energy = 0.0;
    double t490_energy = 0.0;
    double jetson_energy = 0.0;
    double odroid_energy = 0.0;
    double rpi_energy = 0.0;

    double ram = 0.0;
    double disk = 0.0;
    double load = 0.0;
};

vector <string> tokenize(string str,char delimeter,const int& line_num)
{
    vector <string> tokens;
    stringstream check(str);       
    string intermediate; 
    while(getline(check, intermediate, delimeter))     
        tokens.push_back(intermediate); 
    return tokens;
}

//clear vector, the objects are destroyed but the capacity of the vector remains the same
template<typename T>
void clear(std::vector<T>& vec) {
 vec.clear();
}
//swap vector with an empy one
template<typename T>
void empty_swap(std::vector<T>& vec) {
 std::vector<T>().swap(vec);
}

// Why not setup a lambda you can use again & again
auto removeByIndex = []<class T>(std::vector<T> &vec, unsigned int index)
{
  // This is the meat & potatoes
    vec.erase(vec.begin() + index);
};

double rand_float(double a, double b) {
  return ((double)rand() / RAND_MAX) * (b - a) + a;
}

double linear_conversion(double OldMin, double OldMax, double NewMin, double NewMax, double OldValue) {
  return (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin;
}

//percentage range          e.g         1%                      3%              total of 100 nodes
size_t select_random_nodes(const size_t& init_range,const size_t& end_range, const size_t& tot_nodes, size_t &fixed_percentage)
{
    /* initialize random seed: */
  srand(time(NULL));
  /* generate secret number between init_range and end_range: */
  fixed_percentage=rand() % end_range + init_range;
  return round(tot_nodes*fixed_percentage/100.0);
}

size_t select_random_nodes(const size_t& tot_nodes, size_t &fixed_percentage)
{
  return round(tot_nodes*fixed_percentage/100.0);
}

string extractCharacters(string str) 
{ 
    string res;
    for(size_t i = 0; i < str.size(); ++i)
    {
        if (((str[i] >= 'a' && str[i]<='z') || (str[i] >= 'A' && str[i]<='Z')))
        {
            res+=str[i];
        }
    } 
    return res;

}

void ReplaceStringInPlace(std::string& subject, const std::string& search,
  const std::string& replace) {
    size_t pos = 0;
    while((pos = subject.find(search, pos)) != std::string::npos) {
       subject.replace(pos, search.length(), replace);
       pos += replace.length();
   }
}

template<typename T>
void print_equal(T x) {
    std::cout << std::setw(4) << x;
}

template<typename T>
void remove_duplicates(vector<T> &vec)
{
    sort( vec.begin(), vec.end() );
    vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
}


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v) {

  // initialize original index locations
  vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  // using std::stable_sort instead of std::sort
  // to avoid unnecessary index re-orderings
  // when v contains elements of equal values 
  stable_sort(idx.begin(), idx.end(),
     [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

  return idx;
}

bool compareById(const Vertex &a, const Vertex &b)
{
    return a.id < b.id;
}

template<typename T>
typename T::value_type most_frequent_element(T const& v,size_t &maxFrequency)
{    // Precondition: v is not empty
    std::map<typename T::value_type, int> frequencyMap;
    maxFrequency = 0;
    typename T::value_type mostFrequentElement{};
    for (auto&& x : v)
    {
        size_t f = ++frequencyMap[x];
        if (f > maxFrequency)
        {
            maxFrequency = f;
            mostFrequentElement = x;
        }
    }

    return mostFrequentElement;
}

/******************************* Topological Sorting *******************************/

// Data structure to store graph edges
struct Edge {
    size_t src, dest;
};

// Class to represent a graph object
class Graph
{
public:
    // construct a vector of vectors to represent an adjacency list
    vector<vector<size_t>> adjList;

    // stores indegree of a vertex
    vector<size_t> indegree;

    // Graph Constructor
    Graph(vector<Edge> const &edges, size_t N)
    {
        // resize the vector to N elements of type vector<int>
        adjList.resize(N);

        // initialize indegree
        vector<size_t> temp(N, 0);
        indegree = temp;

        // add edges to the directed graph
        for (auto &edge: edges)
        {
            // add an edge from source to destination
            adjList[edge.src].push_back(edge.dest);

            // increment in-degree of destination vertex by 1
            indegree[edge.dest]++;
        }
    }
};

// performs Topological Sort on a given DAG
bool doTopologicalSort(Graph const &graph, vector<size_t> &L, size_t N)
{
    vector<size_t> indegree = graph.indegree;

    // Set of all nodes with no incoming edges
    vector<size_t> S;
    for (size_t i = 0; i < N; i++) {
        if (!indegree[i]) {
            S.push_back(i);
        }
    }

    while (!S.empty())
    {
        // remove a node n from S
        size_t n = S.back();
        S.pop_back();

        // add n to tail of L
        L.push_back(n);

        for (size_t m : graph.adjList[n])
        {
            // remove edge from n to m from the graph
            indegree[m] -= 1;

            // if m has no other incoming edges then
            // insert m into S
            if (!indegree[m]) {
                S.push_back(m);
            }
        }
    }

    // if graph has edges then graph has at least one cycle
    for (size_t i = 0; i < N; i++) {
        if (indegree[i]) {
            return false;
        }
    }

    return true;
}

/******************************* End of Topological Sorting *******************************/


void do_bfs(const vector<vector<size_t>> &g, vector<size_t> &levels, vector<bool> v, size_t u) 
{ 
    queue<size_t> q; 

    q.push(u); 
    v[u] = true; 

    while (!q.empty()) { 

        size_t f = q.front(); 
        q.pop(); 

        // Enqueue all adjacent of f and mark them visited 
        for (auto i = g[f].begin(); i != g[f].end(); i++) { 
            if (!v[*i]) { 
                q.push(*i); 
                v[*i] = true; 
                levels[*i]=levels[f]+1;
            } 
        } 
    } 
}


int main(int argc, char **argv)
{   
    if (argc!=2)
    {
        cout<<red_text<< "You have not defined the TGFF's file path!"<<reset<<endl;
        return -1;
    }
    string input_graph= argv[1];
    if (!std::regex_match (input_graph, std::regex(".*\\b.tgff") ))
    {
        cout<<red_text<< "The file must be a .tgff file!"<<reset<<endl;
        return -1;
    }

    vector<Vertex> vertices;

    //Read input graph file and build graph
    int line_num=0;
    
    //Examples of energy, memory, storage and latency constraints
    vector<size_t> maxEnergy = {100,100,100};//Wh
    vector<size_t> maxMainMemory = {100,100,100};//MB
    vector<size_t> maxSecondaryMemory = {100,100,100};//GB
    size_t latency_constraint=100; //ms
    if (unconstrained)
    {
    	maxEnergy = {100000,100000,100000};
	    maxMainMemory = {10000000000,10000000000,10000000000};
	    maxSecondaryMemory = {10000000000,10000000000,10000000000};
	    latency_constraint=1800000000; //5 hours
    }


    vector<size_t> fixed_edge_allocation_range = {1, 5};
    vector<size_t> fixed_hub_allocation_range = {1, 3};

    //Latency Ranges (Min and Max) between the different computational devices.
    vector<double> rpi_odroid_range_lat{1,10};
    vector<double> odroid_jetson_range_lat{1,10}; 
    vector<double> jetson_t490_range_lat{1,10};
    vector<double> t490_mipro_range_lat{1,10};
    vector<double> mipro_cloud_range_lat{1,10};

    //Energy Ranges (Min and Max) between the different computational devices
    vector<double> rpi_odroid_range_en{1,10};
    vector<double> odroid_jetson_range_en{1,10}; 
    vector<double> jetson_t490_range_en{1,10};
    vector<double> t490_mipro_range_en{1,10};
    vector<double> mipro_cloud_range_en{1,10};

    //Power ranges (Min and Max) between the different computational devices
    vector<double> rpi_actual_range_power{1,10};
    vector<double> odroid_actual_range_power{1,10};
    vector<double> jetson_actual_range_power{1,10};
    vector<double> t490_actual_range_power{1,10};
    vector<double> mipro_actual_range_power{1,10};
    vector<double> cloud_actual_range_power{1,10};


    /************************************************** Reading TGFF File ******************************************/
    ifstream myfile(input_graph);

    /*********** Reading for the first time ************/
    // Empty list that will contain the sorted elements
    vector<size_t> L;
    vector<Edge> tgff_edges;
    vector<size_t> levels;
    vector<bool> v;
    size_t number_of_tasks=0;

    if (myfile.is_open())
    {
        string line;
        while (getline (myfile,line) && (line.find("@TASK_GRAPH") == std::string::npos) ){
            line_num++;
        }
        while (getline (myfile,line) && line.length()!=0){line_num++;}

        //Read tasks/nodes
        while (getline (myfile,line) && line.length()!=0)
        {
            number_of_tasks++;
        }


         //Read Arcs/Edges e.g ARC a0_3528  FROM t0_1604  TO  t0_2204 TYPE 43
        while (getline (myfile,line) && line.length()!=0)
        {
            vector <string> tokens,task_id,ancestor;
            boost::algorithm::trim_all(line);
            line = boost::replace_if(line, boost::is_any_of(WHITESPACE), ' '); 
            tokens = tokenize(line,' ',line_num);
            boost::split(ancestor,tokens[3],boost::is_any_of("_"));
            boost::split(task_id,tokens[5],boost::is_any_of("_"));
            
            /*********** topological sorting ***********/
            Edge tgff_edge;
            tgff_edge.src=stoul(ancestor[1]);
            tgff_edge.dest=stoul(task_id[1]);

            tgff_edges.push_back(tgff_edge);

            /*********** end topological sorting ***********/

        }

        /*********** topological sorting ***********/
         // create a graph from edges
        Graph tgff_graph(tgff_edges, number_of_tasks);

            // Perform Topological Sort
        if (doTopologicalSort(tgff_graph, L, number_of_tasks))
        {
                // print topological order
                /*cout << "\nTopological Sorting: \n";
                for (size_t i=0; i<L.size();i++)
                    print_equal(i);
                cout << endl;
                for (size_t i: L)
                   print_equal(i);
                cout <<endl;*/
            cout << "\nTopological Sorting Procedure has finished!\n";

        } else {
            cout << "Graph has at least one cycle. "
            "Topological sorting is not possible";
            cout <<endl;
        }
        /*********** end topological sorting ***********/

        /*********** End Reading for the first time ************/
        cout << "\nTopological Sorting Run Successfully!\n";

        if (!do_topological_sorting)
        {
        	sort(L.begin(),L.end());
        }



        /*********** BFS ************/

        v.assign(number_of_tasks, false); 
        levels.assign(number_of_tasks, 1); 

        for (size_t i = 0; i < number_of_tasks; i++) { 
            if (!v[i]) 
                do_bfs(tgff_graph.adjList,levels,v,i); 
        }

        cout << "\nBFS Run Successfully!\n";

        /*********** End of BFS ************/

    }
    myfile.close();
    

    /*********** Reading for the second time ************/
    myfile.open(input_graph);
    string rawpath=input_graph.substr(0, input_graph.find_last_of("."));
    string new_file = rawpath+"_edited.tgff";
    if (myfile.is_open())
    {
        vector<size_t> indices;
        for (auto i: sort_indexes(L)) {
            indices.push_back(i);
        }

        std::ofstream ofs(new_file);
        string line;
        while (getline (myfile,line))
        {
            if (line.length()<2){
                ofs << line <<endl;
                continue;
            }
            std::regex task_re("\\bTASK\\b");
            std::regex arc_re("\\bARC\\b");
            std::smatch m;

            string line_to_modify=line;

            vector <string> tokens,task_id, ancestor;

            boost::algorithm::trim_all(line_to_modify);
            line_to_modify = boost::replace_if(line_to_modify, boost::is_any_of(WHITESPACE), ' ');
            tokens = tokenize(line_to_modify,' ',0);

            if (std::regex_search(line_to_modify, m, task_re)) 
            {                
                boost::split(task_id,tokens[1],boost::is_any_of("_"));
                //string search = "t0_"+std::to_string(task_id[1]);
                //string replace = "K0_"+std::to_string(indices[task_id[1]]);

                string search = "t0_"+boost::lexical_cast<std::string>(task_id[1]);
                string replace = "K0_"+boost::lexical_cast<std::string>(indices[stoul(task_id[1])]);
                
                std::regex search_re("\\b"+search+"\\b");
                line=std::regex_replace(line, search_re, replace);

            }
            else if (std::regex_search(line_to_modify, m, arc_re))
            {
                boost::split(ancestor,tokens[3],boost::is_any_of("_"));
                /*string search = "t0_"+std::to_string(ancestor[1]);
                string replace = "K0_"+std::to_string(indices[ancestor[1]]);*/
                string search = "t0_"+boost::lexical_cast<std::string>(ancestor[1]);
                string replace = "K0_"+boost::lexical_cast<std::string>(indices[stoul(ancestor[1])]);

                std::regex search_re("\\b"+search+"\\b");
                line=std::regex_replace(line, search_re, replace);

                boost::split(task_id,tokens[5],boost::is_any_of("_"));
                // search = "t0_"+std::to_string(task_id[1]);
                // replace = "K0_"+std::to_string(indices[task_id[1]]);

                search = "t0_"+boost::lexical_cast<std::string>(task_id[1]);
                replace = "K0_"+boost::lexical_cast<std::string>(indices[stoul(task_id[1])]);
                search_re="\\b"+search+"\\b";
                line=std::regex_replace(line, search_re, replace);
            }
            else
            {
                ofs << line <<endl;
                continue;
            }
            
            ofs << line <<endl;
            line.clear();
        }
        ofs.close();

    }
    myfile.close();
    /*********** End Reading for the second time ************/

    cout << "\nSecond Pass Successfully Run!\n";

    myfile.open(new_file);
    if (myfile.is_open())
    {
        string line;        
        //not found
        while (getline (myfile,line) && (line.find("@TASK_GRAPH") == std::string::npos) ){
            line_num++;
        }
        while (getline (myfile,line) && line.length()!=0){line_num++;}

        //Initial Seed
        srand(time(NULL));
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> random_rpi_lat(30.0, 500.0);
        std::uniform_real_distribution<> random_rpi_power(rpi_actual_range_power[0],rpi_actual_range_power[1]); //Watts
        std::uniform_real_distribution<> random_ram(1.0, 100.0);
        std::uniform_real_distribution<> random_disk(1.0, 200.0);
        std::uniform_real_distribution<> random_load(0.1, 18.0);
        std::uniform_int_distribution<> random_edge_allocation(fixed_edge_allocation_range[0],fixed_edge_allocation_range[1]);
        std::uniform_int_distribution<> random_hub_allocation(fixed_hub_allocation_range[0],fixed_hub_allocation_range[1]);

        cout << "\nReading Nodes!\n";
        //Read Tasks/Nodes
        while (getline (myfile,line) && line.length()!=0)
        {
            vector <string> tokens,task_ids;
            boost::algorithm::trim_all(line);
            line = boost::replace_if(line, boost::is_any_of(WHITESPACE), ' ');
            tokens = tokenize(line,' ',line_num);
            boost::split(task_ids,tokens[1],boost::is_any_of("_"));
            Vertex task;
            task.id=stoul(task_ids[1]);
            /************************************************** Generate Pseudo-Random Values ******************************************/
            
            //Raspberry
            task.rpi_latency = random_rpi_lat(gen);
            task.rpi_power= random_rpi_power(gen);   
            task.rpi_energy = task.rpi_power*task.rpi_latency*2.77777778e-7;  //linear_conversion( 50.0,2000.0,  0.01,3.0,rand_float(0.01,3.0));

            double random_num=0.0,remap_num=0.0;

            //Odroid
            random_num=rand_float( rpi_odroid_range_lat[0],rpi_odroid_range_lat[1]);
            remap_num=linear_conversion( rpi_odroid_range_lat[0],rpi_odroid_range_lat[1],  rpi_odroid_range_en[0],rpi_odroid_range_en[1],random_num);
            task.odroid_latency = task.rpi_latency/random_num;
            task.odroid_power= task.rpi_power * remap_num;
            task.odroid_energy = (task.odroid_power)*task.odroid_latency*2.77777778e-7; //W*ms to Wh

            //Jetson
            random_num=rand_float( odroid_jetson_range_lat[0],odroid_jetson_range_lat[1]);            
            task.jetson_latency = task.odroid_latency/random_num;
            bool lower_bound_change=0, upper_bound_change=0;
            if ((task.odroid_power*odroid_jetson_range_en[0]) < jetson_actual_range_power[0])
                lower_bound_change=1;
            if ((task.odroid_power*odroid_jetson_range_en[1]) > jetson_actual_range_power[1])
                upper_bound_change=1;

            if (lower_bound_change && upper_bound_change) {
                remap_num=linear_conversion( odroid_jetson_range_lat[0],odroid_jetson_range_lat[1],  jetson_actual_range_power[0]/task.odroid_power,jetson_actual_range_power[1]/task.odroid_power,random_num);
            } else if (lower_bound_change) { // no need to test !b here - b==true would be the first case
                remap_num=linear_conversion( odroid_jetson_range_lat[0],odroid_jetson_range_lat[1],  jetson_actual_range_power[0]/task.odroid_power,odroid_jetson_range_en[1],random_num);
            } else if (upper_bound_change) { //no need to test !a here - that would be the first case
                remap_num=linear_conversion( odroid_jetson_range_lat[0],odroid_jetson_range_lat[1],  odroid_jetson_range_en[0],jetson_actual_range_power[1]/task.odroid_power,random_num);
            } else { // !a&&!b - the last remaining
                remap_num=linear_conversion( odroid_jetson_range_lat[0],odroid_jetson_range_lat[1],  odroid_jetson_range_en[0],odroid_jetson_range_en[1],random_num);
            }       

            task.jetson_power= task.odroid_power * remap_num;
            task.jetson_energy = (task.jetson_power)*task.jetson_latency*2.77777778e-7; //W*ms to Wh

            //T490
            random_num=rand_float( jetson_t490_range_lat[0],jetson_t490_range_lat[1]);
            remap_num=linear_conversion( jetson_t490_range_lat[0],jetson_t490_range_lat[1],  jetson_t490_range_en[0],jetson_t490_range_en[1],random_num);
            task.t490_latency = task.jetson_latency/random_num;
            task.t490_power= task.jetson_power * remap_num;
            task.t490_energy = (task.t490_power)*task.t490_latency*2.77777778e-7; //W*ms to Wh

            //Mi Pro

            //Reset flags
            lower_bound_change=0;
            upper_bound_change=0;

            random_num=rand_float( t490_mipro_range_lat[0],t490_mipro_range_lat[1]);           
            task.mi_latency = task.t490_latency/random_num;

            if ((task.t490_power*t490_mipro_range_en[0]) < mipro_actual_range_power[0])
                lower_bound_change=1;
            if ((task.t490_power*t490_mipro_range_en[1]) > mipro_actual_range_power[1])
                upper_bound_change=1;

            if (lower_bound_change && upper_bound_change) {
                remap_num=linear_conversion( t490_mipro_range_lat[0],t490_mipro_range_lat[1],  mipro_actual_range_power[0]/task.t490_power,mipro_actual_range_power[1]/task.t490_power,random_num);
            } else if (lower_bound_change) { // no need to test !b here - b==true would be the first case
                remap_num=linear_conversion( t490_mipro_range_lat[0],t490_mipro_range_lat[1],  mipro_actual_range_power[0]/task.t490_power,t490_mipro_range_en[1],random_num);
            } else if (upper_bound_change) { //no need to test !a here - that would be the first case
                remap_num=linear_conversion( t490_mipro_range_lat[0],t490_mipro_range_lat[1],  t490_mipro_range_en[0],mipro_actual_range_power[1]/task.t490_power,random_num);
            } else { // !a&&!b - the last remaining
                remap_num=linear_conversion( t490_mipro_range_lat[0],t490_mipro_range_lat[1],  t490_mipro_range_en[0],t490_mipro_range_en[1],random_num);
            }

            task.mi_power= task.t490_power * remap_num;
            task.mi_energy = (task.mi_power)*task.mi_latency*2.77777778e-7; //W*ms to Wh

            //Cloud
            random_num=rand_float( mipro_cloud_range_lat[0],mipro_cloud_range_lat[1]);
            remap_num=linear_conversion( mipro_cloud_range_lat[0],mipro_cloud_range_lat[1],  mipro_cloud_range_en[0],mipro_cloud_range_en[1],random_num);
            task.cloud_latency = task.mi_latency/random_num;
            task.cloud_power= task.mi_power * remap_num;
            task.cloud_energy = (task.cloud_power)*task.cloud_latency*2.77777778e-7; //W*ms to Wh


            task.ram = random_ram(gen); //MB
            task.disk = random_disk(gen); //MB
            task.load = random_load(gen); //Mb

            /************************************************** End Generate Pseudo-Random Values ******************************************/

            vertices.push_back(task);
            line_num++;
        }

        cout << "\nVertices Vector size= "<<vertices.size();
        cout << "\nReading Arcs!\n";

         //Read Arcs/Edges e.g ARC a0_3528  FROM t0_1604  TO  t0_2204 TYPE 43
        while (getline (myfile,line) && line.length()!=0)
        {
            vector <string> tokens,task_id,ancestor;
            //line=trim(line,WHITESPACE_EXCLUDING_SPACES);
            //line.erase(boost::remove_if(line, boost::is_any_of(WHITESPACE_EXCLUDING_SPACES)), line.end());
            boost::algorithm::trim_all(line);
            line = boost::replace_if(line, boost::is_any_of(WHITESPACE), ' '); 
            tokens = tokenize(line,' ',line_num);
            boost::split(ancestor,tokens[3],boost::is_any_of("_"));
            boost::split(task_id,tokens[5],boost::is_any_of("_"));
            //cout <<"U: " << ancestor[1] << "V: " << task_id[1] <<endl;

            if (vertices[stoul(task_id[1])].ancestors.size()==1 && vertices[stoul(task_id[1])].ancestors[0]==-1){ //first time
                vertices[stoul(task_id[1])].ancestors[0]=stoul(ancestor[1]);
                vertices[stoul(ancestor[1])].predecessors[0]=stoul(task_id[1]);
            }
            else{
                vertices[stoul(task_id[1])].ancestors.push_back(stoul(ancestor[1]));
                vertices[stoul(ancestor[1])].predecessors.push_back(stoul(task_id[1]));
            }
            line_num++;

        }

        /************************************************** End Reading TGFF File ******************************************/

        /************************************************** Sort Ancestors and predecessors of each task ******************************************/
        cout << "\nSorting Ancestors and Predecessors!\n";
        size_t sum_in_degree=0,sum_out_degree=0,max_in_d=0, max_out_d=0;
        double avg_in_degree=0.0,avg_out_degree=0.0;
        for (size_t i=0;i<vertices.size();i++)
        {
            sort(vertices[i].ancestors.begin(), vertices[i].ancestors.end());
            sum_in_degree+=vertices[i].ancestors.size();
            if (vertices[i].ancestors.size()>max_in_d)
                max_in_d=vertices[i].ancestors.size();
            sort(vertices[i].predecessors.begin(), vertices[i].predecessors.end());
            sum_out_degree+=vertices[i].predecessors.size();
            if (vertices[i].predecessors.size()>max_out_d)
                max_out_d=vertices[i].predecessors.size();
        }

        sort(vertices.begin(),vertices.end(),compareById);
        size_t above_unconnected_nodes=0, below_unconnected_nodes=0;
        for (size_t i=0;i<vertices.size();i++)
        {            
            if (vertices[i].predecessors.size()==1 && vertices[i].predecessors[0]==-1)
                //cout<<vertices[i].id << endl;
                below_unconnected_nodes++;
            if (vertices[i].ancestors.size()==1 && vertices[i].ancestors[0]==-1)
                //cout<<vertices[i].id << endl;
                above_unconnected_nodes++;
        }

        /************************************************** End Sort Ancestors and predecessors of each task ******************************************/

        /************************************************** Randomly Selected Fixed Allocation ******************************************/
        cout << "\nRandomly Selecting Fixed Allocation!\n";
         //copy initial vector
        vector<Vertex> copy_vertices(vertices),fixed_edge_vertices,fixed_hub_vertices,etfg;        

        size_t fixed_edge_percentage,fixed_hub_percentage;
        size_t fixed_edge_nodes = select_random_nodes(vertices.size(),fixed_edge_percentage=random_edge_allocation(gen));
        std::sample(copy_vertices.begin(), copy_vertices.end(), std::back_inserter(fixed_edge_vertices),
            fixed_edge_nodes, std::mt19937{std::random_device{}()});

        for (size_t i = fixed_edge_vertices.size()-1; i < fixed_edge_vertices.size() ; --i)
            removeByIndex(copy_vertices,fixed_edge_vertices[i].id);

        size_t fixed_hub_nodes = select_random_nodes(vertices.size(),fixed_hub_percentage=random_hub_allocation(gen));
        std::sample(copy_vertices.begin(), copy_vertices.end(), std::back_inserter(fixed_hub_vertices),
            fixed_hub_nodes, std::mt19937{std::random_device{}()});

        for (size_t i=0;i<fixed_edge_vertices.size();i++)
            vertices[fixed_edge_vertices[i].id].where=e;
        for (size_t i=0;i<fixed_hub_vertices.size();i++)
            vertices[fixed_hub_vertices[i].id].where=g;


        size_t index=0;
        for (size_t i=0;i<vertices.size();i++)
        {

            if (vertices[i].where==egc)
            {
                empty_swap(vertices[i].etfg_ids);
                for (size_t j=0;j<3;j++)
                {
                    etfg.push_back(vertices[i]);
                    etfg[index].id=index+1;
                    if (j==0)
                        etfg[index].type= std::to_string(i+1)+"e";
                    else if (j==1)
                        etfg[index].type= std::to_string(i+1)+"g";
                    else
                        etfg[index].type= std::to_string(i+1)+"c";
                    vertices[i].etfg_ids.push_back(index);
                    index++;
                }
            }
            else
            {
                etfg.push_back(vertices[i]);
                etfg[index].id=index+1;
                if (etfg[index].where==e)
                    etfg[index].type= std::to_string(i+1)+"e";
                else
                    etfg[index].type= std::to_string(i+1)+"g";
                empty_swap(vertices[i].etfg_ids);
                vertices[i].etfg_ids.push_back(index);
                index++;

            }
        }

        for (size_t i=0;i<etfg.size();i++)
        {
            if (etfg[i].ancestors.size()==1 && etfg[i].ancestors[0]==-1)
                continue;
            else
                empty_swap(etfg[i].etfg_ancestors);

            for (size_t j=0;j<etfg[i].ancestors.size();j++){                
                etfg[i].etfg_ancestors.insert(etfg[i].etfg_ancestors.end(), vertices[etfg[i].ancestors[j]].etfg_ids.begin(), vertices[etfg[i].ancestors[j]].etfg_ids.end());
            }
        }

        //remove duplicates TODO: Fix duplicates
        for (size_t i=0;i<etfg.size();i++)
        {
            remove_duplicates(etfg[i].etfg_ancestors);
        }

        /************************************************** End Randomly Selected Fixed Allocation ******************************************/

        /************************************************** Create Extended Task flow Graphs ******************************************/
        //Create folders, take working directory
        string pwd;
        if (IS_POSIX == 1) { 
            char cwd[1024];           
            getcwd(cwd, sizeof(cwd));
            pwd=cwd;
            //puts("Path info by use getcwd():");
            //printf("\tWorkdir: %s\n", cwd);
            //printf("\tFilepath: %s/%s\n", cwd, __FILE__);
        }
        
        cout <<endl << "PWD= "<< pwd<<endl;
        string energy_folder=pwd+"/coins_journal_graphs/optimize_energy/"+ rawpath.substr(input_graph.find_last_of("/")+1) +"/";
        string performance_folder=pwd+"/coins_journal_graphs/optimize_performance/"+ rawpath.substr(input_graph.find_last_of("/")+1) +"/";
        cout <<endl << "Creating energy_folder "<< energy_folder<<endl; 
        cout <<endl << "Creating performance_folder "<< performance_folder<<endl;  
        fs::create_directories(energy_folder);
        fs::create_directories(performance_folder);

        std::ofstream rpi_mi (performance_folder+"/1rpi_mi.txt", std::ofstream::out);
        std::ofstream rpi_t490 (performance_folder+"/2rpi_t490.txt", std::ofstream::out);
        std::ofstream odroid_mi (performance_folder+"/3odroid_mi.txt", std::ofstream::out);
        std::ofstream odroid_t490 (performance_folder+"/4odroid_t490.txt", std::ofstream::out);
        std::ofstream jetson_mi (performance_folder+"/5jetson_mi.txt", std::ofstream::out);
        std::ofstream jetson_t490 (performance_folder+"/6jetson_t490.txt", std::ofstream::out);

        std::ofstream makefile (performance_folder+"/Makefile", std::ofstream::out);

        std::ofstream rpi_mi_en (energy_folder+"/1rpi_mi.txt", std::ofstream::out);
        std::ofstream rpi_t490_en (energy_folder+"/2rpi_t490.txt", std::ofstream::out);
        std::ofstream odroid_mi_en (energy_folder+"/3odroid_mi.txt", std::ofstream::out);
        std::ofstream odroid_t490_en (energy_folder+"/4odroid_t490.txt", std::ofstream::out);
        std::ofstream jetson_mi_en (energy_folder+"/5jetson_mi.txt", std::ofstream::out);
        std::ofstream jetson_t490_en (energy_folder+"/6jetson_t490.txt", std::ofstream::out);

        std::ofstream makefile_en (energy_folder+"/Makefile", std::ofstream::out);

        makefile << "all:\n";
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/1rpi_mi.txt"<<endl;
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/2rpi_t490.txt"<<endl;
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/3odroid_mi.txt"<<endl;
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/4odroid_t490.txt"<<endl;
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/5jetson_mi.txt"<<endl;
        makefile << "\t~/parser_allocator/linux64/examples/build/pars_allocator_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_performance/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/6jetson_t490.txt";

        makefile_en << "all:\n";
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/1rpi_mi.txt"<<endl;
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/2rpi_t490.txt"<<endl;
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/3odroid_mi.txt"<<endl;
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/4odroid_t490.txt"<<endl;
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/5jetson_mi.txt"<<endl;
        makefile_en << "\t~/parser_allocator/linux64/examples/build/pars_allocator_energy_c++ ~/tgff_to_etfg/coins_journal_graphs/optimize_energy/"<<rawpath.substr(input_graph.find_last_of("/")+1)<<"/6jetson_t490.txt";

        //no scientific notation
        rpi_mi <<std::fixed;
        rpi_t490 <<std::fixed;
        odroid_mi <<std::fixed;
        odroid_t490 <<std::fixed;
        jetson_mi <<std::fixed;
        jetson_t490 <<std::fixed;

        rpi_mi_en <<std::fixed;
        rpi_t490_en <<std::fixed;
        odroid_mi_en <<std::fixed;
        odroid_t490_en <<std::fixed;
        jetson_mi_en <<std::fixed;
        jetson_t490_en <<std::fixed;

        rpi_mi << etfg.size() << " 3\n";
        rpi_t490 << etfg.size() << " 3\n";
        odroid_mi << etfg.size() << " 3\n";
        odroid_t490 << etfg.size() << " 3\n";
        jetson_mi << etfg.size() << " 3\n";
        jetson_t490 << etfg.size() << " 3\n";

        rpi_mi_en << etfg.size() << " 3\n";
        rpi_t490_en << etfg.size() << " 3\n";
        odroid_mi_en << etfg.size() << " 3\n";
        odroid_t490_en << etfg.size() << " 3\n";
        jetson_mi_en << etfg.size() << " 3\n";
        jetson_t490_en << etfg.size() << " 3\n";

        for (size_t j=0;j<maxEnergy.size();j++){
            rpi_mi << maxEnergy[j] << " ";
            rpi_t490 << maxEnergy[j] << " ";
            odroid_mi << maxEnergy[j] << " "; 
            odroid_t490 << maxEnergy[j] << " "; 
            jetson_mi << maxEnergy[j] << " "; 
            jetson_t490 << maxEnergy[j] << " ";

            rpi_mi_en << maxEnergy[j] << " ";
            rpi_t490_en << maxEnergy[j] << " ";
            odroid_mi_en << maxEnergy[j] << " "; 
            odroid_t490_en << maxEnergy[j] << " "; 
            jetson_mi_en << maxEnergy[j] << " "; 
            jetson_t490_en << maxEnergy[j] << " "; 
        }

        rpi_mi << std::endl;
        rpi_t490 << std::endl;
        odroid_mi << std::endl;
        odroid_t490 << std::endl;
        jetson_mi << std::endl;
        jetson_t490 << std::endl;

        rpi_mi_en << std::endl;
        rpi_t490_en << std::endl;
        odroid_mi_en << std::endl;
        odroid_t490_en << std::endl;
        jetson_mi_en << std::endl;
        jetson_t490_en << std::endl;

        for (size_t j=0;j<maxMainMemory.size();j++){
            rpi_mi << maxMainMemory[j] << " ";
            rpi_t490 << maxMainMemory[j] << " ";
            odroid_mi << maxMainMemory[j] << " ";
            odroid_t490 << maxMainMemory[j] << " ";
            jetson_mi << maxMainMemory[j] << " ";
            jetson_t490 << maxMainMemory[j] << " ";

            rpi_mi_en << maxMainMemory[j] << " ";
            rpi_t490_en << maxMainMemory[j] << " ";
            odroid_mi_en << maxMainMemory[j] << " ";
            odroid_t490_en << maxMainMemory[j] << " ";
            jetson_mi_en << maxMainMemory[j] << " ";
            jetson_t490_en << maxMainMemory[j] << " ";
        }
        rpi_mi << std::endl;
        rpi_t490 << std::endl;
        odroid_mi << std::endl;
        odroid_t490 << std::endl;
        jetson_mi << std::endl;
        jetson_t490 << std::endl;

        rpi_mi_en << std::endl;
        rpi_t490_en << std::endl;
        odroid_mi_en << std::endl;
        odroid_t490_en << std::endl;
        jetson_mi_en << std::endl;
        jetson_t490_en << std::endl;


        for (size_t j=0;j<maxSecondaryMemory.size();j++){
            rpi_mi << maxSecondaryMemory[j] << " ";
            rpi_t490 << maxSecondaryMemory[j] << " ";
            odroid_mi << maxSecondaryMemory[j] << " ";
            odroid_t490 << maxSecondaryMemory[j] << " ";
            jetson_mi << maxSecondaryMemory[j] << " ";
            jetson_t490 << maxSecondaryMemory[j] << " ";

            rpi_mi_en << maxSecondaryMemory[j] << " ";
            rpi_t490_en << maxSecondaryMemory[j] << " ";
            odroid_mi_en << maxSecondaryMemory[j] << " ";
            odroid_t490_en << maxSecondaryMemory[j] << " ";
            jetson_mi_en << maxSecondaryMemory[j] << " ";
            jetson_t490_en << maxSecondaryMemory[j] << " ";
        }
        rpi_mi << std::endl;
        rpi_t490 << std::endl;
        odroid_mi << std::endl;
        odroid_t490 << std::endl;
        jetson_mi << std::endl;
        jetson_t490 << std::endl;

        rpi_mi_en << std::endl;
        rpi_t490_en << std::endl;
        odroid_mi_en << std::endl;
        odroid_t490_en << std::endl;
        jetson_mi_en << std::endl;
        jetson_t490_en << std::endl;

        //Latency Constraint
        rpi_mi_en << latency_constraint << std::endl;
        rpi_t490_en<< latency_constraint << std::endl;
        odroid_mi_en << latency_constraint << std::endl;
        odroid_t490_en << latency_constraint << std::endl;
        jetson_mi_en << latency_constraint << std::endl;
        jetson_t490_en << latency_constraint << std::endl;

        for (size_t i=0;i<etfg.size();i++)
        {
            if (extractCharacters(etfg[i].type)=="e")
            {
                rpi_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].rpi_latency << "," << etfg[i].rpi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].rpi_latency << "," << etfg[i].rpi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].odroid_latency << "," << etfg[i].odroid_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].odroid_latency << "," << etfg[i].odroid_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].jetson_latency << "," << etfg[i].jetson_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].jetson_latency << "," << etfg[i].jetson_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";

                rpi_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].rpi_latency << "," << etfg[i].rpi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].rpi_latency << "," << etfg[i].rpi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].odroid_latency << "," << etfg[i].odroid_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].odroid_latency << "," << etfg[i].odroid_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].jetson_latency << "," << etfg[i].jetson_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].jetson_latency << "," << etfg[i].jetson_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";

            }
            else if (extractCharacters(etfg[i].type)=="g")
            {
                rpi_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";

                rpi_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].mi_latency << "," << etfg[i].mi_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].t490_latency << "," << etfg[i].t490_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";

            }
            else if (extractCharacters(etfg[i].type)=="c")
            {
                rpi_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490 << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";

                rpi_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                rpi_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                odroid_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_mi_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                jetson_t490_en << etfg[i].id << ";" << etfg[i].type << ";" << etfg[i].cloud_latency << "," << etfg[i].cloud_energy << "," << etfg[i].ram << "," << etfg[i].disk << ";";
                
            }

            for (size_t k=0;k<etfg[i].etfg_ancestors.size();k++)
            {
                if (k+1 != etfg[i].etfg_ancestors.size())
                {
                    rpi_mi <<etfg[i].etfg_ancestors[k]+1 << ",";
                    rpi_t490 <<etfg[i].etfg_ancestors[k]+1 << ",";
                    odroid_mi <<etfg[i].etfg_ancestors[k]+1 << ",";
                    odroid_t490 <<etfg[i].etfg_ancestors[k]+1 << ",";
                    jetson_mi <<etfg[i].etfg_ancestors[k]+1 << ",";
                    jetson_t490 <<etfg[i].etfg_ancestors[k]+1 << ",";


                    rpi_mi_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                    rpi_t490_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                    odroid_mi_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                    odroid_t490_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                    jetson_mi_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                    jetson_t490_en <<etfg[i].etfg_ancestors[k]+1 << ",";
                }
                else
                {
                    rpi_mi << etfg[i].etfg_ancestors[k]+1 << ";";
                    rpi_t490 << etfg[i].etfg_ancestors[k]+1 << ";";
                    odroid_mi << etfg[i].etfg_ancestors[k]+1 << ";";
                    odroid_t490 << etfg[i].etfg_ancestors[k]+1 << ";";
                    jetson_mi << etfg[i].etfg_ancestors[k]+1 << ";";
                    jetson_t490 << etfg[i].etfg_ancestors[k]+1 << ";";

                    rpi_mi_en << etfg[i].etfg_ancestors[k]+1 << ";";
                    rpi_t490_en << etfg[i].etfg_ancestors[k]+1 << ";";
                    odroid_mi_en << etfg[i].etfg_ancestors[k]+1 << ";";
                    odroid_t490_en << etfg[i].etfg_ancestors[k]+1 << ";";
                    jetson_mi_en << etfg[i].etfg_ancestors[k]+1 << ";";
                    jetson_t490_en << etfg[i].etfg_ancestors[k]+1 << ";";
                }
            }
            rpi_mi << etfg[i].load <<endl;
            rpi_t490 << etfg[i].load <<endl;
            odroid_mi << etfg[i].load <<endl;
            odroid_t490 << etfg[i].load <<endl;
            jetson_mi << etfg[i].load <<endl;
            jetson_t490 << etfg[i].load <<endl;

            rpi_mi_en << etfg[i].load <<endl;
            rpi_t490_en << etfg[i].load <<endl;
            odroid_mi_en << etfg[i].load <<endl;
            odroid_t490_en << etfg[i].load <<endl;
            jetson_mi_en << etfg[i].load <<endl;
            jetson_t490_en << etfg[i].load <<endl;
        }
        rpi_mi.close();
        rpi_t490.close();
        odroid_mi.close();
        odroid_t490.close();
        jetson_mi.close();
        jetson_t490.close();

        rpi_mi_en.close();
        rpi_t490_en.close();
        odroid_mi_en.close();
        odroid_t490_en.close();
        jetson_mi_en.close();
        jetson_t490_en.close();

        /************************************************** End Create Extended Task flow Graphs ******************************************/

        /************************************************** Stats ******************************************/

        size_t freq=0;

        cout <<"\nStatistics for " <<input_graph<<":\n\n";
        cout << "Unconnected Nodes -> Above: "<<above_unconnected_nodes<< " | Below: "<<below_unconnected_nodes<<endl;
        cout <<"Number of Nodes/Edges: "<<vertices.size()<<"/"<<tgff_edges.size()<<endl;
        avg_in_degree=sum_in_degree/double(vertices.size());
        avg_out_degree=sum_out_degree/double(vertices.size());
        cout <<"Avg In/Out Degree: "<<avg_in_degree << "/"<<avg_out_degree<<endl;
        cout <<"Max In/Out Degree: "<<max_in_d << "/"<<max_out_d<<endl;
        cout << "#Levels: " << levels[levels.size()-1]<<endl;
        cout << "Most_frequent_element: " << most_frequent_element(levels,freq) << " seen " << freq << " times." <<endl;
        cout << "Fixed Edge Nodes (%/#) = "<<fixed_edge_percentage<< "% /" << fixed_edge_nodes<<endl;
        cout << "Fixed Hub Nodes (%/#) = "<<fixed_hub_percentage<< "% /" << fixed_hub_nodes<<endl;
        cout << "Size of ETFG = "<< etfg.size()<< " | It should be = "<< (fixed_edge_nodes + fixed_hub_nodes)+ (vertices.size()-(fixed_edge_nodes + fixed_hub_nodes))*3 <<endl;



        fs::create_directories(pwd+"/coins_journal_graphs/graphs_properties/");
        std::ofstream stats_file (pwd+"/coins_journal_graphs/graphs_properties/"+ rawpath.substr(input_graph.find_last_of("/")+1)+"_stats.txt", std::ofstream::out);

        stats_file <<"\nStatistics for " <<input_graph<<":\n\n";
        stats_file << "Unconnected Nodes -> Above: "<<above_unconnected_nodes<< " | Below: "<<below_unconnected_nodes<<endl;
        stats_file <<"Number of Nodes/Edges: "<<vertices.size()<<"/"<<tgff_edges.size()<<endl;
        stats_file <<"Avg In/Out Degree: "<<avg_in_degree << "/"<<avg_out_degree<<endl;
        stats_file <<"Max In/Out Degree: "<<max_in_d << "/"<<max_out_d<<endl;
        stats_file << "#Levels: " << levels[levels.size()-1]<<endl;
        stats_file << "Most_frequent_element: " << most_frequent_element(levels,freq) << " seen " << freq << " times." <<endl;
        stats_file << "Fixed Edge Nodes (%/#) = "<<fixed_edge_percentage<< "% /" << fixed_edge_nodes<<endl;
        stats_file << "Fixed Hub Nodes (%/#) = "<<fixed_hub_percentage<< "% /" << fixed_hub_nodes<<endl;
        stats_file << "Size of ETFG = "<< etfg.size()<< " | It should be = "<< (fixed_edge_nodes + fixed_hub_nodes)+ (vertices.size()-(fixed_edge_nodes + fixed_hub_nodes))*3 <<endl;
        

        /************************************************** End Stats ******************************************/
    }
    myfile.close();

    return 0;
}
