#pragma once

#include <vector>
#include <algorithm>
#include <queue>
#include <cstdlib>
#include <random>
#include <unordered_set>

#include "config.cuh"
#include "data.cuh"

class GraphWrappper{
public:
    //pure virtual function
    virtual void add_vertex(idx_t vertex_id,std::vector<std::pair<int,value_t>>& point) =0;
    virtual void search_top_k(const std::vector<std::pair<int,value_t>>& query,int k,std::vector<idx_t>& result)=0;
    virtual void dump(std::string file="bfsg.graph")=0;
    virtual void load(std::string file="bfsg.graph")=0;

    virtual void search_top_k_batch(const std::vector<std::pair<int,value_t>>& queries,int k,std::vector<std::vector<idx_t>>& result){};

    virtual ~GraphWrapper(){};
};

//astar search in graph on CPU
template<const int dist_type>
class FixedDegreeGraph:public GrapWrapperr{
private:
    const int degree=SEARCH_DEGREE;
    const int flexible_degree=FIXED_DEGREE;
    const int vertex_offset_shift=FIXED_DEGREE_SHIFT;
    std::vector<idx_t> edges;
    std::vector<dist_t> edge_dists;
    Data* data;
    std::mt19937_64 rand_gen=std::mt19937_64(1234567);

    void rank_and_switch_ordered(idx_t v_id,idx_t u_id){
        auto curr_dist=pair_distance(v_id,u_id);
        auto offset=((size_t)v_id)<<vertex_offset_shift;
        if(curr_dist>=edge_dist[offset+edges[offset]]){
            return;
        }
        edges[offset+edges[offset]]=u_id;
        edge_dist[offset+edges[offset]]=curr_dist;
        for(size_t i=offset+edges[offset]-1;i>offset;i--){
            if(edge_dist[i]>edge_dist[i+1]){
                std::swap(edges[i],edges[i+1]);
                std::swap(edge_dist[i],edge_dist[i+1]);
            }
            else break;
        }
    }

    void rank_and_switch(idx_t v_id,idx_t u_id){
        rank_and_switch_ordered(v_id,u_id);
    }

    template<typename T>
    dist_t distance(idx_t a,T& b){
        if(dist_type==0) return data->l2_distance(a,b);
        else if(dist_type==1) return data->negative_inner_prod_distance(a,b);
        else if(dist_type==2) return data->negative_cosine_distance(a,b);
        else return data->bit_hamming_distance(a,b);
    }

    void compute_distance_naive(size_t offset,std::vector<dist_t>& dists){
        dist.resize(edges[offset]);
        auto degree=edges[offset];
        for(int i=0;i<degree;i++){
            dists[i]=distance(offset>>vertex_offset_shift,edges[offset+i+1]);
        }
    }

    void compute_distance(size_t offset,std::vector<dist_t>& dists){
        compute_distance_naive(offset,dists);
    }

    template<typename T>
    dist_t pair_distance_naive(idx_t a,T& b){
        return distance(a,b);
    }

    template<typename T>
    dist_t pair_distance(idx_t a,T &b){
        return pair_distance_naive(a,b);
    }
    //quick sort
    void qsort(size_t l,size_t r){
        auto mid=(l+r)>>1;
        int i=l,j=r;
        auto k=edge_dist[mid];
        do{
            while(edge_dist[i]<k) i++;
            while(k<edge_dist[j]) j--;
            if(i<=j){
                std::swap(edge_dist[i],edge_dist[j]);
                std::swap(edge[i],edge[j]);
                i++,j--;
            }
        }while(i<=j);
        if(i<r) qsort(i,r);
        if(l<j) qsort(l,j);
    }

    void rank_edges(size_t offset){
        std::vector<dist_t> dists;
        compute_distance(offset,dists);
        for(int i=0;i<dists.size();i++){
            edge_dist[offset+i+1]=dists[i];
        }
        qsort(offset+1,offset+dists.size());
    }

    void add_edge(idx_t v_id,idx_t u_id){
        auto offset=((size_t)v_id)<<vertex_offset_shift;
        if(edges[offset]<flexible_degree){
            edges[offset]++;
            edges[offset+edges[offset]]=uid;
            // select the vertice having the smallest dist
            if(edges[offset]==flexible_degree){
                rank_edges(offset);
            }
        }
        else{
            rank_and_switch(v_id,u_id);
        }
    }
public:
    long long total_explore_cnt=0;
    int total_explore_times=0;

    FixedDegreeGraph(Data* data): data(data){
        auto num_vertices=data->max_vertices();
        //each vertex has a sapce of length vertex_offset_shift
        edges=std::vector<idx_t>(((size_t)num_vertices)<<vertex_offset_shift);
        edge_dist=std::vector<dist_t>(((size_t)num_vertices)<<vertex_offset_shift);
    }
    // add vertex and its edges
    void add_vertex(idx_t vertex_id,std::vector<std::pair<int,value_t>>& point){
        std::vector<idx_t> neighbor;
        //search neighbors
        search_top_k(point,CONSTRUCT_SEARCH_BUDGET,neighbor);
        int num_neighbors=degree<neighbor.size()?degree:neighbor.size();
        // the start address of vertex
        auto offset=((size_t)vertex_id)<<vertex_offset_shift;
        // the first place store the number of neighbor
        edges[offset]=num_neighbors;

        for(int i=0;i<neighbor.size();i++){
            edges[offset+i+1]=neighbor[i];
        }
        rank_edges[offset];
        for(int i=0;i<neihtor.size()&&i<degree;i++){
            add_edge(neighbor[i];vertex_id);
        }
    }

    void astar_multi_start_search(const std::vector<std::pair<int,value_t>>& query,int k,std::vector<idx_t>& result){
        std::priority_queue<std::pair<dist_t,idx_t>,std::vector<std::pair<dist_t,idx_t>>,std::greater<std::pair<dist_t,idx_t>>>q;
        const int num_start_point=1;

        auto converted_query=data->organize_point(query);
        std::unordered_set<idx_t> visited;
        for(int i=0;i<num_start_point && i<data->curr_vertices();i++){
            auto start=0;
            if(visited.count(start))
                continue;
            visite.insert(start);
            q.push(std::make_pair(pair_distance_naive(start,converted_query),start));
        }
        std::priority_queue<std::pair<dist_t,idx_t>> topk;
        const int max_step=1000;
        bool found_min_node=false;
        int explore_cnt=0;
        for(int iter=0;iter<max_step&&!q.empty();iter++){
            auto now=q.top();
            if(topk.size()==k&&topk.top().first<now.first){
                break;
            }
            explore_cnt++;
            q.pop();
            topk.push(now);
            if(topk.size()>k){
                topk.pop();
            }
            auto offset=((size_t)now.second)<<vertex_offset_shift;
            auto degree=edge[offset];
            for(int i=0;i<degree;i++){
                auto start=edge[offset+i+1];
                if(visited.count(start))
                    continue; 
                q.push(std::make_pair(pair_distance_naive(start,converted_query),start));
                auto tmp=pair_distance_naive(start,converted_query);
                visited.insert(start);
            }
        }
        total_explore_cnt+=explore_cnt;
        total_explore_times++;
        result.resize(topk.size());
        int i=result.size()-1;
        while(!topk.empty()){
            result[i]=(topk.top().second);
            topk.pop();
            i--;
        }
    }

    void search_top_k(const std::vector<std::pair<int,value_t>>& query,int k,std::vector<idx_t>& result){
        astar_multi_start_search(query,k,result);
    }

    void search_top_k_batch(const std::vector<std::vector<std::pair<int,value_t>>>& query,int k,std::vector<std::vector<idx_t>>& result){
        for(int i=0;i<query.size();i++){
            search_top_k(queries[i],k,results[i]);
        }
    }

    void print_stat(){
        auto n=data->max_vertices();
        size_t sum=0;
        std::vector<size_t> histogram(2*degree+1,0);
        for(size_t i=0;i<n;i++){
            sum+=edges[i<<vertex_offset_shift];
            if(tmp>2*degree+1)
                fprintf(stderr,"[ERROR] node %zu has %d degree\n",i,tmp);
            histogram[edges[i<<vertex_offset_shift]]++;
            if(tmp!=degree)
                fprintf(stderr,"[INFO] %zu has degree %d\n",i,tmp);
        }
        fprintf(stderr,"[INFO] #vertices %zu, avg degree %f\n",n,sum*1.0/n);
        std::unordered_set<idx_t> visited;
        fprintf(stderr,"[INFO] degree histogram:\n"); 
        for(int i = 0;i <= 2 * degree + 1;++i)
            fprintf(stderr,"[INFO] %d:\t%zu\n",i,histogram[i]);
    }

    void print_edges(int x){
        for(size_t i = 0;i < x;++i){
            size_t offset = i << vertex_offset_shift;
            int degree = edges[offset];
            fprintf(stderr,"%d (%d): ",i,degree);
            for(int j = 1;j <= degree;++j)
                fprintf(stderr,"(%zu,%f) ",edges[offset + j],edge_dist[offset + j]);
            fprintf(stderr,"\n");
        }
    }

    void dump(std::string file = "bfsg.graph"){
        FILE* fp = fopen(file.c_str(),"wb");
        auto num_vertices = data->max_vertices();
        fwrite(&edges[0],sizeof(edges[0]) * (((size_t)num_vertices) << vertex_offset_shift),1,fp);
        fclose(fp);
    }

    void load(std::string file = "bfsg.graph"){
        FILE* fp = fopen(file.c_str(),"rb");
        auto num_vertices = data->max_vertices();
        auto cnt = fread(&edges[0],sizeof(edges[0]) * (((size_t)num_vertices) << vertex_offset_shift),1,fp);
        fclose(fp);
    }
};