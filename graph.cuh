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

//??
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
}