#pragma one

#include <random>
#include <vector>

typedef float bithash_t;

class BitHash{
public:
    std::vector<bithash_t> hash_matrix;
    int p,k;
    BitHash(){}
    BitHash(int p,int k,int seed=123):p(p),k(k){
        std::default_random_engine generator (seed+1)；
        std::normal_distribution<bithash_t>
        hash_matrix.resize(p*k);
        for(int i=0;i<hash_matrix.size();i++){
            hash_matrix[i]=distribution(generator);
        }
    }

    std::vector<bool> hash2vecbool(const std::vector<std::pair<int,value_t>>& point){
        std::vector<bool> ret(k);
        for(int i=0;i<k;i++){
            bithash_t sum=0;
            for(const auto& pp:point){
                sum+=hash_matrix[i*p+pp.first]*pp.second;
            }
            ret[i]=(sum>=0);
        }
        //std::move() api
        return std::move(ret);
    }
    //uint8_t is defined by typedef = unsigned char(1 byte = 8 bits)
    uint8_t hash2uint8(const std::vector<pair<int,value_t>>& point){
        uint8_t ret=0;
        for(int i=0;i<k;i++){
            bithash_t sum=0;
            for(const auto& pp:point){
                sum+=hash_matrix[i*p+pp.first]*pp.second;
            }
            ret|=((sum>=0)<<i);
        }
        return std::move(ret);
    }

    std::vector<std::pair<int,value_t>> hash2kv(const std::vector<pair<int,value_t>>& point){
        std::vector<std::pair<int,value_t>> ret(k/sizeof(value_t)/8);
        #ifdef __ENABLE_HASH
            for(int i=0;i<k/sizeof(value_t)/8;i++){
                ret[i].first=i;
                ret[i].second=0;
            }
            for(int i=0;i<k;i++){
                bishash_t sum=0;
                for(const auto& pp: point){
                    sum+=hash_matrix[i*p+pp.first]*pp.second;
                }
                rt[i/sizeof(value_t)/8].second|=(sum>=0)?(1<<i%(sizeof(value_t)*8)):0;
            }
        #endif
        return std::move(ret);
    }
    
};
