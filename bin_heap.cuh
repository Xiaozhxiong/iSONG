//
// Created by xiaoz on 2023/1/6.
//

#ifndef ISONG_BIN_HEAP_H
#define ISONG_BIN_HEAP_H

//[begin,end)

template<typename T>
__device__ void push_heap(T *begin, T *end) {
    T *now = end - 1;
    int parent = (now - begin - 1) / 2;
    while (parent >= 0) {
        if (*(begin + parent) < *now) {
            auto tmp = *now;
            *now = *(begin + parent);
            *(begin + parent) = tmp;
            now = begin + parent;
            parent = (parent - 1) / 2;
        } else {
            break;
        }
    }
}

template<typename T>
__device__ T pop_heap(T* begin,T* end){
    
}
#endif //ISONG_BIN_HEAP_H
