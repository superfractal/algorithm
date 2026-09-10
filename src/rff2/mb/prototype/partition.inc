// Created by GPT-6 on 2026-09-10
#pragma once
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace sensitivity_partition {
template<class Points>
std::vector<size_t> boundaries(const Points&points,size_t count,unsigned long period,
                               unsigned fullbits,double goalbits){
    if(!count||count>points.size()||points.front().n!=0||!period||fullbits<128||!std::isfinite(goalbits))
        throw std::runtime_error("invalid weighted partition input");
    std::vector<double> cumulative(points.size()+1,0);
    for(size_t k=0;k<points.size();++k){
        const auto end=k+1<points.size()?points[k+1].n:period;
        if(end<=points[k].n||end>period||!std::isfinite(points[k].logA))
            throw std::runtime_error("invalid weighted partition pilot");
        const double bits=std::clamp(goalbits-points[k].logA+64.,128.,double(fullbits));
        // Approximate multiplication cost, used only for load balancing.
        const double work=double(end-points[k].n)*bits*std::sqrt(bits);
        cumulative[k+1]=cumulative[k]+work;
        if(!std::isfinite(cumulative[k+1]))throw std::runtime_error("partition weight overflow");
    }
    std::vector<size_t> result;result.reserve(count);result.push_back(0);
    for(size_t j=1;j<count;++j){
        const double target=cumulative.back()*(double(j)/double(count));
        auto it=std::lower_bound(cumulative.begin(),cumulative.end(),target);
        size_t index=size_t(it-cumulative.begin());
        if(index&&index<cumulative.size()&&target-cumulative[index-1]<cumulative[index]-target)--index;
        index=std::clamp(index,result.back()+1,points.size()-(count-j));
        result.push_back(index);
    }
    return result;
}
}
