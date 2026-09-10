// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10
#pragma once
#include "enclosure.hpp"
#include "disks.hpp"
namespace segments {
using namespace inverse_samples;
// Preserve the exact P-bit values while reusing the destination's MPFR storage.
inline void copyDisk(disks::Disk&out,const disks::Disk&in){
    if(mpfr_get_prec(out.value.re.x)!=mpfr_get_prec(in.value.re.x)||
       mpfr_get_prec(out.value.im.x)!=mpfr_get_prec(in.value.im.x))
        throw std::runtime_error("DISK_COPY_PRECISION_MISMATCH");
    mpfr_set(out.value.re.x,in.value.re.x,MPFR_RNDN);
    mpfr_set(out.value.im.x,in.value.im.x,MPFR_RNDN);out.radius=in.radius;
}
struct Part {disks::Disk a,b;C endpoint;upper::Positive radius,critical;uint64_t iterations=0;std::string status="NOT_RUN";};
struct Result {Real scaledResidual;std::string status="INCONCLUSIVE";uint64_t iterations=0,failedLinks=0;double maxLogRatio=-INFINITY,logScaledUpper=NAN,logZoomLower=NAN,logZoomUpper=NAN,seconds=0;bool longestFirst=false;std::string failureDetail;};
inline Result evaluate(const C&c,uint64_t period,double target,const v4::Control&control){
    const auto begin=v4::Clock::now();Result result;
    const auto&seeds=quality::checkpoints;const size_t count=seeds.size();
    if(!count||seeds.front().index!=0||!mpfr_zero_p(seeds.front().value.re.x)||!mpfr_zero_p(seeds.front().value.im.x))return result;
    upper::environment();std::vector<upper::Positive> radii(count);std::vector<Part> parts(count);
    for(size_t j=1;j<count;++j){
        if(seeds[j].index<=seeds[j-1].index||seeds[j].index>=period)return result;
        double exponent=seeds[j].offset-target*std::log2(10.)+std::log2(double(j+1));
        if(!std::isfinite(exponent)||std::abs(exponent)>stms_limits::exponentLimit(v4::precision))return result;
        // The trial radius is this exact binary mantissa/exponent value.
        // Its heuristic construction is not used as an enclosure proof.
        const double whole=std::floor(exponent);
        radii[j]=upper::Positive(std::exp2(exponent-whole),static_cast<int64_t>(whole));
    }
    std::vector<size_t> order;order.reserve(count);
    for(size_t j=0;j<count;++j)order.push_back(j);
    const auto length=[&](size_t j){return (j+1<count?seeds[j+1].index:period)-seeds[j].index;};
    const auto makespan=[&]{
        std::vector<uint64_t> load(std::min<size_t>(15,count),0);
        for(auto j:order)*std::min_element(load.begin(),load.end())+=length(j);
        return *std::max_element(load.begin(),load.end());
    };
    const auto originalSpan=makespan();
    std::sort(order.begin(),order.end(),[&](size_t a,size_t b){
        const auto la=length(a),lb=length(b);return la!=lb?la>lb:a<b;
    });
    result.longestFirst=makespan()<originalSpan;
    if(!result.longestFirst)for(size_t j=0;j<count;++j)order[j]=j;
    std::atomic<size_t> next=0;std::vector<std::jthread> workers;
    for(size_t t=0;t<std::min<size_t>(15,count);++t)workers.emplace_back([&]{
        for(;;){const size_t slot=next.fetch_add(1);if(slot>=count)break;
            const size_t j=order[slot];
            auto&part=parts[j];try{
                upper::environment();enclosure::Orbit orbit;orbit.value=seeds[j].value;orbit.radiusBound=radii[j];
                orbit.iterations=seeds[j].index;
                disks::Ops diskOps;disks::Disk initialCoefficient,suffix,constant,groupA;std::vector<disks::Disk> group(64);size_t used=0;
                std::vector<upper::Positive> groupNorm(group.size());
                mpfr_set_ui(part.a.value.re.x,1,MPFR_RNDN);
                auto flush=[&]{
                    if(!used)return;
                    if(used==1){copyDisk(groupA,group[0]);mpfr_set_ui(constant.value.re.x,1,MPFR_RNDN);mpfr_set_zero(constant.value.im.x,1);constant.radius={};}
                    else{
                        copyDisk(suffix,group[used-1]);copyDisk(constant,suffix);diskOps.addOne(constant);
                        for(size_t k=used-2;k>0;--k){diskOps.multiplyKnownNorms(suffix,suffix,group[k],upper::norm(suffix.value),groupNorm[k]);diskOps.add(constant,constant,suffix);}
                        diskOps.multiplyKnownNorms(groupA,suffix,group[0],upper::norm(suffix.value),groupNorm[0]);
                    }
                    const auto norm=upper::norm(groupA.value);
                    diskOps.multiplyKnownNorms(part.a,groupA,part.a,norm,disks::Ops::l1(part.a.value));
                    diskOps.multiplyKnownNorms(part.b,groupA,part.b,norm,disks::Ops::l1(part.b.value));diskOps.add(part.b,part.b,constant);used=0;
                };
                const uint64_t end=j+1<count?seeds[j+1].index:period;
                while(orbit.iterations<end){
                    if((part.iterations&255)==0&&control.stopped()){orbit.status="TIME_BUDGET";break;}
                    auto&coefficient=orbit.iterations?group[used]:initialCoefficient;
                    const bool grouped=orbit.iterations!=0;
                    mpfr_mul_2ui(coefficient.value.re.x,orbit.value.re.x,1,MPFR_RNDN);
                    mpfr_mul_2ui(coefficient.value.im.x,orbit.value.im.x,1,MPFR_RNDN);coefficient.radius=upper::twice(orbit.radiusBound);
                    if(!orbit.iterations){
                        diskOps.multiply(part.b,coefficient,part.b);diskOps.addOne(part.b);
                    }else{++used;}
                    if(!orbit.step(c))break;++part.iterations;
                    // Exact multiplication by two preserves precision. Scale the
                    // pre-step norm bound rather than extracting it again.
                    if(grouped){groupNorm[used-1]=upper::twice(orbit.previousCenterNorm);if(used==group.size())flush();}
                }
                if(orbit.iterations==end&&orbit.status=="COMPLETE")flush();
                part.endpoint=std::move(orbit.value);part.radius=orbit.radiusBound;part.critical=orbit.criticalBound;part.status=orbit.status;
                if(orbit.iterations!=end&&part.status=="COMPLETE")part.status="INCOMPLETE";
            }catch(const std::exception&e){part.status=e.what();}
        }
        mpfr_free_cache();
    });
    for(auto&t:workers)t.join();
    bool complete=true;auto critical=upper::Positive::power(0);
    for(size_t j=0;j<parts.size();++j){const auto&part=parts[j];result.iterations+=part.iterations;complete=complete&&part.status=="COMPLETE";critical=upper::mul(critical,part.critical);
        if(part.status!="COMPLETE"&&result.failureDetail.empty())result.failureDetail="segment="+std::to_string(j)+" start="+std::to_string(seeds[j].index)+" steps="+std::to_string(part.iterations)+" reason="+part.status+" radiusExponent="+std::to_string(part.radius.e);
    }
    if(complete){
        for(size_t j=0;j+1<count;++j){
            // Component subtraction rounded away from zero bounds its magnitude.
            C delta;mpfr_sub(delta.re.x,parts[j].endpoint.re.x,seeds[j+1].value.re.x,MPFR_RNDA);
            mpfr_sub(delta.im.x,parts[j].endpoint.im.x,seeds[j+1].value.im.x,MPFR_RNDA);
            auto needed=upper::add(upper::norm(delta),parts[j].radius).real();auto available=radii[j+1].real();
            if(mpfr_cmp(needed.x,available.x)>0)++result.failedLinks;
            Real ratio;mpfr_div(ratio.x,needed.x,available.x,MPFR_RNDU);
            if(!mpfr_zero_p(ratio.x))result.maxLogRatio=std::max(result.maxLogRatio,v4::log10value(ratio));
        }
        auto scaled=upper::mul(critical,upper::add(upper::norm(parts.back().endpoint),parts.back().radius)).real();
        result.scaledResidual=scaled;
        result.logScaledUpper=mpfr_zero_p(scaled.x)?-INFINITY:v4::log10value(scaled);
        result.status=result.failedLinks?"CONNECTION_FAILED":mpfr_cmp_d(scaled.x,1e-12)<0?"RESIDUAL_BOUND_PASS":"RESIDUAL_BOUND_FAILED";
        if(result.status=="RESIDUAL_BOUND_PASS"){
            disks::Disk A,D,AD;disks::Ops ops;mpfr_set_ui(A.value.re.x,1,MPFR_RNDN);
            for(const auto&part:parts){
                if(control.stopped()){result.status="TIME_BUDGET";break;}
                ops.multiply(D,part.a,D);ops.add(D,D,part.b);ops.multiply(A,A,part.a);
            }
            if(result.status!="TIME_BUDGET"){
                ops.multiply(AD,A,D);Real lo,hi;auto radius=AD.radius.real();
                mpfr_hypot(lo.x,AD.value.re.x,AD.value.im.x,MPFR_RNDD);mpfr_sub(lo.x,lo.x,radius.x,MPFR_RNDD);
                mpfr_hypot(hi.x,AD.value.re.x,AD.value.im.x,MPFR_RNDU);mpfr_add(hi.x,hi.x,radius.x,MPFR_RNDU);
                if(mpfr_sgn(lo.x)>0){
                    mpfr_log10(lo.x,lo.x,MPFR_RNDD);mpfr_add_ui(lo.x,lo.x,2,MPFR_RNDD);
                    mpfr_log10(hi.x,hi.x,MPFR_RNDU);mpfr_add_ui(hi.x,hi.x,2,MPFR_RNDU);
                    result.logZoomLower=mpfr_get_d(lo.x,MPFR_RNDD);result.logZoomUpper=mpfr_get_d(hi.x,MPFR_RNDU);
                }
                result.status=std::isfinite(result.logZoomLower)&&std::isfinite(result.logZoomUpper)&&result.logZoomLower>0&&result.logZoomUpper-result.logZoomLower<=1e-8?"FULL_BOUND_PASS":"ZOOM_INCONCLUSIVE";
            }
        }
    }else result.status=control.stopped()?"TIME_BUDGET":"SEGMENT_INCONCLUSIVE";
    result.seconds=std::chrono::duration<double>(v4::Clock::now()-begin).count();return result;
}
}
