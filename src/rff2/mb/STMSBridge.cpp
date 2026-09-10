// Created by GPT-6 on 2026-09-10
// Modified by GPT-6 on 2026-09-10
#include "STMSBridge.h"
#include "STMSLimits.hpp"
#include "ParallelRenderState.h"
#include <mutex>
#include <iostream>
#include <chrono>
#include <thread>
#include <stdexcept>
#include "verified_quality/quality.hpp"
#include "verified_quality/segments.hpp"

// Immutable while an invocation and its OpenMP workers run; protected by mutex.
static const merutilm::rff2::ParallelRenderState* activeState=nullptr;
static std::chrono::steady_clock::time_point deadline;
static bool cancelled(){return activeState->interruptRequested()||std::chrono::steady_clock::now()>=deadline;}
static void checkCancellation(){if(cancelled())throw std::runtime_error("STMS cancelled or time budget exceeded");}
#define STMS_LIBRARY
#include "prototype/locate.cpp"

namespace stms_bridge {
std::optional<Proposal> locate(merutilm::rff2::ParallelRenderState&state,
    const std::string&re,const std::string&im,double zoom,uint64_t expectedPeriod,unsigned threads){
    static std::mutex mutex;
    std::unique_lock lock(mutex,std::defer_lock);
    while(!lock.try_lock()){
        if(state.interruptRequested())return std::nullopt;
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    if(state.interruptRequested())return std::nullopt;
    if(!std::isfinite(zoom)||zoom<0||zoom>stms_limits::maxInputLogZoom)
        throw std::runtime_error("STMS input logZoom="+std::to_string(zoom)+
            " outside supported range 0.."+std::to_string(stms_limits::maxInputLogZoom));
    if(!expectedPeriod||expectedPeriod>100000000)
        throw std::runtime_error("STMS FPG period="+std::to_string(expectedPeriod)+
            " outside supported range 1..100000000");
    activeState=&state;deadline=std::chrono::steady_clock::now()+std::chrono::minutes(10);
    persistent_maps.clear();
    struct Cleanup {~Cleanup(){persistent_maps.clear();quality::checkpoints.clear();quality::probes.clear();quality::nextCheckpoint=0;activeState=nullptr;}} cleanup;
    try {
        const auto begin=std::chrono::steady_clock::now();
        const unsigned bits=stms_limits::bitsForZoom(zoom);
        ::C c(bits);
        if(mpf_set_str(c.r,re.c_str(),10)||mpf_set_str(c.i,im.c_str(),10))
            throw std::runtime_error("STMS invalid input coordinate");
        const auto solved=::solve(c,zoom,std::clamp(threads,1u,64u),0,false,expectedPeriod);
        checkCancellation();
        if(solved.period!=expectedPeriod)throw std::runtime_error("STMS period differs from RFF FPG");
        const auto verificationStart=std::chrono::steady_clock::now();
        // Independent, uniform full-precision recurrence: no jets or tapering.
        v4::precision=bits;
        inverse_samples::C candidate;
        mpfr_set_f(candidate.re.x,c.r,MPFR_RNDN);mpfr_set_f(candidate.im.x,c.i,MPFR_RNDN);
        v4::Control control{state,deadline};
        quality::checkpoints.clear();quality::probes.clear();quality::nextCheckpoint=0;
        for(const auto&block:solved.checkpoints){
            quality::Checkpoint seed{block.start,{},block.offset};
            mpfr_set_f(seed.value.re.x,block.z.r,MPFR_RNDN);mpfr_set_f(seed.value.im.x,block.z.i,MPFR_RNDN);
            quality::checkpoints.push_back(std::move(seed));
        }
        auto segmentResult=segments::evaluate(candidate,expectedPeriod,solved.logAB+30,control);
        bool refreshedCheckpoints=false;unsigned fullPrecisionCorrections=0;
        auto refreshCheckpoints=[&]{
            enclosure::Orbit refresh;size_t nextSeed=0;
            while(nextSeed<quality::checkpoints.size()){
                if((refresh.iterations&255)==0)checkCancellation();
                auto&seed=quality::checkpoints[nextSeed];
                if(refresh.iterations==seed.index){
                    mpfr_set(seed.value.re.x,refresh.value.re.x,MPFR_RNDN);
                    mpfr_set(seed.value.im.x,refresh.value.im.x,MPFR_RNDN);
                    ++nextSeed;continue;
                }
                if(refresh.iterations>seed.index||!refresh.step(candidate))
                    throw std::runtime_error("STMS checkpoint refresh failed: "+refresh.status);
            }
            refreshedCheckpoints=true;
        };
        if(segmentResult.status=="CONNECTION_FAILED"){
            refreshCheckpoints();
            segmentResult=segments::evaluate(candidate,expectedPeriod,solved.logAB+30,control);
        }
        while(segmentResult.status=="RESIDUAL_BOUND_FAILED"&&fullPrecisionCorrections<8){
            auto exact=quality::evaluateSerial(candidate,expectedPeriod,control);
            checkCancellation();
            inverse_samples::C correction;
            if(exact.status!="COMPLETE"||exact.iterations!=expectedPeriod||!inverse_samples::divide(exact.value,exact.derivative,correction))throw std::runtime_error("full precision correction failed");
            candidate=inverse_samples::sub(candidate,correction);++fullPrecisionCorrections;
            refreshCheckpoints();
            segmentResult=segments::evaluate(candidate,expectedPeriod,solved.logAB+30,control);
        }
        if(state.interruptRequested())return std::nullopt;
        if(segmentResult.status!="FULL_BOUND_PASS"||segmentResult.iterations!=expectedPeriod)
            throw std::runtime_error("STMS candidate failed bounded parallel verification: "+segmentResult.status+" "+segmentResult.failureDetail);
        Proposal result;
        result.real=v4::text(candidate.re);result.imag=v4::text(candidate.im);
        result.method=solved.method;result.period=expectedPeriod;
        result.logZoomLower=segmentResult.logZoomLower;result.logZoomUpper=segmentResult.logZoomUpper;
        result.logZoom=result.logZoomLower+(result.logZoomUpper-result.logZoomLower)*0.5;
        result.scaledResidual=v4::text(segmentResult.scaledResidual);
        result.solveSeconds=std::chrono::duration<double>(verificationStart-begin).count();
        result.verifySeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-verificationStart).count();
        std::ostringstream profile;profile.precision(17);
        profile<<"search.pilot.seconds="<<solved.pilotSeconds
            <<"\nsearch.kernel.seconds="<<solved.kernelSeconds
            <<"\nsearch.composition.seconds="<<solved.compositionSeconds
            <<"\nsearch.jetBuild.seconds="<<solved.jetBuildSeconds
            <<"\nsearch.jetReplay.seconds="<<solved.jetReplaySeconds
            <<"\nsearch.pilot.points="<<solved.pilotPoints
            <<"\nsearch.final.blocks="<<solved.checkpoints.size()
            <<"\nsearch.iterations="<<solved.iterations<<'\n';
        profile<<"verify.fullPrecisionCorrections="<<fullPrecisionCorrections<<'\n';
        profile<<"verify.checkpoints.refreshed="<<refreshedCheckpoints<<'\n';
        profile<<"verify.schedule="<<(segmentResult.longestFirst?"longest-first":"orbit-order")<<'\n';
        profile<<"verify.segment.lengths=";
        for(size_t j=0;j<solved.checkpoints.size();++j){
            if(j)profile<<',';
            const auto end=j+1<solved.checkpoints.size()?solved.checkpoints[j+1].start:expectedPeriod;
            profile<<end-solved.checkpoints[j].start;
        }
        profile<<'\n';
        for(const auto&stage:solved.stages){
            profile<<"search.stage="<<stage.iteration<<','<<stage.mode<<','<<stage.goal<<','<<stage.target
                <<','<<stage.logStep<<','<<stage.maxState<<','<<stage.stateBits<<','<<stage.derivativeBits
                <<','<<stage.kernelSeconds<<','<<stage.compositionSeconds<<'\n';
        }
        result.searchProfile=profile.str();
        return result;
    } catch(...) {
        if(state.interruptRequested())return std::nullopt;
        throw;
    }
}
}

