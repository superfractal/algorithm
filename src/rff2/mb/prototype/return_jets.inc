// Modified by GPT-6 on 2026-09-10
// Cached Taylor return maps in w=z^2; degree 20. See DESIGN.md and inherited attachment documentation.
constexpr int ORDER=20;
struct Poly {std::vector<C> a;explicit Poly(unsigned bits){for(int k=0;k<=ORDER;k++)a.emplace_back(bits);} };
void pswap(Poly&a,Poly&b){for(int k=0;k<=ORDER;k++){mpf_swap(a.a[k].r,b.a[k].r);mpf_swap(a.a[k].i,b.a[k].i);}}
void pzero(Poly&p){for(auto&x:p.a)zero(x);}
void pmul(Poly&o,const Poly&a,const Poly&b,Ops&op,C&tmp){for(int k=0;k<=ORDER;k++){zero(o.a[k]);for(int j=0;j<=k;j++){op.mul(tmp,a.a[j],b.a[k-j]);add(o.a[k],o.a[k],tmp);}}}
void psquare(Poly&o,const Poly&a,Ops&op,C&tmp){for(int k=0;k<=ORDER;k++){zero(o.a[k]);for(int j=0;j<=k-j;j++){op.mul(tmp,a.a[j],a.a[k-j]);if(j<k-j){mpf_mul_2exp(tmp.r,tmp.r,1);mpf_mul_2exp(tmp.i,tmp.i,1);}add(o.a[k],o.a[k],tmp);}}}
// Reuse thread-private arithmetic storage for independent coefficient products.
struct RawJetProducts {
 int threads;std::vector<std::unique_ptr<Ops>> ops;std::vector<C> scratch;
 explicit RawJetProducts(unsigned bits,int requested):threads(std::clamp(requested,1,ORDER+1)){
  scratch.reserve(threads);for(int j=0;j<threads;++j){ops.push_back(std::make_unique<Ops>(bits));scratch.emplace_back(bits);}
 }
 void multiply(Poly&o,const Poly&a,const Poly&b,bool square=false){
  if(&o==&a||&o==&b)throw std::runtime_error("parallel polynomial output aliases input");
  #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
  for(int k=0;k<=ORDER;++k){
   const auto worker=omp_get_thread_num();auto&op=*ops[worker];auto&tmp=scratch[worker];zero(o.a[k]);
   if(square){
    for(int j=0;j<=k-j;++j){op.mul(tmp,a.a[j],a.a[k-j]);if(j<k-j){mpf_mul_2exp(tmp.r,tmp.r,1);mpf_mul_2exp(tmp.i,tmp.i,1);}add(o.a[k],o.a[k],tmp);}
   }else{
    for(int j=0;j<=k;++j){op.mul(tmp,a.a[j],b.a[k-j]);add(o.a[k],o.a[k],tmp);}
   }
  }
 }
 void compose(Poly&o,const Poly&coeff,const Poly&w,Poly&buffer){
  if(&o==&coeff||&o==&w||&o==&buffer||&buffer==&coeff||&buffer==&w)
   throw std::runtime_error("parallel composition aliases input");
  pzero(o);set(o.a[0],coeff.a[ORDER]);
  for(int j=ORDER-1;j>=0;--j){
   checkCancellation();multiply(buffer,o,w);pswap(o,buffer);add(o.a[0],o.a[0],coeff.a[j]);
  }
 }
 void evaluate(Poly&square,Poly&product,const Poly&s,const Poly&b){
  if(&square==&s||&square==&b||&product==&s||&product==&b||&square==&product)
   throw std::runtime_error("raw jet output aliases input");
  #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
  for(int k=0;k<=ORDER;++k){
   const auto worker=omp_get_thread_num();auto&op=*ops[worker];auto&tmp=scratch[worker];
   zero(square.a[k]);
   for(int j=0;j<=k-j;++j){op.mul(tmp,s.a[j],s.a[k-j]);if(j<k-j){mpf_mul_2exp(tmp.r,tmp.r,1);mpf_mul_2exp(tmp.i,tmp.i,1);}add(square.a[k],square.a[k],tmp);}
   zero(product.a[k]);
   for(int j=0;j<=k;++j){op.mul(tmp,s.a[j],b.a[k-j]);add(product.a[k],product.a[k],tmp);}
  }
 }
};
void peval(Poly&o,const Poly&coeff,const Poly&w,Poly&scratch,Ops&op,C&tmp){pzero(o);set(o.a[0],coeff.a[ORDER]);for(int j=ORDER-1;j>=0;j--){pmul(scratch,o,w,op,tmp);pswap(o,scratch);add(o.a[0],o.a[0],coeff.a[j]);}}
double logadd(double a,double b){if(a<b)std::swap(a,b);if(a<-1e200)return a;return a+std::log10(1+std::pow(10.,std::max(-320.,b-a)));}
double norm10(const C&z){return la(z)/std::log2(10.);}
double power10(double x){return x<-300?0:std::pow(10.,std::min(300.,x));}
struct Map {unsigned long p;Poly s,b;double R,L;Map(unsigned /*bits*/,unsigned long pp,const Poly&ss,const Poly&bb,double r,double l):p(pp),s(ss),b(bb),R(r),L(l){};};
std::vector<Map> persistent_maps;
struct JetResult {unsigned long p=0,raw=0,jumps=0;double logAB=0,phase=0,seconds=0;};
double angle(const C&c){long er=0,ei=0;double r=mpf_get_d_2exp(&er,c.r),i=mpf_get_d_2exp(&ei,c.i);if(r==0&&i==0)return 0;if(r==0)er=ei;if(i==0)ei=er;long e=std::max(er,ei);return std::atan2(std::ldexp(i,int(std::clamp(ei-e,-2000L,0L))),std::ldexp(r,int(std::clamp(er-e,-2000L,0L))));}
JetResult jet_orbit(const C&c,double goal,double inputzoom,unsigned long target,C&oz,C&oA,C&oB,std::vector<Block>*checkpoints=nullptr){
 unsigned bits=int((goal+55)*std::log2(10.));Poly S(bits),B(bits),T(bits),U(bits),W(bits),V(bits),X(bits),Y(bits),H(bits),Q(bits),scratch(bits);Ops op(bits);C tmp(bits),tmp2(bits),w(bits),ev(bits),dv(bits),bv(bits),nextD(bits),cc(bits);set(cc,c);set(S.a[0],cc);one(S.a[1]);one(B.a[0]);
 RawJetProducts rawProducts(bits,omp_get_max_threads());
 unsigned long nextCheckpoint=std::max(256UL,target/180);
 const unsigned long checkpointSpacing=nextCheckpoint;
 if(checkpoints){checkpoints->clear();checkpoints->emplace_back(bits);}
 double logG=std::log10(2.),logH=-1e300,Rmax=0;double R=0,L=std::log10(2.)+R,minz=norm10(S.a[0]),rc=std::log10(32.)-inputzoom,ball=rc;unsigned long n=1;std::vector<Map> maps;JetResult result;auto begin=std::chrono::steady_clock::now();
  auto update_domain=[&](){R=std::min(Rmax,-std::log10(64.)-logH);L=logG+R-std::log10(1-power10(logH+R));};
  auto eval_scalar=[&](const Map&m){op.mul(w,S.a[0],S.a[0]);set(ev,m.s.a[ORDER]);set(bv,m.b.a[ORDER]);mpf_mul_ui(dv.r,m.s.a[ORDER].r,ORDER);mpf_mul_ui(dv.i,m.s.a[ORDER].i,ORDER);for(int j=ORDER-1;j>=0;--j){op.mul(ev,ev,w);add(ev,ev,m.s.a[j]);op.mul(bv,bv,w);add(bv,bv,m.b.a[j]);if(j>=1){op.mul(dv,dv,w);mpf_mul_ui(tmp.r,m.s.a[j].r,j);mpf_mul_ui(tmp.i,m.s.a[j].i,j);add(dv,dv,tmp);}}op.mul(tmp,dv,S.a[0]);mpf_mul_2exp(tmp.r,tmp.r,1);mpf_mul_2exp(tmp.i,tmp.i,1);op.mul(nextD,tmp,B.a[0]);add(nextD,nextD,bv);};
 while(n<(target?target:100000000)){if(cancelled())return {};
  if(checkpoints&&n>=nextCheckpoint){
   checkpoints->emplace_back(bits);auto&seed=checkpoints->back();seed.start=n;
   set(seed.z,S.a[0]);set(seed.fb,B.a[0]);seed.offset=la(S.a[1]);nextCheckpoint=n+checkpointSpacing;
  }
  double lz=norm10(S.a[0]);bool jumped=false;if(!target)ball=logH+rc>=0?1e300:logG+rc-std::log10(1-power10(logH+rc));
  if(lz<-3){for(int mi=int(maps.size())-1;mi>=0;--mi){const Map&m=maps[mi];if(target &&n+m.p>target)continue;double lq=2*lz-m.R;if(lq>std::log10(.1))continue;double q=power10(lq);if(!target &&(2*logadd(lz,ball)>m.R-std::log10(4.) ||rc>m.R-std::log10(4.)))continue;
    eval_scalar(m);double le=m.L+(ORDER+1)*lq-std::log10(1-q);if(le-norm10(nextD)>-goal-12)continue;
    const double logAu=norm10(dv)+1e-12, logBv=norm10(bv)+1e-12;
    const double oldG=logG,oldH=logH;
    const double Gu=logadd(std::log10(2.)+lz+oldG,0.);
    const double Gout=logadd(std::log10(2.)+lz+logAu+oldG,logBv);
    const double Hu=logadd(oldH,2*oldG-Gu);
    const double Hl=logadd(oldH,logAu+2*oldG-Gout);
    const double J=std::log10(2.)-m.R-std::log10(1-q);
    const double Cb=std::log10(8.)+m.L-2*m.R-2*std::log10(1-q);
    const double Hbase=std::max(Hl,logadd(Hu,J+Gu));
    const double Hout=logadd(Hbase,Cb+2*Gu-Gout);
    const double ld=m.R/2+std::log10(1-power10(lz-m.R/2));
    const double allow=ld-logadd(oldG,oldH+ld);
    Rmax=std::min({Rmax,allow,m.R});logG=Gout;logH=Hout;
    rawProducts.multiply(W,S,S,true);rawProducts.compose(V,m.s,W,scratch);rawProducts.compose(X,m.b,W,scratch);
    pzero(H);mpf_mul_ui(H.a[0].r,m.s.a[ORDER].r,ORDER);mpf_mul_ui(H.a[0].i,m.s.a[ORDER].i,ORDER);for(int j=ORDER-1;j>=1;--j){rawProducts.multiply(scratch,H,W);pswap(H,scratch);mpf_mul_ui(tmp.r,m.s.a[j].r,j);mpf_mul_ui(tmp.i,m.s.a[j].i,j);add(H.a[0],H.a[0],tmp);}rawProducts.multiply(Y,H,S);for(auto&x:Y.a){mpf_mul_2exp(x.r,x.r,1);mpf_mul_2exp(x.i,x.i,1);}rawProducts.multiply(Q,Y,B);for(int j=0;j<=ORDER;j++)add(Q.a[j],Q.a[j],X.a[j]);pswap(S,V);pswap(B,Q);update_domain();n+=m.p;result.jumps++;jumped=true;break;
  }}
  if(!jumped){const double e=lz-std::log10(64.);Rmax=std::min(Rmax,e-logadd(logG,e+logH));const double Gn=logadd(std::log10(2.)+lz+logG,0.);logH=logadd(logH,2*logG-Gn);logG=Gn;update_domain();rawProducts.evaluate(T,U,S,B);add(T.a[0],T.a[0],cc);for(auto&x:U.a){mpf_mul_2exp(x.r,x.r,1);mpf_mul_2exp(x.i,x.i,1);}mpf_add_ui(U.a[0].r,U.a[0].r,1);pswap(S,T);pswap(B,U);n++;result.raw++;}
  lz=norm10(S.a[0]);if(lz<minz){minz=lz;if(n>1)maps.emplace_back(bits,n,S,B,R,L);if(std::getenv("JET_TRACE"))std::cerr<<"record "<<n<<" z="<<lz<<" R="<<R<<" L="<<L<<" G="<<logG<<" H="<<logH<<" raw="<<result.raw<<" jumps="<<result.jumps<<"\n";}
  if(!target && lz-norm10(B.a[0])<rc) break;
  if(lz>16)throw std::runtime_error("jet escaped");
  if(result.raw%100000==0&&result.raw>0&&std::getenv("JET_TRACE"))std::cerr<<"progress "<<n<<" raw="<<result.raw<<" jumps="<<result.jumps<<"\n";
 }
 if(!target && n>=100000000 && norm10(S.a[0])-norm10(B.a[0])>=rc) throw std::runtime_error("jet period search limit");
 persistent_maps=std::move(maps);set(oz,S.a[0]);set(oA,S.a[1]);set(oB,B.a[0]);result.p=n;result.logAB=norm10(oA)+norm10(oB);result.phase=angle(oA)+angle(oB);result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();return result;
}
JetResult cached_orbit(const C&base,const C&c,double goal,unsigned long target,C&oz,C&oA,C&oB,std::vector<Block>*checkpoints=nullptr){
 unsigned bits=int((goal+70)*std::log2(10.));C z(bits),A(bits),D(bits),w(bits),v(bits),dv(bits),b(bits),db(bits),dc(bits),t(bits),u(bits),newD(bits),newA(bits);Ops op(bits);sub(dc,c,base);double ldc=norm10(dc);unsigned long n=0;JetResult result;auto begin=std::chrono::steady_clock::now();one(A);
 unsigned long nextCheckpoint=0;const unsigned long spacing=std::max(256UL,target/180);
 if(checkpoints)checkpoints->clear();
 while(n<target){if(cancelled())return {};
  if(checkpoints&&n>=nextCheckpoint){
   checkpoints->emplace_back(bits);auto&seed=checkpoints->back();seed.start=n;set(seed.z,z);seed.offset=la(A);nextCheckpoint=n+spacing;
  }
  double lz=norm10(z);bool jumped=false;
  for(int k=int(persistent_maps.size())-1;k>=0;--k){const auto&m=persistent_maps[k];if(n+m.p>target)continue;double lq=2*lz-m.R,ls=ldc-m.R;if(lq>-.8||ls>-.8)continue;double q=power10(lq),s=power10(ls);op.mul(w,z,z);set(v,m.s.a[ORDER]);set(b,m.b.a[ORDER]);zero(dv);zero(db);
   for(int j=ORDER-1;j>=0;--j){op.mul(dv,dv,w);add(dv,dv,v);op.mul(db,db,w);add(db,db,b);op.mul(v,v,w);add(v,v,m.s.a[j]);op.mul(b,b,w);add(b,b,m.b.a[j]);}
   op.mul(t,db,dc);add(dv,dv,t);op.mul(t,b,dc);add(v,v,t);op.mul(t,dv,z);mpf_mul_2exp(t.r,t.r,1);mpf_mul_2exp(t.i,t.i,1);op.mul(newD,t,D);add(newD,newD,b);if(n==0)set(newA,dv);else op.mul(newA,t,A);
   double le=m.L+logadd((ORDER+1)*lq,2*ls)-std::log10((1-q)*(1-s));if(le-norm10(newD)>-goal-8)continue;
   set(z,v);set(D,newD);set(A,newA);n+=m.p;result.jumps++;jumped=true;break;
  }
  if(!jumped){if(n){op.mul(A,A,z);mpf_mul_2exp(A.r,A.r,1);mpf_mul_2exp(A.i,A.i,1);}op.mul(D,D,z);mpf_mul_2exp(D.r,D.r,1);mpf_mul_2exp(D.i,D.i,1);mpf_add_ui(D.r,D.r,1);op.mul(z,z,z);add(z,z,c);n++;result.raw++;}
  if(norm10(z)>16)throw std::runtime_error("cache escape");
 }
 set(oz,z);set(oA,A);set(oB,D);result.p=n;result.logAB=norm10(A)+norm10(D);result.phase=angle(A)+angle(D);result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();return result;
}


