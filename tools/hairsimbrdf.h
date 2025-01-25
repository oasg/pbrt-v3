#include<vector>
#include<mutex>
#include <iostream>
#include <memory>
#include "def.h"

struct RGB{
  Float r;
  Float g;
  Float b;
};
class hairSimBrdf{
  public:
    hairSimBrdf(const char* file);
    ~hairSimBrdf(){}
    RGB getReflect(Float it, Float ot);
    std::vector<std::vector<RGB>> m_data;
};
class SingBrdf{
  public:
  ~SingBrdf(){
    std::cout<<"sim Brdf destructor!"<<std::endl;
  }
  SingBrdf() = delete;
  SingBrdf& operator=(const SingBrdf&)= delete;
  static std::shared_ptr<hairSimBrdf> get_Instance(){
    if(m_instance_ptr==nullptr){
      std::lock_guard<std::mutex> lk(m_mutex);
      if(m_instance_ptr==nullptr){
        m_instance_ptr = std::shared_ptr<hairSimBrdf>(
          // new hairSimBrdf("../../table/HairDamagedModel/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/")); 
          new hairSimBrdf("../../table/HairModel/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));   
          // new hairSimBrdf("../../table/HairDamagedModelLargeDis/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));   
          // new hairSimBrdf("../../table/HairMultilayerPerlinModel/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));
          // new hairSimBrdf("../../table/HairDamagedModelLargeDisPerlin/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));
          // new hairSimBrdf("../../table/HairRepairedModelLargeDis/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));
          // new hairSimBrdf("../../table/HairRepairedModelLargeDisPerlin/incidenceLayer/(10nm,800cell)/TM/Ns/Reflection/"));
      }
    }
    return m_instance_ptr;
  }
  private:
    static std::shared_ptr<hairSimBrdf> m_instance_ptr;
    static std::mutex m_mutex;
};