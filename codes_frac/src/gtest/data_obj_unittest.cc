#include"gtest.h"
#include"data_obj_new.h"
#include"DPFPWAPoint.h"
#ifndef CONFIG_H
#define CONFIG_H
#include "pwa_conf.h"
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif
#include"whu_constants_and_definitions.h"
TEST(DataObjectTest,phsp_weight)
{
  DPFPWAPoint *pwa_point_phipp = new DPFPWAPoint(mka, mka, mpi, mpi, mpsip);
  DPFPWAPoint *pwa_point_phikk = new DPFPWAPoint(mka, mka, mka, mka, mpsip);
  DataObject phsp_test((std::string)proj_path+"/gtest/gtest_data_obj.dat",pwa_point_phikk);
  DataObject data_test((std::string)proj_path+"/gtest/gtest_data_obj.dat",(std::string)proj_path+"/gtest/gtest_weight_file.dat",pwa_point_phipp);
  ASSERT_EQ(phsp_test.Get_weight()[0],1);
  ASSERT_EQ(data_test.Get_weight()[0],2);
}
