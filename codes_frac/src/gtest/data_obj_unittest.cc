#include"gtest.h"
#include"data_obj_new.h"
#ifndef CONFIG_H
#define CONFIG_H
#include "pwa_conf.h"
#endif

#ifndef STRING_H
#define STRING_H
#include <string>
#endif
TEST(DataObjectTest,phsp_weight)
{
  DataObject phsp_test((std::string)proj_path+"/gtest/gtest_data_obj.dat");
  DataObject data_test((std::string)proj_path+"/gtest/gtest_data_obj.dat",(std::string)proj_path+"/gtest/gtest_weight_file.dat");
  ASSERT_EQ(phsp_test.Get_weight()[0],1);
  ASSERT_EQ(data_test.Get_weight()[0],2);
}
