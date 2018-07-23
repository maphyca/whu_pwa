 #include "kernel.h"
#include"whu_constants_and_definitions.h"
#include<iostream>
kernel::kernel(std::vector<double *> Data, int start, int end, int device_id,int nAmps, int numbers)
{
  d_b2qbv2 = Data[b2qbv2];
}

void kernel::test()
{
  std::cout<<"GPU :"<<d_b2qbv2[0]<<std::endl;
}

