/*****************************************************************************
 * *  PWA GPU part code,head file                                               *
 * *  Copyright (C) 2018 Arthur.Chang  phyzr96@gmail.com.                       *
 * *                                                                            *
 * *  This file is part of PWA.                                                 *
 * *                                                                            *
 * *  This program is free software; you can redistribute it and/or modify      *
 * *  it under the terms of the GNU General Public License version 3 as         *
 * *  published by the Free Software Foundation.                                *
 * *                                                                            *
 * *  You should have received a copy of the GNU General Public License         *
 * *  along with OST. If not, see <http://www.gnu.org/licenses/>.               *
 * *                                                                            *
 * *  Unless required by applicable law or agreed to in writing, software       *
 * *  distributed under the License is distributed on an "AS IS" BASIS,         *
 * *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  *
 * *  See the License for the specific language governing permissions and       *
 * *  limitations under the License.                                            *
 * *                                                                            *
 * *  @file     kernel_calEva.h                                                 *
 * *  @brief    kernel_calEva的头文件                                           *
 * *  声明kernel_calEva各部分，以及cuda_kernel类的定义                          *                                        *
 * *                                                                            *
 * *  @author   Arthur.Chang                                                    *
 * *  @email    phyzr96@gmail.com                                               *
 * *  @version  1.0.0.0(版本号)                                                 *
 * *  @date     2018.4.15                                                       *
 * *  @license  GNU General Public License (GPL)                                *
 * *                                                                            *
 * *----------------------------------------------------------------------------*
 * *  Remark         : Description                                              *
 * *----------------------------------------------------------------------------*
 * *  Change History :                                                          *
 * *  <Date>     | <Version> | <Author>       | <Description>                   *
 * *----------------------------------------------------------------------------*
 * *  2018/04/15 | 1.0.0.0   | Arthur.Chang   | Create file                     *
 * *----------------------------------------------------------------------------*
 * *                                                                            *
 * *****************************************************************************/

#ifndef _KERNEL_CALEVA_H
#define _KERNEL_CALEVA_H
#include "/opt/cuda/include/cuda_runtime.h"
#include "cuComplex.h"
#include <iostream>
#include "cu_PWA_PARAS.h"
#include <vector>
#include <fstream>
#include <math.h>
#include "cu_DPFPropogator.h"
#include <assert.h>
#include <vector>
#include "MultDevice.h"
//double *d_float_pp=NULL;
//#include "cu_PWA_PARAS.h"
    //int initialize_data(std::vector<PWA_PARAS>&, DataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //int data_distribution(DataPointers&, CudaDataPointers&); // 把vector数据和指针数据对应起来，并copy到gpu里面去
    //double calEva(const PWA_PARAS &pp, int idp);
    //double kernel_calEva(const PWA_PARAS &pp,int idp);
    //void func(DataPointers& cpu_data_pointers);
/**
 * @brief 本类是PDF类的成员，包含GPU相关各参数，避免PDF切换带来的各种问题
 * 该类主要包含以d_开头的GPU量在GPU上的全局内存地址，以h_开头的CPU量在CPU端的地址或值
 * 以及用于GPU空间申请，GPU计算的各个成员函数 
*/
class cuda_kernel {
    public:

    cuda_kernel() {};
    double *d_fx[DEVICE_NUM];
    int *d_parameter[DEVICE_NUM];
    double *d_paraList[DEVICE_NUM];
    //double2 * d_complex_para[DEVICE_NUM];
    double *d_mlk[DEVICE_NUM];
    double *h_mlk_pt[DEVICE_NUM];
    double *d_fx_store[DEVICE_NUM];
    double h_fx_store[DEVICE_NUM];

    void cu_malloc_h_pp(double *,double *&,int,int);
    /**
    *  @brief               fx计算以及存储的具体实现，该部分由CPU调用
    *  @param d_float_pp    fx计算所需的pp值的GPU地址
    *  @param h_parameter   fx计算所需参数CPU端地址，用于后续传输进GPU
    *  @param h_paraList    fx计算所需具体值CPU端地址，用于后续传输进GPU
    *  @param para_size     paraList的大小
    *  @param h_fx          fx计算结果的CPU端存储地址
    *  @param h_mlk         mlk求和结果的CPU端存储地址
    *  @param end           事例的结束位置
    *  @param begin         事例的起始位置，实际上该数值必然为0，并无具体作用
    *  @param anaint        外部传入的anaIntegral地址，用于存放fx的求和结果
    *  @note                由于精度问题，以及为了方便后续修改，目然以实现GPU求和功能，目前fx求和部分仍在CPU端完成，大量时间用于fx传输以及求和。
    */

    int host_store_fx(std::vector<double *>,int *,double *,int , double *,double * ,int ,int ,double*);
    //int malloc_mem(int end, int begin, int para_size, int *h_parameter)
    //
    int warp_malloc_mem(int, int, int, int *);
    int malloc_mem(int, int, int, int*);
};
/**
*  @brief               计算fx与mlk的GPU部分
*  @param pp            用于简化计算的pwa_paras部分
*  @param parameter     外部传入的参数，包括paraList中各参数起始位置，事例数等部分
*  @param d_paraList    各个事例属性的具体值，结合parameter定位即可获得所需值
*  @param d_mlk         用于存储计算得到的mlk值，目前通过在共享内存进行原子求和，再以块为单位进行全局内存原子求和的机制，但原子求和存在大量排队时间浪费
*  @param idp           对应线程编号，由于采用原子求和，目前无作用
*  @param offset        偏移量，由于数据集并非均等分割，为了分割Nmc与Nmc_data部分，通过offset控制数值存放位置，该值由global部分判断并传入，Nmc部分offset为0，Nmc_data部
分为1
*  @return              返回计算所得fx值，并利用传入的指针修改对应mlk的值
*  @note                fx求和亦可在此处通过共享内存完成，此时将无需传输fx值， 但目前仍采取传输进入全局内存，并用专用global函数fx_sum完成，此部分依然有优化空间
*  @todo                由于所有线程采用同一fCP值，后续可以将fCP求值部分移出calEva部
分
*/

    __device__ double calEva(const cu_PWA_PARAS *pp, const int * parameter  ,const double * d_paraList,double *d_mlk,int idp,int offset) ;
    /**
    *  @brief               读取计算所需参数，由于线程存在冗余，还需判断是否使用该线程进
    行计算
    *  @warning             由于采用将pp读入寄存器的方式，而cuda单线程寄存器上限255，会有部分寄存器溢出至缓存与全局内存
    *  @todo                采用内存对齐的方式可以减少参数读入的时间，目前尚未对齐
    *  @bug                  paraList与mlk采用动态共享内存的方式存储，当paraList与mlk结构较大时可能存在共享内存不足初始化失败的情况，此时核函数将会启动失败
    *  @param float_pp      一维形式的pp参数，会重新以cu_PWA_PARAS的结构存入寄存器
    *  @param parameter     外部传入的各参数对应全局内存地址
    *  @param d_paraList    各事例具体值的列表对应全局内存地址
    *  @param para_size     paraList的大小
    *  @param d_fx          fx对应全局内存储存地址
    *  @param d_mlk         d_mlk对应全局内存地址，由于求和已写入该部分，mlk分为两个nAmps的部分，分别用于存储Nmc与Nmc_data部分的求和结果
    *  @param end           计算的结束位置
    *  @param begin         计算的开始位置
    */

    __global__ void kernel_store_fx(const double * float_pp,const int *parameter,const double *d_paraList,int para_size,double * d_fx,double *d_mlk,int end,int begin);
/**
*  @brief               在每次操作前将全军内存值初始化为0，用于后续求和
*  @param d_mlk         mlk求和对应全局内存地址，需进行nAmps次操作
*  @param d_fx_store    fx求和对应全局内存地址，只需操作一次
*  @param num           即nAmps，mlk需要初始化的次数
*/

    __global__ void reset_mlk (double *d_mlk,double *d_fx_store,int num);
    /**
    *  @brief               从全局内存读取fx,并利用原子操作进行求和
    *  @note                原子求和以及全局内存操作存在排队以及访问时间浪费，
    可后续以合并进入kernel_store中通过共享内存操作完成
    *  @param d_fx          fx的全局内存地址
    *  @param d_fx_store    fx和的全局内存地址
    *  @param num           由于只需对Nmc部分求和，该参数用于控制fx求和的范围，由外部判断并传入
    */

    __global__ void fx_sum(double *d_fx,double *d_fx_store,int num);
#endif

