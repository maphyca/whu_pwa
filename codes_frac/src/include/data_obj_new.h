/*****************************************************************************
*  OpenST Basic tool library                                                 *
*  Copyright (C) 2018 Ran.Zhang  phyze96@gmail.com.                          *
*                                                                            *
*  This file is part of WHU_PWA.                                             *
*                                                                            *
*  This program is free software; you can redistribute it and/or modify      *
*  it under the terms of the GNU General Public License version 3 as         *
*  published by the Free Software Foundation.                                *
*                                                                            *
*  You should have received a copy of the GNU General Public License         *
*  along with OST. If not, see <http://www.gnu.org/licenses/>.               *
*                                                                            *
*  Unless required by applicable law or agreed to in writing, software       *
*  distributed under the License is distributed on an "AS IS" BASIS,         *
*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.  *
*  See the License for the specific language governing permissions and       *
*  limitations under the License.                                            *
*                                                                            *
*  @file     data_obj.h                                                      *
*  @brief    定义拟合程序使用的数据存取方式                                       *
*  定义初始数据存放结构，以基类DataObject为基础构造data和phsp各自的子类。             *
*                                                                            *
*  @author   Ran.Zhang                                                       *
*  @email    phyzr96@gmail.com                                               *
*  @version  1.0.0.0(版本号)                                                  *
*  @date     2018/09/13                                                       *
*  @license  GNU General Public License (GPL)                                *
*                                                                            *
*----------------------------------------------------------------------------*
*  Remark         : Description                                              *
*----------------------------------------------------------------------------*
*  Change History :                                                          *
*  <Date>     | <Version> | <Author>       | <Description>                   *
*----------------------------------------------------------------------------*
*  2018/09/13 | 1.0.0.1   | Ran.Zhang      | Create file                     *
*----------------------------------------------------------------------------*
*                                                                            *
*****************************************************************************/
#ifndef IO_H
#define IO_H
#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#endif

#ifndef STRING_H
#define STRING_H
#include<string>
#endif

#ifndef POINT_H
#define POINT_H
#include"DPFPWAPoint.h"
#endif
/**
 * @brief 存储事例的四动量，权重系数，并开放接口用于获取转换后的中间参量\n
 * 通过read_weight_file与read_data读取数据和权重，read_trans_result读取转换后的中间参量，利用Get可以获取类中所存储的数据。
 */
class DataObject
{
 public:

  /**
   * @brief 基类默认构造函数\n
   * 一般情况不会被调用
   *
   */
  explicit DataObject();

  /**
   * @brief phsp数据构造函数，无需weight文件
   * @param d_name phsp数据文件名称
   *
   */
  explicit DataObject(std::string d_name, DPFPWAPoint* pwa_point);

  /**
   * @brief data   数据构造函数，需要数据文件与weight文件
   * @param d_name data数据文件名称
   * @param w_name data权重文件名称
   */
  explicit DataObject(std::string d_name,std::string w_name, DPFPWAPoint* pwa_point);

  /**
   * @brief 基类析构函数
   *
   */
   ~DataObject();

  /**
   * @brief 友元函数，用于实现DataObject与输入输出流的交互
   * @param out  用于输出的输出流
   * @param d    被输出的DataObject对象
   */
   friend std::ostream& operator <<(std::ostream &out,DataObject &d);

   /**
    * @brief 友元类，使得计算模块可以操作数据模块数据
    *
    *
    */
   friend class CPUKernel;
   friend class CUDAKernel;

  /**
   * @brief 读取权重参数
   *
   */
  void read_weight();
  /**
   * @brief 读取初始四动量
   *
   */
  void read_data();

  /**
   * @brief 获得初始四动量
   * @return 初始四动量的指针
   */
  double* Get_mcp();

  /**
   * @brief 获得事例数
   * @return 事例数
   */
  int Get_number_of_events();

  /**
   * @brief 获得权重
   * @return 初始权重的指针
   */
  double* Get_weight();

  /**
   * @brief 计算scalar
   * @return 计算结果
   */
  double scalar(double *a1,double *a2) const;

  /**
   * @brief 计算单一事例对应中间变量
   * 
   */
  void calculate0p(int id, double* hp);

  /**
   * @brief pwa中间变量获取接口
   * 
   */
  void read_pwa();

  /**
   * @bref 存储单一事例的中间变量
   *
   */
  void store_parameters(int id,const double *hp);

  /**
   * @bref 获取中间变量指针
   * @return 中间变量集合的指针
   */
  double* Get_hpv();

  std::string  data_file_name; ///数据文件名称
  std::string  weight_file_name; ///权重文件名称

 private:
  double* _mcp; ///数据文件存放数组对应指针
  double* _weight; ///权重文件存放数组对应指针
  int _number_of_events; ///事例数
  double *_hpv; ///中间变量
  DPFPWAPoint* _dp;//point

};

std::ostream &operator <<(std::ostream &out, DataObject &d);

