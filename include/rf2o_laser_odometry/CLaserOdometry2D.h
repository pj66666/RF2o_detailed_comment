#ifndef CLaserOdometry2D_H
#define CLaserOdometry2D_H

// std header
#include <iostream>
#include <fstream>
#include <numeric>

// ROS headers
#include <ros/ros.h>
#include <nav_msgs/Odometry.h>
#include <sensor_msgs/LaserScan.h>

// Eigen headers
#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>

namespace rf2o {

template <typename T> // 小于0返回-1；大于0返回1
inline T sign(const T x) { return x<T(0) ? -1:1; }

template <typename Derived>   // 从旋转矩阵中获取绕 Z 轴的偏航角
inline typename Eigen::MatrixBase<Derived>::Scalar
getYaw(const Eigen::MatrixBase<Derived>& r) //  Eigen::MatrixBase提供了一种方便和通用的方式来操作各种不同类型的矩阵和向量
{// 旋转矩阵的第二列的x分量和第一列的x分量，这两者恰好对应于新坐标系下的X轴和Y轴的投影,表示了新坐标系相对于旧坐标系绕Z轴的旋转角度
  return std::atan2( r(1, 0), r(0, 0) );
}


template<typename T> // 绕定轴rpy旋转   R12  这里可以表示把1系坐标轴转到2系坐标轴，或者2系中点转到1系
inline Eigen::Matrix<T, 3, 3> matrixRollPitchYaw(const T roll,
                                                 const T pitch,
                                                 const T yaw)
{
  const Eigen::AngleAxis<T> ax = Eigen::AngleAxis<T>(roll,  Eigen::Matrix<T, 3, 1>::UnitX());
  const Eigen::AngleAxis<T> ay = Eigen::AngleAxis<T>(pitch, Eigen::Matrix<T, 3, 1>::UnitY());
  const Eigen::AngleAxis<T> az = Eigen::AngleAxis<T>(yaw,   Eigen::Matrix<T, 3, 1>::UnitZ());

  return (az * ay * ax).toRotationMatrix().matrix();
}

template<typename T>// 因为这里使用的是二维，只是绕z轴旋转,   (0, 0, yaw)
inline Eigen::Matrix<T, 3, 3> matrixYaw(const T yaw)
{
  return matrixRollPitchYaw<T>(0, 0, yaw);
}

class CLaserOdometry2D
{
public:

  using Scalar = float;

  using Pose2d = Eigen::Isometry2d;   // 变换矩阵T(3*3)
  using Pose3d = Eigen::Isometry3d;   // 变换矩阵T(4*4)
  using MatrixS31 = Eigen::Matrix<Scalar, 3, 1>;
  using IncrementCov = Eigen::Matrix<Scalar, 3, 3>;

  CLaserOdometry2D();
  virtual ~CLaserOdometry2D() = default;

  void init(const sensor_msgs::LaserScan& scan,
            const geometry_msgs::Pose& initial_robot_pose);

  bool is_initialized();

  bool odometryCalculation(const sensor_msgs::LaserScan& scan);

  void setLaserPose(const Pose3d& laser_pose);

  const Pose3d& getIncrement() const;

  const IncrementCov& getIncrementCovariance() const;

  Pose3d& getPose();
  const Pose3d& getPose() const;

protected:

  bool verbose, module_initialized, first_laser_scan; // 是否输出标志。模块是否初始化，是否第一次激光数据

  // Internal Data
  std::vector<Eigen::MatrixXf> range;       // 存储range_wf滤波后(加权)的每一层金字塔的扫描数据
  std::vector<Eigen::MatrixXf> range_old;   // 上一次测量数据
  std::vector<Eigen::MatrixXf> range_inter; // 存储一些插值后的距离测量数据=0.5*warped+0.5*old
  std::vector<Eigen::MatrixXf> range_warped;// 经过wrap后的数据(双边滤波)

  std::vector<Eigen::MatrixXf> xx;          // 记录每一层的激光点的x坐标(加权后)
  std::vector<Eigen::MatrixXf> xx_inter;
  std::vector<Eigen::MatrixXf> xx_old;
  std::vector<Eigen::MatrixXf> xx_warped;

  std::vector<Eigen::MatrixXf> yy;          // 记录每一层的激光点的y坐标(加权后)
  std::vector<Eigen::MatrixXf> yy_inter;
  std::vector<Eigen::MatrixXf> yy_old;
  std::vector<Eigen::MatrixXf> yy_warped;
  std::vector<Eigen::MatrixXf> transformations;

  Eigen::MatrixXf range_wf;   // 初始激光扫描的数据
  Eigen::MatrixXf dtita;      // 论文公式(19) 导数梯度 R(α)
  Eigen::MatrixXf dt;         // Rt
  Eigen::MatrixXf rtita;      // 论文公式(19)  d(α)
  Eigen::MatrixXf normx, normy, norm_ang;
  Eigen::MatrixXf weights;  // 预加权策略：对应每一个数据的权重  公式14
  Eigen::MatrixXi null;   // 记录当前金字塔层的第i个激光点距离是否小于等于0，即异常，异常为1，正常为0

  Eigen::MatrixXf A,Aw;   // 残差=A*(vx, vy, wz)^T-B   第一个就是线性最小二乘，带w表示加权最小二乘  公式9
  Eigen::MatrixXf B,Bw;   // 残差=A*(vx, vy, wz)^T-B  公式9

  MatrixS31 Var;	//3 unknowns: vx, vy, w
  IncrementCov cov_odo; // 计算协方差矩阵

  //std::string LaserVarName;				//Name of the topic containing the scan lasers \laser_scan
  float fps;								//In Hz
  float fovh;								//Horizontal FOV
  unsigned int cols;        // 采样的激光点数量
  unsigned int cols_i;      // 第i层金字塔采样的激光点数量
  unsigned int width;
  unsigned int ctf_levels;  // 金字塔层数
  unsigned int image_level, level;    // 
  unsigned int num_valid_range;       // 有效点的数量
  unsigned int iter_irls;   // 迭代加权最小二乘问题时执行的迭代次数
  float g_mask[5];          // 权重，处理原始激光数据----处理金字塔数据

  double lin_speed, ang_speed;  // 线性速度 lin_speed 和角速度 ang_speed

  ros::WallDuration	m_runtime;
  ros::Time last_odom_time, current_scan_time;

  MatrixS31 kai_abs_;       //储存线速度, 线速度与角度的乘积 和角度这三个值
  MatrixS31 kai_loc_;       // vx, vy, w
  MatrixS31 kai_loc_old_;   // vx, vy, w
  MatrixS31 kai_loc_level_; // kai_loc_level_ = Var;   vx, vy, w   最小二乘法解出的结果

  Pose3d last_increment_;
  Pose3d laser_pose_on_robot_;
  Pose3d laser_pose_on_robot_inv_;
  Pose3d laser_pose_;
  Pose3d laser_oldpose_;
  Pose3d robot_pose_;
  Pose3d robot_oldpose_;

  bool test;
  std::vector<double> last_m_lin_speeds;
  std::vector<double> last_m_ang_speeds;

  // Methods
  void createImagePyramid();
  void calculateCoord();
  void performWarping();
  void calculaterangeDerivativesSurface();
  void computeNormals();
  void computeWeights();
  void findNullPoints();
  void solveSystemOneLevel();
  void solveSystemNonLinear();
  bool filterLevelSolution();
  void PoseUpdate();
  void Reset(const Pose3d& ini_pose/*, CObservation2DRangeScan scan*/);
};

} /* namespace rf2o */

#endif // CLaserOdometry2D_H
