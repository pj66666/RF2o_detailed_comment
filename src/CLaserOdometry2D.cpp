#include "rf2o_laser_odometry/CLaserOdometry2D.h"

namespace rf2o {

// --------------------------------------------
// CLaserOdometry2D
//---------------------------------------------

CLaserOdometry2D::CLaserOdometry2D() :
  verbose(false),
  module_initialized(false),  // 默认开始时候是没有初始化得
  first_laser_scan(true),
  last_increment_(Pose3d::Identity()),
  laser_pose_on_robot_(Pose3d::Identity()),
  laser_pose_on_robot_inv_(Pose3d::Identity()),
  laser_pose_(Pose3d::Identity()),
  laser_oldpose_(Pose3d::Identity()),
  robot_pose_(Pose3d::Identity()),
  robot_oldpose_(Pose3d::Identity())
{
  //
}

void CLaserOdometry2D::setLaserPose(const Pose3d& laser_pose)
{
  // Set laser pose on the robot

  laser_pose_on_robot_     = laser_pose;
  laser_pose_on_robot_inv_ = laser_pose_on_robot_.inverse();
}

bool CLaserOdometry2D::is_initialized()
{
  return module_initialized;  // 是否被初始化
}


/**
 * @brief 初始化
 * 
 * @param[in] scan  初始化激光数据  
 * @param[in] initial_robot_pose  初始机器人位姿 (四元数) 
 */
void CLaserOdometry2D::init(const sensor_msgs::LaserScan& scan,
                            const geometry_msgs::Pose& initial_robot_pose)
{
  //Got an initial scan laser, obtain its parametes
  ROS_INFO_COND(verbose, "[rf2o] Got first Laser Scan .... Configuring node");
                                 // std::vector<float> 论文中N  激光雷达的激光束数
  width = scan.ranges.size();    // Num of samples (size) of the scan laser 激光采样

  cols = width;						// Max resolution. Should be similar to the width parameter
  fovh = std::abs(scan.angle_max - scan.angle_min); // Horizontal Laser's FOV 水平视野角度 = angle_max-angle_min
  ctf_levels = 5;                     // Coarse-to-Fine levels 金字塔层数
  iter_irls  = 5;                      // Num iterations to solve iterative reweighted least squares
                                      // 迭代加权最小二乘问题时执行的迭代次数
  Pose3d robot_initial_pose = Pose3d::Identity(); 

  robot_initial_pose = Eigen::Quaterniond(initial_robot_pose.orientation.w,
                                          initial_robot_pose.orientation.x,
                                          initial_robot_pose.orientation.y,
                                          initial_robot_pose.orientation.z);
  // 修改了 robot_initial_pose 的平移部分，将其 x 和 y 坐标分别设为 initial_robot_pose 的位置的 x 和 y 坐
  robot_initial_pose.translation()(0) = initial_robot_pose.position.x;
  robot_initial_pose.translation()(1) = initial_robot_pose.position.y;

  ROS_INFO_STREAM_COND(verbose, "[rf2o] Setting origin at:\n"
                       << robot_initial_pose.matrix());

  // Set the initial pose 设置激光雷达的初始位姿，在机器人坐标系下
  laser_pose_    = robot_initial_pose * laser_pose_on_robot_;
  laser_oldpose_ = laser_oldpose_;

  // Init module (internal)
  //---创建了一个大小为 1xN 的矩阵 range_wf，并将所有元素的值设置为 1
  range_wf = Eigen::MatrixXf::Constant(1, width, 1);  // 扫描范围内N个点的值（距离）大小

  // Resize vectors according to levels    std::vector<Eigen::MatrixXf>
  transformations.resize(ctf_levels);   // ctf_levels = 5层金字塔
  for (unsigned int i = 0; i < ctf_levels; i++)
    transformations[i].resize(3, 3);

  // Resize pyramid
  unsigned int s, cols_i;
  // 多分辨率或粗到细策略（Coarse-to-Fine）的级别数 pyr_levels  目前，初始化时候，N=width=cols
  const unsigned int pyr_levels = std::round(std::log2(round(float(width) / float(cols)))) + ctf_levels;
  range.resize(pyr_levels);
  range_old.resize(pyr_levels);
  range_inter.resize(pyr_levels);
  xx.resize(pyr_levels);
  xx_inter.resize(pyr_levels);
  xx_old.resize(pyr_levels);
  yy.resize(pyr_levels);
  yy_inter.resize(pyr_levels);
  yy_old.resize(pyr_levels);
  range_warped.resize(pyr_levels);
  xx_warped.resize(pyr_levels);
  yy_warped.resize(pyr_levels);

  // 初始化金字塔参数
  for (unsigned int i = 0; i < pyr_levels; i++)//6
  {
    s = std::pow(2.f, int(i));    // 1 2 4 8 16 32 , 表示每一层尺度1 0.5 0.25 0.125 ... 0.0625
    cols_i = std::ceil(float(width) / float(s));  // N/s  第i层的cols

    range[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f); // 第i层range激光范围流
    range_old[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f);
    range_inter[i].resize(1, cols_i);

    xx[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f);
    xx_old[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f);

    yy[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f);
    yy_old[i] = Eigen::MatrixXf::Constant(1, cols_i, 0.f);

    xx_inter[i].resize(1, cols_i);
    yy_inter[i].resize(1, cols_i);

    if (cols_i <= cols)
    {
      range_warped[i].resize(1, cols_i);
      xx_warped[i].resize(1, cols_i);
      yy_warped[i].resize(1, cols_i);
    }
  }
  // 每个激光点束的相关参数  
  dt.resize(1, cols);   
  dtita.resize(1, cols);
  normx.resize(1, cols);  // 后面并没有归一化
  normy.resize(1, cols);
  norm_ang.resize(1, cols);
  weights.resize(1, cols);

  null    = Eigen::MatrixXi::Constant(1, cols, 0);
  cov_odo = IncrementCov::Zero();

  fps = 1.f;		//In Hz
  num_valid_range = 0;

  // Compute gaussian mask  // 权重，处理原始激光数据
  g_mask[0] = 1.f / 16.f;
  g_mask[1] = 0.25f;
  g_mask[2] = 6.f / 16.f;
  g_mask[3] = g_mask[1];
  g_mask[4] = g_mask[0];

  kai_abs_     = MatrixS31::Zero();
  kai_loc_old_ = MatrixS31::Zero();

  module_initialized = true;
  last_odom_time = ros::Time::now();
}

const CLaserOdometry2D::Pose3d& CLaserOdometry2D::getIncrement() const
{
  return last_increment_;
}

const Eigen::Matrix<float, 3, 3>& CLaserOdometry2D::getIncrementCovariance() const
{
  return cov_odo; // 协方差
}

CLaserOdometry2D::Pose3d& CLaserOdometry2D::getPose()
{
  return robot_pose_;
}

const CLaserOdometry2D::Pose3d& CLaserOdometry2D::getPose() const
{
  return robot_pose_;
}

bool CLaserOdometry2D::odometryCalculation(const sensor_msgs::LaserScan& scan)
{
  //==================================================================================
  //						DIFERENTIAL  ODOMETRY  MULTILEVEL
  //==================================================================================

  // copy laser scan to internal variable    arg: dataPtr rows = N  cols = 1  
  // Eigen::Map<const Eigen::MatrixXf> 创建了一个 Eigen 矩阵的映射对象   存储初始的激光数据
  range_wf = Eigen::Map<const Eigen::MatrixXf>(scan.ranges.data(), width, 1); 

  ros::WallTime start = ros::WallTime::now();

  createImagePyramid();   // 创建图像金字塔

  // Coarse-to-fine scheme  对每一层金字塔进行下列操作,从金字塔的最顶端开始
  for (unsigned int i=0; i<ctf_levels; i++) // 0 1 2 3 4 5
  {
    // Previous computations
    transformations[i].setIdentity();  // 首先把旋转矩阵初始化

    level = i;// 记录目前计算进度，因为上一层的计算会当作下一层的初值，由粗到细的迭代
    unsigned int s = std::pow(2.f,int(ctf_levels-(i+1)));// 从金字塔的最顶端开始，所以需要做一下处理，计算出最高层的缩放系数
    cols_i = std::ceil(float(cols)/float(s)); // 缩放后的数据长度
    image_level = ctf_levels - i + std::round(std::log2(std::round(float(width)/float(cols)))) - 1;

    //1. Perform warping
    if (i == 0) // 最高层并不会进行数据转换这部分，因此直接计算坐标
    {
      range_warped[image_level] = range[image_level];
      xx_warped[image_level]    = xx[image_level];
      yy_warped[image_level]    = yy[image_level];
    }
    else
      performWarping();

    //2. Calculate inter coords     计算range_inter,xx_inter,yy_inter
    calculateCoord();

    //3. Find null points  找到一些距离为0的点，因为在之前金字塔缩放的时候，距离小与0的不合法激光数据都被置为0了
    findNullPoints();

    //4. Compute derivatives 计算导数
    calculaterangeDerivativesSurface();

    //5. Compute normals
    //computeNormals();

    //6. Compute weights  公式14
    computeWeights();

    //7. Solve odometry
    if (num_valid_range > 3)   // 当前金字塔层有效点的数量-----至少有三个点才能求解方程，因为有三个变量
    {
      solveSystemNonLinear(); 
      //solveSystemOneLevel();    //without robust-function
    }
    else
    {
      /// @todo At initialization something
      /// isn't properly initialized so that
      /// uninitialized values get propagated
      /// from 'filterLevelSolution' first call
      /// Throughout the whole execution. Thus
      /// this 'continue' that surprisingly works.
      continue;
    }

    //8. Filter solution
    if (!filterLevelSolution()) return false;
  }

  m_runtime = ros::WallTime::now() - start;

  ROS_INFO_COND(verbose, "[rf2o] execution time (ms): %f",
                m_runtime.toSec()*double(1000));

  //Update poses
  PoseUpdate();

  return true;
}

void CLaserOdometry2D::createImagePyramid()
{
  const float max_range_dif = 0.3f;

  //Push the frames back   当前的数据备份到旧的容器   交换
  range_old.swap(range);
  xx_old.swap(xx);
  yy_old.swap(yy);

  //The number of levels of the pyramid does not match the number of levels used
  //in the odometry computation (because we sometimes want to finish with lower resolutions)

  unsigned int pyr_levels = std::round(std::log2(std::round(float(width)/float(cols)))) + ctf_levels;

  // Generate levels
  for (unsigned int i = 0; i<pyr_levels; i++)
  {
    unsigned int s = std::pow(2.f,int(i));  // 缩放倍数
    cols_i = std::ceil(float(width)/float(s));  // 每层实际分配的激光点数

    const unsigned int i_1 = i-1;

    // 1 First level -> Filter (not downsampling); 对第一层金字塔的扫描点距离进行加权处理---双边滤波器
    if (i == 0) // i = 0 第一层进行滤波
    {
      for (unsigned int u = 0; u < cols_i; u++)
      {
        const float dcenter = range_wf(u);  // range_wf is  n*1矩阵

        // 1.1 处理内部点
        if ((u>1)&&(u<cols_i-2))  // u从第三列开始到倒数第三列
        {
          if (dcenter > 0.f)      // 当前激光距离大于0,一般都满足
          {
            float sum = 0.f;      // 累加量，用于计算加权平均值   
            float weight = 0.f;   // 权重 

            for (int l=-2; l<3; l++)  // 取附近五个点
            {
              // 该点附近的五个点的距离，与当前点距离的绝对误差
              const float abs_dif = std::abs(range_wf(u+l)-dcenter);
              if (abs_dif < max_range_dif)  // max_range_dif = 0.3 > abs_dif
              { 
                // 认为越靠近该点的 权重应该最大，加权平均值
                const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                weight += aux_w;
                sum += aux_w*range_wf(u+l);
              }
            }
            // 相当于利用u周边5个点，分别给与权重，重新计算u的值
            range[i](u) = sum/weight;
          }
          else
            range[i](u) = 0.f;  // 异常点，距离赋值为0.
        }

        // 1.2 处理边界点(扫描范围边界的点)
        else
        {
          if (dcenter > 0.f)
          {
            float sum = 0.f;
            float weight = 0.f;
            // 本质上依然是加权处理，但边界的点少了几个加权  需要判断是否会越界
            for (int l=-2; l<3; l++)
            {
              const int indu = u+l;
              if ((indu>=0)&&(indu<cols_i))
              {
                const float abs_dif = std::abs(range_wf(indu)-dcenter);
                if (abs_dif < max_range_dif)
                {
                  const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                  weight += aux_w;
                  sum += aux_w*range_wf(indu);
                }
              }
            }
            range[i](u) = sum/weight;
          }
          else
            range[i](u) = 0.f;
        }
      }
    }
    // 2 处理其它层
    // 从第二层开始，做下采样，同样也是分为内部点，和边界点进行处理，原理与第一层相同，计算加权平均值
    // 不过值得注意的是，从第二层开始 使用的数据是上一层下采样后的数据，即range[i_1]，因此每次计算都会缩小2倍。

    //                              Downsampling
    //-----------------------------------------------------------------------------
    else
    {
      for (unsigned int u = 0; u < cols_i; u++)
      {
        const int u2 = 2*u;
        const float dcenter = range[i_1](u2); // 使用的数据是上一层下采样后的数据，即range[i_1]，因此每次计算都会缩小2倍

        // Inner pixels
        if ((u>0)&&(u<cols_i-1))
        {
          if (dcenter > 0.f)
          {
            float sum = 0.f;
            float weight = 0.f;

            for (int l=-2; l<3; l++)
            {
              const float abs_dif = std::abs(range[i_1](u2+l)-dcenter);
              if (abs_dif < max_range_dif)
              {
                const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                weight += aux_w;
                sum += aux_w*range[i_1](u2+l);
              }
            }
            range[i](u) = sum/weight;
          }
          else
            range[i](u) = 0.f;

        }

        //Boundary
        else
        {
          if (dcenter > 0.f)
          {
            float sum = 0.f;
            float weight = 0.f;
            const unsigned int cols_i2 = range[i_1].cols();


            for (int l=-2; l<3; l++)
            {
              const int indu = u2+l;
              if ((indu>=0)&&(indu<cols_i2))
              {
                const float abs_dif = std::abs(range[i_1](indu)-dcenter);
                if (abs_dif < max_range_dif)
                {
                  const float aux_w = g_mask[2+l]*(max_range_dif - abs_dif);
                  weight += aux_w;
                  sum += aux_w*range[i_1](indu);
                }
              }
            }
            range[i](u) = sum/weight;
          }
          else
            range[i](u) = 0.f;

        }
      }
    }

    // 3 计算每一层的激光点坐标"xy" 
    for (unsigned int u = 0; u < cols_i; u++)
    {
      if (range[i](u) > 0.f)
      {
        // 计算第u个点和x轴的夹角θ
        const float tita = -0.5*fovh + float(u)*fovh/float(cols_i-1);
        xx[i](u) = range[i](u)*std::cos(tita);
        yy[i](u) = range[i](u)*std::sin(tita);
      }
      else
      {
        xx[i](u) = 0.f;
        yy[i](u) = 0.f;
      }
    }
  }
}

void CLaserOdometry2D::calculateCoord()
{
  for (unsigned int u = 0; u < cols_i; u++)
  {
     // 如果上一帧该层的数据 或者这一帧数据 为0就不处理
    if ((range_old[image_level](u) == 0.f) || (range_warped[image_level](u) == 0.f))
    {
      range_inter[image_level](u) = 0.f;
      xx_inter[image_level](u)    = 0.f;
      yy_inter[image_level](u)    = 0.f;
    }
    else      // 否则 取两者均值
    {
      range_inter[image_level](u) = 0.5f*(range_old[image_level](u) + range_warped[image_level](u));
      xx_inter[image_level](u)    = 0.5f*(xx_old[image_level](u)    + xx_warped[image_level](u));
      yy_inter[image_level](u)    = 0.5f*(yy_old[image_level](u)    + yy_warped[image_level](u));
    }
  }
}


// 计算导数      计算公式对应于论文公式(19)
void CLaserOdometry2D::calculaterangeDerivativesSurface()
{
  // 我们的算法特别关注 range 梯度的计算。通常，使用固定的离散公式来近似扫描或图像梯度。
  // 在 range 数据的情况下，这种策略会导致对象边界处的梯度值非常高，这并不代表这些对象上的真实梯度。
  // 作为替代方案，我们使用关于环境几何形状的自适应公式。
  // 该公式使用连续观测（点）之间的 2D 距离对扫描中的前向和后向导数进行加权：

  //The gradient size ir reserved at the maximum size (at the constructor)

  //Compute connectivity

  //Defined in a different way now, without inversion
  rtita = Eigen::MatrixXf::Constant(1, cols_i, 1.f);

  for (unsigned int u = 0; u < cols_i-1; u++)
  {
    float dista = xx_inter[image_level](u+1) - xx_inter[image_level](u);// x(α) - x(α-1)
    dista *= dista; // [x(α) - x(α-1)]^2
    float distb = yy_inter[image_level](u+1) - yy_inter[image_level](u);// y(α) - y(α-1)
    distb *= distb; // [y(α) - y(α-1)]^2
    const float dist = dista + distb;   // 论文公式(19)  d(α)^2=[x(α) - x(α-1)]^2 + [y(α) - y(α-1)]^2

    if (dist  > 0.f)
      rtita(u) = std::sqrt(dist);   // 论文公式(19)  d(α)
  }

  // Spatial derivatives  论文公式(19)
  for (unsigned int u = 1; u < cols_i-1; u++)
    // 论文公式 R(α) = dtita(u)
    dtita(u) = (rtita(u-1)*(range_inter[image_level](u+1)-
                range_inter[image_level](u)) + rtita(u)*(range_inter[image_level](u) -
        range_inter[image_level](u-1)))/(rtita(u)+rtita(u-1));

  dtita(0) = dtita(1);
  dtita(cols_i-1) = dtita(cols_i-2);

  //Temporal derivative dt   Rt
  for (unsigned int u = 0; u < cols_i; u++)
    // 范围流R对事件的导数：（当前-上一次）/时间 = （当前-上一次）*频率fps
    dt(u) = fps*(range_warped[image_level](u) - range_old[image_level](u));


  //Apply median filter to the range derivatives
  //MatrixXf dtitamed = dtita, dtmed = dt;
  //vector<float> svector(3);
  //for (unsigned int u=1; u<cols_i-1; u++)
  //{
  //	svector.at(0) = dtita(u-1); svector.at(1) = dtita(u); svector.at(2) = dtita(u+1);
  //	std::sort(svector.begin(), svector.end());
  //	dtitamed(u) = svector.at(1);

  //	svector.at(0) = dt(u-1); svector.at(1) = dt(u); svector.at(2) = dt(u+1);
  //	std::sort(svector.begin(), svector.end());
  //	dtmed(u) = svector.at(1);
  //}

  //dtitamed(0) = dtitamed(1);
  //dtitamed(cols_i-1) = dtitamed(cols_i-2);
  //dtmed(0) = dtmed(1);
  //dtmed(cols_i-1) = dtmed(cols_i-2);

  //dtitamed.swap(dtita);
  //dtmed.swap(dt);
}



void CLaserOdometry2D::computeNormals()
{
  normx.setConstant(1, cols, 0.f);
  normy.setConstant(1, cols, 0.f);
  norm_ang.setConstant(1, cols, 0.f);

  const float incr_tita = fovh/float(cols_i-1);
  for (unsigned int u=0; u<cols_i; u++)
  {
    if (null(u) == 0.f)
    {
      const float tita = -0.5f*fovh + float(u)*incr_tita;
      const float alfa = -std::atan2(2.f*dtita(u), 2.f*range[image_level](u)*incr_tita);
      norm_ang(u) = tita + alfa;
      if (norm_ang(u) < -M_PI)
        norm_ang(u) += 2.f*M_PI;
      else if (norm_ang(u) < 0.f)
        norm_ang(u) += M_PI;
      else if (norm_ang(u) > M_PI)
        norm_ang(u) -= M_PI;

      normx(u) = std::cos(tita + alfa);
      normy(u) = std::sin(tita + alfa);
    }
  }
}

// 论文公式(14):预加权策略来降低范围函数是非线性甚至不可微分的那些点的残差
void CLaserOdometry2D::computeWeights()
{
  //The maximum weight size is reserved at the constructor 最大权重在构造函数中保留
  weights.setConstant(1, cols, 0.f);

  //Parameters for error_linearization
  const float kdtita = 1.f;
  const float kdt = kdtita / (fps*fps); // kdt=Δt^2
  const float k2d = 0.2f;               // Kd

  for (unsigned int u = 1; u < cols_i-1; u++)
    if (null(u) == 0)
    {
      //							Compute derivatives
      //-----------------------------------------------------------------------
      const float ini_dtita = range_old[image_level](u+1) - range_old[image_level](u-1);
      const float final_dtita = range_warped[image_level](u+1) - range_warped[image_level](u-1);

      const float dtitat = ini_dtita - final_dtita; // Rtα 
      const float dtita2 = dtita(u+1) - dtita(u-1); // Rαα

      const float w_der = kdt*(dt(u)*dt(u)) +       // kdt*(dt(u)*dt(u))=Δt^2*Rt^2
          kdtita*(dtita(u)*dtita(u)) +
          k2d*(std::abs(dtitat) + std::abs(dtita2));

      weights(u) = std::sqrt(1.f/w_der);  // 论文公式(14)  w
    }

  const float inv_max = 1.f / weights.maxCoeff();
  weights = inv_max*weights;    // 相当于让所有的权重除以最大权重，归一化
}

// 找到一些距离为0的点，因为在之前金字塔缩放的时候，距离小与0的不合法激光数据都被置为0了
void CLaserOdometry2D::findNullPoints()
{
  //Size of null matrix is set to its maximum size (constructor)
  num_valid_range = 0;

  // 
  for (unsigned int u = 1; u < cols_i-1; u++)
  {
    if (range_inter[image_level](u) == 0.f)
      null(u) = 1;
    else
    {
      num_valid_range++;
      null(u) = 0;
    }
  }
}

// Solves the system without considering any robust-function
void CLaserOdometry2D::solveSystemOneLevel()
{
  A.resize(num_valid_range, 3);
  B.resize(num_valid_range, 1);

  unsigned int cont = 0;
  const float kdtita = (cols_i-1)/fovh;

  //Fill the matrix A and the vector B
  //The order of the variables will be (vx, vy, wz)

  for (unsigned int u = 1; u < cols_i-1; u++)
    if (null(u) == 0)
    {
      // Precomputed expressions
      const float tw = weights(u);
      const float tita = -0.5*fovh + u/kdtita;

      //Fill the matrix A
      A(cont, 0) = tw*(std::cos(tita) + dtita(u)*kdtita*std::sin(tita)/range_inter[image_level](u));
      A(cont, 1) = tw*(std::sin(tita) - dtita(u)*kdtita*std::cos(tita)/range_inter[image_level](u));
      A(cont, 2) = tw*(-yy[image_level](u)*std::cos(tita) + xx[image_level](u)*std::sin(tita) - dtita(u)*kdtita);
      B(cont, 0) = tw*(-dt(u));

      cont++;
    }

  //Solve the linear system of equations using a minimum least squares method
  Eigen::MatrixXf AtA, AtB;
  AtA = A.transpose()*A;
  AtB = A.transpose()*B;
  Var = AtA.ldlt().solve(AtB);

  //Covariance matrix calculation 	Cov Order -> vx,vy,wz
  Eigen::MatrixXf res(num_valid_range,1);
  res = A*Var - B;
  // 计算协方差矩阵
  cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();

  kai_loc_level_ = Var;
}

// Solves the system by considering the Cauchy M-estimator robust-function
void CLaserOdometry2D::solveSystemNonLinear()
{
  // 残差=A*(vx, vy, wz)^T-B   公式9
  A.resize(num_valid_range, 3); Aw.resize(num_valid_range, 3);
  B.resize(num_valid_range, 1); Bw.resize(num_valid_range, 1);
  unsigned int cont = 0;
  const float kdtita = float(cols_i-1)/fovh;  // N-1/FOV  = Kα

  //Fill the matrix A and the vector B
  //The order of the variables will be (vx, vy, wz)

  for (unsigned int u = 1; u < cols_i-1; u++)
    if (null(u) == 0)// 稠密方程，计算每一个激光点的残差系数
    {
      // Precomputed expressions
      const float tw = weights(u);    // 第u个激光点的权重
      const float tita = -0.5*fovh + u/kdtita;  // 第u个激光点与x轴角度θ

      //Fill the matrix A
      A(cont, 0) = tw*(std::cos(tita) + dtita(u)*kdtita*std::sin(tita)/range_inter[image_level](u));
      A(cont, 1) = tw*(std::sin(tita) - dtita(u)*kdtita*std::cos(tita)/range_inter[image_level](u));
      A(cont, 2) = tw*(-yy[image_level](u)*std::cos(tita) + xx[image_level](u)*std::sin(tita) - dtita(u)*kdtita);
      B(cont, 0) = tw*(-dt(u));   // 注意这里是-dt(u)----残差=A*(vx, vy, wz)^T-B

      cont++;
    }

  //Solve the linear system of equations using a minimum least squares method
  Eigen::MatrixXf AtA, AtB;
  AtA = A.transpose()*A;      // 线性最小二乘法通解
  AtB = A.transpose()*B;
  Var = AtA.ldlt().solve(AtB);// Var = AtA.ldlt().solve(AtB);
                              // 使用了 Eigen 库中的 LDLT 分解方法来求解线性方程组 AtA⋅Var=AtB，
                              // 其中 Var 是包含解的向量。

  //Covariance matrix calculation 	Cov Order -> vx,vy,wz
  Eigen::MatrixXf res(num_valid_range,1);
  res = A*Var - B;      // 残差=A*(vx, vy, wz)^T-B
  //cout << endl << "max res: " << res.maxCoeff();
  //cout << endl << "min res: " << res.minCoeff();

  ////Compute the energy
  // Compute the average dt 计算平均aver_dt和平均残差aver_res
  float aver_dt = 0.f, aver_res = 0.f; unsigned int ind = 0;
  for (unsigned int u = 1; u < cols_i-1; u++)
    if (null(u) == 0)
    {
      aver_dt  += std::abs(dt(u));
      aver_res += std::abs(res(ind++));
    }
  aver_dt /= cont; aver_res /= cont;
  //    printf("\n Aver dt = %f, aver res = %f", aver_dt, aver_res);


  const float k = 10.f/aver_dt; //200
  //float energy = 0.f;
  //for (unsigned int i=0; i<res.rows(); i++)
  //	energy += log(1.f + mrpt::math::square(k*res(i)));
  //printf("\n\nEnergy(0) = %f", energy);

  //Solve iterative reweighted least squares迭代加权最小二乘法
  //===================================================================
  for (unsigned int i=1; i<=iter_irls; i++)
  {
    cont = 0;

    for (unsigned int u = 1; u < cols_i-1; u++)
      if (null(u) == 0)
      {
        // 公式12  根据每个激光点的残差，以及平均Rt，赋给每个激光点相应的权重
        // 对应一次计算中残差res大的数据，赋予较小的权重，减弱哪些outlier的影响!!!!!!
        const float res_weight = std::sqrt(1.f/(1.f + ((k*res(cont))*(k*res(cont)))));

        //Fill the matrix Aw 
        Aw(cont,0) = res_weight*A(cont,0);
        Aw(cont,1) = res_weight*A(cont,1);
        Aw(cont,2) = res_weight*A(cont,2);
        Bw(cont)   = res_weight*B(cont);
        cont++;
      }

    // Solve the linear system of equations using a minimum least squares method
    AtA = Aw.transpose()*Aw;
    AtB = Aw.transpose()*Bw;
    Var = AtA.ldlt().solve(AtB);  // 线性最小二乘法通解
    res = A*Var - B;

    ////Compute the energy
    //energy = 0.f;
    //for (unsigned int j=0; j<res.rows(); j++)
    //	energy += log(1.f + mrpt::math::square(k*res(j)));
    //printf("\nEnergy(%d) = %f", i, energy);
  }

  // res.squaredNorm() 计算了一个向量 res 的平方范数 标量
  cov_odo = (1.f/float(num_valid_range-3))*AtA.inverse()*res.squaredNorm();
  kai_loc_level_ = Var;   // 就是我们要求的速度变量

  ROS_INFO_STREAM_COND(verbose && false, "[rf2o] COV_ODO:\n" << cov_odo);
}

void CLaserOdometry2D::Reset(const Pose3d& ini_pose/*, CObservation2DRangeScan scan*/)
{
  //Set the initial pose
  laser_pose_    = ini_pose;
  laser_oldpose_ = ini_pose;

  //readLaser(scan);
  createImagePyramid();
}


/*
1. warping process  // 公式16
2. reprojected onto R1
*/
void CLaserOdometry2D::performWarping()
{
  Eigen::Matrix3f acu_trans;  // 创建了一个 3x3 的浮点型矩阵

  acu_trans.setIdentity();    // 初始化为单位矩阵
  for (unsigned int i=1; i<=level; i++)
    acu_trans = transformations[i-1]*acu_trans;   // 更新矩阵

  Eigen::MatrixXf wacu = Eigen::MatrixXf::Constant(1, cols_i, 0.f);

  range_warped[image_level].setConstant(1, cols_i, 0.f);

  const float cols_lim = float(cols_i-1); // (N-1)
  const float kdtita = cols_lim/fovh;   // (N-1)/FOV 

  for (unsigned int j = 0; j<cols_i; j++)
  {
    if (range[image_level](j) > 0.f)
    {
      // Transform point to the warped reference frame  公式16
      const float x_w = acu_trans(0,0)*xx[image_level](j) + acu_trans(0,1)*yy[image_level](j) + acu_trans(0,2);
      const float y_w = acu_trans(1,0)*xx[image_level](j) + acu_trans(1,1)*yy[image_level](j) + acu_trans(1,2);
      const float tita_w = std::atan2(y_w, x_w);          
      const float range_w = std::sqrt(x_w*x_w + y_w*y_w); // 公式17

      //Calculate warping
      const float uwarp = kdtita*(tita_w + 0.5*fovh);   // 公式18  α

      //The warped pixel (which is not integer in general) contributes to all the surrounding ones
      // 利用滤波后的α来对range流再次进行整理
      if (( uwarp >= 0.f)&&( uwarp < cols_lim))
      {
         // uwarp_l 和 uwarp_r 分别表示了插值过程中所涉及的两个整数坐标 eg  uwarp=2.3  uwarp_l=2  uwarp_r = 3
        const int uwarp_l = uwarp;
        const int uwarp_r = uwarp_l + 1;
               
        const float delta_r = float(uwarp_r) - uwarp;         // rg  delta_r = 0.7 delta_l = 0.3
        const float delta_l = uwarp - float(uwarp_l); 

        //Very close pixel 处理 变换后像素与整数坐标非常接近的情况。
        if (std::abs(std::round(uwarp) - uwarp) < 0.05f)
        {
          range_warped[image_level](round(uwarp)) += range_w;
          wacu(std::round(uwarp)) += 1.f;
        }
        else // 处理变换后像素与整数坐标不接近的情况，相当于把两个边界的点都记录了，根据距离权重
        {
          const float w_r = delta_l*delta_l;
          range_warped[image_level](uwarp_r) += w_r*range_w;
          wacu(uwarp_r) += w_r;

          const float w_l = delta_r*delta_r;
          range_warped[image_level](uwarp_l) += w_l*range_w;
          wacu(uwarp_l) += w_l;
        }
      }
    }
  }

  // Scale the averaged range and compute coordinates
  for (unsigned int u = 0; u<cols_i; u++)
  {
    if (wacu(u) > 0.f)
    {
      // 计算第u个点和x轴的夹角θ
      const float tita = -0.5f*fovh + float(u)/kdtita;
      range_warped[image_level](u) /= wacu(u);
      xx_warped[image_level](u) = range_warped[image_level](u)*std::cos(tita);
      yy_warped[image_level](u) = range_warped[image_level](u)*std::sin(tita);
    }
    else
    {
      range_warped[image_level](u) = 0.f;
      xx_warped[image_level](u) = 0.f;
      yy_warped[image_level](u) = 0.f;
    }
  }
}


// 最常见的情况是激光雷达只观测墙壁。在这种情况下，平行于墙壁的运动是不确定的，
// 因此求解器将为其提供任意解决方案（不仅是我们的方法，而且任何基于纯几何学的方法）。
// 为了缓解这个问题，我们在速度的特征空间中应用了一个低通滤波器，其工作原理如下所述

bool CLaserOdometry2D::filterLevelSolution()
{
  //		Calculate Eigenvalues and Eigenvectors
  //----使用了 Eigen 库中的 SelfAdjointEigenSolver 类来求解协方差矩阵 cov_odo 的特征值和特征向量
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> eigensolver(cov_odo);
  if (eigensolver.info() != Eigen::Success)
  {
    ROS_WARN_COND(verbose, "[rf2o] ERROR: Eigensolver couldn't find a solution. Pose is not updated");
    return false;
  }

  //First, we have to describe both the new linear and angular speeds in the "eigenvector" basis
  //-----------------在“特征向量”中描述新的线速度和角速度
  Eigen::Matrix<float,3,3> Bii;           // 特征向量
  Eigen::Matrix<float,3,1> kai_b;         // 感知速度，即Vx，Vy，w
  Bii = eigensolver.eigenvectors();       // 特征向量

  kai_b = Bii.colPivHouseholderQr().solve(kai_loc_level_);    // kai_loc_level_  vx, vy, w 最小二乘法解出的结果

  assert((kai_loc_level_).isApprox(Bii*kai_b, 1e-5) && "Ax=b has no solution." && __LINE__);

  //Second, we have to describe both the old linear and angular speeds in the "eigenvector" basis too
  //------必须在“特征向量”中描述旧的线速度和角速度
  MatrixS31 kai_loc_sub;

  //Important: we have to substract the solutions from previous levels从以前的层数中减去这个结果
  Eigen::Matrix3f acu_trans;
  acu_trans.setIdentity();
  for (unsigned int i=0; i<level; i++)
    acu_trans = transformations[i]*acu_trans;

  kai_loc_sub(0) = -fps*acu_trans(0,2);
  kai_loc_sub(1) = -fps*acu_trans(1,2);
  if (acu_trans(0,0) > 1.f)
    kai_loc_sub(2) = 0.f;
  else
  {
    kai_loc_sub(2) = -fps*std::acos(acu_trans(0,0))*rf2o::sign(acu_trans(1,0));
  }
  kai_loc_sub += kai_loc_old_;

  Eigen::Matrix<float,3,1> kai_b_old;     // t-1时刻的速度向量
  kai_b_old = Bii.colPivHouseholderQr().solve(kai_loc_sub);

  assert((kai_loc_sub).isApprox(Bii*kai_b_old, 1e-5) && "Ax=b has no solution." && __LINE__);

  // Filter speed  过滤参数：cf=Ke，Kl=df 
  const float cf = 15e3f*std::exp(-float(int(level))),
              df = 0.05f*std::exp(-float(int(level)));

  Eigen::Matrix<float,3,1> kai_b_fil;
  for (unsigned int i=0; i<3; i++)
  {
    // 公式20，这里是把特征值取出来了，所以计算了3次
    kai_b_fil(i) = (kai_b(i) + (cf*eigensolver.eigenvalues()(i,0) + df)*kai_b_old(i))/(1.f + cf*eigensolver.eigenvalues()(i,0) + df);
    //kai_b_fil_f(i,0) = (1.f*kai_b(i,0) + 0.f*kai_b_old_f(i,0))/(1.0f + 0.f);
  }

  // Transform filtered speed to local reference frame and compute transformation
  Eigen::Matrix<float, 3, 1> kai_loc_fil = Bii.inverse().colPivHouseholderQr().solve(kai_b_fil);

  assert((kai_b_fil).isApprox(Bii.inverse()*kai_loc_fil, 1e-5) && "Ax=b has no solution." && __LINE__);

  //transformation
  const float incrx = kai_loc_fil(0)/fps;   // 平移 = 速度*时间 = 速度/频率
  const float incry = kai_loc_fil(1)/fps;
  const float rot   = kai_loc_fil(2)/fps;   // 旋转角度

  transformations[level](0,0) = std::cos(rot);  // 绕z轴旋转的那个矩阵（只是旋转）
  transformations[level](0,1) = -std::sin(rot);
  transformations[level](1,0) = std::sin(rot);
  transformations[level](1,1) = std::cos(rot);
  transformations[level](0,2) = incrx;    // 平移量
  transformations[level](1,2) = incry;

  return true;
}

void CLaserOdometry2D::PoseUpdate()
{
  //  First, compute the overall transformation
  //---------------------------------------------------
  Eigen::Matrix3f acu_trans;
  acu_trans.setIdentity();

  for (unsigned int i=1; i<=ctf_levels; i++)
    acu_trans = transformations[i-1]*acu_trans; // 整体变换--每一层金字塔对应的变换矩阵不一致，逐级累乘

  //				Compute kai_loc and kai_abs
  //--------------------------------------------------------
  kai_loc_(0) = fps*acu_trans(0,2); // 局部速度 kai_loc_          acu_trans(0,2)是平移 
  kai_loc_(1) = fps*acu_trans(1,2); // 

  if (acu_trans(0,0) > 1.f)
    kai_loc_(2) = 0.f;
  else
  {               // acu_trans(0,0)角度
    kai_loc_(2) = fps*std::acos(acu_trans(0,0))*rf2o::sign(acu_trans(1,0));
  }

  //cout << endl << "Arc cos (incr tita): " << kai_loc_(2);

  float phi = rf2o::getYaw(laser_pose_.rotation());

  // 和绝对速度 kai_abs_
  kai_abs_(0) = kai_loc_(0)*std::cos(phi) - kai_loc_(1)*std::sin(phi);  // 论文公式6
  kai_abs_(1) = kai_loc_(0)*std::sin(phi) + kai_loc_(1)*std::cos(phi);  // 论文公式5
  kai_abs_(2) = kai_loc_(2);


  //						Update poses
  //-------------------------------------------------------
  laser_oldpose_ = laser_pose_;

  //  Eigen::Matrix3f aux_acu = acu_trans;
  Pose3d pose_aux_2D = Pose3d::Identity();

  pose_aux_2D = rf2o::matrixYaw(double(kai_loc_(2)/fps));
  pose_aux_2D.translation()(0) = acu_trans(0,2);
  pose_aux_2D.translation()(1) = acu_trans(1,2);

  laser_pose_ = laser_pose_ * pose_aux_2D;

  last_increment_ = pose_aux_2D;

  //				Compute kai_loc_old
  //-------------------------------------------------------
  phi = rf2o::getYaw(laser_pose_.rotation());
  kai_loc_old_(0) =  kai_abs_(0)*std::cos(phi) + kai_abs_(1)*std::sin(phi);
  kai_loc_old_(1) = -kai_abs_(0)*std::sin(phi) + kai_abs_(1)*std::cos(phi);
  kai_loc_old_(2) =  kai_abs_(2);

  ROS_INFO_COND(verbose, "[rf2o] LASERodom = [%f %f %f]",
                laser_pose_.translation()(0),
                laser_pose_.translation()(1),
                rf2o::getYaw(laser_pose_.rotation()));

  //Compose Transformations
  robot_pose_ = laser_pose_ * laser_pose_on_robot_inv_;

  ROS_INFO_COND(verbose, "BASEodom = [%f %f %f]",
                robot_pose_.translation()(0),
                robot_pose_.translation()(1),
                rf2o::getYaw(robot_pose_.rotation()));

  // Estimate linear/angular speeds (mandatory for base_local_planner)
  // last_scan -> the last scan received
  // last_odom_time -> The time of the previous scan lasser used to estimate the pose
  // 估算线性速度 lin_speed 和角速度 ang_speed
  //-------------------------------------------------------------------------------------
  double time_inc_sec = (current_scan_time - last_odom_time).toSec();
  last_odom_time = current_scan_time;
  lin_speed = acu_trans(0,2) / time_inc_sec;
  //double lin_speed = sqrt( mrpt::math::square(robot_oldpose.x()-robot_pose.x()) + mrpt::math::square(robot_oldpose.y()-robot_pose.y()) )/time_inc_sec;

  double ang_inc = rf2o::getYaw(robot_pose_.rotation()) -
      rf2o::getYaw(robot_oldpose_.rotation());

  if (ang_inc > 3.14159)
    ang_inc -= 2*3.14159;
  if (ang_inc < -3.14159)
    ang_inc += 2*3.14159;

  ang_speed = ang_inc/time_inc_sec;
  robot_oldpose_ = robot_pose_;

  //filter speeds
  /*
    last_m_lin_speeds.push_back(lin_speed);
    if (last_m_lin_speeds.size()>4)
        last_m_lin_speeds.erase(last_m_lin_speeds.begin());
    double sum = std::accumulate(last_m_lin_speeds.begin(), last_m_lin_speeds.end(), 0.0);
    lin_speed = sum / last_m_lin_speeds.size();

    last_m_ang_speeds.push_back(ang_speed);
    if (last_m_ang_speeds.size()>4)
        last_m_ang_speeds.erase(last_m_ang_speeds.begin());
    double sum2 = std::accumulate(last_m_ang_speeds.begin(), last_m_ang_speeds.end(), 0.0);
    ang_speed = sum2 / last_m_ang_speeds.size();
    */
}

} /* namespace rf2o */
