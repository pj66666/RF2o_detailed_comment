/** ****************************************************************************************
*  This node presents a fast and precise method to estimate the planar motion of a lidar
*  from consecutive range scans. It is very useful for the estimation of the robot odometry from
*  2D laser range measurements.
*  This module is developed for mobile robots with innacurate or inexistent built-in odometry.
*  It allows the estimation of a precise odometry with low computational cost.
*  For more information, please refer to:
*
*  Planar Odometry from a Radial Laser Scanner. A Range Flow-based Approach. ICRA'16.
*  Available at: http://mapir.isa.uma.es/mapirwebsite/index.php/mapir-downloads/papers/217
*
* Maintainer: Javier G. Monroy
* MAPIR group: http://mapir.isa.uma.es/
*
* Modifications: Jeremie Deray
******************************************************************************************** */

#include "rf2o_laser_odometry/CLaserOdometry2D.h"

#include <tf/transform_broadcaster.h>
#include <tf/transform_listener.h>

namespace rf2o {

class CLaserOdometry2DNode : CLaserOdometry2D
{
public:

  CLaserOdometry2DNode();
  ~CLaserOdometry2DNode() = default;

  void process(const ros::TimerEvent &);
  void publish();

  bool setLaserPoseFromTf();

public:

  bool publish_tf, new_scan_available;  // new_scan_available判断有无新数据发布的标志位，回调函数中调用

  double freq;

  std::string         laser_scan_topic;
  std::string         odom_topic;
  std::string         base_frame_id;
  std::string         odom_frame_id;
  std::string         init_pose_from_topic;

  ros::NodeHandle             n;
  sensor_msgs::LaserScan      last_scan;
  bool                        GT_pose_initialized;
  tf::TransformListener       tf_listener;          // Do not put inside the callback
  tf::TransformBroadcaster    odom_broadcaster;
  nav_msgs::Odometry          initial_robot_pose;

  //Subscriptions & Publishers
  ros::Subscriber laser_sub, initPose_sub;
  ros::Publisher odom_pub;

  bool scan_available();

  //CallBacks
  void LaserCallBack(const sensor_msgs::LaserScan::ConstPtr& new_scan);
  void initPoseCallBack(const nav_msgs::Odometry::ConstPtr& new_initPose);
};

CLaserOdometry2DNode::CLaserOdometry2DNode() :
  CLaserOdometry2D()
{
  ROS_INFO("Initializing RF2O node...");

  //Read Parameters
  //----------------
  ros::NodeHandle pn("~");
  // 从ROS参数服务器中获取名为laser_scan_topic的参数值，将其转化为字符串，并存储在变量laser_scan_topic中。
  // 如果参数不存在，就将默认值"/laser_scan"赋给laser_scan_topic。
  pn.param<std::string>("laser_scan_topic",laser_scan_topic,"/laser_scan");
  pn.param<std::string>("odom_topic", odom_topic, "/odom_rf2o");
  pn.param<std::string>("base_frame_id", base_frame_id, "/base_link");
  pn.param<std::string>("odom_frame_id", odom_frame_id, "/odom");
  pn.param<bool>("publish_tf", publish_tf, true);
  pn.param<std::string>("init_pose_from_topic", init_pose_from_topic, "/base_pose_ground_truth");
  pn.param<double>("freq",freq,10.0);
  pn.param<bool>("verbose", verbose, true);

  //Publishers and Subscribers
  // 发布算法处理后的里程计数据
  odom_pub  = pn.advertise<nav_msgs::Odometry>(odom_topic, 5);
  // 订阅来自laser_scan_topic激光数据
  laser_sub = n.subscribe<sensor_msgs::LaserScan>(laser_scan_topic,1,&CLaserOdometry2DNode::LaserCallBack,this);

  // init pose??  指定任何ROS话题作为初始姿态信息的来源
  if (init_pose_from_topic != "")
  {
    initPose_sub = n.subscribe<nav_msgs::Odometry>(init_pose_from_topic,1,&CLaserOdometry2DNode::initPoseCallBack,this);
    GT_pose_initialized  = false;
  }
  else
  {
    // 这里很关键，不指定任何ROS话题作为初始姿态信息的来源，使用固定点（0，0）. 否则无法进入回调函数
    GT_pose_initialized = true;
    initial_robot_pose.pose.pose.position.x = 0;
    initial_robot_pose.pose.pose.position.y = 0;
    initial_robot_pose.pose.pose.position.z = 0;
    initial_robot_pose.pose.pose.orientation.w = 0;
    initial_robot_pose.pose.pose.orientation.x = 0;
    initial_robot_pose.pose.pose.orientation.y = 0;
    initial_robot_pose.pose.pose.orientation.z = 0;
  }

  setLaserPoseFromTf();

  //Init variables
  module_initialized = false;
  first_laser_scan   = true;

  ROS_INFO_STREAM("Listening laser scan from topic: " << laser_sub.getTopic());
}

bool CLaserOdometry2DNode::setLaserPoseFromTf()
{
  bool retrieved = false;   // 表示尚未成功获取激光雷达的位姿信息

  // Set laser pose on the robot (through tF)
  // This allow estimation of the odometry with respect to the robot base reference system.
  tf::StampedTransform transform;
  transform.setIdentity();
  try
  {
    // 使用tf_listener对象调用lookupTransform函数来查询机器人坐标系 (base_frame_id) 到激光雷达坐标系 (last_scan.header.frame_id) 的变换关系
    tf_listener.lookupTransform(base_frame_id, last_scan.header.frame_id, ros::Time(0), transform);
    retrieved = true;
  }
  catch (tf::TransformException &ex)
  {
    ROS_ERROR("%s",ex.what());
    ros::Duration(1.0).sleep();
    retrieved = false;
  }

  //TF:transform -> Eigen::Isometry3d

  const tf::Matrix3x3 &basis = transform.getBasis();
  Eigen::Matrix3d R;

  for(int r = 0; r < 3; r++)  // 从transform中获取旋转矩阵 basis，并将其转化为一个Eigen的Matrix3d类型的矩阵 R
    for(int c = 0; c < 3; c++)
      R(r,c) = basis[r][c];

  Pose3d laser_tf(R);

  const tf::Vector3 &t = transform.getOrigin(); // 设置变换矩阵T的平移量
  laser_tf.translation()(0) = t[0];
  laser_tf.translation()(1) = t[1];
  laser_tf.translation()(2) = t[2];

  setLaserPose(laser_tf);

  return retrieved;
}

bool CLaserOdometry2DNode::scan_available()
{
  return new_scan_available;
}

void CLaserOdometry2DNode::process(const ros::TimerEvent&)
{
  ROS_INFO_STREAM("is_initialized():" << is_initialized());
  ROS_INFO_STREAM("scan_available():" << scan_available());
  if( is_initialized() && scan_available() )
  {
    //Process odometry estimation
    odometryCalculation(last_scan);   // 里程计真正的处理流程
    publish();                        // 发布算法处理后的结果
    new_scan_available = false; //avoids the possibility to run twice on the same laser scan
  }
  else
  {
    ROS_WARN("Waiting for laser_scans....") ;
  }
}

//-----------------------------------------------------------------------------------
//                                   CALLBACKS
//-----------------------------------------------------------------------------------

void CLaserOdometry2DNode::LaserCallBack(const sensor_msgs::LaserScan::ConstPtr& new_scan)
{
  ROS_INFO("----------------------------------");
  if (GT_pose_initialized)
  {
    // Keep in memory the last received laser_scan
    last_scan = *new_scan;  // 传地址，解引用
    current_scan_time = last_scan.header.stamp;

    //Initialize module on first scan
    if (!first_laser_scan)
    {
      // copy laser scan to internal variable
      for (unsigned int i = 0; i<width; i++)
        range_wf(i) = new_scan->ranges[i];  // 订阅发布的原始激光数据
      new_scan_available = true;            // 接收到新的数据，更改标志位new_scan_available
    }
    else
    {
      // 如果是第一次处理数据，还要初始化
      init(last_scan, initial_robot_pose.pose.pose);
      first_laser_scan = false;
    }
  }
}

// 如果初始没有指定位置，那么这里不会执行
void CLaserOdometry2DNode::initPoseCallBack(const nav_msgs::Odometry::ConstPtr& new_initPose)
{
  //Initialize module on first GT pose. Else do Nothing!
  if (!GT_pose_initialized)
  {
    initial_robot_pose = *new_initPose;
    GT_pose_initialized = true;
  }
}

void CLaserOdometry2DNode::publish()
{
  // first, we'll publish the odometry over tf
  //---------------------------------------
  if (publish_tf)   // 是否发布tf坐标
  {
    // ROS_INFO("[rf2o] Publishing TF: [base_link] to [odom]");
    geometry_msgs::TransformStamped odom_trans;
    odom_trans.header.stamp = ros::Time::now();
    odom_trans.header.frame_id = odom_frame_id;
    odom_trans.child_frame_id = base_frame_id;
    odom_trans.transform.translation.x = robot_pose_.translation()(0);
    odom_trans.transform.translation.y = robot_pose_.translation()(1);
    odom_trans.transform.translation.z = 0.0;
    odom_trans.transform.rotation = tf::createQuaternionMsgFromYaw(rf2o::getYaw(robot_pose_.rotation()));
    //send the transform
    odom_broadcaster.sendTransform(odom_trans);   // 发送里程计信息
  }

  //next, we'll publish the odometry message over ROS
  //-------------------------------------------------
  //ROS_INFO("[rf2o] Publishing Odom Topic");   
  nav_msgs::Odometry odom;    // 发布里程计信息
  odom.header.stamp = ros::Time::now();
  odom.header.frame_id = odom_frame_id;
  //set the position
  odom.pose.pose.position.x = robot_pose_.translation()(0);
  odom.pose.pose.position.y = robot_pose_.translation()(1);
  odom.pose.pose.position.z = 0.0;
  odom.pose.pose.orientation = tf::createQuaternionMsgFromYaw(rf2o::getYaw(robot_pose_.rotation()));
  //set the velocity
  odom.child_frame_id = base_frame_id;
  odom.twist.twist.linear.x = lin_speed;    //linear speed
  odom.twist.twist.linear.y = 0.0;
  odom.twist.twist.angular.z = ang_speed;   //angular speed
  //publish the message
  odom_pub.publish(odom);
}

} /* namespace rf2o */

//-----------------------------------------------------------------------------------
//                                   MAIN
//-----------------------------------------------------------------------------------
int main(int argc, char** argv)
{
  ros::init(argc, argv, "RF2O_LaserOdom");

  rf2o::CLaserOdometry2DNode myLaserOdomNode;

  ros::TimerOptions timer_opt;
  timer_opt.oneshot   = false;
  timer_opt.autostart = true;
  timer_opt.callback_queue = ros::getGlobalCallbackQueue();
  timer_opt.tracked_object = ros::VoidConstPtr();

  timer_opt.callback = boost::bind(&rf2o::CLaserOdometry2DNode::process, &myLaserOdomNode, _1);
  timer_opt.period   = ros::Rate(myLaserOdomNode.freq).expectedCycleTime();

  ros::Timer rf2o_timer = ros::NodeHandle("~").createTimer(timer_opt);

  ros::spin();

  return EXIT_SUCCESS;
}
