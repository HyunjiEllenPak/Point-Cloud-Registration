#ifndef _Functions_H_
#define _Functions_H_
#define _CRT_SECURE_NO_WARNINGS
#define _SCL_SECURE_NO_WARNINGS

#include<iostream>
#include<fstream>
#include<string>
#include<map>
#include<numeric>
#include<Eigen/Core>
#include<Eigen/Dense>
#include <pcl/point_cloud.h>
#include <pcl/io/pcd_io.h>
#include<pcl/point_types.h>
#include<pcl/visualization/cloud_viewer.h>
#include<pcl/visualization/pcl_visualizer.h>
#include <pcl/console/time.h> 
#include <pcl/io/ply_io.h>//ply write
#include<pcl/kdtree/kdtree_flann.h>
#include<pcl/kdtree/kdtree.h>
#define SQ(a) ((a)*(a))
typedef pcl::PointXYZ XYZ;
typedef pcl::PointCloud<XYZ> PC_XYZ;
typedef Eigen::MatrixXf MF;
typedef Eigen::VectorXf VF;

typedef struct neighbor{
	std::vector<float> distances;
	std::vector<int> indices;
}NEIGHBOR;
typedef struct icp_out {
	Eigen::Matrix4f T;
	int iter;
	
}ICP_OUT;
//Read stl file
XYZ read_Point(std::ifstream &s);
void read_stl(std::string fname, PC_XYZ::Ptr points, char* header_info);
void make_Cloud(PC_XYZ::Ptr tris,PC_XYZ::Ptr points);
//PCA registration
XYZ findCentroid(PC_XYZ::Ptr in);
MF covariance(PC_XYZ::Ptr in, XYZ cloud_mean);
MF calculateEigenVector(PC_XYZ::Ptr in, XYZ cloud_mean);
void writePCD(char* fname, PC_XYZ::Ptr in);
void pcaRegistration(PC_XYZ::Ptr ref, PC_XYZ::Ptr floating, PC_XYZ::Ptr tr_floating);
//ICP_registration
void nearest_neighbor(PC_XYZ::Ptr src, PC_XYZ::Ptr dst, NEIGHBOR *neigh);
float distance(XYZ pa, XYZ pb);
Eigen::Matrix4f best_alignment(PC_XYZ::Ptr src, PC_XYZ::Ptr dst);	
void rejectPairs(PC_XYZ::Ptr src, PC_XYZ::Ptr dst);
ICP_OUT icp_registration(PC_XYZ::Ptr src, PC_XYZ::Ptr dst, PC_XYZ::Ptr tf_src, int max_iteration, float epsilon);
void registration();
#endif
