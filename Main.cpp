#include"Functions.h"


void test() {
	Eigen::Vector3f centroid_A(0, 0, 0);
	cout << centroid_A.size() << endl;
	PC_XYZ::Ptr src(new PC_XYZ);
	PC_XYZ::Ptr dst(new PC_XYZ);
	(*src).push_back(XYZ(1, 2, 9));
	
	//std::cout << (*src)[0].x << " " << (*src)[0].y << " " << (*src)[0].z << std::endl;
	//std::cout << (*src)[0].getArray3fMap() << std::endl;
	//std::cout << centroid_A << std::endl;
	Eigen::Matrix3f aa;
	aa << 1, 2, 3, 4, 5, 6, 7, 8, 9;
	Eigen::Matrix3f ab;
	ab=aa;
	aa(0, 0)=10;
	std::cout << ab << std::endl;
	//std::cout << aa.transpose().transpose() << std::endl;
}

void main() {
	registration();
}