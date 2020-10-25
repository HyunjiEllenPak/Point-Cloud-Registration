#include"Functions.h"
#include"Kd_tree.h"
//Read stl file
XYZ read_Point(std::ifstream &s) {
	//save 4bytes one point in triangle
	char* facet=new char[12];//4byte*3:(x,y,z)
	s.read(facet, 12);
	char f1[4] ={ facet[0],
		facet[1],facet[2],facet[3] };

	char f2[4] ={ facet[4],
		facet[5],facet[6],facet[7] };

	char f3[4] ={ facet[8],
		facet[9],facet[10],facet[11] };

	float x = *((float*)f1);
	float y = *((float*)f2);
	float z = *((float*)f3);

	return XYZ(x, y, z);
}
void read_stl(std::string fname, PC_XYZ::Ptr points, char* header_info) {
	std::ifstream stlFile(fname.c_str(), std::ios::binary);
	if (stlFile.fail()) {
		std::cout << "Error:COULD NOT READ FILE" << std::endl;
		assert(false);
	}
	char numTri[4];
	//read data as a block, read 80 byte header+1
	stlFile.read(header_info, 80);
	//80-83th array is about the number of triangles in this mesh
	//convert 4byte char to unsigned long
	stlFile.read(numTri, 4);
	assert(sizeof(long) == 4);
	unsigned long* r = (unsigned long*)numTri;
	unsigned long nTri = *r;
	for (int i=0; i <nTri; i++) {
		auto normal=read_Point(stlFile);
		auto v1=read_Point(stlFile);
		auto v2=read_Point(stlFile);
		auto v3=read_Point(stlFile);
		//tri triangle(normal,v1, v2, v3);
		points->push_back(v1);
		points->push_back(v2);
		points->push_back(v3);
		char dummy[2];
		stlFile.read(dummy, 2);
	}
}
void make_Cloud(PC_XYZ::Ptr tris, PC_XYZ::Ptr points) {
	unsigned long num_points=(*tris).size();
	std::map<std::string, int> table;
	//clock_t begin, end;
	//begin=clock();
	for (int i=0; i < num_points; i++) {
		XYZ p=(*tris)[i];
		std::ostringstream temp_x;
		std::ostringstream temp_y;
		std::ostringstream temp_z;
		temp_x << p.x;
		temp_y << p.y;
		temp_z << p.z;
		std::string temp_string=temp_x.str();
		temp_string.append("#");
		temp_string.append(temp_y.str());
		temp_string.append("#");
		temp_string.append(temp_z.str());
		std::pair<std::string, int> item(temp_string, i);
		table.insert(item);
		//if (i % (num_points / 10) == 0) {
		//	end=clock();
		//	std::cout << "Elapsed time: " << end - begin << " (ms)" << std::endl;
		//}
	}
	std::map<std::string, int>::iterator iter;//x좌표+#+y좌표+#+z좌표 의 string을 만들어서 중복되는 것 제거.
	for (iter=table.begin(); iter != table.end(); iter++) {
		XYZ temp=(*tris)[iter->second];
		points->push_back(temp);
	}
}
//PCA registration
XYZ findCentroid(PC_XYZ::Ptr in) {
	unsigned long len=in->size();
	float sum_x=0;
	float sum_y=0;
	float sum_z=0;
	for (int i=0; i < len; i++) {
		sum_x+=(*in)[i].x;
		sum_y+=(*in)[i].y;
		sum_z+=(*in)[i].z;
	}
	sum_x/=float(len);
	sum_y/=float(len);
	sum_z/= float(len);
	
	return XYZ(sum_x, sum_y, sum_z);
}
MF covariance(PC_XYZ::Ptr in, XYZ cloud_mean) {
	MF cov(3, 3);
	cov << 0, 0, 0,
		0, 0, 0,
		0, 0, 0;
	int num=(*in).size();
	for (int i=0; i < num; i++) {
		XYZ p=(*in)[i];
		p.x-=cloud_mean.x;
		p.y-=cloud_mean.y;
		p.z-=cloud_mean.z;
		cov(0, 0)+=SQ(p.x);
		cov(0, 1)+=(p.x*p.y);
		cov(0, 2)+=(p.x*p.z);
		cov(1, 1)+=SQ(p.y);
		cov(1, 2)+=(p.y*p.z);
		cov(2, 2)+=SQ(p.z);
	}
	cov(1, 0)=cov(0, 1);
	cov(2, 0)=cov(0, 2);
	cov(2, 1)=cov(1, 2);
	cov=cov/(num-1);
	return cov;
}
MF calculateEigenVector(PC_XYZ::Ptr in, XYZ cloud_mean) {
	MF cov=covariance(in, cloud_mean);
	Eigen::EigenSolver<Eigen::MatrixXf> eigensolver(cov);
	Eigen::VectorXcf eigenVal=eigensolver.eigenvalues();
	MF eigenVec=eigensolver.eigenvectors().real();

	return eigenVec;
}

void writePCD(char* fname, PC_XYZ::Ptr in) {
	pcl::io::savePCDFileASCII(fname, *in);
	std::cerr << "Saved " <<(*in).size() << " data points to " << fname << std::endl;
}
void pcaRegistration(PC_XYZ::Ptr ref, PC_XYZ::Ptr floating,PC_XYZ::Ptr tr_floating) {
	std::cout << "PCA starts !" << std::endl;
	NEIGHBOR neigh;
	float error=0;
	float minE=FLT_MAX;
	int num=(*floating).size();
	XYZ ref_mean=findCentroid(ref);
	XYZ float_mean=findCentroid(floating);
	MF ref_eigenVT=calculateEigenVector(ref, ref_mean);
	MF float_eigenVT=calculateEigenVector(floating, float_mean);
	MF ref_eigen_copy=ref_eigenVT;
	MF rot3(3, 3);
	MF rot4(4, 4);
	rot4.setIdentity();
	MF min_rot4(4, 4);
	int sign1[4]={ 1,-1,1,-1 };
	int sign2[4]={ 1,1,-1,-1 };
	for (int nn=0; nn < 4; nn++) {
		ref_eigenVT.col(0) = sign1[nn] * ref_eigen_copy.col(0);
		ref_eigenVT.col(1) = sign2[nn] * ref_eigen_copy.col(1);
		rot3=ref_eigenVT *float_eigenVT.transpose();
		rot4.block<3, 3>(0, 0)=rot3.block<3, 3>(0, 0);
		tr_floating->clear();
		for (int i=0; i < num; i++) {
			VF point(4);
			point << (*floating)[i].x, (*floating)[i].y, (*floating)[i].z, 1;
			VF mul=rot4*point;		
			tr_floating->push_back(XYZ(mul(0), mul(1), mul(2)));
		}
		NEIGHBOR neigh;
		nearest_neighbor(tr_floating, ref, &neigh);

		error=std::accumulate(neigh.distances.begin(), neigh.distances.end(), 0.0) / neigh.distances.size();
		if (nn<minE) {
			minE=error;
			min_rot4=rot4;
		}
	}//nn
	tr_floating->clear();
	for (int i=0; i < num; i++) {
		VF p(4);
		p << (*floating)[i].x, (*floating)[i].y, (*floating)[i].z, 1;
		VF mul=min_rot4*p;
		tr_floating->push_back(XYZ(mul(0), mul(1), mul(2)));
	}

	//(*corr).clear();
	//for (int i=0; i < num; i++) {
	//	(*corr).push_back((*ref)[min_neigh.indices[i]]);
	//}

	//MF world_to_floatCenter(4, 4);
	//world_to_floatCenter << 1, 0, 0, -float_mean.x,
	//			 0, 1, 0, -float_mean.y,
	//			 0, 0, 1, -float_mean.z,
	//			 0, 0, 0, 1;
	//MF floatCenter_to_refCenter(4, 4);
	//floatCenter_to_refCenter << 1, 0, 0, ref_mean.x,
	//			 0, 1, 0, ref_mean.y,
	//			 0, 0, 1, ref_mean.z,
	//			 0, 0, 0, 1;
	//MF init=floatCenter_to_refCenter*rot*world_to_floatCenter;
	//return rot;
	//MF tf= findInitialTransform(ref, floating);
	std::cout << "PCA is finished !" << std::endl;
	
}
//ICP registration
void nearest_neighbor(PC_XYZ::Ptr src, PC_XYZ::Ptr dst, NEIGHBOR *neigh) {
	pcl::KdTreeFLANN<XYZ> kdtree;
	kdtree.setInputCloud(dst);
	std::vector<int> idx(1);
	std::vector<float> dist(1);

	int src_num=(*src).size();
	int dst_num=(*dst).size();
	int index=0;
	(*neigh).distances.clear();
	(*neigh).indices.clear();
	for (int i=0; i < src_num; i++) {
		float min_dist=FLT_MAX;
		XYZ ps=(*src)[i];
		if (kdtree.nearestKSearch(ps, 1, idx, dist) > 0) {
			index=idx[0];
			min_dist=dist[0];
		}
		//for (int j=0; j < dst_num; j++) {
		//	XYZ pd=(*dst)[j];
		//	float dist=distance(ps, pd);
		//	if (dist < min) {
		//		min=dist;
		//		index=j;
		//	}
		//}
		(*neigh).distances.push_back(min_dist);
		(*neigh).indices.push_back(index);
	}	
}
float distance(XYZ pa, XYZ pb) {
	return sqrt(SQ(pa.x - pb.x) + SQ(pa.y - pb.y) + SQ(pa.z - pb.z));
}
//Find the transformation to transform src to dst
Eigen::Matrix4f best_alignment(PC_XYZ::Ptr src, PC_XYZ::Ptr dst) {
	int snum=(*src).size();
	int dnum=(*dst).size();
	//Point number check
	if (snum != dnum) {
		std::cout << "The length is not matched !" << std::endl;
	}
	XYZ centroid_s=findCentroid(src);
	XYZ centroid_d=findCentroid(dst);
	std::cout << "1 =" <<centroid_s << std::endl;
	std::cout << "2 =" << centroid_d << std::endl;
	PC_XYZ trans_src;
	PC_XYZ trans_dst;

	//These matrix should be initialized
	Eigen::Matrix4f T;
	T.setIdentity();
	Eigen::Matrix3f H;
	H.setZero();

	//Translate src and dst axis to the centroid axis
	for (int i=0; i < snum; i++) {
		XYZ p1((*src)[i].x - centroid_s.x, (*src)[i].y - centroid_s.y, (*src)[i].z - centroid_s.z);
		XYZ p2((*dst)[i].x - centroid_d.x, (*dst)[i].y - centroid_d.y, (*dst)[i].z - centroid_d.z);
		//std::cout << p1 << " " << p2 << std::endl;		
		H(0, 0)+=p1.x*p2.x;
		H(0, 1)+=p1.x*p2.y;
		H(0, 2)+=p1.x*p2.z;

		H(1, 0)+=p1.y*p2.x;
		H(1, 1)+=p1.y*p2.y;
		H(1, 2)+=p1.y*p2.z;

		H(2, 0)+=p1.z*p2.x;
		H(2, 1)+=p1.z*p2.y;
		H(2, 2)+=p1.z*p2.z;
		//std::cout << H(2, 0) << std::endl;
	}


	Eigen::JacobiSVD<Eigen::MatrixXf> svd(H, Eigen::ComputeFullU|Eigen::ComputeFullV);
	Eigen::Matrix3f U=svd.matrixU();
	Eigen::Matrix3f Ut=U.transpose();
	Eigen::Matrix3f V=svd.matrixV();
	Eigen::Matrix3f R=V*Ut;
	Eigen::Vector3f t;
	if (R.determinant() < 0 ) {
		std::cout << "det = " << R.determinant()<< std::endl;
		Ut.block<1,3>(2,0)*=-1;
		R=V*Ut;
	}
	t=centroid_d.getVector3fMap() - R*centroid_s.getVector3fMap();
	T.block<3, 3>(0, 0)=R;
	T.block<3, 1>(0, 3)=t;
	std::cout << "rotation matrix" << std::endl;
	std::cout << R << std::endl;
	std::cout << "translation matrix" << std::endl;
	std::cout << t << std::endl;
	return T;

}
void rejectPairs(PC_XYZ::Ptr src, PC_XYZ::Ptr dst) {
	int num=(*src).size();
	float epsilon=0.1;
	for (int i=0; i < num; i++) {
		XYZ temp=(*src)[i];
		for (int j=0; j < num; j++) {
			//if()
		}
		
	}
}
//In:PCA transformed floating points, reference points, icp out points
ICP_OUT icp_registration(PC_XYZ::Ptr src, PC_XYZ::Ptr dst,PC_XYZ::Ptr tf_src, int max_iteration, float epsilon) {
	std::cout << "ICP starts !" << std::endl;
	float pre_error=FLT_MAX;
	float cur_error=0;
	int src_num=(*src).size();
	ICP_OUT result;
	int iter=0;
	NEIGHBOR neigh;
	Eigen::Matrix4f T;
	Eigen::Matrix4f totalT;
	totalT.setIdentity();
	(*tf_src).clear();
	//Deep copy
	for (int i=0; i < src_num; i++) {
		(*tf_src).push_back((*src)[i]);
	}
	PC_XYZ::Ptr corresponding(new PC_XYZ);
	for (int i=0; i < max_iteration; i++) {
		std::cout << "Iteration_number= " << i + 1 << std::endl;
		(*corresponding).clear();
		nearest_neighbor(tf_src, dst, &neigh);

		//Copy corresponding neighbor 
		for (int n1=0; n1 < src_num; n1++) {
			XYZ ss=(*src)[n1];
			XYZ nn=(*dst)[neigh.indices[n1]];
			(*corresponding).push_back((*dst)[neigh.indices[n1]]);
		}
		
		T=best_alignment(tf_src,corresponding);
		totalT=T*totalT;
		for (int n2=0; n2 < src_num; n2++) {
			Eigen::Vector4f temp((*tf_src)[n2].x, (*tf_src)[n2].y, (*tf_src)[n2].z,1);
			temp=T*temp.eval();		
			(*tf_src)[n2].x=temp(0);
			(*tf_src)[n2].y=temp(1);
			(*tf_src)[n2].z=temp(2);
		}
		cur_error=std::accumulate(neigh.distances.begin(), neigh.distances.end(), 0.0) / neigh.distances.size();
		std::cout << "cur_error="<<cur_error << std::endl;
		if (abs(pre_error - cur_error) < epsilon) {
			iter=i + 1;
			break;
		}
		pre_error=cur_error;
	}//iteration
	std::cout << "ICP is finished !" << std::endl;
	for (int i=0; i < src_num; i++) {
		Eigen::Vector4f a((*src)[i].x,(*src)[i].y,(*src)[i].z, 1);
		Eigen::Vector4f temp=totalT*a;
		(*tf_src)[i]=XYZ(temp(0), temp(1), temp(2));
	}
	result.iter=iter;
	result.T=totalT;
	return result;
}

void registration() {
	//Extract points from stl files
	//Reference points
	PC_XYZ::Ptr ref_tris(new PC_XYZ);
	std::string ref_fname="../seminar_data/Reference_image.stl";
	char* ref_header_info=new char[80];
	read_stl(ref_fname, ref_tris, ref_header_info);
	PC_XYZ::Ptr ref_points(new PC_XYZ);
	make_Cloud(ref_tris, ref_points);//907116->(if unified)151633
									 //cout << (*ref_points).size() << endl;
									 //Floating points
	PC_XYZ::Ptr float_tris(new PC_XYZ);
	std::string float_fname="../seminar_data/Floating_image.stl";
	char* float_header_info=new char[80];
	read_stl(float_fname, float_tris, float_header_info);
	PC_XYZ::Ptr float_points(new PC_XYZ);
	make_Cloud(float_tris, float_points);//166232->(if unified)83501
										 //cout << (*float_points).size() << endl;

										 //PCA analysis
	PC_XYZ::Ptr tr_float_points(new PC_XYZ);//83501
	PC_XYZ::Ptr corr(new PC_XYZ);
	std::cout << "Counter start !" << std::endl;
	LARGE_INTEGER Frequency;
	LARGE_INTEGER BeginTime;
	LARGE_INTEGER EndTime;
	QueryPerformanceFrequency(&Frequency);
	QueryPerformanceCounter(&BeginTime);
	pcaRegistration(ref_points, float_points, tr_float_points);
	PC_XYZ::Ptr icp_points(new PC_XYZ);

	//In:PCA transformed floating points, reference points, icp out points
	icp_registration(tr_float_points, ref_points,icp_points,500, 0.000001);
	QueryPerformanceCounter(&EndTime);
	__int64 elapsedTime=(EndTime.QuadPart - BeginTime.QuadPart) / Frequency.QuadPart;
	std::cout << "PCA and ICP elapsed time is " << elapsedTime << "(s)" << std::endl;
	//Write a pcd file
	//Fill in the data

	//PC_XYZ::Ptr icp_points2;
	//pcl::PCLPointCloud2 icp_ply;
	//pcl::toPCLPointCloud2(*icp_points, icp_ply);
	//pcl::io::savePLYFileASCII("icp.ply", *icp_points);
	//writer.write("icp_points.ply", icp_ply, Eigen::Vector4f::Zero(), Eigen::Quaternionf::Identity(), 0, 0);
	//(*icp_points2).is_dense=false;// loat_points.is_dense;	
	//writePCD("Transformed_floating_points.pcd", icp_points2);
	//writePCD("floating_points.pcd", float_points);

	pcl::visualization::PCLVisualizer viewer("ICP demo");
	//Create two vertically seperated viewports
	int v1(0);
	int v2(1);
	viewer.createViewPort(0.0, 0.0, 0.5, 1.0, v1);
	viewer.createViewPort(0.5, 0.0, 1.0, 1.0, v2);

	//The colors will be using
	float bckgr_gray_level=0.0;//Black
	float txt_gray_lvl=1.0 - bckgr_gray_level;

	//Original point cloud is white
	pcl::visualization::PointCloudColorHandlerCustom<XYZ> white(ref_points, (int)255 * txt_gray_lvl, (int)255 * txt_gray_lvl, (int)255 * txt_gray_lvl);
	viewer.addPointCloud(corr, white, "cloud_in_v1", v1);
	viewer.addPointCloud(ref_points, white, "cloud_in_v2", v2);

	//Transformed point cloud is green
	pcl::visualization::PointCloudColorHandlerCustom<XYZ> green(float_points, 20, 180, 20);	
	viewer.addPointCloud(tr_float_points, green, "color_tr_v1", v1);

	//ICP alinged point cloud is red
	pcl::visualization::PointCloudColorHandlerCustom<XYZ> red(icp_points, 180, 20, 20);
	viewer.addPointCloud(icp_points, red, "cloud_icp_v2", v2);
	

	// Adding text descriptions in each viewport
	viewer.addText("White: Original point cloud\nGreen: Matrix transformed point cloud", 10, 15, 16, txt_gray_lvl, txt_gray_lvl, txt_gray_lvl, "icp_info_1", v1);
	viewer.addText("White: Original point cloud\nRed: ICP aligned point cloud", 10, 15, 16, txt_gray_lvl, txt_gray_lvl, txt_gray_lvl, "icp_info_2", v2);

	//std::stringstream ss;
	//ss << iterations;
	//std::string iterations_cnt = "ICP iterations = " + ss.str();
	//viewer.addText(iterations_cnt, 10, 60, 16, txt_gray_lvl, txt_gray_lvl, txt_gray_lvl, "iterations_cnt", v2);

	// Set background color
	viewer.setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, v1);
	viewer.setBackgroundColor(bckgr_gray_level, bckgr_gray_level, bckgr_gray_level, v2);

	// Set camera position and orientation
	viewer.setCameraPosition(-3.68332, 2.94092, 5.71266, 0.289847, 0.921947, -0.256907, 0);
	viewer.setSize(1280, 1024);  // Visualiser window size

								 // Register keyboard callback :
								 //viewer.registerKeyboardCallback(&keyboardEventOccurred, (void*)NULL);

								 // Display the visualiser
	while (!viewer.wasStopped())
	{
		viewer.spinOnce();

	}
}