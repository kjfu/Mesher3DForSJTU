#include "mesher3d_core.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include "mesher3d_io.h"
#include <iostream>
#include <ctime>
#include "mesh.h"



void delaunayTetrahedralization(tetgenio *in, tetgenio *out, REAL size, std::vector<int> &tetMarkers){

	tetgenio tmp_inside_in;
	tetgenio tmp_inside_out;

	tetgenio tmp_outside_in;
	tetgenio tmp_outside_out;

	std::vector<int> insidePointIndices;
	std::vector<int> outsidePointIndices;
	int np = in->numberofpoints;

	for (int i=0; i< in->numberofpoints; i++){
		if (in->pointmarkerlist[i] == 0){
			insidePointIndices.push_back(i);
		}
		else{
			outsidePointIndices.push_back(i);
		}
	}

	/**
	 * @brief copy points  with specific index to another tetgenio
	 * 
	 */
	auto copyPoint = 
	[&in]
	(auto &tetgenContainer, auto &indexContainer, double *center){
		center[0] = 0;
		center[1] = 0;
		center[2] = 0;
		tetgenContainer.numberofpoints = static_cast<int>(indexContainer.size());
		tetgenContainer.pointlist = new REAL[tetgenContainer.numberofpoints*3];
		tetgenContainer.pointmarkerlist = new int[tetgenContainer.numberofpoints];
		for (int i=0; i<tetgenContainer.numberofpoints; i++){
			tetgenContainer.pointlist[i*3+0] = in->pointlist[indexContainer[i]*3+0];
			tetgenContainer.pointlist[i*3+1] = in->pointlist[indexContainer[i]*3+1];
			tetgenContainer.pointlist[i*3+2] = in->pointlist[indexContainer[i]*3+2];
			center[0] += tetgenContainer.pointlist[i*3+0];
			center[1] += tetgenContainer.pointlist[i*3+1];
			center[2] += tetgenContainer.pointlist[i*3+2];
			tetgenContainer.pointmarkerlist[i] = in->pointmarkerlist[indexContainer[i]];
		}
		center[0] /= static_cast<double>(tetgenContainer.numberofpoints);
		center[1] /= static_cast<double>(tetgenContainer.numberofpoints);
		center[2] /= static_cast<double>(tetgenContainer.numberofpoints);
	};

	double center[3];	
	copyPoint(tmp_outside_in, outsidePointIndices, center);
	copyPoint(tmp_inside_in, insidePointIndices, center);


	/**
	 * @brief construct convex hull of points
	 * 
	 */
	tetgenbehavior b0;
	char command0[] = "f";
	b0.parse_commandline(command0);
	fprintf(stdout, "*****Tetrahedralize inside points!\n");
	tetrahedralize(&b0, &tmp_inside_in, &tmp_inside_out);
	fprintf(stdout, "*****Tetrahedralize outside points!\n");
	tetrahedralize(&b0, &tmp_outside_in, &tmp_outside_out);

	/**
	 * @brief Get the max edge in a tethedra
	 * 
	 */
	double h = 0;
	auto maxEdgeSizeInTet = 
	[&tmp_inside_out]
	(auto firstpoint){
		double hmax =0;
		for(int i=0; i<4; i++){
			int p0 = firstpoint[i]*3;
			int p1 = firstpoint[(i+1)%4]*3;
			double *point = tmp_inside_out.pointlist;
			double length = 
			sqrt(pow((point[p0+0] - point[p1+0]),2)
			+ pow((point[p0+1] - point[p1+1]),2)
			+ pow((point[p0+2] - point[p1+2]),2));
			hmax = hmax<length ? length : hmax;
		}
		return hmax;
	}; 

	for (int i=0; i<tmp_inside_out.numberoftetrahedra; i++){
		double hmax = maxEdgeSizeInTet(&tmp_inside_out.tetrahedronlist[4*i]);
		h = h < hmax ? hmax : h; 
	}
	
	size = size > h ? size : h;

	/**
	 * @brief Count boundary triangle facets on a convex hull
	 * 
	 */

	auto travelBoundaryTriangleFacets =
	[]
	(auto &tetgenContainer, auto &pointIndexContainer, auto &facetContainer){				
		std::map<int, int> indexmap;
		for(int i=0; i<tetgenContainer.numberoftrifaces; i++){
			if (tetgenContainer.trifacemarkerlist[i]==1){
				facetContainer.emplace_back();
				for (int j=0; j<3; j++){
					if (std::find(pointIndexContainer.begin(), pointIndexContainer.end(), tetgenContainer.trifacelist[3*i+j]) == pointIndexContainer.end()){
					
						indexmap[tetgenContainer.trifacelist[3*i+j]] = pointIndexContainer.size();
						facetContainer.back().push_back(pointIndexContainer.size());						
						pointIndexContainer.push_back(tetgenContainer.trifacelist[3*i+j]);
						
					}
					else{
						facetContainer.back().push_back(indexmap[tetgenContainer.trifacelist[3*i+j]]);
					}
				}

			}
		}
	};


	tetgenio tmp_shell_in;
	tetgenio tmp_shell_out;
	std::vector<int> outsideShellPointIndices;
	std::vector<int> insideShellPointIndices;
	std::vector<std::vector<int>> outsideFacets;
	std::vector<std::vector<int>> insideFacets;



	travelBoundaryTriangleFacets(tmp_outside_out, outsideShellPointIndices, outsideFacets);

	/**
	 * @brief  Refine outside mesh
	 * 
	 */

	tetgenio tmp_outside_refine_in;
	tetgenio tmp_outside_refine_out;
	tmp_outside_refine_in.numberoffacets = outsideFacets.size();
	tmp_outside_refine_in.facetlist = new tetgenio::facet[tmp_outside_refine_in.numberoffacets];
	tmp_outside_refine_in.numberofholes = 0;
	tmp_outside_refine_in.facetmarkerlist = new int[tmp_outside_refine_in.numberoffacets];
	for(int i=0; i<tmp_outside_refine_in.numberoffacets; i++){
		tmp_outside_refine_in.facetmarkerlist[i]=1;
		tetgenio::facet &f = tmp_outside_refine_in.facetlist[i];		
		f.numberofpolygons = 1;
		f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
		f.numberofholes = 0;
		f.holelist = NULL;
		tetgenio::polygon &p = f.polygonlist[0];
		p.numberofvertices = 3;
		p.vertexlist = new int[p.numberofvertices];
		p.vertexlist[0] = outsideFacets[i][0];
		p.vertexlist[1] = outsideFacets[i][1];
		p.vertexlist[2] = outsideFacets[i][2];	

	}
	tmp_outside_refine_in.numberofpoints = outsideShellPointIndices.size();
	tmp_outside_refine_in.pointlist = new REAL[3*tmp_outside_refine_in.numberofpoints];		
	tmp_outside_refine_in.pointmarkerlist = new int[tmp_outside_refine_in.numberofpoints];
	auto index = 0;		
	for(auto i: outsideShellPointIndices){
		tmp_outside_refine_in.pointlist[(index)*3 + 0] = tmp_outside_out.pointlist[3*i + 0];
		tmp_outside_refine_in.pointlist[(index)*3 + 1] = tmp_outside_out.pointlist[3*i + 1];
		tmp_outside_refine_in.pointlist[(index)*3 + 2] = tmp_outside_out.pointlist[3*i + 2];
		tmp_outside_refine_in.pointmarkerlist[index] = tmp_outside_out.pointmarkerlist[i];
		index++;
	}

	tmp_outside_refine_in.numberofpointmtrs = 1;
	tmp_outside_refine_in.pointmtrlist = new REAL[tmp_outside_refine_in.numberofpoints];
	for (int i=0; i<tmp_outside_refine_in.numberofpoints; i++){
		tmp_outside_refine_in.pointmtrlist[i] = size;
	}
	tetgenbehavior b0_refine;
	char command_refine[] = "pqmDf";
	b0_refine.parse_commandline(command_refine);
	fprintf(stdout, "*****Refine outside!\n");
	tetrahedralize(&b0_refine, &tmp_outside_refine_in, &tmp_outside_refine_out);
	/**
	 * @brief  Build a shell with tetgenio
	 * 
	 */
	outsideShellPointIndices.clear();
	outsideFacets.clear();
	travelBoundaryTriangleFacets(tmp_inside_out, insideShellPointIndices, insideFacets);
	travelBoundaryTriangleFacets(tmp_outside_refine_out, outsideShellPointIndices, outsideFacets);
	
	tmp_shell_in.numberofpoints = insideShellPointIndices.size() + outsideShellPointIndices.size();
	tmp_shell_in.pointlist = new REAL[tmp_shell_in.numberofpoints*3];
	tmp_shell_in.numberofpointmtrs = 1;
	tmp_shell_in.pointmarkerlist = new int[tmp_shell_in.numberofpoints];

	tmp_shell_in.numberoffacets = insideFacets.size() + outsideFacets.size();
	tmp_shell_in.facetlist = new tetgenio::facet[tmp_shell_in.numberoffacets];
	tmp_shell_in.facetmarkerlist = new int[tmp_shell_in.numberoffacets];
	tmp_shell_in.numberofholes = 1;
	tmp_shell_in.holelist = new double[3];
	tmp_shell_in.holelist[0] = center[0];
	tmp_shell_in.holelist[1] = center[1];
	tmp_shell_in.holelist[2] = center[2];

	auto addPoints =
	[&tmp_shell_in]
	(auto &tetgenioContainer, auto &pointIndexContainer, auto base){
		auto index = base;		
		for(auto i: pointIndexContainer){
			tmp_shell_in.pointlist[(index)*3 + 0] = tetgenioContainer.pointlist[3*i + 0];
			tmp_shell_in.pointlist[(index)*3 + 1] = tetgenioContainer.pointlist[3*i + 1];
			tmp_shell_in.pointlist[(index)*3 + 2] = tetgenioContainer.pointlist[3*i + 2];
			tmp_shell_in.pointmarkerlist[index] = tetgenioContainer.pointmarkerlist[i];
			index++;
		}
		// tmp_shell_in.pointmarkerlist[index] = 1;
	};
	addPoints(tmp_inside_out, insideShellPointIndices, 0);
	addPoints(tmp_outside_refine_out, outsideShellPointIndices, insideShellPointIndices.size());



	auto addTriangleFacet = 
	[&tmp_shell_in]
	(auto &tetgenioContainer, auto &facetContainer, int facetBase, int pointBase){
		for(int i=0; i<facetContainer.size(); i++){
			tetgenio::facet &f = tmp_shell_in.facetlist[i+facetBase];		
			f.numberofpolygons = 1;
			f.polygonlist = new tetgenio::polygon[f.numberofpolygons];
			f.numberofholes = 0;
			f.holelist = NULL;
			tetgenio::polygon &p = f.polygonlist[0];
			p.numberofvertices = 3;
			p.vertexlist = new int[p.numberofvertices];
			p.vertexlist[0] = facetContainer[i][0] + pointBase;
			p.vertexlist[1] = facetContainer[i][1] + pointBase;
			p.vertexlist[2] = facetContainer[i][2] + pointBase;
		}
	};
	addTriangleFacet(tmp_inside_out, insideFacets, 0, 0);
	addTriangleFacet(tmp_outside_refine_out, outsideFacets, insideFacets.size(), insideShellPointIndices.size());


	tmp_shell_in.numberofpointmtrs = 1;
	tmp_shell_in.pointmtrlist = new REAL[tmp_shell_in.numberofpoints];
    for(int i=0; i<insideShellPointIndices.size(); i++){

        tmp_shell_in.pointmtrlist[i] = h;
    }
	for(int i=0; i<outsideShellPointIndices.size(); i++){
		tmp_shell_in.pointmtrlist[insideShellPointIndices.size()+i] = size;
	}
	
	tetgenbehavior b1;
	char command1[] = "pqmY";
	b1.parse_commandline(command1);
	fprintf(stdout, "*****Tetrahedralize shell!\n");
	tetrahedralize(&b1, &tmp_shell_in, &tmp_shell_out);
	out->firstnumber = 1;
	out->numberoftetrahedra = tmp_shell_out.numberoftetrahedra + tmp_inside_out.numberoftetrahedra;
	out->tetrahedronlist =  new int[out->numberoftetrahedra*4];
	out->numberofpoints = tmp_shell_out.numberofpoints + tmp_inside_out.numberofpoints - insideShellPointIndices.size();

	out->pointlist = new REAL[out->numberofpoints * 3];
	out->pointmarkerlist = new int[out->numberofpoints];	
	for(int i=0; i<tmp_inside_out.numberofpoints; i++){
		out->pointlist[3*i + 0] = tmp_inside_out.pointlist[3*i + 0];
		out->pointlist[3*i + 1] = tmp_inside_out.pointlist[3*i + 1];
		out->pointlist[3*i + 2] = tmp_inside_out.pointlist[3*i + 2];
		out->pointmarkerlist[i] = 0;
	}	

	int minusPointBase = 0 + tmp_inside_out.numberofpoints - insideShellPointIndices.size(); 
	for(int i=insideShellPointIndices.size(); i<tmp_shell_out.numberofpoints; i++){
	
		out->pointlist[3*(i+minusPointBase) + 0] = tmp_shell_out.pointlist[3*i + 0];
		out->pointlist[3*(i+minusPointBase) + 1] = tmp_shell_out.pointlist[3*i + 1];
		out->pointlist[3*(i+minusPointBase) + 2] = tmp_shell_out.pointlist[3*i + 2];
		out->pointmarkerlist[i+minusPointBase] = 1;
	}

	for(int i=0; i<tmp_inside_out.numberoftetrahedra; i++){
		out->tetrahedronlist[4*i + 0] = tmp_inside_out.tetrahedronlist[4*i + 0]+1;
		out->tetrahedronlist[4*i + 1] = tmp_inside_out.tetrahedronlist[4*i + 1]+1;
		out->tetrahedronlist[4*i + 2] = tmp_inside_out.tetrahedronlist[4*i + 2]+1;
		out->tetrahedronlist[4*i + 3] = tmp_inside_out.tetrahedronlist[4*i + 3]+1;
		tetMarkers.push_back(0);
	}
	
	tetgenio final_convex_out;
	tetgenbehavior b2;
	char command2[] = "f";
	b2.parse_commandline(command2);
	//printf("asdadasdasdads %d\n", tmp_shell_out.numberofpoints);
	fprintf(stdout, "*****Tetrahedralize all points!\n");
	tetrahedralize(&b2, &tmp_shell_out, &final_convex_out);
	for(int i=0; i<final_convex_out.numberoftrifaces; i++){
		if (final_convex_out.trifacemarkerlist[i]==1){
			for(int j=0; j<3; j++){
				out->pointmarkerlist[final_convex_out.trifacelist[3*i+j]+minusPointBase] = 2;
			}
		}
	}



	auto updateIndex=
	[minusPointBase, &insideShellPointIndices]
	(auto pointIndex){
		if (pointIndex<insideShellPointIndices.size()){
			return insideShellPointIndices[pointIndex] + 1;
		}
		else{
			return pointIndex+minusPointBase+1;
		}
	};



	for (int i=0; i<tmp_shell_out.numberoftetrahedra; i++){
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 0]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+0]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 1]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+1]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 2]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+2]);
		out->tetrahedronlist[4*(tmp_inside_out.numberoftetrahedra+i) + 3]	= updateIndex(tmp_shell_out.tetrahedronlist[4*i+3]);
		tetMarkers.push_back(1);
	}	
	
	// int sumsum=0;
	// int queryIndex= 1;
	// for (int i=0; i<out->numberoftetrahedra; i++){
	// 	if (out->tetrahedronlist[4*i+0]==queryIndex || out->tetrahedronlist[4*i+1]==queryIndex  || out->tetrahedronlist[4*i+2]==queryIndex || out->tetrahedronlist[4*i+3]==queryIndex){
	// 		fprintf(stdout, "***%d  %d  %d  %d\n", out->tetrahedronlist[4*i+0], out->tetrahedronlist[4*i+1], out->tetrahedronlist[4*i+2], out->tetrahedronlist[4*i+3]);
	// 		sumsum++;
	// 	}
	// }
	// fprintf(stdout, "********************* %d\n", sumsum);
}






void constrainedTetrahedralization(tetgenio *in, tetgenio *out, REAL size){
	in->numberofpointmtrs = 1;
    in->pointmtrlist = new REAL[in->numberofpoints];
    for(int i=0; i<in->numberofpoints; i++){

        in->pointmtrlist[i] = size;
    }

	tetgenbehavior b;
	
	char command[] = "pqmD";
	b.parse_commandline(command);
    tetrahedralize(&b, in, out);

}

void refineMesh(const std::string &fileInHead, const std::string &fileOutHead){
		Mesh goalMesh;
		Mesh backgroundMesh;
		double timebegin, timeend;
		std::vector<int> refine_elements;
		std::vector<std::array<double,3>> append_points;
		goalMesh.loadMESH(fileInHead+".mesh");		
		loadREMESH(refine_elements, append_points, fileInHead+".remesh");		
		backgroundMesh.clone(goalMesh);
		backgroundMesh.loadNodeValues(fileInHead+".value");	
		std::cout << "[*************] Input nodes: "<<goalMesh.nodes.size() << "; input tets: "<< goalMesh.tetrahedrons.size() << std::endl;
		std::vector<Vector3D> positions;
		for(auto nt:refine_elements){
			Vector3D vec = goalMesh.tetrahedrons[nt-1]->center();
			positions.push_back(vec);
		}
		timebegin = clock_t();
		goalMesh.CavityBasedInsert(positions);
		timeend = clock_t();
		std::cout << "[*************] Refine time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
		std::cout << "[*************] After refine, nodes: "<<goalMesh.nodes.size() << "; tets: "<< goalMesh.tetrahedrons.size() << std::endl;		






		goalMesh.estimateSizing();
		std::vector<Tetrahedron *> tetsToRemove;
		std::vector<Node *> sNodes;
		for (auto &e: goalMesh.tetrahedrons){
			bool removeOut = true;
			for(auto n: e->nodes){
				if (n->label!=0){
					removeOut = false;
					break;
				}

			}
			if (removeOut){
				tetsToRemove.push_back(e);
			}
		}

		Mesh tmpMesh;
		for (auto n: goalMesh.nodes){
			if (n->label==0){
				Node *node = new Node(*n);
				tmpMesh.nodes.push_back(node);
			}
		}

		for(auto p:append_points){
			Node *node = new Node(p.data());
			tmpMesh.nodes.push_back(node);
		}




		tetgenio tmpout;
		tetgenio tmpin;
		Mesh gradMesh;
		std::vector<Node *> grad_nodes;
		std::vector<TriangleFacet> grad_facets; 

		transportNodesToTETGENIO(tmpMesh.nodes, tmpin);
		char cmd[]="f";
		tetrahedralize(cmd, &tmpin, &tmpout);
		instructTetrahedronConnectByTETGENIO(tmpMesh.nodes, tmpout, tmpMesh.tetrahedrons);
		tmpMesh.rebuildIndices();
		tmpMesh.estimateSizing();
		tmpMesh.extractBorder(grad_nodes, grad_facets);
		Vector3D hole;
		for(auto n: grad_nodes){
			hole+=n->pos;
		}
		hole/=grad_nodes.size();



		tmpMesh.readyForSpatialSearch();
		timebegin= clock();
		goalMesh.checkBooleanRemove(tmpMesh);
		timeend = clock();

		std::cout << "\n";
		std::cout <<"[*************]  Boolean remove time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;
		std::cout << "\n";

		std::vector<Tetrahedron *> rmTets;
		for(auto e: goalMesh.tetrahedrons){
			if(e->edit==-1){
				rmTets.push_back(e);
			}
		}
		extratctBorder(rmTets, grad_nodes, grad_facets);


		goalMesh.deleteLargeScaleTetrahedronsPermanently(rmTets);
		
		std::vector<Vector3D> holes;
		holes.emplace_back(0,0,0);
		tetgenio grad_in, grad_out;
		transportFacetsToTETGENIO(grad_nodes, grad_facets, holes, grad_in);
		char cmd0[]="pqmY";
		tetrahedralize(cmd0, &grad_in, &grad_out);
		gradMesh.loadTETGENIO(grad_out);

		std::vector<Node *>mergeNodes0;
		std::cout << "[*************] Start merge!\n";
		goalMesh.mergeMesh(gradMesh, mergeNodes0);
		for(auto n:goalMesh.nodes){
			n->label = 1;
		}
		for(auto e:goalMesh.tetrahedrons){
			e->label = 1;
		}


		for(auto n: tmpMesh.nodes){
			n->label = 0;
		}
		for(auto e:tmpMesh.tetrahedrons){
			e->label = 0;
		}
		std::vector<Node *>mergeNodes1;
		goalMesh.mergeMesh(tmpMesh, mergeNodes1);
		std::vector<Node *> nodes;
		extractBorderNodes(goalMesh.tetrahedrons, nodes);
		for(auto n: nodes){
			n->label = 2;
		}
		std::cout << "[*************] Finish merge!\n";
		timebegin= clock();
		backgroundMesh.interpolateNodeValuesForAnotherMesh(goalMesh);
		timeend = clock();
		std::cout <<"[*************] Interpolation time: "<< (timeend - timebegin)/ CLOCKS_PER_SEC << "s" <<std::endl;

		goalMesh.exportNodeValues(fileOutHead+".value");
		goalMesh.exportMESH(fileOutHead+".mesh");
		std::cout << "Finish Adaption!\n";
}