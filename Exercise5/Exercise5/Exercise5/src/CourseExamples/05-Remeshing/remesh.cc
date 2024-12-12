//== INCLUDES =================================================================

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/Types/TriMesh_ArrayKernelT.hh>
#include <iostream>
#include <set>
//#include <unistd.h> //might be needed for unix (gcc)
#include <float.h>


//== IMPLEMENTATION ===========================================================


typedef OpenMesh::TriMesh_ArrayKernelT<>        Mesh;
Mesh                                            mesh;
OpenMesh::VPropHandleT<Mesh::Point>             update;


// https://www.graphics.rwth-aachen.de/media/openmesh_static/Documentations/OpenMesh-7.0-Documentation/a02238.html#af0bd67ac6c4c3fa8707191dcf28b7649

//-----------------------------------------------------------------------------


void  split_long_edges(float _long);
void  collapse_short_edges(float _short);
void  equalize_valences();
void  tangential_relaxation();


//-----------------------------------------------------------------------------


int main(int argc, char **argv)
{
	if (argc < 4) 
	{
		std::cerr << "Usage: \n" 
			<< argv[0] << " <edge-length>  <input_mesh>  <output_mesh>\n\n";
		exit(1);
	}
	float target_length = (float)atof(argv[1]);
	float low  = 4.0f/5.0f * target_length;
	float high = 4.0f/3.0f * target_length;


	// add required properties
	mesh.request_vertex_status();
	mesh.request_vertex_normals();
	mesh.request_edge_status();
	mesh.request_face_status();
	mesh.request_face_normals();
	mesh.add_property(update);


	// read mesh
	OpenMesh::IO::read_mesh(mesh, argv[2]);
	std::cout << "#vertices: " << mesh.n_vertices() << std::endl;
	mesh.update_normals();


	// main remeshing loop
	for (int i=0; i<5; ++i)
	{
		split_long_edges(high);
		collapse_short_edges(low);
		equalize_valences();
		tangential_relaxation();
	}


	// write mesh
	OpenMesh::IO::write_mesh(mesh, argv[3]);
}


//-----------------------------------------------------------------------------


void  split_long_edges(float _long)
{
	Mesh::EIter     e_it, e_end;
	Mesh::VHandle   v0, v1, vh;
	Mesh::EHandle   eh, e0, e1;
	Mesh::FHandle   f0, f1, f2, f3;
	bool            finished; 
	int             i;



	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
		{

			// Exercise 6.1 ----------------------------------------------
			// INSERT CODE:
			// If the edge is longer than _long
			//  1) add the midpoint to the mesh
			//  2) set the interpolated normal to the vertex
			//  3) split the edge with this vertex (use openMesh function split)
			// Leave the loop running until no splits are done (use the finished variable)
			// -----------------------------------------------------------
			eh = e_it.handle();

			Mesh::HalfedgeHandle heh0 = mesh.halfedge_handle(eh, 0);
			Mesh::HalfedgeHandle heh1 = mesh.halfedge_handle(eh, 1);

			if (!heh0.is_valid()) {
				v0 = mesh.from_vertex_handle(heh1);
				v1 = mesh.to_vertex_handle(heh1);
			}
			else {
				v0 = mesh.from_vertex_handle(heh0);
				v1 = mesh.to_vertex_handle(heh0);
			}

			Mesh::Point p0 = mesh.point(v0);
			Mesh::Point p1 = mesh.point(v1);

			Mesh::Point pmid = (p1 + p0) / 2;
			double len = (p1 - p0).norm();

			if (len > _long) {

				// Add midpoint to mesh and split the edge!
				Mesh::VertexHandle vhmid = mesh.add_vertex(pmid);
				mesh.split(eh, vhmid);

				// Update only affected face & vertex normals, respectively!
				if (!mesh.is_boundary(eh)) {
					Mesh::VertexFaceIter vf_it = mesh.vf_iter(vhmid);
					f0 = mesh.face_handle(vf_it); ++vf_it;
					f1 = mesh.face_handle(vf_it); ++vf_it;
					f2 = mesh.face_handle(vf_it); ++vf_it;
					f3 = mesh.face_handle(vf_it); 
					mesh.set_normal(f0, mesh.calc_face_normal(f0));
					mesh.set_normal(f1, mesh.calc_face_normal(f1));
					mesh.set_normal(f2, mesh.calc_face_normal(f2));
					mesh.set_normal(f3, mesh.calc_face_normal(f3));
					mesh.set_normal(vhmid, mesh.calc_vertex_normal(vhmid));
					for (Mesh::VertexVertexIter vv_it = mesh.vv_iter(vhmid); vv_it; ++vv_it) {
						vh = mesh.vertex_handle(vv_it);
						mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
					}
				}
				else {
					Mesh::VertexFaceIter vf_it = mesh.vf_iter(vhmid);
					f0 = mesh.face_handle(vf_it); ++vf_it;
					f1 = mesh.face_handle(vf_it); 
					mesh.set_normal(f0, mesh.calc_face_normal(f0));
					mesh.set_normal(f1, mesh.calc_face_normal(f1));
					mesh.set_normal(vhmid, mesh.calc_vertex_normal(vhmid));
					for (Mesh::VertexVertexIter vv_it = mesh.vv_iter(vhmid); vv_it; ++vv_it) {
						vh = mesh.vertex_handle(vv_it);
						mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
					}
				}

				// There can still be potentially longer edges => need to loop again to make sure!
				finished = false;
			}
		}
	}
	std::cout << "Done splitting!" << std::endl;
}


//-----------------------------------------------------------------------------


void  collapse_short_edges(float _short)
{
	Mesh::EIter     e_it, e_end;
	Mesh::CVVIter   vv_it;
	Mesh::VHandle   v0, v1;
	Mesh::HHandle   h0, h1, h01, h10;
	bool            finished, b0, b1;
	int             i;
	bool            hcol01, hcol10;

	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
		{
			if (!mesh.status(e_it).deleted()) // might already be deleted
			{

				// Exercise 6.2 ----------------------------------------------
				// INSERT CODE:
				// If the edge is shorter than _short
				//  1) Check if halfedge connects a boundary vertex with a non-boundary vertex. If so, don't collapse. Otherwise
				//  2) Check if halfedges collapsible
				//  3) Select the halfedge to be collapsed if at least one halfedge can be collapsed
				//  4) Collapse the halfedge
				// Leave the loop running until no collapse has been done (use the finished variable)
				// -----------------------------------------------------------
				Mesh::EdgeHandle eh = e_it.handle();

				if (mesh.is_boundary(eh)) continue;

				h0 = mesh.halfedge_handle(eh, 0);
				h1 = mesh.halfedge_handle(eh, 1);

				v0 = mesh.from_vertex_handle(h0);
				v1 = mesh.from_vertex_handle(h1);

				if ((mesh.is_boundary(v0) && !mesh.is_boundary(v1)) || (!mesh.is_boundary(v0) && mesh.is_boundary(v1))) {
					continue;
				}

				Mesh::Point p0 = mesh.point(v0);
				Mesh::Point p1 = mesh.point(v1);

				double len = (p1 - p0).norm();

				if (len < _short) {

					if (mesh.is_collapse_ok(h0) && mesh.is_collapse_ok(h1)) {
						
						if (mesh.valence(v1) >= mesh.valence(v0)) {
							mesh.collapse(h0);

							// Update the affected face normals around a vertex as well as all the vertices of affected regions.
							for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v1); vf_it; ++vf_it) {
								Mesh::FaceHandle fh = vf_it.handle();
								mesh.set_normal(fh, mesh.calc_face_normal(fh));
							}

							for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v1); vf_it; ++vf_it) {
								Mesh::FaceHandle fh = vf_it.handle();

								for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it; ++fv_it) {
									Mesh::VertexHandle vh = fv_it.handle();
									mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
								}
							}

						}
						else {
							mesh.collapse(h1);

							// Update the affected face normals around a vertex as well as all the vertices of affected regions.
							for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it; ++vf_it) {
								Mesh::FaceHandle fh = vf_it.handle();
								mesh.set_normal(fh, mesh.calc_face_normal(fh));
							}

							for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it; ++vf_it) {
								Mesh::FaceHandle fh = vf_it.handle();

								for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it; ++fv_it) {
									Mesh::VertexHandle vh = fv_it.handle();
									mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
								}
							}
						}
						
					}
					else if (mesh.is_collapse_ok(h0) && !mesh.is_collapse_ok(h1)) {
						mesh.collapse(h0);

						// Update the affected face normals around a vertex as well as all the vertices of affected regions.
						for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v1); vf_it; ++vf_it) {
							Mesh::FaceHandle fh = vf_it.handle();
							mesh.set_normal(fh, mesh.calc_face_normal(fh));
						}

						for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v1); vf_it; ++vf_it) {
							Mesh::FaceHandle fh = vf_it.handle();

							for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it; ++fv_it) {
								Mesh::VertexHandle vh = fv_it.handle();
								mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
							}
						}
					}
					else if (!mesh.is_collapse_ok(h0) && mesh.is_collapse_ok(h1)) {
						mesh.collapse(h1);

						// Update the affected face normals around a vertex as well as all the vertices of affected regions.
						for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it; ++vf_it) {
							Mesh::FaceHandle fh = vf_it.handle();
							mesh.set_normal(fh, mesh.calc_face_normal(fh));
						}

						for (Mesh::VertexFaceIter vf_it = mesh.vf_iter(v0); vf_it; ++vf_it) {
							Mesh::FaceHandle fh = vf_it.handle();

							for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh); fv_it; ++fv_it) {
								Mesh::VertexHandle vh = fv_it.handle();
								mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
							}
						}
					}
					else {
						continue;
					}

					finished = false;
				}

			}
		}
	}

	std::cout << "Done Collapsing!" << std::endl;

	mesh.garbage_collection();

	if (i==100) std::cerr << "collapse break\n";
}


//-----------------------------------------------------------------------------


void  equalize_valences()
{
	Mesh::EIter     e_it, e_end;
	Mesh::VHandle   v0, v1, v2, v3;
	Mesh::HHandle   hh;
	int             val0, val1, val2, val3;
	int             val_opt0, val_opt1, val_opt2, val_opt3;
	int             ve0, ve1, ve2, ve3, ve_before, ve_after;
	bool            finished;
	int             i;




	// flip all edges
	for (finished=false, i=0; !finished && i<100; ++i)
	{
		finished = true;

		for (e_it=mesh.edges_begin(), e_end=mesh.edges_end(); e_it!=e_end; ++e_it)
		{
			if (!mesh.is_boundary(e_it.handle()))
			{

				// Reset 
				ve_before = 0;
				ve_after = 0;

				// Exercise 6.3 ----------------------------------------------
				// INSERT CODE:
				//  1) Extract valences of the four vertices involved to an eventual flip.
				//  2) Compute the sum of the squared valence deviances before flip
				//  3) Compute the sum of the squared valence deviances after and eventual flip
				//  4) If valence deviance is decreased and flip is possible, flip the vertex
				// Leave the loop running until no collapse has been done (use the finished variable)
				// -----------------------------------------------------------
				Mesh::EdgeHandle eh = e_it.handle();

				hh = mesh.halfedge_handle(eh, 0);
				Mesh::HHandle hho = mesh.halfedge_handle(eh, 1);

				v0 = mesh.from_vertex_handle(hh);
				v1 = mesh.to_vertex_handle(hh);

				hh = mesh.next_halfedge_handle(hh);
				v2 = mesh.to_vertex_handle(hh);

				hho = mesh.next_halfedge_handle(hho);
				v3 = mesh.to_vertex_handle(hho);

				val0 = mesh.valence(v0);
				val1 = mesh.valence(v1);
				val2 = mesh.valence(v2);
				val3 = mesh.valence(v3);

				if (mesh.is_boundary(v0)) {
					ve_before += std::pow((val0 - 4), 2);
				}
				else {
					ve_before += std::pow((val0 - 6), 2);
				}

				if (mesh.is_boundary(v1)) {
					ve_before += std::pow((val1 - 4), 2);
				}
				else {
					ve_before += std::pow((val1 - 6), 2);
				}

				if (mesh.is_boundary(v2)) {
					ve_before += std::pow((val2 - 4), 2);
				}
				else {
					ve_before += std::pow((val2 - 6), 2);
				}

				if (mesh.is_boundary(v3)) {
					ve_before += std::pow((val3 - 4), 2);
				}
				else {
					ve_before += std::pow((val3 - 6), 2);
				}

				// ve_before = std::pow((val0 - 6), 2) + std::pow((val1 - 6), 2) + std::pow((val2 - 6), 2) + std::pow((val3 - 6), 2);
				
				val_opt0 = val0 - 1;
				val_opt1 = val1 - 1;
				val_opt2 = val2 + 1;
				val_opt3 = val3 + 1;

				if (mesh.is_boundary(v0)) {
					ve_after += std::pow((val_opt0 - 4), 2);
				}
				else {
					ve_after += std::pow((val_opt0 - 6), 2);
				}
				if (mesh.is_boundary(v1)) {
					ve_after += std::pow((val_opt1 - 4), 2);
				}
				else {
					ve_after += std::pow((val_opt1 - 6), 2);
				}
				if (mesh.is_boundary(v2)) {
					ve_after += std::pow((val_opt2 - 4), 2);
				}
				else {
					ve_after += std::pow((val_opt2 - 6), 2);
				}
				if (mesh.is_boundary(v3)) {
					ve_after += std::pow((val_opt3 - 4), 2);
				}
				else {
					ve_after += std::pow((val_opt3 - 6), 2);
				}

				//ve_after = std::pow(val_opt0 - 6, 2) + std::pow(val_opt1 - 6, 2) + std::pow(val_opt2 - 6, 2) + std::pow(val_opt3 - 6, 2);

				// I guess that is_flip_ok wrong in OpenMesh?: https://computergraphics.stackexchange.com/questions/12297/half-edge-criterion-to-check-if-an-edge-flip-is-illegal
				if (mesh.is_flip_ok(eh) && ve_after < ve_before) {

					mesh.flip(eh);
					finished = false;

					// Update face & vertex normals. Note: eh is not necessarily the edge handle for flipped edge
					Mesh::HHandle temp = mesh.find_halfedge(v2, v3);
					eh = mesh.edge_handle(temp);

					hh = mesh.halfedge_handle(eh, 0);
					hho = mesh.halfedge_handle(eh, 1);

					Mesh::FaceHandle fh0 = mesh.face_handle(hh);
					Mesh::FaceHandle fh1 = mesh.face_handle(hho);

					mesh.set_normal(fh0, mesh.calc_face_normal(fh0));
					mesh.set_normal(fh1, mesh.calc_face_normal(fh1));

					for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh0); fv_it; ++fv_it) {
						Mesh::VertexHandle vh = mesh.vertex_handle(fv_it);
						mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
					}

					for (Mesh::FaceVertexIter fv_it = mesh.fv_iter(fh1); fv_it; ++fv_it) {
						Mesh::VertexHandle vh = mesh.vertex_handle(fv_it);
						mesh.set_normal(vh, mesh.calc_vertex_normal(vh));
					}
				} 
			}	
		}
	}
	std::cout << "Done Flipping!" << std::endl;

	if (i==100) std::cerr << "flip break\n";
}


//-----------------------------------------------------------------------------


void  tangential_relaxation()
{
	Mesh::VIter     v_it, v_end(mesh.vertices_end());
	Mesh::CVVIter   vv_it;
	Mesh::Scalar    valence;
	Mesh::Point     u, n;


	// smooth
	for (int iters=0; iters<10; ++iters)
	{
		for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
		{
			if (!mesh.is_boundary(v_it))
			{
				// Exercise 6.4 ----------------------------------------------
				// INSERT CODE:
				//  1) Compute uniform laplacian approximation vector
				//  2) Compute the tangential component of the laplacian vector
				//  3) Store smoothed vertex location in the update vertex property.
				//     (you don't have to use 1/2 attenuation in this case, it's fine without attenuation)
				// -----------------------------------------------------------
				Mesh::VertexHandle vh = v_it.handle();

				Mesh::Point Lu(0.0, 0.0, 0.0);

				int n = 0;
				for (Mesh::VertexVertexIter vv_it = mesh.vv_iter(vh); vv_it; ++vv_it) {

					Mesh::VertexHandle vvh = vv_it.handle();
					Lu += mesh.point(vvh); n++;
				}

				// Perform Tangential Smoothing without the 1/2 attenuation!
				Mesh::Point normal = mesh.normal(vh);

				Lu /= n;
				Lu -= mesh.point(vh);

				// Lu_tangent = Lu - <Lu,n>n 
				Mesh::Point normal_component(Lu[0] * normal[0], Lu[1] * normal[1], Lu[2] * normal[2]);
				Mesh::Point Lu_tangent = Lu - normal * normal_component;
				// v' = v + Lu_tagent(v)
				mesh.property(update, vh) = mesh.point(vh) + Lu_tangent; // defer vertex update
			}
		}

		for (v_it=mesh.vertices_begin(); v_it!=v_end; ++v_it)
			if (!mesh.is_boundary(v_it))
				mesh.point(v_it) = mesh.property(update, v_it);
	}
}


//-----------------------------------------------------------------------------
