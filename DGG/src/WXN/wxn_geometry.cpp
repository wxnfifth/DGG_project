#include "wxn_geometry.h"
#include "wxnMath.h"
#include "wxnTime.h"
#include "ich\Point3D.h"

void TexturedFaces::calculateCoveredFaces()
{
    fprintf(stderr,"coverted vertex %d\n" , planar_coordinates.size());
    set<int> textured_faces_set;
    textured_faces_set.clear();
    for(map<int,Point2D>::const_iterator itr = planar_coordinates.begin(); itr != planar_coordinates.end(); ++itr) {
        const vector<pair<int, double> >& neighbors = model_pt->Neigh(itr->first);
        for(vector<pair<int, double> >::const_iterator itr_neigh = neighbors.begin(); itr_neigh != neighbors.end();++itr_neigh){
            int face_id = model_pt->GetNeighborFaceIndexFromEdge(itr_neigh->first);
            textured_faces_set.insert(face_id);
        }
    }
    for (set<int>::iterator itr = textured_faces_set.begin(); itr != textured_faces_set.end();) {
        const CBaseModel::CFace& face = model_pt->Face(*itr);
        bool delete_flag = false;
        for (int i = 0; i < 3; ++i) {
            if (planar_coordinates.find(face[i]) == planar_coordinates.end()) {
                delete_flag = true;
                break;
            }
        }
        if (delete_flag) {
            itr = textured_faces_set.erase(itr);
        } else {
            itr++;
        }
    }
    box_of_faces.clear();
    textured_faces.clear();
    for (set<int>::iterator itr = textured_faces_set.begin(); itr != textured_faces_set.end(); ++itr) {
        textured_faces.push_back(*itr);
        Box b;
        for(int k = 0; k < 3; ++k) {
            int vert_id = model_pt->Face(*itr)[k];
            b.update( planar_coordinates[vert_id] );
        }
        box_of_faces.push_back(b);
    }
}

void TexturedFaces::mapPointToModel(const Point2D& p , PointOnFace& point_on_face) const {

    vector<int> face_ids;
    //printf("face_ids :\n");
    //for(vector<Box>::iterator itr = box_of_faces.begin(); itr != box_of_faces.end();++itr) {
    ElapasedTime temp_time;
    for(int i = 0; i < box_of_faces.size();++i) {
        if( box_of_faces[i].inBox(p) ) {
            face_ids.push_back( textured_faces[i] );
            //printf("%d " , face_ids[i]);
        }
    }
    //temp_time.printTime("time 1");
    temp_time.start();
    //printf("%d " , face_ids.size());
    //    for(vector<int>::const_iterator face_itr = textured_faces.begin(); face_itr != textured_faces.end(); ++face_itr) {
    for( vector<int>::iterator face_itr = face_ids.begin(); face_itr != face_ids.end(); ++face_itr) {
        const CBaseModel::CFace& face = model_pt->Face(*face_itr);
        Point2D vertex_2d[3];
        CPoint3D vertex_3d[3];
        for (int k = 0; k < 3; ++k) {
            vertex_3d[k] = model_pt->Vert(face[k]);
            map<int,Point2D>::const_iterator texture_iterator = planar_coordinates.find(face[k]);
            assert(texture_iterator != planar_coordinates.end());
            vertex_2d[k] = Point2D(texture_iterator->second.x , texture_iterator->second.y);
        }
        CPoint3D barycentric = pointIn2DTriangleToBarycentric(vertex_2d,p);
        if (barycentric.x >= 0 &&
            barycentric.y >= 0 &&
            barycentric.z >= 0) {
                point_on_face.v = barycentric.x * vertex_3d[0] + barycentric.y * vertex_3d[1] + barycentric.z * vertex_3d[2];
                point_on_face.face_id_ = *face_itr;
                //temp_time.printTime("    time 2\n\n");
                return;
        }
    }
    
    dumpTexturedFaces();
    p.print("Point2D p");
    fprintf(stderr,"face_ids : ");
    for (int i = 0; i < face_ids.size();++i) {
        fprintf(stderr,"%d " , face_ids[i]);
    }
    fprintf(stderr,"\n");
    assert(false);
    //}
}

void TexturedFaces::dumpTexturedFaces() const
{
    string file_name = "textured_faces.obj";
    FILE* file_out = fopen(file_name.c_str() ,"w");
    int cnt = 1;
    for (vector<int>::const_iterator itr = textured_faces.begin(); itr != textured_faces.end(); ++itr) {
        for(int k = 0; k < 3; ++k) {
            fprintf(file_out , "v %lf %lf 0\n" , planar_coordinates.find(model_pt->Face(*itr)[k])->second.x , 
                                                   planar_coordinates.find(model_pt->Face(*itr)[k])->second.y );
        }
        fprintf(file_out , "f %d %d %d\n" , cnt , cnt + 1 , cnt + 2);
        cnt += 3;
    }
    //fprintf(file_out , "v 0 0 0\n" );
    //fprintf(file_out , "v 0 1 0\n" );
    //fprintf(file_out , "v 1 1 0\n" );
    //fprintf(file_out , "v 1 0 0\n" );
    fprintf(file_out , "v -0.5 -0.5 0\n");
    fprintf(file_out , "v -0.5 0.5 0\n");
    fprintf(file_out , "v 0.5 0.5 0\n");
    fprintf(file_out , "v 0.5 -0.5 0\n");
    fprintf(file_out , "f %d %d %d\n" , cnt , cnt + 1 , cnt + 2 );
    fprintf(file_out , "f %d %d %d\n" , cnt , cnt + 2 , cnt + 3 );
    fclose(file_out);
}


void TexturedFaces::mapModelPoint2Plane(const PointOnFace& p_on_face , Point2D& point_2D) const
{
    CPoint3D p_3d[3];
    int v_id[3];
    for(int k = 0; k < 3; ++k) {
        v_id[k] = model_pt->Face(p_on_face.face_id_)[k];
        p_3d[k] = (model_pt->Vert(v_id[k]));
    }

    CPoint3D barycentric = pointIn3DTriangleToBarycentric(p_3d , p_on_face.v);
    point_2D = barycentric.x * planar_coordinates.find(v_id[0])->second + 
               barycentric.y * planar_coordinates.find(v_id[1])->second + 
               barycentric.z * planar_coordinates.find(v_id[2])->second;
}


void TexturedFaces::RenderFaces(const GLenum& mode) const
{
    glPushMatrix();
    GLint shadeModel;
    glGetIntegerv(GL_SHADE_MODEL, &shadeModel);

    const vector<int>& faces = textured_faces;
    for(vector<int>::const_iterator itr = faces.begin(); itr != faces.end();++itr){
        int face_id = *itr;
        int out_boundary_vertex = 0;
        const Point2D* p_texes[3];
        for (int j = 0; j < 3; ++j){
            const Point2D& p_tex = planar_coordinates.find(model_pt->Face(face_id)[j])->second;
            p_texes[j] = &p_tex;
            if (p_tex.x <= 0 || p_tex.y <= 0 || p_tex.x >=1 || p_tex.y >=1 ) {
                out_boundary_vertex++;
                //break;
            }
        }
        if( out_boundary_vertex > 0 ) continue;

        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; ++j)
        {			
            CPoint3D pt = model_pt->Vert(model_pt->Face(face_id)[j]) +  model_pt->Normal(model_pt->Face(face_id)[j]) * RateOfNormalShift ;
            const CPoint3D &normal = model_pt->Normal(model_pt->Face(face_id)[j]);
            glNormal3f((float)normal.x, (float)normal.y, (float)normal.z);
            glTexCoord2d( p_texes[j]->x , p_texes[j]->y);
            glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
        }
        glEnd();
    }		
    glPopMatrix();

}


void TexturedFacesWithControlPoints::addControlPoint(const PointOnFace& hit_point , const bool& flag_mouse_first_point)
{
    TexturedPoint p;
    p.p_3d = hit_point;
    mapModelPoint2Plane(p.p_3d , p.p_2d);
    if (flag_mouse_first_point) {
        temp_control_points.clear();
        temp_control_points.push_back(p);
        temp_control_edges.clear();
    } else {
        temp_control_edges.push_back(ControlEdge(temp_control_points.size() -1 , temp_control_points.size()));
        temp_control_points.push_back(p);
    }
}

void TexturedFacesWithControlPoints::endAddControlPoint()
{
    temp_control_edges.push_back(ControlEdge(temp_control_points.size() -1 , 0));
    //origin_control_points.insert(origin_control_points.end() , temp_control_points.begin() , temp_control_points.end());
    for(int i = 0; i < temp_control_points.size();++i) {
        origin_control_points.push_back(temp_control_points[i].Point_2D());
    }
    for(int i = 0; i < temp_control_edges.size();++i){
        temp_control_edges[i].e1 += control_points.size();
        temp_control_edges[i].e2 += control_points.size();
    }
    origin_control_edges.insert(origin_control_edges.end() , temp_control_edges.begin() , temp_control_edges.end());
    control_points.insert(control_points.end() , temp_control_points.begin() , temp_control_points.end());
    control_edges.insert(control_edges.end() , temp_control_edges.begin() , temp_control_edges.end());
    UpdateControlPoints();

    temp_control_points.clear();
    temp_control_edges.clear();
}

void TexturedFacesWithControlPoints::SetControlPoints(const vector<Point2D>& _control_points,
                                                      const vector<int>& _edges
                                                      )                                                     
{
    if( origin_control_points.size() == 0) {
        origin_control_points = _control_points;
        origin_control_edges.clear();
        origin_control_edges.resize(_edges.size()/2);
        for (int i = 0; i < _edges.size(); i += 2) {
            origin_control_edges[i/2] = ControlEdge(_edges[i],_edges[i+1]);
        }
    }
    //control_points = _control_points;
    //control_points_3d.resize(control_points.size());
    control_points.clear();
    control_points.resize(_control_points.size());
    for (int i = 0; i < control_points.size();++i) {
        control_points[i].p_2d = _control_points[i];
    }
    control_edges = origin_control_edges;
    UpdateControlPoints();
}

void TexturedFacesWithControlPoints::UpdateControlPoints()
{
    for (int i = 0; i < control_points.size(); ++i) {
        mapPointToModel(control_points[i].p_2d , control_points[i].p_3d);
    }
    //control_polygons.clear();
    //control_polygons_3d.clear();
    control_edge_lines.clear();
    //control_polygons.resize(control_points.size());
    control_edge_lines.resize(control_edges.size());
    double min_delta_dis = model_pt->GetMaxEdgeLength();
    for (int i = 0; i < control_edges.size(); ++i) {
        const Point2D& p0 = control_points[control_edges[i].e1].Point_2D();
        const Point2D& p1 = control_points[control_edges[i].e2].Point_2D();
        double dis = (p0-p1).length();
        double delta = min_delta_dis / dis;
        int seg_num = (int)(1 / delta) + 1;
        control_edge_lines[i].resize(seg_num-1);
        for (int j = 1; j < seg_num; ++j) {
            Point2D& p = p1 * j * delta + p0 * ( 1 - j * delta );
            control_edge_lines[i][j-1].p_2d = p;
        }
    }

    //control_polygons_3d.resize(control_polygons.size());
    for (int i = 0; i < control_edge_lines.size(); ++i) { 
        //control_polygons_3d[i].resize(control_polygons[i].size());
        for (int j = 0; j < control_edge_lines[i].size(); ++j) {
            mapPointToModel(control_edge_lines[i][j].p_2d , control_edge_lines[i][j].p_3d);
        }
    }
}

void TexturedFacesWithControlPoints::renderControlPoints()const
{
    
    for (int i = 0; i < control_points.size();++i) {
        if (i == select_control_id) {
            renderSphere(control_points[i].Point_3D() ,
                         color_black , control_point_radius);
        } else {
            renderSphere(control_points[i].Point_3D() , 
                         color_red , control_point_radius);
        }
    }
    for (int i = 0; i < temp_control_points.size();++i) {
        renderSphere(temp_control_points[i].Point_3D() , 
                     color_green ,
                     control_point_radius);
        fprintf(stderr,"draw temp_control_points!\n");
    }
    for (int i = 0; i < control_edges.size(); ++i) {
        const PointOnFace& p0 = control_points[control_edges[i].e1].p_3d;
        const PointOnFace& p1 = control_points[control_edges[i].e2].p_3d;
        //CPoint3D pt = model_pt->Vert(model_pt->Face(face_id)[j]) +  model_pt->Normal(model_pt->Face(face_id)[j]) * RateOfNormalShift ;
        glBegin(GL_LINE_STRIP);
        model_pt->ComputeShiftPoint(p0.v , model_pt->Face(p0.face_id_) ).Show();
        for (int j = 0; j < control_edge_lines[i].size(); ++j) {
            const PointOnFace& p = control_edge_lines[i][j].p_3d;
            model_pt->ComputeShiftPoint(p.v , model_pt->Face(p.face_id_) ).Show();
            
        }
        model_pt->ComputeShiftPoint(p1.v , model_pt->Face(p1.face_id_) ).Show();
        glEnd();
    }
}




void TexturedFacesWithControlPoints::RenderFaces(const GLenum& mode) const
{
    glPushMatrix();
    GLint shadeModel;
    glGetIntegerv(GL_SHADE_MODEL, &shadeModel);

    const vector<int>& faces = textured_faces;
    for(vector<int>::const_iterator itr = faces.begin(); itr != faces.end();++itr){
        int face_id = *itr;
        int out_boundary_vertex = 0;
        const Point2D* p_texes[3];
        const Point2D* p_planar_coord[3];
        for (int j = 0; j < 3; ++j){
            const Point2D& p_tex = deformed_coordinates.find(model_pt->Face(face_id)[j])->second;
            p_texes[j] = &p_tex;
            const Point2D& p_planar = planar_coordinates.find(model_pt->Face(face_id)[j])->second;
            p_planar_coord[j]  = &p_planar;
            //if (p_tex.x <= 0 || p_tex.y <= 0 || p_tex.x >=1 || p_tex.y >=1 ) {
            if( p_planar.x <= 0 || p_planar.y <= 0 || p_planar.x >= 1 || p_planar.y >=1 ) {
                out_boundary_vertex++;
                //break;
            }
        }
        if( out_boundary_vertex > 0 ) continue;

        glBegin(GL_TRIANGLES);
        for (int j = 0; j < 3; ++j)
        {			
            CPoint3D pt = model_pt->Vert(model_pt->Face(face_id)[j]) +  model_pt->Normal(model_pt->Face(face_id)[j]) * RateOfNormalShift ;
            const CPoint3D &normal = model_pt->Normal(model_pt->Face(face_id)[j]);
            glNormal3f((float)normal.x, (float)normal.y, (float)normal.z);
            glTexCoord2d( p_texes[j]->x , p_texes[j]->y);
            glVertex3f((float)pt.x, (float)pt.y, (float)pt.z);
        }
        glEnd();
    }		
    glPopMatrix();

}

void TexturedFacesWithControlPoints::selectControlPoint(const CPoint3D& p)
{
    int current_select_control_id = getClosestControlPoint(p);
    select_control_id = current_select_control_id;
}

int TexturedFacesWithControlPoints::getClosestControlPoint(const CPoint3D&p)
{
    int closest_control_point = -1;
    double min_dis_sqr = 1e10;
    for (int i = 0; i < control_points.size(); ++i) {
        double temp_dis_sqr = (p - control_points[i].Point_3D()).LenSqr();
        if( temp_dis_sqr < control_point_radius * control_point_radius
            && temp_dis_sqr < min_dis_sqr) {
            closest_control_point = i;
            min_dis_sqr = temp_dis_sqr;
        }
    }
    return closest_control_point;
}


void  TexturedFacesWithControlPoints::moveSelection(const PointOnFace& first_p ,
                                                    const PointOnFace& last_p)
{
    Point2D first_point_2D;
    Point2D last_point_2D;
    mapModelPoint2Plane(first_p,first_point_2D);
    mapModelPoint2Plane(last_p,last_point_2D);
    Point2D delta = last_point_2D - first_point_2D;
    if( select_control_id >= 0 ){
        assert(select_control_id < control_points.size());
        control_points[select_control_id].p_2d += delta;
    }
    UpdateControlPoints();
}
