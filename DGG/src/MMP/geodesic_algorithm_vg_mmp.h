#ifndef GEODESIC_ALGORITHM_EXACT
#define GEODESIC_ALGORITHM_EXACT

#include "geodesic_algorithm_mmp_basement.h"
#include <set>
#include <algorithm>

namespace geodesic{

class GeodesicAlgorithmVGMMP : public GeodesicAlgorithmMMPBasement
{
public:
	GeodesicAlgorithmVGMMP(geodesic::Mesh* mesh):
	  	GeodesicAlgorithmMMPBasement(mesh)
	{
	};	

	~GeodesicAlgorithmVGMMP(){};

	void propagate(std::vector<SurfacePoint>& sources,
   				   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
				   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

	void propagate_local(std::vector<SurfacePoint>& sources,
    int fixed_k, std::map<int,double>& dest_verts,
   				   double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
				   std::vector<SurfacePoint>* stop_points = NULL); //or after ensuring that all the stop_points are covered

  void propagate_vg(std::vector<SurfacePoint>& sources,
    double eps_vg, std::map<int,double>& dest_verts,
    double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
    std::vector<SurfacePoint>* stop_points = NULL
    ); //or after ensuring that all the stop_points are covered

  void propagate_dgg(std::vector<SurfacePoint>& sources,
    double eps_vg, std::map<int,double>& dest_verts,
    double max_propagation_distance = GEODESIC_INF,			//propagation algorithm stops after reaching the certain distance from the source
    std::vector<SurfacePoint>* stop_points = NULL
    ); //or after ensuring that all the stop_points are covered


private:
	typedef std::set<interval_pointer, Interval> IntervalQueue;

	void update_list_and_queue(list_pointer list,
							   IntervalWithStop* candidates,	//up to two candidates
							   unsigned num_candidates);

	bool check_stop_conditions(unsigned& index);

	void initialize_propagation_data();		

	bool erase_from_queue(interval_pointer p);

	IntervalQueue m_queue;	//interval queue
};

inline bool GeodesicAlgorithmVGMMP::erase_from_queue(interval_pointer p)
{
	if(p->min() < GEODESIC_INF/10.0)// && p->min >= queue->begin()->first)
	{
		assert(m_queue.count(p)<=1);			//the set is unique

		IntervalQueue::iterator it = m_queue.find(p);

		if(it != m_queue.end())
		{
			m_queue.erase(it);
			return true;
		}
	}

	return false;
}

inline void GeodesicAlgorithmVGMMP::initialize_propagation_data()
{
	clear();

	IntervalWithStop candidate;
	std::vector<edge_pointer> edges_visible_from_source;
	for(unsigned i=0; i<m_sources.size(); ++i)		//for all edges adjacent to the starting vertex			
	{
		SurfacePoint* source = &m_sources[i];
		
		edges_visible_from_source.clear();
		list_edges_visible_from_source(source->base_element(), 
									   edges_visible_from_source);
		
		for(unsigned j=0; j<edges_visible_from_source.size(); ++j)
		{
			edge_pointer e = edges_visible_from_source[j];
			candidate.initialize(e, source, i);
            candidate.stop() = e->length();
			candidate.compute_min_distance(candidate.stop());
			candidate.direction() = Interval::FROM_SOURCE;

			update_list_and_queue(interval_list(e), &candidate, 1);
		}
	}
}

inline void GeodesicAlgorithmVGMMP::propagate_vg(std::vector<SurfacePoint>& sources,
                                                 double eps_vg, std::map<int,double>& dest_verts,
                                                 double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
                                                 std::vector<SurfacePoint>* stop_points
                                                 )  //or after ensuring that all the stop_points are covered
{
  set_stop_conditions(stop_points, max_propagation_distance);
  set_sources(sources);
  initialize_propagation_data();

  clock_t start = clock();

  unsigned satisfied_index = 0;

  m_iterations = 0;		//for statistics
  m_queue_max_size = 0;

  IntervalWithStop candidates[2];

  dest_verts.clear();
  double d_max = GEODESIC_INF;
  //double d_max_origin = GEODESIC_INF;
  double d_max_current = GEODESIC_INF;
  
  double average_2_part1 = 0;
  double average_2_part2 = 0;
  double average_1 = 0;
  double average_2 = 0;
  int average_cnt = 0;

  while(!m_queue.empty())
  {
    m_queue_max_size = std::max((int)m_queue.size(), (int)m_queue_max_size);

    m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();

    unsigned const check_period = 10;
    if(++m_iterations % check_period == 0)		//check if we covered all required vertices
    {
      if (check_stop_conditions(satisfied_index))
      {
        break;
      }
    }

    interval_pointer min_interval = *m_queue.begin();
    m_queue.erase(m_queue.begin());
    edge_pointer edge = min_interval->edge();
    list_pointer list = interval_list(edge);

    assert(min_interval->d() < GEODESIC_INF);



    int v0 = edge->v0()->id();
    int v1 = edge->v1()->id();
    //printf("v0 %d v1 %d\n" , v0 , v1);
    if (false) {
    geodesic::SurfacePoint dest_v0 = geodesic::SurfacePoint(&mesh()->vertices()[v0]);
    double dis_v0 = 0;
    best_source(dest_v0, dis_v0);
    geodesic::SurfacePoint dest_v1 = geodesic::SurfacePoint(&mesh()->vertices()[v1]);
    double dis_v1 = 0;
    best_source(dest_v0, dis_v1);
    }
    //printf("v0 best dis %lf v1 best dis %lf\n" , dis_v0, dis_v1);
    double distance_v0 = interval_list(edge)->signal(0.0);
    double distance_v1 = interval_list(edge)->signal(edge->length());
    //if (distance_v0 < min_interval->min()) {
    //  dest_verts[v0] = distance_v0;
    //}
    //if (distance_v1 < min_interval->min()) {
    //  dest_verts[v1] = distance_v1;
    //}
    if (distance_v0 < GEODESIC_INF) {
      dest_verts[v0] = distance_v0;
    }
    if (distance_v1 < GEODESIC_INF) {
    dest_verts[v1] = distance_v1;
    }

    //printf("min_interval %lf d_max %lf\n" , min_interval->min(), d_max);
    //update d_max

    //double tmp_d_max = edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min()) + min_interval->min();
      
    //d_max_origin = min(d_max_origin, edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min()));
    if (false) {
    double x1 = min_interval->start();
    double x2 = min_interval->stop();
    double a = x1;
    double b = 1-x2;
    double tmp_t = (1 - fabs(a-b) ) / 2.0;
    //double tmp_t = 0.5;
    if (edge->length() / min_interval->min() < 4 * sqrt(eps_vg))
    {
      double tmp_edge_dist_new = fabs(min_interval->pseudo_y() * tmp_t * edge->length() / sqrt( (min_interval->pseudo_x() - tmp_t * edge->length() ) * (min_interval->pseudo_x() - tmp_t * edge->length() ) + min_interval->pseudo_y() * min_interval->pseudo_y() ) ); 

      double tmp_d_max = tmp_edge_dist_new * tmp_edge_dist_new / (2.0 * eps_vg * min_interval->min() ) + min_interval->min() + edge->length() * 0.5;
      //tmp_d_max = tmp_edge_dist_new * 4.0/ (3.0 * eps_vg);
      d_max_current = min(d_max_current, tmp_d_max); 
    }
    average_cnt++;
    if (min_interval->min() > d_max_current - 1e-6) continue;
    //cout << "tmp_d_max " << d_max_current << "origin_xx " << d_max_origin << " interval length " << min_interval->min() << " d " << min_interval->d() << "\n";
    } 

    if (true) {
      double tmp_e  = 0.5 * edge->length();
      //tmp_e = min(tmp_e, min_interval->start());
      //tmp_e = min(tmp_e, edge->length() - min_interval->stop());

      double tmp_d_max = GEODESIC_INF;
      if (min_interval->min() - tmp_e > 0) {
        tmp_d_max = tmp_e * tmp_e / (2.0 * eps_vg * (min_interval->min() - tmp_e)) + tmp_e;
      }
      d_max_current = min(d_max_current, tmp_d_max); 
      // f fprintf(stderr,stderr,"e %lf d_max %lf current %lf\n" , tmp_e, tmp_d_max, min_interval->min());
      if (min_interval->min() > d_max_current - 1e-6) continue;
    }

    if (false) {
      double tmp_e = 0;
      double b0 = min_interval->start();
      double b1 = min_interval->stop();
      double M = 0.5 * edge->length();
      Point2D A(0,0);
      Point2D B(edge->length(),0);
      Point2D S(min_interval->pseudo_x(),min_interval->pseudo_y());
      Point2D Pb0(b0,0);
      Point2D Pb1(b1,0);
      Point2D PM(M,0);
      double d0 = S.distance(Pb0.xy());
      double d1 = S.distance(Pb1.xy());
      double Ab1 = A.distance(Pb1.xy());
      double b0B = B.distance(Pb0.xy());
      if (M > b1 && M < edge->length()) {
        double a = S.distance(A.xy());
        double b = d1;//S.distance(Pb1.xy());
        double c = Ab1;//A.distance(Pb1.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      } else if (M > 0 && M < b0) {
        double a = S.distance(B.xy());
        double b = S.distance(Pb0.xy());
        double c = b0B;//B.distance(Pb0.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      } else {//if (M > b0 && M < b1) {
        double a = S.distance(A.xy());
        double b = S.distance(PM.xy());
        double c = A.distance(PM.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      }
      double tmp_d_max = GEODESIC_INF;
      if (min_interval->min() - M > 0) {
        tmp_d_max = tmp_e * tmp_e / (2.0 * eps_vg * (min_interval->min() - M)) + max(d0,d1) + min(min(Ab1,b0B),M);
      }
      d_max_current = min(d_max_current, tmp_d_max); 
      // fprintf(stderr,"e %lf d_max %lf current %lf\n" , tmp_e, tmp_d_max, min_interval->min());
      if (min_interval->min() > d_max_current - 1e-6) continue;
    }
  
#if 0
    if ( fabs(min_interval->min()) > 1e-10) {


      double tmp_d_max = edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min()) + min_interval->min();
      //double tmp_d_max = edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min());
      average_1 += edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min());
      average_2 += tmp_d_max;
      average_2_part1 += edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min());
      average_2_part2 += min_interval->min();
      average_cnt ++;

      //double e = min_interval->stop() - min_interval->start();
      //double tmp_d_max = e * e / (2.0 * eps_vg * min_interval->min());
      d_max = min(tmp_d_max, d_max);
      if (min_interval->min() > d_max - 1e-6) continue;
 
    }

#endif
    //vertex_pointer v = m_stop_vertices[index].first;
    //if (dest_verts.size() >= fixed_k) {
    //  break;
    //}

    bool const first_interval = min_interval->start() == 0.0;
    //bool const last_interval = min_interval->stop() == edge->length();
    bool const last_interval = min_interval->next() == NULL;

    //bool const turn_left = edge->v0()->saddle_or_boundary();
    //bool const turn_right = edge->v1()->saddle_or_boundary();
    bool const turn_left = false;
    bool const turn_right = false;

    for(unsigned i=0; i<edge->adjacent_faces().size(); ++i)		//two possible faces to propagate
    {
      if(!edge->is_boundary())		//just in case, always propagate boundary edges
      {
        if((i == 0 && min_interval->direction() == Interval::FROM_FACE_0) ||
          (i == 1 && min_interval->direction() == Interval::FROM_FACE_1))
        {
          continue;
        }
      }

      face_pointer face = edge->adjacent_faces()[i];			//if we come from 1, go to 2
      edge_pointer next_edge = face->next_edge(edge,edge->v0());

      unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), 
        min_interval->pseudo_y(), 
        min_interval->d(),		//parameters of the interval
        min_interval->start(), 
        min_interval->stop(),		//start/end of the interval
        face->vertex_angle(edge->v0()),	//corner angle
        next_edge->length(),		//length of the new edge
        first_interval,		//if it is the first interval on the edge
        last_interval,
        turn_left,
        turn_right,
        candidates);		//if it is the last interval on the edge
      bool propagate_to_right = true;

      if(num_propagated)
      {
        if(candidates[num_propagated-1].stop() != next_edge->length()) 
        {
          propagate_to_right = false;
        }

        bool const invert = next_edge->v0()->id() != edge->v0()->id(); //if the origins coinside, do not invert intervals

        construct_propagated_intervals(invert,		//do not inverse 
          next_edge, 
          face,
          candidates,
          num_propagated,
          min_interval);

        update_list_and_queue(interval_list(next_edge), 
          candidates, 
          num_propagated);
      }

      if(propagate_to_right)
      {
        //propogation to the right edge
        double length = edge->length();
        next_edge = face->next_edge(edge,edge->v1());

        num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), 
          min_interval->pseudo_y(), 
          min_interval->d(),		//parameters of the interval
          length - min_interval->stop(), 
          length - min_interval->start(),		//start/end of the interval
          face->vertex_angle(edge->v1()),	//corner angle
          next_edge->length(),		//length of the new edge
          last_interval,		//if it is the first interval on the edge
          first_interval,
          turn_right,
          turn_left,
          candidates);		//if it is the last interval on the edge

        if(num_propagated)
        {
          bool const invert = next_edge->v0()->id() != edge->v1()->id();		//if the origins coinside, do not invert intervals

          construct_propagated_intervals(invert,		//do not inverse 
            next_edge, 
            face,
            candidates,
            num_propagated,
            min_interval);

          update_list_and_queue(interval_list(next_edge), 
            candidates, 
            num_propagated);
        }
      }
    } 
  } 

  m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
  clock_t stop = clock();
  m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;

  //printf("average_cnt %d\n" , average_cnt);
  average_1 /= average_cnt;
  average_2 /= average_cnt;
  average_2_part1 /= average_cnt;
  average_2_part2 /= average_cnt;

  //printf("average_1 %lf average_2 %lf average_2_part1 %lf average_2_part2 %lf\n" , average_1, average_2, average_2_part1, average_2_part2);

}


inline void GeodesicAlgorithmVGMMP::propagate_dgg(std::vector<SurfacePoint>& sources,
                                                 double eps_vg, std::map<int,double>& dest_verts,
                                                 double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
                                                 std::vector<SurfacePoint>* stop_points
                                                 )  //or after ensuring that all the stop_points are covered
{
  set_stop_conditions(stop_points, max_propagation_distance);
  set_sources(sources);
  initialize_propagation_data();

  clock_t start = clock();

  unsigned satisfied_index = 0;

  m_iterations = 0;		//for statistics
  m_queue_max_size = 0;

  IntervalWithStop candidates[2];

  dest_verts.clear();
  double d_max = GEODESIC_INF;
  //double d_max_origin = GEODESIC_INF;
  double d_max_current = GEODESIC_INF;
  
  double average_2_part1 = 0;
  double average_2_part2 = 0;
  double average_1 = 0;
  double average_2 = 0;
  int average_cnt = 0;

  while(!m_queue.empty())
  {
    m_queue_max_size = std::max((int)m_queue.size(), (int)m_queue_max_size);

    m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();

    unsigned const check_period = 10;
    if(++m_iterations % check_period == 0)		//check if we covered all required vertices
    {
      if (check_stop_conditions(satisfied_index))
      {
        break;
      }
    }

    interval_pointer min_interval = *m_queue.begin();
    m_queue.erase(m_queue.begin());
    edge_pointer edge = min_interval->edge();
    list_pointer list = interval_list(edge);

    assert(min_interval->d() < GEODESIC_INF);



    int v0 = edge->v0()->id();
    int v1 = edge->v1()->id();
    //printf("v0 %d v1 %d\n" , v0 , v1);
    geodesic::SurfacePoint dest_v0 = geodesic::SurfacePoint(&mesh()->vertices()[v0]);
    double dis_v0 = 0;
    best_source(dest_v0, dis_v0);
    geodesic::SurfacePoint dest_v1 = geodesic::SurfacePoint(&mesh()->vertices()[v1]);
    double dis_v1 = 0;
    best_source(dest_v0, dis_v1);
    //printf("v0 best dis %lf v1 best dis %lf\n" , dis_v0, dis_v1);
    double distance_v0 = interval_list(edge)->signal(0.0);
    double distance_v1 = interval_list(edge)->signal(edge->length());
    //if (distance_v0 < min_interval->min()) {
    //  dest_verts[v0] = distance_v0;
    //}
    //if (distance_v1 < min_interval->min()) {
    //  dest_verts[v1] = distance_v1;
    //}
    if (distance_v0 < GEODESIC_INF) {
      dest_verts[v0] = distance_v0;
    }
    if (distance_v1 < GEODESIC_INF) {
    dest_verts[v1] = distance_v1;
    }




    //printf("min_interval %lf d_max %lf\n" , min_interval->min(), d_max);
    //update d_max

    //double tmp_d_max = edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min()) + min_interval->min();
      
    //d_max_origin = min(d_max_origin, edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min()));
    if (false) {
    double x1 = min_interval->start();
    double x2 = min_interval->stop();
    double a = x1;
    double b = 1-x2;
    double tmp_t = (1 - fabs(a-b) ) / 2.0;
    //double tmp_t = 0.5;
    if (edge->length() / min_interval->min() < 4 * sqrt(eps_vg))
    {
      double tmp_edge_dist_new = fabs(min_interval->pseudo_y() * tmp_t * edge->length() / sqrt( (min_interval->pseudo_x() - tmp_t * edge->length() ) * (min_interval->pseudo_x() - tmp_t * edge->length() ) + min_interval->pseudo_y() * min_interval->pseudo_y() ) ); 

      double tmp_d_max = tmp_edge_dist_new * tmp_edge_dist_new / (2.0 * eps_vg * min_interval->min() ) + min_interval->min() + edge->length() * 0.5;
      //tmp_d_max = tmp_edge_dist_new * 4.0/ (3.0 * eps_vg);
      d_max_current = min(d_max_current, tmp_d_max); 
    }
    average_cnt++;
    if (min_interval->min() > d_max_current - 1e-6) continue;
    //cout << "tmp_d_max " << d_max_current << "origin_xx " << d_max_origin << " interval length " << min_interval->min() << " d " << min_interval->d() << "\n";
    } 

    if (true) {
      double tmp_e  = 0.5 * edge->length();
      //tmp_e = min(tmp_e, min_interval->start());
      //tmp_e = min(tmp_e, edge->length() - min_interval->stop());

      double tmp_d_max = GEODESIC_INF;
      if (min_interval->min() - tmp_e > 0) {
        tmp_d_max = tmp_e * tmp_e / (2.0 * eps_vg * (min_interval->min() - tmp_e)) +tmp_e;
      }
      d_max_current = min(d_max_current, tmp_d_max); 
      // fprintf(stderr,"e %lf d_max %lf current %lf\n" , tmp_e, tmp_d_max, min_interval->min());
      if (min_interval->min() > d_max_current - 1e-6) continue;
    }

    if (false) {
      double tmp_e = 0;
      double b0 = min_interval->start();
      double b1 = min_interval->stop();
      double M = 0.5 * edge->length();
      Point2D A(0,0);
      Point2D B(edge->length(),0);
      Point2D S(min_interval->pseudo_x(),min_interval->pseudo_y());
      Point2D Pb0(b0,0);
      Point2D Pb1(b1,0);
      Point2D PM(M,0);
      double d0 = S.distance(Pb0.xy());
      double d1 = S.distance(Pb1.xy());
      double Ab1 = A.distance(Pb1.xy());
      double b0B = B.distance(Pb0.xy());
      if (M > b1 && M < edge->length()) {
        double a = S.distance(A.xy());
        double b = d1;//S.distance(Pb1.xy());
        double c = Ab1;//A.distance(Pb1.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      } else if (M > 0 && M < b0) {
        double a = S.distance(B.xy());
        double b = S.distance(Pb0.xy());
        double c = b0B;//B.distance(Pb0.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      } else {//if (M > b0 && M < b1) {
        double a = S.distance(A.xy());
        double b = S.distance(PM.xy());
        double c = A.distance(PM.xy());
        double costheta = (a*a+b*b-c*c)/(2.0*a*b);
        if (costheta > 1) costheta = 1;
        if (costheta < -1) costheta = -1;
        double sintheta = sqrt(1-costheta*costheta);
        tmp_e = a * sintheta;
      }
      double tmp_d_max = GEODESIC_INF;
      if (min_interval->min() - M > 0) {
        tmp_d_max = tmp_e * tmp_e / (2.0 * eps_vg * (min_interval->min() - M)) + max(d0,d1) + min(min(Ab1,b0B),M);
      }
      d_max_current = min(d_max_current, tmp_d_max); 
      // fprintf(stderr,"e %lf d_max %lf current %lf\n" , tmp_e, tmp_d_max, min_interval->min());
      if (min_interval->min() > d_max_current - 1e-6) continue;
    }


    
     
#if 0
    if ( fabs(min_interval->min()) > 1e-10) {


      double tmp_d_max = edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min()) + min_interval->min();
      //double tmp_d_max = edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min());
      average_1 += edge->length() * edge->length() / (2.0 * eps_vg * min_interval->min());
      average_2 += tmp_d_max;
      average_2_part1 += edge->length() * edge->length() / (8.0 * eps_vg * min_interval->min());
      average_2_part2 += min_interval->min();
      average_cnt ++;

      //double e = min_interval->stop() - min_interval->start();
      //double tmp_d_max = e * e / (2.0 * eps_vg * min_interval->min());
      d_max = min(tmp_d_max, d_max);
      if (min_interval->min() > d_max - 1e-6) continue;
 
    }

#endif
    //vertex_pointer v = m_stop_vertices[index].first;



    //if (dest_verts.size() >= fixed_k) {
    //  break;
    //}



    bool const first_interval = min_interval->start() == 0.0;
    //bool const last_interval = min_interval->stop() == edge->length();
    bool const last_interval = min_interval->next() == NULL;

    //bool const turn_left = edge->v0()->saddle_or_boundary();
    //bool const turn_right = edge->v1()->saddle_or_boundary();
    bool const turn_left = false;
    bool const turn_right = false;


    for(unsigned i=0; i<edge->adjacent_faces().size(); ++i)		//two possible faces to propagate
    {
      if(!edge->is_boundary())		//just in case, always propagate boundary edges
      {
        if((i == 0 && min_interval->direction() == Interval::FROM_FACE_0) ||
          (i == 1 && min_interval->direction() == Interval::FROM_FACE_1))
        {
          continue;
        }
      }

      face_pointer face = edge->adjacent_faces()[i];			//if we come from 1, go to 2
      edge_pointer next_edge = face->next_edge(edge,edge->v0());

      unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), 
        min_interval->pseudo_y(), 
        min_interval->d(),		//parameters of the interval
        min_interval->start(), 
        min_interval->stop(),		//start/end of the interval
        face->vertex_angle(edge->v0()),	//corner angle
        next_edge->length(),		//length of the new edge
        first_interval,		//if it is the first interval on the edge
        last_interval,
        turn_left,
        turn_right,
        candidates);		//if it is the last interval on the edge
      bool propagate_to_right = true;

      if(num_propagated)
      {
        if(candidates[num_propagated-1].stop() != next_edge->length()) 
        {
          propagate_to_right = false;
        }

        bool const invert = next_edge->v0()->id() != edge->v0()->id(); //if the origins coinside, do not invert intervals

        construct_propagated_intervals(invert,		//do not inverse 
          next_edge, 
          face,
          candidates,
          num_propagated,
          min_interval);

        update_list_and_queue(interval_list(next_edge), 
          candidates, 
          num_propagated);
      }

      if(propagate_to_right)
      {
        //propogation to the right edge
        double length = edge->length();
        next_edge = face->next_edge(edge,edge->v1());

        num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), 
          min_interval->pseudo_y(), 
          min_interval->d(),		//parameters of the interval
          length - min_interval->stop(), 
          length - min_interval->start(),		//start/end of the interval
          face->vertex_angle(edge->v1()),	//corner angle
          next_edge->length(),		//length of the new edge
          last_interval,		//if it is the first interval on the edge
          first_interval,
          turn_right,
          turn_left,
          candidates);		//if it is the last interval on the edge

        if(num_propagated)
        {
          bool const invert = next_edge->v0()->id() != edge->v1()->id();		//if the origins coinside, do not invert intervals

          construct_propagated_intervals(invert,		//do not inverse 
            next_edge, 
            face,
            candidates,
            num_propagated,
            min_interval);

          update_list_and_queue(interval_list(next_edge), 
            candidates, 
            num_propagated);
        }
      }
    } 
  } 

  m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
  clock_t stop = clock();
  m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;

  //printf("average_cnt %d\n" , average_cnt);
  average_1 /= average_cnt;
  average_2 /= average_cnt;
  average_2_part1 /= average_cnt;
  average_2_part2 /= average_cnt;

  //printf("average_1 %lf average_2 %lf average_2_part1 %lf average_2_part2 %lf\n" , average_1, average_2, average_2_part1, average_2_part2);

}

inline void GeodesicAlgorithmVGMMP::propagate_local(std::vector<SurfacePoint>& sources,
                                                    int fixed_k, std::map<int,double>& dest_verts,
                                                    double max_propagation_distance ,		
                                                    //propagation algorithm stops after reaching the certain distance from the source
                                                    std::vector<SurfacePoint>* stop_points ) 
                                                    //or after ensuring that all the stop_points are covered
{
  set_stop_conditions(stop_points, max_propagation_distance);
  set_sources(sources);
  initialize_propagation_data();

  clock_t start = clock();

  unsigned satisfied_index = 0;

  m_iterations = 0;		//for statistics
  m_queue_max_size = 0;

  IntervalWithStop candidates[2];

  dest_verts.clear();

  while(!m_queue.empty())
  {
    m_queue_max_size = std::max((int)m_queue.size(), (int)m_queue_max_size);

    m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();

    unsigned const check_period = 10;
    if(++m_iterations % check_period == 0)		//check if we covered all required vertices
    {
      if (check_stop_conditions(satisfied_index))
      {
        break;
      }
    }

    interval_pointer min_interval = *m_queue.begin();



    m_queue.erase(m_queue.begin());
    edge_pointer edge = min_interval->edge();
    list_pointer list = interval_list(edge);

    assert(min_interval->d() < GEODESIC_INF);

    //vertex_pointer v = m_stop_vertices[index].first;
    

    int v0 = edge->v0()->id();
    int v1 = edge->v1()->id();
    double distance_v0 = interval_list(edge)->signal(0.0);
    double distance_v1 = interval_list(edge)->signal(edge->length());
    if (distance_v0 < min_interval->min()) {
      dest_verts[v0] = distance_v0;
    }
    if (distance_v1 < min_interval->min()) {
      dest_verts[v1] = distance_v1;
    }
    if (dest_verts.size() >= fixed_k) {
      break;
    }



    bool const first_interval = min_interval->start() == 0.0;
    //bool const last_interval = min_interval->stop() == edge->length();
    bool const last_interval = min_interval->next() == NULL;

    //bool const turn_left = edge->v0()->saddle_or_boundary();
    //bool const turn_right = edge->v1()->saddle_or_boundary();
    bool const turn_left = false;
    bool const turn_right = false;


    for(unsigned i=0; i<edge->adjacent_faces().size(); ++i)		//two possible faces to propagate
    {
      if(!edge->is_boundary())		//just in case, always propagate boundary edges
      {
        if((i == 0 && min_interval->direction() == Interval::FROM_FACE_0) ||
          (i == 1 && min_interval->direction() == Interval::FROM_FACE_1))
        {
          continue;
        }
      }

      face_pointer face = edge->adjacent_faces()[i];			//if we come from 1, go to 2
      edge_pointer next_edge = face->next_edge(edge,edge->v0());

      unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), 
        min_interval->pseudo_y(), 
        min_interval->d(),		//parameters of the interval
        min_interval->start(), 
        min_interval->stop(),		//start/end of the interval
        face->vertex_angle(edge->v0()),	//corner angle
        next_edge->length(),		//length of the new edge
        first_interval,		//if it is the first interval on the edge
        last_interval,
        turn_left,
        turn_right,
        candidates);		//if it is the last interval on the edge
      bool propagate_to_right = true;

      if(num_propagated)
      {
        if(candidates[num_propagated-1].stop() != next_edge->length()) 
        {
          propagate_to_right = false;
        }

        bool const invert = next_edge->v0()->id() != edge->v0()->id(); //if the origins coinside, do not invert intervals

        construct_propagated_intervals(invert,		//do not inverse 
          next_edge, 
          face,
          candidates,
          num_propagated,
          min_interval);

        update_list_and_queue(interval_list(next_edge), 
          candidates, 
          num_propagated);
      }

      if(propagate_to_right)
      {
        //propogation to the right edge
        double length = edge->length();
        next_edge = face->next_edge(edge,edge->v1());

        num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), 
          min_interval->pseudo_y(), 
          min_interval->d(),		//parameters of the interval
          length - min_interval->stop(), 
          length - min_interval->start(),		//start/end of the interval
          face->vertex_angle(edge->v1()),	//corner angle
          next_edge->length(),		//length of the new edge
          last_interval,		//if it is the first interval on the edge
          first_interval,
          turn_right,
          turn_left,
          candidates);		//if it is the last interval on the edge

        if(num_propagated)
        {
          bool const invert = next_edge->v0()->id() != edge->v1()->id();		//if the origins coinside, do not invert intervals

          construct_propagated_intervals(invert,		//do not inverse 
            next_edge, 
            face,
            candidates,
            num_propagated,
            min_interval);

          update_list_and_queue(interval_list(next_edge), 
            candidates, 
            num_propagated);
        }
      }
    } 
  } 

  m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
  clock_t stop = clock();
  m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;


}

inline void GeodesicAlgorithmVGMMP::propagate(std::vector<SurfacePoint>& sources,
   									   double max_propagation_distance,			//propagation algorithm stops after reaching the certain distance from the source
									   std::vector<SurfacePoint>* stop_points)
{
	set_stop_conditions(stop_points, max_propagation_distance);
	set_sources(sources);
	initialize_propagation_data();

	clock_t start = clock();

	unsigned satisfied_index = 0;

	m_iterations = 0;		//for statistics
	m_queue_max_size = 0;

	IntervalWithStop candidates[2];

  
	while(!m_queue.empty())
	{
		m_queue_max_size = std::max((int)m_queue.size(), (int)m_queue_max_size);

    unsigned const check_period = 10;
    if(++m_iterations % check_period == 0)		//check if we covered all required vertices
    {
      if (check_stop_conditions(satisfied_index))
      {
        break;
      }
    }

		interval_pointer min_interval = *m_queue.begin();
		m_queue.erase(m_queue.begin());
		edge_pointer edge = min_interval->edge();
		list_pointer list = interval_list(edge);

		assert(min_interval->d() < GEODESIC_INF);




		bool const first_interval = min_interval->start() == 0.0;
		//bool const last_interval = min_interval->stop() == edge->length();
		bool const last_interval = min_interval->next() == NULL;

		//bool const turn_left = edge->v0()->saddle_or_boundary();
		//bool const turn_right = edge->v1()->saddle_or_boundary();
		bool const turn_left = false;
		bool const turn_right = false;


		for(unsigned i=0; i<edge->adjacent_faces().size(); ++i)		//two possible faces to propagate
		{
			if(!edge->is_boundary())		//just in case, always propagate boundary edges
			{
				if((i == 0 && min_interval->direction() == Interval::FROM_FACE_0) ||
					(i == 1 && min_interval->direction() == Interval::FROM_FACE_1))
				{
					continue;
				}
			}

			face_pointer face = edge->adjacent_faces()[i];			//if we come from 1, go to 2
			edge_pointer next_edge = face->next_edge(edge,edge->v0());

			unsigned num_propagated = compute_propagated_parameters(min_interval->pseudo_x(), 
																	 min_interval->pseudo_y(), 
																	 min_interval->d(),		//parameters of the interval
																	 min_interval->start(), 
																	 min_interval->stop(),		//start/end of the interval
																	 face->vertex_angle(edge->v0()),	//corner angle
																	 next_edge->length(),		//length of the new edge
																	 first_interval,		//if it is the first interval on the edge
																	 last_interval,
																	 turn_left,
																	 turn_right,
																	 candidates);		//if it is the last interval on the edge
			bool propagate_to_right = true;

			if(num_propagated)
			{
				if(candidates[num_propagated-1].stop() != next_edge->length()) 
				{
					propagate_to_right = false;
				}
				
				bool const invert = next_edge->v0()->id() != edge->v0()->id(); //if the origins coinside, do not invert intervals

				construct_propagated_intervals(invert,		//do not inverse 
											 next_edge, 
											 face,
											 candidates,
											 num_propagated,
											 min_interval);
				
				update_list_and_queue(interval_list(next_edge), 
									  candidates, 
									  num_propagated);
			}

			if(propagate_to_right)
			{
									//propogation to the right edge
				double length = edge->length();
				next_edge = face->next_edge(edge,edge->v1());

				num_propagated = compute_propagated_parameters(length - min_interval->pseudo_x(), 
															 min_interval->pseudo_y(), 
															 min_interval->d(),		//parameters of the interval
															 length - min_interval->stop(), 
															 length - min_interval->start(),		//start/end of the interval
															 face->vertex_angle(edge->v1()),	//corner angle
															 next_edge->length(),		//length of the new edge
															 last_interval,		//if it is the first interval on the edge
															 first_interval,
															 turn_right,
															 turn_left,
															 candidates);		//if it is the last interval on the edge

				if(num_propagated)
				{
					bool const invert = next_edge->v0()->id() != edge->v1()->id();		//if the origins coinside, do not invert intervals

					construct_propagated_intervals(invert,		//do not inverse 
												 next_edge, 
												 face,
												 candidates,
												 num_propagated,
												 min_interval);

					update_list_and_queue(interval_list(next_edge), 
									      candidates, 
										  num_propagated);
				}
			}
		} 
	} 

	m_propagation_distance_stopped = m_queue.empty() ? GEODESIC_INF : (*m_queue.begin())->min();
	clock_t stop = clock();
	m_time_consumed = (static_cast<double>(stop)-static_cast<double>(start))/CLOCKS_PER_SEC;

/*	for(unsigned i=0; i<m_edge_interval_lists.size(); ++i)
	{
		list_pointer list = &m_edge_interval_lists[i];
		interval_pointer p = list->first();
		assert(p->start() == 0.0);
		while(p->next())
		{
			assert(p->stop() == p->next()->start());
			assert(p->d() < GEODESIC_INF);
			p = p->next();
		}
	}*/
}


inline bool GeodesicAlgorithmVGMMP::check_stop_conditions(unsigned& index)
{
	double queue_distance = (*m_queue.begin())->min();
	if(queue_distance < stop_distance())
	{
		return false;
	}

	while(index < m_stop_vertices.size())
	{
		vertex_pointer v = m_stop_vertices[index].first;
		edge_pointer edge = v->adjacent_edges()[0];				//take any edge

		double distance = edge->v0()->id() == v->id() ? 
						  interval_list(edge)->signal(0.0) :
						  interval_list(edge)->signal(edge->length());

		if(queue_distance < distance + m_stop_vertices[index].second)
		{
			return false;
		}

		++index;
	}
	return true;
}


inline void GeodesicAlgorithmVGMMP::update_list_and_queue(list_pointer list,
												IntervalWithStop* candidates,	//up to two candidates
												unsigned num_candidates)
{
	assert(num_candidates <= 2);
	//assert(list->first() != NULL);
	edge_pointer edge = list->edge();
	double const local_epsilon = SMALLEST_INTERVAL_RATIO * edge->length(); 

	if(list->first() == NULL) 
	{
		//winNum += num_candidates;

		interval_pointer* p = &list->first();
		IntervalWithStop* first;
		IntervalWithStop* second; 

		if(num_candidates == 1)
		{
			first = candidates;
			second = candidates;
			first->compute_min_distance(first->stop());
		}
		else 
		{	
			if(candidates->start() <= (candidates+1)->start())
			{
				first = candidates;
				second = candidates+1;
			}
			else 
			{
				first = candidates+1;
				second = candidates;
			}
			assert(first->stop() == second->start());

			first->compute_min_distance(first->stop());
			second->compute_min_distance(second->stop());
		}

		if(first->start() > 0.0)
		{
			*p = m_memory_allocator.allocate();
			(*p)->initialize(edge);
			p = &(*p)->next();
		}

		*p = m_memory_allocator.allocate();
		memcpy(*p,first,sizeof(Interval));
		m_queue.insert(*p);
		++winNum;

		if(num_candidates == 2)
		{
			p = &(*p)->next();
			*p = m_memory_allocator.allocate();
			memcpy(*p,second,sizeof(Interval));
			m_queue.insert(*p);
			++winNum;
		}

		if(second->stop() < edge->length())
		{
			p = &(*p)->next();
			*p = m_memory_allocator.allocate();
			(*p)->initialize(edge);
			(*p)->start() = second->stop();
		}
		else
		{
			(*p)->next() = NULL;
		}
		return;
	}

	bool propagate_flag;

	for(unsigned i=0; i<num_candidates; ++i)				//for all new intervals
	{
		IntervalWithStop* q = &candidates[i];
	
		interval_pointer previous = NULL;

		interval_pointer p = list->first();
		assert(p->start() == 0.0);

		while(p != NULL && p->stop() - local_epsilon < q->start())
		{
			p = p->next(); 
		}

		while(p != NULL && p->start() < q->stop() - local_epsilon)			//go through all old intervals
		{
			unsigned const N = intersect_intervals(p, q);								//interset two intervals

			if(N == 1)			
			{
				if(map[0]==OLD)	//if "p" is always better, we do not need to update anything)
				{
					if(previous)		//close previous interval and put in into the queue
					{
						previous->next() = p;
						previous->compute_min_distance(p->start());
						m_queue.insert(previous);
						++winNum;
						previous = NULL;
					}

					p = p->next(); 
					
				}
				else if(previous)	//extend previous interval to cover everything; remove p
				{
					previous->next() = p->next(); 
					erase_from_queue(p);
					m_memory_allocator.deallocate(p);

					p = previous->next();
				}
				else				//p becomes "previous"
				{
					previous = p;
					interval_pointer next = p->next();
					erase_from_queue(p);

					memcpy(previous,q,sizeof(Interval));

					previous->start() = start[0];
					previous->next() = next;

					p = next; 
				}
				continue;
			}

			//update_flag = true;

			Interval swap(*p);							//used for swapping information
			propagate_flag = erase_from_queue(p);

			for(unsigned j=1; j<N; ++j)				//no memory is needed for the first one
			{
				i_new[j] = m_memory_allocator.allocate();	//create new intervals
			}

			if(map[0]==OLD)	//finish previous, if any
			{
				if(previous)
				{
					previous->next() = p;
					previous->compute_min_distance(previous->stop());
					m_queue.insert(previous);
					++winNum;
					previous = NULL;
				}
				i_new[0] = p;
				p->next() = i_new[1];
				p->start() = start[0];
			}
			else if(previous)	//extend previous interval to cover everything; remove p
			{
				i_new[0] = previous;
				previous->next() = i_new[1]; 
				m_memory_allocator.deallocate(p);
				previous = NULL;
			}
			else				//p becomes "previous"
			{
				i_new[0] = p;
				memcpy(p,q,sizeof(Interval));

				p->next() = i_new[1];
				p->start() = start[0];
			}

			assert(!previous);

			for(unsigned j=1; j<N; ++j)					
			{
				interval_pointer current_interval = i_new[j];

				if(map[j] == OLD)	
				{
					memcpy(current_interval,&swap,sizeof(Interval));
				}
				else
				{
					memcpy(current_interval,q,sizeof(Interval));
				}
				
				if(j == N-1)	
				{
					current_interval->next() = swap.next();
				}
				else			
				{
					current_interval->next() = i_new[j+1];
				}

				current_interval->start() = start[j];
			}

			for(unsigned j=0; j<N; ++j)								//find "min" and add the intervals to the queue
			{
				if(j==N-1 && map[j]==NEW)
				{
					previous = i_new[j];
				}
				else
				{
					interval_pointer current_interval = i_new[j];

					current_interval->compute_min_distance(current_interval->stop());					//compute minimal distance

					if(map[j]==NEW || (map[j]==OLD && propagate_flag))
					{
						m_queue.insert(current_interval);
						++winNum;
					}
				}
			}

			p = swap.next();
		}

		if(previous)		//close previous interval and put in into the queue
		{
			previous->compute_min_distance(previous->stop());
			m_queue.insert(previous);
			++winNum;
			previous = NULL;
		}
	}
}

}		//geodesic

#endif //GEODESIC_ALGORITHM_EXACT
