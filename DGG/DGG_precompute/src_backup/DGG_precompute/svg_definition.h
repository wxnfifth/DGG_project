#ifndef _SVG_DEFINITION_H_

struct HeadOfSVG {
    int begin_vertex_index;
    int end_vertex_index;
    int num_of_vertex;
    HeadOfSVG(int _begin_vertex_index , 
        int _end_vertex_index , 
        int _num_of_vertex)
        :
    begin_vertex_index(_begin_vertex_index) ,
        end_vertex_index(_end_vertex_index) , 
        num_of_vertex(_num_of_vertex){}
    HeadOfSVG(){}
};

struct BodyHeadOfSVG{
    int source_index;
    int neighbor_num;
    BodyHeadOfSVG(int _source_index , int _neighbor_num)
    {
        source_index = _source_index;
        neighbor_num = _neighbor_num;
    }
    BodyHeadOfSVG(){}
};
struct BodyPartOfSVG{
    int dest_index;
    float dest_dis;
    BodyPartOfSVG(){}
    BodyPartOfSVG(int _dest_index , float _dest_dis ):
        dest_index(_dest_index),
        dest_dis(_dest_dis)
    {
        dest_index = _dest_index;
        dest_dis = _dest_dis;
    }
};

struct BodyPartOfSVGWithAngle{
    int dest_index;
    float dest_dis;
    float angle;
    int begin_pos; 
    int end_pos;
    BodyPartOfSVGWithAngle(){}
    BodyPartOfSVGWithAngle(int _dest_index , float _dest_dis, float _angle, int _begin_pos, int _end_pos):
        dest_index(_dest_index),
        dest_dis(_dest_dis),
        angle(_angle),
        begin_pos(_begin_pos),
        end_pos(_end_pos)
    {
    }
    bool operator<(const BodyPartOfSVGWithAngle& o) const{
      return angle < o.angle;
    }
};


struct BodyPartOfSVGWithK : BodyPartOfSVG{
    int rank_k;
    BodyPartOfSVGWithK(){}
    BodyPartOfSVGWithK(int _dest_index , float _dest_dis , int _rank_k ):
        BodyPartOfSVG(_dest_index , _dest_dis),
        rank_k(_rank_k){}
    bool operator<(const BodyPartOfSVGWithK& other)const{
        return rank_k < other.rank_k;
    }
};

#endif
