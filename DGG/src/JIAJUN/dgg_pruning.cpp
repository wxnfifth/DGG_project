// LocalGeodesics.cpp : Defines the entry point for the console application.
//
#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<algorithm>
#include<queue>
#include<vector>
#include <windows.h>
#include<map>
#include<assert.h>
#include<iostream>
#include<string>
#include<sstream>
#include "wxnTime.h"
#include "svg_definition.h"
using namespace std;

namespace JIAJUN_DGG_PRUNING{

#define _printToFile 0
#define _debug 0
  const double maxDist = 1e5;
  //const double eps = 1e-3;

  class SVGEdge{
  public:
    int v;
    bool deleted;
    double dis;
    short int begin_pos;
    short int end_pos;
    short int pos;
    SVGEdge(): deleted(false){}

  };


  class SVG{
  public:
    int N; 
    int K;
    int * degree;
    SVGEdge * edge;
    SVG() : degree(NULL), edge(NULL){	}
    ~SVG(){ if(degree != NULL) delete [] degree; if(edge != NULL) delete [] edge; }
    SVG(const SVG & svg){
      N = svg.N;
      K = svg.K;
      degree = new int[sizeof(*svg.degree) / sizeof(int)];
      for (int i = 0; i != sizeof(*svg.degree) / sizeof(int); ++i)
        degree[i] = svg.degree[i];
      edge = new SVGEdge[sizeof(*svg.edge) / sizeof(SVGEdge)];
      for (int i = 0; i != sizeof(*svg.edge) / sizeof(SVGEdge); ++i)
        edge[i] = svg.edge[i];
    }
    SVG& operator=(const SVG & svg)
    {
      N = svg.N;
      K = svg.K;
      delete [] degree;
      delete [] edge;
      degree = new int[sizeof(*svg.degree) / sizeof(int)];
      for (int i = 0; i != sizeof(*svg.degree) / sizeof(int); ++i)
        degree[i] = svg.degree[i];
      edge = new SVGEdge[sizeof(*svg.edge) / sizeof(SVGEdge)];
      for (int i = 0; i != sizeof(*svg.edge) / sizeof(SVGEdge); ++i)
        edge[i] = svg.edge[i];
    }
    virtual void dijkstra(int src, double * dis, bool * mark, double) = 0;
    //virtual void read(const char *);

  };

  class DGG : public SVG{
  public:
    std::vector<int> origin_neigh_num;
    double eps;
    //void read(const char *);
    void readWxnBinary(const char *);
    void readWxnBinary_new(const char *);
    void dijkstra(int src, double * dis, bool * mark, double);
    void pruning();
    void removeDeletedEdges();
    void write(const char *);
    void writeSVGBinary(const char *);
    void writeSVGBinary_new(const char *);
  };

  void DGG::readWxnBinary(const char* svg_filename)
  {
    {
      std::ifstream input_file (svg_filename, std::ios::in | std::ios::binary);
      HeadOfSVG head_of_svg;
      input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
      N = head_of_svg.num_of_vertex;
      origin_neigh_num.resize(N);
      degree = new int[N+1];
      degree[0] = 0;
      for (int i = 0; i < head_of_svg.num_of_vertex; ++i) {
        BodyHeadOfSVG body_head;
        input_file.read( (char*)&body_head , sizeof(body_head));
        origin_neigh_num[i] = body_head.neighbor_num;
        degree[i+1] = degree[i] + body_head.neighbor_num;
        for(int j = 0; j < body_head.neighbor_num;++j){ 
          BodyPartOfSVGWithAngle body_part;
          input_file.read((char*)&body_part , sizeof(body_part));
        }
      }
      input_file.close();
    }
    {
      std::ifstream input_file(svg_filename, std::ios::in | std::ios::binary);
      edge = new SVGEdge[degree[N]+100];
      HeadOfSVG head_of_svg;
      input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
      //printf("%d", degree[N]);
      for(int i = 0; i < N; ++i) {
        BodyHeadOfSVG body_head;
        input_file.read( (char*)&body_head , sizeof(body_head));
        for(int j = degree[i]; j < degree[i+1]; ++j){
          BodyPartOfSVGWithAngle body_part;
          input_file.read((char*)&body_part , sizeof(body_part));
          //printf("%d\t", __temp__tempSVGEdge__.v);
          edge[j].v = body_part.dest_index;
          edge[j].dis = body_part.dest_dis;
          edge[j].deleted = 0;
          edge[j].begin_pos = body_part.begin_pos;
          edge[j].end_pos = body_part.end_pos;
          edge[j].pos = j - degree[i];
        }
      }
    }


  }
  void DGG::readWxnBinary_new(const char* svg_filename)
  {
    {
      std::ifstream input_file (svg_filename, std::ios::in | std::ios::binary);
      HeadOfSVG head_of_svg;
      input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
      N = head_of_svg.num_of_vertex;
      origin_neigh_num.resize(N);
      degree = new int[N+1];
      degree[0] = 0;
      for (int i = 0; i < head_of_svg.num_of_vertex; ++i) {
        BodyHeadOfSVG body_head;
        input_file.read( (char*)&body_head , sizeof(body_head));
        origin_neigh_num[i] = body_head.neighbor_num;
        degree[i+1] = degree[i] + body_head.neighbor_num;
        for(int j = 0; j < body_head.neighbor_num;++j){ 
          BodyPartOfSVGWithRange body_part;
          input_file.read((char*)&body_part , sizeof(body_part));
        }
      }
      input_file.close();
    }
    {
      std::ifstream input_file(svg_filename, std::ios::in | std::ios::binary);
      edge = new SVGEdge[degree[N]+100];
      HeadOfSVG head_of_svg;
      input_file.read( (char*)&head_of_svg , sizeof(head_of_svg));
      //printf("%d", degree[N]);
      for(int i = 0; i < N; ++i) {
        BodyHeadOfSVG body_head;
        input_file.read( (char*)&body_head , sizeof(body_head));
        for(int j = degree[i]; j < degree[i+1]; ++j){
          BodyPartOfSVGWithRange body_part;
          input_file.read((char*)&body_part , sizeof(body_part));
          //printf("%d\t", __temp__tempSVGEdge__.v);
          edge[j].v = body_part.dest_index;
          edge[j].dis = body_part.dest_dis;
          edge[j].deleted = 0;
          edge[j].begin_pos = body_part.begin_pos;
          edge[j].end_pos = body_part.end_pos;
          edge[j].pos = j - degree[i];
        }
      }
    }


  }

#if _debug
  void SVG::dijkstra(int src, double * dis, bool * mark)
  {
  }
#endif
  void DGG::dijkstra(int src, double * dis, bool * mark, double eps){
    //printf("%d", eps);
    DGG & dgg = * this;
    typedef SVGEdge node;
    const double maxError = 1e-5;
    
    struct CMP{
      node a, b;
      bool operator()(node a, node b){ return a.dis > b.dis; }
    };
    std::priority_queue <node, std::vector<node>, CMP > q;

    dis[src] = 0;
    mark[src] = 1;
    std::map<int, int> nodeMap;
    for (int i = dgg.degree[src]; i != dgg.degree[src + 1]; ++i)
    {
      dis[dgg.edge[i].v] = dgg.edge[i].dis;
      q.push(dgg.edge[i]);
      nodeMap[dgg.edge[i].v] = i;
    }
    int cnt = 0;
    while (!q.empty()) {
      node a = q.top();
      cnt++;
      q.pop();
      if(mark[a.v])continue;
      bool found = 0;
      for(int i = dgg.degree[a.v]; i < dgg.degree[a.v + 1]; ++i){
        if(abs(dis[dgg.edge[i].v] - maxDist) < maxError || a.v == dgg.edge[i].v) continue;
        if(a.dis + dgg.edge[i].dis < dis[dgg.edge[i].v] * (1 + eps)){
          node b;
          b.v = dgg.edge[i].v;
          
          b.dis = min((double)a.dis + (double)dgg.edge[i].dis, dis[dgg.edge[i].v] );
          dis[b.v] = b.dis;
          assert(nodeMap.count(dgg.edge[i].v));
          dgg.edge[nodeMap[dgg.edge[i].v]].deleted = 1;
#if _debug
          fprintf(stderr,"Found:: src: %d a.v: %d dgg.edge[i]: %d              %lf   %lf   %lf\n",src, a.v, dgg.edge[i].v, dis[dgg.edge[i].v], a.dis, dgg.edge[i].dis);
          system("pause");
#endif
#if _debug
          if (!found){

            fprintf(stderr,"Not found!! %d %d %d %lf %lf %lf \n", src, a.v, dgg.edge[i].v, dis[dgg.edge[i].v], a.dis, dgg.edge[i].dis);
            fprintf(stderr,"%d\n",src);

            for (int j = dgg.degree[dgg.degree[src] ]; j < dgg.degree[dgg.degree[src] + 1]; ++j)
              fprintf(stderr,"%d  ", dgg.edge[j].v);
            fprintf(stderr,"\n");

            system("pause");
#endif      //++++++dgg.edge[i].deleted = 1;
#if _debug				
            fprintf(stderr,"%d %d %d %lf   %lf   %lf %d\n", src, a.v, i, dis[dgg.edge[i].v], a.dis, dgg.edge[i].dis, finalDeg);
            if (i % 5 == 0) system("pause");
#endif
            q.push(b);
          }
        }
        mark[a.v] = 1;
      }
      //fprintf(stderr,"cnt %d\n" , cnt);
      for (int i = dgg.degree[src]; i != dgg.degree[src + 1]; ++i)
      {
        dis[dgg.edge[i].v] = maxDist;
        mark[dgg.edge[i].v] = 0;
      }
      dis[src] = maxDist;
      mark[src] = 0;
#if _debug
      fprintf(stderr,"%d  %d\n", dgg.degree[src + 1] - dgg.degree[src], finalDeg);
#endif
    }


    void DGG::writeSVGBinary(const char * output_filename)
    {
      ofstream output_file (output_filename , ios::out | ios::binary);
      HeadOfSVG head_of_svg(0 , N-1 , N );     
      output_file.write((char*)&head_of_svg , sizeof(head_of_svg));


			double ave_degree = 0;
      for (int source_index = 0; source_index < N; ++source_index) {
        ave_degree += degree[source_index+1] - degree[source_index];
				BodyHeadOfSVG body_header(source_index , degree[source_index+1] - degree[source_index]);
        output_file.write((char*)&body_header , sizeof(body_header));
        //printf("%d %d \n" , source_index,  degree[source_index+1] - degree[source_index]);
        //vector<BodyPartOfSVGWithK> body_parts(degree[source_index+1] - degree[source_index]);

        vector<int> origin2current(origin_neigh_num[source_index],-1);
        for (int j = degree[source_index]; j < degree[source_index+1];j++) {
          origin2current[edge[j].pos] = j - degree[source_index];
        }

        int start_pos = 0;
        for (int j = 0; j < origin2current.size(); ++j) {
          if (origin2current[j] == -1) {
            origin2current[j] = start_pos;
          }else{
            start_pos = origin2current[j];
          }
        }

        for (int j = degree[source_index]; j < degree[source_index+1]; j++) {
          BodyPartOfSVGWithAngle b;
          b.angle = 0;
          if (edge[j].begin_pos == -1) {
            b.begin_pos = -1;
            b.end_pos = -1;
          }else{
            b.begin_pos = origin2current[edge[j].begin_pos];
            b.end_pos = origin2current[edge[j].end_pos];
          }
          b.dest_dis = edge[j].dis;
          b.dest_index = edge[j].v;
          //printf("angle %lf pos %d dis %lf index %d end_pos %d\n" , b.angle, b.begin_pos, b.dest_dis, b.dest_index, b.end_pos);
          output_file.write((char*)&b , sizeof(b));
        }
      }
      output_file.close();
			printf("ave degree after pruning %lf\n" , ave_degree);
		}

    void DGG::writeSVGBinary_new(const char * output_filename)
    {
      ofstream output_file (output_filename , ios::out | ios::binary);
      HeadOfSVG head_of_svg(0 , N-1 , N );     
      output_file.write((char*)&head_of_svg , sizeof(head_of_svg));

      for (int source_index = 0; source_index < N; ++source_index) {
        BodyHeadOfSVG body_header(source_index , degree[source_index+1] - degree[source_index]);
        output_file.write((char*)&body_header , sizeof(body_header));

        vector<int> origin2current(origin_neigh_num[source_index],-1);
        for (int j = degree[source_index]; j < degree[source_index+1];j++) {
          origin2current[edge[j].pos] = j - degree[source_index];
        }

        int start_pos = 0;
        for (int j = 0; j < origin2current.size(); ++j) {
          if (origin2current[j] == -1) {
            origin2current[j] = start_pos;
          }else{
            start_pos = origin2current[j];
          }
        }

        for (int j = degree[source_index]; j < degree[source_index+1]; j++) {
          BodyPartOfSVGWithRange b;
          if (edge[j].begin_pos == -1) {
            b.begin_pos = -1;
            b.end_pos = -1;
          }else{
            b.begin_pos = origin2current[edge[j].begin_pos];
            b.end_pos = origin2current[edge[j].end_pos];
          }
          b.dest_dis = edge[j].dis;
          b.dest_index = edge[j].v;
          //printf("angle %lf pos %d dis %lf index %d end_pos %d\n" , b.angle, b.begin_pos, b.dest_dis, b.dest_index, b.end_pos);
          output_file.write((char*)&b , sizeof(b));
        }
      }
      output_file.close();
    }

    void DGG::write(const char * filename)
    {
      //printf("__filename %s\n" , filename);
      //assert(strlen(filename) >= 70);
      fprintf(stderr,"Writing start");
      char buf[1024];
      sprintf(buf, "%s_new.edge", filename);
      fprintf(stderr,"%s", buf);
      FILE * fp = fopen(buf, "wb");
      fwrite(edge, sizeof(SVGEdge), degree[N], fp);
      fclose(fp);
      fprintf(stderr,"Writing edge done!\n");
      sprintf(buf, "%s_new_degree.txt", filename);
      fp = fopen(buf, "w");	
      fprintf(fp, "%d\n" , N);
      for(int i = 0; i < N; ++i){
        fprintf(fp, "%d\n", degree[i+1] - degree[i]);
      }
      fclose(fp);
      fprintf(stderr,"Written to file!\n");

    }

#if _debug
    void dgg_dijkstra_2(const SVG & dgg, int src, double * dis, bool * mark){

      int iniDeg = 0, finalDeg = 0;
      typedef SVGEdge node;
      node a, b;
      std::map <int, int> nodeToLabel;
      for (int i = dgg.degree[src]; i != dgg.degree[src + 1]; ++i)
      {
        nodeToLabel[dgg.edge[i].v] = i;
      }
      struct CMP{
        bool operator()(node a, node b){ return a.dis > b.dis; }
      };

      std::priority_queue <node, std::vector<node>, CMP> q;
      int n = dgg.degree[src + 1] - dgg.degree[src];
      dis[src] = 0;
      mark[src] = 1;

      for (int i = dgg.degree[src]; i != dgg.degree[src + 1]; ++i)
      {
        a.v = dgg.edge[i].v;
        a.dis = dgg.edge[i].dis;
        dis[a.v] = a.dis;
#ifdef _debug
        fprintf(stderr,"%d %lf\n", a.v, a.dis);
#endif
        q.push(a);
      }
      while (!q.empty()){
        a = q.top();
        q.pop();
        if(mark[a.v])continue;
        bool found = 0;
        for(int i = dgg.degree[a.v]; i < dgg.degree[a.v + 1]; ++i){
          if(abs(dis[dgg.edge[i].v] - maxDist) < 1e-5 || a.v == dgg.edge[i].v) continue;
#ifdef _debug
          fprintf(stderr,"src: %d a.v: %d dgg.edge[i]: %d              %lf   %lf   %lf\n",src, a.v, dgg.edge[i].v, dis[dgg.edge[i].v], a.dis, dgg.edge[i].dis);
#endif
          if(a.dis + dgg.edge[i].dis < dis[dgg.edge[i].v] * (1 + eps)){
            b.v = dgg.edge[i].v;
            b.dis = min(a.dis + dgg.edge[i].dis, dis[dgg.edge[i].v] );
            dis[b.v] = b.dis;
            found = 0;

            int j = nodeToLabel[b.v];
            if (!dgg.edge[j].deleted)finalDeg ++;
            dgg.edge[j].deleted = 1;

#ifdef _debug				
            fprintf(stderr,"%d %d %d %lf   %lf   %lf %d\n", src, a.v, i, dis[dgg.edge[i].v], a.dis, dgg.edge[i].dis, finalDeg);
            if (i % 5 == 0) system("pause");
#endif
            q.push(b);
          }
        }
        mark[a.v] = 1;
      }

      for (int i = dgg.degree[src]; i != dgg.degree[src + 1]; ++i)
      {
        dis[dgg.edge[i].v] = maxDist;
        mark[dgg.edge[i].v] = 0;
      }
      dis[src] = maxDist;
      mark[src] = 0;
#ifdef _debug
      fprintf(stderr,"%d  %d\n", dgg.degree[src + 1] - dgg.degree[src], finalDeg);
#endif
    }
#endif

    void DGG::pruning(){
      //read file
      fprintf(stderr,"%lf\n", eps);
      DGG & dgg = * this;
      int Count = 0;
      //pruning using dijkstra
      //initialization
      double * dis = new double[dgg.degree[dgg.N] + 100];
      bool * mark = new bool[dgg.degree[dgg.N] + 100];
      for (int i = 0; i != dgg.degree[dgg.N] + 100; ++i)
      {
        dis[i] = maxDist;
        mark[i] = 0;
      }

      for (int i = 0; i != dgg.N; ++i)
      {
        dijkstra(i, dis, mark, eps);
#if _debug
        for (int j = 0; j != dgg.degree[dgg.N]; ++ j)
          if (abs(dis[j] - maxDist) > eps || mark[j])
          {
            fprintf(stderr,"Error!! %d", j );
            system("pause");
          }
#endif
      }
      removeDeletedEdges();
      delete [] dis;
      delete [] mark;
    }

    void DGG::removeDeletedEdges()
    {
      DGG & dgg = * this;
      //write new dgg
      DGG newdgg;
      newdgg.N = dgg.N;
      newdgg.K = dgg.K;
      newdgg.degree = new int[newdgg.N+1];
      newdgg.degree[0] = 0;
      int numOfEdges = 0;
      for (int i = 0; i < dgg.N; ++i){
        numOfEdges += std::count_if(dgg.edge + dgg.degree[i], dgg.edge + dgg.degree[i + 1], 
          [](SVGEdge e){return !e.deleted;});
        newdgg.degree[i + 1] = numOfEdges;
      }
      fprintf(stderr,"Total number of edges: %d\n", numOfEdges);
      newdgg.edge = new SVGEdge[numOfEdges + 100];

      int _numOfEdges = 0;
      for (int i=0; i < dgg.N; ++i){
        for (int j = dgg.degree[i]; j < dgg.degree[i+1]; ++j){
          if (dgg.edge[j].deleted == 0)
          {
            newdgg.edge[_numOfEdges] = dgg.edge[j];
            _numOfEdges++;
          }
        }
      }
      //printf("Total number of edges: %d\n", _numOfEdges);
      assert(numOfEdges == _numOfEdges);
      delete [] dgg.degree;
      delete [] dgg.edge;
      dgg.N = newdgg.N;
      dgg.K = newdgg.K;
      dgg.edge = newdgg.edge;
      newdgg.edge = NULL;
      dgg.degree = newdgg.degree;
      newdgg.degree = NULL;
      fprintf(stderr,"Removing edge#edge: %d\n", dgg.degree[dgg.N]);
    }


    void dgg_pruning(const std::string& input_file_name, double eps, std::string& output_filename, double& prune_time)
    {
      DGG dgg;
      dgg.readWxnBinary(input_file_name.c_str());
      dgg.eps = eps;
			double average_degree_before = dgg.degree[dgg.N] / dgg.N;
			printf("Average_degree_before %lf\n", average_degree_before);
      ElapasedTime t;
      dgg.pruning();
      prune_time = t.getTime();
			double average_degree_after = dgg.degree[dgg.N] / dgg.N;
			printf("Average_degree_after %lf , percenter %lf\n", average_degree_after, (double)average_degree_after / average_degree_before );
      output_filename = input_file_name.substr(0,input_file_name.length()-7) + "_pruning.binary";
      dgg.writeSVGBinary(output_filename.c_str());
    }
    void dgg_pruning_new(const std::string& input_file_name, double eps, std::string& output_filename, double& prune_time)
    {
      DGG dgg;
      dgg.readWxnBinary_new(input_file_name.c_str());
      dgg.eps = eps;
      fprintf(stderr,"Before: %d\n", dgg.degree[dgg.N]);
      ElapasedTime t;
      dgg.pruning();
      prune_time = t.getTime();
      fprintf(stderr,"After: %d\n", dgg.degree[dgg.N]);
      output_filename = input_file_name.substr(0,input_file_name.length()-7) + "_pruning.binary";
      dgg.writeSVGBinary_new(output_filename.c_str());
    }

  }
