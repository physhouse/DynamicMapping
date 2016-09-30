#ifndef _FG_ATOMS_H
#define _FG_ATOMS_H

#include <fstream>
#include "pointers.h"

struct FrameSource {
  std::ifstream trajectory;
  char trajName[1024];
  int header_size;
  int x_pos;
  int v_pos;
  int type_pos;
  std::string* elements;
};

class Fg_atoms : protected Pointers
{
public:
   Fg_atoms(class Mapping* map);
   void init(int argc, char** argv);
   void readInitialFrame();
   void readNextFrame();
   void readFinalFrame();
  
   //helper functions
  
   void readLammpsHeader();
   void readLammpsBody();
   int  stringSplit(std::string, const char*, std::string*);
   void finishReading();
   void cleanup();
   void output();

   double	L;
   int		fg_num;
   int		nframes; //Number of frames in the trajectory
   int		nSkipFrames;
   double**	r;
   double*	v;
   int		currentStep;
   FrameSource* input;
   std::ofstream fgtrj;
};

#endif
