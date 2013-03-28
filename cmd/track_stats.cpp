/*
Written by Richard Beare, Murdoch Childrens Research Institute.
Derived from filter_tracks
*/
#include <omp.h>

#include "app.h"
#include "image/interp.h"
#include "math/vector.h"
#include "point.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/roi.h"
#include "dwi/tractography/tracker/base.h"


using namespace std;
using namespace MR;
using namespace MR::DWI;
using namespace MR::DWI::Tractography;

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "Produce a set of track statistics (length,..) \n ",
  NULL
};

ARGUMENTS = {
  Argument ("input",  "input tracks file",  "the input file containing the tracks to be filtered.").type_file(),
  Argument ("output", "output stats file", "the output file containing the tracks selected.")     .type_file(),
  Argument::End
};

OPTIONS = {  
  Option::End
};






class Stats_filter {

  public:
  //Stats_filter (){}
  float measure(const std::vector<Point>& tck)
  {
    // just measure length
    float total=0.0;
    for (unsigned k=1;k<tck.size();k++)
      {
	Point Pd=tck[k]-tck[k-1];
	total+=Pd.norm();
      }
    return(total);
  }


};




EXECUTE
{

  Reader reader;
  Properties properties;
  Writer writer;

  reader.open (argument[0].get_string(), properties);
  const float progress_multiplier = properties["count"].empty() ? 0.0 : 100.0 / to<float> (properties["count"]);



  Stats_filter filter;
  //writer.create (argument[1].get_string(), properties);
  std::ofstream out(argument[1].get_string());
  if (!out) throw Exception ("cannot open statistics file \"" + std::string(argument[1].get_string()) + "\": " + Glib::strerror (errno));

  out << "Index,Length" << std::endl;
  std::vector<Point> tck;
  unsigned count=0;
  while (reader.next (tck)) 
    {
    float len = filter.measure (tck);
    out << count << "," << len << std::endl;
    ++count;
    fprintf (stderr, "\r%8u read,  [%3d%%]",
         count, int(progress_multiplier * count));
    }
  reader.close();

}
