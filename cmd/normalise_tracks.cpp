/*
    Copyright 2008 Brain Research Institute, Melbourne, Australia

    Written by J-Donald Tournier, 27/06/08.

    This file is part of MRtrix.

    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.


    31-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * various optimisations to improve performance

    03-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * skip points in the tracks file if they are outside the supplied warp

*/

#include <fstream>
#include <glibmm/stringutils.h>

#include "app.h"
#include "get_set.h"
#include "dwi/tractography/file.h"
#include "dwi/tractography/properties.h"

using namespace MR; 
using namespace MR::DWI; 
using namespace std; 

SET_VERSION_DEFAULT;

DESCRIPTION = {
  "apply a normalisation map to a tracks file.",
  NULL
};

ARGUMENTS = {
  Argument ("tracks", "track file", "the input track file.").type_file (),
  Argument ("transform", "transform image", "the image containing the transform.").type_image_in(),
  Argument ("output", "output file", "the output track file").type_file(),
  Argument::End
};



OPTIONS = { Option::End };




EXECUTE {
  Tractography::Properties properties;
  Tractography::Reader file;
  file.open (argument[0].get_string(), properties);

  Image::Object& tranform_image (*argument[1].get_image());

  Tractography::Writer writer;
  writer.create (argument[2].get_string(), properties);

  std::vector<Point> tck, out;

  Image::Interp interp (tranform_image);
  ProgressBar::init (0, "normalising tracks...");

  while (file.next (tck)) {
    out.clear();
    for (std::vector<Point>::iterator i = tck.begin(); i != tck.end(); ++i) {
      interp.R (*i);
      if (!interp) continue;
      Point p;
      interp.set(3,0); p[0] = interp.value(); if (!gsl_finite (p[0])) continue;
      interp.set(3,1); p[1] = interp.value(); if (!gsl_finite (p[1])) continue;
      interp.set(3,2); p[2] = interp.value(); if (!gsl_finite (p[2])) continue;
      out.push_back (p);
    }
    writer.append (out);
    writer.total_count++;
    ProgressBar::inc();
  }

  ProgressBar::done();
}

