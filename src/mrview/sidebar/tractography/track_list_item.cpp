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

    15-12-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * a few bug fixes + memory performance improvements for the depth blend option
    
    12-03-2010 J-Donald Tournier <d.tournier@brain.org.au>
    * a few bug fixes for the colour handling and support for depth blend on
    * 64 bit systems
    
*/

#include <gtk/gtkgl.h>
#include <GL/gl.h>
#include <GL/glext.h>
#include <glibmm/stringutils.h>

#include "use_gl.h"
#include "mrview/sidebar/tractography/track_list_item.h"
#include "dwi/tractography/file.h"

namespace MR {
  namespace Viewer {
    namespace SideBar {

      Point TrackListItem::Track::Point::normal;
      const TrackListItem::Track::Point* TrackListItem::Track::Point::root = NULL;

      void TrackListItem::load (const String& filename) 
      {
        file = filename;
        tracks.clear();
        alloc.clear(); 
        struct_stat64 S;
        if (stat64 (file.c_str(), &S)) throw Exception ("error accessing tracks file \"" + file + "\": " + Glib::strerror (errno));
        mtime = S.st_mtime;

        DWI::Tractography::Reader reader;
        reader.open (file, properties); 
        GLubyte default_colour[] = { 255, 255, 255 };
        GLubyte A = (GLubyte) (255.0*get_alpha());

        std::vector<Point> tck;
        while (reader.next (tck)) {
          if (tck.size()) {
            Track dest (alloc, tck.size());

            if (colour_by_dir) {
              dest[0].set_pos (tck[0]);
              if (dest.size() > 1) {
                Point dir = tck[1] - tck[0];
                dest[0].set_colour (dir.normalise(), A);

                guint i;
                for (i = 1; i < dest.size()-1; i++) {
                  dest[i].set_pos (tck[i]);
                  dir = tck[i+1] - tck[i-1];
                  dest[i].set_colour (dir.normalise(), A);
                }

                dest[i].set_pos (tck[i]);
                dir = tck[i] - tck[i-1];
                dest[i].set_colour (dir.normalise(), A);
              }
              else dest[0].set_colour (default_colour, A);
            }
            else {
              for (guint i = 0; i < dest.size(); i++) {
                dest[i].set_pos (tck[i]);
                dest[i].set_colour (colour, A);
              }
            }

            tracks.push_back (dest);
          }
        }

      }








      void TrackListItem::draw () 
      {
        glBlendColor (1.0, 1.0, 1.0, get_alpha());
        for (std::list<Track>::iterator i = tracks.begin(); i != tracks.end(); ++i) {
          const Track& tck (*i);
          glVertexPointer (3, GL_FLOAT, sizeof (Track::Point), tck[0].pos);
          glColorPointer (3, GL_UNSIGNED_BYTE, sizeof (Track::Point), tck[0].C);
          glDrawArrays (GL_LINE_STRIP, 0, tck.size());
        }
      }


    }
  }
}

