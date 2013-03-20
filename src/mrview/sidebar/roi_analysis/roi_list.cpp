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

*/

#include <gtkmm/treerowreference.h>
#include <gtkmm/colorselection.h>
#include <gtkmm/stock.h>
#include <gdkmm/pixbuf.h>

#include "mrview/sidebar/roi_analysis/roi_list.h"
#include "mrview/sidebar/roi_analysis.h"
#include "mrview/window.h"
#include "dialog/file.h"

namespace MR {
  namespace Viewer {
    namespace SideBar {

      DP_ROIList::DP_ROIList (const ROIAnalysis& sidebar) : parent (sidebar), set (true), editing (false)
      {
        using namespace Gtk::Menu_Helpers;
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::OPEN, sigc::mem_fun(*this, &DP_ROIList::on_open) ) );
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::NEW, sigc::mem_fun(*this, &DP_ROIList::on_new) ) );
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLOSE, sigc::mem_fun(*this, &DP_ROIList::on_close) ) );
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (MenuElem("_Change colour", sigc::mem_fun(*this, &DP_ROIList::on_set_colour)));
        popup_menu.items().push_back (SeparatorElem());
        popup_menu.items().push_back (StockMenuElem(Gtk::Stock::CLEAR, sigc::mem_fun(*this, &DP_ROIList::on_clear) ) );
        popup_menu.accelerate (*this);

        model = Gtk::ListStore::create (columns);
        set_model (model);

        int tick_column_index = append_column_editable ("", columns.show) - 1;
        append_column ("", columns.pix);
        append_column ("file", columns.name);

        set_tooltip_text ("right-click for more options");

        set_headers_visible (false);

        Gtk::CellRendererToggle* tick = dynamic_cast<Gtk::CellRendererToggle*> (get_column_cell_renderer (tick_column_index));
        tick->signal_toggled().connect (sigc::mem_fun (*this, &DP_ROIList::on_tick));
      }






      DP_ROIList::~DP_ROIList()  { }





      void DP_ROIList::draw (int transparency)
      {
        Gtk::TreeModel::Children rois = model->children();
        if (rois.size() == 0) return;

        Pane& pane (Window::Main->pane());
        Slice::Current slice (pane);

        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glEnable (GL_TEXTURE_2D);
        glDisable (GL_DEPTH_TEST);
        glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
        glDepthMask (GL_FALSE);
        glColorMask (GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);

        for (Gtk::TreeModel::Children::iterator iter = rois.begin(); iter != rois.end(); ++iter) {
          bool show = (*iter)[columns.show];
          if (show) {
            RefPtr<ROI> roi = (*iter)[columns.roi];

            Slice::Info S;
            S.image = roi->mask;
            S.focus = slice.focus;

            bool transforms_differ = false;
            if (slice.orientation) 
              transforms_differ = true;
            else {
              for (guint i = 0; i < 4; ++i) {
                for (guint j = 0; j < 3; ++j) {
                  if (fabs (slice.image->image->header().transform()(i,j) - roi->mask->image->header().transform()(i,j)) > 1e-5) {
                    transforms_differ = true;
                    break;
                  }
                }
              }
            }

            if (transforms_differ) {
              const GLdouble* M (pane.get_modelview());
              float matrix[] = { 
                float(M[0]), float(M[1]), float(M[2]),
                float(M[4]), float(M[5]), float(M[6]),
                float(M[8]), float(M[9]), float(M[10])
              };
              S.orientation.from_matrix (matrix);
              S.projection = 2;
            }
            else 
              S.projection = slice.projection;

            S.interpolate = false;

            Slice::Current current (S);

            roi->render.update (current);
            glColor4ub ((roi->colour & 0xFF000000) >> 24, (roi->colour & 0x00FF0000) >> 16, (roi->colour & 0x0000FF00) >> 8, transparency);
            roi->render.draw();
          }
        }
        glDepthMask (GL_TRUE);
        glDisable (GL_TEXTURE_2D);
      }







      void DP_ROIList::load (RefPtr<MR::Image::Object> image)
      { 
        guint32 colour = 0xFFFF00FF;
        Glib::RefPtr<Gdk::Pixbuf> pix = Gdk::Pixbuf::create  (Gdk::COLORSPACE_RGB, false, 8, 16, 16);
        pix->fill (colour);
        Gtk::TreeModel::Row row = *(model->append());
        row[columns.show] = true;
        row[columns.pix] = pix;
        row[columns.name] = image->name();
        row[columns.roi] = RefPtr<ROI> (new ROI (image, colour));;
        get_selection()->select (row);
      }







      bool DP_ROIList::on_button_press_event (GdkEventButton* event)
      {
        if ((event->type == GDK_BUTTON_PRESS) && (event->button == 3)) {
          Gtk::TreeModel::Path path;
          Gtk::TreeViewColumn* col;
          int x, y;
          bool is_row = get_path_at_pos ((int) event->x, (int) event->y, path, col, x, y);
          if (is_row) {
            if (!get_selection()->is_selected (path)) 
              TreeView::on_button_press_event(event);
          }
          else get_selection()->unselect_all();

          popup_menu.items()[2].set_sensitive (is_row);
          popup_menu.items()[4].set_sensitive (is_row);
          popup_menu.popup (event->button, event->time);
          return (true);
        }
        return (TreeView::on_button_press_event(event));
      }




      void DP_ROIList::on_set_colour () 
      {
        Gtk::ColorSelectionDialog dialog ("Choose colour for ROI");
        if (dialog.run() == Gtk::RESPONSE_OK) {
          Gtk::TreeModel::iterator iter = get_selection()->get_selected();
          RefPtr<ROI> roi = (*iter)[columns.roi];
          Gdk::Color colour (dialog.get_colorsel()->get_current_color());
          GLubyte C[] = { 
            GLubyte (colour.get_red() >> 8),
            GLubyte (colour.get_green() >> 8),
            GLubyte (colour.get_blue() >> 8)
          };
          roi->colour = C[0] << 24 | C[1] << 16 | C[2] << 8 | 255;
          Glib::RefPtr<Gdk::Pixbuf> pix = (*iter)[columns.pix];
          pix->fill (roi->colour);
          Window::Main->update (&parent);
        }
      }









      void DP_ROIList::on_open () 
      {
        Dialog::File dialog ("Open mask image", true, true);

        if (dialog.run() == Gtk::RESPONSE_OK) {
          std::vector<RefPtr<MR::Image::Object> > selection = dialog.get_images();
          if (selection.size()) {
            for (guint n = 0; n < selection.size(); n++) 
              load (selection[n]);
            Window::Main->update (&parent);
          }
        }
      }







      void DP_ROIList::on_new () 
      {
        Dialog::File dialog ("Create mask image", false, false);

        if (dialog.run() == Gtk::RESPONSE_OK) {
          std::vector<String> selection = dialog.get_selection();
          if (selection.size()) {
            MR::Image::Header header (*Window::Main->image->image);
            header.data_type = DataType::Bit;
            header.axes.set_ndim (3);
            RefPtr<MR::Image::Object> obj (new MR::Image::Object);
            try {
              obj->create (selection[0], header); 
              load (obj);
            }
            catch (...) { error ("error creating mask image \"" + selection[0] + "\": " + Glib::strerror (errno)); }
          }
        }
      }








      void DP_ROIList::on_close ()
      {
        std::list<Gtk::TreeModel::Path> paths = get_selection()->get_selected_rows();
        std::list<Gtk::TreeModel::RowReference> rows;

        for (std::list<Gtk::TreeModel::Path>::iterator row = paths.begin(); row != paths.end(); row++) 
          rows.push_back (Gtk::TreeModel::RowReference (model, *row));

        for (std::list<Gtk::TreeModel::RowReference>::iterator row = rows.begin(); row != rows.end(); row++) 
          model->erase (model->get_iter (row->get_path()));

        Window::Main->update (&parent);
      }





      void DP_ROIList::on_clear () 
      {
        model->clear(); 
        Window::Main->update (&parent);
      }




      void DP_ROIList::on_tick (const String& path) { Window::Main->update (&parent); }
      
      bool DP_ROIList::on_button_press (GdkEventButton* event, float brush) 
      {
        Gtk::TreeModel::iterator iter = get_selection()->get_selected();
        if (!iter) return (false);
	guint state = (event->state & MODIFIERS) & (~GDK_BUTTON1_MASK);
        if (state != GDK_SHIFT_MASK && 
            state != ( GDK_SHIFT_MASK | CTRL_CMD_MASK ))
          return (false);
        row = *iter;
        bool show = row[columns.show];
        if (!show) return (false);

        set = state == GDK_SHIFT_MASK;
        editing = true;

        process (event->x, event->y, brush);
        return (true);
      }







      void DP_ROIList::process (gdouble x, gdouble y, float brush)
      {
        RefPtr<ROI> roi = row[columns.roi];
        Point pos (roi->mask->interp->R2P (position (x, y)));
        MR::Image::Position ima (*roi->mask->image);
        int p[] = { round (pos[0]), round(pos[1]), round(pos[2]) };
        int e = ceil(brush/2.0);
        float dist = (brush*brush)/4.0;

        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);
#if 0
	// this choice should be selectable via the dialog box
        for (ima.set (2, p[2]-e); ima[2] <= p[2]+e; ima.inc(2)) {
          if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
          for (ima.set (1, p[1]-e); ima[1] <= p[1]+e; ima.inc(1)) {
            if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
            for (ima.set (0, p[0]-e); ima[0] <= p[0]+e; ima.inc(0)) {
              if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
              if ((ima[0]-p[0])*(ima[0]-p[0]) + (ima[1]-p[1])*(ima[1]-p[1]) + (ima[2]-p[2])*(ima[2]-p[2]) < dist)
                ima.value (set ? 1.0 : 0.0);
            }
          }
        }
#else
	const unsigned projection(S.projection);
	switch (projection) {
	case 0:
	  // sagittal
	  ima.set (0, p[0]);
	  for (ima.set (2, p[2]-e); ima[2] <= p[2]+e; ima.inc(2)) 
	    {
	    if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
	    for (ima.set (1, p[1]-e); ima[1] <= p[1]+e; ima.inc(1)) 
	      {
	      if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
	      if ((ima[1]-p[1])*(ima[1]-p[1]) + (ima[2]-p[2])*(ima[2]-p[2]) < dist)
		ima.value (set ? 1.0 : 0.0);
	      }
	    }
	  break;
	case 1:
	  // coronal
	  ima.set (1, p[1]);
	  for (ima.set (2, p[2]-e); ima[2] <= p[2]+e; ima.inc(2)) 
	    {
	    if (ima[2] < 0 || ima[2] >= ima.dim(2)) continue;
	    for (ima.set (0, p[0]-e); ima[0] <= p[0]+e; ima.inc(0)) 
	      {
	      if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
	      if ((ima[0]-p[0])*(ima[0]-p[0]) + (ima[2]-p[2])*(ima[2]-p[2]) < dist)
                ima.value (set ? 1.0 : 0.0);
	      }
	    }
	  break;
	case 2:
	  // axial
	  ima.set (2, p[2]);
	  for (ima.set (1, p[1]-e); ima[1] <= p[1]+e; ima.inc(1)) 
	    {
	    if (ima[1] < 0 || ima[1] >= ima.dim(1)) continue;
	    for (ima.set (0, p[0]-e); ima[0] <= p[0]+e; ima.inc(0)) 
	      {
	      if (ima[0] < 0 || ima[0] >= ima.dim(0)) continue;
	      if ((ima[0]-p[0])*(ima[0]-p[0]) + (ima[1]-p[1])*(ima[1]-p[1]) < dist)
		ima.value (set ? 1.0 : 0.0);
	      }
	    }
	  break;
	default:
	  std::cerr << "Unrecognised projection " << projection << " RJBs fault" << std::endl;
	}
#endif
        Window::Main->update (&parent);
      }



      Point DP_ROIList::position (gdouble x, gdouble y)
      {
        Pane& pane (Window::Main->pane());
        const Slice::Current S (pane);
        Point f = pane.model_to_screen (S.focus);
        f[0] = x; 
        f[1] = pane.height() - y;
        return (pane.screen_to_model (f));
      }

    }
  }
}


