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


    17-10-2008 J-Donald Tournier <d.tournier@brain.org.au>
    * modify version information printed out by "-version" option
      to take account of new command version, copyright and author variables

*/

#include <glibmm/stringutils.h>

#include "app.h"
#include "file/config.h"


#define NUM_DEFAULT_OPTIONS 5


namespace MR {

  namespace ProgressBar {
    namespace {

      const gchar* busy[] = {
        ".    ",
        " .   ",
        "  .  ",
        "   . ",
        "    .",
        "   . ",
        "  .  ",
        " .   "
      };



      void init_func_cmdline () 
      { 
      }


      void display_func_cmdline ()
      {
        if (gsl_isnan (multiplier)) 
          fprintf (stderr, "\r%s: %s %s", Glib::get_application_name().c_str(), message.c_str(), busy[percent%8]);
        else fprintf (stderr, "\r%s: %s %3u%%", Glib::get_application_name().c_str(), message.c_str(), percent);
      }


      void done_func_cmdline () 
      { 
        if (gsl_isnan (multiplier)) fprintf (stderr, "\r%s: %s  - ok\n", Glib::get_application_name().c_str(), message.c_str()); 
        else fprintf (stderr, "\r%s: %s %3u%%\n", Glib::get_application_name().c_str(), message.c_str(), 100); 
      }
    }
  }





  void cmdline_print (const String& msg) { std::cout << msg; }

  void cmdline_error (const String& msg) 
  {
    if (App::log_level) std::cerr << Glib::get_application_name() << ": " << msg << "\n"; 
  }

  void cmdline_info  (const String& msg) 
  { 
    if (App::log_level > 1) std::cerr << Glib::get_application_name() << " [INFO]: " <<  msg << "\n"; 
  }

  void cmdline_debug (const String& msg)
  { 
    if (App::log_level > 2) std::cerr << Glib::get_application_name() << " [DEBUG]: " <<  msg << "\n"; 
  }





  namespace {

    void print_formatted_paragraph (const String& header, const String& text, int header_indent, int indent, int width)
    {
      int current = fprintf (stderr, "%-*s%-*s ", header_indent, "", indent-header_indent-2, header.c_str());

      String::size_type start = 0, end;
      do {
        end = start;
        while (!g_ascii_isspace(text [end]) && end < text.size()) end++;
        String token (text.substr (start, end-start));
        if (current + (int) token.size() + 1 >= width) 
          current = fprintf (stderr, "\n%*s%s", indent, "", token.c_str()) - 1;
        else current += fprintf (stderr, " %s", token.c_str());
        start = end + 1;
      } while (end < text.size());
      fprintf (stderr, "\n");
    }

  }


#define HELP_INDENT 10
#define HELP_WIDTH  80


#define HELP_PURPOSE_INDENT 0, HELP_INDENT, HELP_WIDTH
#define HELP_ARG_INDENT 12, 24, HELP_WIDTH
#define HELP_OPTION_INDENT 2, 16, HELP_WIDTH



  int App::log_level = 1;
  const gchar**   App::command_description = NULL;
  const Argument* App::command_arguments = NULL;
  const Option*   App::command_options = NULL;
  const guint*    App::version = NULL;
  const gchar*    App::copyright = NULL;
  const gchar*    App::author = NULL;

  const Option App::default_options[] = {
    Option ("info", "display information", "display information messages."),
    Option ("quiet", "suppress reporting", "do not display information messages or progress status."),
    Option ("debug", "display debug messages", "display debugging messages."),
    Option ("help", "show help page", "display this information page and exit."),
    Option ("version", "show version", "display version information and exit.")
  };




  guint App::match_option (const gchar* stub) const
  {
    std::vector<guint> candidates;
    String s (stub);

    for (guint n = 0; command_options[n].is_valid(); n++) 
      if (s.compare (0, s.size(), command_options[n].sname, s.size()) == 0)
        candidates.push_back (n);

    for (guint n = 0; n < NUM_DEFAULT_OPTIONS; n++) 
      if (s.compare (0, s.size(), default_options[n].sname, s.size()) == 0)
        candidates.push_back (n + DEFAULT_OPTIONS_OFFSET);

    if (candidates.size() == 0) return (G_MAXUINT);
    if (candidates.size() == 1) return (candidates[0]);

    s = "several matches possible for option \"" + s + "\": \"" + option_name (candidates[0]) + "\", \"" + option_name (candidates[1]) + "\"";
    for (guint n = 2; n < candidates.size(); n++) { s += ", "; s += option_name (candidates[n]); s += "\""; }
    throw Exception (s);
  }




  App::App (int argc, gchar** argv, const gchar** cmd_desc, const MR::Argument* cmd_args, const MR::Option* cmd_opts, 
          const guint* cmd_version, const char* cmd_author, const char* cmd_copyright)
  {
#ifdef G_OS_WIN32
    // force stderr to be unbuffered, and stdout to be line-buffered:
    setvbuf (stderr, NULL, _IONBF, 0);
    setvbuf (stdout, NULL, _IOLBF, 0);
#endif

    command_description = cmd_desc;
    command_arguments = cmd_args;
    command_options = cmd_opts;
    author = cmd_author;
    version = cmd_version;
    copyright = cmd_copyright;

    if (argc == 2) {
      if (strcmp (argv[1], "__print_full_usage__") == 0) {
        print_full_usage ();
        throw (0);
      }
    }

    String application_name (Glib::path_get_basename(argv[0]));
#ifdef G_OS_WIN32
    if (Glib::str_has_suffix (application_name, ".exe")) 
      application_name.erase (application_name.size()-4);
#endif
    Glib::set_application_name (application_name);

    log_level = 1;

    ProgressBar::init_func = ProgressBar::init_func_cmdline;
    ProgressBar::display_func = ProgressBar::display_func_cmdline;
    ProgressBar::done_func = ProgressBar::done_func_cmdline;

    print = cmdline_print;
    error = cmdline_error;
    info = cmdline_info;
    debug = cmdline_debug;


    sort_arguments (argc, argv); 

    srand (time (NULL));
      
    File::Config::init ();
  }


  App::~App () { }

  void App::run (int argc, gchar** argv)
  {
    parse_arguments ();
    execute ();
  }





  void App::sort_arguments (int argc, gchar** argv)
  {
    for (int n = 1; n < argc; n++) {
      const gchar* arg = argv[n];
      if (arg[0] == '-' && arg[1]) {

        while (*arg == '-') arg++;
        guint opt = match_option (arg);

        if (opt == G_MAXUINT) {
          throw Exception (String ("unknown option \"-") + arg + "\"");
        }
        else if (opt == DEFAULT_OPTIONS_OFFSET) {
          if (log_level < 2) log_level = 2;
        }
        else if (opt == DEFAULT_OPTIONS_OFFSET+1) {
          log_level = 0;
          ProgressBar::display = false;
        }
        else if (opt == DEFAULT_OPTIONS_OFFSET+2) {
          log_level = 3;
        }
        else if (opt == DEFAULT_OPTIONS_OFFSET+3) {
          print_help ();
          throw 0;
        }
        else if (opt == DEFAULT_OPTIONS_OFFSET+4) {
          std::printf ("%s %d.%d.%d\n  Author: %s\n  %s\n  using MRtrix %d.%d.%d, glib %d.%d.%d, GSL %s (build " __DATE__ ")\n",
              Glib::get_application_name().c_str(), 
              version[0], version[1], version[2], author, copyright,
              mrtrix_major_version, mrtrix_minor_version, mrtrix_micro_version, 
              glib_major_version, glib_minor_version, glib_micro_version, gsl_version);
          throw 0;
        }
        else {
          if (n + (int) command_options[opt].size() >= argc) {
            throw Exception (String ("not enough parameters to option \"-") + command_options[opt].sname + "\"");
          }

          parsed_options.push_back (ParsedOption());
          parsed_options.back().index = opt;
          while (parsed_options.back().args.size() < command_options[opt].size()) 
            parsed_options.back().args.push_back (argv[++n]);
        }
      }
      else parsed_arguments.push_back (argv[n]);
    }
  }






  void App::parse_arguments ()
  {
    guint num_args_required = 0, num_command_arguments = 0;
    bool has_optional_arguments = false;
    for (const Argument* arg = App::command_arguments; arg->is_valid(); arg++) {
      num_command_arguments++;
      if (arg->mandatory) num_args_required++; 
      if (arg->allow_multiple) has_optional_arguments = true;
    }

    if (has_optional_arguments && num_args_required > parsed_arguments.size()) 
      throw Exception ("expected at least " + str (num_args_required) + " arguments (" + str(parsed_arguments.size()) + " supplied)");
    
    if (!has_optional_arguments && num_args_required != parsed_arguments.size()) 
      throw Exception ("expected exactly " + str (num_args_required) + " arguments (" + str (parsed_arguments.size()) + " supplied)");

    guint optional_argument = G_MAXUINT;
    for (guint n = 0; n < parsed_arguments.size(); n++) {

      if (n < optional_argument) 
        if (!command_arguments[n].mandatory || command_arguments[n].allow_multiple) optional_argument = n;

      guint index = n;
      if (n >= optional_argument) {
        if ((int) (num_args_required - optional_argument) < (int) (parsed_arguments.size() - n)) index = optional_argument;
        else index = num_args_required - parsed_arguments.size() + n + (command_arguments[optional_argument].mandatory ? 0 : 1);
      }

      if (index >= num_command_arguments) throw Exception ("too many arguments");

      argument.push_back (ArgBase (command_arguments[index], parsed_arguments[n]));
      if (argument.back().type() == Undefined) 
        throw Exception (String ("error parsing argument \"") + command_arguments[index].sname + 
              "\" (specified as \"" + parsed_arguments[n] + "\")"); 
    }

    for (guint n = 0; n < parsed_options.size(); n++) {
      option.push_back (OptBase());
      option.back().index = parsed_options[n].index;
      for (guint a = 0; a < parsed_options[n].args.size(); a++) {
        ArgBase arg (command_options[parsed_options[n].index][a], parsed_options[n].args[a]);
        if (arg.type() == Undefined) 
          throw Exception (String ("error parsing argument \"") + command_options[parsed_options[n].index][a].sname 
                + "\" of option \"-" + command_options[parsed_options[n].index].sname
                + "\" (specified as \"" + parsed_options[n].args[a] + "\")");
        option.back().push_back (arg);
      }
    }

    for (guint index = 0; command_options[index].is_valid(); index++) {
      guint count = 0;
      for (guint n = 0; n < option.size(); n++)
        if (option[n].index == index)
          count++;

      if (command_options[index].mandatory && count < 1) 
        throw Exception (String ("mandatory option \"") + command_options[index].sname + "\" must be specified");

      if (!command_options[index].allow_multiple && count > 1) 
        throw Exception (String ("multiple instances of option \"") +  command_options[index].sname + "\" are not allowed");
    }

  }






  void App::print_help () const
  {
    fprintf (stderr, "%s: part of the MRtrix package\n\n", Glib::get_application_name().c_str());
    if (command_description[0]) {
      print_formatted_paragraph ("PURPOSE:", command_description[0], HELP_PURPOSE_INDENT);
      fprintf (stderr, "\n");
      for (const gchar** p = command_description+1; *p; p++) {
        print_formatted_paragraph ("", *p, HELP_PURPOSE_INDENT);
        fprintf (stderr, "\n");
      }
    }
    else fprintf (stderr, "(no description available)\n\n");

    fprintf (stderr, "%-*s%s [ options ]", HELP_INDENT, "SYNTAX:", Glib::get_application_name().c_str());
    for (const Argument* arg = command_arguments; arg->is_valid(); arg++) {
      if (!arg->mandatory) fprintf (stderr, " [");
      fprintf (stderr, " %s", arg->sname);
      if (arg->allow_multiple) {
        if (arg->mandatory) fprintf (stderr, " [ %s", arg->sname);
        fprintf (stderr, " ...");
      }
      if (!arg->mandatory || arg->allow_multiple) fprintf (stderr, " ]");
    }
    fprintf (stderr, "\n\n");



    for (const Argument* arg = command_arguments; arg->is_valid(); arg++) {
      print_formatted_paragraph (arg->sname, arg->desc, HELP_ARG_INDENT);
      fprintf (stderr, "\n");
    }


    fprintf (stderr, "OPTIONS:\n\n");
    for (const Option* opt = command_options; opt->is_valid(); opt++) {
      String text ("-");
      text += opt->sname;
      for (guint n = 0; n < opt->size(); n++) { text += " "; text += (*opt)[n].sname; }
      print_formatted_paragraph (text, opt->desc, HELP_OPTION_INDENT);
      for (guint n = 0; n < opt->size(); n++) {
        fprintf (stderr, "\n");
        print_formatted_paragraph ("", String ((*opt)[n].sname) + ": " + (*opt)[n].desc, HELP_OPTION_INDENT);
      }
      fprintf (stderr, "\n");
    }

    for (guint n = 0; n < NUM_DEFAULT_OPTIONS; n++) {
      String text ("-");
      text += default_options[n].sname;
      print_formatted_paragraph (text, default_options[n].desc, HELP_OPTION_INDENT);
      fprintf (stderr, "\n");
    }
  }




  void App::print_full_argument_usage (const Argument& arg) const
  {
    std::cout << "ARGUMENT " << arg.sname << " " << (arg.mandatory ? '1' : '0') << " " << (arg.allow_multiple ? '1' : '0') << " ";
    switch (arg.type) {
      case Integer: std::cout << "INT " << arg.extra_info.i.min << " " << arg.extra_info.i.max << " " << arg.extra_info.i.def; break;
      case Float: std::cout << "FLOAT " << arg.extra_info.f.min << " " << arg.extra_info.f.max << " " << arg.extra_info.f.def; break;
      case Text: std::cout << "TEXT"; if (arg.extra_info.string) std::cout << " " << arg.extra_info.string; break;
      case ArgFile: std::cout << "FILE"; break;
      case Choice: std::cout << "CHOICE"; for (const char** p = arg.extra_info.choice; *p; p++) std::cout << " " << *p; break;
      case ImageIn: std::cout << "IMAGEIN"; break;
      case ImageOut: std::cout << "IMAGEOUT"; break;
      case IntSeq: std::cout << "ISEQ"; break;
      case FloatSeq: std::cout << "FSEQ"; break;
      default:
                     throw (1);
    }
    std::cout << "\n" << arg.lname << "\n" << arg.desc << "\n";
  }






  void App::print_full_option_usage (const Option& opt) const
  {
    std::cout << "OPTION " << opt.sname << " " << (opt.mandatory ? '1' : '0') << " " << (opt.allow_multiple ? '1' : '0') << "\n";
    std::cout << opt.lname << "\n" << opt.desc << "\n";

    for (std::vector<Argument>::const_iterator i = opt.begin(); i != opt.end(); ++i) 
      print_full_argument_usage (*i);
  }






  void App::print_full_usage () const
  {
    for (const gchar** p = command_description; *p; p++) 
      std::cout << *p << "\n";

    for (const Argument* arg = command_arguments; arg->is_valid(); arg++) 
      print_full_argument_usage (*arg);

    for (const Option* opt = command_options; opt->is_valid(); opt++) 
      print_full_option_usage (*opt);

    for (guint n = 0; n < NUM_DEFAULT_OPTIONS; n++) 
      print_full_option_usage (default_options[n]);
  }






  std::ostream& operator<< (std::ostream& stream, const App& app)
  {
    stream 
      << "----------------------------------\n  COMMAND: " 
      << Glib::get_application_name()
      << "\n----------------------------------\n\n";

    const gchar** c = App::command_description;
    while (*c) { stream << *c << "\n\n"; c++; }

    stream << "ARGUMENTS:\n\n";
    for (guint n = 0; App::command_arguments[n].is_valid(); n++)
      stream << "[" << n << "] " << App::command_arguments[n] << "\n\n";

    stream << "OPTIONS:\n\n";
    for (guint n = 0; App::command_options[n].is_valid(); n++)
      stream << App::command_options[n] << "\n";

    return (stream);
  }



}


