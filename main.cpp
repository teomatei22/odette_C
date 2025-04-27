#include <iostream>
#include "TDMParser.h"
#include "interpret.h"
#include "orbmath.h"
#include "popl.h"
// --help
// --method [sgp4 | gauss|
// --tle-json [TLE file .json for json, otherwise plaintext]
// --tdm [for gauss, format example in data folder]
// --output [output file]
// --format [csv|json|etc]
// --visualize [starts a 3d window to see the orbit

#define HELP_ME \
"Odette - Orbit Determination Technologies\n" \
"-----------------------------------------\n" \
"\n" \
"Usage:\n" \
"    odette [OPTIONS]\n" \
"\n" \
"Description:\n" \
"    Odette is a command-line utility for orbit determination and propagation.\n" \
"    It supports both SGP4 (TLE-based) and Gauss (angle-only) methods, accepts\n" \
"    various input formats, and can optionally visualize computed orbits.\n" \
"\n" \
"Options:\n" \
"    --method <sgp4|gauss>\n" \
"        Orbit determination method:\n" \
"          sgp4   - Uses Two-Line Element sets (SGP4 propagation)\n" \
"          gauss  - Computes orbit from RA/Dec observations (Gauss's method)\n" \
"\n" \
"    --tle-json <file>\n" \
"        Path to TLE input file. If the file ends in '.json', it will be parsed\n" \
"        as a JSON array. Otherwise, it is assumed to be plaintext TLE format.\n" \
"\n" \
"    --tdm <file>\n" \
"        Path to a TDM (Tracking Data Message) file. Required for Gauss method.\n" \
"\n" \
"    --output <file>\n" \
"        Output file for orbit results (ephemeris or state vectors).\n" \
"\n" \
"    --format <csv|json|...>\n" \
"        Format for the output file. Options include:\n" \
"          csv   - Comma-separated values\n" \
"          json  - JSON structured format\n" \
"\n" \
"    --visualize\n" \
"        Launches an interactive 3D visualization window of the orbit.\n" \
"\n" \
"    -h, --help\n" \
"        Display this help message and exit.\n" \
"\n" \
"Examples:\n" \
"    odette --method sgp4 --tle-json tle.txt --output orbit.csv --format csv\n" \
"    odette --method gauss --tdm obs.tdm --output orbit.json --format json --visualize\n" \
"\n" \
"For more information, visit the project repository or documentation.\n"


int main(int argc, char* argv[]) {
    using namespace popl;

    OptionParser op("Allowed options");
    auto help_option   = op.add<Switch>("h", "help", HELP_ME);
    op.parse(argc, argv);

    // print auto-generated help message
    if (help_option->is_set())
        std::cout << op << "\n";


}