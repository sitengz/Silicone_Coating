#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

constexpr double kAvogadroScale = 0.602; // Converts (g/mol)/(g/cm^3) to A^3.

struct Settings {
    int n1 = 384, m1 = 306; // Bifunctional network strands.
    int n2 = 32,  m2 = 153; // Tetrafunctional cross-linkers.
    int n3 = 32,  m3 = 425; // Neutral linear oil fillers.
    int n4 = 5,   m4 = 6;   // Five-bead star-like moderators.
    int crosslinker_functionality = 4;
    std::string crosslink_distribution = "random";
    std::uint32_t crosslink_seed = 20260722u;
    double bead_mass = 74.0;
    double density = 0.1;
    double target_density = 0.8;
    double bond_length = 2.801;
    double spacing = 7.0;
    double thickness = -1.0; // Negative selects the original cubic bulk box.
    std::uint32_t seed = 5489u;
    std::string output = "data.V35_PDMS_N32_10wt";
    bool output_explicit = false;
    bool filler_length_explicit = false;
    double filler_weight_percent = -1.0;
};

struct Atom { int id, molecule, type; double charge, x, y, z; };
struct Bond { int id, type, a, b; };
struct Angle { int id, type, a, b, c; };
struct Dihedral { int id, type, a, b, c, d; };
struct Box { double lx, ly, lz; };

struct System {
    std::vector<Atom> atoms;
    std::vector<Bond> bonds;
    std::vector<Angle> angles;
    std::vector<Dihedral> dihedrals;
};

int parse_int(const std::string& value, const std::string& option) {
    try { size_t used = 0; int result = std::stoi(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid integer for " + option + ": " + value); }
}

double parse_double(const std::string& value, const std::string& option) {
    try { size_t used = 0; double result = std::stod(value, &used); if (used != value.size()) throw std::invalid_argument(""); return result; }
    catch (...) { throw std::runtime_error("Invalid number for " + option + ": " + value); }
}

void print_help(const char* program) {
    std::cout
        << "Usage: " << program << " [options]\n\n"
        << "V35 LAMMPS data-file generator. N is beads per molecule; M is molecule count.\n\n"
        << "  --n1 N --m1 M    bifunctional network strands (defaults: 384, 306)\n"
        << "  --n2 N           cross-linker length (default: 32); M2 is stoichiometric\n"
        << "  --n3 N --m3 M    neutral linear oil fillers (defaults: 32, 425)\n"
        << "  --n4 N --m4 M    star-moderator size/count (N4 must be 5; defaults: 5, 6)\n"
        << "  --filler-length N  PDMS filler chain length (alias for --n3)\n"
        << "  --filler-wt X      PDMS weight percent; computes M3 from fixed V35 base\n"
        << "  --functionality F cross-linker reactive sites, 3-16 (default: 4)\n"
        << "  --crosslink-distribution MODE\n"
        << "                     reactive-site placement: regular or random (default: random)\n"
        << "  --crosslink-seed N random-site seed (default: 20260722)\n"
        << "  --mass X          bead mass in g/mol (default: 74)\n"
        << "  --density X       initial packing density in g/cm^3 (default: 0.1)\n"
        << "  --target-density X density after scripted compression (default: 0.8)\n"
        << "  --bond-length X   initial bond length in angstrom (default: 2.801)\n"
        << "  --spacing X       spacing between placed molecules (default: 7.0)\n"
        << "  --thickness X     fixed film thickness Lz in angstrom; omit for bulk\n"
        << "  --seed N          reproducible random seed (default: 5489)\n"
        << "  --output FILE     override the automatically generated data filename\n"
        << "                    companion input, submit, and .info files share its directory\n"
        << "  --help             show this help\n";
}

Settings parse_args(int argc, char** argv) {
    Settings s;
    for (int i = 1; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help") { print_help(argv[0]); std::exit(0); }
        if (i + 1 >= argc) throw std::runtime_error("Missing value after " + option);
        const std::string value = argv[++i];
        if      (option == "--n1") s.n1 = parse_int(value, option);
        else if (option == "--m1") s.m1 = parse_int(value, option);
        else if (option == "--n2") s.n2 = parse_int(value, option);
        else if (option == "--m2")
            throw std::runtime_error("M2 is determined by M2=2*M1/functionality; do not supply --m2");
        else if (option == "--n3") s.n3 = parse_int(value, option);
        else if (option == "--m3") s.m3 = parse_int(value, option);
        else if (option == "--n4") s.n4 = parse_int(value, option);
        else if (option == "--m4") s.m4 = parse_int(value, option);
        else if (option == "--functionality") s.crosslinker_functionality = parse_int(value, option);
        else if (option == "--crosslink-distribution") s.crosslink_distribution = value;
        else if (option == "--crosslink-seed") s.crosslink_seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--mass") s.bead_mass = parse_double(value, option);
        else if (option == "--density") s.density = parse_double(value, option);
        else if (option == "--target-density") s.target_density = parse_double(value, option);
        else if (option == "--bond-length") s.bond_length = parse_double(value, option);
        else if (option == "--spacing") s.spacing = parse_double(value, option);
        else if (option == "--thickness") s.thickness = parse_double(value, option);
        else if (option == "--seed") s.seed = static_cast<std::uint32_t>(parse_int(value, option));
        else if (option == "--output") { s.output = value; s.output_explicit = true; }
        else if (option == "--filler-length") { s.n3 = parse_int(value, option); s.filler_length_explicit = true; }
        else if (option == "--filler-wt") s.filler_weight_percent = parse_double(value, option);
        else throw std::runtime_error("Unknown option: " + option);
    }
    return s;
}

std::string filename_number(double value) {
    std::ostringstream text;
    text << std::fixed << std::setprecision(4) << value;
    std::string result = text.str();
    while (!result.empty() && result.back() == '0') result.pop_back();
    if (!result.empty() && result.back() == '.') result.pop_back();
    std::replace(result.begin(), result.end(), '.', 'p');
    return result;
}

void apply_filler_composition(Settings& s) {
    if (s.filler_length_explicit && s.filler_weight_percent < 0.0)
        s.filler_weight_percent = 10.0;
    if (s.filler_weight_percent < 0.0) return;
    if (s.filler_weight_percent >= 100.0)
        throw std::runtime_error("--filler-wt must be at least 0 and less than 100");
    if (s.n3 <= 0 && s.filler_weight_percent > 0.0)
        throw std::runtime_error("A positive --filler-wt requires a positive --filler-length");

    // Components 1, 2, and 4 form the fixed V35 base. With equal bead masses,
    // w = filler_beads / (base_beads + filler_beads).
    const long long base_beads = 1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n4*s.m4;
    if (base_beads <= 0) throw std::runtime_error("The fixed V35 base must contain at least one bead");
    if (s.filler_weight_percent == 0.0) {
        s.m3 = 0;
    } else {
        const double fraction = s.filler_weight_percent / 100.0;
        const double desired_filler_beads = base_beads * fraction / (1.0 - fraction);
        s.m3 = std::max(1, static_cast<int>(std::lround(desired_filler_beads / s.n3)));
    }
    if (!s.output_explicit)
        s.output = "data.V35_PDMS_N" + std::to_string(s.n3) + "_" + filename_number(s.filler_weight_percent) + "wt";
}

void apply_crosslinker_stoichiometry(Settings& s) {
    if (s.crosslinker_functionality < 3)
        throw std::runtime_error("Cross-linker functionality must be at least 3 to form a network");
    const long long reactive_ends = 2LL * s.m1;
    if (reactive_ends % s.crosslinker_functionality != 0) {
        std::ostringstream message;
        message << "Exact stoichiometry is impossible: 2*M1=" << reactive_ends
                << " is not divisible by functionality=" << s.crosslinker_functionality;
        throw std::runtime_error(message.str());
    }
    const long long molecule_count = reactive_ends / s.crosslinker_functionality;
    if (molecule_count > 100000000LL)
        throw std::runtime_error("Stoichiometric M2 is too large");
    s.m2 = static_cast<int>(molecule_count);
}

void apply_geometry_filename(Settings& s) {
    if (s.thickness > 0.0 && !s.output_explicit)
        s.output += "_film_Lz" + filename_number(s.thickness);
}

void report_composition(const Settings& s) {
    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    long long total_beads = 0;
    long long total_molecules = 0;
    for (size_t i = 0; i < 4; ++i) total_beads += 1LL * n[i] * m[i];
    for (int count : m) total_molecules += count;
    std::cerr << "Composition (equal bead masses):\n";
    for (size_t i = 0; i < 4; ++i) {
        const long long component_beads = 1LL * n[i] * m[i];
        const double realized = total_beads == 0 ? 0.0 : 100.0 * component_beads / total_beads;
        const double mole_fraction = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        std::cerr << "  component " << i + 1 << ": N=" << n[i] << ", M=" << m[i]
                  << ", mole%=" << std::fixed << std::setprecision(4) << mole_fraction
                  << ", realized wt%=" << realized;
        if (i == 2 && s.filler_weight_percent >= 0.0)
            std::cerr << ", requested filler wt%=" << s.filler_weight_percent;
        std::cerr << '\n';
    }
    std::cerr << "  total beads=" << total_beads << '\n';
}

void validate(const Settings& s) {
    const int values[] = {s.n1, s.n2, s.n3, s.n4, s.m1, s.m2, s.m3, s.m4};
    if (std::any_of(std::begin(values), std::end(values), [](int x) { return x < 0; }))
        throw std::runtime_error("N and M values cannot be negative");
    if (s.n4 != 5 && s.m4 != 0)
        throw std::runtime_error("The implemented moderator is a five-bead star; use --n4 5");
    if (s.crosslinker_functionality < 3 || s.crosslinker_functionality > 16 || s.crosslinker_functionality > s.n2)
        throw std::runtime_error("Cross-linker functionality must be between 3 and min(16, N2)");
    if (s.crosslink_distribution != "regular" && s.crosslink_distribution != "random")
        throw std::runtime_error("--crosslink-distribution must be regular or random");
    if (s.crosslink_distribution == "regular" && s.crosslinker_functionality > (s.n2 + 1) / 2)
        throw std::runtime_error("Regular placement at 1,3,5,... requires functionality <= ceil(N2/2)");
    if (s.bead_mass <= 0 || s.density <= 0 || s.target_density <= 0 ||
        s.bond_length <= 0 || s.spacing <= 0)
        throw std::runtime_error("Mass, densities, bond length, and spacing must be positive");
    if (s.thickness == 0.0)
        throw std::runtime_error("--thickness must be positive; omit it for the cubic bulk system");
    const long long total_beads = 1LL*s.n1*s.m1 + 1LL*s.n2*s.m2 + 1LL*s.n3*s.m3 + 1LL*s.n4*s.m4;
    if (total_beads <= 0) throw std::runtime_error("The formulation must contain at least one bead");
    if (s.output.empty()) throw std::runtime_error("Output filename cannot be empty");
}

std::set<int> regular_sites(int functionality) {
    std::set<int> sites;
    for (int i = 0; i < functionality; ++i) sites.insert(1 + 2*i);
    return sites;
}

std::set<int> random_sites(int bead_count, int functionality, std::uint32_t seed) {
    std::vector<int> candidates(static_cast<size_t>(bead_count));
    for (int i = 0; i < bead_count; ++i) candidates[static_cast<size_t>(i)] = i + 1;
    std::mt19937 rng(seed);
    std::shuffle(candidates.begin(), candidates.end(), rng);
    return std::set<int>(candidates.begin(), candidates.begin() + functionality);
}

std::set<int> crosslinker_sites(const Settings& s) {
    if (s.crosslink_distribution == "regular")
        return regular_sites(s.crosslinker_functionality);
    return random_sites(s.n2, s.crosslinker_functionality, s.crosslink_seed);
}

void add_linear_topology(System& sys, int first_atom, int bead_count) {
    for (int i = 0; i + 1 < bead_count; ++i)
        sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 1, first_atom + i, first_atom + i + 1});
    for (int i = 0; i + 2 < bead_count; ++i)
        sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2});
    for (int i = 0; i + 3 < bead_count; ++i)
        sys.dihedrals.push_back({static_cast<int>(sys.dihedrals.size()) + 1, 1, first_atom + i, first_atom + i + 1, first_atom + i + 2, first_atom + i + 3});
}


struct PathPoint { double x, y; };

std::vector<PathPoint> folded_strand_path(int bead_count, const Box& box,
                                          double bond_length, double spacing) {
    std::vector<PathPoint> path;
    if (bead_count <= 0) return path;
    const double available_x = box.lx - 2.0 * spacing;
    if (available_x < bond_length)
        throw std::runtime_error("Lx is too small to place a folded network strand");
    const int horizontal_steps_max =
        std::max(1, static_cast<int>(std::floor(available_x / bond_length)));
    const int connector_steps_max =
        std::max(1, static_cast<int>(std::ceil(spacing / bond_length)));

    path.reserve(static_cast<size_t>(bead_count));
    double x = 0.0, y = 0.0;
    int direction = 1;
    int horizontal_steps = 0;
    int connector_steps = 0;
    bool connecting_rows = false;

    for (int bead = 0; bead < bead_count; ++bead) {
        path.push_back({x, y});
        if (bead + 1 == bead_count) break;
        if (!connecting_rows && horizontal_steps < horizontal_steps_max) {
            x += direction * bond_length;
            ++horizontal_steps;
        } else if (!connecting_rows) {
            y += bond_length;
            connecting_rows = true;
            connector_steps = 1;
        } else if (connector_steps < connector_steps_max) {
            y += bond_length;
            ++connector_steps;
        } else {
            direction = -direction;
            x += direction * bond_length;
            horizontal_steps = 1;
            connector_steps = 0;
            connecting_rows = false;
        }
    }
    return path;
}

void add_network_strands(System& sys, const Settings& s, const Box& box) {
    if (s.m1 == 0) return;
    if (s.n1 <= 0) throw std::runtime_error("A nonzero M1 requires N1 > 0");
    if (box.ly < 2.0*s.spacing || box.lz < 2.0*s.spacing)
        throw std::runtime_error("Ly and Lz must each be at least twice the placement spacing");

    const std::vector<PathPoint> path =
        folded_strand_path(s.n1, box, s.bond_length, s.spacing);
    double footprint_y = 0.0;
    for (const PathPoint& point : path) footprint_y = std::max(footprint_y, point.y);
    if (footprint_y > box.ly - 2.0*s.spacing)
        throw std::runtime_error("Ly is too small for the folded network-strand footprint");

    const int ny = std::max(1, static_cast<int>(
        std::floor((box.ly - 2.0*s.spacing) / (footprint_y + s.spacing))) + 1);
    const int nz = std::max(1, static_cast<int>(box.lz / s.spacing));
    const long long capacity = 1LL * ny * nz;
    if (capacity < s.m1) {
        std::ostringstream msg;
        msg << "Folded network-strand placement capacity (" << capacity
            << ") is smaller than M1 (" << s.m1 << "). Reduce spacing.";
        throw std::runtime_error(msg.str());
    }

    for (int molecule_index = 0; molecule_index < s.m1; ++molecule_index) {
        const int iy = molecule_index % ny;
        const int iz = molecule_index / ny;
        const bool mirror_x = molecule_index % 2 != 0;
        const double y0 = -box.ly/2 + s.spacing +
                          iy * (footprint_y + s.spacing);
        const double z = -box.lz/2 + s.spacing + iz * s.spacing;
        const int molecule_id =
            static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;

        for (int bead = 1; bead <= s.n1; ++bead) {
            const PathPoint& point = path[static_cast<size_t>(bead - 1)];
            const double x = mirror_x
                ? box.lx/2 - s.spacing - point.x
                : -box.lx/2 + s.spacing + point.x;
            const int atom_type = (bead == 1 || bead == s.n1) ? 2 : 1;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1,
                                 molecule_id, atom_type, 0.0, x, y0 + point.y, z});
        }
        add_linear_topology(sys, first_atom, s.n1);
    }
}

void add_linear_component(System& sys, int bead_count, int molecule_count, int component,
                          const Box& box, double bond_length, double spacing,
                          const std::set<int>& reactive_sites = {}) {
    if (bead_count == 0 && molecule_count != 0) throw std::runtime_error("A nonzero molecule count requires N > 0");
    if (molecule_count == 0) return;

    const double span = std::max(0, bead_count - 1) * bond_length;
    if (box.lx < span + 2.0*spacing)
        throw std::runtime_error("Lx is too small for a straight chain plus placement margins");
    if (box.ly < 2.0*spacing || box.lz < 2.0*spacing)
        throw std::runtime_error("Ly and Lz must each be at least twice the placement spacing");
    const int nx = std::max(1, static_cast<int>(box.lx / (span + spacing)));
    const int ny = std::max(1, static_cast<int>(box.ly / spacing));
    const int nz = component == 3
        ? std::max(1, static_cast<int>((box.lz / 2.0 - spacing) / spacing) + 1)
        : std::max(1, static_cast<int>(box.lz / spacing));
    const long long capacity = 1LL * nx * ny * nz;
    if (capacity < molecule_count) {
        std::ostringstream msg;
        msg << "Placement grid capacity (" << capacity << ") is smaller than M" << component
            << " (" << molecule_count << "). Reduce spacing or molecule length.";
        throw std::runtime_error(msg.str());
    }

    for (int molecule_index = 0; molecule_index < molecule_count; ++molecule_index) {
        const int ix = molecule_index % nx;
        const int iy = (molecule_index / nx) % ny;
        const int iz = molecule_index / (nx * ny);
        const bool reverse = component != 1;
        double x = reverse ? box.lx / 2 - spacing - ix * (span + spacing)
                           : -box.lx / 2 + spacing + ix * (span + spacing);
        const double y = (reverse ? box.ly / 2 - spacing : -box.ly / 2 + spacing) + (reverse ? -1 : 1) * iy * spacing;
        double z = 0.0;
        if (component == 2) z = box.lz / 2 - spacing - iz * spacing;
        else if (component == 3) z = iz * spacing; // Filler starts at the midplane and grows toward +z.
        else z = -box.lz / 2 + spacing + iz * spacing;

        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int first_atom = static_cast<int>(sys.atoms.size()) + 1;
        for (int bead = 1; bead <= bead_count; ++bead) {
            int atom_type = 1;
            if (component == 1 && (bead == 1 || bead == bead_count)) atom_type = 2;
            if (component == 2 && reactive_sites.count(bead)) atom_type = 3;
            sys.atoms.push_back({static_cast<int>(sys.atoms.size()) + 1, molecule_id, atom_type, 0.0, x, y, z});
            x += reverse ? -bond_length : bond_length;
        }
        add_linear_topology(sys, first_atom, bead_count);
    }
}

void add_star_moderators(System& sys, const Settings& s, const Box& box, std::mt19937& rng) {
    const double margin = s.bond_length;
    if (box.lx <= 2.0*margin || box.ly <= 2.0*margin)
        throw std::runtime_error("Lateral box dimensions are too small for a star moderator");
    std::uniform_real_distribution<double> xdist(-box.lx/2 + margin, box.lx/2 - margin);
    std::uniform_real_distribution<double> ydist(-box.ly/2 + margin, box.ly/2 - margin);
    std::uniform_real_distribution<double> zdist(0.02 * box.lz, 0.12 * box.lz); // Preserves V35's original slab.
    for (int i = 0; i < s.m4; ++i) {
        const int molecule_id = static_cast<int>(sys.atoms.empty() ? 1 : sys.atoms.back().molecule + 1);
        const int center = static_cast<int>(sys.atoms.size()) + 1;
        const double x = xdist(rng), y = ydist(rng), z = zdist(rng), b = s.bond_length;
        // Preserve the legacy force-field mapping: ordinary center, type-2 arms.
        sys.atoms.push_back({center,     molecule_id, 1, 0.0, x,     y,     z});
        sys.atoms.push_back({center + 1, molecule_id, 2, 0.0, x + b, y,     z});
        sys.atoms.push_back({center + 2, molecule_id, 2, 0.0, x - b, y,     z});
        sys.atoms.push_back({center + 3, molecule_id, 2, 0.0, x,     y + b, z});
        sys.atoms.push_back({center + 4, molecule_id, 2, 0.0, x,     y - b, z});
        for (int arm = 1; arm <= 4; ++arm)
            sys.bonds.push_back({static_cast<int>(sys.bonds.size()) + 1, 2, center, center + arm});
        for (int a = 1; a <= 4; ++a)
            for (int c = a + 1; c <= 4; ++c)
                sys.angles.push_back({static_cast<int>(sys.angles.size()) + 1, 1, center + a, center, center + c});
    }
}

System build_system(const Settings& s, const Box& box) {
    System sys;
    const long long expected_atoms = 1LL * s.n1*s.m1 + 1LL * s.n2*s.m2 + 1LL * s.n3*s.m3 + 1LL * s.n4*s.m4;
    if (expected_atoms > 100000000LL) throw std::runtime_error("Requested system is too large");
    sys.atoms.reserve(static_cast<size_t>(expected_atoms));
    add_network_strands(sys, s, box);
    const std::set<int> reactive_sites = crosslinker_sites(s);
    std::cerr << "Cross-linker type-3 sites (" << s.crosslink_distribution << "):";
    for (int site : reactive_sites) std::cerr << ' ' << site;
    std::cerr << '\n';
    add_linear_component(sys, s.n2, s.m2, 2, box, s.bond_length, s.spacing,
                         reactive_sites);
    add_linear_component(sys, s.n3, s.m3, 3, box, s.bond_length, s.spacing);
    std::mt19937 rng(s.seed);
    add_star_moderators(sys, s, box, rng);
    if (static_cast<long long>(sys.atoms.size()) != expected_atoms)
        throw std::logic_error("Internal error: generated atom count does not match requested composition");
    return sys;
}

void write_data(const Settings& s, const System& sys, const Box& box) {
    std::ofstream out(s.output);
    if (!out) throw std::runtime_error("Cannot open output file: " + s.output);
    out << "LAMMPS data file for the V35 four-component formulation\n\n"
        << sys.atoms.size() << " atoms\n4 atom types\n"
        << sys.bonds.size() << " bonds\n2 bond types\n"
        << sys.angles.size() << " angles\n1 angle types\n"
        << sys.dihedrals.size() << " dihedrals\n1 dihedral types\n\n"
        << std::fixed << std::setprecision(6)
        << -box.lx/2 << ' ' << box.lx/2 << " xlo xhi\n"
        << -box.ly/2 << ' ' << box.ly/2 << " ylo yhi\n"
        << -box.lz/2 << ' ' << box.lz/2 << " zlo zhi\n\nMasses\n\n";
    for (int type = 1; type <= 4; ++type) out << type << ' ' << s.bead_mass << '\n';
    out << "\nAtoms # full\n\n" << std::setprecision(3);
    for (const auto& a : sys.atoms)
        out << a.id << ' ' << a.molecule << ' ' << a.type << ' ' << a.charge << ' '
            << a.x << ' ' << a.y << ' ' << a.z << '\n';
    if (!sys.bonds.empty()) {
        out << "\nBonds\n\n";
        for (const auto& b : sys.bonds) out << b.id << ' ' << b.type << ' ' << b.a << ' ' << b.b << '\n';
    }
    if (!sys.angles.empty()) {
        out << "\nAngles\n\n";
        for (const auto& a : sys.angles) out << a.id << ' ' << a.type << ' ' << a.a << ' ' << a.b << ' ' << a.c << '\n';
    }
    if (!sys.dihedrals.empty()) {
        out << "\nDihedrals\n\n";
        for (const auto& d : sys.dihedrals) out << d.id << ' ' << d.type << ' ' << d.a << ' ' << d.b << ' ' << d.c << ' ' << d.d << '\n';
    }
    if (!out) throw std::runtime_error("Failed while writing output file: " + s.output);
}

struct OutputFiles {
    std::string data;
    std::string data_basename;
    std::string input;
    std::string input_basename;
    std::string submit;
    std::string submit_basename;
    std::string info;
    std::string info_basename;
    std::string case_name;
};

struct LjParameters {
    double epsilon;
    double sigma;
    double cutoff;
};

std::string basename_of(const std::string& path) {
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? path : path.substr(slash + 1);
}

std::string directory_of(const std::string& path) {
    const size_t slash = path.find_last_of("/\\");
    return slash == std::string::npos ? "" : path.substr(0, slash + 1);
}

OutputFiles output_files(const Settings& s) {
    OutputFiles files;
    files.data = s.output;
    files.data_basename = basename_of(s.output);
    if (files.data_basename.empty())
        throw std::runtime_error("Output filename cannot end with a directory separator");
    files.case_name = files.data_basename.rfind("data.", 0) == 0
        ? files.data_basename.substr(5) : files.data_basename;
    if (files.case_name.empty())
        throw std::runtime_error("Cannot derive a case name from the output filename");
    const std::string directory = directory_of(s.output);
    files.input_basename = "in." + files.case_name;
    files.submit_basename = "submit." + files.case_name + ".sh";
    files.info_basename = files.case_name + ".info";
    files.input = directory + files.input_basename;
    files.submit = directory + files.submit_basename;
    files.info = directory + files.info_basename;
    return files;
}

LjParameters lj_parameters(double temperature) {
    const double a = temperature - 186.04682;
    const double b = 0.00758 * a;
    const double c = 1.0 + std::exp(b);
    const double epsilon = (4.77795 / c + 1.47169) * 0.350646;
    const double sigma = ((7.86548e-05) * temperature + 1.27856) * 4.95013;
    return {epsilon, sigma, sigma * std::pow(2.0, 1.0 / 6.0)};
}

std::string sanitize_job_name(const std::string& value) {
    std::string result;
    for (char c : value) {
        const unsigned char uc = static_cast<unsigned char>(c);
        result.push_back(std::isalnum(uc) || c == '_' || c == '-' ? c : '_');
    }
    if (result.size() > 100) result.resize(100);
    return result.empty() ? "V35" : result;
}

std::string shell_single_quote(const std::string& value) {
    std::string result = "'";
    for (char c : value) {
        if (c == '\'') result += "'\\''";
        else result.push_back(c);
    }
    result += "'";
    return result;
}

std::string json_escape(const std::string& value) {
    std::ostringstream escaped;
    for (unsigned char c : value) {
        switch (c) {
            case '"': escaped << "\\\""; break;
            case '\\': escaped << "\\\\"; break;
            case '\b': escaped << "\\b"; break;
            case '\f': escaped << "\\f"; break;
            case '\n': escaped << "\\n"; break;
            case '\r': escaped << "\\r"; break;
            case '\t': escaped << "\\t"; break;
            default:
                if (c < 0x20)
                    escaped << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                            << static_cast<int>(c) << std::dec << std::setfill(' ');
                else
                    escaped << static_cast<char>(c);
        }
    }
    return escaped.str();
}

void write_lammps_input(const Settings& s, const OutputFiles& files) {
    std::ofstream out(files.input);
    if (!out) throw std::runtime_error("Cannot open LAMMPS input file: " + files.input);

    const bool film = s.thickness > 0.0;
    const double compression_scale = film
        ? std::sqrt(s.density / s.target_density)
        : std::cbrt(s.density / s.target_density);
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const std::string suffix = files.case_name;

    out << "# Generated by the V35 model generator\n"
        << "# Geometry: " << (film ? "film with fixed Lz" : "periodic bulk") << "\n\n"
        << "units           real\n"
        << "boundary        p p " << (film ? "f" : "p") << "\n"
        << "atom_style      full\n"
        << "bond_style      harmonic\n"
        << "angle_style     harmonic\n"
        << "dihedral_style  nharmonic\n"
        << std::fixed << std::setprecision(8)
        << "pair_style      lj/cut " << hot.cutoff << "\n"
        << "comm_modify     cutoff 15\n"
        << "read_data       " << files.data_basename
        << " extra/bond/per/atom 4 extra/angle/per/atom 10"
        << " extra/dihedral/per/atom 10 extra/special/per/atom 30\n\n"
        << "mass            1 " << s.bead_mass << "\n"
        << "mass            2 " << s.bead_mass << "\n"
        << "mass            3 " << s.bead_mass << "\n"
        << "mass            4 " << s.bead_mass << "\n\n"
        << "bond_coeff      1 115.4086 2.801\n"
        << "bond_coeff      2 115.4086 7.235\n"
        << "angle_coeff     1 64.62431 111.623\n"
        << "dihedral_coeff  1 4 3.280141429 -0.59019769 1.991530534 3.31026047\n\n"
        << "# Calibrated 800 K repulsive pair interaction\n"
        << "pair_coeff      * * " << hot.epsilon << ' ' << hot.sigma << "\n\n"
        << "neighbor        2 bin\n"
        << "neigh_modify    delay 5 every 1\n"
        << "timestep        5\n"
        << "thermo          1000\n"
        << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n"
        << "restart         100000 restart." << suffix << ".1 restart." << suffix << ".2\n"
        << "dump            traj all custom 100000 dump." << suffix
        << ".lammpstrj id mol type q x y z ix iy iz\n"
        << "dump_modify     traj format line \"%d %d %d %.1f %.3f %.3f %.3f %d %d %d\" sort id\n\n";

    if (film) {
        out << "# Repulsive walls follow the fixed lower and upper z box edges.\n"
            << "fix             zwalls all wall/lj126 zlo EDGE zhi EDGE "
            << hot.epsilon << ' ' << hot.sigma << ' ' << hot.cutoff << " units box\n\n";
    }

    out << "# 800 K relaxation\n"
        << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
        << "run             1000000\n"
        << "write_data      data." << suffix << ".rep_800 nocoeff\n\n"
        << "unfix           integrate\n\n";

    if (film) {
        out << "# Lateral compression with crosslinking at fixed film thickness\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob 0.1 348154\n"
            << "thermo_style    custom step temp density etotal epair ebond eangle"
            << " edihed f_xlink[1] f_xlink[2] bonds\n"
            << "fix             compress all deform 1 x scale " << compression_scale
            << " y scale " << compression_scale << " units box\n"
            << "run             1000000\n"
            << "unfix           compress\n\n"
            << "# Relax at the compressed dimensions with continued crosslinking\n"
            << "run             1000000\n\n"
            << "# Continue crosslinking at fixed dimensions\n"
            << "run             2000000\n"
            << "unfix           integrate\n"
            << "unfix           xlink\n";
    } else {
        out << "# Isotropic compression with crosslinking at 800 K\n"
            << "fix             integrate all nvt temp 800.0 800.0 50.0\n"
            << "fix             xlink all bond/create 1 2 3 " << hot.cutoff
            << " 2 iparam 1 1 jparam 1 1 prob 0.1 348154\n"
            << "thermo_style    custom step temp density etotal epair ebond eangle"
            << " edihed f_xlink[1] f_xlink[2] bonds\n"
            << "fix             compress all deform 1 x scale " << compression_scale
            << " y scale " << compression_scale << " z scale " << compression_scale
            << " units box\n"
            << "run             1000000\n"
            << "unfix           compress\n\n"
            << "# Relax at the compressed dimensions with continued crosslinking\n"
            << "run             1000000\n\n"
            << "# Continue crosslinking at fixed dimensions\n"
            << "run             2000000\n"
            << "unfix           integrate\n"
            << "unfix           xlink\n";
    }

    out << "write_data      data." << suffix << ".xlink_800 nocoeff\n\n"
        << "# Restore the reacted-bond equilibrium length and 300 K pair interaction\n"
        << "bond_coeff      2 115.4086 2.801\n"
        << "pair_style      lj/gromacs 12 15\n"
        << "pair_coeff      * * " << cold.epsilon << ' ' << cold.sigma << "\n";

    if (film) {
        out << "unfix           zwalls\n"
            << "fix             zwalls all wall/lj126 zlo EDGE zhi EDGE "
            << cold.epsilon << ' ' << cold.sigma << ' ' << cold.cutoff << " units box\n";
    }

    out << "thermo_style    custom step temp density lx ly lz pxx pyy pzz"
        << " etotal epair ebond eangle edihed\n\n"
        << "# Cool from 800 K to 300 K\n";

    if (film) {
        out << "fix             integrate all npt temp 800.0 300.0 50.0"
            << " x 1.0 1.0 500.0 y 1.0 1.0 500.0 couple xy\n";
    } else {
        out << "fix             integrate all npt temp 800.0 300.0 50.0"
            << " iso 1.0 1.0 500.0\n";
    }

    out << "run             1000000\n"
        << "write_data      data." << suffix << ".300 nocoeff\n"
        << "unfix           integrate\n\n"
        << "# Final 300 K equilibration\n";

    if (film) {
        out << "fix             integrate all npt temp 300.0 300.0 50.0"
            << " x 1.0 1.0 500.0 y 1.0 1.0 500.0 couple xy\n";
    } else {
        out << "fix             integrate all npt temp 300.0 300.0 50.0"
            << " iso 1.0 1.0 500.0\n";
    }

    out << "run             1000000\n"
        << "write_data      data." << suffix << ".npt_eq nocoeff\n";

    if (!out) throw std::runtime_error("Failed while writing LAMMPS input file: " + files.input);
}

void write_submit_script(const OutputFiles& files) {
    std::ofstream out(files.submit);
    if (!out) throw std::runtime_error("Cannot open Slurm submit file: " + files.submit);
    out << "#!/bin/bash\n"
        << "#SBATCH --job-name=" << sanitize_job_name(files.case_name) << "_xlink\n"
        << "#SBATCH --time=48:00:00\n"
        << "#SBATCH --nodes=1\n"
        << "#SBATCH --ntasks-per-node=48\n"
        << "#SBATCH --mem=200G\n"
        << "#SBATCH --partition=nova\n"
        << "#SBATCH --mail-user=siteng@iastate.edu\n"
        << "#SBATCH --mail-type=END,FAIL\n"
        << "#SBATCH --output=slurm-%j.out\n"
        << "#SBATCH --error=slurm-%j.err\n\n"
        << "set -euo pipefail\n"
        << "cd -- \"${SLURM_SUBMIT_DIR:?SLURM_SUBMIT_DIR is not set}\"\n\n"
        << "module purge\n"
        << "module load intel/22.3.1\n"
        << "module load mpi/2021.7.1\n"
        << "module load lammps/20230802.2-py310-openmpi4-ezoqd7f\n\n"
        << "export OMP_NUM_THREADS=1\n\n"
        << "INPUT=" << shell_single_quote(files.input_basename) << "\n"
        << "OUTPUT=" << shell_single_quote("out." + files.case_name) << "\n"
        << "srun lmp -in \"$INPUT\" > \"$OUTPUT\"\n";
    if (!out) throw std::runtime_error("Failed while writing Slurm submit file: " + files.submit);
}

void write_info(const Settings& s, const System& sys, const Box& box,
                const OutputFiles& files) {
    std::ofstream out(files.info);
    if (!out) throw std::runtime_error("Cannot open model info file: " + files.info);

    const int n[] = {s.n1, s.n2, s.n3, s.n4};
    const int m[] = {s.m1, s.m2, s.m3, s.m4};
    const char* names[] = {"network_strands", "crosslinkers", "neutral_pdms_filler",
                           "star_moderators"};
    long long total_beads = 0;
    long long total_molecules = 0;
    for (int i = 0; i < 4; ++i) {
        total_beads += 1LL * n[i] * m[i];
        total_molecules += m[i];
    }
    const double volume = box.lx * box.ly * box.lz;
    const double realized_filler_wt = total_beads == 0
        ? 0.0 : 100.0 * (1LL * s.n3 * s.m3) / total_beads;
    const double compression_scale = s.thickness > 0.0
        ? std::sqrt(s.density / s.target_density)
        : std::cbrt(s.density / s.target_density);
    const LjParameters hot = lj_parameters(800.0);
    const LjParameters cold = lj_parameters(300.0);
    const std::set<int> sites = crosslinker_sites(s);

    out << std::fixed << std::setprecision(8)
        << "{\n"
        << "  \"format\": \"V35-model-info\",\n"
        << "  \"format_version\": 1,\n"
        << "  \"formulation\": \"V35\",\n"
        << "  \"case_name\": \"" << json_escape(files.case_name) << "\",\n"
        << "  \"geometry\": \"" << (s.thickness > 0.0 ? "film" : "bulk") << "\",\n"
        << "  \"film_thickness_angstrom\": ";
    if (s.thickness > 0.0) out << s.thickness;
    else out << "null";
    out << ",\n  \"film_thickness_source\": ";
    if (s.thickness > 0.0)
        out << "\"command-line value captured from the equilibrated bulk result\"";
    else
        out << "null";
    out << ",\n"
        << "  \"files\": {\n"
        << "    \"data\": \"" << json_escape(files.data_basename) << "\",\n"
        << "    \"lammps_input\": \"" << json_escape(files.input_basename) << "\",\n"
        << "    \"slurm_submit\": \"" << json_escape(files.submit_basename) << "\",\n"
        << "    \"model_info\": \"" << json_escape(files.info_basename) << "\"\n"
        << "  },\n"
        << "  \"components\": {\n";
    for (int i = 0; i < 4; ++i) {
        const long long beads = 1LL * n[i] * m[i];
        const double weight_percent = total_beads == 0 ? 0.0 : 100.0 * beads / total_beads;
        const double molecule_percent = total_molecules == 0 ? 0.0 : 100.0 * m[i] / total_molecules;
        out << "    \"" << names[i] << "\": {\"N\": " << n[i] << ", \"M\": " << m[i]
            << ", \"beads\": " << beads << ", \"weight_percent\": " << weight_percent
            << ", \"molecule_percent\": " << molecule_percent << "}"
            << (i == 3 ? "\n" : ",\n");
    }
    out << "  },\n"
        << "  \"composition\": {\n"
        << "    \"total_beads\": " << total_beads << ",\n"
        << "    \"total_molecules\": " << total_molecules << ",\n"
        << "    \"requested_filler_weight_percent\": ";
    if (s.filler_weight_percent >= 0.0) out << s.filler_weight_percent;
    else out << "null";
    out << ",\n"
        << "    \"realized_filler_weight_percent\": " << realized_filler_wt << "\n"
        << "  },\n"
        << "  \"crosslinker\": {\n"
        << "    \"functionality\": " << s.crosslinker_functionality << ",\n"
        << "    \"distribution\": \"" << json_escape(s.crosslink_distribution) << "\",\n"
        << "    \"reactive_bead_sites\": [";
    bool first = true;
    for (int site : sites) {
        if (!first) out << ", ";
        out << site;
        first = false;
    }
    out << "],\n"
        << "    \"type2_sites_total_including_moderators\": " << 2LL*s.m1 + 4LL*s.m4 << ",\n"
        << "    \"type3_sites_total\": " << 1LL*s.m2*s.crosslinker_functionality << "\n"
        << "  },\n"
        << "  \"initial_state\": {\n"
        << "    \"bead_mass_g_per_mol\": " << s.bead_mass << ",\n"
        << "    \"density_g_per_cm3\": " << s.density << ",\n"
        << "    \"box_angstrom\": {\"Lx\": " << box.lx << ", \"Ly\": " << box.ly
        << ", \"Lz\": " << box.lz << "},\n"
        << "    \"volume_angstrom3\": " << volume << ",\n"
        << "    \"bond_length_angstrom\": " << s.bond_length << ",\n"
        << "    \"placement_spacing_angstrom\": " << s.spacing << ",\n"
        << "    \"network_strand_initial_shape\": \"folded_serpentine\"\n"
        << "  },\n"
        << "  \"topology\": {\n"
        << "    \"atoms\": " << sys.atoms.size() << ",\n"
        << "    \"bonds\": " << sys.bonds.size() << ",\n"
        << "    \"angles\": " << sys.angles.size() << ",\n"
        << "    \"dihedrals\": " << sys.dihedrals.size() << "\n"
        << "  },\n"
        << "  \"random_seeds\": {\n"
        << "    \"crosslink_site_seed\": " << s.crosslink_seed << ",\n"
        << "    \"star_moderator_seed\": " << s.seed << ",\n"
        << "    \"bond_creation_seed\": 348154\n"
        << "  },\n"
        << "  \"simulation_template\": {\n"
        << "    \"target_compressed_density_g_per_cm3\": " << s.target_density << ",\n"
        << "    \"compression_scale_per_deformed_dimension\": " << compression_scale << ",\n"
        << "    \"hot_temperature_K\": 800.0,\n"
        << "    \"final_temperature_K\": 300.0,\n"
        << "    \"wall_style\": ";
    if (s.thickness > 0.0) out << "\"lj126 at zlo/zhi EDGE\"";
    else out << "null";
    out << ",\n"
        << "    \"hot_lj\": {\"epsilon\": " << hot.epsilon << ", \"sigma\": " << hot.sigma
        << ", \"cutoff\": " << hot.cutoff << "},\n"
        << "    \"cold_lj\": {\"epsilon\": " << cold.epsilon << ", \"sigma\": " << cold.sigma
        << ", \"cutoff\": " << cold.cutoff << "},\n"
        << "    \"timestep_fs\": 5.0,\n"
        << "    \"total_run_steps\": 7000000,\n"
        << "    \"bond_creation_active_steps\": 4000000\n"
        << "  }\n"
        << "}\n";
    if (!out) throw std::runtime_error("Failed while writing model info file: " + files.info);
}

} // namespace

int main(int argc, char** argv) {
    try {
        Settings settings = parse_args(argc, argv);
        apply_crosslinker_stoichiometry(settings);
        apply_filler_composition(settings);
        apply_geometry_filename(settings);
        validate(settings);
        report_composition(settings);
        const long long beads = 1LL*settings.n1*settings.m1 + 1LL*settings.n2*settings.m2
                              + 1LL*settings.n3*settings.m3 + 1LL*settings.n4*settings.m4;
        const double volume = beads * settings.bead_mass / (settings.density * kAvogadroScale);
        Box box{};
        if (settings.thickness > 0.0) {
            box.lz = settings.thickness;
            box.lx = box.ly = std::sqrt(volume / box.lz);
        } else {
            box.lx = box.ly = box.lz = std::cbrt(volume);
        }
        const System system = build_system(settings, box);
        const OutputFiles files = output_files(settings);
        write_data(settings, system, box);
        write_lammps_input(settings, files);
        write_submit_script(files);
        write_info(settings, system, box, files);
        std::cerr << "Wrote model package:\n"
                  << "  " << files.data << "\n"
                  << "  " << files.input << "\n"
                  << "  " << files.submit << "\n"
                  << "  " << files.info << "\n"
                  << "System: " << system.atoms.size() << " atoms, "
                  << system.bonds.size() << " bonds, " << system.angles.size() << " angles, "
                  << system.dihedrals.size() << " dihedrals; box "
                  << box.lx << " x " << box.ly << " x " << box.lz << " A\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}

