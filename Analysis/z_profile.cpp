#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::size_t kComponentCount = 4;

enum ComponentIndex : std::size_t {
    kNetwork = 0,
    kCrosslinker = 1,
    kFiller = 2,
    kModerator = 3
};

struct ComponentInfo {
    std::string key;
    long long molecules = 0;
    long long expected_beads = 0;
};

struct ModelInfo {
    std::string case_name;
    std::string geometry;
    std::array<ComponentInfo, kComponentCount> components;
};

struct DataHeader {
    long long atoms = 0;
    int atom_types = 0;
    double xlo = 0.0;
    double xhi = 0.0;
    double ylo = 0.0;
    double yhi = 0.0;
    double zlo = 0.0;
    double zhi = 0.0;
    bool have_x = false;
    bool have_y = false;
    bool have_z = false;
    std::vector<double> masses;
};

struct WindowData {
    long long total_count = 0;
    double total_mass = 0.0;
    std::array<long long, kComponentCount> component_counts{{0, 0, 0, 0}};
    std::array<double, kComponentCount> component_masses{{0.0, 0.0, 0.0, 0.0}};
};

struct AtomRecord {
    double z = 0.0;
    double mass = 0.0;
    std::size_t component = 0;
};

struct Options {
    std::string data_file;
    std::string info_file;
    std::string output_file;
    double requested_sample_spacing = 1.0;
    double requested_window_width = 1.0;
};

std::string trim(const std::string &input) {
    const std::string whitespace = " \t\r\n";
    const std::size_t first = input.find_first_not_of(whitespace);
    if (first == std::string::npos) {
        return "";
    }
    const std::size_t last = input.find_last_not_of(whitespace);
    return input.substr(first, last - first + 1);
}

bool begins_with(const std::string &text, const std::string &prefix) {
    return text.size() >= prefix.size() &&
           text.compare(0, prefix.size(), prefix) == 0;
}

std::string read_text_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open info file: " + path);
    }
    std::ostringstream buffer;
    buffer << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading info file: " + path);
    }
    return buffer.str();
}

std::size_t matching_brace(const std::string &text, std::size_t open_position) {
    if (open_position >= text.size() || text[open_position] != '{') {
        throw std::runtime_error("internal JSON object parsing error");
    }

    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t i = open_position; i < text.size(); ++i) {
        const char c = text[i];
        if (in_string) {
            if (escaped) {
                escaped = false;
            } else if (c == '\\') {
                escaped = true;
            } else if (c == '"') {
                in_string = false;
            }
            continue;
        }

        if (c == '"') {
            in_string = true;
        } else if (c == '{') {
            ++depth;
        } else if (c == '}') {
            --depth;
            if (depth == 0) {
                return i;
            }
        }
    }
    throw std::runtime_error("unterminated JSON object");
}

std::string json_object_for_key(const std::string &text, const std::string &key) {
    const std::string quoted_key = "\"" + key + "\"";
    const std::size_t key_position = text.find(quoted_key);
    if (key_position == std::string::npos) {
        throw std::runtime_error("missing JSON object: " + key);
    }
    const std::size_t colon = text.find(':', key_position + quoted_key.size());
    const std::size_t open = text.find('{', colon);
    if (colon == std::string::npos || open == std::string::npos) {
        throw std::runtime_error("invalid JSON object: " + key);
    }
    const std::size_t close = matching_brace(text, open);
    return text.substr(open, close - open + 1);
}

std::vector<std::pair<std::string, std::string>>
json_child_objects(const std::string &object_text) {
    std::vector<std::pair<std::string, std::string>> children;
    if (object_text.size() < 2 || object_text.front() != '{') {
        throw std::runtime_error("invalid components JSON object");
    }

    std::size_t i = 1;
    while (i + 1 < object_text.size()) {
        while (i < object_text.size() &&
               (std::isspace(static_cast<unsigned char>(object_text[i])) ||
                object_text[i] == ',')) {
            ++i;
        }
        if (i >= object_text.size() || object_text[i] == '}') {
            break;
        }
        if (object_text[i] != '"') {
            throw std::runtime_error("invalid key in components JSON object");
        }

        const std::size_t key_begin = ++i;
        bool escaped = false;
        while (i < object_text.size()) {
            if (!escaped && object_text[i] == '"') {
                break;
            }
            if (!escaped && object_text[i] == '\\') {
                escaped = true;
            } else {
                escaped = false;
            }
            ++i;
        }
        if (i >= object_text.size()) {
            throw std::runtime_error("unterminated key in components JSON object");
        }
        const std::string key = object_text.substr(key_begin, i - key_begin);
        ++i;

        while (i < object_text.size() &&
               std::isspace(static_cast<unsigned char>(object_text[i]))) {
            ++i;
        }
        if (i >= object_text.size() || object_text[i] != ':') {
            throw std::runtime_error("missing colon in components JSON object");
        }
        ++i;
        while (i < object_text.size() &&
               std::isspace(static_cast<unsigned char>(object_text[i]))) {
            ++i;
        }
        if (i >= object_text.size() || object_text[i] != '{') {
            throw std::runtime_error("component entry is not a JSON object: " + key);
        }

        const std::size_t close = matching_brace(object_text, i);
        children.emplace_back(key, object_text.substr(i, close - i + 1));
        i = close + 1;
    }
    return children;
}

long long json_integer(const std::string &object_text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*([0-9]+)");
    std::smatch match;
    if (!std::regex_search(object_text, match, pattern)) {
        throw std::runtime_error("missing integer field in model info: " + key);
    }
    return std::stoll(match[1].str());
}

std::string json_string(const std::string &text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error("missing string field in model info: " + key);
    }
    return match[1].str();
}

ComponentInfo make_component(const std::pair<std::string, std::string> &entry) {
    ComponentInfo component;
    component.key = entry.first;
    component.molecules = json_integer(entry.second, "M");
    component.expected_beads = json_integer(entry.second, "beads");
    return component;
}

ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    ModelInfo model;
    model.case_name = json_string(text, "case_name");
    model.geometry = json_string(text, "geometry");
    if (model.geometry != "bulk" && model.geometry != "film") {
        throw std::runtime_error(
            "unsupported geometry in model info: " + model.geometry);
    }

    const std::string components_object = json_object_for_key(text, "components");
    const auto entries = json_child_objects(components_object);
    if (entries.size() != kComponentCount) {
        throw std::runtime_error(
            "the analyzer requires exactly four component entries in the info file");
    }

    bool have_network = false;
    bool have_crosslinker = false;
    bool have_moderator = false;
    std::vector<ComponentInfo> unassigned;

    for (const auto &entry : entries) {
        const ComponentInfo component = make_component(entry);
        if (entry.first == "network_strands") {
            model.components[kNetwork] = component;
            have_network = true;
        } else if (entry.first == "crosslinkers") {
            model.components[kCrosslinker] = component;
            have_crosslinker = true;
        } else if (entry.first == "star_moderators") {
            model.components[kModerator] = component;
            have_moderator = true;
        } else {
            unassigned.push_back(component);
        }
    }

    if (!have_network || !have_crosslinker || !have_moderator ||
        unassigned.size() != 1) {
        throw std::runtime_error(
            "could not identify network, crosslinker, filler, and moderator components");
    }
    model.components[kFiller] = unassigned.front();
    return model;
}

DataHeader parse_data_header(const std::string &path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open LAMMPS data file: " + path);
    }

    DataHeader header;
    bool in_masses = false;
    int masses_read = 0;
    std::string line;
    while (std::getline(input, line)) {
        const std::string clean = trim(line);
        if (clean.empty()) {
            continue;
        }

        if (begins_with(clean, "Masses")) {
            if (header.atom_types <= 0) {
                throw std::runtime_error(
                    "atom type count must appear before the Masses section");
            }
            header.masses.assign(static_cast<std::size_t>(header.atom_types + 1), 0.0);
            in_masses = true;
            continue;
        }

        if (in_masses && masses_read < header.atom_types) {
            std::istringstream mass_line(clean);
            int type = 0;
            double mass = 0.0;
            if (mass_line >> type >> mass) {
                if (type < 1 || type > header.atom_types || mass <= 0.0) {
                    throw std::runtime_error("invalid entry in Masses section");
                }
                if (header.masses[static_cast<std::size_t>(type)] != 0.0) {
                    throw std::runtime_error("duplicate atom type in Masses section");
                }
                header.masses[static_cast<std::size_t>(type)] = mass;
                ++masses_read;
                if (masses_read == header.atom_types) {
                    in_masses = false;
                }
                continue;
            }
        }

        std::istringstream fields(clean);
        long long integer_value = 0;
        std::string word1;
        std::string word2;
        if (fields >> integer_value >> word1) {
            if (word1 == "atoms") {
                header.atoms = integer_value;
                continue;
            }
            if (word1 == "atom" && fields >> word2 && word2 == "types") {
                header.atom_types = static_cast<int>(integer_value);
                continue;
            }
        }

        std::istringstream bounds(clean);
        double lo = 0.0;
        double hi = 0.0;
        std::string lower_label;
        std::string upper_label;
        if (bounds >> lo >> hi >> lower_label >> upper_label) {
            if (lower_label == "xlo" && upper_label == "xhi") {
                header.xlo = lo;
                header.xhi = hi;
                header.have_x = true;
            } else if (lower_label == "ylo" && upper_label == "yhi") {
                header.ylo = lo;
                header.yhi = hi;
                header.have_y = true;
            } else if (lower_label == "zlo" && upper_label == "zhi") {
                header.zlo = lo;
                header.zhi = hi;
                header.have_z = true;
            }
        }
    }

    if (header.atoms <= 0) {
        throw std::runtime_error("invalid or missing atom count in data file");
    }
    if (header.atom_types <= 0 ||
        masses_read != header.atom_types ||
        header.masses.size() != static_cast<std::size_t>(header.atom_types + 1)) {
        throw std::runtime_error("incomplete Masses section in data file");
    }
    if (!header.have_x || !header.have_y || !header.have_z ||
        header.xhi <= header.xlo || header.yhi <= header.ylo ||
        header.zhi <= header.zlo) {
        throw std::runtime_error("invalid or missing orthogonal box bounds");
    }
    return header;
}

std::size_t component_for_molecule(
    long long molecule,
    const std::array<long long, kComponentCount> &cumulative_molecules) {
    if (molecule < 1 || molecule > cumulative_molecules.back()) {
        throw std::runtime_error(
            "atom has molecule ID outside the four component ranges: " +
            std::to_string(molecule));
    }
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        if (molecule <= cumulative_molecules[component]) {
            return component;
        }
    }
    throw std::runtime_error("internal molecule-range classification error");
}

std::string sanitized_case_name(std::string name) {
    for (char &character : name) {
        const unsigned char c = static_cast<unsigned char>(character);
        if (!std::isalnum(c) && character != '-' && character != '_') {
            character = '_';
        }
    }
    if (name.empty()) {
        return "model";
    }
    return name;
}

double parse_positive_double(const std::string &text, const std::string &option) {
    std::size_t consumed = 0;
    double value = 0.0;
    try {
        value = std::stod(text, &consumed);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid value for " + option + ": " + text);
    }
    if (consumed != text.size() || !std::isfinite(value) || value <= 0.0) {
        throw std::runtime_error("value for " + option + " must be positive");
    }
    return value;
}

void print_help(const char *program) {
    std::cout
        << "Usage:\n"
        << "  " << program
        << " DATA_FILE INFO_FILE [--sample-spacing X] [--window-width X]\n"
        << "      [--output FILE]\n\n"
        << "Analyze one final LAMMPS atom_style full data file along z.\n"
        << "The four components are identified from molecule-ID ranges in the\n"
        << "V22/V35 JSON .info file. Output is a numeric-only table for MATLAB.\n\n"
        << "Options:\n"
        << "  --sample-spacing X  Distance between profile points in angstrom\n"
        << "                      (default: 1.0)\n"
        << "  --window-width X    Width of the counting window in angstrom\n"
        << "                      (default: 1.0)\n"
        << "  --bin-width X       Alias for --window-width\n"
        << "  --output FILE       Output path (default: z_profile.<case>.dat)\n"
        << "  --help              Show this command reference\n";
}

Options parse_options(int argc, char **argv) {
    Options options;
    std::vector<std::string> positional;

    for (int i = 1; i < argc; ++i) {
        const std::string argument = argv[i];
        if (argument == "--help" || argument == "-h") {
            print_help(argv[0]);
            std::exit(0);
        }
        if (argument == "--sample-spacing") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--sample-spacing requires a value");
            }
            options.requested_sample_spacing =
                parse_positive_double(argv[++i], "--sample-spacing");
        } else if (argument == "--window-width" || argument == "--bin-width") {
            if (i + 1 >= argc) {
                throw std::runtime_error(argument + " requires a value");
            }
            options.requested_window_width =
                parse_positive_double(argv[++i], argument);
        } else if (argument == "--output") {
            if (i + 1 >= argc) {
                throw std::runtime_error("--output requires a file path");
            }
            options.output_file = argv[++i];
            if (options.output_file.empty()) {
                throw std::runtime_error("--output cannot be empty");
            }
        } else if (begins_with(argument, "--")) {
            throw std::runtime_error("unknown option: " + argument);
        } else {
            positional.push_back(argument);
        }
    }

    if (positional.size() != 2) {
        throw std::runtime_error(
            "expected DATA_FILE and INFO_FILE; use --help for usage");
    }
    options.data_file = positional[0];
    options.info_file = positional[1];
    return options;
}

void analyze(const Options &input_options) {
    Options options = input_options;
    const ModelInfo model = parse_model_info(options.info_file);
    const DataHeader header = parse_data_header(options.data_file);

    if (options.output_file.empty()) {
        options.output_file =
            "z_profile." + sanitized_case_name(model.case_name) + ".dat";
    }

    std::array<long long, kComponentCount> cumulative_molecules{{0, 0, 0, 0}};
    long long molecule_total = 0;
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        molecule_total += model.components[component].molecules;
        cumulative_molecules[component] = molecule_total;
    }
    if (molecule_total <= 0) {
        throw std::runtime_error("model info contains no molecules");
    }

    const double lz = header.zhi - header.zlo;
    if (options.requested_window_width > lz) {
        throw std::runtime_error(
            "window width cannot exceed the final box thickness");
    }
    const long long rounded_sample_points = static_cast<long long>(
        std::llround(lz / options.requested_sample_spacing));
    const std::size_t number_of_sample_points =
        static_cast<std::size_t>(std::max(1LL, rounded_sample_points));
    const double sample_spacing =
        lz / static_cast<double>(number_of_sample_points);
    const long long rounded_window_points = static_cast<long long>(
        std::llround(options.requested_window_width / sample_spacing));
    const std::size_t window_points = static_cast<std::size_t>(
        std::max(1LL, std::min(
            static_cast<long long>(number_of_sample_points),
            rounded_window_points)));
    const double window_width =
        sample_spacing * static_cast<double>(window_points);

    std::ifstream data(options.data_file);
    if (!data) {
        throw std::runtime_error("cannot reopen LAMMPS data file: " + options.data_file);
    }

    bool found_atoms_section = false;
    long long atoms_read = 0;
    double system_mass = 0.0;
    std::array<long long, kComponentCount> observed_counts{{0, 0, 0, 0}};
    std::vector<AtomRecord> atoms;
    atoms.reserve(static_cast<std::size_t>(header.atoms));
    std::string line;

    while (std::getline(data, line)) {
        const std::string clean = trim(line);
        if (!found_atoms_section) {
            if (begins_with(clean, "Atoms")) {
                found_atoms_section = true;
            }
            continue;
        }
        if (atoms_read == header.atoms) {
            break;
        }
        if (clean.empty() || clean.front() == '#') {
            continue;
        }

        std::istringstream atom_line(clean);
        long long atom_id = 0;
        long long molecule = 0;
        int type = 0;
        double charge = 0.0;
        double x = 0.0;
        double y = 0.0;
        double z = 0.0;
        if (!(atom_line >> atom_id >> molecule >> type >> charge >> x >> y >> z)) {
            throw std::runtime_error(
                "could not parse atom record after reading " +
                std::to_string(atoms_read) + " atoms");
        }
        (void)atom_id;
        (void)charge;
        (void)x;
        (void)y;

        if (type < 1 || type > header.atom_types) {
            throw std::runtime_error("atom type is outside the Masses table");
        }
        const double mass = header.masses[static_cast<std::size_t>(type)];
        const std::size_t component =
            component_for_molecule(molecule, cumulative_molecules);

        const double tolerance = 1.0e-8 * std::max(1.0, lz);
        if (z < header.zlo - tolerance || z > header.zhi + tolerance) {
            throw std::runtime_error(
                "atom z coordinate is outside the primary box: " + std::to_string(z));
        }
        const double clamped_z = std::min(header.zhi, std::max(header.zlo, z));
        atoms.push_back(AtomRecord{clamped_z, mass, component});
        ++observed_counts[component];
        system_mass += mass;
        ++atoms_read;
    }

    if (!found_atoms_section) {
        throw std::runtime_error("missing Atoms section in data file");
    }
    if (atoms_read != header.atoms) {
        throw std::runtime_error(
            "Atoms section ended after " + std::to_string(atoms_read) +
            " records; header declares " + std::to_string(header.atoms));
    }
    if (system_mass <= 0.0) {
        throw std::runtime_error("calculated system mass is not positive");
    }
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        if (observed_counts[component] != model.components[component].expected_beads) {
            throw std::runtime_error(
                "component bead count does not match info file for " +
                model.components[component].key + ": observed " +
                std::to_string(observed_counts[component]) + ", expected " +
                std::to_string(model.components[component].expected_beads));
        }
    }

    std::vector<WindowData> windows(number_of_sample_points);
    std::vector<double> window_lower(number_of_sample_points, 0.0);
    std::vector<double> window_upper(number_of_sample_points, 0.0);
    std::vector<double> effective_window_width(number_of_sample_points, 0.0);
    const double half_window = 0.5 * window_width;
    const bool periodic_z = model.geometry == "bulk";

    for (std::size_t point = 0; point < number_of_sample_points; ++point) {
        const double center =
            header.zlo + (static_cast<double>(point) + 0.5) * sample_spacing;
        if (periodic_z) {
            window_lower[point] = center - half_window;
            window_upper[point] = center + half_window;
            effective_window_width[point] = window_width;
        } else {
            window_lower[point] = std::max(header.zlo, center - half_window);
            window_upper[point] = std::min(header.zhi, center + half_window);
            effective_window_width[point] =
                window_upper[point] - window_lower[point];
        }

        WindowData &values = windows[point];
        for (const AtomRecord &atom : atoms) {
            bool included = false;
            if (periodic_z) {
                double distance = std::abs(atom.z - center);
                distance = std::min(distance, lz - distance);
                included = distance < half_window;
            } else {
                included =
                    atom.z >= window_lower[point] && atom.z < window_upper[point];
                if (atom.z == header.zhi &&
                    window_upper[point] == header.zhi) {
                    included = true;
                }
            }
            if (!included) {
                continue;
            }

            ++values.total_count;
            values.total_mass += atom.mass;
            ++values.component_counts[atom.component];
            values.component_masses[atom.component] += atom.mass;
        }
    }

    std::ofstream output(options.output_file);
    if (!output) {
        throw std::runtime_error("cannot create output file: " + options.output_file);
    }
    output << std::fixed << std::setprecision(12);

    for (std::size_t point = 0; point < number_of_sample_points; ++point) {
        const double z_center =
            header.zlo + (static_cast<double>(point) + 0.5) * sample_spacing;
        const WindowData &values = windows[point];
        const double relative_density =
            values.total_mass * lz /
            (system_mass * effective_window_width[point]);

        output << (point + 1) << ' '
               << window_lower[point] << ' '
               << z_center << ' '
               << window_upper[point] << ' '
               << values.total_count << ' '
               << relative_density;

        for (const long long count : values.component_counts) {
            output << ' ' << count;
        }

        if (values.total_mass > 0.0) {
            for (const double component_mass : values.component_masses) {
                output << ' ' << (component_mass / values.total_mass);
            }
        } else {
            for (std::size_t component = 0; component < kComponentCount; ++component) {
                output << " NaN";
            }
        }
        output << '\n';
    }
    if (!output) {
        throw std::runtime_error("failed while writing output file");
    }

    std::cout << "Analyzed case: " << model.case_name << '\n'
              << "Atoms: " << atoms_read << '\n'
              << "Geometry: " << model.geometry << '\n'
              << "Box z range: [" << header.zlo << ", " << header.zhi << "] A\n"
              << "Requested sample spacing: "
              << options.requested_sample_spacing << " A\n"
              << "Number of sample points: " << number_of_sample_points << '\n'
              << "Actual sample spacing: " << sample_spacing << " A\n"
              << "Requested window width: "
              << options.requested_window_width << " A\n"
              << "Window spans: " << window_points << " sample spacing(s)\n"
              << "Actual window width: " << window_width << " A\n"
              << "Component molecule ranges:\n"
              << "  network:     1-" << cumulative_molecules[kNetwork] << '\n'
              << "  crosslinker: " << (cumulative_molecules[kNetwork] + 1)
              << '-' << cumulative_molecules[kCrosslinker] << '\n'
              << "  filler:      " << (cumulative_molecules[kCrosslinker] + 1)
              << '-' << cumulative_molecules[kFiller] << '\n'
              << "  moderator:   " << (cumulative_molecules[kFiller] + 1)
              << '-' << cumulative_molecules[kModerator] << '\n'
              << "Output: " << options.output_file << '\n';
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        analyze(options);
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
