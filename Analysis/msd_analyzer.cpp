#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

enum ComponentIndex {
    kNetwork = 0,
    kCrosslinker = 1,
    kFiller = 2,
    kModerator = 3,
    kComponentCount = 4
};

enum class Selection { kAll, kFiller };
enum class ParticleMode { kBeads, kMoleculeCom };
enum class AveragingMode { kRaw, kTimeAveraged };

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

Vec3 operator-(const Vec3 &left, const Vec3 &right) {
    return {left.x - right.x, left.y - right.y, left.z - right.z};
}

struct ComponentInfo {
    std::string key;
    long long molecules = 0;
    long long beads = 0;
};

struct ModelInfo {
    std::string format;
    long long format_version = 0;
    std::string formulation;
    std::string case_name;
    std::string geometry;
    std::array<ComponentInfo, kComponentCount> components;
    long long total_beads = 0;
    long long total_molecules = 0;
    double total_mass = 0.0;
    std::vector<double> atom_type_masses;
    double timestep_fs = 0.0;
    long long production_steps = 0;
    double production_duration_ns = 0.0;
    long long dump_every_steps = 0;
    long long expected_frames = 0;
};

struct Options {
    std::string trajectory_file;
    std::string info_file;
    Selection selection = Selection::kFiller;
    ParticleMode particle_mode = ParticleMode::kMoleculeCom;
    AveragingMode averaging_mode = AveragingMode::kTimeAveraged;
    bool remove_system_com = true;
    long long origin_stride = 1;
    long long max_lag_frames = -1;
    double fit_start_ns = std::numeric_limits<double>::quiet_NaN();
    double fit_end_ns = std::numeric_limits<double>::quiet_NaN();
    std::string output_file;
    std::string report_file;
};

struct Box {
    std::array<double, 3> lo{{0.0, 0.0, 0.0}};
    std::array<double, 3> hi{{0.0, 0.0, 0.0}};
    std::array<bool, 3> periodic{{false, false, false}};

    double length(std::size_t dimension) const {
        return hi[dimension] - lo[dimension];
    }
};

struct TrajectoryData {
    std::vector<long long> timesteps;
    std::vector<Vec3> system_com;
    std::vector<Vec3> stored_positions;
    std::vector<Vec3> reference_positions;
    std::vector<struct MsdPoint> raw_points;
    long long entity_count = 0;
    long long selected_bead_count = 0;
    long long selected_molecule_count = 0;
    long long frame_atom_count = 0;
    Box box;
    double system_mass = 0.0;
};

struct MsdPoint {
    long long lag_index = 0;
    long long lag_steps = 0;
    double time_ps = 0.0;
    long long origins = 0;
    long long entities = 0;
    double msd_x = 0.0;
    double msd_y = 0.0;
    double msd_z = 0.0;
    double msd_xy = 0.0;
    double msd_3d = 0.0;
};

struct LinearFit {
    bool valid = false;
    long long points = 0;
    double slope = 0.0;
    double intercept = 0.0;
    double r_squared = std::numeric_limits<double>::quiet_NaN();
};

struct DiffusionFits {
    double requested_start_ns = 0.0;
    double requested_end_ns = 0.0;
    double actual_start_ns = 0.0;
    double actual_end_ns = 0.0;
    bool automatic_range = false;
    LinearFit x;
    LinearFit y;
    LinearFit z;
    LinearFit xy;
    LinearFit three_d;
};

std::string trim(const std::string &input) {
    const std::string whitespace = " \t\r\n";
    const std::size_t first = input.find_first_not_of(whitespace);
    if (first == std::string::npos) return "";
    const std::size_t last = input.find_last_not_of(whitespace);
    return input.substr(first, last - first + 1);
}

std::vector<std::string> split_words(const std::string &line) {
    std::istringstream input(line);
    std::vector<std::string> words;
    std::string word;
    while (input >> word) words.push_back(word);
    return words;
}

std::string read_text_file(const std::string &path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open info file: " + path);
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading info file: " + path);
    }
    return contents.str();
}

std::size_t matching_brace(const std::string &text, std::size_t open) {
    if (open >= text.size() || text[open] != '{') {
        throw std::runtime_error("internal JSON parsing error");
    }
    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t index = open; index < text.size(); ++index) {
        const char value = text[index];
        if (in_string) {
            if (escaped) escaped = false;
            else if (value == '\\') escaped = true;
            else if (value == '"') in_string = false;
            continue;
        }
        if (value == '"') in_string = true;
        else if (value == '{') ++depth;
        else if (value == '}' && --depth == 0) return index;
    }
    throw std::runtime_error("unterminated JSON object");
}

std::string json_object_for_key(const std::string &text, const std::string &key) {
    const std::string quoted = "\"" + key + "\"";
    const std::size_t key_position = text.find(quoted);
    if (key_position == std::string::npos) {
        throw std::runtime_error("missing JSON object: " + key);
    }
    const std::size_t colon = text.find(':', key_position + quoted.size());
    const std::size_t open = colon == std::string::npos
        ? std::string::npos : text.find('{', colon + 1);
    if (open == std::string::npos) {
        throw std::runtime_error("JSON field is not an object: " + key);
    }
    const std::size_t close = matching_brace(text, open);
    return text.substr(open, close - open + 1);
}

std::vector<std::pair<std::string, std::string>>
json_child_objects(const std::string &object_text) {
    std::vector<std::pair<std::string, std::string>> children;
    std::size_t position = 1;
    while (position + 1 < object_text.size()) {
        while (position < object_text.size() &&
               (std::isspace(static_cast<unsigned char>(object_text[position])) ||
                object_text[position] == ',')) {
            ++position;
        }
        if (position >= object_text.size() || object_text[position] == '}') break;
        if (object_text[position] != '"') {
            throw std::runtime_error("invalid key in JSON object");
        }
        const std::size_t key_begin = ++position;
        while (position < object_text.size() && object_text[position] != '"') {
            if (object_text[position] == '\\') ++position;
            ++position;
        }
        if (position >= object_text.size()) {
            throw std::runtime_error("unterminated JSON key");
        }
        const std::string key = object_text.substr(key_begin, position - key_begin);
        const std::size_t colon = object_text.find(':', ++position);
        const std::size_t open = colon == std::string::npos
            ? std::string::npos
            : object_text.find_first_not_of(" \t\r\n", colon + 1);
        if (open == std::string::npos || object_text[open] != '{') {
            throw std::runtime_error("JSON child is not an object: " + key);
        }
        const std::size_t close = matching_brace(object_text, open);
        children.emplace_back(key, object_text.substr(open, close - open + 1));
        position = close + 1;
    }
    return children;
}

double json_number(const std::string &text, const std::string &key) {
    const std::regex pattern(
        "\"" + key +
        "\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error("missing numeric info field: " + key);
    }
    return std::stod(match[1].str());
}

long long json_integer(const std::string &text, const std::string &key) {
    const double value = json_number(text, key);
    if (!std::isfinite(value) || value < 0.0 ||
        std::fabs(value - std::round(value)) > 1.0e-8) {
        throw std::runtime_error("info field is not a nonnegative integer: " + key);
    }
    return static_cast<long long>(std::llround(value));
}

std::string json_string(const std::string &text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(text, match, pattern)) {
        throw std::runtime_error("missing string info field: " + key);
    }
    return match[1].str();
}

ComponentInfo component_from_entry(
    const std::pair<std::string, std::string> &entry) {
    ComponentInfo component;
    component.key = entry.first;
    component.molecules = json_integer(entry.second, "M");
    component.beads = json_integer(entry.second, "beads");
    return component;
}

ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    ModelInfo info;
    info.format = json_string(text, "format");
    info.format_version = json_integer(text, "format_version");
    info.formulation = json_string(text, "formulation");
    info.case_name = json_string(text, "case_name");
    info.geometry = json_string(text, "geometry");
    if (info.format_version != 2 ||
        (info.format != "V22-model-info" && info.format != "V35-model-info")) {
        throw std::runtime_error(
            "MSD analysis requires a current V22/V35 model-info version-2 file");
    }
    if (info.geometry != "bulk" && info.geometry != "film") {
        throw std::runtime_error("unsupported geometry: " + info.geometry);
    }

    const auto component_entries =
        json_child_objects(json_object_for_key(text, "components"));
    bool have_network = false;
    bool have_crosslinker = false;
    bool have_filler = false;
    bool have_moderator = false;
    for (const auto &entry : component_entries) {
        if (entry.first == "network_strands") {
            info.components[kNetwork] = component_from_entry(entry);
            have_network = true;
        } else if (entry.first == "crosslinkers") {
            info.components[kCrosslinker] = component_from_entry(entry);
            have_crosslinker = true;
        } else if (entry.first == "star_moderators") {
            info.components[kModerator] = component_from_entry(entry);
            have_moderator = true;
        } else {
            if (have_filler) {
                throw std::runtime_error("more than one filler component in info file");
            }
            info.components[kFiller] = component_from_entry(entry);
            have_filler = true;
        }
    }
    if (!have_network || !have_crosslinker || !have_filler || !have_moderator) {
        throw std::runtime_error("could not identify all four formulation components");
    }

    const std::string composition = json_object_for_key(text, "composition");
    info.total_beads = json_integer(composition, "total_beads");
    info.total_molecules = json_integer(composition, "total_molecules");
    info.total_mass = json_number(composition, "total_mass_g_per_mol_equivalent");

    long long molecule_sum = 0;
    long long bead_sum = 0;
    for (const ComponentInfo &component : info.components) {
        molecule_sum += component.molecules;
        bead_sum += component.beads;
    }
    if (molecule_sum != info.total_molecules || bead_sum != info.total_beads) {
        throw std::runtime_error("component totals do not match composition totals");
    }

    const std::string force_field = json_object_for_key(text, "force_field");
    const std::string atom_types = json_object_for_key(force_field, "atom_types");
    const auto atom_type_entries = json_child_objects(atom_types);
    int maximum_type = 0;
    for (const auto &entry : atom_type_entries) {
        maximum_type = std::max(maximum_type, std::stoi(entry.first));
    }
    info.atom_type_masses.assign(static_cast<std::size_t>(maximum_type + 1), 0.0);
    for (const auto &entry : atom_type_entries) {
        const int type = std::stoi(entry.first);
        info.atom_type_masses[static_cast<std::size_t>(type)] =
            json_number(entry.second, "mass");
    }

    const std::string simulation = json_object_for_key(text, "simulation_template");
    info.timestep_fs = json_number(simulation, "timestep_fs");
    const std::string production = json_object_for_key(simulation, "msd_production");
    info.production_steps = json_integer(production, "steps");
    info.production_duration_ns = json_number(production, "duration_ns");
    info.dump_every_steps = json_integer(production, "dump_every_steps");
    info.expected_frames = json_integer(production, "expected_frames");
    if (info.timestep_fs <= 0.0 || info.dump_every_steps <= 0) {
        throw std::runtime_error("invalid MSD timing fields in info file");
    }
    return info;
}

std::array<long long, kComponentCount> component_ends(const ModelInfo &info) {
    std::array<long long, kComponentCount> ends{{0, 0, 0, 0}};
    long long total = 0;
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        total += info.components[component].molecules;
        ends[component] = total;
    }
    return ends;
}

int component_for_molecule(
    long long molecule,
    const std::array<long long, kComponentCount> &ends) {
    if (molecule < 1 || molecule > ends.back()) {
        throw std::runtime_error(
            "trajectory molecule ID outside info ranges: " +
            std::to_string(molecule));
    }
    for (std::size_t component = 0; component < kComponentCount; ++component) {
        if (molecule <= ends[component]) return static_cast<int>(component);
    }
    throw std::runtime_error("internal molecule classification error");
}

long long parse_nonnegative_integer(
    const std::string &value, const std::string &option, bool allow_zero) {
    std::size_t used = 0;
    long long result = 0;
    try {
        result = std::stoll(value, &used);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid integer for " + option + ": " + value);
    }
    if (used != value.size() || result < (allow_zero ? 0 : 1)) {
        throw std::runtime_error("invalid integer for " + option + ": " + value);
    }
    return result;
}

double parse_nonnegative_double(
    const std::string &value, const std::string &option) {
    std::size_t used = 0;
    double result = 0.0;
    try {
        result = std::stod(value, &used);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid number for " + option + ": " + value);
    }
    if (used != value.size() || !std::isfinite(result) || result < 0.0) {
        throw std::runtime_error("invalid number for " + option + ": " + value);
    }
    return result;
}

bool parse_yes_no(const std::string &value, const std::string &option) {
    if (value == "yes") return true;
    if (value == "no") return false;
    throw std::runtime_error(option + " must be yes or no");
}

std::string selection_name(Selection selection) {
    return selection == Selection::kAll ? "all" : "filler";
}

std::string particle_name(ParticleMode mode) {
    return mode == ParticleMode::kBeads ? "beads" : "molecule_com";
}

std::string averaging_name(AveragingMode mode) {
    return mode == AveragingMode::kRaw ? "raw" : "time_averaged";
}

void print_help(const char *program) {
    std::cout
        << "Usage:\n  " << program << " TRAJECTORY INFO_FILE [options]\n\n"
        << "Calculate bead or molecular-COM MSD from the generated MSD trajectory.\n\n"
        << "Options:\n"
        << "  --selection MODE       filler (default) or all\n"
        << "  --particle MODE        molecule-com (default) or beads\n"
        << "  --averaging MODE       time-averaged (default) or raw\n"
        << "  --remove-system-com X  yes (default) or no\n"
        << "  --origin-stride N      time-origin stride (default: 1)\n"
        << "  --max-lag-frames N     largest time-averaged lag; default: all\n"
        << "  --fit-start-ns X       diffusion-fit lower time\n"
        << "  --fit-end-ns X         diffusion-fit upper time\n"
        << "  --output FILE          numeric MSD table path\n"
        << "  --report-output FILE   text report path\n"
        << "  --help                 show this reference\n\n"
        << "The default is the recommended oil self-diffusion observable: exact\n"
        << "time-averaged, mass-weighted filler-chain COM MSD with whole-system\n"
        << "COM drift removed. If no fit range is given, the last half of the\n"
        << "available lag-time range is used as a preliminary fit.\n";
}

Options parse_options(int argc, char **argv) {
    Options options;
    std::vector<std::string> positional;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--help" || argument == "-h") {
            print_help(argv[0]);
            std::exit(0);
        }
        if (argument.rfind("--", 0) == 0) {
            if (index + 1 >= argc) {
                throw std::runtime_error("missing value after " + argument);
            }
            const std::string value = argv[++index];
            if (argument == "--selection") {
                if (value == "all") options.selection = Selection::kAll;
                else if (value == "filler") options.selection = Selection::kFiller;
                else throw std::runtime_error("--selection must be all or filler");
            } else if (argument == "--particle") {
                if (value == "beads") options.particle_mode = ParticleMode::kBeads;
                else if (value == "molecule-com") {
                    options.particle_mode = ParticleMode::kMoleculeCom;
                } else {
                    throw std::runtime_error(
                        "--particle must be beads or molecule-com");
                }
            } else if (argument == "--averaging") {
                if (value == "raw") options.averaging_mode = AveragingMode::kRaw;
                else if (value == "time-averaged") {
                    options.averaging_mode = AveragingMode::kTimeAveraged;
                } else {
                    throw std::runtime_error(
                        "--averaging must be raw or time-averaged");
                }
            } else if (argument == "--remove-system-com") {
                options.remove_system_com = parse_yes_no(value, argument);
            } else if (argument == "--origin-stride") {
                options.origin_stride =
                    parse_nonnegative_integer(value, argument, false);
            } else if (argument == "--max-lag-frames") {
                options.max_lag_frames =
                    parse_nonnegative_integer(value, argument, true);
            } else if (argument == "--fit-start-ns") {
                options.fit_start_ns = parse_nonnegative_double(value, argument);
            } else if (argument == "--fit-end-ns") {
                options.fit_end_ns = parse_nonnegative_double(value, argument);
            } else if (argument == "--output") {
                options.output_file = value;
            } else if (argument == "--report-output") {
                options.report_file = value;
            } else {
                throw std::runtime_error("unknown option: " + argument);
            }
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() != 2) {
        throw std::runtime_error(
            "expected TRAJECTORY and INFO_FILE; use --help");
    }
    options.trajectory_file = positional[0];
    options.info_file = positional[1];
    if (std::isfinite(options.fit_start_ns) &&
        std::isfinite(options.fit_end_ns) &&
        options.fit_end_ns <= options.fit_start_ns) {
        throw std::runtime_error("--fit-end-ns must exceed --fit-start-ns");
    }
    return options;
}

std::string sanitized_case_name(std::string name) {
    for (char &character : name) {
        const unsigned char value = static_cast<unsigned char>(character);
        if (!std::isalnum(value) && character != '-' && character != '_') {
            character = '_';
        }
    }
    return name.empty() ? "model" : name;
}

void resolve_output_paths(Options &options, const ModelInfo &info) {
    if (!options.output_file.empty() && !options.report_file.empty()) return;
    const std::string sample = sanitized_case_name(info.case_name);
    const std::filesystem::path directory(sample);
    std::error_code error;
    std::filesystem::create_directories(directory, error);
    if (error) {
        throw std::runtime_error(
            "cannot create sample output directory " + directory.string() +
            ": " + error.message());
    }
    const std::string tag = selection_name(options.selection) + "_" +
        particle_name(options.particle_mode) + "_" +
        averaging_name(options.averaging_mode);
    if (options.output_file.empty()) {
        options.output_file =
            (directory / ("msd_" + tag + "." + sample + ".dat")).string();
    }
    if (options.report_file.empty()) {
        options.report_file =
            (directory / ("msd_report_" + tag + "." + sample + ".txt")).string();
    }
}

long long parse_integer_line(const std::string &line, const std::string &label) {
    const std::string clean = trim(line);
    std::size_t used = 0;
    long long value = 0;
    try {
        value = std::stoll(clean, &used);
    } catch (const std::exception &) {
        throw std::runtime_error("invalid " + label + ": " + line);
    }
    if (used != clean.size()) {
        throw std::runtime_error("invalid " + label + ": " + line);
    }
    return value;
}

std::array<double, 2> parse_bound_line(const std::string &line) {
    const char *position = line.c_str();
    char *end = nullptr;
    std::array<double, 2> values{{0.0, 0.0}};
    for (double &value : values) {
        value = std::strtod(position, &end);
        if (end == position || !std::isfinite(value)) {
            throw std::runtime_error("invalid BOX BOUNDS line: " + line);
        }
        position = end;
    }
    return values;
}

void parse_numeric_fields(
    const std::string &line,
    std::size_t count,
    std::array<double, 16> &values) {
    if (count > values.size()) {
        throw std::runtime_error("too many columns in trajectory atom record");
    }
    const char *position = line.c_str();
    char *end = nullptr;
    for (std::size_t index = 0; index < count; ++index) {
        values[index] = std::strtod(position, &end);
        if (end == position || !std::isfinite(values[index])) {
            throw std::runtime_error("invalid trajectory atom record: " + line);
        }
        position = end;
    }
    while (*position != '\0' &&
           std::isspace(static_cast<unsigned char>(*position))) {
        ++position;
    }
    if (*position != '\0') {
        throw std::runtime_error("extra values in trajectory atom record");
    }
}

long long exact_integer(double value, const std::string &field) {
    if (std::fabs(value - std::round(value)) > 1.0e-8) {
        throw std::runtime_error("noninteger trajectory field: " + field);
    }
    return static_cast<long long>(std::llround(value));
}

std::size_t required_column(
    const std::map<std::string, std::size_t> &columns,
    const std::string &name) {
    const auto found = columns.find(name);
    if (found == columns.end()) {
        throw std::runtime_error("trajectory is missing atom column: " + name);
    }
    return found->second;
}

Vec3 calculate_raw_msd_sum(
    const std::vector<Vec3> &current,
    const std::vector<Vec3> &reference,
    const Vec3 &drift) {
    if (current.size() != reference.size()) {
        throw std::runtime_error("internal raw-MSD entity count mismatch");
    }
    Vec3 sum;
    for (std::size_t index = 0; index < current.size(); ++index) {
        const Vec3 displacement{
            current[index].x - reference[index].x - drift.x,
            current[index].y - reference[index].y - drift.y,
            current[index].z - reference[index].z - drift.z
        };
        sum.x += displacement.x * displacement.x;
        sum.y += displacement.y * displacement.y;
        sum.z += displacement.z * displacement.z;
    }
    return sum;
}

MsdPoint make_msd_point(
    long long lag_index,
    long long lag_steps,
    double timestep_fs,
    long long origins,
    long long entities,
    const Vec3 &squared_sum) {
    if (origins <= 0 || entities <= 0) {
        throw std::runtime_error("MSD average has no observations");
    }
    const double denominator =
        static_cast<double>(origins) * static_cast<double>(entities);
    MsdPoint point;
    point.lag_index = lag_index;
    point.lag_steps = lag_steps;
    point.time_ps = static_cast<double>(lag_steps) * timestep_fs / 1000.0;
    point.origins = origins;
    point.entities = entities;
    point.msd_x = squared_sum.x / denominator;
    point.msd_y = squared_sum.y / denominator;
    point.msd_z = squared_sum.z / denominator;
    point.msd_xy = point.msd_x + point.msd_y;
    point.msd_3d = point.msd_xy + point.msd_z;
    return point;
}

TrajectoryData read_trajectory(
    const Options &options,
    const ModelInfo &info,
    std::vector<std::string> &warnings) {
    std::ifstream input(options.trajectory_file);
    if (!input) {
        throw std::runtime_error(
            "cannot open trajectory: " + options.trajectory_file);
    }

    const auto ends = component_ends(info);
    const long long filler_first_molecule = ends[kCrosslinker] + 1;
    const long long filler_last_molecule = ends[kFiller];
    if (options.selection == Selection::kFiller &&
        info.components[kFiller].molecules == 0) {
        throw std::runtime_error("filler selection requested for a no-oil system");
    }

    const long long expected_selected_beads =
        options.selection == Selection::kAll
            ? info.total_beads : info.components[kFiller].beads;
    const long long expected_selected_molecules =
        options.selection == Selection::kAll
            ? info.total_molecules : info.components[kFiller].molecules;
    const long long expected_entities =
        options.particle_mode == ParticleMode::kBeads
            ? expected_selected_beads : expected_selected_molecules;

    TrajectoryData trajectory;
    trajectory.entity_count = expected_entities;
    trajectory.selected_bead_count = expected_selected_beads;
    trajectory.selected_molecule_count = expected_selected_molecules;
    trajectory.frame_atom_count = info.total_beads;
    if (options.averaging_mode == AveragingMode::kTimeAveraged) {
        const long double reserve_bytes =
            static_cast<long double>(expected_entities) *
            static_cast<long double>(info.expected_frames) *
            static_cast<long double>(sizeof(Vec3));
        if (reserve_bytes <= 256.0L * 1024.0L * 1024.0L) {
            trajectory.stored_positions.reserve(
                static_cast<std::size_t>(expected_entities * info.expected_frames));
        }
    }

    std::vector<long long> id_to_entity(
        static_cast<std::size_t>(info.total_beads + 1), -1);
    std::vector<long long> molecule_to_entity(
        static_cast<std::size_t>(info.total_molecules + 1), -1);
    if (options.particle_mode == ParticleMode::kMoleculeCom) {
        long long entity = 0;
        for (long long molecule = 1; molecule <= info.total_molecules; ++molecule) {
            const bool selected = options.selection == Selection::kAll ||
                (molecule >= filler_first_molecule &&
                 molecule <= filler_last_molecule);
            if (selected) {
                molecule_to_entity[static_cast<std::size_t>(molecule)] = entity++;
            }
        }
        if (entity != expected_entities) {
            throw std::runtime_error("internal selected-molecule count mismatch");
        }
    }

    std::vector<double> reference_entity_masses;
    std::vector<int> atom_seen_stamp(
        static_cast<std::size_t>(info.total_beads + 1), -1);
    std::string line;
    long long frame_index = 0;
    Box reference_box;
    bool have_reference_box = false;

    while (std::getline(input, line)) {
        if (trim(line).empty()) continue;
        if (trim(line) != "ITEM: TIMESTEP") {
            throw std::runtime_error(
                "expected ITEM: TIMESTEP before frame " +
                std::to_string(frame_index + 1));
        }
        if (!std::getline(input, line)) {
            throw std::runtime_error("missing trajectory timestep value");
        }
        const long long timestep = parse_integer_line(line, "trajectory timestep");
        if (!trajectory.timesteps.empty() &&
            timestep <= trajectory.timesteps.back()) {
            throw std::runtime_error("trajectory timesteps are not increasing");
        }

        if (!std::getline(input, line) || trim(line) != "ITEM: NUMBER OF ATOMS" ||
            !std::getline(input, line)) {
            throw std::runtime_error("invalid NUMBER OF ATOMS block");
        }
        const long long atom_count = parse_integer_line(line, "atom count");
        if (atom_count != info.total_beads) {
            throw std::runtime_error(
                "trajectory atom count does not match info total_beads");
        }

        if (!std::getline(input, line)) {
            throw std::runtime_error("missing BOX BOUNDS header");
        }
        const std::vector<std::string> box_header = split_words(line);
        if (box_header.size() != 6 || box_header[0] != "ITEM:" ||
            box_header[1] != "BOX" || box_header[2] != "BOUNDS") {
            throw std::runtime_error(
                "only orthogonal ITEM: BOX BOUNDS trajectories are supported");
        }
        Box box;
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            box.periodic[dimension] = !box_header[dimension + 3].empty() &&
                box_header[dimension + 3][0] == 'p';
            if (!std::getline(input, line)) {
                throw std::runtime_error("incomplete BOX BOUNDS block");
            }
            const auto bounds = parse_bound_line(line);
            box.lo[dimension] = bounds[0];
            box.hi[dimension] = bounds[1];
            if (box.length(dimension) <= 0.0) {
                throw std::runtime_error("nonpositive trajectory box length");
            }
        }
        if (info.geometry == "bulk" &&
            (!box.periodic[0] || !box.periodic[1] || !box.periodic[2])) {
            throw std::runtime_error("bulk info file does not match dump boundaries");
        }
        if (info.geometry == "film" &&
            (!box.periodic[0] || !box.periodic[1] || box.periodic[2])) {
            throw std::runtime_error("film info file does not match dump boundaries");
        }
        if (!have_reference_box) {
            reference_box = box;
            trajectory.box = box;
            have_reference_box = true;
        } else {
            for (std::size_t dimension = 0; dimension < 3; ++dimension) {
                const double scale = std::max(1.0, reference_box.length(dimension));
                if (std::fabs(box.lo[dimension] - reference_box.lo[dimension]) >
                        1.0e-10 * scale ||
                    std::fabs(box.hi[dimension] - reference_box.hi[dimension]) >
                        1.0e-10 * scale ||
                    box.periodic[dimension] != reference_box.periodic[dimension]) {
                    throw std::runtime_error(
                        "MSD production box changed between frames; the current "
                        "unwrapping assumes the generated fixed-box NVT stage");
                }
            }
        }

        if (!std::getline(input, line)) {
            throw std::runtime_error("missing ATOMS header");
        }
        const std::vector<std::string> atom_header = split_words(line);
        if (atom_header.size() < 3 || atom_header[0] != "ITEM:" ||
            atom_header[1] != "ATOMS") {
            throw std::runtime_error("invalid ITEM: ATOMS header");
        }
        const std::size_t field_count = atom_header.size() - 2;
        std::map<std::string, std::size_t> columns;
        for (std::size_t index = 0; index < field_count; ++index) {
            columns[atom_header[index + 2]] = index;
        }
        const std::size_t id_column = required_column(columns, "id");
        const std::size_t molecule_column = required_column(columns, "mol");
        const std::size_t type_column = required_column(columns, "type");
        const std::size_t x_column = required_column(columns, "x");
        const std::size_t y_column = required_column(columns, "y");
        const std::size_t z_column = required_column(columns, "z");
        const std::size_t ix_column = required_column(columns, "ix");
        const std::size_t iy_column = required_column(columns, "iy");
        const std::size_t iz_column = required_column(columns, "iz");

        std::vector<Vec3> current_positions(
            static_cast<std::size_t>(expected_entities));
        std::vector<double> current_entity_masses;
        if (options.particle_mode == ParticleMode::kMoleculeCom) {
            current_entity_masses.assign(
                static_cast<std::size_t>(expected_entities), 0.0);
        }
        long long selected_beads_seen = 0;
        double system_mass = 0.0;
        Vec3 weighted_system_position;
        std::array<double, 16> fields{};

        for (long long atom_line = 0; atom_line < atom_count; ++atom_line) {
            if (!std::getline(input, line)) {
                throw std::runtime_error("trajectory ends inside an atom block");
            }
            parse_numeric_fields(line, field_count, fields);
            const long long atom_id = exact_integer(fields[id_column], "id");
            const long long molecule =
                exact_integer(fields[molecule_column], "mol");
            const long long type = exact_integer(fields[type_column], "type");
            const long long ix = exact_integer(fields[ix_column], "ix");
            const long long iy = exact_integer(fields[iy_column], "iy");
            const long long iz = exact_integer(fields[iz_column], "iz");
            if (atom_id < 1 || atom_id > info.total_beads ||
                atom_seen_stamp[static_cast<std::size_t>(atom_id)] == frame_index) {
                throw std::runtime_error("invalid or duplicate atom ID in frame");
            }
            atom_seen_stamp[static_cast<std::size_t>(atom_id)] =
                static_cast<int>(frame_index);
            if (molecule < 1 || molecule > info.total_molecules) {
                throw std::runtime_error("invalid molecule ID in trajectory");
            }
            if (type < 1 ||
                static_cast<std::size_t>(type) >= info.atom_type_masses.size() ||
                info.atom_type_masses[static_cast<std::size_t>(type)] <= 0.0) {
                throw std::runtime_error("invalid atom type in trajectory");
            }
            const double mass = info.atom_type_masses[static_cast<std::size_t>(type)];
            const Vec3 position{
                fields[x_column] +
                    (box.periodic[0] ? static_cast<double>(ix) * box.length(0) : 0.0),
                fields[y_column] +
                    (box.periodic[1] ? static_cast<double>(iy) * box.length(1) : 0.0),
                fields[z_column] +
                    (box.periodic[2] ? static_cast<double>(iz) * box.length(2) : 0.0)
            };
            system_mass += mass;
            weighted_system_position.x += mass * position.x;
            weighted_system_position.y += mass * position.y;
            weighted_system_position.z += mass * position.z;

            const int component = component_for_molecule(molecule, ends);
            const bool selected = options.selection == Selection::kAll ||
                component == kFiller;
            if (!selected) continue;
            ++selected_beads_seen;

            if (options.particle_mode == ParticleMode::kBeads) {
                long long entity = id_to_entity[static_cast<std::size_t>(atom_id)];
                if (frame_index == 0) {
                    entity = selected_beads_seen - 1;
                    id_to_entity[static_cast<std::size_t>(atom_id)] = entity;
                } else if (entity < 0) {
                    throw std::runtime_error(
                        "selected atom IDs changed between trajectory frames");
                }
                current_positions[static_cast<std::size_t>(entity)] = position;
            } else {
                const long long entity =
                    molecule_to_entity[static_cast<std::size_t>(molecule)];
                if (entity < 0) {
                    throw std::runtime_error("internal molecule selection error");
                }
                Vec3 &sum = current_positions[static_cast<std::size_t>(entity)];
                sum.x += mass * position.x;
                sum.y += mass * position.y;
                sum.z += mass * position.z;
                current_entity_masses[static_cast<std::size_t>(entity)] += mass;
            }
        }

        if (selected_beads_seen != expected_selected_beads) {
            throw std::runtime_error(
                "selected bead count does not match the model-info component range");
        }
        if (std::fabs(system_mass - info.total_mass) >
            1.0e-8 * std::max(1.0, info.total_mass)) {
            throw std::runtime_error(
                "trajectory mass calculated from atom types differs from info mass");
        }
        const Vec3 current_system_com{
            weighted_system_position.x / system_mass,
            weighted_system_position.y / system_mass,
            weighted_system_position.z / system_mass
        };
        if (options.particle_mode == ParticleMode::kMoleculeCom) {
            if (frame_index == 0) {
                reference_entity_masses = current_entity_masses;
            }
            for (std::size_t entity = 0; entity < current_positions.size(); ++entity) {
                const double mass = current_entity_masses[entity];
                if (mass <= 0.0 ||
                    (frame_index > 0 &&
                     std::fabs(mass - reference_entity_masses[entity]) >
                         1.0e-10 * std::max(1.0, reference_entity_masses[entity]))) {
                    throw std::runtime_error(
                        "selected molecule mass changed or molecule is incomplete");
                }
                current_positions[entity].x /= mass;
                current_positions[entity].y /= mass;
                current_positions[entity].z /= mass;
            }
        }

        trajectory.timesteps.push_back(timestep);
        trajectory.system_com.push_back(current_system_com);
        if (frame_index == 0) {
            trajectory.reference_positions = current_positions;
            trajectory.system_mass = system_mass;
        }

        if (options.averaging_mode == AveragingMode::kRaw) {
            const Vec3 drift = options.remove_system_com
                ? current_system_com - trajectory.system_com.front() : Vec3{};
            const Vec3 squared_sum = calculate_raw_msd_sum(
                current_positions, trajectory.reference_positions, drift);
            trajectory.raw_points.push_back(make_msd_point(
                frame_index,
                timestep - trajectory.timesteps.front(),
                info.timestep_fs,
                1,
                expected_entities,
                squared_sum));
        } else {
            try {
                trajectory.stored_positions.insert(
                    trajectory.stored_positions.end(),
                    current_positions.begin(), current_positions.end());
            } catch (const std::bad_alloc &) {
                throw std::runtime_error(
                    "insufficient memory for exact time-averaged MSD; use raw "
                    "averaging, molecule-com particles, or fewer frames");
            }
        }

        ++frame_index;
    }
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading trajectory");
    }
    if (trajectory.timesteps.size() < 2) {
        throw std::runtime_error("MSD analysis requires at least two frames");
    }
    if (trajectory.timesteps.size() !=
        static_cast<std::size_t>(info.expected_frames)) {
        warnings.push_back(
            "trajectory contains " + std::to_string(trajectory.timesteps.size()) +
            " frames, while the info file expects " +
            std::to_string(info.expected_frames) +
            "; diffusion from a short excerpt is not production-quality");
    }
    const long long interval = trajectory.timesteps[1] - trajectory.timesteps[0];
    if (interval <= 0) throw std::runtime_error("invalid trajectory frame interval");
    for (std::size_t frame = 2; frame < trajectory.timesteps.size(); ++frame) {
        if (trajectory.timesteps[frame] - trajectory.timesteps[frame - 1] != interval) {
            if (options.averaging_mode == AveragingMode::kTimeAveraged) {
                throw std::runtime_error(
                    "time-averaged MSD requires uniformly spaced frames");
            }
            warnings.push_back("raw MSD trajectory has nonuniform frame spacing");
            break;
        }
    }
    if (interval != info.dump_every_steps) {
        warnings.push_back(
            "observed frame interval differs from info dump_every_steps");
    }
    return trajectory;
}

std::vector<MsdPoint> calculate_time_averaged_msd(
    const Options &options,
    const ModelInfo &info,
    const TrajectoryData &trajectory) {
    const long long frames = static_cast<long long>(trajectory.timesteps.size());
    const long long interval = trajectory.timesteps[1] - trajectory.timesteps[0];
    const long long max_lag = options.max_lag_frames < 0
        ? frames - 1 : std::min(options.max_lag_frames, frames - 1);
    std::vector<MsdPoint> points(static_cast<std::size_t>(max_lag + 1));
    const long long entities = trajectory.entity_count;

    long double estimated_pairs = 0.0L;
    for (long long lag = 0; lag <= max_lag; ++lag) {
        const long long origins =
            (frames - lag + options.origin_stride - 1) / options.origin_stride;
        estimated_pairs +=
            static_cast<long double>(origins) * static_cast<long double>(entities);
    }
    std::cout << "Time-averaged displacement evaluations: "
              << std::fixed << std::setprecision(0)
              << estimated_pairs << "\n";
    if (estimated_pairs > 1.0e9L) {
        std::cout
            << "Warning: this exact calculation is expensive. "
            << "Consider --particle molecule-com, --origin-stride N, or "
            << "--max-lag-frames N.\n";
    }

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)
#endif
    for (long long lag = 0; lag <= max_lag; ++lag) {
        Vec3 squared_sum;
        long long origin_count = 0;
        for (long long origin = 0; origin + lag < frames;
             origin += options.origin_stride) {
            const long long later = origin + lag;
            const Vec3 drift = options.remove_system_com
                ? trajectory.system_com[static_cast<std::size_t>(later)] -
                    trajectory.system_com[static_cast<std::size_t>(origin)]
                : Vec3{};
            const std::size_t first_offset =
                static_cast<std::size_t>(origin * entities);
            const std::size_t later_offset =
                static_cast<std::size_t>(later * entities);
            for (long long entity = 0; entity < entities; ++entity) {
                const Vec3 &first = trajectory.stored_positions[
                    first_offset + static_cast<std::size_t>(entity)];
                const Vec3 &second = trajectory.stored_positions[
                    later_offset + static_cast<std::size_t>(entity)];
                const double dx = second.x - first.x - drift.x;
                const double dy = second.y - first.y - drift.y;
                const double dz = second.z - first.z - drift.z;
                squared_sum.x += dx * dx;
                squared_sum.y += dy * dy;
                squared_sum.z += dz * dz;
            }
            ++origin_count;
        }
        points[static_cast<std::size_t>(lag)] = make_msd_point(
            lag, lag * interval, info.timestep_fs,
            origin_count, entities, squared_sum);
    }
    return points;
}

LinearFit linear_fit(
    const std::vector<MsdPoint> &points,
    double start_ns,
    double end_ns,
    double MsdPoint::*member) {
    std::vector<std::pair<double, double>> selected;
    for (const MsdPoint &point : points) {
        const double time_ns = point.time_ps / 1000.0;
        if (point.time_ps > 0.0 && time_ns >= start_ns && time_ns <= end_ns) {
            selected.emplace_back(point.time_ps, point.*member);
        }
    }
    LinearFit fit;
    fit.points = static_cast<long long>(selected.size());
    if (selected.size() < 3) return fit;
    double mean_x = 0.0;
    double mean_y = 0.0;
    for (const auto &value : selected) {
        mean_x += value.first;
        mean_y += value.second;
    }
    mean_x /= static_cast<double>(selected.size());
    mean_y /= static_cast<double>(selected.size());
    double covariance = 0.0;
    double variance_x = 0.0;
    double total_y = 0.0;
    for (const auto &value : selected) {
        const double dx = value.first - mean_x;
        const double dy = value.second - mean_y;
        covariance += dx * dy;
        variance_x += dx * dx;
        total_y += dy * dy;
    }
    if (variance_x <= 0.0) return fit;
    fit.slope = covariance / variance_x;
    fit.intercept = mean_y - fit.slope * mean_x;
    double residual = 0.0;
    for (const auto &value : selected) {
        const double predicted = fit.intercept + fit.slope * value.first;
        const double difference = value.second - predicted;
        residual += difference * difference;
    }
    fit.r_squared = total_y > 0.0 ? 1.0 - residual / total_y
                                  : std::numeric_limits<double>::quiet_NaN();
    fit.valid = true;
    return fit;
}

DiffusionFits calculate_fits(
    const Options &options,
    const std::vector<MsdPoint> &points) {
    if (points.size() < 2) {
        throw std::runtime_error("not enough MSD points for diffusion fitting");
    }
    const double maximum_ns = points.back().time_ps / 1000.0;
    DiffusionFits fits;
    fits.automatic_range = !std::isfinite(options.fit_start_ns) &&
                           !std::isfinite(options.fit_end_ns);
    fits.requested_start_ns = std::isfinite(options.fit_start_ns)
        ? options.fit_start_ns : 0.5 * maximum_ns;
    fits.requested_end_ns = std::isfinite(options.fit_end_ns)
        ? options.fit_end_ns : maximum_ns;
    if (fits.requested_end_ns <= fits.requested_start_ns ||
        fits.requested_start_ns > maximum_ns) {
        throw std::runtime_error("diffusion fit range does not overlap MSD times");
    }
    bool found = false;
    for (const MsdPoint &point : points) {
        const double time_ns = point.time_ps / 1000.0;
        if (point.time_ps > 0.0 &&
            time_ns >= fits.requested_start_ns &&
            time_ns <= fits.requested_end_ns) {
            if (!found) {
                fits.actual_start_ns = time_ns;
                found = true;
            }
            fits.actual_end_ns = time_ns;
        }
    }
    fits.x = linear_fit(points, fits.requested_start_ns,
                        fits.requested_end_ns, &MsdPoint::msd_x);
    fits.y = linear_fit(points, fits.requested_start_ns,
                        fits.requested_end_ns, &MsdPoint::msd_y);
    fits.z = linear_fit(points, fits.requested_start_ns,
                        fits.requested_end_ns, &MsdPoint::msd_z);
    fits.xy = linear_fit(points, fits.requested_start_ns,
                         fits.requested_end_ns, &MsdPoint::msd_xy);
    fits.three_d = linear_fit(points, fits.requested_start_ns,
                              fits.requested_end_ns, &MsdPoint::msd_3d);
    return fits;
}

void write_numeric_output(
    const std::string &path,
    const std::vector<MsdPoint> &points) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open MSD output: " + path);
    output << std::scientific << std::setprecision(10);
    for (const MsdPoint &point : points) {
        const double time_ns = point.time_ps / 1000.0;
        const double d3 = point.time_ps > 0.0
            ? point.msd_3d / (6.0 * point.time_ps)
            : std::numeric_limits<double>::quiet_NaN();
        const double dxy = point.time_ps > 0.0
            ? point.msd_xy / (4.0 * point.time_ps)
            : std::numeric_limits<double>::quiet_NaN();
        const double dz = point.time_ps > 0.0
            ? point.msd_z / (2.0 * point.time_ps)
            : std::numeric_limits<double>::quiet_NaN();
        output
            << point.lag_index << " "
            << point.lag_steps << " "
            << point.time_ps << " "
            << time_ns << " "
            << point.origins << " "
            << point.entities << " "
            << point.msd_x << " "
            << point.msd_y << " "
            << point.msd_z << " "
            << point.msd_xy << " "
            << point.msd_3d << " "
            << d3 << " "
            << d3 * 1.0e-4 << " "
            << dxy << " "
            << dz << "\n";
    }
    if (!output) throw std::runtime_error("failed while writing MSD output");
}

void write_fit_line(
    std::ostream &output,
    const std::string &label,
    const LinearFit &fit,
    double divisor) {
    output << std::left << std::setw(28) << label << std::right << ": ";
    if (!fit.valid) {
        output << "unavailable (" << fit.points << " fit points)\n";
        return;
    }
    const double diffusion_a2_ps = fit.slope / divisor;
    output << std::scientific << std::setprecision(8)
           << "slope=" << fit.slope << " A^2/ps"
           << ", D=" << diffusion_a2_ps << " A^2/ps"
           << ", D=" << diffusion_a2_ps * 1.0e-4 << " cm^2/s"
           << ", R^2=" << fit.r_squared
           << ", points=" << fit.points << "\n"
           << std::fixed << std::setprecision(8);
}

void write_report(
    const std::string &path,
    const Options &options,
    const ModelInfo &info,
    const TrajectoryData &trajectory,
    const std::vector<MsdPoint> &points,
    const DiffusionFits &fits,
    const std::vector<std::string> &warnings) {
    std::ofstream output(path);
    if (!output) throw std::runtime_error("cannot open MSD report: " + path);
    const long long observed_interval =
        trajectory.timesteps[1] - trajectory.timesteps[0];
    const double observed_duration_ns =
        static_cast<double>(trajectory.timesteps.back() -
                            trajectory.timesteps.front()) *
        info.timestep_fs / 1.0e6;
    const double stored_memory_mib =
        static_cast<double>(trajectory.stored_positions.size() * sizeof(Vec3)) /
        (1024.0 * 1024.0);
    const auto ends = component_ends(info);
    output << std::fixed << std::setprecision(8);
    output << "MSD AND DIFFUSION ANALYSIS\n"
           << "==========================\n"
           << "Case                        : " << info.case_name << "\n"
           << "Formulation                 : " << info.formulation << "\n"
           << "Geometry                    : " << info.geometry << "\n"
           << "Trajectory                  : " << options.trajectory_file << "\n"
           << "Model info                  : " << options.info_file << "\n"
           << "Numeric output              : " << options.output_file << "\n\n"
           << "ANALYSIS DEFINITION\n"
           << "Selection                   : " << selection_name(options.selection) << "\n"
           << "Particle                    : " << particle_name(options.particle_mode) << "\n"
           << "Averaging                   : " << averaging_name(options.averaging_mode) << "\n"
           << "Remove whole-system COM     : "
           << (options.remove_system_com ? "yes" : "no") << "\n"
           << "Time-origin stride          : " << options.origin_stride << " frame(s)\n"
           << "Entities averaged           : " << trajectory.entity_count << "\n"
           << "Selected beads              : " << trajectory.selected_bead_count << "\n"
           << "Selected molecules          : " << trajectory.selected_molecule_count << "\n";
    output << "Filler molecule ID range    : ";
    if (info.components[kFiller].molecules > 0) {
        output << ends[kCrosslinker] + 1 << " to " << ends[kFiller] << "\n";
    } else {
        output << "none\n";
    }
    output << "Stored coordinate memory    : " << stored_memory_mib << " MiB\n\n"
           << "TRAJECTORY\n"
           << "Frames read                 : " << trajectory.timesteps.size()
           << " / " << info.expected_frames << " expected\n"
           << "First/last timestep         : " << trajectory.timesteps.front()
           << " / " << trajectory.timesteps.back() << "\n"
           << "Observed frame interval     : " << observed_interval << " steps\n"
           << "Timestep                    : " << info.timestep_fs << " fs\n"
           << "Observed lag-time span      : " << observed_duration_ns << " ns\n"
           << "Configured production       : " << info.production_duration_ns << " ns\n"
           << "Box lengths (A)             : "
           << trajectory.box.length(0) << " "
           << trajectory.box.length(1) << " "
           << trajectory.box.length(2) << "\n"
           << "MSD points                  : " << points.size() << "\n\n"
           << "DIFFUSION FIT\n"
           << "Fit range source            : "
           << (fits.automatic_range ? "automatic last half" : "command line/default bound")
           << "\n"
           << "Requested fit range (ns)    : "
           << fits.requested_start_ns << " to " << fits.requested_end_ns << "\n"
           << "Actual fit points span (ns) : "
           << fits.actual_start_ns << " to " << fits.actual_end_ns << "\n";
    write_fit_line(output, "x: D_x = slope/2", fits.x, 2.0);
    write_fit_line(output, "y: D_y = slope/2", fits.y, 2.0);
    write_fit_line(output, "z: D_z = slope/2", fits.z, 2.0);
    write_fit_line(output, "xy: D_xy = slope/4", fits.xy, 4.0);
    write_fit_line(output, "3D: D = slope/6", fits.three_d, 6.0);
    output << "\nINTERPRETATION\n"
           << "The fitted slope, not a single final MSD/(2*d*t) value, is the\n"
           << "recommended diffusion estimate. Inspect linearity and R^2, then rerun\n"
           << "with an explicit fit range inside the long-time diffusive regime.\n";
    if (options.particle_mode == ParticleMode::kBeads) {
        output << "Bead MSD includes internal molecular motion. Filler molecular-COM\n"
               << "MSD is preferred for silicone-oil translational diffusion.\n";
    }
    if (info.geometry == "film") {
        output << "For a confined film, xy diffusion is normally more meaningful than\n"
               << "a long-time 3D diffusion coefficient because z motion is bounded.\n";
    }
    output << "\nCHECKS AND WARNINGS\n";
    if (warnings.empty()) output << "PASS: no consistency warnings.\n";
    else {
        for (const std::string &warning : warnings) {
            output << "WARNING: " << warning << "\n";
        }
    }
    if (!output) throw std::runtime_error("failed while writing MSD report");
}

}  // namespace

int main(int argc, char **argv) {
    try {
        Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        resolve_output_paths(options, info);
        std::vector<std::string> warnings;
        const TrajectoryData trajectory =
            read_trajectory(options, info, warnings);
        std::vector<MsdPoint> points;
        if (options.averaging_mode == AveragingMode::kRaw) {
            points = trajectory.raw_points;
        } else {
            points = calculate_time_averaged_msd(options, info, trajectory);
        }
        const DiffusionFits fits = calculate_fits(options, points);
        if (!fits.three_d.valid) {
            warnings.push_back(
                "fewer than three points lie in the diffusion-fit range");
        } else {
            if (fits.three_d.slope <= 0.0) {
                warnings.push_back(
                    "the fitted 3D MSD slope is nonpositive; this interval is "
                    "not a valid diffusive regime");
            }
            if (std::isfinite(fits.three_d.r_squared) &&
                fits.three_d.r_squared < 0.90) {
                warnings.push_back(
                    "the 3D diffusion fit has R^2 below 0.90; inspect the MSD "
                    "and choose a better linear fit range");
            }
        }
        write_numeric_output(options.output_file, points);
        write_report(
            options.report_file, options, info, trajectory, points, fits, warnings);
        std::cout << "Wrote MSD data: " << options.output_file << "\n"
                  << "Wrote MSD report: " << options.report_file << "\n";
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
