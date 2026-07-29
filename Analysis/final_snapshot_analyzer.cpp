#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
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

constexpr double kAvogadro = 6.02214076e23;
constexpr double kAngstromCubedToCmCubed = 1.0e-24;
constexpr double kBoltzmannKcalPerMolK = 0.00198720425864083;
constexpr double kLammpsRealMvv2e = 2390.05736153349;

enum ComponentIndex {
    kNetwork = 0,
    kCrosslinker = 1,
    kFiller = 2,
    kModerator = 3,
    kComponentCount = 4
};

struct ComponentInfo {
    std::string key;
    std::string label;
    long long molecules = 0;
    long long beads = 0;
    double weight_percent = 0.0;
};

struct ModelInfo {
    std::string format;
    long long format_version = 0;
    std::string formulation;
    std::string case_name;
    std::string geometry;
    std::vector<ComponentInfo> components;
    long long total_beads = 0;
    long long total_molecules = 0;
    double total_mass_g_per_mol = 0.0;
    double requested_filler_weight_percent = 0.0;
    double realized_filler_weight_percent = 0.0;
    long long initial_bonds = 0;
    long long type2_sites_total = 0;
    long long type3_sites_total = 0;
    int crosslink_bond_type = 2;
    double initial_density = 0.0;
    double target_compressed_density = 0.0;
    double timestep_fs = 0.0;
    long long total_run_steps = 0;
    long long bond_creation_active_steps = 0;
    double final_temperature_k = 0.0;
};

struct AtomRecord {
    bool seen = false;
    long long molecule = 0;
    int type = 0;
    double mass = 0.0;
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    long long ix = 0;
    long long iy = 0;
    long long iz = 0;
    bool have_images = false;
};

struct MoleculeRecord {
    long long beads = 0;
    double mass = 0.0;
    long long inter_crosslink_degree = 0;
    std::vector<long long> crosslink_attachment_atoms;
};

struct Bounds {
    double lo = 0.0;
    double hi = 0.0;
    bool seen = false;
};

struct DataSummary {
    long long timestep = -1;
    long long declared_atoms = 0;
    long long declared_bonds = 0;
    long long declared_angles = 0;
    long long declared_dihedrals = 0;
    int atom_types = 0;
    int bond_types = 0;
    int angle_types = 0;
    int dihedral_types = 0;
    Bounds x;
    Bounds y;
    Bounds z;
    std::vector<double> masses;

    long long atoms_read = 0;
    long long velocities_read = 0;
    long long bonds_read = 0;
    long long angles_read = 0;
    long long dihedrals_read = 0;
    long long nonzero_image_atoms = 0;
    double total_mass_g_per_mol = 0.0;
    double charge_sum = 0.0;
    double x_min = std::numeric_limits<double>::infinity();
    double x_max = -std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity();
    double y_max = -std::numeric_limits<double>::infinity();
    double z_min = std::numeric_limits<double>::infinity();
    double z_max = -std::numeric_limits<double>::infinity();

    std::vector<long long> atom_type_counts;
    std::map<int, long long> bond_type_counts;
    std::map<int, long long> angle_type_counts;
    std::map<int, long long> dihedral_type_counts;
    std::vector<long long> component_beads;
    std::vector<double> component_masses;
    std::vector<AtomRecord> atoms;
    std::vector<MoleculeRecord> molecules;

    long long crosslink_bonds_inter = 0;
    long long crosslink_bonds_intra = 0;
    bool pair_coeffs_present = false;
    bool bond_coeffs_present = false;
    bool angle_coeffs_present = false;
    bool dihedral_coeffs_present = false;

    double velocity_mass = 0.0;
    double momentum_x = 0.0;
    double momentum_y = 0.0;
    double momentum_z = 0.0;
    double mass_velocity_squared = 0.0;
};

struct Options {
    std::string data_file;
    std::string info_file;
    std::string output_file;
    std::string strand_output_file;
};

class DisjointSet {
  public:
    explicit DisjointSet(std::size_t size) : parent_(size), rank_(size, 0) {
        for (std::size_t i = 0; i < size; ++i) {
            parent_[i] = static_cast<long long>(i);
        }
    }

    long long find(long long value) {
        const std::size_t index = checked_index(value);
        if (parent_[index] != value) {
            parent_[index] = find(parent_[index]);
        }
        return parent_[index];
    }

    void unite(long long first, long long second) {
        long long root_first = find(first);
        long long root_second = find(second);
        if (root_first == root_second) {
            return;
        }
        std::size_t first_index = checked_index(root_first);
        std::size_t second_index = checked_index(root_second);
        if (rank_[first_index] < rank_[second_index]) {
            std::swap(root_first, root_second);
            std::swap(first_index, second_index);
        }
        parent_[second_index] = root_first;
        if (rank_[first_index] == rank_[second_index]) {
            ++rank_[first_index];
        }
    }

  private:
    std::size_t checked_index(long long value) const {
        if (value < 0 || static_cast<std::size_t>(value) >= parent_.size()) {
            throw std::runtime_error(
                "molecule ID outside the range declared by the info file: " +
                std::to_string(value));
        }
        return static_cast<std::size_t>(value);
    }

    std::vector<long long> parent_;
    std::vector<unsigned char> rank_;
};

struct NetworkComponent {
    long long root = 0;
    long long precursor_molecules = 0;
    long long network_molecules = 0;
    long long crosslinker_molecules = 0;
    long long precursor_beads = 0;
    double precursor_mass = 0.0;
};

struct EndToEndRecord {
    long long molecule = 0;
    long long first_atom = 0;
    long long second_atom = 0;
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;
    double distance = 0.0;
    bool used_unwrapped_coordinates = false;
    bool in_largest_component = false;
};

struct DistanceStatistics {
    long long count = 0;
    double mean = 0.0;
    double mean_squared = 0.0;
    double rms = 0.0;
    double standard_deviation = 0.0;
    double minimum = 0.0;
    double first_quartile = 0.0;
    double median = 0.0;
    double third_quartile = 0.0;
    double maximum = 0.0;
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
    std::ostringstream contents;
    contents << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading info file: " + path);
    }
    return contents.str();
}

std::size_t matching_brace(const std::string &text, std::size_t open_position) {
    if (open_position >= text.size() || text[open_position] != '{') {
        throw std::runtime_error("internal JSON object parsing error");
    }
    int depth = 0;
    bool in_string = false;
    bool escaped = false;
    for (std::size_t i = open_position; i < text.size(); ++i) {
        const char character = text[i];
        if (in_string) {
            if (escaped) {
                escaped = false;
            } else if (character == '\\') {
                escaped = true;
            } else if (character == '"') {
                in_string = false;
            }
            continue;
        }
        if (character == '"') {
            in_string = true;
        } else if (character == '{') {
            ++depth;
        } else if (character == '}') {
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
    if (colon == std::string::npos) {
        throw std::runtime_error("invalid JSON field: " + key);
    }
    const std::size_t open = text.find('{', colon + 1);
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
        if (position >= object_text.size() || object_text[position] == '}') {
            break;
        }
        if (object_text[position] != '"') {
            throw std::runtime_error("invalid key in JSON object");
        }
        const std::size_t key_begin = ++position;
        while (position < object_text.size() && object_text[position] != '"') {
            if (object_text[position] == '\\') {
                ++position;
            }
            ++position;
        }
        if (position >= object_text.size()) {
            throw std::runtime_error("unterminated JSON key");
        }
        const std::string key =
            object_text.substr(key_begin, position - key_begin);
        ++position;
        const std::size_t colon = object_text.find(':', position);
        if (colon == std::string::npos) {
            throw std::runtime_error("missing colon after JSON key: " + key);
        }
        const std::size_t open = object_text.find_first_not_of(" \t\r\n", colon + 1);
        if (open == std::string::npos || object_text[open] != '{') {
            throw std::runtime_error("JSON child is not an object: " + key);
        }
        const std::size_t close = matching_brace(object_text, open);
        children.emplace_back(key, object_text.substr(open, close - open + 1));
        position = close + 1;
    }
    return children;
}

double json_number(const std::string &object_text, const std::string &key) {
    const std::regex pattern(
        "\"" + key +
        "\"\\s*:\\s*(-?(?:[0-9]+(?:\\.[0-9]*)?|\\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::smatch match;
    if (!std::regex_search(object_text, match, pattern)) {
        throw std::runtime_error("missing numeric field in info file: " + key);
    }
    return std::stod(match[1].str());
}

long long json_integer(const std::string &object_text, const std::string &key) {
    const double value = json_number(object_text, key);
    if (!std::isfinite(value) || value < 0.0 ||
        std::fabs(value - std::round(value)) > 1.0e-8) {
        throw std::runtime_error("info field is not a nonnegative integer: " + key);
    }
    return static_cast<long long>(std::llround(value));
}

std::string json_string(const std::string &object_text, const std::string &key) {
    const std::regex pattern("\"" + key + "\"\\s*:\\s*\"([^\"]*)\"");
    std::smatch match;
    if (!std::regex_search(object_text, match, pattern)) {
        throw std::runtime_error("missing string field in info file: " + key);
    }
    return match[1].str();
}

ComponentInfo component_from_entry(
    const std::pair<std::string, std::string> &entry,
    const std::string &label) {
    ComponentInfo component;
    component.key = entry.first;
    component.label = label;
    component.molecules = json_integer(entry.second, "M");
    component.beads = json_integer(entry.second, "beads");
    component.weight_percent = json_number(entry.second, "weight_percent");
    return component;
}

ModelInfo parse_model_info(const std::string &path) {
    const std::string text = read_text_file(path);
    ModelInfo info;
    info.format = json_string(text, "format");
    info.format_version = json_integer(text, "format_version");
    info.formulation = json_string(text, "formulation");
    if (info.format_version != 2 ||
        (info.format != "V22-model-info" &&
         info.format != "V35-model-info")) {
        throw std::runtime_error(
            "unsupported model-info format: " + info.format +
            " version " + std::to_string(info.format_version) +
            " (supported: V22/V35 model-info version 2)");
    }
    info.case_name = json_string(text, "case_name");
    info.geometry = json_string(text, "geometry");
    info.components.resize(kComponentCount);

    const auto component_entries =
        json_child_objects(json_object_for_key(text, "components"));
    bool have_network = false;
    bool have_crosslinker = false;
    bool have_filler = false;
    bool have_moderator = false;
    for (const auto &entry : component_entries) {
        if (entry.first == "network_strands") {
            info.components[kNetwork] =
                component_from_entry(entry, "network strands");
            have_network = true;
        } else if (entry.first == "crosslinkers") {
            info.components[kCrosslinker] =
                component_from_entry(entry, "crosslinkers");
            have_crosslinker = true;
        } else if (entry.first == "star_moderators") {
            info.components[kModerator] =
                component_from_entry(entry, "star moderators");
            have_moderator = true;
        } else {
            if (have_filler) {
                throw std::runtime_error(
                    "more than one unrecognized filler component in info file");
            }
            info.components[kFiller] =
                component_from_entry(entry, entry.first);
            have_filler = true;
        }
    }
    if (!have_network || !have_crosslinker || !have_filler || !have_moderator) {
        throw std::runtime_error(
            "could not identify the four expected components in the info file");
    }

    const std::string composition = json_object_for_key(text, "composition");
    info.total_beads = json_integer(composition, "total_beads");
    info.total_molecules = json_integer(composition, "total_molecules");
    info.total_mass_g_per_mol =
        json_number(composition, "total_mass_g_per_mol_equivalent");
    info.requested_filler_weight_percent =
        json_number(composition, "requested_filler_weight_percent");
    info.realized_filler_weight_percent =
        json_number(composition, "realized_filler_weight_percent");

    const std::string topology = json_object_for_key(text, "topology");
    info.initial_bonds = json_integer(topology, "bonds");

    const std::string crosslinker = json_object_for_key(text, "crosslinker");
    info.type2_sites_total =
        json_integer(crosslinker, "type2_sites_total_including_moderators");
    info.type3_sites_total = json_integer(crosslinker, "type3_sites_total");

    const std::string force_field = json_object_for_key(text, "force_field");
    const std::string bond_remap =
        json_object_for_key(force_field, "oil_bond_type_remap");
    info.crosslink_bond_type =
        static_cast<int>(json_integer(bond_remap, "formulation_crosslink"));

    const std::string initial_state = json_object_for_key(text, "initial_state");
    info.initial_density = json_number(initial_state, "density_g_per_cm3");

    const std::string simulation =
        json_object_for_key(text, "simulation_template");
    info.target_compressed_density =
        json_number(simulation, "target_compressed_density_g_per_cm3");
    info.timestep_fs = json_number(simulation, "timestep_fs");
    info.total_run_steps = json_integer(simulation, "total_run_steps");
    info.bond_creation_active_steps =
        json_integer(simulation, "bond_creation_active_steps");
    info.final_temperature_k =
        json_number(simulation, "final_temperature_K");

    long long molecule_sum = 0;
    long long bead_sum = 0;
    for (const ComponentInfo &component : info.components) {
        molecule_sum += component.molecules;
        bead_sum += component.beads;
    }
    if (molecule_sum != info.total_molecules || bead_sum != info.total_beads) {
        throw std::runtime_error(
            "component totals do not match composition totals in info file");
    }
    return info;
}

std::vector<long long> component_ends(const ModelInfo &info) {
    std::vector<long long> ends(kComponentCount, 0);
    long long total = 0;
    for (int component = 0; component < kComponentCount; ++component) {
        total += info.components[static_cast<std::size_t>(component)].molecules;
        ends[static_cast<std::size_t>(component)] = total;
    }
    return ends;
}

int component_for_molecule(
    long long molecule, const std::vector<long long> &ends) {
    if (molecule < 1 || molecule > ends.back()) {
        throw std::runtime_error(
            "molecule ID outside the ranges declared by the info file: " +
            std::to_string(molecule));
    }
    for (int component = 0; component < kComponentCount; ++component) {
        if (molecule <= ends[static_cast<std::size_t>(component)]) {
            return component;
        }
    }
    throw std::runtime_error("internal component classification error");
}

enum class Section {
    kHeader,
    kMasses,
    kAtoms,
    kVelocities,
    kBonds,
    kAngles,
    kDihedrals,
    kOther
};

bool section_line(
    const std::string &line,
    const std::string &name) {
    if (!begins_with(line, name)) {
        return false;
    }
    return line.size() == name.size() ||
           std::isspace(static_cast<unsigned char>(line[name.size()])) ||
           line[name.size()] == '#';
}

bool update_section(
    const std::string &line, Section &section, DataSummary &summary) {
    if (section_line(line, "Masses")) {
        section = Section::kMasses;
    } else if (section_line(line, "Pair Coeffs") ||
               section_line(line, "PairIJ Coeffs")) {
        summary.pair_coeffs_present = true;
        section = Section::kOther;
    } else if (section_line(line, "Bond Coeffs")) {
        summary.bond_coeffs_present = true;
        section = Section::kOther;
    } else if (section_line(line, "Angle Coeffs")) {
        summary.angle_coeffs_present = true;
        section = Section::kOther;
    } else if (section_line(line, "Dihedral Coeffs")) {
        summary.dihedral_coeffs_present = true;
        section = Section::kOther;
    } else if (section_line(line, "Atoms")) {
        section = Section::kAtoms;
    } else if (section_line(line, "Velocities")) {
        section = Section::kVelocities;
    } else if (section_line(line, "Bonds")) {
        section = Section::kBonds;
    } else if (section_line(line, "Angles")) {
        section = Section::kAngles;
    } else if (section_line(line, "Dihedrals")) {
        section = Section::kDihedrals;
    } else if (section_line(line, "Impropers") ||
               section_line(line, "Improper Coeffs")) {
        section = Section::kOther;
    } else {
        return false;
    }
    return true;
}

void parse_header_line(const std::string &line, DataSummary &summary) {
    const std::regex timestep_pattern("timestep\\s*=\\s*([0-9]+)");
    std::smatch timestep_match;
    if (std::regex_search(line, timestep_match, timestep_pattern)) {
        summary.timestep = std::stoll(timestep_match[1].str());
    }

    long long count = 0;
    std::string word1;
    std::string word2;
    std::istringstream count_line(line);
    if (count_line >> count >> word1) {
        if (word1 == "atoms") {
            summary.declared_atoms = count;
            return;
        }
        if (word1 == "bonds") {
            summary.declared_bonds = count;
            return;
        }
        if (word1 == "angles") {
            summary.declared_angles = count;
            return;
        }
        if (word1 == "dihedrals") {
            summary.declared_dihedrals = count;
            return;
        }
        if (count_line >> word2 && word2 == "types") {
            if (word1 == "atom") {
                summary.atom_types = static_cast<int>(count);
            } else if (word1 == "bond") {
                summary.bond_types = static_cast<int>(count);
            } else if (word1 == "angle") {
                summary.angle_types = static_cast<int>(count);
            } else if (word1 == "dihedral") {
                summary.dihedral_types = static_cast<int>(count);
            }
            return;
        }
    }

    double lo = 0.0;
    double hi = 0.0;
    std::string lower_label;
    std::string upper_label;
    std::istringstream bounds_line(line);
    if (bounds_line >> lo >> hi >> lower_label >> upper_label) {
        Bounds *bounds = nullptr;
        if (lower_label == "xlo" && upper_label == "xhi") {
            bounds = &summary.x;
        } else if (lower_label == "ylo" && upper_label == "yhi") {
            bounds = &summary.y;
        } else if (lower_label == "zlo" && upper_label == "zhi") {
            bounds = &summary.z;
        }
        if (bounds != nullptr) {
            bounds->lo = lo;
            bounds->hi = hi;
            bounds->seen = true;
        }
    }
}

DataSummary parse_data_file(
    const std::string &path,
    const ModelInfo &info,
    DisjointSet &sets) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open LAMMPS data file: " + path);
    }

    DataSummary summary;
    summary.component_beads.assign(kComponentCount, 0);
    summary.component_masses.assign(kComponentCount, 0.0);
    summary.molecules.resize(
        static_cast<std::size_t>(info.total_molecules + 1));
    const std::vector<long long> ends = component_ends(info);

    Section section = Section::kHeader;
    std::string line;
    long long line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        const std::string clean = trim(line);
        if (clean.empty() || clean.front() == '#') {
            continue;
        }
        if (update_section(clean, section, summary)) {
            if (section == Section::kAtoms) {
                if (summary.declared_atoms <= 0 || summary.atom_types <= 0) {
                    throw std::runtime_error(
                        "atom and atom-type counts must precede the Atoms section");
                }
                if (summary.masses.size() !=
                    static_cast<std::size_t>(summary.atom_types + 1)) {
                    throw std::runtime_error(
                        "incomplete Masses section before Atoms");
                }
                summary.atoms.resize(
                    static_cast<std::size_t>(summary.declared_atoms + 1));
                summary.atom_type_counts.assign(
                    static_cast<std::size_t>(summary.atom_types + 1), 0);
            }
            continue;
        }

        if (section == Section::kHeader) {
            parse_header_line(clean, summary);
            continue;
        }
        if (section == Section::kMasses) {
            int type = 0;
            double mass = 0.0;
            std::istringstream fields(clean);
            if (!(fields >> type >> mass)) {
                throw std::runtime_error(
                    "invalid Masses entry at line " + std::to_string(line_number));
            }
            if (summary.atom_types <= 0 || type < 1 ||
                type > summary.atom_types || mass <= 0.0) {
                throw std::runtime_error(
                    "invalid atom type or mass at line " +
                    std::to_string(line_number));
            }
            if (summary.masses.empty()) {
                summary.masses.assign(
                    static_cast<std::size_t>(summary.atom_types + 1), 0.0);
            }
            if (summary.masses[static_cast<std::size_t>(type)] != 0.0) {
                throw std::runtime_error(
                    "duplicate mass for atom type " + std::to_string(type));
            }
            summary.masses[static_cast<std::size_t>(type)] = mass;
            continue;
        }
        if (section == Section::kAtoms) {
            long long id = 0;
            long long molecule = 0;
            int type = 0;
            double charge = 0.0;
            double x = 0.0;
            double y = 0.0;
            double z = 0.0;
            long long ix = 0;
            long long iy = 0;
            long long iz = 0;
            std::istringstream fields(clean);
            if (!(fields >> id >> molecule >> type >> charge >> x >> y >> z)) {
                throw std::runtime_error(
                    "invalid atom record at line " + std::to_string(line_number));
            }
            const bool have_images = static_cast<bool>(fields >> ix >> iy >> iz);
            if (id < 1 || id > summary.declared_atoms ||
                summary.atoms[static_cast<std::size_t>(id)].seen) {
                throw std::runtime_error(
                    "invalid or duplicate atom ID at line " +
                    std::to_string(line_number));
            }
            if (type < 1 || type > summary.atom_types ||
                summary.masses[static_cast<std::size_t>(type)] <= 0.0) {
                throw std::runtime_error(
                    "atom has invalid type at line " +
                    std::to_string(line_number));
            }
            const int component = component_for_molecule(molecule, ends);
            const double mass = summary.masses[static_cast<std::size_t>(type)];
            AtomRecord &atom = summary.atoms[static_cast<std::size_t>(id)];
            atom.seen = true;
            atom.molecule = molecule;
            atom.type = type;
            atom.mass = mass;
            atom.x = x;
            atom.y = y;
            atom.z = z;
            atom.ix = ix;
            atom.iy = iy;
            atom.iz = iz;
            atom.have_images = have_images;

            MoleculeRecord &molecule_record =
                summary.molecules[static_cast<std::size_t>(molecule)];
            ++molecule_record.beads;
            molecule_record.mass += mass;
            ++summary.component_beads[static_cast<std::size_t>(component)];
            summary.component_masses[static_cast<std::size_t>(component)] += mass;
            ++summary.atom_type_counts[static_cast<std::size_t>(type)];
            summary.total_mass_g_per_mol += mass;
            summary.charge_sum += charge;
            summary.x_min = std::min(summary.x_min, x);
            summary.x_max = std::max(summary.x_max, x);
            summary.y_min = std::min(summary.y_min, y);
            summary.y_max = std::max(summary.y_max, y);
            summary.z_min = std::min(summary.z_min, z);
            summary.z_max = std::max(summary.z_max, z);
            if (have_images && (ix != 0 || iy != 0 || iz != 0)) {
                ++summary.nonzero_image_atoms;
            }
            ++summary.atoms_read;
            continue;
        }
        if (section == Section::kVelocities) {
            long long id = 0;
            double vx = 0.0;
            double vy = 0.0;
            double vz = 0.0;
            std::istringstream fields(clean);
            if (!(fields >> id >> vx >> vy >> vz) ||
                id < 1 || id > summary.declared_atoms ||
                !summary.atoms[static_cast<std::size_t>(id)].seen) {
                throw std::runtime_error(
                    "invalid velocity record at line " +
                    std::to_string(line_number));
            }
            const double mass =
                summary.atoms[static_cast<std::size_t>(id)].mass;
            summary.velocity_mass += mass;
            summary.momentum_x += mass * vx;
            summary.momentum_y += mass * vy;
            summary.momentum_z += mass * vz;
            summary.mass_velocity_squared +=
                mass * (vx * vx + vy * vy + vz * vz);
            ++summary.velocities_read;
            continue;
        }
        if (section == Section::kBonds) {
            long long id = 0;
            int type = 0;
            long long first_atom = 0;
            long long second_atom = 0;
            std::istringstream fields(clean);
            if (!(fields >> id >> type >> first_atom >> second_atom) ||
                first_atom < 1 || first_atom > summary.declared_atoms ||
                second_atom < 1 || second_atom > summary.declared_atoms ||
                !summary.atoms[static_cast<std::size_t>(first_atom)].seen ||
                !summary.atoms[static_cast<std::size_t>(second_atom)].seen) {
                throw std::runtime_error(
                    "invalid bond record at line " + std::to_string(line_number));
            }
            ++summary.bond_type_counts[type];
            ++summary.bonds_read;
            if (type == info.crosslink_bond_type) {
                const long long first_molecule =
                    summary.atoms[static_cast<std::size_t>(first_atom)].molecule;
                const long long second_molecule =
                    summary.atoms[static_cast<std::size_t>(second_atom)].molecule;
                if (first_molecule == second_molecule) {
                    ++summary.crosslink_bonds_intra;
                } else {
                    ++summary.crosslink_bonds_inter;
                    ++summary.molecules[
                        static_cast<std::size_t>(first_molecule)]
                          .inter_crosslink_degree;
                    summary.molecules[
                        static_cast<std::size_t>(first_molecule)]
                        .crosslink_attachment_atoms.push_back(first_atom);
                    ++summary.molecules[
                        static_cast<std::size_t>(second_molecule)]
                          .inter_crosslink_degree;
                    summary.molecules[
                        static_cast<std::size_t>(second_molecule)]
                        .crosslink_attachment_atoms.push_back(second_atom);
                    sets.unite(first_molecule, second_molecule);
                }
            }
            continue;
        }
        if (section == Section::kAngles) {
            long long id = 0;
            int type = 0;
            std::istringstream fields(clean);
            if (!(fields >> id >> type)) {
                throw std::runtime_error(
                    "invalid angle record at line " + std::to_string(line_number));
            }
            ++summary.angle_type_counts[type];
            ++summary.angles_read;
            continue;
        }
        if (section == Section::kDihedrals) {
            long long id = 0;
            int type = 0;
            std::istringstream fields(clean);
            if (!(fields >> id >> type)) {
                throw std::runtime_error(
                    "invalid dihedral record at line " +
                    std::to_string(line_number));
            }
            ++summary.dihedral_type_counts[type];
            ++summary.dihedrals_read;
        }
    }

    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed while reading LAMMPS data file");
    }
    return summary;
}

std::vector<NetworkComponent> network_components(
    const ModelInfo &info,
    const DataSummary &data,
    DisjointSet &sets) {
    std::map<long long, NetworkComponent> by_root;
    const std::vector<long long> ends = component_ends(info);
    for (long long molecule = 1;
         molecule <= ends[static_cast<std::size_t>(kCrosslinker)];
         ++molecule) {
        const int component = component_for_molecule(molecule, ends);
        const long long root = sets.find(molecule);
        NetworkComponent &entry = by_root[root];
        entry.root = root;
        ++entry.precursor_molecules;
        entry.precursor_beads +=
            data.molecules[static_cast<std::size_t>(molecule)].beads;
        entry.precursor_mass +=
            data.molecules[static_cast<std::size_t>(molecule)].mass;
        if (component == kNetwork) {
            ++entry.network_molecules;
        } else if (component == kCrosslinker) {
            ++entry.crosslinker_molecules;
        }
    }
    std::vector<NetworkComponent> result;
    for (const auto &item : by_root) {
        result.push_back(item.second);
    }
    std::sort(
        result.begin(), result.end(),
        [](const NetworkComponent &left, const NetworkComponent &right) {
            if (left.precursor_molecules != right.precursor_molecules) {
                return left.precursor_molecules > right.precursor_molecules;
            }
            return left.precursor_beads > right.precursor_beads;
        });
    return result;
}

std::map<long long, long long> degree_distribution(
    const ModelInfo &info,
    const DataSummary &data,
    int requested_component) {
    std::map<long long, long long> distribution;
    const std::vector<long long> ends = component_ends(info);
    const long long first =
        requested_component == 0
            ? 1
            : ends[static_cast<std::size_t>(requested_component - 1)] + 1;
    const long long last = ends[static_cast<std::size_t>(requested_component)];
    for (long long molecule = first; molecule <= last; ++molecule) {
        ++distribution[data.molecules[static_cast<std::size_t>(molecule)]
                           .inter_crosslink_degree];
    }
    return distribution;
}

double minimum_image_displacement(double displacement, double box_length) {
    return displacement - std::round(displacement / box_length) * box_length;
}

std::vector<EndToEndRecord> strand_end_to_end_records(
    const ModelInfo &info,
    const DataSummary &data,
    DisjointSet &sets,
    long long largest_component_root,
    long long &non_distinct_attachment_count) {
    std::vector<EndToEndRecord> records;
    non_distinct_attachment_count = 0;
    const std::vector<long long> ends = component_ends(info);
    const double lx = data.x.hi - data.x.lo;
    const double ly = data.y.hi - data.y.lo;
    const double lz = data.z.hi - data.z.lo;

    for (long long molecule = 1;
         molecule <= ends[static_cast<std::size_t>(kNetwork)];
         ++molecule) {
        const MoleculeRecord &molecule_record =
            data.molecules[static_cast<std::size_t>(molecule)];
        if (molecule_record.inter_crosslink_degree != 2) {
            continue;
        }
        std::vector<long long> attachment_atoms =
            molecule_record.crosslink_attachment_atoms;
        std::sort(attachment_atoms.begin(), attachment_atoms.end());
        attachment_atoms.erase(
            std::unique(attachment_atoms.begin(), attachment_atoms.end()),
            attachment_atoms.end());
        if (attachment_atoms.size() != 2) {
            ++non_distinct_attachment_count;
            continue;
        }

        const AtomRecord &first =
            data.atoms[static_cast<std::size_t>(attachment_atoms[0])];
        const AtomRecord &second =
            data.atoms[static_cast<std::size_t>(attachment_atoms[1])];
        EndToEndRecord record;
        record.molecule = molecule;
        record.first_atom = attachment_atoms[0];
        record.second_atom = attachment_atoms[1];
        record.used_unwrapped_coordinates =
            first.have_images && second.have_images;
        if (record.used_unwrapped_coordinates) {
            record.dx =
                (second.x + static_cast<double>(second.ix) * lx) -
                (first.x + static_cast<double>(first.ix) * lx);
            record.dy =
                (second.y + static_cast<double>(second.iy) * ly) -
                (first.y + static_cast<double>(first.iy) * ly);
            record.dz =
                (second.z + static_cast<double>(second.iz) * lz) -
                (first.z + static_cast<double>(first.iz) * lz);
        } else {
            record.dx = minimum_image_displacement(second.x - first.x, lx);
            record.dy = minimum_image_displacement(second.y - first.y, ly);
            record.dz = minimum_image_displacement(second.z - first.z, lz);
        }
        record.distance = std::sqrt(
            record.dx * record.dx +
            record.dy * record.dy +
            record.dz * record.dz);
        record.in_largest_component =
            sets.find(molecule) == largest_component_root;
        records.push_back(record);
    }
    return records;
}

double percentile(const std::vector<double> &sorted_values, double fraction) {
    if (sorted_values.empty()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double position =
        fraction * static_cast<double>(sorted_values.size() - 1);
    const std::size_t lower = static_cast<std::size_t>(std::floor(position));
    const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
    const double weight = position - static_cast<double>(lower);
    return sorted_values[lower] * (1.0 - weight) +
           sorted_values[upper] * weight;
}

DistanceStatistics distance_statistics(
    const std::vector<EndToEndRecord> &records,
    bool largest_component_only) {
    DistanceStatistics statistics;
    std::vector<double> distances;
    for (const EndToEndRecord &record : records) {
        if (largest_component_only && !record.in_largest_component) {
            continue;
        }
        distances.push_back(record.distance);
        statistics.mean += record.distance;
        statistics.mean_squared += record.distance * record.distance;
    }
    statistics.count = static_cast<long long>(distances.size());
    if (distances.empty()) {
        return statistics;
    }
    statistics.mean /= static_cast<double>(statistics.count);
    statistics.mean_squared /= static_cast<double>(statistics.count);
    statistics.rms = std::sqrt(statistics.mean_squared);
    statistics.standard_deviation = std::sqrt(std::max(
        0.0,
        statistics.mean_squared - statistics.mean * statistics.mean));
    std::sort(distances.begin(), distances.end());
    statistics.minimum = distances.front();
    statistics.first_quartile = percentile(distances, 0.25);
    statistics.median = percentile(distances, 0.50);
    statistics.third_quartile = percentile(distances, 0.75);
    statistics.maximum = distances.back();
    return statistics;
}

void write_strand_output(
    const std::string &path,
    const std::vector<EndToEndRecord> &records) {
    std::ofstream output(path);
    if (!output) {
        throw std::runtime_error(
            "cannot open strand end-to-end output file: " + path);
    }
    output
        << "# molecule_id attachment_atom_1 attachment_atom_2 "
        << "dx_A dy_A dz_A end_to_end_A in_largest_component "
        << "used_unwrapped_coordinates\n";
    output << std::scientific << std::setprecision(10);
    for (const EndToEndRecord &record : records) {
        output
            << record.molecule << " "
            << record.first_atom << " "
            << record.second_atom << " "
            << record.dx << " "
            << record.dy << " "
            << record.dz << " "
            << record.distance << " "
            << (record.in_largest_component ? 1 : 0) << " "
            << (record.used_unwrapped_coordinates ? 1 : 0) << "\n";
    }
    if (!output) {
        throw std::runtime_error(
            "failed while writing strand end-to-end output: " + path);
    }
}

std::string yes_no(bool value) {
    return value ? "yes" : "no";
}

template <typename Key>
std::string distribution_text(const std::map<Key, long long> &values) {
    std::ostringstream output;
    bool first = true;
    for (const auto &entry : values) {
        if (!first) {
            output << ", ";
        }
        output << entry.first << ":" << entry.second;
        first = false;
    }
    return output.str();
}

void validate(
    const ModelInfo &info,
    const DataSummary &data,
    std::vector<std::string> &warnings) {
    if (data.declared_atoms != info.total_beads) {
        warnings.push_back(
            "data atom count does not match info total_beads");
    }
    if (data.atoms_read != data.declared_atoms) {
        throw std::runtime_error(
            "Atoms section count differs from the data-file header");
    }
    if (data.velocities_read != 0 &&
        data.velocities_read != data.declared_atoms) {
        warnings.push_back(
            "velocity count differs from atom count; temperature is omitted");
    }
    if (data.bonds_read != data.declared_bonds ||
        data.angles_read != data.declared_angles ||
        data.dihedrals_read != data.declared_dihedrals) {
        throw std::runtime_error(
            "one or more topology section counts differ from the header");
    }
    if (!data.x.seen || !data.y.seen || !data.z.seen ||
        data.x.hi <= data.x.lo || data.y.hi <= data.y.lo ||
        data.z.hi <= data.z.lo) {
        throw std::runtime_error("invalid or missing orthogonal box bounds");
    }
    for (int component = 0; component < kComponentCount; ++component) {
        if (data.component_beads[static_cast<std::size_t>(component)] !=
            info.components[static_cast<std::size_t>(component)].beads) {
            warnings.push_back(
                "observed bead count differs from info for " +
                info.components[static_cast<std::size_t>(component)].label);
        }
    }
    const double mass_scale =
        std::max(1.0, std::fabs(info.total_mass_g_per_mol));
    if (std::fabs(data.total_mass_g_per_mol - info.total_mass_g_per_mol) >
        1.0e-8 * mass_scale) {
        warnings.push_back(
            "mass calculated from data differs from info total mass");
    }
    if (data.timestep >= 0 && data.timestep != info.total_run_steps) {
        warnings.push_back(
            "snapshot timestep differs from info total_run_steps");
    }
    const long long bonds_created =
        data.declared_bonds - info.initial_bonds;
    if (bonds_created != data.crosslink_bonds_inter) {
        warnings.push_back(
            "net bond increase differs from the number of intermolecular "
            "formulation-crosslink bonds; conversion estimates need review");
    }
    if (data.atom_type_counts.size() > 3) {
        if (data.atom_type_counts[2] + bonds_created !=
            info.type2_sites_total) {
            warnings.push_back(
                "final type-2 atoms plus created bonds differ from the "
                "initial type-2 site count");
        }
        if (data.atom_type_counts[3] + bonds_created !=
            info.type3_sites_total) {
            warnings.push_back(
                "final type-3 atoms plus created bonds differ from the "
                "initial type-3 site count");
        }
    }
    if (!data.pair_coeffs_present) {
        warnings.push_back(
            "Pair Coeffs are absent; use the matching LAMMPS input to recover "
            "the nonbonded force field");
    }
}

void write_report(
    std::ostream &output,
    const Options &options,
    const ModelInfo &info,
    const DataSummary &data,
    const std::vector<NetworkComponent> &components,
    const std::vector<EndToEndRecord> &end_to_end_records,
    long long non_distinct_attachment_count,
    const std::vector<std::string> &warnings) {
    const double lx = data.x.hi - data.x.lo;
    const double ly = data.y.hi - data.y.lo;
    const double lz = data.z.hi - data.z.lo;
    const double volume_a3 = lx * ly * lz;
    const double density =
        (data.total_mass_g_per_mol / kAvogadro) /
        (volume_a3 * kAngstromCubedToCmCubed);
    const double duration_ns =
        static_cast<double>(info.total_run_steps) * info.timestep_fs / 1.0e6;
    const double reaction_duration_ns =
        static_cast<double>(info.bond_creation_active_steps) *
        info.timestep_fs / 1.0e6;
    const long long bonds_created =
        data.declared_bonds - info.initial_bonds;

    double temperature = std::numeric_limits<double>::quiet_NaN();
    double com_vx = std::numeric_limits<double>::quiet_NaN();
    double com_vy = std::numeric_limits<double>::quiet_NaN();
    double com_vz = std::numeric_limits<double>::quiet_NaN();
    if (data.velocities_read == data.declared_atoms &&
        data.velocity_mass > 0.0 && data.velocities_read > 1) {
        com_vx = data.momentum_x / data.velocity_mass;
        com_vy = data.momentum_y / data.velocity_mass;
        com_vz = data.momentum_z / data.velocity_mass;
        const double corrected_mass_velocity_squared =
            data.mass_velocity_squared -
            data.velocity_mass *
                (com_vx * com_vx + com_vy * com_vy + com_vz * com_vz);
        const double kinetic_energy =
            0.5 * kLammpsRealMvv2e * corrected_mass_velocity_squared;
        const long long degrees_of_freedom = 3 * data.velocities_read - 3;
        temperature =
            2.0 * kinetic_energy /
            (static_cast<double>(degrees_of_freedom) *
             kBoltzmannKcalPerMolK);
    }

    const NetworkComponent &largest = components.front();
    const long long precursor_molecules =
        info.components[kNetwork].molecules +
        info.components[kCrosslinker].molecules;
    const long long precursor_beads =
        info.components[kNetwork].beads +
        info.components[kCrosslinker].beads;
    const double precursor_mass =
        data.component_masses[kNetwork] +
        data.component_masses[kCrosslinker];
    const DistanceStatistics all_end_to_end =
        distance_statistics(end_to_end_records, false);
    const DistanceStatistics largest_end_to_end =
        distance_statistics(end_to_end_records, true);
    long long unwrapped_distance_count = 0;
    for (const EndToEndRecord &record : end_to_end_records) {
        if (record.used_unwrapped_coordinates) {
            ++unwrapped_distance_count;
        }
    }

    long long linked_component_count = 0;
    long long isolated_precursor_count = 0;
    for (const NetworkComponent &component : components) {
        if (component.precursor_molecules > 1) {
            ++linked_component_count;
        } else {
            ++isolated_precursor_count;
        }
    }

    output << std::fixed << std::setprecision(6);
    output << "FINAL LAMMPS SNAPSHOT ANALYSIS\n";
    output << "================================\n";
    output << "Case                         : " << info.case_name << "\n";
    output << "Formulation                  : " << info.formulation << "\n";
    output << "Info format/version          : "
           << info.format << " / " << info.format_version << "\n";
    output << "Geometry                     : " << info.geometry << "\n";
    output << "Data file                    : " << options.data_file << "\n";
    output << "Info file                    : " << options.info_file << "\n\n";

    output << "RUN AND FINAL STATE\n";
    output << "Snapshot timestep            : " << data.timestep << "\n";
    output << "Configured total steps       : " << info.total_run_steps << "\n";
    output << "Timestep (fs)                : " << info.timestep_fs << "\n";
    output << "Configured duration (ns)     : " << duration_ns << "\n";
    output << "Bond-creation duration (ns)  : " << reaction_duration_ns << "\n";
    output << "Box lengths (A)              : "
           << lx << " " << ly << " " << lz << "\n";
    output << "Volume (A^3)                 : " << volume_a3 << "\n";
    output << "Initial density (g/cm^3)     : " << info.initial_density << "\n";
    output << "Compression target (g/cm^3) : "
           << info.target_compressed_density << "\n";
    output << "Final density (g/cm^3)       : " << density << "\n";
    if (std::isfinite(temperature)) {
        output << "Velocity temperature (K)     : " << temperature << "\n";
        output << "Target final temperature (K) : "
               << info.final_temperature_k << "\n";
        output << "COM velocity (A/fs)          : "
               << std::scientific << com_vx << " " << com_vy << " " << com_vz
               << std::fixed << "\n";
    } else {
        output << "Velocity temperature (K)     : unavailable\n";
    }
    output << "\n";

    output << "COMPOSITION\n";
    output << "Atoms                        : " << data.declared_atoms << "\n";
    output << "Molecules                    : " << info.total_molecules << "\n";
    output << "Mass (g/mol equivalent)      : "
           << std::setprecision(4) << data.total_mass_g_per_mol
           << std::setprecision(6) << "\n";
    output << "Requested filler wt%         : "
           << info.requested_filler_weight_percent << "\n";
    output << "Info realized filler wt%     : "
           << info.realized_filler_weight_percent << "\n";
    output << "Observed components:\n";
    for (int component = 0; component < kComponentCount; ++component) {
        const ComponentInfo &component_info =
            info.components[static_cast<std::size_t>(component)];
        const double observed_weight =
            100.0 *
            data.component_masses[static_cast<std::size_t>(component)] /
            data.total_mass_g_per_mol;
        output << "  " << std::left << std::setw(22) << component_info.label
               << std::right
               << " molecules=" << std::setw(6) << component_info.molecules
               << " beads=" << std::setw(8)
               << data.component_beads[static_cast<std::size_t>(component)]
               << " wt%=" << std::setw(10) << observed_weight << "\n";
    }
    output << "Atom type populations:\n";
    for (int type = 1; type <= data.atom_types; ++type) {
        output << "  type " << type
               << "  mass=" << data.masses[static_cast<std::size_t>(type)]
               << "  count="
               << data.atom_type_counts[static_cast<std::size_t>(type)] << "\n";
    }
    output << "\n";

    output << "TOPOLOGY AND REACTION\n";
    output << "Initial bonds from info      : " << info.initial_bonds << "\n";
    output << "Final bonds                  : " << data.declared_bonds << "\n";
    output << "Net bonds created            : " << bonds_created << "\n";
    output << "Crosslink bond type          : "
           << info.crosslink_bond_type << "\n";
    output << "Intermolecular type-"
           << info.crosslink_bond_type << " bonds : "
           << data.crosslink_bonds_inter << "\n";
    output << "Intramolecular type-"
           << info.crosslink_bond_type << " bonds : "
           << data.crosslink_bonds_intra << "\n";
    if (info.type2_sites_total > 0) {
        output << "Type-2 site conversion (%)   : "
               << 100.0 * static_cast<double>(bonds_created) /
                      static_cast<double>(info.type2_sites_total)
               << "\n";
    }
    if (info.type3_sites_total > 0) {
        output << "Type-3 site conversion (%)   : "
               << 100.0 * static_cast<double>(bonds_created) /
                      static_cast<double>(info.type3_sites_total)
               << "\n";
    }
    if (data.atom_type_counts.size() > 3) {
        output << "Final reactive type-2 atoms  : "
               << data.atom_type_counts[2] << "\n";
        output << "Final reactive type-3 atoms  : "
               << data.atom_type_counts[3] << "\n";
    }
    output << "Network strand degrees       : "
           << distribution_text(degree_distribution(
                  info, data, kNetwork))
           << "  (degree:molecules)\n";
    output << "Crosslinker degrees          : "
           << distribution_text(degree_distribution(
                  info, data, kCrosslinker))
           << "  (degree:molecules)\n";
    output << "Moderator degrees            : "
           << distribution_text(degree_distribution(
                  info, data, kModerator))
           << "  (degree:molecules)\n";
    output << "Bond type populations        : "
           << distribution_text(data.bond_type_counts) << "\n";
    output << "Angle type populations       : "
           << distribution_text(data.angle_type_counts) << "\n";
    output << "Dihedral type populations    : "
           << distribution_text(data.dihedral_type_counts) << "\n\n";

    output << "CROSSLINK NETWORK CONNECTIVITY\n";
    output << "Precursor molecules          : " << precursor_molecules << "\n";
    output << "Linked components (>1 mol)   : " << linked_component_count << "\n";
    output << "Isolated precursor molecules : " << isolated_precursor_count << "\n";
    output << "Largest component molecules  : "
           << largest.precursor_molecules << " / " << precursor_molecules
           << " ("
           << 100.0 * static_cast<double>(largest.precursor_molecules) /
                  static_cast<double>(precursor_molecules)
           << "%)\n";
    output << "  network/crosslinker         : "
           << largest.network_molecules << " / "
           << largest.crosslinker_molecules << "\n";
    output << "Largest component beads      : "
           << largest.precursor_beads << " / " << precursor_beads
           << " ("
           << 100.0 * static_cast<double>(largest.precursor_beads) /
                  static_cast<double>(precursor_beads)
           << "%)\n";
    output << "Largest component mass frac  : "
           << 100.0 * largest.precursor_mass / precursor_mass << "%\n";
    if (components.size() > 1) {
        output << "Other precursor components   :\n";
        for (std::size_t index = 1; index < components.size(); ++index) {
            const NetworkComponent &component = components[index];
            output << "  molecules=" << component.precursor_molecules
                   << " network=" << component.network_molecules
                   << " crosslinkers=" << component.crosslinker_molecules
                   << " beads=" << component.precursor_beads << "\n";
        }
    }
    output << "\n";

    output << "DOUBLY CONNECTED STRAND END-TO-END DISTANCE\n";
    output
        << "Definition                    : distance between the two distinct strand\n"
        << "                                atoms carrying intermolecular crosslinks\n";
    output << "Eligible degree-2 strands    : "
           << all_end_to_end.count + non_distinct_attachment_count << "\n";
    output << "Distances calculated         : "
           << all_end_to_end.count << "\n";
    output << "In largest network component : "
           << largest_end_to_end.count << "\n";
    output << "Using unwrapped coordinates  : "
           << unwrapped_distance_count << " / " << all_end_to_end.count << "\n";
    output << "Rejected non-distinct sites  : "
           << non_distinct_attachment_count << "\n";
    if (all_end_to_end.count > 0) {
        output << "Mean R_ee (A)                : "
               << all_end_to_end.mean << "\n";
        output << "Mean R_ee^2 (A^2)            : "
               << all_end_to_end.mean_squared << "\n";
        output << "RMS R_ee (A)                 : "
               << all_end_to_end.rms << "\n";
        output << "Population std. dev. (A)     : "
               << all_end_to_end.standard_deviation << "\n";
        output << "Min / Q1 / median (A)        : "
               << all_end_to_end.minimum << " / "
               << all_end_to_end.first_quartile << " / "
               << all_end_to_end.median << "\n";
        output << "Q3 / max (A)                 : "
               << all_end_to_end.third_quartile << " / "
               << all_end_to_end.maximum << "\n";
    }
    if (largest_end_to_end.count > 0 &&
        largest_end_to_end.count != all_end_to_end.count) {
        output << "Largest-component mean (A)   : "
               << largest_end_to_end.mean << "\n";
        output << "Largest-component RMS (A)    : "
               << largest_end_to_end.rms << "\n";
    }
    if (!options.strand_output_file.empty()) {
        output << "Per-strand output            : "
               << options.strand_output_file << "\n";
    }
    output << "\n";

    output << "DATA COMPLETENESS\n";
    output << "Atoms parsed                 : " << data.atoms_read
           << " / " << data.declared_atoms << "\n";
    output << "Velocities parsed            : " << data.velocities_read
           << " / " << data.declared_atoms << "\n";
    output << "Bonds/angles/dihedrals       : "
           << data.bonds_read << " / " << data.angles_read << " / "
           << data.dihedrals_read << "\n";
    output << "Atoms with nonzero images    : "
           << data.nonzero_image_atoms << "\n";
    output << "Coordinate ranges x/y/z (A)  : ["
           << data.x_min << ", " << data.x_max << "] ["
           << data.y_min << ", " << data.y_max << "] ["
           << data.z_min << ", " << data.z_max << "]\n";
    output << "Pair Coeffs present          : "
           << yes_no(data.pair_coeffs_present) << "\n";
    output << "Bond Coeffs present          : "
           << yes_no(data.bond_coeffs_present) << "\n";
    output << "Angle Coeffs present         : "
           << yes_no(data.angle_coeffs_present) << "\n";
    output << "Dihedral Coeffs present      : "
           << yes_no(data.dihedral_coeffs_present) << "\n";
    output << "Charge sum                   : " << data.charge_sum << "\n";
    output << "\n";

    output << "CHECKS AND LIMITATIONS\n";
    if (warnings.empty()) {
        output << "PASS: no consistency warnings.\n";
    } else {
        for (const std::string &warning : warnings) {
            output << "WARNING: " << warning << "\n";
        }
    }
    output
        << "NOTE: one final snapshot supports static composition, density, and\n"
        << "      connectivity analysis. Equilibration, pressure/energy stability,\n"
        << "      diffusion, and reaction timing require log or trajectory files.\n";
}

void print_help(const char *program) {
    std::cout
        << "Usage:\n"
        << "  " << program
        << " FINAL_DATA_FILE MODEL_INFO_FILE [--output REPORT_FILE]\n"
        << "      [--strand-output TABLE_FILE]\n\n"
        << "Analyze a final LAMMPS atom_style full data snapshot using its matching\n"
        << "V22/V35 model-info version-2 JSON file. The report covers final density,\n"
        << "velocity-derived temperature, reaction conversion, crosslink-network\n"
        << "connectivity, doubly connected strand end-to-end distances, topology\n"
        << "counts, and data completeness.\n\n"
        << "Options:\n"
        << "  --output FILE         Write the report instead of using standard output\n"
        << "  --strand-output FILE  Write one end-to-end record per eligible strand\n"
        << "  --help                Show this command reference\n\n"
        << "Example:\n"
        << "  " << program
        << " data.case.npt_eq case.info --output case.analysis.txt \\\n"
        << "      --strand-output case.strand_ree.dat\n";
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
        if (argument == "--output") {
            if (index + 1 >= argc) {
                throw std::runtime_error("--output requires a file path");
            }
            options.output_file = argv[++index];
        } else if (argument == "--strand-output") {
            if (index + 1 >= argc) {
                throw std::runtime_error(
                    "--strand-output requires a file path");
            }
            options.strand_output_file = argv[++index];
        } else if (begins_with(argument, "--")) {
            throw std::runtime_error("unknown option: " + argument);
        } else {
            positional.push_back(argument);
        }
    }
    if (positional.size() != 2) {
        throw std::runtime_error(
            "expected FINAL_DATA_FILE and MODEL_INFO_FILE; use --help");
    }
    options.data_file = positional[0];
    options.info_file = positional[1];
    return options;
}

}  // namespace

int main(int argc, char **argv) {
    try {
        const Options options = parse_options(argc, argv);
        const ModelInfo info = parse_model_info(options.info_file);
        DisjointSet sets(
            static_cast<std::size_t>(info.total_molecules + 1));
        const DataSummary data =
            parse_data_file(options.data_file, info, sets);
        std::vector<std::string> warnings;
        validate(info, data, warnings);
        const std::vector<NetworkComponent> components =
            network_components(info, data, sets);
        if (components.empty()) {
            throw std::runtime_error(
                "no network-strand or crosslinker molecules were found");
        }
        long long non_distinct_attachment_count = 0;
        const std::vector<EndToEndRecord> end_to_end_records =
            strand_end_to_end_records(
                info,
                data,
                sets,
                components.front().root,
                non_distinct_attachment_count);
        if (non_distinct_attachment_count > 0) {
            warnings.push_back(
                std::to_string(non_distinct_attachment_count) +
                " degree-2 strands did not have two distinct crosslink "
                "attachment atoms");
        }
        const long long fallback_distance_count =
            static_cast<long long>(std::count_if(
                end_to_end_records.begin(),
                end_to_end_records.end(),
                [](const EndToEndRecord &record) {
                    return !record.used_unwrapped_coordinates;
                }));
        if (fallback_distance_count > 0) {
            warnings.push_back(
                std::to_string(fallback_distance_count) +
                " strand distances used minimum-image coordinates because "
                "image flags were unavailable");
        }
        if (!options.strand_output_file.empty()) {
            write_strand_output(
                options.strand_output_file, end_to_end_records);
            std::cout << "Wrote strand end-to-end table: "
                      << options.strand_output_file << "\n";
        }

        if (options.output_file.empty()) {
            write_report(
                std::cout,
                options,
                info,
                data,
                components,
                end_to_end_records,
                non_distinct_attachment_count,
                warnings);
        } else {
            std::ofstream output(options.output_file);
            if (!output) {
                throw std::runtime_error(
                    "cannot open report file: " + options.output_file);
            }
            write_report(
                output,
                options,
                info,
                data,
                components,
                end_to_end_records,
                non_distinct_attachment_count,
                warnings);
            if (!output) {
                throw std::runtime_error(
                    "failed while writing report: " + options.output_file);
            }
            std::cout << "Wrote analysis report: "
                      << options.output_file << "\n";
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "Error: " << error.what() << "\n";
        return 1;
    }
}
