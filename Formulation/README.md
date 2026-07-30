# Shared formulation oil component

`silicone_oil_component.hpp` is the shared implementation used by both
`V22/v22_generator.cpp` and `V35/v35_generator.cpp`.

It contains the tested standalone Oil model's:

- PDMS, PMPS, and PDMS–PMPS sequence construction;
- non-crossing local oil conformations and pendant placement;
- atom types 1, 4, and 5;
- oil bond, angle, and dihedral classification;
- 800 K and 300 K nonbonded matrices;
- `0.579966 ×` geometric-epsilon and arithmetic-sigma DMS–MPS mixing rule;
- mass constants and mass-aware chain composition;
- whole-chain placement inside orthogonal bulk or film boxes.

In the formulation generators, oil chains are restricted to the central 40%
of `Lz`. Their placement is checked against the bottom-up network and top-down
cross-linkers. Star moderators occupy a separate upper-middle band and are
rejected if any moderator bead is too close to a previously placed bead.

The formulation reserves bond type 2 for network crosslinks and moderator
bonds. Consequently, the standalone Oil bond mapping is remapped during
integration:

```text
integrated type 1 = DMS or mixed backbone
integrated type 2 = formulation crosslink/moderator
integrated type 3 = MPS-MPS backbone
integrated type 4 = MPS pendant
```

This header is not compiled by itself. Compile either formulation wrapper:

```bash
g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    V22/v22_generator.cpp -o V22/v22_generator

g++ -std=c++17 -O2 -Wall -Wextra -Wpedantic \
    V35/v35_generator.cpp -o V35/v35_generator
```

Both wrappers use `special_bonds lj 0 0 0.5` for no-oil and oil-containing
systems.
