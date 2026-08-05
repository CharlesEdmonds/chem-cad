# ChemCAD

A ChemDraw-style chemical structure sketcher with automatic bond angles, a
periodic table element picker, structure/name conversion in both directions, and
a reaction planner that suggests reagents, multi-step routes and predicted side
products.

C++20 + Dear ImGui (docking) + GLFW/OpenGL 3, with RDKit doing the chemistry.
The UI layer is plain ImGui with no framework of its own, so it can later be
hosted inside another ImGui application.

## Build

RDKit is not packaged by Fedora, so it is built once into a user-local prefix.

```bash
sudo dnf install -y boost-devel eigen3-devel freetype-devel libcurl-devel \
                    cmake ninja-build gcc-c++ java-latest-openjdk-headless
./scripts/setup_deps.sh          # builds RDKit + downloads the OPSIN jar (~20 min, cached)
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
./build/chemcad
```

`scripts/setup_deps.sh` is idempotent: it installs RDKit to
`~/.local/share/chemcad-deps/rdkit` and OPSIN to
`~/.local/share/chemcad-deps/share/opsin/opsin.jar`, and skips work already done.
Override the location with `CHEMCAD_DEPS_PREFIX`.

## Tests

```bash
ctest --test-dir build --output-on-failure
```

Every suite is hermetic — no network required. `test_ui_interaction` drives the
real canvas and panel code through a null ImGui backend, so sketching gestures
are verified without a display.

## Layout

| Area | What it does |
| --- | --- |
| Tool column (left) | Select, eraser, bond, chain, ring template, atom, charge tools, plus bond-order / stereo / ring pickers |
| Sketch tab | The drawing canvas |
| Reaction Planner tab | Starting-material boxes → product box, route suggestions |
| Periodic Table (right) | All 118 elements; click one to draw with it |
| Properties (right) | Formula, MW, cLogP, rings, canonical SMILES, IUPAC name, name→structure |

### Sketching

Bond tool: click empty canvas for a new fragment, click an atom to sprout the
next bond at the correct angle, drag for a snapped bond, and drag onto an
existing atom to close a ring. Click a bond to cycle single → double → triple.

Shortcuts: `C N O S P F B I` retype the hovered atom, `L`=Cl, `R`=Br; `1/2/3` set
the hovered bond order; `W` wedge, `Shift+W` hash; `Del` deletes the selection;
`Esc` returns to Select; `Ctrl+Z` / `Ctrl+Shift+Z` undo/redo; `Ctrl+C` / `Ctrl+V`
copy and paste SMILES. Wheel zooms, middle-drag (or `Space`+drag) pans.
**Structure → Clean Up Structure** relays the whole sketch with CoordGen.

## Naming

`name → structure` uses OPSIN locally (offline, systematic names) and falls back
to PubChem for trivial names like *aspirin*. `structure → name` uses PubChem,
since no good offline IUPAC namer exists. Both directions are cached in
`~/.cache/chemcad/namecache.json`; without a network the sketcher is unaffected
and naming simply reports that it is offline.

## Reaction planner

Fill one or more starting-material boxes and the product box (by SMILES, by name
lookup, or from the current sketch), then **Suggest Routes**. A breadth-first
search applies 82 curated reaction templates
(`data/reactions/*.json`) and reports each step's reagents, conditions,
**predicted side products**, and a mechanism note.

Routes found from the knowledge base are badged `KB`. Setting an API key enables
an `AI` fallback lane for targets the templates cannot reach:

```bash
export CHEMCAD_LLM_API_KEY=sk-...
export CHEMCAD_LLM_BASE_URL=https://api.openai.com   # optional
export CHEMCAD_LLM_MODEL=gpt-4o-mini                 # optional
```

Every LLM-proposed structure is validated through RDKit before it is shown, and
any LLM failure silently degrades to knowledge-base-only results.

### Adding reactions

Drop another object into any `data/reactions/*.json` array:

```json
{
  "id": "fischer_esterification",
  "name": "Fischer esterification",
  "smarts": "[C:1](=[O:2])[OX2H1].[OX2H1:3][C:4]>>[C:1](=[O:2])[O:3][C:4]",
  "arity": 2,
  "reagents": ["H2SO4 (cat.)"],
  "conditions": "reflux",
  "byproducts": ["O"],
  "priority": 5,
  "notes": "Equilibrium; driven by excess alcohol or removal of water.",
  "tags": ["esterification", "carbonyl"]
}
```

`byproducts` are the co-products the atom mapping cannot express — this is what
the planner surfaces as side products. `test_kb` validates every template
against RDKit, so a malformed SMARTS fails the build's test run rather than
silently never matching.

## Files

```
src/core/     molecule graph, undo, bond-angle geometry (STL only)
src/chem/     the single RDKit boundary
src/naming/   OPSIN subprocess + PubChem client + cache
src/rxn/      reaction knowledge base, route search, LLM client
src/ui/       ImGui panels: canvas, periodic table, properties, planner
src/app/      entry point, worker pool, project files, PNG capture
data/         periodic table + reaction knowledge base
```
