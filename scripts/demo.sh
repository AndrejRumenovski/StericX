#!/usr/bin/env bash
# Offline rehearsal using checked-in inputs. No study artifacts are regenerated.
set -Eeuo pipefail

usage() {
    cat <<'EOF'
Usage: scripts/demo.sh [--pause] [--screen]

Show descriptors, compare two ligands, search the shipped database, and inspect
the compact Ni-hDA model. --screen adds the retrospective S001 fit/deck workflow.
--pause waits for Enter before each step when run in an interactive terminal.

Requires Bash and a release binary built with: cargo build --release
Every run writes a fresh .stericx/demo-* directory inside the repository.
No network, Python, quantum calculations, or package installation is used.
EOF
}

pause=false
screen=false
for arg in "$@"; do
    case "$arg" in
        --pause) pause=true ;;
        --screen) screen=true ;;
        -h|--help) usage; exit 0 ;;
        *) printf 'Unknown option: %s\n' "$arg" >&2; usage >&2; exit 2 ;;
    esac
done

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd -P)"
binary="$repo_root/target/release/stericx"
if [[ ! -x "$binary" ]]; then
    printf 'Release binary missing: %s\nBuild it first:\n  cd %q\n  cargo build --release\n' \
        "$binary" "$repo_root" >&2
    exit 1
fi
if $pause && [[ ! -t 0 ]]; then
    printf 'No interactive terminal; continuing without pauses.\n'
    pause=false
fi

required=(
    data/xyz/SIG-NIHDA-723_1db84e66.xyz
    data/xyz/SIG-NIHDA-724_0d5c5a07.xyz
    data/ligand_db/kraken_phosphines.csv
    docs/study_001/stericx_portable_model.json
)
if $screen; then
    required+=(data/reactions.sigpack data/reactions_raw.csv
        docs/examples/ni_hda_screening_split.csv
        docs/examples/ni_hda_screening_tested.csv)
fi
for input in "${required[@]}"; do
    if [[ ! -r "$repo_root/$input" ]]; then
        printf 'Required demo input is missing or unreadable: %s\n' "$repo_root/$input" >&2
        exit 1
    fi
done

mkdir -p -- "$repo_root/.stericx"
output_dir="$(mktemp -d "$repo_root/.stericx/demo-$(date -u +%Y%m%dT%H%M%SZ)-XXXXXX")"
step='preparation'
trap 'printf "\nDemo failed during %s. Review files in: %s\n" "$step" "$output_dir" >&2' ERR
cd -- "$repo_root"

database="$repo_root/data/ligand_db/kraken_phosphines.csv"
saved_model="$repo_root/docs/study_001/stericx_portable_model.json"
ligand_723="$repo_root/data/xyz/SIG-NIHDA-723_1db84e66.xyz"
ligand_724="$repo_root/data/xyz/SIG-NIHDA-724_0d5c5a07.xyz"

run_step() {
    step="$1"
    local log="$2"
    shift 2
    printf '\n--- %s ---\n' "$step"
    if $pause; then
        read -r -p 'Press Enter to run this step: ' demo_pause_input
    fi
    {
        printf '\n# %s\n' "$step"
        printf '%q ' "$@"
        printf '\n'
    } >> "$output_dir/commands.txt"
    "$@" 2>&1 | tee "$output_dir/$log"
}

cat <<EOF
StericX offline research demo
Run folder: $output_dir
Guide: $repo_root/docs/DEMO.md

Frame this as descriptor tooling and an auditable research workflow.
The strongest benchmark is for max_delta_qvbur_min on matched DFT geometries
(R-squared 0.9852); it does not establish accuracy for every descriptor.
Sterimol B1 is approximate. Geometric similarity does not establish reactivity.
The compact Ni-hDA model has near-zero validation Q-squared.
The source response is ddG_abs: magnitude, without an R/S assignment.
No prospective measured validation is demonstrated here.
EOF

run_step 'Binary version' 00-version.txt "$binary" --version
run_step '1. Featurize ligand 723 (one representative conformer)' 01-descriptors.txt \
    "$binary" descriptors "$ligand_723"
run_step '2. Compare representative geometries for ligands 723 and 724' 02-compare.txt \
    "$binary" compare "$ligand_723" "$ligand_724" \
    --database "$database"
run_step '3. Search 1,541 precomputed ligands by geometric constraints' 03-search.txt \
    "$binary" search --database "$database" --vbur 30:35 --b5-max 8 --top 5
run_step '4. Check the saved model document' 04-model-validation.txt \
    "$binary" model validate "$saved_model"
printf '\nDocument validity checks structure; predictive quality requires the next diagnostics.\n'
printf 'Read ddG_abs as selectivity magnitude only; see docs/SCIENTIFIC_NOTES.md for response provenance.\n'
run_step '5. Inspect the compact Ni-hDA model and its weak validation' 05-model-inspect.txt \
    "$binary" model inspect "$saved_model"

if $screen; then
    cat <<'EOF'

Optional retrospective screening demonstration
Study 011 top-1 recovery: 0.158 versus 0.333 random; 10 of 57 fits failed.
Training uses ensemble Sterimol descriptors; this CSV screening path recomputes
one representative conformer. The resulting mismatch is a recorded limitation.
This S001 example was predefined, not selected for favorable performance.
Its group-LOO Q-squared is about -192.56: it is unsuitable for directing experiments.
The source CSV includes known responses. Screen ignores those response columns;
this remains a retrospective demonstration, and the exported deck is illustrative.
EOF
    run_step '6. Fit the predefined eight-train / three-held-out S001 split' 06-fit.txt \
        "$binary" fit \
        --data "$repo_root/data/reactions.sigpack" \
        --metadata "$repo_root/docs/examples/ni_hda_screening_split.csv" \
        --output "$output_dir/fit-report.json" \
        --predictions "$output_dir/frozen-predictions.csv" \
        --portable-model "$output_dir/model.json" \
        --model-id ni-hda-professor-demo-s001 \
        --reaction-family 'Ni-catalyzed homo-Diels-Alder' \
        --catalyst-metal Ni \
        --ligand-class 'monodentate phosphorus(III)' \
        --source-url 'https://raw.githubusercontent.com/SigmanGroup/Ni-Catalyzed-hDA/main/data/kraken.csv' \
        --response-temp-k 298.15 --optimize maximize \
        --response-sign-convention 'Magnitude |ddG| from ddG_abs; larger values mean greater enantioselectivity; no R/S assignment.' \
        --bootstrap 1000 --permutations 500 --seed 20260904
    run_step '7. Validate the newly fitted model document' 07-fit-validation.txt \
        "$binary" model validate "$output_dir/model.json"
    run_step '8. Inspect the scaffold-validation failure before screening' 08-fit-inspect.txt \
        "$binary" model inspect "$output_dir/model.json"
    run_step '9. Export an illustrative two-ligand deck after excluding training ligands' 09-screen.txt \
        "$binary" screen "$output_dir/model.json" \
        --library "$repo_root/data/reactions_raw.csv" \
        --exclude-tested "$repo_root/docs/examples/ni_hda_screening_tested.csv" \
        --top 2 --diverse --diversity-weight 0.5 \
        --export-deck "$output_dir/illustrative-candidate-deck.csv"
    printf '\nThe deck demonstrates export and traceability. This model has not earned experimental deployment.\n'
fi

printf '\nDemo completed. Commands and output logs: %s\n' "$output_dir"
