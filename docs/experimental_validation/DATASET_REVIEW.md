# Primary-literature and dataset review

Reviewed 2026-11-24. Accessibility means actual retrieved material, not a paper's
mention of modeling. `source_receipts.json` binds the initial downloads (including
403 failures), and raw bytes remain in `.stericx/experimental_validation/sources/`.
Additional structure and conformer-data receipts are stored beside their inputs.

## Selected primary study: nickel-catalyzed homo-Diels–Alder selectivity

**Cadge, Lozano, Merriman, Oblad, Sigman, Reisman**, *A Data Science-Guided Approach
for the Development of Nickel-Catalyzed Homo-Diels–Alder Reactions*, JACS (2025),
[DOI 10.1021/jacs.5c09948](https://doi.org/10.1021/jacs.5c09948).
The [Caltech primary SI](https://authors.library.caltech.edu/records/ef9d0-y2c15)
was actually downloaded and Table S3 (S22) rendered/visually inspected. The authors'
[GitHub repository](https://github.com/SigmanGroup/Ni-Catalyzed-hDA), immutable commit
`6b5a406dbd5af32ad0bc7739f686ab982d3aabef`, supplies the modeling table/notebook.

1. **Reaction:** Ni(cod)2/monophosphine hDA of norbornadiene and methyl vinyl ketone.
2. **Measured response:** product ee by chiral GC-FID; yield separately by GC-FID
   with p-xylene internal standard (reported average of two runs). Table also gives
   an absolute apparent selectivity ΔΔG, which is derived rather than directly measured.
3. **Observations:** Table S3 has thirteen entries: twelve quantified ee results and
   one unquantified L9 result. Eleven of the twelve measured ligand IDs occur in the
   author CSV, whose 1,566 total rows are mostly unlabeled and never count as experiments.
4. **Ligands:** 723 ((S)-L6),1057(L7),1058(L8),724(L10),401(L11),785(L12),906(L13),
   2063,2064,2067,2062,498(L14); L9 has trace yield and ee not determined. The SI
   Figure S12 supplies structures; graph identities are linked by Kraken ID, not names alone.
5. **Structures:** original MolSSI Kraken DFT-library V2000 exports were retrieved
   for eleven ligands, **135 conformers**, with explicit atoms/bonds and 4-decimal
   coordinates. ID906 is absent from the author CSV and its original structure
   metadata endpoint returns HTTP500; no identity or structure is invented to fill it.
6. **Conditions:** Ni(cod)2 6.5 mol%, ligand 13 mol%, 1,2-DCE 0.1 M, 80°C/353.15 K,
   16 h; norbornadiene 0.1 mmol and MVK 2 equivalents, as reported in Table S3.
7. **Conformers/energies:** exported free-ligand geometries are available; their
   historical selection/energy identities are not established. All 135 valid
   `/conformers/data/{id}` responses have `data:null`. Initial guessed endpoint
   `/conformers/{id}` produced 404 and is preserved as an acquisition error, not
   evidence of missing chemistry. The valid endpoint confirms the energy gap.
8. **Missing/conflicting data:** one unquantified ee; one unavailable molecular
   identity/geometry (906); no experimentally quantified conformer populations;
   no replicate ee uncertainties; incomplete stereochemical encoding for some axes
   or trivalent P centers. ID2064 has 3% ee but 0.0281 kcal/mol, inconsistent at
   353.15 K. ID723's SI enantiopure-ligand yield (63%) differs from the author CSV's
   96% screening yield; yields are not silently joined as the same condition.
9. **Provenance:** SI table and original author CSV retained separately. Published
   geometries are not replaced by generated geometries. The primary protocol uses
   the measured ee and reported reaction temperature; original ΔΔG is a sensitivity.
10. **Usage:** public author repository/SI access does not itself establish an open
    license. No LICENSE was identified in the captured Ni-hDA tree. Original files
    are retained locally for reproducibility; this campaign does not confer
    redistribution rights or publish a dataset/release.

The article's [Zenodo record](https://doi.org/10.5281/zenodo.15338599) was inspected.
It contains ORCA mechanistic intermediate/transition-state outputs, XYZ, Hessian
and density files (plus PPh3 and reactants). Those are not a matched conformer-energy
table for the twelve chiral ligand-screen entries. They cannot supply invented
Boltzmann energies for the exported Kraken ensembles.

## Other candidate studies

| Paper and DOI | Experiment / response / count | Identities and structures | Conditions, conformers and energies | Missing data, provenance and usage decision |
|---|---|---|---|---|
| Reid & Sigman, *Holistic prediction of enantioselectivity in asymmetric catalysis*, Nature (2019), [10.1038/s41586-019-1384-z](https://doi.org/10.1038/s41586-019-1384-z) | BINOL phosphoric-acid-catalyzed nucleophilic additions to imines; enantiomeric ratios / ΔΔG. Primary article describes 367 modeling entries. | Catalysts, imines and nucleophiles vary. Primary article identifies descriptor tables in an accompanying spreadsheet and acquisition in SI. No complete structure-to-row export was independently retrieved here. | Conditions vary across literature reactions; a single temperature assumption would be invalid. Common conformer-energy linkage not verified. | [Author manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6641578/) readable through web tool; direct archive download returns403 and publisher redirect fails. Count is article-level, not an independently enumerated clean matrix. Multi-component and condition confounding requires a separate mapping campaign. Academic-use terms remain those of the source. |
| Zahrt et al., *Prediction of higher-selectivity catalysts by computer-driven workflow and machine learning*, Science (2019), [10.1126/science.aau5631](https://doi.org/10.1126/science.aau5631) | Chiral phosphoric-acid-catalyzed thiol addition to imines; **1,075** catalyst/substrate combinations, 43 catalysts ×25 products, each run twice (2,150 experiments); duplicate mean selectivity. | Article identifies 43 catalysts and 5×5 substrate combinations. Structure drawings do not by themselves establish a reusable stereochemically mapped conformer set. | Conditions and temperature must be read from original SI per experiment; not verified here. No complete conformer-energy mapping retrieved. | [Primary manuscript](https://pmc.ncbi.nlm.nih.gov/articles/PMC6417887/) searchable/readable via web source, direct retrieval403. Valuable larger future test, but no fabricated structures, temperatures or independent conformer rows. AAAS/source terms apply. |
| Paton/Fletcher collaborators, *Enantioselective Conjugate Addition Catalyzed by a Copper Phosphoramidite Complex: Computational and Experimental Exploration of Asymmetric Induction*, ACS Catalysis (2017), [10.1021/acscatal.7b01453](https://doi.org/10.1021/acscatal.7b01453) | Copper/phosphoramidite conjugate addition; enantioselectivity and quantitative structure–selectivity relationships. Observation count not independently enumerated from original SI here. | Phosphoramidite identities/structures appear in the study; a complete row-to-coordinate map was not obtained. | Reaction-specific conditions, conformer selection and energy availability remain unverified. | Relevant Paton primary research, not an admitted numerical dataset. Missing reusable measurements/geometry mapping prevents a trustworthy reproduction in this campaign; no synthetic labels substituted. ACS/source terms apply. |
| Gensch et al., *A Comprehensive Discovery Platform for Organophosphorus Ligands for Catalysis*, JACS (2022), [10.1021/jacs.1c09718](https://doi.org/10.1021/jacs.1c09718) | Kraken molecular descriptor library; **not itself a catalytic-selectivity experiment for every ligand**. Historical repository audit covers 1,566 IDs, 1,546 geometry-bearing ligands and 31,721 conformers. | SMILES, graph/SDF and conformer IDs available through original MolSSI service; selected primary-study subset freshly retrieved. | Free-ligand DFT ensembles; source convention and selection provenance distinct from reaction-bound species. Current per-conformer energy data unavailable in checked subset. | Used as structure/reference source, never as 1,566 measured reactions. Original SI/code, API receipt and historical energy limitations are cited in the existing Kraken audit. Public access does not waive source terms. |
| Newman-Stonebraker et al., *Univariate classification of phosphine ligation state and reactivity in cross-coupling catalysis*, Science (2021), [10.1126/science.abj4213](https://doi.org/10.1126/science.abj4213) | Ni/Pd cross-coupling yield/activity and ligation thresholds, 12 reaction panels. Existing Pd replay contains 267 ligand–reaction rows, not 267 independent ligands. | Phosphine IDs and Kraken geometries; yield tables in primary SI. Existing studies retain scope and threshold disagreements. | Multiple reaction-specific cutoffs/conditions; conformer min %Vbur, not direct measured conformer populations. | Primary [author-hosted paper](https://doyle.chem.ucla.edu/wp-content/uploads/2021/10/62.-Univariate-classification-of-phosphine-ligation-state-and-reactivity-in-cross-coupling-catalysis.pdf) located. Original SI is absent from current `data/external`; earlier comparisons are historical evidence, not a new held-out experiment. AAAS source terms apply. |
| *Bisphosphine ligand conformer selection to enhance descriptor database representation: improving statistical modelling outcomes*, Chemical Science (2025), [10.1039/D5SC04691B](https://doi.org/10.1039/D5SC04691B) | Pd Hayashi–Heck regioselectivity and aryl-carbonylation enantioselectivity. Retrieved modeling matrices contain **63 rows** (Heck) and **42 rows** (SIA); these counts describe files, not verified unique experimental replicates. | [Original repository](https://github.com/SigmanGroup/BisphosphineConformerSelection), tree 95f29e73a66847fe8d4284a5458ead0f19dc94a8; 7,367 tree entries, structures and modeling data. Retrieved matrices have ligand IDs/selected published predictors/ΔΔG. | Bound bisphosphine complexes and conformer selection require separate two-donor mapping, metal removal and energy verification. Primary [SI](https://www.rsc.org/suppdata/d5/sc/d5sc04691b/d5sc04691b1.pdf) located; not assumed to provide automatically interchangeable free-ligand geometry. | Strong follow-up candidate with accessible original assets; selected-feature tables are not leakage-free raw descriptor panels. Not substituted post hoc for the locked primary study. Repository/article licenses must be checked for redistribution. |

## Selection rationale

Ni-hDA supplies a direct original measured-selectivity table, a fixed reaction,
identifiable monophosphines and reproducible geometry exports compatible with the
existing authoritative Rust engine. It is small and chemically clustered; that
limitation is fundamental, not repaired by the much larger unlabeled library or by
treating conformers as experimental rows. Its primary hypothesis, endpoint, feature
comparison, grouping and uncertainty analysis were locked before the new predictions
in `STATISTICAL_PROTOCOL.md` / `protocol_lock.json`. Larger, differently coordinated
systems are useful independent next studies, not opportunities to select a favorable
endpoint after looking at new results.

## Published correction checked

The [published correction, DOI 10.1021/jacs.5c18763](https://pubs.acs.org/doi/10.1021/jacs.5c18763)
was posted 22 December 2025 (2026 issue). It corrects other reaction conditions,
product assignments and characterization and provides a revised SI. We retrieved
[the corrected SI](https://doi.org/10.1021/jacs.5c18763.s001) through Figshare article
30933433, file 60556465, and visually inspected page S22. Its full Table S3 text,
including notes, equals the original after whitespace normalization. Thus all selected
screen ee/yield entries and the 2064 inconsistency remain unchanged. The corrected
scope/EDBO+ products are not interchangeable endpoints for this study.

The corrected SI metadata declares CC BY-NC 4.0. Original publisher, API and author
repository materials retain their own terms; this campaign does not relicense them.
The corrected PDF, metadata, HTTP receipt, rendered table and comparison hashes are
retained in `.stericx/experimental_validation/sources/`. The late discovery of the
correction and unchanged endpoint check are disclosed; this is not a claim that all
original paper measurements or assignments remain valid.
