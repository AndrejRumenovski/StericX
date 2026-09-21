"""Finish independent model evidence; retain frozen SUT observations unchanged.

Uses the corrected primary-table extraction from models_reactions.py. All
reference fits enumerate Gini cutpoints independently of the production studies.
"""
from models_reactions import *


def main():
    rows = pd.read_json(M / 'inputs/crosscoupling_preprint_rows_v2.json')
    old = pd.read_json(M / 'inputs/crosscoupling_preprint_rows.json')
    old_keys = set(zip(old.reaction, old.id))
    added = [r for r in rows.to_dict('records') if (r['reaction'], r['id']) not in old_keys]
    save(M / 'results/extraction_correction.json', {
        'prior_rows': len(old), 'corrected_rows': len(rows), 'restored_rows': added,
        'cause': 'A numeric line can contain ID plus values when the ligand name wraps above or below it. The old audit parser incorrectly required a name on the numeric line.',
        'verification': 'Rendered original SI pages S90, S98, S99 and S100 inspected. Counts now match all twelve frozen study dataset sizes.',
        'retained_previous_results': 'results/prior_extraction_v1/',
    })
    # Transcribed from the primary preprint SI Tables S7 (Ni) and S8 (Pd),
    # visually checked in rendered original pages S85 and S86.
    published = {
        'I': [32.42, 'Left', .79, .62], 'II': [32.74, 'Left', .70, .53],
        'III': [31.55, 'Left', .67, .50], 'IV': [31.89, 'Left', .64, .45],
        'V': [51.53, 'Left', .66, .36], 'RS1': [31.89, 'Left', .70, .54],
        'VII': [32.43, 'Left', .67, .43], 'VIII': [28.87, 'Right', .86, .75],
        'IX': [28.65, 'Right', .93, .84], 'X': [30.53, 'Right', .93, .86],
        'XI': [28.82, 'Right', .90, .79], 'XII': [29.58, 'Right', .68, .38],
    }
    freeze(M / 'inputs/crosscoupling_primary_summary_tables.json',
           (json.dumps({'source': 'crosscoupling_preprint_si.pdf, Tables S7/S8, pages S85/S86',
                        'columns': ['threshold', 'direction', 'accuracy', 'mcc'],
                        'rows': published}, indent=2) + '\n').encode())
    all_models = json.loads((M / 'results/crosscoupling_reference_models.json').read_text())
    comparisons = []
    for rx, p in published.items():
        models = {f"{m['descriptor']}:{m['label_convention']}": m for m in all_models if m['reaction'] == rx}
        comparisons.append({'reaction': rx, 'primary_summary': dict(zip(['threshold', 'direction', 'accuracy', 'mcc'], p)), 'fits': models})
    save(M / 'results/crosscoupling_primary_comparison.json', comparisons)
    cuts = {'I': 10, 'II': 10, 'III': 5, 'IV': 5, 'V': 20, 'RS1': 10, 'VII': 10, 'VIII': 10, 'IX': 10, 'X': 5, 'XI': 30, 'XII': 10}
    boundary = rows[rows.apply(lambda r: r['yield'] == cuts[r.reaction], axis=1)]
    boundary.to_csv(M / 'results/crosscoupling_boundary_rows.csv', index=False, float_format='%.17g')
    native = pd.read_csv(M / 'inputs/current_native_percent_buried_volume_min.csv', index_col='molecule_id')
    by_id = native.sut.to_dict()
    assert not set(rows.id) - set(by_id), 'Retain and explain all missing SUT geometries.'
    rows['native_percent_buried_volume_min'] = rows.id.map(by_id)
    fidelities = []
    for name, subset in [('Ni', rows[rows.reaction.isin(list(cuts)[:6])]), ('Pd', rows[~rows.reaction.isin(list(cuts)[:6])])]:
        a = subset.vbur_min_published_rounded.to_numpy(); b = subset.native_percent_buried_volume_min.to_numpy()
        e = b-a; slope, intercept = np.polyfit(a,b,1)
        fidelities.append({'family':name, **metrics(a,b), 'median_absolute_error':float(np.median(abs(e))), 'maximum_absolute_error':float(max(abs(e))), 'slope':float(slope), 'intercept':float(intercept), 'note':'Repeated ligand IDs across reactions are retained. Primary values are printed to 0.1 percentage point; not full-precision reference geometries.'})
    save(M / 'results/crosscoupling_descriptor_fidelity.json', fidelities)
    rows.to_csv(M / 'results/crosscoupling_complete_inputs.csv', index=False, float_format='%.17g')

    # Independently reproduce leave-one-reaction-out transfer and explicitly
    # quantify chemical identity overlap, which is distinct from target leakage.
    ni = rows[rows.reaction.isin(list(cuts)[:6])].copy()
    transfers = []
    for convention in ['official_ge', 'stericx_study_gt']:
        ni['active'] = ni.apply(lambda r: int(r['yield'] >= cuts[r.reaction] if convention == 'official_ge' else r['yield'] > cuts[r.reaction]), axis=1)
        for held in list(cuts)[:6]:
            train = ni[ni.reaction != held]; test = ni[ni.reaction == held]
            fit = stump(train.native_percent_buried_volume_min.to_numpy(), train.active.to_numpy(), False)
            prediction = (test.native_percent_buried_volume_min.to_numpy() <= fit['threshold']).astype(int)
            if fit['direction'] == 'Right': prediction = 1-prediction
            assert fit['direction'] in ['Left','Right']
            overlap = sorted(set(train.id) & set(test.id))
            transfers.append({'held_reaction':held, 'label_convention':convention, 'n_test':len(test), 'n_train':len(train), 'threshold':fit['threshold'], 'direction':fit['direction'], 'accuracy':float(np.mean(prediction==test.active.to_numpy())), 'mcc':float(matthews_corrcoef(test.active.to_numpy(),prediction)), 'shared_ligand_ids':overlap, 'number_test_ligands_also_in_training':len(overlap), 'number_unique_test_ligands':int(len(set(test.id)-set(train.id))), 'interpretation':'Valid reaction-level withholding for this declared task; does not withhold all ligand identities or scaffold families.'})
    save(M / 'results/crosscoupling_transfer.json', transfers)

    # Check the updated official supporting information's retained target rows,
    # without claiming that the inaccessible correction letter was reviewed.
    correction_ee = {723:88,1057:79,1058:68,724:91,401:4,785:68,2063:18,2064:3,2067:32,2062:24,498:8}
    correction_yield = {723:63,1057:48,1058:16,724:95,401:69,785:23,2063:70,2064:61,2067:12,2062:43,498:35}
    src=pd.read_csv(M/'sources/ni_hda_reaction_data.csv',index_col=0)
    checked=[]
    for ident,ee in correction_ee.items():
        target=float(src.loc[ident,'ddG_abs']); expected=.00198720425864083*353.15*np.log((100+ee)/(100-ee))
        checked.append({'id':ident,'corrected_SI_ee_percent':ee,'corrected_SI_yield_percent':correction_yield[ident],'primary_csv_target_kcal_mol':target,'independent_target_at_353_15K':float(expected),'absolute_target_error':abs(target-expected)})
    freeze(M/'inputs/corrected_ni_hda_tableS3_transcription_v2.json',(json.dumps({'source':'correction_si.pdf, printed page S22 / PDF page 22, Table S3; only the eleven retained target IDs', 'rows':checked, 'amendment':'The initial transcription named the printed page incorrectly as 21. Numeric transcription is unchanged; the first artifact is retained.'},indent=2)+'\n').encode())
    save(M/'results/corrected_ni_hda_targets.json',{'rows':checked,'maximum_target_error':max(r['absolute_target_error'] for r in checked),'scope':'Retained ligand IDs and stored DDG targets reproduce the corrected SI table at its printed DDG precision. However ID2064 is internally inconsistent: SI and official CSV give ee=3%, while DDG=0.0281 corresponds to approximately 2% at353.15K. The correction letter was inaccessible (HTTP403); no assertion is made about all changes in that letter or original chromatograms.'})

    save(M/'manifest_finalize.json', {'purpose':'Post-analysis consolidation of model audit artifacts; this is not a claim that this manifest predated calculation. Immutable SUT observations are hashed, not rewritten.', 'files':{str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(M.rglob('*')) if p.is_file() and p.name not in ['manifest_finalize.json']}, 'reference_scripts':{str(p.relative_to(A)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted((A/'scripts').glob('models_*.py'))}})


if __name__ == '__main__':
    main()
