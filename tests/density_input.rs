use steric_x::density::*;
fn declaration() -> DensityDeclaration {
    DensityDeclaration {
        field: DensityField::TotalElectronDensity,
        unit: DensityUnit::ElectronsPerBohrCubed,
        provenance: "analytic scalar fixture".into(),
    }
}
fn cube() -> String {
    "title\ntotal density\n1 -1 -2 -3 1\n2 1 0 0\n2 0 2 0\n2 0 0 3\n6 6 0.5 1 2\n0 1D-3 2 3 4 5 6 7\n".into()
}
#[test]
fn preserves_field_geometry_and_order() {
    let g = DensityGrid::from_cube_str("fixture", &cube(), declaration()).unwrap();
    let m = g.metadata();
    assert_eq!(m.dimensions, [2, 2, 2]);
    assert_eq!(m.vectors_bohr, [[1., 0., 0.], [0., 2., 0.], [0., 0., 3.]]);
    assert_eq!(m.origin_bohr, [-1., -2., -3.]);
    assert_eq!(
        m.origin_angstrom,
        [-BOHR_ANGSTROM, -2. * BOHR_ANGSTROM, -3. * BOHR_ANGSTROM]
    );
    assert_eq!(m.atoms[0].position_bohr, [0.5, 1., 2.]);
    assert_eq!(m.comments, ["title", "total density"]);
    assert_eq!(g.value([0, 0, 1]), Some(0.001));
    assert_eq!(g.value([1, 0, 0]), Some(4.));
    assert_eq!(g.value([1, 1, 1]), Some(7.));
    assert_eq!(g.value([2, 0, 0]), None);
}
#[test]
fn rejects_incompatible_fields_units_and_bad_values() {
    let c = cube();
    for changed in [
        c.replace("1 -1 -2 -3 1", "-1 -1 -2 -3 1"),
        c.replace("1 -1 -2 -3 1", "1 -1 -2 -3 2"),
        c.replace("2 1 0 0", "-2 1 0 0"),
        c.replace("2 1 0 0", "2 0 0 0"),
        c.replace("1D-3", "NaN"),
        c.replace("1D-3", "-0.01"),
        c.replace("6 7\n", "6\n"),
        format!("{c}8\n"),
        c.replace("6 6 0.5", "119 6 0.5"),
    ] {
        assert!(
            DensityGrid::from_cube_str("fixture", &changed, declaration()).is_err(),
            "accepted {changed}"
        );
    }
    let mut missing = declaration();
    missing.provenance.clear();
    assert!(DensityGrid::from_cube_str("fixture", &c, missing).is_err());
}
#[test]
fn oblique_grid_is_preserved_without_reinterpretation() {
    let g = DensityGrid::from_cube_str(
        "fixture",
        &cube().replace("2 1 0 0", "2 1 0.25 0"),
        declaration(),
    )
    .unwrap();
    assert_eq!(g.metadata().vectors_bohr[0], [1., 0.25, 0.]);
}
