triangleRing = method(Options => {Degrees => {{0,2},{1,0},{-2,2}}})
triangleRing(Ring, Symbol, Symbol, Symbol) := Ring => o -> (R, dvar, evar, fvar) -> (
    S := R[{dvar, evar, fvar}, Degrees => o.Degrees];
    (d, e, f) := (S_0, S_1, S_2);
    S / ideal(d^2, e^3)
    )

coupleRing = method(Options => {Degrees => {{1,0},{-2,2}}})
coupleRing(Ring, ZZ, Symbol, Symbol) := Ring => o -> (R, r, e, f) -> (
    Q := R[e_r, f_r, Degrees=>o.Degrees];
    expectCoupleRing Q; -- installs Q.Page, Q.isEvenDegree, and Q.isOddDegree
    Q
    )

-- TODO: if Q already has its derived couple ring, return it.  Otherwise, stash
derivedCoupleRing = method()
derivedCoupleRing(Ring) := Ring => Q -> (
    expectCoupleRing Q;
    bne := baseName Q_0;
    bnf := baseName Q_1;
    r := bne#1;
    se := bne#0;
    sf := bnf#0;
    R := coefficientRing Q;
    external := externalDegreeIndices Q;
    de := (degree Q_0)_external;
    df := (degree Q_1)_external;
    hdf := first entries sub(matrix ({df/2}), ZZ);
    ret := R[se_(r+1), sf_(r+1), Degrees => {de - hdf, df}];
    expectCoupleRing ret; --installs ret.Page, ret.isEvenDegree, ret.isOddDegree
    ret
    );
