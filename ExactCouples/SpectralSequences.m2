
-- E should be an exact couple
-- TODO: let user supply output ring
pageModule = method()
pageModule(ZZ, Symbol, Module) := Module => (r, dSymbol, E) -> (
    R := coefficientRing ring E;
    R[dSymbol_r];
    pageModule(r, value((baseName(dSymbol_r))#0), E)
    );
pageModule(ZZ, IndexedVariableTable, Module) := Module => (r, dVariable, E) -> (
    Q := ring E;
    expectCoupleRing Q;
    curpage := Q.Page;
    M := derivedCouple(r - curpage, E);
    Q = ring M; -- switch to this ring
    (e,f) := (Q_0,Q_1);


    external := externalDegreeIndices Q;
    internal := internalDegreeIndices Q;

    R := coefficientRing Q;
    DCover := R[dVariable_r, Degrees => {(degree e)_external}];
    D := DCover / (DCover_0)^2;
    d := D_0;

    degreeLaw := deg -> (
        if Q.isEvenDegree(deg_external) then (
            hd := ((deg_external) // 2) | deg_internal;
            D^{-hd}
            ) else (
            shiftDeg := deg + (degree e);
            hsd := ((shiftDeg_external) // 2) | shiftDeg_internal;
            coker oneEntry(hsd,,d)
            )
        );
    entryLaw := m -> (
        row := first degrees target m;
        col := first degrees source m;
        ent := m_(0,0);
        nent := if Q.isEvenDegree(row_external) then (
            if Q.isEvenDegree(col_external) then (
                ent
                ) else (
                e*ent
                )
            ) else (
            if Q.isEvenDegree(col_external) then (
                ent // e
                ) else (
                ent
                )
            );
        phi := map(D, Q, {0, 0});
        d * phi(nent // (e^2)) + phi(nent % (e^2))
        );
    pres := presentation M;
    rows := toList(0..<numgens target pres);
    cols := toList(0..<numgens source pres);
    --  D^{} appears to avoid empty direct sum
    tar := directSum ({D^{}} | apply(degrees target pres, d -> degreeLaw(d)));
    src := directSum ({D^{}} | apply(degrees source pres, d -> degreeLaw(d)));
    if #rows == 0 or #cols == 0 then (
        return map(tar, src, {});
        );
    rawEntries := matrix for r in rows list for c in cols list entryLaw(submatrix(pres, {r}, {c}));
    coker map(tar, src, rawEntries)
    )


-- ranges is (ps,qs,rs) all three are tuples
-- but rs can be a single integer.
-- E should be an exact couple
plotPages = method()
plotPages(Sequence, Function, Module) := Nothing => (ranges, dispfunc, E) -> (
    Q := ring E;
    expectCoupleRing Q;
    external := externalDegreeIndices Q;
    page := Q.Page;
    ps := toList ranges_0;
    qs := toList ranges_1;
    rs := if instance(ranges_2, ZZ) then {ranges_2} else toList ranges_2;
    local C;
    for r in rs do (
        dr := r - page;
        hdf := (degree Q_1) // 2;
        diffDeg := ((degree Q_0) - dr * hdf)_external;
        print("page " | r | ", with differential of degree " | toString(diffDeg) | ":");
        C = derivedCouple(dr, E);
        ents := reverse table(qs,ps,(q,p)->dispfunc({2*p,2*q},C));
        entsp := ents | {apply(ps,p->"p="|toString(p))};
        entspq := transpose({reverse(apply(qs,q->"q="|toString(q))) | {""}} |
                  {apply(qs,q->"") | {""}} |
                  transpose(entsp));
        print(netList entspq);
        print "";
        );
    )
