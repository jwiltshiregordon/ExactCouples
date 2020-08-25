
chainModule = method()
chainModule(Ring, ChainComplex) := Module => (Q, X) -> (
    expectChainRing(Q);
    R := ring X;
    supp := (min X)..(max X);
    t := local t;
    external := externalDegreeIndices Q;
    F := R[t, Degrees=>{(degree Q_0)_external}];
    XX := prune X;
    sh := (max X) * (degree Q_0);
    map(Q, F, {Q_0}) ** sequenceModule(F, reverse toList apply(supp, k -> XX.dd_k)) ** Q^{sh}
    )

chainModule(ChainComplex) := Module => X -> (
    R := ring X;
    d := getSymbol "d";
    A := R[d];
    Q := A / (A_0)^2;
    chainModule(Q, X)
    )

-- origin is always at zero
toChainComplex = method()
toChainComplex(Module) := ChainComplex => M -> (
    if not (isHomogeneous M) then (
        error "chain module must be homogeneous";
        );
    Q := ring M;
    expectChainRing Q;
    d := Q_0;
    degd := degree d;
    external := externalDegreeIndices Q;
    dm := oneEntry(0 * degd,, d);
    h := (heft Q)_external;
    degHefts := apply(degrees M, deg -> sum apply(deg_external, h, (i,j)->i*j));
    maxHeft := max degHefts;
    minHeft := min degHefts;
    dHeft := sum apply(degd_external, h, (i,j)->i*j);
    divDegs := select(minHeft..(maxHeft+dHeft), i -> (i % dHeft) == 0);
    ch := new ChainComplex; ch.ring = Q;
    diffs := reverse apply(divDegs, i -> evaluateInDegree({i}, M ** dm));
    i := local i;
    ldeg := -(last divDegs) // dHeft;
    ch#ldeg = target first diffs;
    for i from ldeg + 1 to ldeg - 1 + #diffs do (
        j := i - ldeg - 1;
        ch#i = source diffs#j;
        ch.dd#i = diffs#j;
        );
    ch
    )

chainModuleHomology = method()
chainModuleHomology(Module) := Module => C -> (
    R := ring C;
    expectChainRing R;
    d := R_0;
    HZ := coker matrix {{d}};
    Hom(R^{degree d}, Ext^1(HZ,C))
    )

chainModuleHomology(ZZ, Module) := Module => (k, C) -> (
    evaluateInDegree({k}, chainModuleHomology(C))
    )

sequenceModule = method()
sequenceModule(Ring, List) := (Q, L) -> (
    l := #L;
    if l == 0 then error "emtpy list of maps; for a single module, use extensionInDegree";
    external := externalDegreeIndices Q;
    dt := (degree Q_0)_external;
    i := local i;
    rows := for i from 0 to l list (
        M := if i < l then source L#i else target last L;
        extensionInDegree(i * dt, Q, M)
        );
    j := local j;
    cols := for j from 1 to l list extensionInDegree(j * dt, Q, source L#(j-1));
    ents := for i from 0 to l list (
        for j from 0 to l - 1 list (
            if i == j then cover(Q_0 * id_(rows#i))
            else if i == j + 1 then map(Q, ring L#j) ** -cover(L#j)
            else cover map(rows#i, cols#j, {})
            )
        );
    tar := directSum({Q^{}} | rows);
    src := directSum({Q^{}} | cols);
    -- why does matrix only work on maps between free modules?
    if not isHomogeneous map(tar, src, matrix ents) then (
        error "sequenceModule failed because the input was not homogeneous";
        );
    coker map(tar, src, matrix ents)
    )

sequenceModule(List) := L -> (
    if #L == 0 then error "emtpy list of maps; for a single module, use extensionInDegree";
    R := ring first L;
    t := getSymbol "t";
    Q := R[t];
    sequenceModule(Q, L)
    )


-- filtration will come from action of Q_0
-- so Q must have at least one generator
-- Start in degree zero, and proceed by degree Q_0
-- TODO: omit output ring
filtrationModule = method()
filtrationModule(Ring, List) := (Q, L) -> (
    l := #L;
    expectSequenceRing Q;
    dz := 0 * (degree Q_0);
    expectFiltrationList L;
    -- TODO: is l == 0 check already in expectFiltrationList L?
    if l == 0 then error "filtration must have at least one module";
    if l == 1 then return extensionInDegree(dz, Q, first L);
    sequenceModule(Q, apply(l - 1, i -> inducedMap(L#(i+1), L#i)))
    );
