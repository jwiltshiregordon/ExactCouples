internalDegreeIndices = method()
internalDegreeIndices(Ring) := List => R -> (
    S := coefficientRing R;
    dlr := degreeLength R;
    dls := degreeLength S;
    toList((dlr - dls)..(-1 + dlr))
    )


externalDegreeIndices = method()
externalDegreeIndices(Ring) := List => R -> (
    S := coefficientRing R;
    dlr := degreeLength R;
    dls := degreeLength S;
    toList(0..(-1 + dlr - dls))
    )

declareGenerators = method()
declareGenerators(Ring, List) := Module => (R, genList) -> (
    h := hashTable genList;
    degs := values h;
    varnames := keys h;
    free := R^(-degs);
    apply(numgens free, i -> globalAssign(varnames_i, free_i));
    free
    );

cospan = method(Dispatch => Thing)
cospan Sequence := Module => rels -> (
    if #rels == 0 then (
        error "cospan must have at least one relation";
        );
    tar := class rels_0;
    coker map(tar, , matrix {rels})
    );
cospan Thing := Module => (rel) -> (
    tar := class rel;
    coker map(tar, , matrix {rel})
    );

oneEntry = method()
oneEntry(List, List, RingElement) := Matrix => (rowdeg, coldeg, x) -> (
    R := ring x;
    map(R^{-rowdeg}, R^{-coldeg}, {{x}})
    )

oneEntry(List, Nothing, RingElement) := Matrix => (rowdeg, null, x) -> oneEntry(rowdeg, rowdeg + degree x, x)
oneEntry(Nothing, List, RingElement) := Matrix => (null, coldeg, x) -> oneEntry(coldeg - degree x, coldeg, x)
oneEntry(ZZ, ZZ, RingElement) := Matrix => (rowdeg, coldeg, x) -> oneEntry(rowdeg*(degree x),coldeg*(degree x),x)
oneEntry(ZZ, Nothing, RingElement) := Matrix => (rowdeg, null, x) -> oneEntry(rowdeg*(degree x),,x)
oneEntry(Nothing, ZZ, RingElement) := Matrix => (null, coldeg, x) -> oneEntry(,coldeg*(degree x),x)
