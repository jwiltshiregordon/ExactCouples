
-- flat thing goes first
tensorFlat = method()
tensorFlat(RingMap, Module) := Module => (phi, M) -> (
    law := x -> phi ** x;
    applyLawToModule(law, M)
    )

tensorFlat(RingMap, Matrix) := Matrix => (phi, f) -> (
    src := tensorFlat(phi, source f);
    tar := tensorFlat(phi, target f);
    map(tar, src, phi ** (cover f))
    )

tensorFlat(Module, Module) := Module => (F, M) -> (
    law := x -> F ** x;
    applyLawToModule(law, M)
    )


-- Returns a matrix over the coefficient ring of the ring of m
-- degreeLaw takes a multidegree to a list of multidegrees
-- entryLaw( 1 x 1 matrix ) returns an appropriate entry matrix
-- degreeLaw operates on the true degree, not just the external degree.

-- The output ring must be included as an argument.
-- Only apply this function on matrices that give a map of free modules.
-- TODO: allow non-free modules to be returned
applyEntrywise = method()
applyEntrywise(Ring, FunctionClosure, FunctionClosure, Matrix) := Matrix => (S, degreeLaw, entryLaw, m) -> (
    if not ((isFreeModule source m) and (isFreeModule target m)) then (
        error "applyEntrywise only operates on maps of free modules";
        );
    rows := toList(0..(-1 + numgens target m));
    cols := toList(0..(-1 + numgens source m));
    rowdegs := flatten apply(degrees target m, d -> degreeLaw(d));
    coldegs := flatten apply(degrees source m, d -> degreeLaw(d));
    if #rows == 0 or #cols == 0 then (
        return map(S^(-rowdegs), S^(-coldegs), {});
        );
    rawEntries := matrix for r in rows list for c in cols list entryLaw(submatrix(m, {r}, {c}));
    map(S^(-rowdegs), S^(-coldegs), sub(rawEntries, S))
    )

applyLawToModule = method()
applyLawToModule(FunctionClosure, Module) := Module => (law, M) -> (
    if M.?generators then (
        if M.?relations then (
            subquotient(law(gens M), law(relations M))
            ) else (
            image(law(gens M))
            )
        ) else (
        if M.?relations then (
            coker(law(relations M))
            ) else (
            ambient(target law(gens M))
            )
        )
    )


applyLawToMatrix = method()
-- Use this version when applyLawToModule has already been computed on the source and
-- target of f.  Pass these values as tar and src.
applyLawToMatrix(FunctionClosure, Module, Matrix, Module) := Matrix => (law, tar, f, src) -> (
    map(tar, src, law(cover f))
    )

applyLawToMatrix(FunctionClosure, Matrix) := Matrix => (law, f) -> (
    tar := applyLawToModule(law, target f);
    src := applyLawToModule(law, source f);
    applyLawToMatrix(law, tar, f, src)
    )

applyLawToChainComplex = method()
applyLawToChainComplex(Ring, FunctionClosure, ChainComplex) := ChainComplex => (Q, law, C) -> (
    i := local i;
    ret := new ChainComplex;
    for i from min C to max C do (
        ret#i = applyLawToModule(law, C_i);
        );
    for i from 1 + min C to max C do (
        ret.dd#i = applyLawToMatrix(law, ret_(i-1), C.dd_i, ret_i);
        );
    ret.ring = Q;
    ret
    )
