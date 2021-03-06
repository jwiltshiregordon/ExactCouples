-- This file, LawUtility.m2, is about exact functors that take R-modules
-- to S-modules.  Such a functor is determined by its operation on free
-- S-modules.

-- If F : (R-mod) --> (S-mod) is exact, then let "law" be a FunctionClosure
-- giving the operation of F on maps of free modules.  Specifically, given
-- a Matrix m describing a map A --> B between frees, let law(m) be a Matrix
-- describing FA --> FB (here, FA and FB need not be free modules).

-- With such a law in hand, it is easy to apply F to a general M2 module:

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

-- As an application, we can tensor with flat modules without computing
-- a gb.  Of course, it is up to the programmer to be sure of flatness.

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

-- Since any exact functor F : (R-mod) --> (S-mod) as above is additive,
-- its law is actually determined by its operation on 1x1 matrices that
-- give maps between free modules of rank one.  This restricted function
-- is called "entryLaw".

--The following function computes law from entryLaw.

--generateLaw = method()
--generateLaw(Ring, Ring, FunctionClosure) := FunctionClosure => (R, S, entryLaw) -> (
--    m -> ( -- m should be a map of free modules
--        if not ((isFreeModule source m) and (isFreeModule target m)) then (
--            error "law only operates on maps of free modules";
--            );
--        rows := toList(0..<numgens target m);
--        cols := toList(0..<numgens source m);
--        degreeLaw := deg -> target entryLaw(id_(R^({-deg})));
--        src := directSum({S^{}} | apply(degrees source m, degreeLaw));
--        tar := directSum({S^{}} | apply(degrees target m, degreeLaw));
--        if #rows == 0 or #cols == 0 then (
--            return map(tar, src, {});
--            );
--        -- TODO: add parallelization here
--        rawEntries := matrix for r in rows list for c in cols list entryLaw(submatrix(m, {r}, {c}));
--        map(tar, src, sub(rawEntries, S))
--        )
--    )

-- TODO: faster potentially to keep values on generators and compute products.

-- Returns a matrix over the coefficient ring of the ring of m
-- degreeLaw takes a multidegree to a list of multidegrees
-- entryLaw( 1 x 1 matrix ) returns an appropriate entry matrix
-- degreeLaw operates on the true degree, not just the external degree.

-- The output ring must be included as an argument.
-- Only apply this function on matrices that give a map of free modules.
applyEntrywise = method()
applyEntrywise(Ring, FunctionClosure, FunctionClosure, Matrix) := Matrix => (S, degreeLaw, entryLaw, m) -> (
    if not ((isFreeModule source m) and (isFreeModule target m)) then (
        error "applyEntrywise only operates on maps of free modules";
        );
    rows := toList(0..(-1 + numgens target m));
    cols := toList(0..(-1 + numgens source m));
    rowbatches := apply(degrees target m, d -> degreeLaw(d));
    colbatches := apply(degrees source m, d -> degreeLaw(d));
    rowdegs := flatten rowbatches;
    coldegs := flatten colbatches;
    if #rows == 0 or #cols == 0 then (
        return map(S^(-rowdegs), S^(-coldegs), {});
        );
    ret := mutableMatrix(S, #rowdegs, #coldegs, Dense=>true);  -- already tried false
    rowwidths := apply(rowbatches, e->#e);
    colwidths := apply(colbatches, e->#e);
    --i := 0;
    --j := 0; -- keep track of current upper left entry of block
    --for r in rows do (
    --    for c in cols do (
    --       if m_(r,c) != 0 then (
    --            blk := entryLaw(submatrix(m, {r}, {c}));
    --            for i' from 0 to rowwidths#r - 1 do (
    --                for j' from 0 to colwidths#c - 1 do (
    --                    ret_(i + i', j + j') = sub(blk_(i',j'),S);
    --                    );
    --                );
    --            );
    --        j = j + colwidths#c;
    --        );
    --    j = 0;
    --    i = i + rowwidths#r;
    --    );
    rawEntries := matrix for r in rows list for c in cols list 
        if m_(r,c) != 0 then entryLaw(submatrix(m, {r}, {c})) else map(S^(-rowbatches#r),S^(-colbatches#c),{});
    map(S^(-rowdegs), S^(-coldegs), sub(rawEntries, S))
    --map(S^(-rowdegs), S^(-coldegs), matrix ret)
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
