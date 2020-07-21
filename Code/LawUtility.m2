-- If R and S are rings, and if F is an exact functor
-- F : (R-mod) --> (S-mod)
-- with the property that F sends free R-modules to free S-modules,
-- then it is easy to apply F to a general M2 module if we have the 
-- following two functions:
-- 
-- degreeLaw : (R-degrees) --> (lists of S-degrees)
-- 
--             (1x1 matrices over R  )     (matrices over S)
-- entryLaw  : (encoding maps between) --> (encoding maps  )
--             (rank-1 free modules  )     (of free modules)
-- 
-- Here, degreeLaw says what happens to free modules:
--
-- F(R^{-deg}) = S^(degreeLaw deg)
--
-- and entryLaw says what happens to a 1x1 matrix m
--
-- F(m) = entryLaw(m)
-- 
-- We require that entryLaw be compatible with degreeLaw so that
-- the source and target of entryLaw(m) match degreeLaw applied
-- to the source and target.  For this reason, it is possible to
-- infer degreeLaw from R and entryLaw:

--inferDegreeLaw = method()
--inferDegreeLaw(Ring, FunctionClosure) := FunctionClosure => (R, entryLaw) -> (
--    deg -> degrees target entryLaw(id_(R^{-deg}))
--    )


-- The following code actually works for any exact F with F(free) = f.g.
-- We may compute F(Macaulay2 module) by the following 
-- easy code (which uses only entryLaw since degreeLaw is determined):
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
-- a gb.  This is much faster, but it is up to the programmer to be
-- sure that the module or map is actually flat!

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

-- TODO: general functoriality for adjointable laws
-- TODO: complete adjunction from degreeLaws, one entryLaw, and the unit or counit