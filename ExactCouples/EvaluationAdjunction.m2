-- The adjunction for evaluation in a particular degree.  As the right adjoint, we have
-- r : (internal degrees) --> (full degrees)
--          ideg          -->   d | ideg
--
-- The left adjoint takes the form
-- l : (full degrees) --> (internal degrees)
--     (edeg | ideg)  --> ideg + ... + ideg, where the number of summands is the number of
--                                           external monomials of degree (d - edeg).
-- We will solve for l as a functor by using naturality of the adjunction unit.
-- The unit is a natural transformation 1 ==> rl.  In this case, given a full degree x,
-- the component at x will be a single row whose degree is x, and columns for every external
-- monomial of degree (d - x_external).

evaluateInDegreeLaw = method()
evaluateInDegreeLaw(List, Ring) := FunctionClosure => (d, R) -> (
    internal := internalDegreeIndices R;
    external := externalDegreeIndices R;
    mR := presentation module R;
    fastBasis := deg -> (
        degz := deg | (0 * internal);
        map(R^1,,rawBasis(raw mR, degz, degz, heft R, 0..(-1 + numgens R), false, -1))
        );
    unitLaw := x -> fastBasis(d - x_external) ** R^{-x};
    unitLaw = memoize unitLaw;
    degreeLaw := x -> apply(degrees source unitLaw x, deg -> deg_internal);
    entryLaw := m -> m * unitLaw((first degrees source m)) // unitLaw((first degrees target m));
    m -> applyEntrywise(coefficientRing R, degreeLaw, entryLaw, m)
    --generateLaw(R, coefficientRing R, entryLaw)
    )

evaluateInDegree = method()
evaluateInDegree(List, Module) := Module => (d, M) -> (
    external := externalDegreeIndices ring M;
    if #d != #external then (
        error "d must be an external degree for ring M";
        );
    law := evaluateInDegreeLaw(d, ring M);
    applyLawToModule(law, M)
    )

evaluateInDegree(List, Matrix) := Matrix => (d, f) -> (
    external := externalDegreeIndices ring f;
    if #d != #external then (
        error "d must be an external degree for ring f";
        );
    law := evaluateInDegreeLaw(d, ring f);
    applyLawToMatrix(law, f)
    )

-- Convenience methods for evaluating at a ring element
structureMap = method()
structureMap(List, List, RingElement, Module) := Matrix => (rowdeg, coldeg, x, M) -> (
    external := externalDegreeIndices ring M;
    if #rowdeg != #external then (
        error "rowdeg must be an external degree for ring M";
        );
    if #coldeg != #external then (
        error "coldeg must be an external degree for ring M";
        );
    zed := 0 * (internalDegreeIndices ring M);
    evaluateInDegree(0*rowdeg, Hom(oneEntry(rowdeg | zed, coldeg | zed, sub(x,ring M)), M))
    )

structureMap(List, Nothing, RingElement, Module) := Matrix => (rowdeg, null, x, M) -> (
    zed := 0 * (internalDegreeIndices ring M);
    evaluateInDegree(0*rowdeg, Hom(oneEntry(rowdeg | zed, , sub(x,ring M)), M))
    )

structureMap(Nothing, List, RingElement, Module) := Matrix => (null, coldeg, x, M) -> (
    zed := 0 * (internalDegreeIndices ring M);
    evaluateInDegree(0*coldeg, Hom(oneEntry(, coldeg | zed, sub(x,ring M)), M))
    )



evaluateInDegree(List, ChainComplex) := ChainComplex => (d, C) -> (
    law := evaluateInDegreeLaw(d, ring C);
    applyLawToChainComplex(coefficientRing ring C, law, C)
    )



extensionInDegreeLaw = method()
extensionInDegreeLaw(List, Ring) := FunctionClosure => (d, R) -> (
    degreeLaw := x -> {d | x};
    entryLaw := m -> map(R^(-degreeLaw(first degrees target m)),
                   R^(-degreeLaw(first degrees source m)), sub(m, R));
    m -> applyEntrywise(R, degreeLaw, entryLaw, m)
    )

-- Given a ring, a multidegree, and a module for its coefficient ring,
-- build the left kan extension along the multidegree
extensionInDegree = method()
extensionInDegree(List, Ring, Module) := Module => (d, E, M) -> (
    law := extensionInDegreeLaw(d, E);
    applyLawToModule(law, M)
    )

extensionInDegree(List, Ring, Matrix) := Module => (d, E, f) -> (
    law := extensionInDegreeLaw(d, E);
    applyLawToMatrix(law, f)
    )

