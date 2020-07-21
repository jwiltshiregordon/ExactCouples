
-- Q should be R[t,d]/d^2
-- ring M should be R[d]/d^2
-- TODO: omit output ring
canonicalFiltration = method()
canonicalFiltration(Ring, Module) := Module => (Q, M) -> (
    expectChainSequenceRing Q;
    expectChainRing(ring M);
    if coefficientRing Q =!= coefficientRing ring M then (
        error "coefficient rings do not match";
        );
    degMap := deg -> deg_{0} | deg;
    phi := map(Q, ring M, {Q_0 * Q_1}, DegreeMap => degMap);
    phi ** M
    )


forceAction = method()
forceAction(Matrix, Matrix) := Module => (f, m) -> (
    coker matrix {{f ** (id_(source m))}, {-(id_(source f)) ** m}}
    )
