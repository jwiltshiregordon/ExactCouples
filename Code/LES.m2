
longExactSequence = method()
longExactSequence(Matrix) := Module => m -> (
    C := ring m;
    expectChainRing C;
    external := externalDegreeIndices C;
    degd := (degree C_0)_external;
    R := coefficientRing C;
    d := local d;
    f := local f;
    F := R[d, f, Degrees => {{0}|degd,{-1}|degd}] / d^2;
    --F := R[d, f, Degrees => {{0,1},{-1,1}}]/ideal(d^2);
    --phi := map(F,C,{d},DegreeMap => (deg -> {0} | deg));
    --zeds := toList((degreeLength R):0);
    --evencols := forceAction(oneEntry({1,0}|zeds,,F_1),phi**m);
    sm := sequenceModule((ring m)[f,Degrees=>{{-1}}],{m});
    sm' := map(F, ring sm, {F_1, F_0}, DegreeMap => deg -> deg) ** sm;
    --Q := coupleRing(R, 1, getSymbol "e", getSymbol "f");
    exactCouple(sm' ** F^{degree(F_1)}) -- banana
    )

longExactSequence(Ring, Matrix) := Module => (Q, m) -> (
    expectCoupleRing Q;
    M := longExactSequence(m);
    toQ := map(Q, ring M, {Q_0, Q_1});
    toQ ** M
    )

longExactSequenceToChainComplex = method()
longExactSequenceToChainComplex(List, ZZ, Module) := ChainComplex => (startDeg, rotations, M) -> (
    Q := ring M;
    expectCoupleRing Q;
    (e, f) := (Q_0, Q_1);
    (de, df) := degree \ (e, f);
    zeds := toList((degreeLength coefficientRing Q): 0);
    sd := startDeg | zeds;
    diffs := {map(Q^{sd}, Q^{sd - df}, {{f}})};
    for i from 1 to rotations do (
        diffs = diffs | {map(Q^{sd + de}, Q^{sd}, {{e}}),
                         map(Q^{sd + 2*de}, Q^{sd + de}, {{e}}),
                         map(Q^{sd + 2*de + df}, Q^{sd + 2*de}, {{f}})};
        sd = sd + 2*de + df;
        );
    external := externalDegreeIndices Q;
    evaluateInDegree(0 * external, chainComplex(reverse apply(diffs, ddd -> ddd ** M)))
    )

arrowAbove = method()
arrowAbove(Matrix) := Net => m -> (
    n := net m;
    w := width n;
    r := max(3, 2 + floor(w / 2));
    ar := concatenate(r: " -") | "> "; -- #ar == (r + 1) * 2
    sp := concatenate(floor(1 + r - w/2): " ");
    ar || (sp | n)
    )

arrowBelow = method()
arrowBelow(Matrix) := Net => m -> (
    n := net m;
    w := width n;
    r := max(3, 2 + floor(w / 2));
    ar := concatenate(r: " -") | "> ";
    sp := concatenate(floor(1 + r - w/2): " ");
    (sp | n)^(1 + depth n) || ar
    )

excerptCouple = method()
excerptCouple(List, ZZ, Module) := Net => (startDegree, rotations, M) -> (
    Q := ring M;
    expectCoupleRing Q;
    external := externalDegreeIndices Q;
    assert(#startDegree == #external);
    if not Q.isEvenDegree(startDegree) then (
        error("Start degree must be even");
        );
    -- longExactSequenceToChainComplex is not recommended
    ch := prune longExactSequenceToChainComplex(startDegree - (degree Q_0)_external, rotations, M);
    leftArt := " .- -> " || "(      " || " \\     ";
    rightArt := ("    \\ " || " - -' ")^1;
    o := {{(net ch.dd_1)^(-1), leftArt, ch_0, , , , , }};
    o = o | for i from 0 to rotations - 1 list (
        {(net ch.dd_(4 + 3 * i))^(-1), leftArt, ch_(3 + 3 * i),
          arrowAbove(ch.dd_(3 + 3 * i)), ch_(2 + 3 * i),
          arrowBelow(ch.dd_(2 + 3 * i)), ch_(1 + 3 * i), rightArt}
        );
    netList(o | {{ , , , , , , ch_(1 + 3 * (rotations + 0*1)), rightArt}},
            Boxes => false, Alignment => Center, VerticalSpace => 1)
    )

excerptLES = method()
-- Here, we assume that M is an exact couple produced by longExactSequence
excerptLES(ZZ, ZZ, Module) := Net => (k, l, M) -> excerptCouple({0,2*k-2},l-k+1,M);

excerptLES(ZZ, Module) := Net => (k, M) -> excerptLES(k, k, M)
