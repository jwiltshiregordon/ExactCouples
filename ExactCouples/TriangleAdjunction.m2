mapToTriangleRing = method()
mapToTriangleRing(Ring) := RingMap => S -> (
    expectChainSequenceRing S;
    R := coefficientRing S;
    degd := degree S_0;
    degf := degree S_1;
    d := baseName S_0;
    e := getSymbol "e";
    f := baseName S_1;
    external := externalDegreeIndices S;
    internal := internalDegreeIndices S;

    degd' := 2 * degd_external;
    degf' := 2 * degf_external;
    dege' := (degd - degf)_external;

    T := triangleRing(R, d, e, f, Degrees => {degd', dege', degf'});
    degMap := deg -> (2 * (deg_external)) | deg_internal;
    map(T, S, DegreeMap => degMap)
    )

distinguishedTriangleLaw = method()
distinguishedTriangleLaw(Ring) := FunctionClosure => R -> (
    expectTriangleRing R;
    (d, e, f) := (R_0, R_1, R_2);
    external := externalDegreeIndices R;
    (degd, dege, degf) := (d, e, f) / (v -> (degree v)_external);
    evenDegrees := image matrix transpose {degd, degf};
    spanDegrees := image matrix transpose {degd, dege, degf};
    ged := gens evenDegrees;
    projDeg := deg -> ((matrix transpose {deg_external}) // ged);
    -- we use projDeg to alternate signs in a consistent way,
    -- even away from spanDegrees.
    sgnDeg := deg -> (
        v := projDeg(deg);
	      nds := v_(0,0);
        if nds % 2 == 0 then 1 else -1
        );
    degreeLaw := x -> {x - 2*(degree e)};
    entryLaw := m -> (
        src := first degrees source m;
        tar := first degrees target m;
        sign := sgnDeg(src_external);
        dsub := -d + sign * e^2 * f;
        map(R^(degreeLaw tar), R^(degreeLaw src), sub(m, {d => dsub}))
        );
    m -> applyEntrywise(R, degreeLaw, entryLaw, m)
    )

distinguishedTriangle = method()
distinguishedTriangle(Ring, Module) := Module => (Q, M) -> (
    expectTriangleRing Q;
    S := ring M;
    expectChainSequenceRing S;
    external := externalDegreeIndices S;
    internal := internalDegreeIndices S;
    degMap := deg -> 2 * (deg_external) | deg_internal;
    title := "chain sequence ring " | (toString S)
             | "not compatible with triangle ring " | (toString Q) | ": ";
    if coefficientRing Q =!= coefficientRing S then (
        error(title | "coefficient rings do not match");
        );
    if degMap(degree S_0) != degree Q_0 then (
        error(title | "external degrees must double");
        );
    if degMap(degree S_1) != degree Q_2 then (
        error(title | "external degrees must double");
        );
    law := distinguishedTriangleLaw Q;
    sq := map(Q, S, {Q_0, Q_2}, DegreeMap => degMap);
    applyLawToModule(law, tensorFlat(sq, M))
    )

distinguishedTriangle(Module) := Module => M -> (
    S := ring M;
    expectChainSequenceRing S;
    phi := mapToTriangleRing S;
    Q := target phi;
    distinguishedTriangle(Q, M)
    )
