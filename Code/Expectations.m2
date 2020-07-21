expectSequenceRing = method()
expectSequenceRing(Ring) := Nothing => Q -> (
    title := "expectSequenceRing " | (toString Q) | ": ";
    if numgens Q != 1 then (
        error(title | "expected exactly one generator over coefficient ring");
        );
    )

expectChainRing = method()
expectChainRing(Ring) := Nothing => Q -> (
    title := "expectChainRing " | (toString Q) | ": ";
    if numgens Q != 1 then (
        error(title | "expected exactly one generator over coefficient ring");
        );
    R := coefficientRing Q;
    d := Q_0;
    if d * d != 0 then (
        error "generator must square to zero";
        );
    )

expectTriangleRing = method()
expectTriangleRing(Ring) := Nothing => Q -> (
    title := "expectTriangleRing " | (toString Q) | " =?= R[d,e,f]/(d^2, e^3): ";
    if numgens Q != 3 then (
        error "expected exactly three generators over coefficient ring";
        );
    (d, e, f) := (Q_0, Q_1, Q_2);
    if degree d != degree (e^2*f) then (
	      error(title | "monomials d and e^2*f must have the same degree");
	      );
    if d * d != 0 then (
        error(title | "variable d must square to zero");
        );
    if e^3 != 0 then (
        error(title | "variable e must cube to zero");
        );
    )

expectChainSequenceRing = method()
expectChainSequenceRing(Ring) := Nothing => Q -> (
    title := "expectChainSequenceRing " | (toString Q) | " =?= R[d,f]/d^2: ";
    if numgens Q != 2 then (
        error(title | "must have exactly two generators over coefficient ring");
        );
    (d, f) := (Q_0, Q_1);
    if d * d != 0 then (
        error(title | "variable d must square to zero");
        );
    )

expectCoupleRing = method()
expectCoupleRing(Ring) := Nothing => Q -> (
    title := "expectCoupleRing " | (toString Q) | " =?= R[e_r, f_r]: ";
    if numgens Q != 2 then (
        error(title | "must have exactly two generators over coefficient ring");
        );
    (e, f) := (Q_0, Q_1);
    baseNames := baseName \ (e, f);
    indexed := all(baseNames, s -> instance(s, IndexedVariable));
    if not indexed then (
        error(title | "generators must be indexed variables");
        );
    (re, rf) := last \ baseNames;
    if re != rf then (
        error(title | "generators must have same index");
        );
    r := re;
    Q.Page = r; -- store the page number in the couple ring
    (de, df) := degree \ (e, f);
    -- if other grading groups become available, this part should change
    -- the user should supply an element that doubles to (degree f).
    if df % 2 != 0 * df then (
        error(title | "degree of f must be divisible by two");
        );
    external := externalDegreeIndices Q;
    evenDegrees := image matrix transpose {2*de_external, df_external};

    gsd := matrix transpose {de_external, df_external};
    projDeg := deg -> ((matrix transpose {deg_external}) // gsd);
    deglen := #external;
    Q.isEvenDegree = deg -> (
        if #deg != deglen then (
            error("deg has wrong length; should be " | toString(deglen));
            );
        imDeg := image matrix (gsd * projDeg(deg));
	      isSubset(imDeg, evenDegrees)
        );
    Q.isOddDegree = deg -> (
        if #deg != deglen then (
            error("deg has wrong length; should be " | toString(deglen));
            );
        not Q.isEvenDegree deg
        );
    );


expectFiltrationList = method()
expectFiltrationList(List) := Nothing => L -> (
    if #L == 0 then (
        error "a filtration list must contain at least one module";
        );
    if not all(#L - 1, q -> isSubset(L#q, L#(q+1))) then (
        error "expected a list of submodules, each contained in the next";
        );
    );
