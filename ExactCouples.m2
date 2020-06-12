-- -*- coding: utf-8 -*-
-- installPackage("ExactCouples",FileName => "/Users/jwiltshiregordon/Dropbox/Programming/Macaulay2/ExactCouples.m2")
newPackage(
  "ExactCouples",
      Version => "0.5",
      Date => "May 20, 2020",
      Authors => {
       {
       Name => "John D. Wiltshire-Gordon",
       Email => "jwiltshiregordon@gmail.com"}
       },
      HomePage => "https://sites.google.com/wisc.edu/jwg/home",
      Headline => "spectral sequences by Massey's method of exact couples",
  AuxiliaryFiles => false, -- set to true if package comes with auxiliary files
      DebuggingMode => true,     -- set to true only during development
  PackageImports => {"Elimination"}
      )

-- Any symbols or functions that the user is to have access to
-- must be placed in one of the following two lists
export {"applyEntrywise",
  "declareGenerators", "cospan",
  "internalDegreeIndices", "externalDegreeIndices",
  "evaluateInDegree",
  "extensionInDegree",
  "evaluateInDegreeLaw", "extensionInDegreeLaw", "distinguishedTriangleLaw",
  "expectChainRing", "expectCoupleRing", "expectTriangleRing",
  "chainModuleHomology",
  "toChainComplex",
  "triangleRing", "distinguishedTriangle",
  "longExactSequence", "excerptLES", --"longExactSequenceToChainComplex",
  "restackRing", "restackModule",
  --"arrowAbove", "arrowBelow",
  --"tensorFlat",
  "oneEntry",
  "exactCouple", "expectExactCouple", "derivedCouple", "pageModule",
  "coupleRing", "declareCouple",
  "isEvenDegree", "isOddDegree", "Page",
  "plotPages",
  "derivedCoupleRing", "enforceCoupleRelations", "excerptCouple",
  "sequenceModule", "filtrationModule",
  "canonicalFiltration", "expectFiltrationList",
  "contravariantExtCouple", "contravariantExtLES", "covariantExtCouple", "covariantExtLES",
  "TorCouple", "TorLES",
  "chainModule", "mapToTriangleRing"}
exportMutable {}

-- we need to use rawBasis
debug Core

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
    )

evaluateInDegree = method()
evaluateInDegree(List, Module) := Module => (d, M) -> (
    law := evaluateInDegreeLaw(d, ring M);
    applyLawToModule(law, M)
    )

evaluateInDegree(List, Matrix) := Matrix => (d, f) -> (
    law := evaluateInDegreeLaw(d, ring f);
    applyLawToMatrix(law, f)
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


-- Set up rings that encode various common setups in homological algebra

-- Ring builders

triangleRing = method(Options => {Degrees => {{0,2},{1,0},{-2,2}}})
triangleRing(Ring, Symbol, Symbol, Symbol) := Ring => o -> (R, dvar, evar, fvar) -> (
    S := R[{dvar, evar, fvar}, Degrees => o.Degrees];
    (d, e, f) := (S_0, S_1, S_2);
    S / ideal(d^2, e^3)
    )

coupleRing = method(Options => {Degrees => {{1,0},{-2,2}}})
coupleRing(Ring, ZZ, Symbol, Symbol) := Ring => o -> (R, r, e, f) -> (
    Q := R[e_r, f_r, Degrees=>o.Degrees];
    expectCoupleRing Q; -- installs Q.Page, Q.isEvenDegree, and Q.isOddDegree
    Q
    )

-- expectos
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

-- Set up common conversions

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

forceAction = method()
forceAction(Matrix, Matrix) := Module => (f, m) -> (
    coker matrix {{f ** (id_(source m))}, {-(id_(source f)) ** m}}
    )

longExactSequence = method()
longExactSequence(Matrix) := Module => m -> (
    C := ring m;
    expectChainRing C;
    R := coefficientRing C;
    d := local d;
    f := local f;
    F := R[d, f, Degrees => {{0,1},{-1,1}}]/ideal(d^2);
    phi := map(F,C,{d},DegreeMap => (deg -> {0} | deg));
    zeds := toList((degreeLength R):0);
    evencols := forceAction(oneEntry({1,0}|zeds,,F_1),phi**m);
    g := getSymbol "e";
    h := getSymbol "f";
    Q := coupleRing(R, 1, g, h);
    exactCouple(Q, evencols)
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
excerptLES(ZZ, ZZ, Module) := Net => (k, l, M) -> excerptCouple({0,0},l-k+1,M);

excerptLES(ZZ, Module) := Net => (k, M) -> excerptLES(k, k, M)


restackRing = method()  -- this code can be tightened
restackRing(List, Ring) := RingMap => (l, R) -> (
    q := #l;
    S := R;
    {F, psi} := flattenRing S;
    FF := ambient F;
    tI := ideal presentation F;
    dls := {};  --record the degree lengths along the tower
    for i to q - 1 do (
        try coefficientRing S then (
            dls = {degreeLength S} | dls;
            S = coefficientRing S;
            ) else (
                error ("Input ring must be a tower of height at least #l = " | q);
            );
        );
    dls = {0} | dls;
    idls := reverse apply(-1 + #dls, i -> -dls_i + dls_(i+1));
    S = R;
    batches := {};
    for i to q - 1 do (
        external := toList(0..<(idls_i));
        batches = {{gens S, apply(degrees S, deg -> deg_external)}} | batches;
        S = coefficientRing S;
        );
    rels := {};
    for i to q - 1 do (
        ovars := flatten (apply(i,j->batches_j_0) | apply(toList((i+1)..<q),j->batches_j_0));
        ovars =  apply(ovars, v -> (
            psi := map(FF, ring v);
            psi(v)
            ));
        rels = rels | {eliminate(ovars, tI)};
        );
    ll := apply(l, i -> i - 1);
    batches = batches_ll;
    rels = rels_ll;
    T := S;
    for i to q - 1 do (
        T = T[batches_i_0, Degrees => batches_i_1];
        psi := map(T,ring (rels_i));
        T = T / ideal(psi ** (gens rels_i));
        );
    ie := apply(q, i -> toList((dls_i)..<(dls_(i+1))));
    perm := inversePermutation flatten ie_ll;
    regrade := deg -> deg_perm;
    phi := map(R, T);
    I := ker phi;
    map((source phi) / I, R, DegreeMap => regrade)
    )

restackModule = method()
restackModule(List, Module) := Module => (l, M) -> tensorFlat(restackRing(l, ring M), M)


enforceCoupleRelations = method()
enforceCoupleRelations(Module) := Module => M -> (
    Q := ring M;
    expectCoupleRing Q;
    (e, f) := (Q_0, Q_1);
    pageRels := pdeg -> map(Q^{-pdeg},, {{e^3, f}});
    auxRels := adeg -> map(Q^{-adeg},, {{e^2, e*f}});
    external := externalDegreeIndices Q;
    rels := deg -> if Q.isEvenDegree(deg_external) then pageRels(deg) else auxRels(deg);
    P := presentation M;
    -- directSum fails if (degrees target P) == {}; that's the reason for id_(Q^{})
    coker(P | directSum({id_(Q^{})} | apply(degrees target P, rels)))
    );


-- TODO: I think this code ruins the variable R_1
exactCouple = method()
exactCouple(Module) := Module => (M) -> (
    R := ring M;
    expectChainSequenceRing R;
    external := externalDegreeIndices R;
    internal := internalDegreeIndices R;
    e := getSymbol "e";
    f := baseName R_1;
    S := coefficientRing R;
    dege := ((degree R_0) - (degree R_1))_external;
    degf := 2 * (degree R_1)_external;
    Q := S[e_1,f_1, Degrees => {dege, degf}]; -- TODO: use coupleRing
    exactCouple(Q, M)
    )

-- TODO: fix error messages
-- no way to omit Q since it contains page information
-- exactCouple(coupleRing(?), M) is pretty easy, though
-- TODO: do examples in documentation.
exactCouple(Ring, Module) := Module => (Q, M) -> (
    expectCoupleRing Q;
    R := ring M;
    expectChainSequenceRing R;
    T := target mapToTriangleRing R;
    if degree Q_0 != degree T_1 then (
	error "degree of first generator should be (degree d) - (degree f)";
	);

    if degree Q_1 != degree T_2 then (
	error "external degree of f should double";
	);
    C := distinguishedTriangle(prune M);
    T = ring C; -- update T to match C
    Hd := coker oneEntry(-degree T_0,,T_0);
    -- although the following ring map is not well-defined, the subsequent
    -- use of enforceCoupleRelations corrects the problem.
    phi := map(Q,T,{0,Q_0,Q_1});
    enforceCoupleRelations(phi ** Ext^1(Hd,C))
    )

expectExactCouple = method()
expectExactCouple(Module) := Nothing => M -> (
    Q := ring M;
    expectCoupleRing Q;
    e := Q_0;
    f := Q_1;
    dz := 0 * (degree e);
    external := externalDegreeIndices Q;
    gensEven := N -> all(degrees N, deg -> Q.isEvenDegree(deg_external));
    gensOdd := N ->  all(degrees N, deg -> Q.isOddDegree(deg_external));
    if not gensEven(prune image(M ** oneEntry(dz,,e^2))) then (
        error "e^2 fails to annihilate aux";
        );
    if not gensOdd(prune image(M ** oneEntry(dz,, f))) then (
        error "f fails to annihilate page";
        );
    kerfmodime := prune(kernel(M ** oneEntry(,dz,f))/image(M ** oneEntry(dz,, e)));
    if not gensEven(kerfmodime) then (
        error "failure of exactness at aux: ker f != im e.";
        );
    keremodimf := prune(kernel(M ** oneEntry(,dz,e))/image(M ** oneEntry(dz,,f)));
    if not gensEven(keremodimf) then (
        error "failure of exactness at aux: ker e != im f.";
        );
    -- plan to check exactness at page:
    -- build the maps M -e-> M -e-> M/e^2 and take homology;
    -- on the page, this gives A --> E --> A/e^2 = A
    e1 := M ** oneEntry(dz,,e);
    e2 := M ** oneEntry(,dz,e);
    e2' := map((target e2)/e^2,target e2,id_(target e2)) * e2;
    keremodime := prune((kernel e2')/(image e1));
    if not gensOdd(keremodime) then (
        error "failure of exactness at page: ker e != im e.";
        );
    );

-- TODO: if Q already has its derived couple ring, return it.  Otherwise, stash
derivedCoupleRing = method()
derivedCoupleRing(Ring) := Ring => Q -> (
    expectCoupleRing Q;
    bne := baseName Q_0;
    bnf := baseName Q_1;
    r := bne#1;
    se := bne#0;
    sf := bnf#0;
    R := coefficientRing Q;
    external := externalDegreeIndices Q;
    de := (degree Q_0)_external;
    df := (degree Q_1)_external;
    hdf := first entries sub(matrix ({df/2}), ZZ);
    ret := R[se_(r+1), sf_(r+1), Degrees => {de - hdf, df}];
    expectCoupleRing ret; --installs ret.Page, ret.isEvenDegree, ret.isOddDegree
    ret
    );

derivedCouple = method()
derivedCouple(Module) := Module => M -> (
    if M.cache.?derivedCouple then return M.cache.derivedCouple;
    Q := ring M;
    expectCoupleRing Q;
    external := externalDegreeIndices Q;
    internal := internalDegreeIndices Q;

    e := Q_0;
    f := Q_1;

    de := (degree e);
    df := (degree f);
    hdf := df // 2;
    dz := 0 * de;

    ff := f ** id_M;
    iso := (inducedMap(image ff, coimage ff, ff)) ** Q^{hdf};
    A' := source iso;
    imf := target iso;

    keree := kernel(oneEntry(,dz,e^2) ** id_M);
    pushEven := deg -> if Q.isEvenDegree(deg_external) then 1_Q else e;
    toEven := diagonalMatrix(Q, apply(degrees cover keree, pushEven));

    cocycles := image map(keree,, toEven);
    zinc := inducedMap(M, cocycles);
    coboundaries := image((e^2 ** id_M) // zinc);
    E' := cocycles / coboundaries;
    ztoe := inducedMap(E', cocycles);

    m := map(M,, oneEntry(-df,, e) ** id_(cover A'));
    A'toE' := prune map(E', A' ** Q^{-de+hdf}, ztoe * (m // zinc));
    n := (oneEntry(dz,,e) ** id_(ambient E')) * (gens E');
    nn := map(coker(relations imf), cover E', n);
    slv := (nn // inducedMap(target nn, imf));
    E'toA' := prune map(A',E' ** Q^{-de+hdf},iso^-1 * slv);

    Q' := derivedCoupleRing Q;

    eDrop := map(Q', Q, {0, Q'_1});
    E'Lan := eDrop ** (prune E');
    A'Lan := eDrop ** (prune A');

    QtoQ' := map(Q', Q, {Q'_0, Q'_1});

    eE'Lan := Q'_0 * id_E'Lan;
    eA'Lan := Q'_0 * id_A'Lan;

    A'toE'Lan := map(E'Lan,, QtoQ' ** A'toE');
    E'toA'Lan := map(A'Lan,, QtoQ' ** E'toA');

    ret := coker map(E'Lan ++ A'Lan,, matrix {{-eE'Lan,A'toE'Lan},{E'toA'Lan,-eA'Lan}});
    ret = prune ret;
    r := 1;

    M.cache.derivedCouple = ret;
    ret
    )

derivedCouple(ZZ, Module) := Module => (n, M) -> (
    if n < 0 then error "cannot form derived couples a negative number of times";
    if n == 0 then M else derivedCouple derivedCouple(n - 1, M)
    )

-- E should be an exact couple
-- Let user supply output ring?
-- TODO: current way is bad mix.
pageModule = method()
pageModule(ZZ, Symbol, Module) := Module => (r, dSymbol, E) -> (
    R := coefficientRing ring E;
    R[dSymbol_r];
    pageModule(r, value((baseName(dSymbol_r))#0), E)
    );
pageModule(ZZ, IndexedVariableTable, Module) := Module => (r, dVariable, E) -> (
    Q := ring E;
    expectCoupleRing Q;
    curpage := Q.Page;
    M := derivedCouple(r - curpage, E);
    Q = ring M; -- switch to this ring
    (e,f) := (Q_0,Q_1);


    external := externalDegreeIndices Q;
    internal := internalDegreeIndices Q;

    R := coefficientRing Q;
    DCover := R[dVariable_r, Degrees => {(degree e)_external}];
    D := DCover / (DCover_0)^2;
    d := D_0;

    degreeLaw := deg -> (
        if Q.isEvenDegree(deg_external) then (
            hd := ((deg_external) // 2) | deg_internal;
            D^{-hd}
            ) else (
            shiftDeg := deg + (degree e);
            hsd := ((shiftDeg_external) // 2) | shiftDeg_internal;
            coker oneEntry(hsd,,d)
            )
        );
    entryLaw := m -> (
        row := first degrees target m;
        col := first degrees source m;
        ent := m_(0,0);
        nent := if Q.isEvenDegree(row_external) then (
            if Q.isEvenDegree(col_external) then (
                ent
                ) else (
                e*ent
                )
            ) else (
            if Q.isEvenDegree(col_external) then (
                ent // e
                ) else (
                ent
                )
            );
        phi := map(D, Q, {0, 0});
        d * phi(nent // (e^2)) + phi(nent % (e^2))
        );
    pres := presentation M;
    rows := toList(0..<numgens target pres);
    cols := toList(0..<numgens source pres);
    --  D^{} appears to avoid empty direct sum
    tar := directSum ({D^{}} | apply(degrees target pres, d -> degreeLaw(d)));
    src := directSum ({D^{}} | apply(degrees source pres, d -> degreeLaw(d)));
    if #rows == 0 or #cols == 0 then (
        return map(tar, src, {});
        );
    rawEntries := matrix for r in rows list for c in cols list entryLaw(submatrix(pres, {r}, {c}));
    coker map(tar, src, rawEntries)
    )


-- ranges is (ps,qs,rs) all three are tuples
-- but rs can be a single integer.
-- E should be an exact couple
plotPages = method()
plotPages(Sequence, Function, Module) := Nothing => (ranges, dispfunc, E) -> (
    Q := ring E;
    expectCoupleRing Q;
    external := externalDegreeIndices Q;
    page := Q.Page;
    ps := toList ranges_0;
    qs := toList ranges_1;
    rs := if instance(ranges_2, ZZ) then {ranges_2} else toList ranges_2;
    local C;
    for r in rs do (
        dr := r - page;
        hdf := (degree Q_1) // 2;
        diffDeg := ((degree Q_0) - dr * hdf)_external;
        print("page " | r | ", with differential of degree " | toString(diffDeg) | ":");
        C = derivedCouple(dr, E);
        ents := reverse table(qs,ps,(q,p)->dispfunc({2*p,2*q},C));
        entsp := ents | {apply(ps,p->"p="|toString(p))};
        entspq := transpose({reverse(apply(qs,q->"q="|toString(q))) | {""}} | 
                  {apply(qs,q->"") | {""}} |
                  transpose(entsp));
        print(netList entspq);
        print "";
        );
    )

declareCouple = method()
declareCouple(Ring, List, List) := Module => (Q, pageGens, auxGens) -> (
    expectCoupleRing Q;
    tabA := hashTable auxGens;
    tabE := hashTable pageGens;
    internal := internalDegreeIndices Q;
    external := externalDegreeIndices Q;
    if not all(values tabE, deg -> Q.isEvenDegree(deg_external)) then (
        error "page generators must be even";
        );
    if not all(values tabA, deg -> Q.isOddDegree(deg_external)) then (
        error "aux generators must be odd";
        );
    degs := values tabA | values tabE;
    varnames := keys tabA | keys tabE;
    ret := enforceCoupleRelations Q^(-degs);
    apply(numgens ret, i -> globalAssign(varnames_i, ret_i));
    ret
    );


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

-- TODO: let user supply output ring
-- algorithm:
-- Left Kan extend Y to the sequence ring, and then right Kan to chain ring.
-- call the resulting thing Y'.  
-- Take the cofiltration of (submods | {0}), and resolve by frees.  
-- Call the resulting chain sequence module M.
-- claim: exactCouple Hom(M, Y') is the correct couple
-- Pf: Hom(M, Y') = Hom(M, Ran Lan Y) = Hom(Res M, Lan Y).
-- Since M is a resolution, Res M is free, say Res M = Lan G.
-- Hom(Lan G, Lan Y) = Hom(G, Res Lan Y).  And this thing is
-- the chain sequence module that defines the couple.
-- It would be mathematically similar to start with the cofiltration instead;
-- this could be faster at least some of the time.
contravariantExtCouple = method()
contravariantExtCouple(List, Module) := Module => (submods, Y) -> (
    R := ring last submods;
    expectFiltrationList submods;
    f := local f;
    F := R[f, Degrees=>{{-1}}];
    l := #submods;
    -- the next commented line leads to many correct entries, but converges to zero
    -- so it is incorrect.
    -- It could be a prototype for implementing edge homomorphisms, though.
    --submods' := submods | {(last submods)/(last submods)};
    submods' := {last submods} | apply(submods, m->(last submods)/m);
    fm := sequenceModule(F, apply(l, k -> inducedMap(submods'#(k+1), submods'#k)));
    -- only use res over flat rings due to M2 bug
    flattenModule := m -> ((flattenRing ring m)#1) ** m;
    fm = flattenModule fm;
    d := local d;
    ringfm := ring fm;
    ch := ringfm[d]/d^2;
    rfm := flattenModule chainModule(ch, res fm);
    S := R[d,f,Degrees=>{{1,0},{0,-1}}]/d^2;
    why := map(ring rfm, ring Y,DegreeMap=>(deg->{0,0}|deg));
    Y' := why ** Y;
    hm := Hom(rfm,Y');
    pres := presentation hm;
    -- Must remove certain rows and columns from pres
    -- because our chain sequence module is only correct in certain degrees;
    -- we had truncated to keep it fg
    dtp := degrees target pres;
    rowselect := select(#dtp,k->((dtp#k)#1)>0);
    dsp := degrees source pres;
    colselect := select(#dsp,k->((dsp#k)#1)>0);
    hm = coker(pres_colselect^rowselect);
    couple := exactCouple((map(S,ring hm)) ** hm);
    Q := ring couple;
    sh := Q^{4*degree(Q_0)+1*degree(Q_1)};
    couple ** sh
    )

TEST ///
    R = (ZZ/17)[x,y,z]
    randomFilteredModule = () -> (
        genmod := R^(apply(5,p->-random(3)));
        relmod := R^(apply(2,p->-random(4)));
        pres := random(genmod, relmod);
        X := coker pres;
        A0 := image map(X,,(id_genmod)_{0});
        A1 := image map(X,,(id_genmod)_{0,1,2});
        {A0,A1,X}
        );
    submods = randomFilteredModule();
    expectFiltrationList submods
    Y = prune(submods#1);
    couple = prune contravariantExtCouple(submods,Y);
    expectExactCouple couple
    flattenModule = m -> ((flattenRing ring m)#1) ** m;
    C = flattenModule couple;
    C' = flattenModule derivedCouple(5,couple);
    X = last submods;
    A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
    E1 = (p,q) -> Ext^p(A(q)/A(q-1),Y);
    A1 = (p,q) -> Ext^p(X/A(q-1),Y);
    proj = q -> inducedMap(X/A(q),X);
    filt = (p,q) -> image Ext^p(proj q,Y);
    Einfty = (p,q) -> prune(filt(p,q-1)/filt(p,q));
    entropy = 0;
    expectEqual = (x,y) -> (
        assert(x == y);
        if x != 0 then (
            entropy = entropy + size2(x);
            );
        );
    hf = hilbertFunction;
    for p from 0 to 2 do (
        for q from 0 to 3 do (
            page1 = E1(p,q);
            aux1 = A1(p,q);
            pageInfty = Einfty(p,q);
            for n from 0 to 20 do (
                expectEqual(hf({n},page1),     hf({2*p,   2*q,   n}, C));
                expectEqual(hf({n},aux1),      hf({2*p-1, 2*q-1, n}, C));
                expectEqual(hf({n},pageInfty), hf({2*p,   2*q,   n}, C'));
                );
            );
        );
    print("total assertion entropy " | (toString entropy));
///

contravariantExtLES = method()
contravariantExtLES(ZZ, Module, Module, Module) := Net => (k, X, A, Y) -> (
    expectFiltrationList {A,X};
    E1 := contravariantExtCouple({A,X},Y);
    excerptCouple({-2,0},k,E1)
    )

-- TODO: allow user to supply output ring
covariantExtCouple = method()
covariantExtCouple(Symbol, Module, Module) := Module => (eSymbol, W, seqmod) -> (
    R := ring W;
    F := ring seqmod;
    expectSequenceRing F;
    if not (coefficientRing ring seqmod) === ring W then (
        error("last two arguments incompatible: covariantExtCouple(eSymbol, W, seqmod) " + 
              "requires (coefficientRing ring seqmod) === ring W");
        );
    t := baseName(F_0);
    external := externalDegreeIndices F;
    degt := (degree F_0)_external;
    Q := coupleRing(R,1,eSymbol,t,Degrees=>{{1}|(-degt),{0}|(2*degt)});
    d := local d;
    Co := R[d]/d^2;
    {R', theta} := flattenRing R;
    cW := chainModule(Co, res(theta ** W));
    M := Hom(cW, Co^{degree(Co_0)});
    
    S := R[d,t,Degrees=>{{1}|(0*degt),{0}|degt}]/d^2;
    phi := map(S,Co,DegreeMap=>deg->{deg_0}|(0*degt)|(deg_{1..<#deg}));
    psi := map(S,F,DegreeMap=>deg->{0}|deg);
    C := (phi ** M) ** (psi ** seqmod);
    exactCouple(Q, C)
    )

covariantExtCouple(Module, List) := Module => (W, submods) -> (
    expectFiltrationList submods;
    R := ring last submods;
    F := R[getSymbol "f"];
    fm := filtrationModule(F, submods);
    covariantExtCouple(getSymbol "e", W, fm)
    )

TEST ///
    R = (ZZ/17)[x,y,z]
    randomFilteredModule = () -> (
        genmod := R^(apply(5,p->-random(3)));
        relmod := R^(apply(2,p->-random(4)));
        pres := random(genmod, relmod);
        X := coker pres;
        A0 := image map(X,,(id_genmod)_{0});
        A1 := image map(X,,(id_genmod)_{0,1,2});
        {A0,A1,X}
        );
    submods = randomFilteredModule();
    expectFiltrationList submods
    W = prune(submods#1);
    couple = prune covariantExtCouple(W,submods);
    expectExactCouple couple
    flattenModule = m -> ((flattenRing ring m)#1) ** m;
    C = flattenModule couple;
    C' = flattenModule derivedCouple(5,couple);
    X = last submods;
    A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
    E1 = (p,q) -> Ext^p(W,A(q)/A(q-1));
    A1 = (p,q) -> Ext^p(W,A(q));
    inc = q -> inducedMap(X,A(q));
    filt = (p,q) -> image Ext^p(W, inc q);
    Einfty = (p,q) -> prune(filt(p,q)/filt(p,q-1));
    entropy = 0;
    expectEqual = (x,y) -> (
        assert(x == y);
        if x != 0 then (
            entropy = entropy + size2(x);
            );
        );
    hf = hilbertFunction;
    for p from 0 to 2 do (
        for q from 0 to 3 do (
            page1 = E1(p,q);
            aux1 = A1(p,q);
            pageInfty = Einfty(p,q);
            for n from 0 to 20 do (
                expectEqual(hf({n},page1),     hf({2*p,   2*q,   n}, C));
                expectEqual(hf({n},aux1),      hf({2*p-1, 2*q+1, n}, C));
                expectEqual(hf({n},pageInfty), hf({2*p,   2*q,   n}, C'));
                );
            );
        );
    print("total assertion entropy " | (toString entropy));
///

covariantExtLES = method()
covariantExtLES(ZZ, Module, Module, Module) := Net => (k, W, X, A) -> (
    expectFiltrationList {A,X};
    E1 := covariantExtCouple(W,{A,X});
    excerptCouple({-2,2},k,E1)
    );

-- TODO: allow user to supply output ring
TorCouple = method()
TorCouple(Symbol, Module, Module) := Module => (eSymbol, W, seqmod) -> (
    R := ring W;
    F := ring seqmod;
    expectSequenceRing F;
    if not (coefficientRing ring seqmod) === ring W then (
        error("last two arguments incompatible: TorCouple(eSymbol, W, seqmod) " + 
              "requires (coefficientRing ring seqmod) === ring W");
        );
    t := baseName(F_0);
    external := externalDegreeIndices F;
    degt := (degree F_0)_external;
    Q := coupleRing(R,1,eSymbol,t,Degrees=>{{-1}|(-degt),{0}|(2*degt)});
    d := local d;
    Ch := R[d,Degrees=>{-1}]/d^2;
    
    -- work around res bug for ring towers
    {R', theta} := flattenRing R;
    rW := (theta^(-1))(res(theta ** W));
    M := chainModule(Ch, rW);
    
    S := R[d,t,Degrees=>{{-1}|(0*degt),{0}|degt}]/d^2;
    phi := map(S,Ch,DegreeMap=>deg->{deg_0}|(0*degt)|(deg_{1..<#deg}));
    psi := map(S,F,DegreeMap=>deg->{0}|deg);
    C := (phi ** M) ** (psi ** seqmod);
    exactCouple(Q, C)
    )
    
    
TorCouple(Module, List) := Module => (W, submods) -> (
    expectFiltrationList submods;
    R := ring last submods;
    F := R[getSymbol "f"];
    fm := filtrationModule(F, submods);
    TorCouple(getSymbol "e", W, fm)
    )

TEST ///
    R = (ZZ/17)[x,y,z]
    randomFilteredModule = () -> (
        genmod := R^(apply(5,p->-random(3)));
        relmod := R^(apply(2,p->-random(4)));
        pres := random(genmod, relmod);
        X := coker pres;
        A0 := image map(X,,(id_genmod)_{0});
        A1 := image map(X,,(id_genmod)_{0,1,2});
        {A0,A1,X}
        );
    submods = randomFilteredModule();
    expectFiltrationList submods
    W = prune(submods#1);
    couple = prune TorCouple(W,submods);
    expectExactCouple couple
    flattenModule = m -> ((flattenRing ring m)#1) ** m;
    C = flattenModule couple;
    C' = flattenModule derivedCouple(5,couple);
    X = last submods;
    A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
    E1 = (p,q) -> Tor_p(W,A(q)/A(q-1));
    A1 = (p,q) -> Tor_p(W,A(q));
    inc = q -> inducedMap(X,A(q));
    filt = (p,q) -> image Tor_p(W, inc q);
    Einfty = (p,q) -> prune(filt(p,q)/filt(p,q-1));
    entropy = 0;
    expectEqual = (x,y) -> (
        assert(x == y);
        if x != 0 then (
            entropy = entropy + size2(x);
            );
        );
    hf = hilbertFunction;
    for p from 0 to 2 do (
        for q from 0 to 3 do (
            page1 = E1(p,q);
            aux1 = A1(p,q);
            --pageInfty = Einfty(p,q); -- No direct check available because M2 is missing
                                       -- the method Tor(Module,Matrix).
            for n from 0 to 20 do (
                expectEqual(hf({n},page1),     hf({2*p,   2*q,   n}, C));
                expectEqual(hf({n},aux1),      hf({2*p+1, 2*q+1, n}, C));
                --expectEqual(hf({n},pageInfty), hf({2*p,   2*q,   n}, C'));
                );
            );
        );
    print("total assertion entropy " | (toString entropy));
///


TorLES = method()
TorLES(ZZ, Module, Module, Module) := Net => (k, W, X, A) -> (
    expectFiltrationList {A,X};
    E1 := TorCouple(W,{A,X});
    Q := ring E1;
    external := externalDegreeIndices Q;
    rotDeg := (2 * (degree Q_0) + (degree Q_1))_external;
    startDeg := {2,2} - (k - 2) * rotDeg;
    excerptCouple(startDeg,k,E1)
    );

beginDocumentation()

doc ///
    Key
        ExactCouples
    Headline
        spectral sequences by Massey's method of exact couples
    SeeAlso
        "Conventions and first examples"
        "Bockstein spectral sequence"
        "Serre spectral sequence in homology"
///


doc ///
    Key
        evaluateInDegree
        (evaluateInDegree,List,Module)
        (evaluateInDegree,List,Matrix)
        (evaluateInDegree,List,ChainComplex)
    Headline
        evaluates a module in a particular degree
    Usage
        evaluateInDegree(deg, M)
    Inputs
        M:Module
           over some ring R
        deg:List
            an external degree for R
    Outputs
        :Module
            the part of M sitting in degree deg, considered as a module for
            the coefficient ring of R
    Description
        Text
        Example
            S = QQ[s, t, u]; R = S[x, y]; m = matrix {{s*x^2+t*x*y+u*y^2}}
            N = evaluateInDegree({4}, coker matrix {{s*x^2+t*x*y+u*y^2}})
            apply(10, i -> hilbertFunction({i}, N))
            (F, f) = flattenRing R; apply(10, i -> hilbertFunction({4, i}, coker (f m)))
    SeeAlso
        extensionInDegree
///


 document {
     Key => {applyEntrywise, (applyEntrywise, Ring, FunctionClosure, FunctionClosure, Matrix)},
     Headline => "apply a matrix-valued function to the entries of a matrix",
     Usage => "applyEntrywise(S, degreeLaw, entryLaw, m)",
     Inputs => {
   "m" => Matrix => {"a homogeneous matrix over some ring R"},
   "S" => Ring => {"the output ring"},
   "degreeLaw" => FunctionClosure => {"a function that converts a multidegree for R ",
                                      "to a list of multidegrees for S"},
   "entryLaw" => FunctionClosure => {"a function that accepts a homogeneous one-by-one matrix over R and ",
                                     "returns a matrix over S whose degrees match those determined by degreeLaw.",
     "This function will be called on every one-by-one submatrix of m"}

        },
     Outputs => {"the block matrix obtained by replacing each entry e with the matrix entryLaw(e)"}
     }


doc ///
    Key
        chainModule
        (chainModule,Ring,ChainComplex)
        (chainModule,ChainComplex)
    Headline
        writes a chain complex of R-modules as an R[d]/d^2-module
    Usage
        M = chainModule(Q, X)
    Inputs
        X:ChainComplex
            bounded, and whose terms are modules for some ring R
        Q:Ring
            of the form R[d]/d^2
    Outputs
        M:Module
            a graded Q-module encoding the chain complex
    Description
        Text
            If $\alpha$ is the degree of the variable d, then M is supported
            in degrees that are multiples of $\alpha$.  The part of M sitting
            in degree $d \cdot \alpha$ matches X_{-d}.  If $\alpha = -1$, then we have
            X_d = M_d for all d.
        Example
            R = QQ[x, y, z]; M = coker vars R; C = res M -- a Koszul complex
            Q = R[d, Degrees => {-1}] / ideal(d^2); m = chainModule(Q, C)
            (F, f) = flattenRing Q;
            matrix table(10, 10, (i, j) -> hilbertFunction({j,i}, f ** m))
            matrix table(10, 10, (i, j) -> hilbertFunction(i, C_j))
    SeeAlso
        (toChainComplex,Module)
///


doc ///
    Key
        extensionInDegree
        (extensionInDegree,List,Ring,Module)
        (extensionInDegree,List,Ring,Matrix)
    Headline
        places a copy of a module in a certain degree
    Usage
        E = extensionInDegree(deg,Q,M)
    Inputs
        Q:Ring
        M:Module
            over the coefficient ring of Q
        deg:List
            a multidegree for the ring Q
    Outputs
        E:Module
            freely generated by a copy of M sitting in degree deg
///

doc ///
    Key
        distinguishedTriangle
        (distinguishedTriangle,Ring,Module)
        (distinguishedTriangle,Module)
    Headline
        constructs the mapping cone along with its inclusion and projection
    Usage
        distinguishedTriangle(Q, M)
    Inputs
        Q:Ring
            satisfying the criteria in expectTriangleRing; if omitted, a suitable ring is generated
        M:Module
            for a ring of the form R[d,f]/d^2, a subring of Q
    Outputs
        :Module
            over the ring Q, encoding a distinguished triangle
    Description
        Text
            There are a lot of details to record.  (The following is actually describing 
                distinguishedTriangleLaw, and I've since decided that it is too technical
                to be useful to ordinary users)

            Suppose $R$ is a ring of coefficients, and that the triangle ring Q takes
            the form $Q = R[d,e,f]/(d^2,e^3)$.  Write $Even$ for the subgroup generated
            by $(deg d)$ and $(deg f)$.  Since Q is a triangle ring, we know $(deg e)$ is
            not an element of $Even$, but $2 * (deg e) = (deg d) - (deg f)$, and so
            $2 * (deg e) \in Even$.

            We also know that $Even$ admits a group homomorphism $sgn : Even \to \{1,-1\}$ to
            the multiplicative unit group of $\ZZ$ so that


            $sgn(d) = -1$


            $sgn(f) = 1$.


            Write $All$ for the subgroup generated by $Even$ and the element $(deg e)$ so
            that $|All : Even| = 2$.

            The true grading group of $Q$ might be much more than just $Even$.  In this case,
            a general degree $\alpha$ can be written $\alpha = (\alpha // All) + (\alpha % All)$.
            In what follows, we assume all degrees are in $All$ by replacing $\alpha$
            by $(\alpha // All)$ if necessary.

            We define an additive functor $r : R[d,e,f]/(d^2,e^3) \to R[d,f]/d^2$ by the
            following formulas.


            Set $\delta, \epsilon, \phi$ to be the degrees of $d, e, f$.  A general even degree
            takes the form $\alpha = a * \delta + b * \phi$, and $\alpha + \epsilon$ is then odd.


            {\bf The right adjoint}


            On degrees:


            $r(\alpha) = \alpha \oplus (\alpha + 2\epsilon)$


            $r(\alpha + \epsilon) = \alpha + 2\epsilon$


            On entries with even row-degree {\tt alpha}, $d, e, f$ become

        CannedExample
            | d  0 |         |     0     |         | f  0 |
            | f -d |         | sgn alpha |         | 0  f |
        Text
            On entries with odd row-degree, $d, e, f$ become
        CannedExample
             | d |              | 1  0 |            | f |
        Text
            These matrices give the usual construction of the mapping cone and its related maps.
            It is easy to check that these six matrices commute when composable.  The most interesting
            case is
        CannedExample
            | d  0 | * |        0        |  ==  |     0     | * | d  0 |
            | f -d |   | sgn (d * alpha) |      | sgn alpha |   | f -d |
        Text
            which holds because $sgn (d * \alpha) = sgn(d) * sgn(\alpha) = -sgn(\alpha)$.
            It must also be checked that these matrices satisfy the relations $d^2=e^3=0$.

            {\bf Constructing the left adjoint}

            We want to find a left adjoint to $r$ so that we can compute $r^*$ as $l_!$.
            The key fact is that $r^*$ takes free modules to free modules.  Specifically,
            a free module with generator in degree $\alpha$ becomes a free module with
            generator in degree $\alpha - \epsilon$.  The new basis vector
            is encoded in the adjunction unit, whose component at an even degree $\alpha$
            takes the form
        CannedExample
            | 1 |
        Example
            2*2
    Caveat
    SeeAlso
        expectTriangleRing
///


doc ///
    Key
        longExactSequence
        (longExactSequence,Matrix)
        (longExactSequence,Ring,Matrix)
    Headline
        finds the long exact sequence associated to a map of R[d]/d^2-modules
    Usage
        M = longExactSequence m
    Inputs
        m:Matrix
            whose ring is of the form R[d]/d^2 for some coefficient ring R, and where d has degree +1
    Outputs
        LES:Module
            an exact couple encoding the sequence, which lives in bidegrees
            $\{-1,q\} \to \{0,q\} \to \{1,q} \to \{-1,q+1\} \to ...$
    Description
        Text
            Suppose $m:A \to B$, where A and B are considered cochain complexes by virtue
            of the square-zero action of d, and note that m commutes with this action since
            it is a map of modules.  Writing C(m) for the mapping cone of the map m, the long
            exact sequence in cohomology takes the form
            $\cdots \to H^p A \to H^p B \to H^p C(m) \to H^{p+1} A \to \cdots$
            The output module M encodes this sequence as follows:


            The cohomology of A appears in the first (zero-indexed) column: M_{1,2p} = $H^p A$;


            The cohomology of B appears in the minus-first column, with a shift: M_{-1,2p} = $H^{p-1} B$;


            The cohomology of C(m) appears in the even rows of the zeroth column;


            Moreover, setting e = (ring LES)_0 and f = (ring LES)_1, we have:

            Multiplication by e induces the horizontal maps $H^p B \to H^p C(m) \to H^{p+1}A$ and

            Multiplication by f induces diagonal maps $H^p A \to H^p B$ of degree (-2, 1).
            In fact, under the isomorphisms M_{2,p} = $H^p A$ and M_{0,p} = $H^{p-1} B$, this map coincides with $H^p m : H^p A \to H^p B$.

        Example
            R = QQ[x]; S = R[d] / ideal(d^2); declareGenerators(S, {a => {0,0}}); A = cospan(x^2*a, d*x*a)
            declareGenerators(S, {b => {0,0}}); B = cospan(x^2*b, d*b)
            m = map(B, A, matrix {b}); LES = longExactSequence m;
            excerptLES(0,2,LES)
    Caveat
        The source and target of m must be homogeneous, and m itself must be degree-preserving
    SeeAlso
        excerptLES
///


doc ///
    Key
        declareGenerators
        (declareGenerators,Ring,List)
    Headline
        builds a free module and names its generators
    Usage
        declareGenerators(R,genList)
    Inputs
        R:Ring
            the coefficients
        genList:List
            of the form \{ ... , x => deg, ... \} where x is a symbol and deg is a degree
    Outputs
        :Module
            free on the generating set of variables in their various degrees, and with coefficients drawn from R
    Consequences
        Item
          The symbols in the list variableNames bind to the generators of the output module.
    Description
        Text
            Useful for declaring many variables at once.
        Example
          declareGenerators(ZZ[x],{a => 1,b => 2,c => 3})
          cospan(x*a-2*b,x*b-2*c)
    SeeAlso
        cospan
///


doc ///
    Key
        cospan
        (cospan,Sequence)
        (cospan,Thing)
    Headline
        mods out by a collection of module elements
    Usage
        cospan(elements)
    Inputs
        elements:Sequence
            all members of the same ambient module
    Outputs
        :Module
            the quotient of the ambient module by the listed elements
    Description
        Text
            Useful for building modules by generators and relations
        Example
            declareGenerators(ZZ[x],{a => 1, b => 2, c => 3})
            cospan(x*a-2*b,x*b-2*c)
    Caveat
        In M2, multiplication of ring elements by module elements happens on the left, so use
        x*a, not a*x.
    SeeAlso
        declareGenerators
///


doc ///
    Key
        excerptLES
        (excerptLES,ZZ,Module)
        (excerptLES,ZZ,ZZ,Module)
    Headline
        displays a few entries of a long exact sequence
    Usage
        excerptLES(k, l, M)
    Inputs
        M:Module
            over a couple ring, encoding a long exact sequence
        k:ZZ
            the first cohomological degree to be examined
        l:ZZ
            the last (and if omitted, l=k)
    Outputs
        :Net
            depicting rows k..l of the long exact sequence, as well as the entries just before
            and after this row.  There is also some ASCII art showing the snaking paths
            of the connecting maps.
    Description
        Text
            The encoding of a long exact sequence as a module is described in detail in @ TO "longExactSequence" @
        Example
            R = QQ[d]/d^2;
            M = coker map(R^{-1,-1,-1,-1,-1,0,0,0,0},R^{-1,-1,-1,-2,-2,-2,-2,-2},
                {{0,0,0,d,0,0,0,0},{0,0,0,0,d,0,0,0},{0,0,0,0,0,d,0,0},
                 {0,0,0,0,0,0,d,0},{0,0,0,0,0,0,0,d},{-4*d,4*d,4*d,0,0,0,0,0},
                 {5*d,0,0,0,0,0,0,0},{0,5*d,0,0,0,0,0,0},{0,0,5*d,0,0,0,0,0}});
            N = coker map(R^{-1,-1,0,0,0,1,1},R^{0,-1,-1,-2,-2},
                {{0,0,0,d,0},{0,0,0,0,d},{0,0,-4*d,0,0},{0,d,0,0,0},
                 {0,0,3*d,0,0},{0,0,0,0,0},{d,0,0,0,0}});
            f = map(N, M, {{1,1,1,1,1,0,0,0,0},{1,1,1,1,1,0,0,0,0},{d,d,d,d,d,-14,-12,4,4},
                           {0,0,0,0,0,1,1,1,1},{0,0,0,0,0,3,3,3,3},{0,0,0,0,0,d,d,d,0},
                           {0,0,0,0,0,0,0,0,0}});
            LES = longExactSequence f;
            expectExactCouple LES;
            excerptLES(0,2,LES)
    SeeAlso
        longExactSequence
///

doc ///
    Key
        restackRing
        (restackRing,List,Ring)
        restackModule
        (restackModule,List,Module)
    Headline
        changes the order in which variables were adjoined
    Usage
        restackRing(p,R)
    Inputs
        R:Ring
            whose coefficient ring is a quotient ring, and whose coefficient ring's coefficient ring is a quotient ring, etc. for at least n levels
        p:List
            the desired reordering of these levels as a permutation of the list \{1..n\}
    Outputs
        :RingMap
            from R to a new ring where the variables are adjoined in the order determined by p
    Description
        Text
            Here's an example restacking a ring that is four levels deep.
        Example
            A=QQ[x,y, Degrees => {{1,2},{1,2}}]/(x^2+y^2);
            B=A[b];
            C=B[p,q]/(p^3-2*q^3);
            D=C[d];
            restackRing({2,3,4,1}, D)
        Text
            Let $R=QQ[x,y,z]$, and suppose X is a cochain
            complex of R-modules, expressed as a module for the ring R[d]/d^2.  To obtain the module in
            cohomological degree k, evaluateInDegree({k},X) suffices.  But what if we want to evaluate
            in *polynomial* degree k?

            In order to evaluate X in polynomial
            degree g, however, we must restack the ring.
        Example
            R = QQ[x,y,z]; E = R[d]/d^2;
            C = coker map(E^{{0,0},{0,-2},{-1,-3}},E^{{0,-2},{-1,-6}},{{x^2+y^2+z^2,d},{-1,0},{0,-x*y*z}})
            apply(4,k->evaluateInDegree({k},C))
            phi = restackRing({2,1},E);
            C' = phi**C; apply(4,g->prune evaluateInDegree({g},C'))
    Caveat
        Each stage of R may only introduce relations among the most-recent variables.  So, in the last
        example, C=B[p,q]/(p^3-2*q^3) was allowed, but C=B[p,q]/(x*p^3-2*y*q^3) would not be.
    SeeAlso
///

doc ///
    Key
        chainModuleHomology
        (chainModuleHomology, ZZ, Module)
        (chainModuleHomology, Module)
    Headline
        computes the d-cohomology of an R[d]/d^2-module
    Usage
        chainModuleHomology(k,M)
        chainModuleHomology M
    Inputs
        M:Module
            whose ring is of the form R[d]/d^2 for some coefficient ring R
        k:ZZ
            the cohomological degree; if omitted, then the cohomology is
            computed in all degrees and returned as an R[d]/d^2-module
            where d acts by zero.
    Outputs
        :Module
            the d-cohomology of M
    Description
        Text
            We build the cochain complex for the simplicial complex with
            vertices \{a,b,c\}  and facets \{ab,ac,bc\}.  Topologically, this
            is a circle, so the cohomology is QQ^1 in degrees 0 and 1.
        Example
            C = QQ[d]/d^2;
            declareGenerators(C,{a=>0,b=>0,c=>0,ab=>1,ac=>1,bc=>1});
            M = cospan(d*a+ab+ac, d*b-ab+bc, d*c-ac-bc, d*ab, d*ac, d*bc);
            apply(5,i->prune evaluateInDegree({i},M))
            H = chainModuleHomology(M);
            apply(5,i->prune evaluateInDegree({i},H))
            apply(5,i->prune chainModuleHomology(i,M))
    Caveat
    SeeAlso
///

doc ///
    Key
        exactCouple
        (exactCouple,Ring,Module)
        (exactCouple,Module)
    Headline
        builds an exact couple from a R[d,f]/d^2-module
    Usage
        exactCouple(Q, M)
    Inputs
        M:Module
            over a ring of the form R[d,f]/d^2 for
            some coefficient ring R
        Q:Ring
            a couple ring with the same coefficient ring R.  If this argument is omitted,
            a suitable ring is constructed automatically.
    Outputs
        :Module
            over Q, and encoding an exact couple as explained below
    Description
        Text
            Any map of cochain complexes gives a long exact sequence in cohomology.
            Considering M to be a sequence of cochain complexes connected by maps---
            the cochain complex structure comes from the action of d, and the maps
            come from the action of f---we obtain an interlocking sequence of long
            exact sequences.

            The output module encodes the exact couple by placing the page data
            in degrees of the form (2*p, 2*q) and the auxilliary data at the
            midpoints of the differentials.  (The 2s ensure that
            these midpoints are still valid bidegrees.)

            We build the cochain complex for the simplicial complex with vertices
            \{a,b,c\} and facets \{ab,ac,bc\}, placing it in row 0.  In row 1,
            we mod out by (bc); in row 2, by (ac,bc), continuing until every simplex
            is annihilated in row 7.
        Example
            R = QQ[d,t,Degrees=>{{0,1},{1,0}}]/d^2;
            declareGenerators(R,{a=>{0,0},b=>{0,0},c=>{0,0},ab=>{0,1},ac=>{0,1},bc=>{0,1}});
            M = cospan(d*a+ab+ac, d*b-ab+bc, d*c-ac-bc, d*ab, d*ac, d*bc,
                       t*bc, t^2*ac, t^3*ab, t^4*c, t^5*b, t^6*a);

            netList table(7,4,(i,j)->hilbertFunction({6-i,j},M)) -- each row is a cochain complex
            Q = QQ[e_1,f_1,Degrees=>{{-1,1},{2,0}}];
            E1 = exactCouple(Q, M)
            for r from 1 to 7 do (
                print("page " | r |": ");
                print prune pageModule(r,D,E1);
                print " ";
                );
            plotPages((0..7,-2..2,1..7),prune @@ evaluateInDegree,E1)
    Caveat
    SeeAlso
        "Conventions and first examples"
        expectExactCouple
        derivedCouple
///

doc ///
    Key
        derivedCouple
        (derivedCouple,Module)
        (derivedCouple,ZZ,Module)
    Headline
        builds the derived couple of an exact couple
    Usage
        derivedCouple M
    Inputs
        M:Module
            over a couple ring, and encoding an exact couple
    Outputs
        :Module
            over the derived couple ring, and encoding the derived couple
    Description
        Text
            Suppose the ring of M is the couple ring R[e, f].

            Let S be the subring R[e^2, f].  Homogeneous elements of S are restricted to
            an index-two subgroup of the bidegrees of M; as an S-module, M splits as a
            direct sum of its even part and its odd part.  We write A for the odd part
            and E for the even part.  Multiplication by e induces maps from E to A and
            back again.  Since M encodes an exact couple, we have

            image(f : A --> A) = kernel(e : A --> E)

            image(e : A --> E) = kernel(e : E --> A)

            image(e : E --> A) = kernel(f : A --> A).

            The derived couple then replaces A with image(f : A --> A) and E with
            ker(e^2 : E --> E) / im(e^2 : E --> E).  Our grading convention is logical,
            if nonstandard.  Since the differential on E is e^2, which the the composite
            E -e-> A -e-> E, we place A at the midpoint of the differential.  Moreover,
            since the construction of the derived couple makes use of the first isomorphism
            theorem between the image and coimage of f, which therefore have equal claim
            to being A', we place A' at the midpoint of f.  We keep E' in the same
            degrees.

            This all works well, except that the degree of f must be divisible by 2
            so that A' can live at its midpoint.  The current implementation doubles
            all degrees of the input module so that these midpoints exist and are
            unique.  In a future version of M2 that allows grading by a general abelian group,
            the user would be expected to supply a degree that doubles to the degree of f.
        Example
            R = QQ[d,t,Degrees=>{{0,1},{1,0}}]/d^2;
            declareGenerators(R,{a=>{0,0},b=>{0,0},c=>{0,0},ab=>{0,1},ac=>{0,1},bc=>{0,1}});
            M = cospan(d*a+ab+ac, d*b-ab+bc, d*c-ac-bc, d*ab, d*ac, d*bc,
                       t*bc, t^2*ac, t^3*ab, t^4*c, t^5*b, t^6*a);
            Q = QQ[e_1,f_1,Degrees=>{{-1,1},{2,0}}];
            E1 = exactCouple(Q, M)
            expectExactCouple E1;
            E2 = derivedCouple E1
            expectExactCouple E2;
    Caveat
    SeeAlso
        exactCouple
        expectExactCouple
///

doc ///
    Key
        expectExactCouple
        (expectExactCouple,Module)
    Headline
        accepts a module if it encodes an exact couple
    Usage
        expectExactCouple M
    Inputs
        M:Module
            over a couple ring R[e,f]
    Consequences
        Item
            Causes an error if M is not exact.
    Description
        Text
            Let S be the subring R[e^2, f].  Homogeneous elements of S are restricted to
            an index-two subgroup of the bidegrees of M; as an S-module, M splits as a
            direct sum of its even part and its odd part.  We write E for the odd part
            and A for the even part.  Multiplication by e induces maps from E to A and
            back again.  We say that M is exact if

            image(f : A --> A) = kernel(e : A --> E)

            image(e : A --> E) = kernel(e : E --> A)

            image(e : E --> A) = kernel(f : A --> A).

        Example
            R = QQ[d,t,Degrees=>{{0,1},{1,0}}]/d^2;
            declareGenerators(R,{a=>{0,0},b=>{0,0},c=>{0,0},ab=>{0,1},ac=>{0,1},bc=>{0,1}});
            M = cospan(d*a+ab+ac, d*b-ab+bc, d*c-ac-bc, d*ab, d*ac, d*bc,
                       t*bc, t^2*ac, t^3*ab, t^4*c, t^5*b, t^6*a);
            netList table(7,4,(i,j)->hilbertFunction({6-i,j},M))
            Q = QQ[e_1,f_1,Degrees=>{{-1,1},{2,0}}];
            E1 = exactCouple(Q,M);
            -- TODO use CannedExample
            expectExactCouple E1; -- No error
            E1' = E1 / E1_0; -- but expectExactCouple E1' would give the error "failure of exactness at page: ker e != im e."
    Caveat
    SeeAlso
        exactCouple
        derivedCouple
///

doc ///
    Key
        "Conventions and first examples"
    Headline
        specifics on encoding exact couples as modules for a ring
    Description
        Text
            {\bf Encoding an exact couple as a module for a ring}

            An exact couple is a pair of R-modules E and A together with maps A --> E --> A and A --> A
            with the conditions that

            im(A --> E) = ker(E --> A)

            im(E --> A) = ker(A --> A)

            im(A --> A) = ker(A --> E).

            In this package, exact couples are encoded by an action of R[e,f] on the direct sum
            $M = A \oplus E$, where e acts by the maps  A --> E --> A, f acts on A by the map A --> A,
            and on E by 0.  A ring of the form R[e,f]
            form is called a 'couple ring' and the degrees of its variables must satisfy
            the conditions in @ TO expectCoupleRing @.

            The module E is the "page" and the module A is the "auxiliary".  Typically, E and A occupy
            disjoint sets of degrees, termed "even" and "odd".  The terms of the spectral sequence live
            in even degrees, since they come from the page.  The auxiliary data lives in odd degree.

            {\bf Constructing exact couples}

            Exact couples arise in algebra any time a chain complex $C=(C,d)$ carries a self map
            $f : C \to C$.  Setting $A = H_* C$ to be the homology of $C$, and $E=H_* C(f)$ to be
            the homology of the mapping cone of $f$, we have a long exact sequence

            ... ----> A ----> E ----> A --f--> A ----> ...

            defining the maps in the exact couple.

            The all-purpose method
            @ TO exactCouple @ accepts as input C considered as an $R[d,f]/d^2$ module.
            Note that the action of $d$ and $f$ is encoded by the ring action.


            To build exact couples by hand, use @ TO declareCouple @, and check your work with
            @ TO isHomogeneous @ and @ TO expectExactCouple @.  For the usual long
            exact sequence induced by a map of complexes, use @ TO longExactSequence @.
            Once you have an exact couple, @ TO derivedCouple @ produces another, as we now explain.

            {\bf Derived couples}

            Every exact couple determines a new exact couple by the formulas $A' = im(f)$ and
            $E' = ker(e^2)/im(e^2)$.  The maps $A' \to E' \to A'$ are induced by $e$, and the
            map $A' \to A'$ is induced by $f$.  There is a wrinkle, however: the map
            $A' \to E'$ factors $A' \to A / ker(f) \to E'$, making use of the first
            isomorphism theorem.  Some degree confusion can result.  Indeed, the inclusion
            $A' \subseteq A$, which is induced by the identity, is a degree 0 map; in contrast,
            the projection $A / (ker f) \to A'$, which is induced by $f$, has the same degree
            as $f$ itself.  This forces the inverse map $A' \to A / ker(f)$ to have degree $-deg(f)$.

            If we had taken A' to be the coimage of f instead of the image, 
            this would have resulted in different degrees.  So either of these conventions
            would have an assymetry: an unjustified preference for image or coimage.

            Our convention avoids this assymetry by averaging: we place A' exactly halfway between
            (image f) and (coimage f).  This has the effect of giving the comparison
            maps $A \to A' \to A$ the same degree, namely, $(degree f)/2$.  (This division-by-two
            accounts for the requirement in @ TO expectCoupleRing @ that the degree of $f$ be even.)

            The new degrees can be expressed in terms of the old:


            $ degree f' = degree f$


            $ degree e' = (degree e) - (degree f) / 2$.


            The resulting ring $R[e',f']$ is the derived couple ring of $R[e,f]$; it acts on the
            derived couple, and can be obtained using @ TO derivedCoupleRing @.
        Text
            {\bf How to determine $deg(e)$ and $deg(f)$ in your example}
            
            Suppose you have in mind a particular spectral sequence of $R$-modules,
            starting on page k,
            and you want to program
            its exact couple.  What degrees should you use for the couple ring $R[e_k,f_k]$?
            
            Think about the degree of $D_k$, the differential on the starting page.  This tells
            you the degree of $e_k$:
            
            $deg(e_k) = deg(D_k)$.
            
            Now think about the degree of $D_{k+1}$, the differential on the next page.  This
            lets you compute the degree of $f_k$:
            
            $deg(f_k) = 2 * deg(D_k) - 2 * deg(D_{k+1})$.
    SeeAlso
        exactCouple
        expectExactCouple
        derivedCouple
        contravariantExtCouple
        covariantExtCouple
        TorCouple
        "Bockstein spectral sequence"
        "Serre spectral sequence in homology"
///

doc ///
    Key
        "Bockstein spectral sequence"
    Headline
        a singly-graded spectral sequence built from the chain self-map "multiplication by p"
    Description
        Text
            {\bf Bockstein Spectral Sequence}

            Let p be a prime number, and suppose C is a chain complex over the integers; then
            multiplication by p induces a chain map $C --> C$, and so we have a chain complex
            with a self-map, placing us in the algebraic context to obtain an exact couple.

            For example, let C be the cellular cochain complex for the real projective space
            $\mathbb{R}P^3$ and its usual cell structure with a single cell in each degree:

            Z --0--> Z --2--> Z --0--> Z ----> 0

            Name the classes p0, p1, p2, and p3, specify the differential by imposing relations of
            the form d*pk = d(pk), and set t to act by 2
            by tensoring with R^1/(t-2) (this is a convenient way to impose the relation that
            every generator g has t*g = 2*g):
        Example
            Q = ZZ[d, f, Degrees => {1,0}]/d^2;
            declareGenerators(Q, {p0 => 0, p1 => 1, p2 => 2, p3 => 3});
            C = cospan(d*p0, d*p1-2*p2, d*p2, d*p3) ** Q^1/(f-2); C
            isHomogeneous C
        Text
            This C is the right sort of module to give to exactCouple since it carries an action
            of a ring of the form R[d,f]/d^2
        Example
            bock = exactCouple C
            expectExactCouple bock
            P1 = prune pageModule(1,D,bock)
        Text
            Since the generators of the E_1-page are annihilated by 2, the same will be true on subsequent pages.
        Example
            P2 = prune pageModule(2,D,bock)
            P3 = prune pageModule(3,D,bock)
        Text
            It is always the case that the the pages
            of the Bockstein spectral sequence are defined over the field ZZ/p; indeed this is its
            main useful property.
        Example
            P1' = prune(map((ZZ/2)[D_1],ring P1) ** P1)
            P2' = prune(map((ZZ/2)[D_2],ring P1) ** P1)
            P3' = prune(map((ZZ/2)[D_3],ring P1) ** P1)
    SeeAlso
        pageModule
        derivedCouple
        expectExactCouple
        "Conventions and first examples"
        "Serre spectral sequence in homology"
///

doc ///
    Key
        "Serre spectral sequence in homology"
    Headline
        exact couple associated to a fibration
    Description
        Text
            {\bf Homology Serre spectral sequence for the Hopf fibration}

            The Serre spectral sequence can be constructed by choosing a cell structure
            on the base space, and looking at the induced filtration on the total
            space.  In the case of the Hopf fibration $p: S^3 \to S^2$, and making use
            of the easy cell structure on S^2 with only two cells, we obtain the
            following filtration:

            $X_k = \emptyset$ for $k \leq 0$

            $X_0 = X_1 = S^1$

            $X_k = S^3$ for $k \geq 2$

            --The homology exact couple of a filtered space has $A = \oplus_{p,q} A_{p,q}$ with
            --$A_{p,q} = H_p X_q$ and $E_{p,q} = H_p(X_q , X_{q-1})$.  The maps are then
            The homology exact couple of a filtered space is built from the various long exact
            sequences associated to the pairs $(X_p , X_{p-1})$.

            Its first page is the relative homology, and its first auxiliary is the
            (absolute) homology.  The variables act as follows:

            $f : H_n X_p \to H_n X_{p+1}$

            $e : H_n X_p \to H_n(X_p , X_{p-1})$

            $e : H_n(X_p , X_{p-1}) \to  H_{n-1} X_{p-1}$,

            where this last map is the usual connecting homomorphism.  In our example, the most
            interesting long exact sequence concerns the pair $(S^3,S^1)$.  Since the homology of
            $S^3$ and $S^1$ are pretty sparse, the homology of the pair is determined by isomorphisms

            $e : H_3 S^3 \to H_3(S^3,S^1)$ since $H_3$ and $H_2$ vanish on $S^1$

            $e : H_2(S^3,S^1) \to H_1 S^1$ since $H_2$ and $H_1$ vanish on $S^3$

            $H_1(S^3,S^1) = 0$ since $H_1 S^3 = 0$ and $H_0 S^1 \to H_0 S^3$ is an injection

            $H_0(S^3,S^1) = 0$ since $H_0 S^1 \to H_0 S^3$ is a surjection.

            This analysis shows that the four groups $H_0 X_0$, $H_1 X_0$, $H_2(X_2, X_1)$, and
            $H_3 X_2$ are free abelian of rank one.
            Let $x \in H_0 X_0$, $y \in H_1 X_0$, $z \in H_2(X_2,X_1)$, and $w \in H_3 X_2$
            be generating classes for these groups.

            By the general properties of exact couples, we have that $e^2$ annihilates auxiliary
            generators and $f$ annihilates page generators.  It turns out that
            the only remaining relation in the exact couple is $e*z - f*y$, which says that
            the connecting map $H_2(X_2,X_1) \to H_1(X_1)$ sends $z$ to the same place as the
            filtation map $H_1(X_0) \to H_1(X_1)$ sends $y$.  (Depending on your choices for
                generators, this relation may read $e*z + f*y$.)

            It is straightforward to give all this information to M2, but determining
            degrees does take a bit of thought the first time you do it.  We continue the example
            by {\bf analyzing degrees ...}

            Since we have specific degree preferences (exactly
            double the usual ones for a Serre spectral sequence) we will
            do this by hand rather than relying on the default degrees provided by @ TO coupleRing @.

            The usual Serre spectral sequence has $E^1_{pq} = H_{p+q}(X_p , X_{p-1})$, so we place
            this module in degree \{2p,2q}.  The module $A^1_{pq} = H_{p+q} X_p$ sits
            halfway along the map $E^1_{pq} \to E^1_{(p-1)q}$, so it has degree \{2p-1,2q}.

            In usual indexing, the differential has degree \{-1,0}, so our differential has degree
            \{-2,0}.  On the other hand, since the differential is given by $e^2$, this means that
            the degree of $e$ itself is \{-1,0}.  (And in general, the degree of e in our conventions
            should be the degree of the differential in the usual conventions.)

            Everything so far has been about the first page, $E^1$.  However, the easiest way to
            determine the degree of $f$ is to consider the second page.  The usual conventions give
            that the differential on $E^2$ has degree \{-2,1}, and so (by the same argument as above)
            this is the degree of $e_2$.  We have generally $deg e_{k+1} = deg e_k - (deg f_k)/2$.
            Solving, we see that the degree of $f = f_1$ is 2 * (\{-1,0\} - \{-2,1\}) = \{2,-2}.

            We conclude the general rule that{\bf ... the degree of e_k in our convention equals the
                degree of the page-k differential in the standard convention; the degree of f_k is then
                given by $2 * (deg e_k - deg e_{k+1})$.}

            We are now able to set up the couple ring.
        Example
            R = QQ;
            Q = R[e_1,f_1,Degrees=>{{-1,0},{2,-2}}];
        Text
            Now it is time to build the couple.  We give our couple ring to the function
             @ TO declareCouple @, together with names and degrees for our generating classes.  The
            page generators are always given first, and the auxiliary generators second.
        Example
            declareCouple(Q, {z => {4,0}}, {x => {1,0}, y => {1,2}, w => {5,2}})
        Text
            The couple we have built is "free" in the sense that the only relations imposed
            are the tautologous ones that hold in every exact couple (f acts by zero on the page,
                and e^2 acts by zero on the auxiliary).  To obtain our couple, we must impose the
            relation $e*z-f*y$.
        Example
            C = cospan(e_1*z-f_1*y)
            isHomogeneous C
            expectExactCouple C
        Text
            Since @ TO expectExactCouple @ accepts C without an error, we truly have an exact couple,
            and are ready to compute its spectral sequence.

            The next line displays the some entries in the first four pages of the
            spectral sequence determined by C.
        Example
            plotPages((-2..4,-2..3,1..4), prune @@ evaluateInDegree,C)
    SeeAlso
        exactCouple
        expectExactCouple
        derivedCouple
        pageModule
        derivedCouple
        expectExactCouple
        "Conventions and first examples"
        "Bockstein spectral sequence"
///

doc ///
    Key
        contravariantExtLES
        (contravariantExtLES, ZZ, Module, Module, Module)
    Headline
        the long exact sequence in Ext induced by an inclusion in the first coordinate of Hom
    Usage
        contravariantExtLES(k,X,A,Y)
    Inputs
        k:ZZ
            number of rotations of the long exact sequence to display
        X:Module
        A:Module
            a submodule of X
        Y:Module
            giving a hom-functor Hom(-,Y)
    Outputs
        :Net
            showing the long exact sequence in Ext
    Description
        Text
            The long exact sequence is returned as a Net with the following general format:
        CannedExample
            |   .- ->        ...    (k-3) more rows appearing
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->  Ext^1(X,Y)  - - ->  Ext^1(A,Y)  - - ->  Ext^2(X/A,Y)  - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->   Hom(X,Y)   - - ->   Hom(A,Y)   - - ->  Ext^1(X/A,Y)  - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->      0       - - ->       0      - - ->   Hom(X/A,Y)   - -'
            |  (
            |   \
            |
            |                                                                   \
            |                                                        0       - -'
        Text
            The next example gives a typical use.
        Example
            R = QQ[x]
            X = R^1 / x^9
            A = image map(X,,{{x^7}})
            Y = coker map(R^1,,{{x^3}})
            contravariantExtLES(3,X,A,Y)
            apply(2, p -> prune Ext^p(X,Y))
            apply(2, p -> prune Ext^p(A,Y))
            apply(2, p -> prune Ext^p(X/A,Y))
    Caveat
        For computational access to the maps in the sequence, use contravariantExtCouple instead.
    SeeAlso
        excerptLES
        contravariantExtCouple
///

doc ///
    Key
        covariantExtLES
        (covariantExtLES, ZZ, Module, Module, Module)
    Headline
        the long exact sequence in Ext induced by an inclusion in the last coordinate of Hom
    Usage
        covariantExtLES(k,W,X,A)
    Inputs
        k:ZZ
            number of rotations of the long exact sequence to display
        W:Module
            giving a hom-functor Hom(W,-)
        X:Module
        A:Module
            a submodule of X
    Outputs
        :Net
            showing the long exact sequence in Ext
    Description
        Text
            The long exact sequence is returned as a Net with the following general format:
        CannedExample
            |   .- ->        ...    (k-3) more rows appearing
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->  Ext^1(W,X)  - - ->  Ext^1(W,X/A)  - - ->  Ext^2(W,A)  - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->   Hom(W,X)   - - ->   Hom(W,X/A)   - - ->  Ext^1(W,A)  - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->      0       - - ->        0       - - ->   Hom(W,A)   - -'
            |  (
            |   \
            |
            |                                                                   \
            |                                                        0       - -'
        Text
            The next example gives a typical use.
        Example
            R = QQ[x]
            X = R^1 / x^9
            A = image map(X,,{{x^7}})
            W = coker map(R^1,,{{x^3}})
            covariantExtLES(3,W,X,A)
            apply(2, p -> prune Ext^p(W,X))
            apply(2, p -> prune Ext^p(W,X/A))
            apply(2, p -> prune Ext^p(W,A))
    Caveat
        For computational access to the maps in the sequence, use covariantExtCouple instead.
    SeeAlso
        excerptLES
        covariantExtCouple
///

doc ///
    Key
        covariantExtCouple
        (covariantExtCouple,Module,List)
    Headline
        the exact couple obtained by applying Ext(W,-) to a filtered module
    Usage
        covariantExtCouple(W, submods)
    Inputs
        W:Module
            for some ring R, giving a functor Hom(W,-)
        submods:List
            of R-modules \{A_0, A_1, ..., A_m} with each A_i inside A_{i+1}
    Outputs
        M:Module
            an exact couple
    Description
        Text
            For notational convenience, set $X = A_m$, and
            extend the sequence $A_i$ to all $i \in \ZZ$
            by setting $A_i = 0$ for $i < 0$, and $A_i = X$ for $i > m$.  
            
            The returned couple $M$ is a module for  
            the ring R[e_1,f_1,Degrees=>\{\{1,-1},\{0,2}}].
            We describe the module $M$ in every bidegree $\{s,t}$.  The description
            depends on the parity of $s$ and $t$.
            
            
            If $s$ and $t$ are both even, say $\{s,t} = \{2p 2q}$, then
    
            $M_{s,t} = Ext^p(W, A_q/A_{q-1})$;
            
            if $s$ and $t$ are both odd, say $\{s,t} = \{2p-1,2q+1}$, then
    
            $M_{s,t} = Ext^p(W, A_q)$;
            
            and otherwise, if $s$ and $t$ sum to an odd number, then $M_{s,t} = 0$.  
            
            The variables $e_1$ and $f_1$ act by the maps in the 
            various long exact sequences
            
            $Ext^p(W, A_{q-1}) \to Ext^p(W, A_q) \to 
             Ext^p(W, A_q / A_{q-1}) \to Ext^{p+1}(W, A_{q-1})$.
            
            {\bf Associated spectral sequence}
            
            The spectral sequence associated to this couple converges to $Ext^p(W,X)$.
            The differential on page $r$ has bidegree \{1,-r}.  The first page has
            
            $E^{pq}_1 = Ext^p(W,A_q/A_{q-1})$.
            
            Setting $F^p_q = image(Ext^p(W,A_q) \to Ext^p(W,X))$, the infinity page has
            
            $E^{p,q}_{\infty} = F^p_{q} / F^p_{q-1}$.
        Example
            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}));
            for m in submods do print m;
            W = coker map(R^1,,{{x^3}})
            couple = prune covariantExtCouple(W,submods)
            expectExactCouple couple
            plotPages((-1..2,-1..5,1..3), prune @@ evaluateInDegree, couple)
            A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
            E1 = (q,p) -> prune Ext^p(W,A(q)/A(q-1));
            netList reverse table(5,2,E1)
            inc = q -> inducedMap(X,A(q));
            filt = (p,q) -> image Ext^p(W,inc q);
            Einfty = (q,p) -> prune(filt(p,q)/filt(p,q-1));
            netList reverse table(5,2,Einfty)
    SeeAlso
        contravariantExtCouple
///

doc ///
    Key
        contravariantExtCouple
        (contravariantExtCouple,List,Module)
    Headline
        the exact couple obtained by applying Ext(-,Y) to a filtered module
    Usage
        contravariantExtCouple(submods, Y)
    Inputs
        Y:Module
            for some ring R, giving a functor Hom(-,Y)
        submods:List
            of R-modules \{A_0, A_1, ..., A_m} with each A_i inside A_{i+1}
    Outputs
        M:Module
            an exact couple
    Description
        Text
            For notational convenience, set $X = A_m$, and
            extend the sequence $A_i$ to all $i \in \ZZ$
            by setting $A_i = 0$ for $i < 0$, and $A_i = X$ for $i > m$.  
            
            The returned couple $M$ is a module for  
            the ring R[e_1,f_1,Degrees=>\{\{1,1},\{0,-2}}].  
            We describe the module $M$ in every bidegree $\{s,t}$.  The description
            depends on the parity of $s$ and $t$.
            
            
            If $s$ and $t$ are both even, say $\{s,t} = \{2p 2q}$, then
    
            $M_{s,t} = Ext^p(A_q / A_{q-1}, Y)$;
            
            if $s$ and $t$ are both odd, say $\{s,t} = \{2p-1,2q-1}$, then
    
            $M_{s,t} = Ext^p(X / A_{q-1}, Y)$;
            
            and otherwise, if $s$ and $t$ sum to an odd number, then $M_{s,t} = 0$.  
            
            The variables $e_1$ and $f_1$ act by the maps in the 
            various long exact sequences
            
            $Ext^p((X / A_q), Y) \to Ext^p((X / A_{q-1}), Y) \to 
             Ext^p((A_q / A_{q-1}), Y) \to Ext^{p+1}((X / A_q), Y)$.
            
            {\bf Associated spectral sequence}
            
            The spectral sequence associated to this couple converges to $Ext^p(X,Y)$.
            The differential on page $r$ has bidegree \{1,r}.  The first page has
            
            $E^{p,q}_1 = Ext^p(A_q/A_{q-1},Y)$.
            
            Setting $F^p_q = image(Ext^p((X/A_q),Y) \to Ext^p(X,Y))$, the infinity page has
            
            $E^{p,q}_{\infty} = F^p_{q-1} / F^p_q$.
        Example
            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}));
            for m in submods do print m;
            Y = coker map(R^1,,{{x^3}})
            couple = prune contravariantExtCouple(submods,Y)
            expectExactCouple couple
            plotPages((-1..2,-1..5,1..3), prune @@ evaluateInDegree, couple)
            A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
            E1 = (q,p) -> prune Ext^p(A(q)/A(q-1),Y);
            netList reverse table(5,2,E1)
            proj = q -> inducedMap(X/A(q),X);
            filt = (p,q) -> image Ext^p(proj q,Y);
            Einfty = (q,p) -> prune(filt(p,q-1)/filt(p,q));
            netList reverse table(5,2,Einfty)
        Text
            It seems to me that this is the same spectral sequence as the one you would get
            from the couple
            
            $Ext^p((A_q), Y) \to Ext^p((A_{q-1}), Y) \to 
             Ext^p((A_q / A_{q-1}), Y) \to Ext^{p+1}((A_q), Y)$;
            
            If I learn of a proof of this fact, then I will put the reference here.
    SeeAlso
        covariantExtCouple
///

doc ///
    Key
        filtrationModule
        (filtrationModule,Ring,List)
    Headline
        converts a filtered module to an R[t]-module
    Usage
        filtrationModule(Q, L)
    Inputs
        Q:Ring
            of the form R[t] for some coefficient ring R and variable t
        L:List
            of the form \{ M_0, M_1, ..., M_k \}  for some k > 0
            with each M_i $\subseteq$ M_{i+1} an inclusion of R-modules
    Outputs
        :Module
            graded, with M_0 sitting in degree zero, and where t acts by inclusions
    Description
        Example
            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}))
            Q = R[t]
            filtrationModule(Q, submods)
    Caveat
        The ring Q should be valid for expectSequenceRing.
        The list L should be valid for expectFiltrationList.
    SeeAlso
        sequenceModule
///

doc ///
    Key
        expectChainRing
        (expectChainRing,Ring)
    Headline
        accepts rings of the form R[d]/d^2
    Usage
        expectChainRing Q
    Inputs
        Q:Ring
    Consequences
        Item
            causes an error if Q is not of the form R[d]/d^2
    Description
        Text
            Specifically, Q should have exactly one generator over
            its coefficient ring, and that generator should square to zero.
            There is no restriction on the degree of the generator.
        Example
            Q = QQ[d]/d^2;
            expectChainRing Q
///


doc ///
    Key
        expectTriangleRing
        (expectTriangleRing,Ring)
    Headline
        accepts certain rings of the form R[d,e,f]/(d^2, e^3)
    Usage
        expectTriangleRing Q
    Inputs
        Q:Ring
    Consequences
        Item
            causes an error if Q is not of the form R[d,e,f]/(d^2, e^3)
            with $deg(d) = deg(e^2f)$
    Description
        Text
            Specifically, Q should have exactly three generators over
            its coefficient ring, with $d^2 = e^3 = 0$.
            The degree of $d$ must equal the degree of $e^2f$.
        Example
            Q = triangleRing(QQ, d, e, f)
            expectTriangleRing Q
    SeeAlso
        triangleRing
///

doc ///
    Key
        expectCoupleRing
        (expectCoupleRing,Ring)
    Headline
        accepts certain rings of the form R[e_r,f_r], and installs Page, isEvenDegree, and isOddDegree
    Usage
        expectCoupleRing Q
    Inputs
        Q:Ring
    Consequences
        Item
            causes an error if Q is not of the form R[e_r,f_r]
            with $deg(f)$ divisible by two
        Item
            Sets Q.Page, Q.isEvenDegree, and Q.isOddDegree
    Description
        Text
            Specifically, Q should have exactly two generators over
            its coefficient ring, and these must be indexed variables
            with the same subscript.  Moreover, the degree of $f$ must be
            divisible by two.

            Q.Page = r

            Q.isEvenDegree is set to be a boolean-valued function so that


            Q.isEvenDegree(deg) = not Q.isOddDegree(deg)


            Q.isEvenDegree(deg + degree e_r) = Q.isOddDegree(deg)


            Q.isEvenDegree(deg + degree f_r) = Q.isEvenDegree(deg)

        Example
            Q = coupleRing(QQ, 7, e, f)
            describe Q
            expectCoupleRing Q
            Q.Page
            netList table(5,10,(i,j)->Q.isEvenDegree({i,j}))
            netList table(5,10,(i,j)->Q.isOddDegree({i,j}))
    SeeAlso
        coupleRing
///


doc ///
    Key
        triangleRing
        (triangleRing,Ring,Symbol,Symbol,Symbol)
    Headline
        builds a triangle ring
    Usage
        triangleRing(R,d,e,f)
    Inputs
        R:Ring
            the coefficients
        d:Symbol
        e:Symbol
        f:Symbol
    Outputs
        :Ring
            R[d,e,f,Degrees=>\{\{0,2},\{1,0},\{-2,2}}]/(d^2,e^3)
    Description
        Text
           These default degrees are chosen to meet the conditions of expectTriangleRing.
           Other degrees are possible using the optional argument Degrees or by building
           the triangle ring by hand directly.
        Example
           Q = triangleRing(QQ,d,e,f)
           describe Q
           expectTriangleRing Q
    SeeAlso
        expectTriangleRing
///

doc ///
    Key
        coupleRing
        (coupleRing,Ring,ZZ,Symbol,Symbol)
    Headline
        builds a couple ring
    Usage
        coupleRing(R,r,e,f)
    Inputs
        R:Ring
            the coefficients
        r:ZZ
            the page number
        e:Symbol
        f:Symbol
    Outputs
        :Ring
            R[e_r,f_r,Degrees=>\{\{1,0},\{-2,2}}]
    Description
        Text
           These default degrees are chosen to meet the conditions of expectCoupleRing.
           Other degrees are possible using the optional argument Degrees or by building
           the couple ring by hand directly.
        Example
           Q = coupleRing(QQ,7,e,f)
           describe Q
           expectCoupleRing Q
    SeeAlso
        expectCoupleRing
///

doc ///
    Key
        excerptCouple
        (excerptCouple,List,ZZ,Module)
    Headline
        displays one of the long exact sequences in an exact couple
    Usage
        excerptCouple(deg, k, C)
    Inputs
        C:Module
            an exact couple given as a module for a couple ring Q
        deg:List
            an exterior multidegree of C, the starting degree; should be an even degree
            in the sense of Q.isEvenDegree(deg).  This degree appears at the bottom of
            the middle column in the displayed long exact sequence
        k:ZZ
            the number of winds of the long exact sequence to display
    Outputs
        :Net
            the long exact sequence starting in degree deg and continuing for
            3*k terms
    Description
        Example
            R = QQ[x]
            X = R^1 / x^9
            A = apply(5,k->image map(X,,{{x^(8-2*k)}}))
            W = coker map(R^1,,{{x^3}})
        Text
            We build the exact couple coming from applying Hom(W,-) to a
            filtered module
        Example
            couple = prune covariantExtCouple(W,A)
            expectExactCouple couple
        Text
            We check that \{0,4} is an even degree, and then use excerptCouple
        Example
            Q = ring couple
            Q.isEvenDegree({0,4})
            excerptCouple({0,4},2,couple)
        Text
            The middle column of the displayed long exact sequence is in even degrees; in other
            words, its entries appear on the E_1 page of the couple.  Specifically, 
            $E_1^{pq} = Ext^p(W,A_q/A_{q-1})$
            can be found in degree \{2p,2q}.  The bottom element of the middle
            column is degree \{0,4}, which is then $Ext^0(W,A_2/A_1)$.  
            The top element of the middle column is $Ext^1(W,A_2/A_1)$.
        Example
            prune Ext^0(W,A_2/A_1)
            prune Ext^1(W,A_2/A_1)
    SeeAlso
        covariantExtCouple
///

doc ///
    Key
        pageModule
        (pageModule,ZZ,Symbol,Module)
        (pageModule,ZZ,IndexedVariableTable,Module)
    Headline
        gives a page of a spectral sequence as a module for R[d]/d^2 where d is the differential
    Usage
        pageModule(r, D, C)
    Inputs
        C:Module
            an exact couple, acted on by a couple ring Q = R[e,f]
        D:Symbol
            representing the differential
        r:ZZ
            the page number; must have r \ge Q.Page
    Outputs
        :Module
            The page E_r of the exact couple C considered as a module for R[D_r]/D_r^2 where D_r is the
            differential on page r.
    Description
        Text
            We show how to use pageModule to study the differentials in a spectral sequence.  The following
            lines construct the homological Serre spectral sequence for the Hopf fibration S^3 \to S^2.
        Example
            Q = coupleRing(ZZ,1,e,f,Degrees=>{{-1,0},{2,-2}})
            declareCouple(Q, {z => {4,0}}, {x => {1,0}, y => {1,2}, w => {5,2}})
            C = cospan(e_1*z-f_1*y)
            isHomogeneous C
            expectExactCouple C
        Text
            We extract the E^1 page as a module for D_1
        Example
            prune pageModule(1,D,C)
        Text
            Note that the differential annihilates all four generators.  We now extract the
            E^1 page with its differential, D_2:
        Example
            prune pageModule(2,D,C)
            degree D_2
        Text
            This time, the generator in degree \{2,0} maps via D_2 to a nontrivial element in degree \{0,1}.
            Since the module has no additional generators in that degree, the differential is an isomorphism
            between these two degrees.  Computing the next page shows the cancellation:
        Example
            prune pageModule(3,D,C)
    SeeAlso
        plotPages
///

doc ///
    Key
        plotPages
        (plotPages,Sequence,Function,Module)
    Headline
        displays a few pages of a spectral sequence
    Usage
        plotPages(ranges,f,C)
    Inputs
        C:Module
            an exact couple, acted on by a couple ring Q = R[e,f]
        f:Function
            accepting arguments f(deg, M) with deg an exterior degree of Q and M an iterated
            derived couple of C
        ranges:Sequence
            with three terms (ps,qs,rs), each itself a sequence; rs indicates the page numbers to be displayed;
            ps the entries on the p axis, and qs the entries on the q axis.  We allow rs to be a single integer
            if only one page is to be displayed.
    Outputs
        :Net
    Description
        Text
            The following
            lines construct the homological Serre spectral sequence for the Hopf fibration S^3 \to S^2.
        Example
            Q = coupleRing(ZZ,1,e,f,Degrees=>{{-1,0},{2,-2}})
            declareCouple(Q, {z => {4,0}}, {x => {1,0}, y => {1,2}, w => {5,2}})
            C = cospan(e_1*z-f_1*y)
            isHomogeneous C
            expectExactCouple C
        Text
            Use plotPages to show the first three pages
        Example
            plotPages((0..3,0..2,1..3), prune @@ evaluateInDegree, C)
        Text
            and the tenth page
        Example
            plotPages((0..3,0..2,10), prune @@ evaluateInDegree, C)
        Text
            The usual choices for f are evaluateInDegree or hilbertFunction:
        Example
            plotPages((0..3,0..2,1..3), hilbertFunction, C)
    Caveat
        This function assumes that the couple C is bi-graded so that $deg(e)$ and $deg(f)$ are lists
        of length two.  If this is not the case, then you can still form derived couples and probe them
        using evaluateInDegree, but plotPages will not work.
    SeeAlso
        derivedCouple
        evaluateInDegree

///

doc ///
    Key
        declareCouple
        (declareCouple,Ring,List,List)
    Headline
        initializes generating classes for an exact couple
    Usage
        declareCouple(Q, pageClasses, auxClasses)
    Inputs
        Q:Ring
            a couple ring
        pageClasses:List
            of the form \{... , generatorName => degree, ...} with generatorName a Symbol and degree a List.
            The degrees should be external degrees for the ring Q.  Since these classes appear on the page
            of the spectral sequence, their degrees should be even, i.e., Q.isEvenDegree(degree) should return true.
        auxClasses:List
            also of the form \{... , generatorName => degree, ...}; the degrees should be odd external degrees
            for Q.
    Outputs
        :Module
            a free exact couple with the specified generators
    Consequences
        Item
            assigns appropriate values to the symbols supplied as generator names
    Description
        Text
            We build the homological Serre spectral sequence for the Hopf fibration S^3 \to S^2
        Example
            Q = coupleRing(ZZ,1,e,f,Degrees=>{{-1,0},{2,-2}})
            declareCouple(Q, {z => {4,0}}, {x => {1,0}, y => {1,2}, w => {5,2}})
            C = cospan(e_1*z-f_1*y)
            isHomogeneous C
            expectExactCouple C
    Caveat
        The output of declareCouple is a free exact couple, but it is not a free module for Q; the
        tautological couple relations are enforced.  If Q = R[e,f], then we have that f annihilates
        every page generator and e^2 annihilates every aux generator.  We also have that e*f and e^3
        act by zero on all generators.
    SeeAlso
        enforceCoupleRelations
///

doc ///
    Key
        derivedCoupleRing
        (derivedCoupleRing,Ring)
    Headline
        forms the ring that acts on a derived couple
    Usage
        derivedCoupleRing Q
    Inputs
        Q:Ring
            a couple ring.  (see @ TO expectCoupleRing @)
    Outputs
        :Ring
            a couple ring with incremented subscripts, and modified degrees
    Description
        Text
            Suppose Q = R[e_r,f_r].  The derived couple ring of Q will be R[e_{r+1},f_{r+1}].
            The degree of f will not change, but the degree of e does transvect against the
            direction of f.  Specicially, in our convention,
            the degree of f is assumed to be even, and the new degree of e is given by the
            formula


            $deg e_{r+1} = deg(e_r) - deg(f)/2$.
        Example
            Q = coupleRing(ZZ,1,e,f,Degrees=>{{-1,0},{2,-2}})
            Q' = derivedCoupleRing Q
            degree Q'_0
            degree Q'_1
    SeeAlso
        derivedCouple
///

doc ///
    Key
        enforceCoupleRelations
        (enforceCoupleRelations,Module)
    Headline
        mods out by tautological relations satisfied by every exact couple
    Usage
        enforceCoupleRelations M
    Inputs
        M:Module
            for a couple ring Q (see @ TO expectCoupleRing @)
    Outputs
        :Module
            M / (tautological couple relations)
    Description
        Text
            If Q = R[e,f] and M is an exact couple, then f annihilates the even degrees of M
            and e^2 annihilates the odd degrees of M.  (In this context, even and odd are determined
            by the functions Q.isEvenDegree and Q.isOddDegree).
        Example
            Q = coupleRing(ZZ,1,e,f,Degrees=>{{-1,0},{2,-2}})
        Text
            A free Q module is not a free couple because the tautological couple relations do not hold
        CannedExample
            i2 : expectExactCouple Q^{{0,0},{-1,0},{-2,0}}
                 error: e^2 fails to annihilate aux
        Text
            To obtain the free couple as a quotient of the free module, use enforceCoupleRelations
        Example
            C = enforceCoupleRelations Q^{{0,0},{-1,0},{-2,0}}
            expectExactCouple C
    SeeAlso
        expectExactCouple
///


doc ///
    Key
        expectFiltrationList
        (expectFiltrationList,List)
    Headline
        accepts a list of modules if each includes in the next
    Usage
        expectFiltrationList L
    Inputs
        L:List
            of the form \{..., M_i, M_{i+1}, ...}; M_i modules for the same ring, and M_i \subseteq M_{i+1}.
    Consequences
        Item
            causes an error if L is not of the indicated form, or if L is empty
    Description
        Text
            The following code is run by expectFiltrationList, and must return true
            to avoid an error.
        CannedExample
            all(#L - 1, q -> isSubset(L#q, L#(q+1)))
    SeeAlso
        filtrationModule
///


doc ///
    Key
        externalDegreeIndices
        (externalDegreeIndices,Ring)
    Headline
        for a ring Q, returns the degree-coordinates present in Q but not in its coefficient ring
    Usage
        externalDegreeIndices Q
    Inputs
        Q:Ring
            with some coefficient ring R
    Outputs
        :List
            given by \{0, ..., (degreeLength Q) - (degreeLength R) - 1}
    Description
        Text
            We build a polynomial ring with coefficients in a polynomial ring
        Example
            R = QQ[x,y]
            Q = R[s,t]
            degree s
            degree t
        Text
            Note that the generators of Q have two degrees: an "internal" degree coming from R and
            an "external" degree coming from Q proper.
        Example
            external = externalDegreeIndices Q
            (degree s)_external
    SeeAlso
        internalDegreeIndices
///


doc ///
    Key
        internalDegreeIndices
        (internalDegreeIndices,Ring)
    Headline
        for a ring, returns the degree-coordinates of its coefficient ring
    Usage
        internalDegreeIndices Q
    Inputs
        Q:Ring
            with some coefficient ring R
    Outputs
        :List
            given by \{(degreeLength Q) - (degreeLength R), ..., (degreeLength Q) - 1}
    Description
        Text
            We build a polynomial ring with coefficients in a polynomial ring
        Example
            R = QQ[x,y]
            Q = R[s,t]
            degree s
            degree x
        Text
            Note that the generators of Q have two degrees: an "internal" degree coming from R and
            an "external" degree coming from Q proper.
        Example
            internal = internalDegreeIndices Q
            (degree s)_internal
            (degree x_Q)_internal
    SeeAlso
        externalDegreeIndices
///


doc ///
    Key
        oneEntry
        (oneEntry,List,List,RingElement)
        (oneEntry,List,Nothing,RingElement)
        (oneEntry,Nothing,List,RingElement)
        (oneEntry,Nothing,ZZ,RingElement)
        (oneEntry,ZZ,Nothing,RingElement)
        (oneEntry,ZZ,ZZ,RingElement)
    Headline
        builds a one-by-one matrix
    Usage
        oneEntry(row,col,ent)
    Inputs
        ent:RingElement
            drawn from some ring R
        row:List
            a degree for the ring R
        col:List
            a degree for the ring R
    Outputs
        :Matrix
            a single-entry matrix with row degree row, column degree col, and entry ent
    Description
        Example
            R = QQ[x]
            oneEntry({3},{5},x^2)
            oneEntry(,{5},x^2)
            oneEntry({3},,x^2)
        Text
            If integers are used in place of either row or col, they are multiplied by
            degree(ent).
        Example
            S = QQ[y,Degrees=>{{4,-5}}]
            oneEntry(3,,y)
    Caveat
        the ring element ent should be homogeneous
///

doc ///
    Key
        mapToTriangleRing
        (mapToTriangleRing,Ring)
    Headline
        embeds a ring of the form R[d,f]/d^2 in its triangle ring R[d,e,f]/(d^2,e^3)
    Usage
        mapToTriangleRing Q
    Inputs
        Q:Ring
            of the form R[d,f]/d^2
    Outputs
        :RingMap
            from Q to its triangle ring R[d,e,f]/(d^2,e^3)
    Description
        Text
            The external degrees of d and f double, and the external degree of the new variable "e"
            is the difference of external degrees (degree d) - (degree f).
        Example
            R = QQ[x];
            Q = R[d,f,Degrees=>{{1,2,3},{10,11,12}}]/d^2
            phi = mapToTriangleRing Q
            T = target phi
            degree \ {Q_0, Q_1}
            degree \ {T_0, T_1, T_2}
    Caveat
        The input ring Q must be valid for @ TO expectChainSequenceRing @.
    SeeAlso
        expectChainSequenceRing
///


doc ///
    Key
        toChainComplex
        (toChainComplex,Module)
    Headline
        converts a module for R[d]/d^2 to a chain complex
    Usage
        toChainComplex M
    Inputs
        M:Module
            for a ring of the form R[d]/d^2
    Outputs
        C:ChainComplex
            consisting of the the degrees of M that are multiples of (degree d),
            together with a differential reflecting the action of d.
    Description
        Text
            Suppose d has degree v.  The ouput chain complex C has C_0 = M_{0*v},
            and since the differential in a chain complex has degree -1, it has
            generally


            $C_i = M_{-iv}$.
        Example
            R = ZZ[d,Degrees=>{2}]/d^2;
            M = cokernel map(R^(-{{0},{1},{2},{3}}),,{{4,0,d,0},{0,6,0,d},{0,0,8,0},{0,0,0,10}})
            isHomogeneous M
            prune toChainComplex M
            apply(10,d->prune evaluateInDegree({d},M))
    Caveat
        M must be homogeneous
    SeeAlso
        chainModule
///

doc ///
    Key
        TorCouple
        (TorCouple,Module,List)
    Headline
        the exact couple obtained by applying Tor(W,-) to a filtered module
    Usage
        TorCouple(W,submods)
    Inputs
        W:Module
            for some ring R, giving a functor Tor(W,-)
        submods:List
            of R-modules \{A_0, A_1, ..., A_m} with each A_i inside A_{i+1}
    Outputs
        M:Module
            an exact couple
    Description
        Text
            For notational convenience, set $X = A_m$, and
            extend the sequence $A_i$ to all $i \in \ZZ$
            by setting $A_i = 0$ for $i < 0$, and $A_i = X$ for $i > m$.  
            
            The returned couple $M$ is a module for  
            the ring R[e_1,f_1,Degrees=>\{\{-1,-1},\{0,2}}].
            We describe the module $M$ in every bidegree $\{s,t}$.  The description
            depends on the parity of $s$ and $t$.
            
            
            If $s$ and $t$ are both even, say $\{s,t} = \{2p, 2q}$, then
    
            $M_{s,t} = Tor_p(W,A_q/A_{q-1})$;
            
            if $s$ and $t$ are both odd, say $\{s,t} = \{2p+1,2q+1}$, then
    
            $M_{s,t} = Tor_p(W,A_q)$;
            
            and otherwise, if $s$ and $t$ sum to an odd number, then $M_{s,t} = 0$.  
            
            The variables $e_1$ and $f_1$ act by the maps in the 
            various long exact sequences
            
            $Tor_p(W, A_{q-1}) \to Tor_p(W, A_q) \to 
             Tor_p(W, A_q / A_{q-1}) \to Tor_{p-1}(W, A_{q-1})$.
            
            {\bf Associated spectral sequence}
            
            The spectral sequence associated to this couple converges to $Tor_p(W,X)$.
            The differential on page $r$ has bidegree \{-1,-r}.  The first page has
            
            $E^{pq}_1 = Tor_p(W,A_q/A_{q-1})$.
            
            Setting $F^p_q = image(Tor_p(W,A_q) \to Tor_p(W,X))$, the infinity page has
            
            $E^{p,q}_{\infty} = F^p_{q} / F^p_{q-1}$.
        Example
            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}));
            for m in submods do print m;
            W = coker map(R^1,,{{x^3}})
            couple = prune TorCouple(W,submods)
            expectExactCouple couple
            plotPages((-1..2,-1..5,1..3), prune @@ evaluateInDegree, couple)
            A = i -> if i < 0 then image(0*id_X) else if i >= #submods then X else submods#i;
            E1 = (q,p) -> prune Tor_p(W,A(q)/A(q-1));
            netList reverse table(5,2,E1)
            inc = q -> inducedMap(X,A(q));
            filt = (p,q) -> image Tor_p(W,inc q); --no method for this?
            Einfty = (q,p) -> prune(filt(p,q)/filt(p,q-1));
            --netList reverse table(5,2,Einfty)
    SeeAlso
        covariantExtCouple
        contravariantExtCouple
///

doc ///
    Key
        TorLES
        (TorLES, ZZ, Module, Module, Module)
    Headline
        the long exact sequence in Tor induced by an inclusion in the second coordinate
    Usage
        TorLES(k,W,X,A)
    Inputs
        k:ZZ
            number of rotations of the long exact sequence to display
        W:Module
            giving a functor Tor(W,-)
        X:Module
        A:Module
            a submodule of X
    Outputs
        :Net
            showing the long exact sequence in Tor
    Description
        Text
            The long exact sequence is returned as a Net with the following general format:
        CannedExample
            |   .- ->      0
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->    W ** X    - - ->    W ** X/A   - - ->       0       - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->  Tor_1(W,X)  - - ->  Tor_1(W,X/A)  - - ->    W ** A    - -'
            |  (
            |   \
            |
            |                                                                   \
            |   .- ->  Tor_2(W,X)  - - ->  Tor_2(W,X/A)  - - ->  Tor_1(W,A)  - -'
            |  (
            |   \
            |
            |                                                                   \
            |                          (k-3) more rows appearing    ...      - -'
        Text
            The next example gives a typical use.
        Example
            R = QQ[x]
            X = R^1 / x^9
            A = image map(X,,{{x^7}})
            W = coker map(R^1,,{{x^3}})
            TorLES(3,W,X,A)
            apply(2, p -> prune Tor_p(W,X))
            apply(2, p -> prune Tor_p(W,X/A))
            apply(2, p -> prune Tor_p(W,A))
    Caveat
        For computational access to the maps in the sequence, use TorCouple instead.
    SeeAlso
        excerptLES
        TorCouple
        covariantExtCouple
        contravariantExtCouple
        covariantExtLES
        contravariantExtLES
///

doc ///
    Key
        canonicalFiltration
        (canonicalFiltration,Ring,Module)
    Headline
        filters a complex by its truncations
    Usage
        canonicalFiltration(Q,C)
    Inputs
        Q:Ring
            a ring of the form R[d,f]/d^2
        C:Module
            with an action of R[d]/d^2 considered as a subring of Q
    Outputs
        :Module
            with an action of Q; the variable f acts by including truncations
    Description
        Text
            The truncation of a chain complex is engineered to have homology only in certain degrees.
            There are several possible conventions.  We use the "coker" convention so that a complex

            ... - - -> C_{i-1} - - -> C_i - - -> C_{i+1} - - -> ...

            when truncated at position i, becomes

            ... - - -> C_{i-1} - - -> C_i - - -> coker d - - -> 0 - - -> 0 ...

            This convention has computational benefits (it does not require any gb calculation) but
            has the conceptual drawback that the natural map from the truncation into the original
            complex is not an inclusion.  At the level of the filtered derived category, however,
            all is well.  Some more details would be good, but this function is not yet being
            used for anything, so we stop here.
        Example
            R = QQ[x,y,z];
            C = chainModule res coker vars R;
            phi = map(R[d,f,Degrees=>{{1,0},{0,1}}]/d^2, ring C)
            canonicalFiltration(target phi, C)
///

doc ///
    Key
        sequenceModule
        (sequenceModule,List)
        (sequenceModule,Ring,List)
    Headline
        builds a graded R[t]-module from a sequence of maps
    Usage
        sequenceModule(Q, L)
    Inputs
        Q:Ring
            of the form R[t]; if omitted, a suitable ring is supplied
        L:List
            of maps $\{A_0 \to A_1, A_1 \to A_2, ..., A_{k-1} \to A_k\}$
    Outputs
        M:Module
            with an action of Q
    Description
        Text
            The module M is concentrated in nonnegative degrees.  For $0 \le n \le k$, we have

            $M_n \cong A_n$.

            For $n \ge k$, the module is constant

            $M_n \cong A_k$.

            The variable t acts by the supplied maps in the range $0 \le n \le k$, and by the identity
            past degree k.
        Example
            R = QQ[x]
            Q = R[t]
            m1 = random(R^{-3,-4,-5},R^{-6,-7,-8})
            m2 = random(R^{0,-1,-2},R^{-3,-4,-5})
            M = sequenceModule(Q,{m1,m2})
            isHomogeneous M
    Caveat
        If the degree of the variable t is not 1, the above still applies, but M_n should be
        interpreted as the degree $n * (deg t)$ part of M.
    SeeAlso
        filtrationModule
///

doc ///
    Key
        Page
    Headline
        for a couple ring Q, Q.Page is the page number
    Caveat
        this value must be installed by hand or by expectCoupleRing(Q)
    SeeAlso
        expectCoupleRing
///

doc ///
    Key
        isEvenDegree
    Headline
        for a couple ring Q, Q.isEvenDegree returns true on page-degrees of Q
    Caveat
        this value must be installed by hand or by expectCoupleRing(Q)
    SeeAlso
        expectCoupleRing
///

doc ///
    Key
        isOddDegree
    Headline
        for a couple ring Q, Q.isOddDegree returns true on auxiliary-degrees of Q
    Caveat
        this value must be installed by hand or by expectCoupleRing(Q)
    SeeAlso
        expectCoupleRing
///


TEST ///
S = QQ[s, t, u]; -- TODO: this won't work over ZZ
R = S[x, y, z];
-- TODO: the following two lines don't work; they always return zero
m1 = random(R^{{0,0},{-1,-1}},R^{{-2,-2},{-3,-3}});
m2 = random(R^{{-2,-2},{-3,-3}},R^{{-4,-4},{-5,-5}});
law = evaluateInDegreeLaw({9},R);
-- Check functoriality of evaluateInDegreeLaw
assert(law(m1) * law(m2) == law(m1 * m2));
m1 = random(S^{{0},{-1}},S^{{-2},{-3}});
m2 = random(S^{{-2},{-3}},S^{{-4},{-5}});
law = extensionInDegreeLaw({9},R);
-- Check functoriality of extensionInDegreeLaw
assert(law(m1) * law(m2) == law(m1 * m2));
-- TODO: assert unit/counit equations
///

TEST ///
S = QQ[s, t, u];
R = triangleRing(S, d, e, f);
m1 = oneEntry(,degree d,d);
m2 = oneEntry(degree d,,f);
n1 = map(R^{-{0,2,0}}, , {{f}});
n2 = map(source n1, , {{d}});
law = distinguishedTriangleLaw(R);
assert(law(m1) * law(m2) == law(m1 * m2))
assert(law(n1) * law(n2) == law(n1 * n2))
assert(law(oneEntry(,{0,0,0},d)) * law(oneEntry({0,0,0},,d)) == 0)
///


end--

-- Here place M2 code that you find useful while developing this
-- package.  None of it will be executed when the file is loaded,
-- because loading stops when the symbol "end" is encountered.


check ExactCouples

installPackage("ExactCouples",FileName => "/Users/jwiltshiregordon/Dropbox/Programming/Macaulay2/ExactCouples/ExactCouples.m2")
installPackage("ExactCouples",FileName => "/Users/jwiltshiregordon/Dropbox/Programming/Macaulay2/ExactCouples/ExactCouples.m2", RerunExamples=>true)

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/packages PACKAGES=ExactCouples pre-install"
-- End:


-- TODO: explain pruning lemma
-- TODO: explain algorithms for (co/contra)variantExtCouple and TorCouple
-- 
-- TODO: long exact sequence of a triple documentation example
-- TODO: map of filtered modules gives induced map on LES "functoriality" page for docs
-- TODO: better links in docs
--
-- TODO: clean up couple code
-- TODO: highlight Ext and Tor couples in docs.
-- TODO: longExactSequence probably has doc errors, and anyhow should be rewritten
-- TODO: test unit-counit formulas for the adjunctions

-- surj of exact couples has exact kernel
-- at chain level, any map of chain sequence modules has a cone, also a chain sequence module
-- if this thing has no homology, then the original map is q.i. and gives identical couples
-- so you can prove two chain sequence modules give the same couple by comparing with a map,
-- taking the cone, and proving that couple is zero.
-- What if you only want isos on the page? Proof idea relies on new lemma:
-- if A->B->C->A[1]->B[1] is an exact sequence of chain complexes, then it induces a long 
-- exact sequence in homology.

-- TODO: sequence modules allow any degrees. So... must allow throughout.  Currently failing
-- for triangle rings somewhere.

restart
needsPackage "ExactCouples"
            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}));
            for m in submods do print m;
            W = coker map(R^1,,{{x^3}})
            --Q = coupleRing(R,1,ee,ff,Degrees=>{{-1, -1}, {0, 2}})
            seqmod = filtrationModule(R[t,Degrees=>{{-2}}],submods)
            couple = prune TorCouple(ee, W, seqmod)
            --couple = prune TorCouple(W,submods)
            expectExactCouple couple
            plotPages((-1..2,-1..5,1..3), prune @@ evaluateInDegree, couple)
            
            
            
            

            R = QQ[x]
            X = R^1 / x^9
            submods = apply(5,k->image map(X,,{{x^(8-2*k)}}));
            for m in submods do print m;
            W = coker map(R^1,,{{x^3}})
            couple = prune TorCouple(W,submods)
            expectExactCouple couple
            plotPages((-1..2,-1..5,1..3), prune @@ evaluateInDegree, couple)

-- What refactoring can I do to my couple code?
-- submods -> chainModule is filtrationModule
-- safe res
-- TorCouple(Module, Module) would be ok? second module would be a sequence module
-- fm in the code.  Then one method could call the other.
-- variable t is never used by name.
-- use coupleRing so that you don't need to name the variables local e or local f
-- Use getSymbol "e" and getSymbol "f".
-- For contravariantExtCouple, second module should correspond to the cofiltration

-- Conf(2,X) is space of possible edge attachments


S = QQ[s]
A = S[a,b]
R = A[f]
declareGenerators(R, {x=>{0,1,0},y=>{1,0,0}})
M = cospan(f*x-(a+b)*y,(f*a-f*b)*y,s^3*x,s^4*y)
isHomogeneous M
xyplot := (a,b,mm)->netList reverse table(toList b,toList a,(j,i)->prune evaluateInDegree({i,j},mm));
flattenModule := m -> ((flattenRing ring m)#1) ** m;
phi = map(S[f,a,b,Degrees=>{{1,0},{0,1},{0,1}}],R,DegreeMap=>deg->deg)
xyplot(0..2,0..1,phi **  M)
M0 = evaluateInDegree({0},restackModule({2,1,3},M))
M1 = evaluateInDegree({1},restackModule({2,1,3},M))
W = 



doc ///
    Key
    Headline
    Usage
    Inputs
    Outputs
    Description
        Text
        Example
    Caveat
    SeeAlso
///

@ TO expectCoupleRing @

(see @ TO expectCoupleRing @)

-- Unused doc snippets here:
        Text
            {\bf Making the spectral sequence equivariant for an easy action}
            
            Thinking of $S^3$ as the unit vectors in the complex Hilbert space $\mathbb{C}^2$,
            the Hopf fibration becomes the map to the Riemann sphere $S^2$.
            Complex conjugation acts on this construction of the Hopf fibration.  Let us redo
            the construction to carry the action of complex conjugation.  Our cell structure on
            $S^2$ must be preserved by the action, and the basepoint must be fixed.  This is easily
            accomplished by choosing a real point to be the basepoint of $S^2 = \mathbb{C}P^1$.

            So the calculation proceeds exactly as before, except we must add relations to enforce
            the action of conjugation, which is encoded as an action of a new variable c.
        Example
            erase(symbol w); erase(symbol x); erase(symbol y); erase(symbol z);
            R = QQ[c,s,Degrees=>{1,1}]/(c^2-s^2); -- An unknown M2 issue prevents me from using s=1
            Q = R[e_1,f_1,Degrees=>{{-1,0},{2,-2}}];
            declareCouple(Q, {z => {4,0,0}}, {x => {1,0,0}, y => {1,2,0}, w => {5,2,0}})
            C = cospan(e_1*z-f_1*y,c*x-s*x,c*y+s*y,c*z+s*z,c*w-s*w)
            isHomogeneous C
            expectExactCouple C
            plotPages((-2..4,-2..3,1..4), prune @@ evaluateInDegree,C)