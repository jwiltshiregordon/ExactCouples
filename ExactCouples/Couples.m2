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
    Q := S[e_1,f_1, Degrees => {dege, degf}];
    exactCouple(Q, M)
    )

-- no way to omit Q argument since it contains page information
-- exactCouple(coupleRing(?), M) is pretty easy, though
exactCouple(Ring, Module) := Module => (Q, M) -> (
    expectCoupleRing Q;
    R := ring M;
    expectChainSequenceRing R;
    T := target mapToTriangleRing R;
    title := "computing exactCouple(Q, M) where M is an S[d,f]-module: ";
    if degree Q_0 != degree T_1 then (
	error(title | "degree of Q_0 should be (degree d) - (degree f)");
	);

    if degree Q_1 != degree T_2 then (
	error(title | "degree of Q_1 should be 2 * (degree f)");
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
    if n < 0 then error("cannot form derived couples a negative number of times n=" | toString(n));
    if n == 0 then M else derivedCouple derivedCouple(n - 1, M)
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
