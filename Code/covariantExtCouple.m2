-- TODO: allow user to supply output ring
covariantExtCouple = method()
covariantExtCouple(Module, Module) := Module => (W, seqmod) -> (
    covariantExtCouple(getSymbol "e", W, seqmod)
    )

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
