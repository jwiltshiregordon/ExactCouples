-- TODO: allow user to supply output ring
TorCouple = method()
TorCouple(Module, Module) := Module => (W, seqmod) -> (
    TorCouple(getSymbol "e", W, seqmod)
    )
TorCouple(Symbol, Module, Module) := Module => (eSymbol, W, seqmod) -> (
    R := ring W;
    F := ring seqmod;
    expectSequenceRing F;
    if not (coefficientRing ring seqmod) === ring W then (
        error("last two arguments incompatible: TorCouple(eSymbol, W, seqmod) " |
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
