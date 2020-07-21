
-- TODO: let user supply output ring
-- algorithm:
-- Left Kan extend Y to the sequence ring, and then right Kan to chain ring.
-- call the resulting thing Y'.
-- Take the cofiltration of ({0} | submods), and resolve by frees.
-- Call the resulting chain sequence module M.
-- claim: exactCouple Hom(M, Y') is the correct couple (after truncating)
-- Pf sketch: Hom(M, Y') = Hom(M, Ran Lan Y) = Hom(Res M, Lan Y).
-- Since M is a resolution, Res M is free, say Res M = Lan G.
-- Hom(Lan G, Lan Y) = Hom(G, Res Lan Y).  And this thing is
-- the chain sequence module that defines the couple.
-- It would be mathematically similar to start with the cofiltration instead;
-- this could be faster at least some of the time.
contravariantExtCouple = method()
contravariantExtCouple(Module, Module) := (seqmod, Y) -> (
    contravariantExtCouple(getSymbol "e", seqmod, Y)
    )

contravariantExtCouple(Symbol, Module, Module) := (eSymbol, seqmod, Y) -> (
    R := ring Y;
    if not (coefficientRing ring seqmod) === ring Y then (
        error("last two arguments incompatible: contravariantExtCouple(eSymbol, seqmod, Y) " +
              "requires (coefficientRing ring seqmod) === ring Y");
        );
    F := ring seqmod;
    expectSequenceRing F;
    f := F_0;
    external := externalDegreeIndices F;
    degf := (degree f)_external;
    flattenModule := m -> ((flattenRing ring m)#1) ** m;
    fm := flattenModule seqmod;
    d := local d;
    ringfm := ring fm;
    ch := ringfm[d]/d^2;
    rfm := flattenModule chainModule(ch, res fm);
    S := R[d,f,Degrees=>{{1}|(0*degf),{0}|degf}]/d^2;
    phi := map(ring rfm, ring Y, DegreeMap=>(deg->{0}|(0*degf)|deg));
    Y' := phi ** Y;
    hm := Hom(rfm,Y');
    pres := presentation hm;
    -- Must remove certain rows and columns from pres
    -- because our chain sequence module is only correct in certain degrees;
    -- we had truncated to keep it fg
    if (degree f) == 0*(degree f) then (
        error("contravariantExtCouple(eSymbol, seqmod, Y) relies on a nonzero degree for " +
              "(ring seqmod)_0");
        );
    edF := externalDegreeIndices F;
    eheft := (heft F)_edF;
    poscheck := deg -> (sum(apply(eheft, (deg_{1..<#deg})_edF, (p,q)->p*q)) < 0);
    dtp := degrees target pres;
    rowselect := select(#dtp,k->poscheck(dtp#k));
    dsp := degrees source pres;
    colselect := select(#dsp,k->poscheck(dsp#k));
    hm = coker(pres_colselect^rowselect);
    couple := exactCouple((map(S,ring hm)) ** hm);
    Q := ring couple;
    sh := Q^{4*degree(Q_0)+1*degree(Q_1)};
    couple ** sh
    )

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
    seqmod := sequenceModule(F, apply(l, k -> inducedMap(submods'#(k+1), submods'#k)));
    contravariantExtCouple(getSymbol "e", seqmod, Y)
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
