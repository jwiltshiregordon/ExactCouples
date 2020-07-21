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
