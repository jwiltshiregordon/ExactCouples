restackRing = method()
restackRing(List, Ring) := RingMap => (fn, R) -> (
    n := #fn;
    if n == 0 then return id_R;
    m := max fn;
    R0 := R;
    filt := reverse(for i from 1 to n list R0 do R0 = coefficientRing R0);
    coefs := coefficientRing first filt;
    ringInfo := S -> (
        variables := first entries vars S;
        external := externalDegreeIndices S;
        vardegs := apply(variables, v -> (degree v)_external);
        (variables, vardegs, ideal S)
        );
    ri := apply(filt, ringInfo);
    dls := apply(ri, z -> #(first z#1));
    S0 := coefs;
    for i from 1 to m do (
        variables := {};
        vardegs := {};
        ideals := {};
        for j from 1 to n do (
            if fn#(j-1) == i then (
                {vs, vds, I} := ri#(j-1);
                variables = variables | vs;
                zedsA := if #vardegs == 0 then {} else 0 * (first vardegs);
                zedsB := if #vds == 0 then {} else 0 * (first vds);
                vardegs = apply(vardegs, deg -> deg | zedsB) | apply(vds, deg -> zedsA | deg);
                ideals = ideals | {I};
                );
            );
        S0 = S0[variables, Degrees=>vardegs];
        S0 = S0 / (sum apply(ideals, I -> sub(I, S0)));
        );
    batchcount := wts -> (
        awts := {0} | accumulate(plus, {0} | wts);
        apply(#wts, i->toList(awts#i..<(awts#(i+1))))
        );
    fiberflip := wts -> (
        flatten apply(batchcount wts, reverse)
        );
    sumfibers := (wts, surj) -> (
        apply(1 + max surj, v -> sum apply(wts, surj, (w,x) -> if x == v then w else 0))
        );
    surjcable := (wts, surj) -> (
        m := 1 + max surj;
        counts := new MutableList from m:0;
        bc := batchcount(sumfibers(wts, surj));
        cable := new MutableList from flatten apply(wts, surj, (w,v) -> toList(w:v));
        for i in 0 ..< #cable do (
            k := cable#i;
            cable#i = bc#k#(counts#k);
            counts#k = counts#k + 1;
            );
        toList cable
        );
    flipsurj := surj -> (
        pi0 := reverse toList(0..<#surj);
        pi1 := reverse toList(0..(max surj));
        pi1_(surj_pi0)
        );
    surj := flipsurj apply(fn,v->v-1);
    wts := reverse dls;
    perm0 := fiberflip(wts);
    perm1 := surjcable(wts, surj);
    perm2 := fiberflip(sumfibers(wts, surj));
    perm := inversePermutation(perm2_(perm1_perm0)) | toList(#perm0..<(degreeLength R));
    map(S0,R,DegreeMap=>(deg -> deg_perm))
    );

restackModule = method()
restackModule(List, Module) := Module => (l, M) -> tensorFlat(restackRing(l, ring M), M)

TEST ///
R = ((ZZ[a,Degrees=>{{1,2,3}}])[b,c,Degrees=>{2,4}])[d,e,f,Degrees=>{{1,1},{1,2},{2,1}}];
surjs = {{1,2,3},{1,3,2},{2,1,3},{2,3,1},{3,1,2},{3,2,1},{1,1,2},{1,2,1},{2,1,1},{1,2,2},{2,1,2},
         {2,2,1},{1,1,1},{1,2},{2,1},{1,1},{1},{}};
for surj in surjs do (
    assert(isHomogeneous restackRing(surj, R));
    );
///
