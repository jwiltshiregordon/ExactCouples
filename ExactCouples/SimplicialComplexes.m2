filteredSimplicialComplexCouple = method()
filteredSimplicialComplexCouple(List, Function) := Module => (facets, filt) -> (
    Xdim := -1 + max apply(facets, f -> #f);
    faces := apply(1 + Xdim, k -> unique flatten apply(facets, f -> subsets(f, k + 1)));
    t := getSymbol "t";
    R := ZZ[t];
    omega := (a,b)->if isSubset(a,b) then 
                      (-1)^(position(b, v->not member(v,a))) * (R_0)^(filt(b)-filt(a)) else 0;
    frees := apply(faces, fs -> R^(-apply(fs, filt)));
    diffs := apply(Xdim, i -> map(frees#i, frees#(i+1), matrix table(faces#i, faces#(i+1), omega)));
    diffs = {map(R^{},frees#0,{})} | diffs;
    d := getSymbol "d";
    S := R[d,Degrees=>{{-1}}];
    sm := sequenceModule(S/(S_0)^2, reverse diffs);
    sm = sm ** (ring sm)^{{-Xdim,0}};
    smm := restackModule({2,1},sm);
    M := restackModule({1,1},smm);
    couple := prune exactCouple M;
    return couple;
    )
