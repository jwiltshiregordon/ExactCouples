-- -*- coding: utf-8 -*-
newPackage(
  "ExactCouples",
      Version => "0.6",
      Date => "July 20, 2020",
      Authors => {
       {
       Name => "John D. Wiltshire-Gordon",
       Email => "jwiltshiregordon@gmail.com"}
       },
      HomePage => "https://jwiltshiregordon.github.io",
      Headline => "spectral sequences by Massey's method of exact couples",
  AuxiliaryFiles => true, -- set to true if package comes with auxiliary files
      DebuggingMode => true,     -- set to true only during development
  PackageImports => {"Elimination"}
      )

-- Any symbols or functions that the user is to have access to
-- must be placed in one of the following two lists
export {--"applyEntrywise",
  "declareGenerators", "cospan",
  "internalDegreeIndices", "externalDegreeIndices",
  "evaluateInDegree", "structureMap", "eid",
  "extensionInDegree",
  "evaluateInDegreeLaw", "extensionInDegreeLaw", "distinguishedTriangleLaw",
  "expectChainRing", "expectCoupleRing", "expectTriangleRing","expectSequenceRing",
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
  "chainModule", "mapToTriangleRing",
  "filteredSimplicialComplexCouple"
  --,"generateLaw"
  }
exportMutable {}

-- we need to use rawBasis
debug Core

load "./ExactCouples/Utility.m2"
load "./ExactCouples/LawUtility.m2"
load "./ExactCouples/EvaluationAdjunction.m2"
load "./ExactCouples/RingBuilders.m2"
load "./ExactCouples/Expectations.m2"
load "./ExactCouples/Conversions.m2"
load "./ExactCouples/TriangleAdjunction.m2"
load "./ExactCouples/LES.m2"
load "./ExactCouples/Restack.m2"
load "./ExactCouples/Couples.m2"
load "./ExactCouples/SpectralSequences.m2"
load "./ExactCouples/contravariantExtCouple.m2"
load "./ExactCouples/covariantExtCouple.m2"
load "./ExactCouples/TorCouple.m2"
load "./ExactCouples/SimplicialComplexes.m2"
load "./ExactCouples/Deprecated.m2"


eid = prune @@ evaluateInDegree;

beginDocumentation()

load "./ExactCouples/Documentation.m2"
load "./ExactCouples/AdjunctionTests.m2"
end--

check ExactCouples

installPackage("ExactCouples",FileName => "/Users/jwg/Dropbox/Programming/Macaulay2/ExactCouples/ExactCouples.m2")
installPackage("ExactCouples",FileName => "/Users/jwg/Dropbox/Programming/Macaulay2/ExactCouples/ExactCouples.m2", RerunExamples=>true)

-- Scratch:

-- Local Variables:
-- compile-command: "make -C $M2BUILDDIR/Macaulay2/packages PACKAGES=ExactCouples pre-install"
-- End:


-- TODO: explain pruning lemma
-- TODO: explain algorithms for (co/contra)variantExtCouple and TorCouple
--
-- TODO: long exact sequence of a triple documentation example
-- TODO: snake lemma doc using ker = Tor_1 and a couple
-- TODO: better links in docs
-- 
-- TODO: clean up couple code
-- TODO: test unit-counit formulas for the adjunctions
-- TODO: spell check docs

-- TODO: toChainComplex can sometimes return modules over inconsistent base rings?

-- TODO: restackModule is more important than restackRing.  Change docs to reflect this.

-- TODO: plotPages should default to prune @@ evaluateInDegree

-- surj of exact couples has exact kernel
-- at chain level, any map of chain sequence modules has a cone, also a chain sequence module
-- if this thing has no homology, then the original map is q.i. and gives identical couples
-- so you can prove two chain sequence modules give the same couple by comparing with a map,
-- taking the cone, and proving that couple is zero.
-- What if you only want isos on the page? Proof idea relies on new lemma:
-- if A->B->C->A[1]->B[1] is an exact sequence of chain complexes, then it induces a long
-- exact sequence in homology.

-- from original chain sequence module, every entry should be Tor, and the full page should be
-- Tor against a bimodule.  Then, this should generalize with koszul duality.

-- Comments on law utility
-- TODO: explain how this setup makes writing tests easier
-- TODO: general functoriality for adjointable laws
-- TODO: complete adjunction from degreeLaws, one entryLaw, and the unit or counit (degreeLaw doesn't
-- make sense until we are in adjoint scenario)  We can go direct to law without stopping at entryLaw.
-- choose a good size of matrix, and parallelize.

-- Generate adjoint pair from unit

-- sheaf cohomology example
-- Leray/MVSS example: Conf(4,R^2) --> Conf(2,R^2)xConf(2,R^2)
-- Cech cover of target, preimage, then refine.
--
-- 
-- First column: sum of cech of each open
-- Second column: sum of cech of each overlap...
-- so write each p-fold intersection of opens as a good union
-- and include all.
--
-- So in original good cover of X, draw poset of all overlaps.  
-- Start with an intersection-closed good open cover of X.
-- Now pick several subsets of these... partial unions.
-- So each tiny open is in at least one of these unions.
-- New simplicial cx? (,,,several original opens,,,;,,several collections,,)
-- every open in first part is in every collection.
-- 
-- 
-- For example, given a good cover of a genus 3 surface
-- we can use three open sets.
-- 





restart
needsPackage "ExactCouples"
Xfacets = {{1,2,3},{1,3,4},{1,4,5},{1,2,5},{2,3,6},{3,4,6},{4,5,6},{2,5,6},{1,6}}
Usets = {{1,2,3,4,5},{2,3,4,5,6}};
-- every vertex must appear in at least one of the Usets.
Xfaces = k -> unique flatten apply(Xfacets, f -> subsets(f, k+1));
Xdim = max apply(Xfacets, f -> #f - 1);
Xfacelist = flatten apply(Xdim + 1, Xfaces)

MVfacets = flatten for f in Xfacelist list {f|select(Usets, us->isSubset(f,us))}
MVdim = max apply(MVfacets, f -> #f - 1)
rowfilt = f -> #select(f,e->instance(e,ZZ)) - 1;
colfilt = f -> #select(f,e->not instance(e,ZZ)) - 1;
MVfaces = k -> if k == -1 then {{}} else select(unique flatten apply(MVfacets, f -> subsets(f, k+2)), (
        f -> (rowfilt(f) >= 0 and colfilt(f) >= 0)
        ));
R = ZZ[t];
chains = apply(1+MVdim,k->R^(-apply(MVfaces k, colfilt)))
omega = (a,b)->if isSubset(a,b) then (-1)^(position(b, v->not member(v,a))) * t^(colfilt(b)-colfilt(a)) else 0;
diffs = apply(MVdim-1,k->map(chains#k, chains#(k+1), matrix table(MVfaces k, MVfaces (k+1),omega)));

ch = chainComplex(diffs)
prune homology eid({0},ch)
prune homology eid({1},ch)
prune homology eid({2},ch)

sm = sequenceModule(R[D,Degrees=>{{-1}}]/D^2, (reverse diffs) | {map(R^{}, first chains, {})});
sm = sm ** (ring sm)^{{1-MVdim,0}};
smm = restackModule({2,1},sm);
print("pruning module");
M = prune restackModule({1,1},smm);
print("computing couple");
couple = prune exactCouple M;
plotPages((-1..(MVdim+1),-3..3,1..2), eid, couple)

-- given simplicial cx
-- and several subsets of variables with full union
-- build new complex with faces
-- (v_1,...,v_k;V_1,...,V_l)
-- with v_i \in V_j for all i,j
-- Then filter this space by V stuff.






R = QQ[x_1..x_4,y_1..y_4] -- must cut down number of vars.
point = ideal(x_2-x_1,x_3-x_2,x_4-x_3,y_2-y_1,y_3-y_2,y_4-y_3)
parallelogram = ideal(x_2-x_1+x_4-x_3,x_3-x_2+x_1-x_4,
                      y_2-y_1+y_4-y_3,y_3-y_2+y_1-y_4)
square = parallelogram + ideal(y_4-y_1+x_2-x_1,x_1-x_4+y_2-y_1)
rectangle = parallelogram + ideal((x_1-x_3)^2+(y_1-y_3)^2-(x_2-x_4)^2-(y_2-y_4)^2)
filt = module \ {rectangle,square,point}
k = max(regularity \ filt)
W = module ideal(x_1^k,x_2^k,x_3^k,x_4^k,y_1^k,y_2^k,y_3^k,y_4^k)
couple covariantExtCouple(W, filt)





loadPackage "PruneComplex"
facets = {{1,2},{1,3},{2,3},{2,4},{3,4}};
n = 3;
d = -1 + max apply(facets, f -> #f);
faces = k -> unique flatten apply(facets, f -> subsets(f, k));
ff = flatten apply(2 + d, faces);
nextsteps = f -> select(unique flatten select(facets, g -> isSubset(f, g)), v -> #f == 0 or v >= last f);

pe = mat -> (
    q := numrows mat;
    nxts := apply(q, r -> nextsteps unique first entries mat^{r});
    prod := apply(fold((a,b)->a**b,nxts), flatten);
    ret := {mat};
    for p in prod do (
        col := map(ZZ^q,ZZ^1,apply(p,v->{v}));
        if numcols mat == 0 or col != mat_{-1 + numcols mat} then (
            ret = ret | pe(mat | col);
            );
        );
    ret
    );

cs = m -> apply(numcols m, c -> m_{c});
prod = apply(pe(map(ZZ^n,ZZ^0,{})),cs);
pd = -1 + max apply(prod, f -> #f);
--filt = f -> #f - 1
filt = f -> (
    m = fold((a,b)->a|b,f);
    urows = unique apply(numrows m, r -> m^{r});
    #urows - 1
    );
R = ZZ[t];
print("generating prodfaces");
prodfaces = apply(1 + pd, k -> select(prod, f -> #f - 1 == k));
chains = apply(1+pd,k->R^(-apply(prodfaces#k, filt)));
omega = (a,b)->if isSubset(a,b) then (-1)^(position(b, v->not member(v,a))) * t^(filt(b)-filt(a)) else 0;
print("computing differentials");
diffs = apply(pd,k->map(chains#k, chains#(k+1), matrix table(prodfaces#k,prodfaces#(k+1),omega)));
cx = chainComplex diffs;
sm = sequenceModule(R[D,Degrees=>{{-1}}]/D^2,(reverse diffs) | {map(R^{}, first chains, {})});
sm = sm ** (ring sm)^{{-pd,0}};
smm = restackModule({2,1},sm);
print("pruning module");
M = prune restackModule({1,1},smm);
print("computing couple");
couple = prune exactCouple M;
plotPages((-1..(pd+1),-1..n,1..2), prune @@ evaluateInDegree, couple)




n = 4
ordring = (n,x) -> (
    ijs = apply(flatten(apply(subsets(1..n,2),ij->{ij,reverse ij})), toSequence);
    xs = apply(ijs,ij->x_ij);
    R = ZZ[xs];
    I = ideal(flatten(apply(ijs,ij->{((x_ij)_R)^2,(x_ij)_R * (x_(reverse ij))_R})));
    R/I
    )

ordring(4,u)



us = flatten(apply(subsets(1..n,2),ij->{u_ij,u_(reverse ij)}))
vs = flatten(apply(subsets(1..n,2),ij->{v_ij,v_(reverse ij)}))
cov = ZZ[us|vs]
cov/ideal(apply(us|vs,x->(x_cov)^2) | apply(subsets(1..n,2),ij->u_ij * u_(reverse ij)) | 











eid = prune @@ evaluateInDegree
face2 = {{1,2,4},{1,3,4},{3,4,6},{3,5,6},{2,4,8},{4,6,8},{5,6,7},{6,7,8},
                      {1,2,8},{1,7,8},{1,3,7},{3,5,7}};
face1 = unique(flatten(apply(subsets({0,1,2},2), s-> apply(face2, f->f_s))));
face0 = unique(flatten(apply(subsets({0,1,2},1), s-> apply(face2, f->f_s))));
diags = {{1,4},{3,6},{4,8},{6,7},{1,8},{3,7}};
skel = f -> if #f != 2 then #f-1 else if member(f,diags) then 2 else 1;
R = ZZ[t]
free0 = R^(-apply(face0,skel))
free1 = R^(-apply(face1,skel))
free2 = R^(-apply(face2,skel))
omega = (a,b)->if isSubset(a,b) then (-1)^(position(b, v->not member(v,a))) * t^(skel(b)-skel(a)) else 0
d12 = map(free1, free2, matrix table(face1,face2,omega))
d01 = map(free0, free1, matrix table(face0,face1,omega))
sm = sequenceModule(R[d,Degrees=>{{-1}}]/d^2,{d12,d01,map(R^{},free0,{})})
sm = sm ** (ring sm)^{{-2,0}}
smm = restackModule({2,1},sm)
M = restackModule({1,1},smm)
couple = prune exactCouple M
plotPages((-1..3,-1..3,1..2), prune @@ evaluateInDegree, couple)
prune pageModule(1,D,couple)


R = ZZ[d,t,Degrees=>{{1,0},{0,1}}]
free = R^(-apply(faces, f->{#f-1,skel(f)}))



apply(face2 | face1 | face0, skel)

apply(facets, f-> d




R = QQ[x,y,z];
p = y^2*z-x^3+17*z^3;
X = Proj(R/p)


R = (ZZ/11)[x,y,z];
p = y^2*z-x^3+17*z^3;
F = sheaf(R^1/p)
apply(3, k-> HH^k F)


            

            R = ZZ[x,y,z];
            p = y^2*z-x^3+17*z^3;
            filt = {module ideal(9*p), module ideal(3*p), module ideal(p)};
            k = max({0} | apply(filt,regularity));
            W = module ideal(x^k,y^k,z^k);
            couple = prune covariantExtCouple(W, filt)
            couple' = prune evaluateInDegree({0},restackModule({2,1}, couple));
            plotPages((-1..3,-1..3,1..2),prune@@evaluateInDegree,couple');


R = ZZ[x,y,z]
p = y^2*z-x^3+17*z^3
q = y^2*z-x^3-2*x*z^2
filt = {module ideal p, module ideal(p,q), R^1}
for d from -3 to 3 do (
    print("d = " | toString(d));
    filtd = apply(filt,m->m**R^({d}));
    k = max({0} | apply(filtd,regularity));
    use R;
    W = module ideal(x^k,y^k,z^k);
    couple = prune covariantExtCouple(W, filtd);
    couple' = prune evaluateInDegree({0},restackModule({2,1}, couple));
    plotPages((-1..3,-1..3,1..2),prune@@evaluateInDegree,couple');
    );



Q=ZZ[r,Degrees => {{1,4,5}}];
A=Q[x,y, Degrees => {{1,2},{1,2}}]/(x^2+y^2);
B=A[b];
C=B[p,q]/(p^3-2*q^3);
D=C[d];
isHomogeneous restackRing({1,2,3,4},D)
isHomogeneous restackRing({2,1,4,3},D)
isHomogeneous restackRing({1,1,1,1},D)


-- D = QQ[x,y][b][p,q][d]
-- flatten would be QQ[x,y,b,p,q,d]
-- and this would be the function 1111
-- what would be 1122?
-- QQ[x,y,b][p,q,d]
-- what would be 2211?
-- QQ[p,q,d][x,y,b]
-- 2231 should give QQ[d][x,y,b][p,q]




wts = {2,1,1,1}
surj = {1,1,2,0}
surjcable(wts,flipsurj(surj))


netList apply(6,i->{R_i,degree R_i})
netList apply(6,i->{S0_i,degree S0_i})

-- Record original degree vector bunches.
-- Start at end.  1 for d: {0}
-- then 1 for p,q: {1}
-- then 1 for b: {2}
-- then 2 for x,y: {3,4}
-- 




variables = apply(filt, S -> first entries vars S);
exdegs = apply(filt, externalDegreeIndices);
vardegs = apply(filt, variables, exdegs, (S, V, E) -> apply(V, v -> v_E));

apply(filt, S -> apply(first entries vars S, v -> (degree v)_(externalDegreeIndices S)))
apply    






restackRing({2,3,4,1}, D)

-- find coefficients underneath
-- can sub() ideal into whatever ring.  






R = (ZZ/3)[x,y,z]
p = x^4+y^4+z^4
q = p*(x^3+y^3)
filt = {module ideal p, module ideal(p,q), R^1}
for d from -3 to 3 do (
    print("d = " | toString(d));
    filtd = apply(filt,m->m**R^({d}));
    k = max({0} | apply(filtd,regularity));
    use R;
    W = module ideal(x^k,y^k,z^k);
    couple = prune covariantExtCouple(W, filtd);
    couple' = prune evaluateInDegree({0},restackModule({2,1}, couple));
    plotPages((-1..3,-1..3,1..2),prune@@evaluateInDegree,couple');
    );





excerptLES(0,3,couple')

plotPages


regularity M



restart
needsPackage "ExactCouples"
needsPackage "PushForward"
R = QQ[x,y, Degrees=>{{1,0},{0,1}}]
m = matrix {{1_R,x*y^3,x^2*y^8,x^3*y,x^4*y^2,x^5*y^4,x^6*y^7,x^7*y^5,x^8*y^6}}
sparegens = M -> image map(M,, (prune coker cover (prune M).cache.pruningMap).cache.pruningMap)
M = image m
viennot = reverse for i from 0 to 4 list M do M = sparegens M
couple = prune TorCouple(coker vars R, viennot)
disp = (deg,E) -> degrees minimalPresentation evaluateInDegree(deg,E);
plotPages((-1..3,-1..6,1..5), disp, couple);



M = prune image m
pF = prune cokernel cover M.cache.pruningMap
image map(image M, , pF.cache.pruningMap)
image map(image m,, map(cover M, , cover pM.cache.pruningMap))


restart
needsPackage "ExactCouples"
needsPackage "PushForward"
needsPackage "SpechtModule"
lambda = {3,2}
R = QQ[x,y]
I = ideal apply(lambda | {0}, 0..#lambda, (k,l)->x^l*y^k)
boxes = flatten apply(#lambda,r->apply(lambda#r,c->(r,c)))
p = new Partition from lambda
tableaux = toListOfTableaux(standardTableaux p)
filts = apply(tableaux, t -> reverse apply(1 + sum lambda, k -> 
        module(I + ideal apply(boxes, b -> (if t_b >= k then x^(b#0)*y^(b#1) else 0_R)))))
for i in 0..<#tableaux do (
    print(tableaux#i);
    filt = filts#i;
    couple = prune TorCouple(coker vars R, filt);
    --disp = (deg,E) -> prune pushFwd(map(R,QQ),evaluateInDegree(deg,E))
    disp = (deg,E) -> sort flatten degrees minimalPresentation evaluateInDegree(deg,E);
    plotPages((-1..3,-1..6,1..5), disp, couple);
    )



restart
needsPackage "ExactCouples"
needsPackage "PushForward"
R = QQ[x,y]
ideals = {ideal(0_R), ideal(x^4,y^4,x^2*y), ideal(x^2,y^2,x^2*y),ideal(x,y)}
M = coker map(R^{0,-1},, {{x^3,y^4,0,0,x^2,y^2},{0,0,x^2,y^2,-x-y,2*y+7*x}})
quos = apply(ideals, I->M /(I*M))
inds = apply(-1+#ideals,k->inducedMap(quos#(k+1),quos#(k)))
seq = sequenceModule(R[t],inds)
couple = prune TorCouple(coker vars R, seq)
--disp = (deg,E) -> prune pushFwd(map(R,QQ),evaluateInDegree(deg,E))
disp = (deg,E) -> sort flatten degrees minimalPresentation evaluateInDegree(deg,E)
plotPages((-1..4,-1..4,1..5), disp, couple)





restart
needsPackage "ExactCouples"
needsPackage "PushForward"
R = QQ[x,y]
lambdas = {{1},{4,1,1,1},{4,4,2,1},{4,4,4,3}}
f = lambda -> module ideal apply(lambda | {0}, 0..#lambda, (k,l)->x^k*y^l)
filt = reverse apply(lambdas,f)
couple = prune TorCouple(coker vars R, filt)
--disp = (deg,E) -> prune pushFwd(map(R,QQ),evaluateInDegree(deg,E))
disp = (deg,E) -> sort flatten degrees minimalPresentation evaluateInDegree(deg,E)
plotPages((-1..3,-1..4,1..4), disp, couple)






restart
needsPackage "ExactCouples"
needsPackage "PushForward"
R = QQ[x,y,z,t]
m = matrix {{t^2,t*x,x^2},{t^2,t*y,y^2},{t^2,t*z,z^2}}
filt = reverse apply({1,2,3},i->module minors(i,m))
couple = prune TorCouple(coker matrix {{t}}, filt)
disp = (deg,E) -> prune pushFwd(map(R,QQ[x,y,z]),evaluateInDegree(deg,E))
plotPages((-1..2,-1..4,1..2), disp, couple)



restart
needsPackage "ExactCouples"
needsPackage "PushForward"
R = QQ[x,y,z,t]
m = matrix {{x,y,z},{y,x,t},{z,t,x}}
I = minors(2,m)
filt = reverse apply(3,k->module intersect(I,ideal(x^k)))
couple = prune TorCouple(coker vars R, filt)
disp = (deg,E) -> prune pushFwd(map(R,QQ),evaluateInDegree(deg,E))
plotPages((-1..4,-1..3,1..2), disp, couple)




S = R / det(m)
M = coker(sub(m,S) - t * id_(S^3))
prune Hom(coker S_3, M)

Is = reverse apply(4,i->module minors(i,sub(m,S)))
couple = prune covariantExtCouple(M ** coker t, Is)
plotPages((0..4,0..4,1..3),prune @@ evaluateInDegree, couple)





R = QQ[x,y,z]
m = matrix {{-y,-z,0},{x,0,-z},{0,x,y}}
Is = reverse apply(3,i->module minors(i,m))
couple = prune covariantExtCouple(coker vars R, Is)
plotPages((0..4,0..4,1..3),prune @@ evaluateInDegree, couple)


minors(0,p)
minors(1,p)
minors(2,p)
minors(3,p)
minors(4,p)

W = coker matrix {{t}}
prune Hom(W,M)
prune Ext^1(W,M)


M = coker(m - t * id_(R^3))























R = ZZ/101[a..d];
I = ideal(a*b-c*d, (a*c-b*d)^2);
pd = primaryDecomposition I
fpd = apply(reverse accumulate(intersect, pd), module)
couple = prune contravariantExtCouple(fpd,R^1)
plotPages((0..5,-3..3,1..5), prune @@ evaluateInDegree, couple)




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









needsPackage "NoetherNormalization"
loadPackage "Dmodules"
needsPackage "TestIdeals"
needsPackage "PushForward"

R=QQ[a,b,c,d,e,f]
I = ideal(a*c-1,b*e-1,d*f-1,a*e*d-1,b*c*f-1)
isCohenMacaulay(R/I)
M = R^1/I
(aut, J, vv) = noetherNormalization I
S = QQ[vv]
phi = aut^(-1) * map(R,S)
acts = apply(6,i->pushFwd(phi,map(M,M,matrix {{R_i}})))
for a in acts do (print(a); print(""))
(A,B,C,D,E,F)=toSequence(acts)
gd = m -> matrix(degree \ (first entries transpose m_{0}))
gd(A^50)


R=QQ[x,y]
I=ideal(x*y-1)
M=R^1/I
(aut,J,vv)=noetherNormalization I
S=QQ[vv]
phi=aut^(-1) * map(R,S)
apply(2,i->pushFwd(phi,map(M,M,matrix {{R_i}})))




R=QQ[x,y]
I=ideal(x^2-y^2)
S=QQ[x]
phi=map(R,S)
isCohenMacaulay(R/I)
unit=lift(basis(R/(I+ideal(R_0))),R)
k=numcols unit

--find action of each var, starting with R_0
cosetRep=((matrix {{R_0}})*unit)//unit
gbU = groebnerBasis((kernel unit)/I)
Q=QQ[q,x,y]
qs=matrix {apply(k,i->q^i)}
enc = first entries (qs*(map(Q,R)**cosetRep))
egb = qs*(map(Q,R)**gbU)
enm = matrix {apply(enc,r->first first entries gens eliminate({y},ideal(r)+ideal(egb)+ideal(map(Q,R)**(gens I))))}


map(Q,R)**cosetRep

J=ideal(x*y,y-q,x^2-y^2)
eliminate({y},J)
J=ideal(x,y-q,x^2-y^2)
eliminate({y},J)

uu = (unit || (0*unit))
(x*uu) % gbU
(y*uu) % gbU



R=QQ[a,b,c,d,e,f,t]
I = ideal(a*c-t^2,b*e-t^2,d*f-t^2,a*e*d-t^3,b*c*f-t^3)
isCohenMacaulay(R/I)
(phi, J, vv) = noetherNormalization I
S=QQ[vv]
psi = map(R,S)
pushForward(psi,R^1/J)



R=QQ[x,y]
noetherNormalization ideal(x*y-1)

restart
S=QQ[x,y,z,t]
M=S^1/ideal(x*y-t^2,x-2*y-z)
R=QQ[z,t]
phi=map(S,R)
pushForward(phi,M)



W = QQ[X, dX, WeylAlgebra=>{X=>dX}]
M = W^1/(dX)
L = Dlocalize(M,X)
H = DHom(M,M)
H_0












W = QQ[X, dX, Y, dY, WeylAlgebra=>{X=>dX, Y=>dY}]
M = W^1/(X*Y-1,X*dX+2)
isHolonomic M
apply(10,k->hilbertFunction(k,M))







declareGenerators(W,{g=>{0},h=>{0}})
M = cospan(Y^2*g-h,X^2*h-g,Y*g-X*h,(Y*dX+X*dY)*g,(Y*dX+X*dY)*h)
(g,h)=(M_0,M_1)

apply(10,k->hilbertFunction(k,M))




--M = coker matrix {{X*Y-1,0},{0,X*Y-1}}
M = W^1/(X*Y-1,dX*Y+dY*X)
isHolonomic M
apply(10,k->hilbertFunction(k,M))



M = W^1/ideal(X*Y-1,dX,dY)
isHolonomic M


I = ideal (X*Y-1,dX*Y+Y^2,dY*X+X^2)


W = QQ[X, dX, WeylAlgebra=>{X=>dX}]
M = W^1/(ideal(dX))
Dlocalize(M, X)
DlocalizeMap(M, X)


presentation I




W = QQ[X, dX, Y, dY, WeylAlgebra=>{X=>dX, Y=>dY}]
I = ideal (X*Y-1)
h = localCohom(I, W^1 / ideal{dX,dY})







R = QQ[z][g,Degrees=>{-1}]/g^100
declareGenerators(R,{x=>{0,5},y=>{-1,0}})
X = cospan(z^6*x,z^10*y,g*x-z^5*y)
isHomogeneous X
declareGenerators(coefficientRing R,{w=>{0}})
Y = cospan(z^10*w)
Y' = extensionInDegree({99},R,Y)
eid = prune @@ evaluateInDegree

(prune Ext^0(eid({0},X),Y),eid({0},Ext^0(X,Y')))
(prune Ext^0(eid({-1},X),Y),eid({1},Ext^0(X,Y')))

(prune Ext^1(eid({0},X),Y),eid({0},Ext^1(X,Y')))
(prune Ext^1(eid({-1},X),Y),eid({1},Ext^1(X,Y')))

contravariantExtCouple(X,Y')





            R = QQ[z]; S = R[g][t, Degrees=>{{-1}}]; declareGenerators(S,{x=>{0,0,5},y=>{0,1,0}});
            M = cospan(z^6*x,z^3*t*x,z^10*y,z^7*t*y,g*x-z^5*y,t^2*x,t^2*y); isHomogeneous M

            eid = prune @@ evaluateInDegree; (dt, dg) = degree \ (S_0, S_1);
            {A,X,B,Y} = (deg -> prune eid({deg#1},eid({deg#0},M))) \ ({0,0},dt,dg,dt+dg);
            netList {{A, X}, {B, Y}}
            Y = R^1 / (R_0^10);
            Y' = extensionInDegree({0}, coefficientRing S, Y)
            couple = prune contravariantExtCouple(M,Y')
            expectExactCouple couple
            excerptCouple({-2,0},4,couple)

phi = map(R[g,t,Degrees=>{{0,1},{-1,0}}],S)
xyplot := (a,b,mm)->netList reverse table(toList b,toList a,(j,i)->prune evaluateInDegree({i,j},mm));
xyplot(-2..2,-2..2,phi**M)
xyplot(-5..5,-5..5,couple)

-- covariantExtCouple functoriality
            R = QQ[z]; S = R[g][t]; declareGenerators(S,{a=>{0,0,3},x=>{1,0,1},b=>{0,1,2},y=>{1,1,0}});
            M = cospan(z^13*a,z^15*x,z^6*b,z^8*y,g*a-z*b,g*x-z*y,t*a-z^2*x,t*b-z^2*y); isHomogeneous M
            eid = prune @@ evaluateInDegree; (dt, dg) = degree \ (S_0, S_1);
            {A,X,B,Y} = (deg -> prune eid({deg#1},eid({deg#0},M))) \ ({0,0},dt,dg,dt+dg);
            netList {{A, X}, {B, Y}}
            W = R^1 / (R_0^10);
            W' = extensionInDegree({0}, coefficientRing S, W)
            couple = prune covariantExtCouple(W',M)
            expectExactCouple couple










R = QQ[x]
S = R[t][g]
declareGenerators(S,{a=>{0,0,3},b=>{0,1,1},c=>{1,0,2},d=>{1,1,0}})
M = cospan(x^13*a,x^15*b,x^6*c,x^8*d,g*a-x*c,g*b-x*d,t*a-x^2*b,t*c-x^2*d)
isHomogeneous M
{A,B,C,D} = (deg -> prune evaluateInDegree({deg#0},evaluateInDegree({deg#1},M))) \ ({0,0},{1,0},{0,1},{1,1})
W = R^1/(R_0^10)
AB = prune evaluateInDegree({0},M)
coupleAB = covariantExtCouple(W,AB)
CD = prune evaluateInDegree({1},M)
coupleCD = covariantExtCouple(W,CD)
--M' = restackModule({1,3,2},M) -- fails. Finally time to fix restackings
-- gives non-homogenoeus; hand-prepare below
phi = map(R[g][t],S,DegreeMap=>deg->deg_{1,0,2})
M' = phi ** M
isHomogeneous M'
couple = covariantExtCouple(extensionInDegree({0},coefficientRing ring M',W),M')
excerptCouple({-2,2},4,coupleAB)
excerptCouple({-2,2},4,coupleCD)
excerptCouple({-2,2},4,couple)








couple' = restackModule({1,3,2},couple)
coupleAB' = evaluateInDegree({0},couple')
excerptCouple({-2,2},4,coupleAB')
excerptCouple({-2,2},4,coupleAB)


AB = prune evaluateInDegree({0},restackModule({1,3,2},M))
CD = prune evaluateInDegree({1},restackModule({1,3,2},M))

ABcouple = covariantExtCouple(W,map(R[t],ring AB)**AB)
excerptCouple({0,2},2,ABcouple)



table({1,2,3},{1,2,3},(x,y)->x*y)


phi = map((QQ[x])[g,t,Degrees=>{{0,1},{1,0}}],R)
xyplot := (a,b,mm)->netList reverse table(toList b,toList a,(j,i)->prune evaluateInDegree({i,j},mm));
xyplot(0..1,0..1,phi**M)

W' = extensionInDegree({0},coefficientRing R,W)
Ext^1(W',evaluateInDegree({0},M))




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
