
doc ///
    Key
        ExactCouples
    Headline
        spectral sequences by Massey's method of exact couples
    Description
        Text
            An exact couple is a pair of R-modules E and A together with maps A --> E --> A and A --> A
            with the conditions that

            im(A --> E) = ker(E --> A)

            im(E --> A) = ker(A --> A)

            im(A --> A) = ker(A --> E).
            
            Every exact couple gives rise to a spectral sequence, and most spectral sequences
            can be constructed by means of an exact couple.  Exact couples were introduced in
            
            Massey, W. (1952). Exact Couples in Algebraic Topology (Parts I and II). 
            Annals of Mathematics, 56(2), second series, 363-396.
            
            This package constructs exact couples arising from a bounded chain
            complex of finitely generated modules $C_*$ and an endomorphism $f : C_* \to C_*$.  Then, 
            we take A to be the homology of $C_*$, and E to be the homology of the mapping cone
            of $f$.  The maps of the exact couple then come from the usual long exact sequence.  
            The resulting spectral sequence converges to the
            homology of $C_*$.
            For example, several standard exact sequences for Tor and Ext are of this form; see 
            @ TO "Exact couples for Tor and Ext" @.
            
            Since an exact couple of R-modules is itself an R-module with some extra action maps, it
            can be considered a module for a larger ring with new variables acting by these extra maps.
            This observation allows us to encode exact couples as modules
            and thereby compute with them using Groebner bases and other standard Macaulay2 methods.
            
            This encoding strategy works generally for commuting diagrams; 
            see @ TO "Encoding diagrams as modules" @.  This style may feel unfamilar at first.  One 
            benefit comes in our approach to functoriality; see @ TO "Functoriality for Tor and Ext couples" @.
            
            If this encoding makes sense, then @ TO "Conventions and first examples" @ is a good place
            to look.
            
            {\bf A common complaint}
            
            It may seem strange or useless to compute a spectral sequence converging to
            $H_* C_*$, since in our context, it would be much easier to compute these modules directly.
            However, it is a misconception that the purpose of a spectral sequence
            is to compute its abutment.
            
            Rather, a spectral sequence is a tool for proving, and this proving is easier if you can look at 
            examples.
            --For instance, a skillful filtration of $C_*$ could clarify its
            --homology.  Or perhaps the map $f$ is
            --the mathematical focus, and so the intricacies of the spectral sequence tell you about
            --the way iterates of $f$ interact with homology.  In either case, a spectral sequence can 
            --lead to direct, understandable proofs.
            
            For the serious user of this package, there is at least one common mathematical difficulty:
            
            Even if a spectral sequence with finitely generated entries can be constructed
            by means of a couple, the auxilliary data is often not finitely generated.  The author
            of this package would be interested in any systematic method of replacing such couples
            with finitely generated versions that recover the same spectral sequence.  If you
            have information or ideas, please send an email!
           
    SeeAlso
        "Elemetary introduction: solving linear equations in abelian groups"
        "Encoding diagrams as modules"
        "Conventions and first examples"
        "Bockstein spectral sequence"
        "Serre spectral sequence in homology"
        "Exact couples for Tor and Ext"
        "Functoriality for Tor and Ext couples"
        "An exact couple associated to a Young tableau"
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
            an external R-degree
    Outputs
        :Module
            the part of M sitting in degree deg, considered as a module for
            the coefficient ring of R
    Description
        Text
            An "external" degree deg is one so that, give a degree deg' for the coefficient ring of
            R, the concatenation (deg | deg') is a valid degree for R.  So an external degree is a list
            of integers of length (degreeLength R) - (degreeLength coefficientRing R).
        Example
            S = QQ[s, t, u]; R = S[x, y]; m = matrix {{s*x^2+t*x*y+u*y^2}}; M = coker m
            N = evaluateInDegree({4}, M)
            apply(10, i -> hilbertFunction({i}, N))
            (F, f) = flattenRing R; apply(10, i -> hilbertFunction({4, i}, f ** M))
    SeeAlso
        structureMap
        extensionInDegree
///

doc ///
    Key
        structureMap
        (structureMap,List,List,RingElement,Module)
        (structureMap,List,Nothing,RingElement,Module)
        (structureMap,Nothing,List,RingElement,Module)
    Headline
        computes the action of a ring element on a particular degree
    Usage
        structureMap(x,y,r,M)
    Inputs
        M:Module
           over some ring R with coefficient ring $k$
        x:List
            an external R-degree
        y:List
            an external R-degree
        r:RingElement
            of degree y-x
    Outputs
        :Matrix
            the $k$-linear map $M_x \to M_y$, $m \mapsto r \cdot m$
    Description
        Text
            An "external" R-degree is a list of integers
            of length equal to (degreeLength R) - (degreeLength k).
            See @ TO "Encoding diagrams as modules" @
            for more discussion about internal and external degrees for a ring.
            
            A graded R-module has a component $R_x$ in every external degree $x$, and this
            component is a $k$-module.
            If $r \in R$ has degree
            $y-x$, then multiplication by $r$ gives a $k$-linear map $R_x \to R_y$.  This is the
            map that is returned by structureMap.
            
            Only one of $x$ and $y$ needs to be supplied if $r$ is nonzero, since then the other
            degree can be inferred.
        Example
            k = QQ[s, t, u]; R = k[x, y]; m = matrix {{s*x^2+t*x*y+u*y^2}}; M = coker m
            phi = structureMap({4}, {7}, x^2*y, M)
            source phi
            target phi
    SeeAlso
        evaluateInDegree
        extensionInDegree
///


-- document {
--     Key => {applyEntrywise, (applyEntrywise, Ring, FunctionClosure, FunctionClosure, Matrix)},
--     Headline => "apply a matrix-valued function to the entries of a matrix",
--     Usage => "applyEntrywise(S, degreeLaw, entryLaw, m)",
--     Inputs => {
--   "m" => Matrix => {"a homogeneous matrix over some ring R"},
--   "S" => Ring => {"the output ring"},
--   "degreeLaw" => FunctionClosure => {"a function that converts a multidegree for R ",
--                                      "to a list of multidegrees for S"},
--   "entryLaw" => FunctionClosure => {"a function that accepts a homogeneous one-by-one matrix over R and ",
--                                     "returns a matrix over S whose degrees match those determined by degreeLaw.",
--     "This function will be called on every one-by-one submatrix of m"}
--
--        },
--     Outputs => {"the block matrix obtained by replacing each entry e with the matrix entryLaw(e)"}
--     }


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
        longExactSequence
        (longExactSequence,Matrix)
        (longExactSequence,Ring,Matrix)
    Headline
        finds the long exact sequence associated to a map of R[d]/d^2-modules
    Usage
        M = longExactSequence m
    Inputs
        m:Matrix
            whose ring is of the form R[d]/d^2 for some coefficient ring R
    Outputs
        M:Module
            for a ring of the form R[e,f,Degrees=>\{\{1,0},\{-2,2}}]
    Description
        Text
            Suppose $m:A \to B$ where $A$ and $B$ are both $R[d]/d^2$-modules considered as cochain complexes.
            Writing C(m) for the mapping cone of the map m, the long
            exact sequence in cohomology takes the form

            $\cdots \to H^p A \to H^p B \to H^p C(m) \to H^{p+1} A \to \cdots$

            The output module M is bigraded, and encodes this sequence as follows:

            $\cdots \to M_{\{1,2p\}} \to M_{\{-1,2p+2\}} \to M_{\{0,2p+2\}} \to M_{\{1,2p+2\}} \to \cdots $

            where the maps are multiplication by $e$ or $f$, in the repeating pattern $...,e,f,e,e,f,e,e,f,e...$.
            In other words, for all $p \in \ZZ$,

            $M_{\{1,2p\}} = H^p A$

            $M_{\{-1,2p\}} = H^{p-1} B$

            $M_{\{0,2p\}} = H^{p-1} B$;

            Multiplication by e induces the maps $H^{p-1} B \to H^{p-1} C(m) \to H^{p}A$ and

            Multiplication by f induces diagonal maps $H^p A \to H^p B$.

            Note that this setup leaves the odd rows empty.  The reasoning behind these conventions
            is explained in @ TO "Conventions and first examples" @.
        Example
            R = QQ[x]; S = R[d] / ideal(d^2); declareGenerators(S, {a => {0,0}}); A = cospan(x^2*a, d*x*a)
            declareGenerators(S, {b => {0,0}}); B = cospan(x^2*b, d*b)
            m = map(B, A, matrix {b});
            LES = longExactSequence m;
            excerptLES(0,2,LES)
    Caveat
        The source and target of m must be homogeneous, and m itself must be degree-preserving.  If the variable
        d has degree different from 1, then the other degrees are adjusted to match d.

        The output module usually has nonzero entries outside of the three columns of interest indicated above.
        This is because we compute M by building a natural exact couple associated to m, and other nonzero
        entries appear organically.
    SeeAlso
        excerptLES
        "Conventions and first examples"
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
    Caveat
        This function gives an error if any of the variables already have values.
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
            declareGenerators(ZZ[x], {a => 1, b => 2, c => 3})
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
            -- build two R-modules
            M = coker map(R^{-1,-1,-1,-1,-1,0,0,0,0},R^{-1,-1,-1,-2,-2,-2,-2,-2},
                {{0,0,0,d,0,0,0,0},{0,0,0,0,d,0,0,0},{0,0,0,0,0,d,0,0},
                 {0,0,0,0,0,0,d,0},{0,0,0,0,0,0,0,d},{-4*d,4*d,4*d,0,0,0,0,0},
                 {5*d,0,0,0,0,0,0,0},{0,5*d,0,0,0,0,0,0},{0,0,5*d,0,0,0,0,0}});
            N = coker map(R^{-1,-1,0,0,0,1,1},R^{0,-1,-1,-2,-2},
                {{0,0,0,d,0},{0,0,0,0,d},{0,0,-4*d,0,0},{0,d,0,0,0},
                 {0,0,3*d,0,0},{0,0,0,0,0},{d,0,0,0,0}});
            -- and an interesting map N <-- M
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

            We use neither convention.  Instead, we place A' exactly halfway between
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
        expectSequenceRing
        (expectSequenceRing,Ring)
    Headline
        accepts rings of the form R[t]
    Usage
        expectSequenceRing Q
    Inputs
        Q:Ring
    Consequences
        Item
            causes an error if Q is not of the form R[t]
    Description
        Text
            Specifically, Q should have exactly one generator over
            its coefficient ring.
        Example
            Q = QQ[t];
            expectSequenceRing Q
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
        Text
            There is also a more-general method available where the argument submods is
            replaced by a module for a ring of the form R[f].  If f acts by inclusions, we
            recover the case above, but if not, then we still obtain an exact couple from
            relative homology.
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

doc ///
    Key
        "Exact couples for Tor and Ext"
        (TorCouple, Module, Module)
        (TorCouple, Symbol, Module, Module)
        (contravariantExtCouple, Module, Module)
        (contravariantExtCouple, Symbol, Module, Module)
        (covariantExtCouple, Module, Module)
        (covariantExtCouple, Symbol, Module, Module)
    Headline
        building couples by applying Tor or Ext to a filtered module or a graded R[t]-module
    Usage
        TorCouple(W, seqmod)
        TorCouple(eSymbol, W, seqmod)
    Inputs
        W:Module
            for some ring R
        seqmod:Module
            for a ring of the form R[t]
        eSymbol:Symbol
            naming the variable usually called "e"; if omitted, this argument is replaced by
            getSymbol "e"
    Outputs
        C:Module
            for the couple ring R[e_1,t_1], an exact couple
    Description
        Text
            The
            filtered versions are described in detail at their individual
            pages @ TO TorCouple @, @ TO contravariantExtCouple @, and
            @ TO covariantExtCouple @.  We invite the reader to check these first,
            since the corresponding spectral sequences are easier to describe.

            The remainder of this page describes an alternative, and more-general way of using these
            three functions.  Specifically, we
            explain how to use a graded R[t]-module in place of a filtered module.

            The idea is to replace a sequence of inclusions with a general sequence of maps.
            The first page then consists of the homology of mapping cones, which play the role
            of the associated graded.
            If the variable t acts by inclusions, the mapping cone is quasi-isomorphic to the
            associated graded, and we recover the filtered case.

            {\bf How to use the R[t] versions}

            The schematic for using these generalized versions: replace the "submods" argument, which is
            an increasing list of submodules of a fixed R-module, with a "seqmod" argument which
            is a module for a ring of the form R[t].

            Such a seqmod can be built directly, or from a sequence of maps using
            @ TO sequenceModule @.

            In the case of contravariantExtCouple, the appropriate sequence module (the one generalizing
                the filtered case) is not the one corresponding to the filtration, but rather the
            cofiltration $X/A_i$.  So, in recovering the filtered case, the variable t acts by surjections.

            These generalized versions also allow a symbol as the first argument, which will be used
            to name the variable "e" in the resulting couple ring.  (The variable "f" in the couple
            ring will be named after the variable (ring seqmod)_0).
        Example
            R = QQ[x,y,z]
            A = coker matrix {{x^2+y^2+z^2}};
            B = coker vars R;
            toSeqMod = cm -> ( -- This short function converts chain modules to sequence modules
                amb := ambient ring cm; -- by forgetting that the action of d squares to zero.
                pres := lift(presentation cm, amb); -- It can be replaced by pushForward once that
                coker(pres | (id_(target pres) ** (matrix {{(amb_0)^2}}))) -- function is fixed.
                );
            B' = toSeqMod chainModule(res B)
            torCouple = prune TorCouple(A,B')
            plotPages((-3..3,-4..2,1..3), prune @@ evaluateInDegree, torCouple)
            covExtCouple = prune covariantExtCouple(A,B')
            plotPages((-3..3,-4..2,1..3), prune @@ evaluateInDegree, covExtCouple)
            contraExtCouple = prune contravariantExtCouple(B' ** (ring B')^{{-4,0}},A) -- see Caveat
            plotPages((-3..1,-5..1,1..3), prune @@ evaluateInDegree, contraExtCouple)
        Text
            It bears mentioning that there are many other spectral sequences arising from Tor and Ext.
            For example, both coordinates could vary instead of just one, or these functors could be composed
            in various
            ways.  For such applications, work directly with the function @ TO exactCouple @.  If you
            obtain anything compelling, the author of this package would appreciate hearing
            about it!
    Caveat
        In the case of @ TO contravariantExtCouple @, the module seqmod should be concentrated in degrees
        that are positive multiples of $deg(t)$.
    SeeAlso
        TorCouple
///

doc ///
    Key
        "Functoriality for Tor and Ext couples"
    Headline
        induced maps between couples and spectral sequences
    Description
        Text
            The exact couples described in @ TO "Exact couples for Tor and Ext" @ are functorial in
            a sense that we explain here.

            As a basic first example, suppose we have two R-modules, $X$ and $Y$, two submodules
            $A \subseteq X$ and $B \subseteq Y$, and a map $g : X \to Y$ with $g(A) \subseteq B$.

            Fix some other R-module W.  Applying Hom(W,-), we obtain two long exact sequences in Ext

            $0 \to Hom(W,A) \to Hom(W,X) \to Hom(W,X/A) \to Ext^1(W,A) \to \cdots $

            $0 \to Hom(W,B) \to Hom(W,Y) \to Hom(W,Y/B) \to Ext^1(W,B) \to \cdots $

            and downwards maps induced by $g$ that give a commuting ladder.

            The pages @ TO "Exact couples for Tor and Ext" @ and @ TO covariantExtCouple @ explain how
            to build each row of this ladder individually, but how can we obtain the downward maps?

            In brief, we encode the action of g by introducing a ring R[g] and reformulate the ladder
            as a long exact sequence over this larger ring.  Each vertical rung becomes a snippet of a graded
            R[g]-module (since there are two rows, we only care about g-degrees 0 and 1).  So the final
            output is a long exact sequence of R[g]-modules recording the two sequences in degrees 0 and
            1 and the map as multiplication by g.
            
            Shapiro's lemma allows us to compute this new sequence in terms of a single covariant Ext couple 
            over R[g].  Specifically, we replace W with an appropriate R[g]-module W', and
            replace X, Y, A, B with a single R[g][t]-module M encoding their comparison maps, and then
            call covariantExtCouple(W',M).

            {\bf A small example}

            Set R = $\QQ[z]$, and build the commuting square
        CannedExample
            |   A - - t - -> X
            |   |            |
            |   g            g
            |   |            |
            |   v            v
            |   B - - t - -> Y
        Text
            by giving a presentation over the ring R[g][t], placing A in bidegree \{0,0}.  The example square
            we have in mind:
        CannedExample
            |  cokernel {3} | z13 | - z^2 -> cokernel {1} | z15 |
            |           |                             |
            |           z                             z
            |           |                             |
            |           v                             v
            |  cokernel {2} | z6 |  - x^2 -> cokernel {0} | z8 |
        Text
            Building it is a matter of naming generators, specifying degrees, and imposing relations.
        Example
            R = QQ[z]; S = R[g][t]; declareGenerators(S,{a=>{0,0,3},x=>{1,0,1},b=>{0,1,2},y=>{1,1,0}});
            M = cospan(z^13*a,z^15*x,z^6*b,z^8*y,g*a-z*b,g*x-z*y,t*a-z^2*x,t*b-z^2*y); isHomogeneous M
            eid = prune @@ evaluateInDegree; (dt, dg) = degree \ (S_0, S_1);
            {A,X,B,Y} = (deg -> prune eid({deg#1},eid({deg#0},M))) \ ({0,0},dt,dg,dt+dg);
            netList {{A, X}, {B, Y}}
        Text
            We also need another module W.
        Example
            W = R^1 / (R_0^10);
        Text
            Since R[g] is projective as an R-module, Shapiro's lemma applies, and Ext can be computed
            over R[g] if we replace W with its extension.
        Example
            W' = extensionInDegree({0}, coefficientRing S, W)
        Text
            And now we can build the couple.
        Example
            couple = prune covariantExtCouple(W',M)
            expectExactCouple couple
        Text
            To view the long exact sequence, we use @ TO excerptCouple @.
        Example
            excerptCouple({-2,2},4,couple)
        Text
            Note that each entry in this long exact sequence is still an R[g]-module.
            This is by design!  The long exact ladder is now encoded as an exact sequence
            of R[g]-modules whose degree 0 part and degree 1 part make up the two exact
            rows,
            and where the action of g gives the ladder's rungs.

            {\bf Interpreting the output}

            As an example, we focus on the bottom nonzero entry in the middle column:
        Example
            relHom = eid({0,2},couple)
        Text
            The degree \{0,2} comes from the conventions explained at @ TO covariantExtCouple @, which
            also explains that relHom, when expanded along
            its g action, should look like this:
        CannedExample
            |   Hom(W,X/A)
            |       |
            |       g
            |       |
            |       v
            |   Hom(W,Y/B)
        Text
            By evaluating in degrees 0 and 1, we verify that relHom has these expected values.
        Example
            (eid({0},relHom), Hom(W,coker map(X,A,{{z^2}})))
            (eid({1},relHom), Hom(W,coker map(Y,B,{{z^2}})))
        Text
            We can also compute the action of g on relHom in these degrees:
        Example
            eid({1}, relHom ** matrix {{g}})
        Text
            {\bf Contravariant functoriality}

            We give a similar example for @ TO contravariantExtCouple @, evaluating on a cofiltration
            (surjections instead of inclusions).
            See @ TO "Exact couples for Tor and Ext" @ for more details on this way of
            using contravariantExtCouple.

            Start with a commuting square of R-modules
        CannedExample
            |   X - - t - ->> C
            |   |             |
            |   g             g
            |   |             |
            |   v             v
            |   Y - - t - ->> D
        Text
            For convenience, set A = ker(t : X - -> C) and B = ker(t : Y - -> D).  If Z is some other
            module, we have a ladder
        CannedExample
            |   0 - - -> Hom(D,Z) - - -> Hom(Y,Z) - - -> Hom(B,Z) - - -> Ext^1(D,Z) - - -> ...
            |   |            |               |               |                 |
            |   g            g               g               g                 g
            |   |            |               |               |                 |
            |   v            v               v               v                 v
            |   0 - - -> Hom(C,Z) - - -> Hom(X,Z) - - -> Hom(A,Z) - - -> Ext^1(C,Z) - - -> ...
        Text
            To use Shapiro's lemma in this setting, we would like to replace the R-module Z with an
            R[g]-module Z' that captures the action of g on these Ext groups.  Such a Z' exists, but
            it is usually not finitely generated.  In this case, Z' would have a copy of Z in
            every non-positive degree.

            In order to maintain finite generation, we replace $R[g]$ with a quotient ring $R[g]/g^n$
            for some $n$ large enough to accommodate the maps we care about.  (In our case, we have
            a two-row ladder, so $R[g]/g^2$ suffices.  If we had other action maps, a similar truncation
            trick would work for larger rings $R[g,h,...]$.)

            Using the ring $R[g]/g^n$, we may take
        CannedExample
            Z' = extensionInDegree({-n+1},Z)
        Text
            since this module is isomorphic to the right Kan extension of Z to the ring $R[g]/g^n$; this
            module is also known as a "coinduced module".

            {\bf A small example}

            Set $R=\QQ[z]$, and build the commuting square
        CannedExample
            |  cokernel {5} | z6 |  - - ->> cokernel {5} | z3 |
            |           |                             |
            |           z5                            z5
            |           |                             |
            |           v                             v
            |  cokernel {0} | z10 | - - ->> cokernel {0} | z7 |
        Text
            As indicated above, since we only care about the map induced by the matrix
        CannedExample
            {0} | g |
        Text
            and not the action of g in other degrees, we may take n=2 since this number exceeds
            the row- and column-degrees of this matrix.
        Example
            erase(symbol x); erase(symbol y);
            n = 2;
            R = QQ[z]; S = (R[g]/g^n)[t]; declareGenerators(S,{x=>{0,0,5},y=>{0,1,0}});
            M = cospan(z^6*x,z^3*t*x,z^10*y,z^7*t*y,g*x-z^5*y,t^2*x,t^2*y); isHomogeneous M
            eid = prune @@ evaluateInDegree; (dt, dg) = degree \ (S_0, S_1);
            {X,C,Y,D} = (deg -> prune eid({deg#1},eid({deg#0},M))) \ ({0,0},dt,dg,dt+dg);
            netList {{X, C}, {Y, D}}
        Text
            Set Z to be $\QQ[z]/z^7$, and build the coinduced/Hom/right-Kan-extension module Z'
        Example
            Z = R^1 / (R_0^7);
            Z' = extensionInDegree({-n+1}, coefficientRing S, Z)
        Text
            We are now ready to build the couple using the larger ring.
        Example
            couple = prune contravariantExtCouple(M,Z')
            excerptCouple({-2,0},4,couple)
        Text
            It may be alarming to see non-zero entries in positions that appear at the point in
            the sequence that should hold Ext^3.  Fortunately, these entries disappear after
            evaluation at the two degrees \{0} and \{-1}, which are the two long exact sequences
            we are trying to relate.  (The degrees are negated because of
            contravariance.  Alternatively, we could have set the degree of g to be -1)

            To check this, we restack and evaluate
            to see both rows.
        Example
            --restack by hand
            Q = (R[j_1,k_1,Degrees=>{{1,-1},{0,2}}])[g]/g^n
            phi = map(Q,ring couple,{j_1,k_1,g,z},DegreeMap=>deg->deg_{2,0,1,3})
            C = phi ** couple
            C0 = evaluateInDegree({0},C)
            expectExactCouple C0
            C1 = evaluateInDegree({-1},C)
            expectExactCouple C1
            excerptCouple({-2,0},4,C0)
            excerptCouple({-2,0},4,C1)
        Text
            We can verify by recomputing these two sequences directly.
        Example
            A = image map(X,,{{z^3}});
            B = image map(Y,,{{z^7}});
            contravariantExtLES(4,X,A,Z)
            contravariantExtLES(4,Y,B,Z)
    SeeAlso
        "Exact couples for Tor and Ext"
///

doc ///
    Key
        "Encoding diagrams as modules"
    Headline
        building graded modules with specified modules in certain degrees, and with specified action maps
    Description
        Text
            Many algorithms in computer algebra accept as input a finite commuting diagram of modules.  
            This can pose a challenge to programmers, since specifying such a diagram can be unwieldy.
            For example, to input a commuting cube, a user could specify twelve maps in a list... but...
            these maps come without any preferred ordering, which makes any convention hard to remember,
            and moreover, the programmer will be forced to include
            many compatibility checks in order to supply useful error messages.
            
            This package takes the following approach, which we illustrate in the case of a commuting
            square (a cube would be similar).
            
            Let R be the base ring, and build a new ring
        CannedExample
            S = R[f,g,Degrees=>{{1,0},{0,1}}]
        Text
            A graded S-module is an infinite grid of R-modules connected by maps induced by multiplication
            by the variables f and g.  In particular, it encodes an infinite number of commuting squares!
            To specify a single commuting square, we restrict attention to the four bidegrees \{0,0\}, \{1,0\}, 
            \{0,1\}, and \{1,1\}.  (We could ask that the module vanish away from these degrees, but in
            practice it is more efficient to just say "we don't care what happens in other degrees".)
            
            {\bf Internal and external degrees}
            
            Some terminology, since the ring R may itself have some grading.  We call this grading 
            "internal" since it happens inside the coefficients.  The
            variables f and g have internal degree zero, even though their external degrees are \{1,0} and
            \{0,1} respectively.  When building a ring, the Degrees option specifies external degrees.  Suppose
            that R has degree length three so that deg(1_R) = \{0,0,0}, for example.  We have then
        CannedExample
            deg(f) = {1,0,0,0,0}
            deg(g) = {0,1,0,0,0}
        Text
            The first two coordinates are the external degree, and the last three are internal.  To obtain
            this information about a ring, you can use @ TO internalDegreeIndices @ and 
            @ TO externalDegreeIndices @.  For example:
        Example
            R = QQ[x,y,Degrees=>{{1,2,3},{4,5,6}}]
            S = R[f,g,Degrees=>{{1,0},{0,1}}]
            internal = internalDegreeIndices S
            external = externalDegreeIndices S
        Text
            Later, given a multidegree, it is easy to find the internal and external degrees
        Example
            deg = {2,3,4,5,6}
            deg_internal
            deg_external
        Text
            We generally wish to accommodate a wide range of coefficient rings R, which in particular means
            we accommodate any number of internal degrees, including towers of rings where the coefficients
            themselves have coefficients, etc.  In such cases, all degrees that are not external count as
            internal.
            
            {\bf Example: encoding a commuting square}
            
            In this example, we take R=QQ[z].
        Example
            R = QQ[z]
            S = R[f,g,Degrees=>{{1,0},{0,1}}]
        Text
            We wish to encode a commuting square with the
            general layout
        CannedExample
            |   A - -g- -> B
            |   |          |
            |   f          f
            |   |          |
            |   v          v
            |   C - -g- -> D
        Text
            The downward maps will be encoded by the action of f, and the rightward maps by g.  Here
            is a particular example.
        CannedExample
            |  cokernel {3} | z13 | - z^2 -> cokernel {1} | z15 |
            |           |                             |
            |           z                             z
            |           |                             |
            |           v                             v
            |  cokernel {2} | z6 |  - x^2 -> cokernel {0} | z8 |
        Text
            We aim to encode this commuting square as an S-module.  
            The external degree will give
            the position in the commuting square, and the internal degree will record the R-degree.  Then,
            in a multidegree \{r,c,d\}, r gives the row (0 or 1 in our case), c gives the column (also
                0 or 1), and d gives the internal degree (usual grading on R=QQ[z]).
            
            {\bf Declaring generators}
            
            We begin the process in the upper left corner, with the module
        CannedExample
            |  A = cokernel {3} | z13 |
        Text
            Let us name the generator of this module $a$.  Since we
            are in the upper left corner, the external degree is \{0,0\}.  And since the generator appears
            in R-degree 3, the internal degree is \{3\}.  In total, then, the degree of $a$ is \{0,0,3\}.
            
            Similarly, let $b$, $c$, and $d$ be generators with external degrees \{0,1}, \{1,0}, and \{1,1},
            and with interal degrees 1, 2, and 0.  This information can be given to M2 using the function
            @ TO declareGenerators @.
        Example
            declareGenerators(S, {a => {0,0,3}, b => {0,1,1}, c => {1,0,2}, d => {1,1,0}})
        Text
            We must now impose relations on these four generators so that the four modules match our intent,
            and same for the maps.
            
            The first four relations come from the original descriptions of A, B, C, and D:
        CannedExample
            z^13*a
            z^15*b
             z^6*c
             z^8*d
        Text
            The next four relations come from the descriptions of the four maps:
        CannedExample
            g*a - z^2*b
            g*c - z^2*d
            f*a  -  z*c
            f*b  -  z*d
        Text
            The first of these, for example, forces ga = z^2b, and this is what we want since g is 
            supposed to act by the horizontal map, which sends the generator for A to z^2 times the
            generator for B.   With these four relations, the action of f and g is determined in 
            the four degrees of interest.
        Example
            M = cospan(z^13*a, z^15*b, z^6*c, z^8*d,
                       g*a - z^2*b, g*c - z^2*d, f*a - z*c, f*b - z*d)
        Text
            The module M now contains a complete description of the commuting square.
            
            {\bf Evaluating a module at various external degrees}
            
            In order to check that M is correct, we can use the function @ TO evaluateInDegree @ to
            make sure the proper R-module lives in each of the four external degrees.
        Example
            netList apply(2, r -> apply(2, c -> prune evaluateInDegree({r,c}, M)))
        Text
            
            {\bf Evaluating a module at a structure map}
            
            In order to check the action of f and g, we use another form of evaluateInDegree.
        Example
            prune structureMap({0,0},,g,M)
            prune structureMap({1,0},,g,M)
            prune structureMap({0,0},,f,M)
            prune structureMap({0,1},,f,M)
        Text
            
            {\bf Example calculation: computing kernels of cokernels}
            
            In order to perform calculations on a diagram encoded as above, one main strategy involves
            changing which variables are internal and which are external.  In this example, we take
            the cokernel of the downward maps, and then take the kernel of the induced rightward map,
            resulting in a single R-module 
        CannedExample
            ker( coker(A - -g- -> B) - - -> coker(C - -g- -> D) )
        Text
            We want to take the cokernel of the g action map in a way that retains the action of f.  So
            build a ring where f is an internal variable, and only g is external:
        Example
            S' = R[f][g]
            phi = map(S',S,DegreeMap=>deg->deg_{1,0,2})
            isHomogeneous phi
            M' = phi ** M
        Text
            Since only g is external, we may evaluate to obtain a map of R[f]-modules, and then take its
            cokernel:
        Example
            cokerg = coker structureMap({0},,g,M')
        Text
            Now g is gone, and f is an external variable.  We may evaluate to obtain the map on cokernels,
            and take the kernel:
        Example
            ker structureMap({0},,f,cokerg)
    SeeAlso
        evaluateInDegree
        cospan
        declareGenerators
///

doc ///
    Key
        "Elemetary introduction: solving linear equations in abelian groups"
    Headline
        a discussion of the equation 3x+6y=0.
    Description
        Text
            Suppose we are interested in solving the equation $3a_1 + 6a_2 = 0$ for pairs of elements 
            $a_1, a_2$ drawn from an abelian group $A$.  Introduce notation for the solution set
            
            $GA := \{ (a_1, a_2) \in A^2 | 3a_1 + 6a_2 = 0 \}$
            
            We make three observations about G:
            
            1) There is a natural abelian group structure on $GA$
            
            2) If $f: A \to B$ is a map, then any solution $(a_1, a_2)$
            for $A$ pushes along $f$ to a solution $(f(a_1), f(a_2))$ for $B$.  We call the resulting
            induced map $Gf : GA \to GB$.  
            In fact, $Gf$ is a group homomorphism.
            
            3) If $i : A \subseteq X$ is an inclusion of abelian groups, then any solution in $A$ is
            also a solution in $X$; in other words, the map $Gi : GA \subseteq GX$ is an inclusion. 
            In fact, the maps from 2) give an exact sequence $0 \to GA \to GX \to G(X/A)$ so that 
            $G$ is left exact.
            
            Observation 3) leads to a divide-and-conquor strategy for understanding $GX$.  If $X$ is 
            a complicated abelian group,
            and if $A \subseteq X$ is an easier subgroup, then the exact sequence from 3) says that every
            solution in $X$ is patched from solutions in $GA$ and $G(X/A)$.  One must bear in mind that 
            the map $GX \to G(X/A)$ may 
            not be surjective, so not all solutions in $X/A$ lift to solutions in $X$.  However, if we are able to
            find generators for $I = image(G(X) \to G(X/A))$, then we can lift each of these to a solution in
            $X$.  From these, every solution can be obtained by adding elements of $GA$.
            
            To assist in our understanding of $I$, we enlist the extended sequence (called (*) below)
            
            $0 \to GA \to GX \to G(X/A) \to (R^1G)A \to (R^1G)X \to (R^1G)(X/A) \to 0$    {\bf (*) }
            
            which we will construct shortly.  (Here, $R^1G$ is an operation like $G$ in that it satisfies
                properties 1) and 2) above.  It is known as the "first right derived functor of G".  The
                zeroth right derived functor of $G$ is $G$ itself: $R^0G = G$)
            
            Then, $I = kernel(\delta : G(X/A) \to (R^1G)A)$, the kernel of the 
            connecting homomorphism.  
            --One theoretical benefit of the connecting homomorphism is that its
            --source and target depend only on $X/A$ and $A$.  (Of course, the values of $\delta$
            --depend in an essential way on the way $X$ is built from these two groups!)
        
            {\bf Description of the derived functor $R^1G$ and connecting homomorphism}
            
            By a short homological argument, we shall prove an isomorphism
            
            $(R^1G)(A) = A / \{ 3a_1 + 6a_2 | a_1, a_2 \in A \}$ for any abelian group $A$.
            
            Proof:
            
            Let $W$ be the abelian group $\ZZ^2 / (3,6)$.  
            The functor $G$ is naturally isomorphic to $Hom(W,-)$, so $R^1G$ is naturally 
            isomorphic to $Ext^1(W,-)$, which can be computed from a free resolution $0 \to \ZZ \to \ZZ^2 \to W$.
            This resolution also shows that $R^pG = 0$ for $p > 1$.
            
            We now describe the connecting homomorphism.
            
            Let $(x_1,x_2) \in X^2$ be a solution in $X/A$ so that $3x_1 + 6x_2 \in A$.  Then, the connecting 
            homomorphism is given by the formula
            
            $\delta(x_1,x_2) = 3x_1 + 6x_2 \in (R^1G)A$.
            
            We must check that $\delta$ does not depend on the choice of representatives.  Suppose $(x_1', x_2')$
            represents the same solution in $X/A$.  Then $x_1' - x_1 \in A$, $x_2' - x_2 \in A$, and so
            
            $\delta(x_1',x_2') - \delta(x_1,x_2) = 3(x_1'-x_1) + 6(x_2'-x_2)$
            
            is of the form $3a_1 + 6a_2$, thereby vanishing in $(R^1G)A$.
            
            {\bf Multi-step filtrations}
            
            If $X$ is quite complicated, then subgroup $A$ may not be simple enough
            to analyze directly.  It is natural to iterate the procedure above, choosing 
            subgroups $A' \subseteq A$,
            and $A'' \subseteq X/A$ to subdivide the filtration dyadically as needed.  
            It is equivalent, and better-
            organized, to allow ourselves a filtration by a chain of subgroups
            
            $A_0 \subseteq A_1 \subseteq \cdots \subseteq A_n = X$.
            
            We want to choose enough subgroups so that the quotients $A_q / A_{q-1}$ are not too
            difficult.  Specifically, we assume that each $A_q / A_{q-1}$ is easy enough that we can obtain 
            good information about the $(R^pG)(A_q / A_{q-1})$.
            
            Of course, we would be happier to understand the quotients $G(A_q) / G(A_{a-1})$, since from these
            the procedure for computing $GX$ is quite clear: simply pick representatives for each group 
            $G(A_q) / G(A_{q-1})$; then, each element of $GX$ expands uniquely as a sum of repesentatives, one 
            drawn from each successive quotient.
            
            In other words, we start with information about $(R^pG)(A_q/A_{q-1})$ with the goal of understanding
            each $G(A_q) / G(A_{q-1})$.
            
            {\bf The spectral sequence}
            
            --Sometimes, we can obtain enough control on the $(R^pG)(A_q/A_{q-1})$ to tell us what we
            --want to know about $G(A_q) / G(A_{q-1})$.  
            The proper tool for this is a spectral sequence, which is a sequence
            of two-dimensional arrays ("pages") of abelian groups and maps in a certain pattern, and with certain 
            compatibilities.  Spectral sequences can be notationally dense, but our example will remain 
            light.
            
            The first page of our spectral sequence looks like this:
            
        CannedExample
            page 1, with differential of degree {1, -1}:
            +----++----+-----------------+-----------------+---+
            |q=4 ||0   |       ...       |       ...       |0  |
            +----++----+-----------------+-----------------+---+
            |q=3 ||0   |    G A_3/A_2    |  R^1G A_3/A_2   |0  |
            +----++----+-----------------+-----------------+---+
            |q=2 ||0   |    G A_2/A_1    |  R^1G A_2/A_1   |0  |
            +----++----+-----------------+-----------------+---+
            |q=1 ||0   |    G A_1/A_0    |  R^1G A_1/A_0   |0  |
            +----++----+-----------------+-----------------+---+
            |q=0 ||0   |      G A_0      |    R^1G A_0     |0  |
            +----++----+-----------------+-----------------+---+
            |q=-1||0   |0                |0                |0  |
            +----++----+-----------------+-----------------+---+
            |    ||p=-1|p=0              |p=1              |p=2|
            +----++----+-----------------+-----------------+---+
        Text
            This is a table of (R^pG)(A_q / A_{q-1}), essentially.  There are maps (not pictured) of degree \{1,-1}
            relating the groups in the p=0 column to those in the p=1 column.  They are called differentials,
            and their degree changes depending on the page number.  The differentials on page one are given
            by the composite
            
            $D_1 : G(A_{q+1} / A_q) \to (R^1G)(A_q) \to (R^1G)(A_q / A_{q-1})$
            
            where the first map is the connecting homomorphism in the sequence (*) for the pair
            $A_q \subseteq A_{q+1}$, and the second map comes from property 2) for $R^1G$.
            
            Note that $D_1$ squares to zero because our diagram has only two nonzero columns.  
            The modules on page two can be computed as the homology of the differential $D_1$.
            
            We think of the p=0 column as containing potential
            solutions, and the p=1 column as containing potential obstructions to solutions.  The differential's
            job is to cancel out potential solutions with potential obstructions as appropriate, so that 
            these elements do not appear on the next page.
            
            After turning the page enough times, (for an $n$-step filtration, $n$ turns is enough), the answer 
            stabilizes to the "infinity" page, which looks like this
        CannedExample
            page r >> 0, with differential of degree {1, -r}:
            +----++----+-----------------+------------------------+---+
            |q=4 ||0   |       ...       |       ...              |0  |
            +----++----+-----------------+------------------------+---+
            |q=3 ||0   |  G A_3 / G A_2  |  R^1G A_3 / R^1G A_2   |0  |
            +----++----+-----------------+------------------------+---+
            |q=2 ||0   |  G A_2 / G A_1  |  R^1G A_2 / R^1G A_1   |0  |
            +----++----+-----------------+------------------------+---+
            |q=1 ||0   |  G A_1 / G A_0  |  R^1G A_1 / R^1G A_0   |0  |
            +----++----+-----------------+------------------------+---+
            |q=0 ||0   |      G A_0      |        R^1G A_0        |0  |
            +----++----+-----------------+------------------------+---+
            |q=-1||0   |0                |0                       |0  |
            +----++----+-----------------+------------------------+---+
            |    ||p=-1|p=0              |p=1                     |p=2|
            +----++----+-----------------+------------------------+---+
        Text
            Now the p=0 column consists of actual solutions, and the p=1 column consists of hibernating 
            obstructions that might
            reawaken if $X$ maps somewhere else.
            
            To see how this might be useful, we mention some general facts about spectral sequences.  
            First, the differential on 
            each page always squares to zero.  Second, the next page is isomorphic to to the homology of
            the differential on the previous page.  Consequently, any zero entry persists 
            forever, and more generally, entries only get smaller (in the sense of subquotients) as 
            the pages turn.
            
            In our example, the r=1 page is concentrated in the columns p=0,1, so every subsequent page
            is similarly bounded; and the 
            differential on page $r$ has degree \{1,-r}, so degree considerations ensure the square-zero
            property.  
            Taking homology replaces the p=0 groups with kernels and the p=1 groups with
            cokernels.  This operation is therefore responsible for removing the extraneous solutions that 
            may be present on the first page.  Let's try it!
        
            {\bf Running the spectral sequence on an example}
            
            Let $X$ be the abelian group with three generators $a,b,c$ and defining relations $4a=b$,
            $4b=c$, and $4c=a$.  Let $A_0 = span(a+b+c)$, $A_1 = span(a+b+c,9a)$, and $A_2=X$.  Using the natural
            isomorphism of $G$ and $Hom(W,-)$ mentioned earlier, we build the spectral sequence as follows:
        Example
            declareGenerators(ZZ, {a=>{}, b=>{}, c=>{}}); -- degrees in X are empty
            X = cospan(4*a-b,4*b-c,4*c-a)
            A0 = image map(X,,matrix {a+b+c})
            A1 = image map(X,,matrix {a+b+c,9*a})
            W = coker matrix {{3},{6}}
            M = prune covariantExtCouple(W,{A0,A1,X});
            plotPages((-1..2,-1..4,1..3), prune @@ evaluateInDegree, M)
        Text
            Some things we can conclude just from the first page:
            
            1) The differential on page 1 vanishes.  The only one with any potential maps from $\ZZ/7$ to
               $\ZZ/3$, and is therefore zero.
            
            2) On later pages, there is only one potential nonzero differential, and it will happen 
               (if it happens) on page 2.
               Moreover, since $R^1GA_0$ is a simple group, this differential is either zero or surjective.
               In either case, we have enough information to compute the order of $GX$: it will 
               have size either $3^4 \cdot 7$ or $3^3 \cdot 7$.
            
            3) For degree reasons, the entries in positions (0,0) and (0,1) survive to the infinity page.
               Practically speaking, this means that every potential solution in $G(A_1/A_0)$ lifts to a true
               solution.
            
            4) We also know that $R^1GX$ is either $\ZZ/3$ or $\ZZ/3 \oplus \ZZ/3$, and this information would
               also tell us the order of $GX$.  ($\ZZ/9$ is not possible because $R^1G$ always has exponent 3)
            
            5) We know that $GX$ is not cyclic because its subgroup $G(A_0)$ is not cyclic.
            
            Let's try to understand $GX$ by hand using this information.
            
            We start by finding the nine solutions in $A_0$.  This is easy, since $4(a+b+c)=a+b+c$ and so
            $3(a+b+c)=0$.  Therefore, $A_0 \cong \ZZ/3$ and so any pair in $A_0$ solves the equation.  Explicitly,
            setting $w=a+b+c$, we have nine solutions $(k*w,l*w) \in GA_0$ for $k,l \in \{0,1,2\}$.
            
            Next, by observation 3), we know that solutions in $G(A_1/A_0)$ will lift and that 
            $G(A_1)/G(A_0) = \ZZ / 7$.  Since $A_1$ is 
            generated by $A_0$ together with one new element, we know $A_1/A_0$ is cyclic with generator $9a$.
            Its order must be prime to 3 because $R^1G$ vanishes on this group.  We may therefore divide 
            the equation $3x+6y=0$ by three, obtaining $x+2y=0$ and $x=-2y$.  Taking $y=9a$ gives the solution
            $(-18a,9a)$ in $A_1/A_0$.  The group $A_1/A_0$ is nonzero since $G(A_1/A_0) = \ZZ/7$, and thus
            this solution is nonzero, since its second coordinate generates the nonzero cyclic group $A_1/A_0$.
            This proves that we have found a generator for $G(A_1)/G(A_0)$, since any nonzero element generates
            $\ZZ/7$.  
            Observing that the solution $(-18a,9a)$
            also works in $GA_1$, we find that we have fully computed $GA_1$: it is isomorphic
            to $(\ZZ/3)^2 \oplus \ZZ/7$, generated by $(a+b+c,0)$, $(0,a+b+c)$, and $(-18a,9a)$.
            
            We could now contend with $G(X/A_1)$, which is isomorphic to $\ZZ/3 \oplus \ZZ/3$.  We know that
            either all nine or only three of these lift to solutions in GX.
            
            Instead, let's compute $R^1GX$.  Since the p>1 derived functors vanish, $R^1G$ is right exact.
            Moreover, by the general calculation, it sends $\ZZ$ to $\ZZ/3$.  This proves that $R^1G$ 
                coincides with $- \otimes \ZZ/3$.
            By right exactness, we may
            apply $R^1G$ directly to
            the presentation for $X$ to obtain a matrix over $\ZZ/3$:
        Example
            prune coker(map(ZZ/3,ZZ) ** (presentation X))
        Text
            According to this calculation, the total $\ZZ/3$-dimension of the p=1 column at the infinity 
            page must
            be 1.  The only way this could happen is if the differential on page 2 is surjective.
            
            And indeed, that is what we see on pages 2 and 3.
            
            (In fact, the author was aided by knowledge of pages 2 and 3 in the construction of the argument
                above.  Knowing the answer often makes proving it easier.)
            
            Using this package, it is also possible to compute directly the differentials on these pages.
            First we compute the pages as modules
        Example
            E1 = prune pageModule(1, D, M)
            E2 = prune pageModule(2, D, M)
            E3 = prune pageModule(3, D, M)
        Text
            Each of these modules is over a different ring.  Specifically, page $r$ is a module for 
            $\ZZ[D_r]/(D_r)^2$ where the
            degree of $D_r$ is \{1,-r}.
            
            We may compute the interesting differential on page two using the following command.
        Example
            prune structureMap({0,2},{1,0},D_2,E2)
        Text
            So we see directly that this differential is a surjection.
            
            {\bf A peek at the exact couple}
            
            One of the main headaches of a spectral sequence is that the modules and differentials on one
            page compute the next page of modules, but not the differentials.
            
            However, there is another object, called an exact couple, that determines a spectral sequence 
            but also has the property that each of its "pages" fully determines the next.  It may not be
            true that every spectral sequence is constructed by means of an exact couple, but certainly
            most of the common ones can.  
            
            An exact couple is a 
            generalization of the exact sequence (*) only in that it consists in this sequence for every pair 
            $A_{q-1} \subseteq A_q$.  In order to acccommodate so many long exact sequences, and because these
            sequences overlap, we adjust the layout of (*) to a stackable N-shaped zig-zag
        CannedExample
            |
            |      GX        R^1GX
            |      ^ \        ^ \
            |      |  \       |  \
            |      |   \      |   \
            | 0    |  G(X/A)  |   R^1G(X/A)
            |  \   |     \    |     \
            |   \  |      \   |      \
            |    \ |       \  |       \
            |     GA       R^1GA       0
        Text
            resulting in a larger, more complicated-looking diagram:
        CannedExample
            |      .        .
            |      :        :
            |     GA2     R1GA2     
            |      ^\       ^\
            |      | \      | \
            |      |  \     |  \
            | 0    | GA2/A1 | R1GA2/A1
            |  \   |    \   |    \
            |   \  |     \  |     \
            |    \ |      \ |      \
            |     GA1     R1GA1     0
            |      ^\       ^\
            |      | \      | \
            |      |  \     |  \
            | 0    | GA1/A0 | R1GA1/A0
            |  \   |    \   |    \
            |   \  |     \  |     \
            |    \ |      \ |      \
            |     GA0     R1GA0     0
            |      ^\       ^\
            |      | \      | \
            |      |  \     |  \
            | 0    | *GA0*  |  R1GA0
            |  \   |    \   |    \
            |   \  |     \  |     \
            |    \ |      \ |      \
            |     0        0        0
        Text
            Let's put coordinates on this diagram.  The entry marked *GA0* is in position (0,0).
            The diagonal maps, which go downward and rightward, have degree (1,-1), and the upward maps
            have degree (0,2).  So the diagram is supported in degrees (x,y) where x+y is even.
            --If x and y are even, say, (x,y) = (2p,2q), we have the group $(R^pG)(A_q/A_{q-1})$.  And if they
            --are odd, say (x,y) = (2p-1,2q+1), then we have $R^pA_q$.
            
            If we introduce a ring R = \ZZ[e,f, Degrees=>\{\{1,-1},\{0,2}}], then
            the diagram becomes a graded R-module where e acts by the diagonal maps, 
            and f acts by the upwards maps.  
            
            In the example above, we named this module M.  Observe that M encodes an exact couple
            by setting $E$ to be the odd part and $A$ to be the even part.
            The purpose of this package is to make easy the calculation of couples like $M$, their derived
            couples, and the resulting spectral sequences.
    SeeAlso
        "Encoding diagrams as modules"
        "Conventions and first examples"
        "Bockstein spectral sequence"
        "Serre spectral sequence in homology"
        "Exact couples for Tor and Ext"
        "Functoriality for Tor and Ext couples"
///

doc ///
    Key
        "An exact couple associated to a Young tableau"
    Headline
        a standard filling provides a filtration by monomial ideals
    Description
        Text
            Let $\lambda$ be a partition, and let $I \subseteq R=\QQ[x,y]$ be its
            corresponding monomial ideal.  Every standard filling of $\lambda$ provides a 
            filtration by monomial ideals of the interval
            
            $I \subseteq ... \subseteq R$
            
            We compute the resulting homology spectral sequence for every standard tableau
            of shape \{3,2\}, where homology means Tor(coker vars R, -).
        Example
            needsPackage "PushForward"
            needsPackage "SpechtModule"
            lambda = {3,2}
            R = QQ[x,y,Degrees=>{{1,0},{0,1}}]
            I = ideal apply(lambda | {0}, 0..#lambda, (k,l)->x^l*y^k)
            boxes = flatten apply(#lambda,r->apply(lambda#r,c->(r,c)))
            p = new Partition from lambda
            tableaux = toListOfTableaux(standardTableaux p)
            filts = apply(tableaux, t -> reverse apply(1 + sum lambda, k -> 
                    module(I + ideal apply(boxes, b -> (if t_b >= k then x^(b#0)*y^(b#1) else 0_R)))))
            disp = (deg,E) -> prune pushFwd(map(R,QQ),evaluateInDegree(deg,E))
            for i in 0..<#tableaux do (
                print(tableaux#i);
                filt = filts#i;
                couple = prune TorCouple(coker vars R, filt);
                plotPages((-1..3,-1..6,1..(sum lambda)), disp, couple);
                )
        Text
            By changing the display function, we can view degree information
        Example
            disp = (deg,E) -> degrees minimalPresentation evaluateInDegree(deg,E);
            for i in 0..<#tableaux do (
                print(tableaux#i);
                filt = filts#i;
                couple = prune TorCouple(coker vars R, filt);
                plotPages((-1..3,-1..6,1..(sum lambda)), disp, couple);
                )
///
