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
