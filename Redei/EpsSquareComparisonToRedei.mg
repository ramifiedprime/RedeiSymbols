load "../Redei/Redei.mg";
MAX:=10000;
K:=QuadraticField(-3);
OK:=RingOfIntegers(K);
gal, aut, galIdent:=AutomorphismGroup(K);

interestingprimes:=[p: p in PrimesUpTo(MAX)| p mod 24 eq 1];
for p in interestingprimes do
	_,eps:=NormEquation(K, p);
	eps:=eps[1];
	epsbar:=[galIdent(g)(eps): g in gal| galIdent(g)(eps) ne eps][1];
	K_eps, m_eps:=Completion(K,eps*OK);
	print IsSquare(m_eps(epsbar)) eq (RedeiSymbol(-3,p,p) eq 1);
end for;