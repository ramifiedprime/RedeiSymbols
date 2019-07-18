// This will be an implementation of computing Redei symbols.
// This functionality has not yet been implemented in magma
// so is new and likely wrong.  The ideas will be taken from
// Stevenhagens paper on Redeis quadratic artin symbol.

function TestHilbert(a,b) // 0 for all were 1, 1 for something was not trivial
	divs:={p : p in PrimeDivisors(a) cat PrimeDivisors(b)};
	for p in divs do
		if HilbertSymbol(a,b,p) ne 1 then return 1; end if;
	end for;
	return 0;
end function;

function TestDefined(a,b,c) // returns 1 if fails Hilbert testing, 2 if fails gcd testing, 0 is pass
	if 1 in {TestHilbert(a,b), TestHilbert(b,c), TestHilbert(a,c)} then return 1; end if;
	discrims:=[Discriminant(QuadraticField(x)): x in [a,b,c]];
	if GCD(discrims) ne 1 then return 2; end if;
	return 0;
end function;

function DefinedErrorCodes(x)
	if x eq 1 then return "UNDEFINED: Hilbert symbols of inputs failed to be locally trivial everywhere.";
	elif x eq 2 then return "UNDEFINED: The GCD of the associated discriminants of the inputs is not 1.";
	else "Unknown error code";end if;
end function;

function BuildMinimallyRamifiedField(a,b)// returns F, K
	K:=QuadraticField(a*b);
	P<x>:=PolynomialRing(K);
	//To be completed.
end function;

function Redei(a,b,c)
	T:=TestDefined(a,b,c);
	if T ne 0 then print DefinedErrorCodes(T); return 0; end if;
	if IsSquare(a) or IsSquare(b) or IsSquare(c) then return 1 ; end if;

end function;