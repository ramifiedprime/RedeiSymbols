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

function ErrorCodesRedei(x)
	if x eq 1 then return "UNDEFINED: Hilbert symbols of inputs failed to be locally trivial everywhere.";
	elif x eq 2 then return "UNDEFINED: The GCD of the associated discriminants of the inputs is not 1.";
	elif x eq 3 then return "UNDEFINED: a,b does not give a solution to x^2=ay^2+bz^2";
	elif x eq 9 then return "UNDEFINED: Something unexpected happened";
	else "Unknown error code";end if;
end function;

function FFieldConstructor(a,b: abs:=false) //Constructs F/E/K_a/QQ in this order. Abs returns absolute of F.
	if IsSquare(a) or IsSquare(b) or IsSquare(a/b) then return 9; end if;
	Ka:=QuadraticField(a);
	isnorm,beta:=NormEquation(RingOfIntegers(Ka), b);
	if not isnorm then return 3; end if;
	beta:=beta[1];
	P<x>:=PolynomialRing(RationalField());
	E:=ext<Ka|x^2-b:Abs:=true>;
	P<x>:=PolynomialRing(E);
	F:=ext<E|x^2-beta:Abs:=true>;
	if abs then F:=AbsoluteField(F); end if;

	return true, F;
end function;

primesForMinimalRamification:=PrimeFactors(Discriminant(AbsoluteField(F)));
Delta_a:=PrimeFactors(Discriminant(Ka));
Delta_b:=PrimeFactors(Discriminant(Kb));

//AbsoluteField(F) will reground to \mathbb{Q}.

// function ClassicalRedei

// function BuildMinimallyRamifiedField(a,b)// returns F, K
// 	K:=QuadraticField(a*b);
// 	P<x>:=PolynomialRing(K);
// 	//To be completed.
// end function;

function Redei(a,b,c)
	T:=TestDefined(a,b,c);
	if T ne 0 then print ErrorCodesRedei(T); return 0; end if;
	if IsSquare(a) or IsSquare(b) or IsSquare(c) then return 1 ; end if;
	print "Valid!";
	return 0;
end function;