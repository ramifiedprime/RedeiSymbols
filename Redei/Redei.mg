// This will be an implementation of computing Redei symbols.
// This functionality has not yet been implemented in magma
// so is new and likely wrong.  The ideas will be taken from
// Stevenhagens paper on Redeis quadratic artin symbol.

///////////////////////////////////////
// Useful additional functions.
///////////////////////////////////////
//Returns the quadratic discriminant associated to \mathbb{Q}(\sqrt{a})
function QuadraticDiscriminant(a)
	if IsSquare(a) then return 1; end if;
	a:=&*[p: p in PrimeFactors(a)| IsOdd(Valuation(a,p))]; //redefine as squarefree part
	if a mod 4 eq 1 then return a;
	else return 4*a; 
	end if;
end function;

//Returns a list of rational primes which have ramification in the field $K$ which
//is assumed to be Galois.
function RamifiedRationalPrimes(K)
	return [p: p in PrimeFactors(Discriminant(AbsoluteField(K)))];
end function;

// returns the ramification index of the rational prime $p$ in the GALOIS extension K/QQ
function MyRamificationIndex(K,p)
	Factorisation(7*RingOfIntegers(E))[1][2];
end function;

///////////////////////////////////////
// Tests for errors catching.
///////////////////////////////////////
function TestHilbert(a,b) // 0 for all were 1, 1 for something was not trivial
	divs:={p : p in PrimeDivisors(a) cat PrimeDivisors(b)};
	for p in divs do
		if HilbertSymbol(a,b,p) ne 1 then return 1; end if;
	end for;
	return 0;
end function;

function TestDefined(a,b,c) // returns 1 if fails Hilbert testing, 2 if fails gcd testing, 0 is pass
	if 1 in {TestHilbert(a,b), TestHilbert(b,c), TestHilbert(a,c)} then return 1; end if;
	discrims:=[QuadraticDiscriminant(x): x in [a,b,c]];
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

///////////////////////////////////////
// Redei Symbol specific functions.
///////////////////////////////////////

//Constructs E/K_a/QQ and the beta needed for finding some choice of F
function EFieldBetaConstructor(a,b: abs:=false, minram:=false) 
	if IsSquare(a) or IsSquare(b) or IsSquare(a/b) then return 9; end if;
	Ka:=QuadraticField(a);
	isnorm,beta:=NormEquation(RingOfIntegers(Ka), b);
	if not isnorm then return 3; end if;
	beta:=beta[1];
	P<x>:=PolynomialRing(RationalField());
	E:=ext<Ka|x^2-b:Abs:=true>;
	return true, E, beta;
end function;

function MinimallyRamifiedFConstructor(a,b,E,beta) // Uses 7.1 of Stevenhagen a lot!
	P<x>:=PolynomialRing(RationalField());
	K:=NumberField(x^2-a*b);
	E:=AbsoluteField(E);
	P<x>:=PolynomialRing(E);
	F:=AbsoluteField(NumberField(x^2-beta));
	Delta_a:=QuadraticDiscriminant(a);
	Delta_b:=QuadraticDiscriminant(b);

	ramprimes:=[p: p in RamifiedRationalPrimes(F)| MyRamificationIndex(F,p) ne MyRamificationIndex(K,p)];
	oddprimes:=[p:p in ramprimes| IsOdd(p)];
	avoidablyramifiedat2:= 2 in ramprimes and (IsOdd(Delta_a) or IsOdd(Delta_b));//the second statement is for forced ramification
	for p in oddprimes do
		if Delta_a mod p eq 0 and Delta_b mod p eq 0 then continue; // Forced Ramification
		else beta:=p*beta;end if;
	end for;
	F:=AbsoluteField(NumberField(x^2-beta));
	//2-minimal ramification
	if avoidablyramifiedat2 then
		if (IsOdd(Delta_a) and IsOdd(Delta_b))//Lemma 7.1 part 1
		or (IsEven(Delta_a) and Delta_b mod 8 eq 1) //Lemma 7.1 part 3
		or (IsEven(Delta_b) and Delta_a mod 8 eq 1) 
		then
			for t in [-1,2,-2] do
				beta_t:=beta*t;
				if not 2 in [p: p in RamifiedRationalPrimes(AbsoluteField(NumberField(x^2-beta_t)))| MyRamificationIndex(F,p) ne MyRamificationIndex(K,p)] then
					beta:=beta_t;
					F:=AbsoluteField(NumberField(x^2-beta));
				end if;
			end for;
		else
			continue; // need to handle case 4!!
		end if;
	end if;
	return F;
end function;

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