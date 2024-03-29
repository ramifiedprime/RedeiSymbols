// This will be an implementation of computing Redei symbols.
// This functionality has not yet been implemented in magma
// so is new and likely wrong.  The ideas will be taken from
// Stevenhagens paper on Redeis quadratic artin symbol.

///////////////////////////////////////
// Useful additional functions.
///////////////////////////////////////
//Returns the squarefree part of the input integer.
function SquarefreeIntegerPart(a)
	if IsSquare(a) then return 1; 
	elif IsSquare(-a) then return -1;end if;
	return Sign(a)*(&*[p: p in PrimeFactors(a)| IsOdd(Valuation(a,p))]); //redefine as squarefree part
end function;

//Returns the quadratic discriminant associated to \mathbb{Q}(\sqrt{a})
function QuadraticDiscriminant(a)
	a:=SquarefreeIntegerPart(a);
	if a mod 4 eq 1 then return a;
	else return 4*a; 
	end if;
end function;

//Returns a list of rational primes which have ramification in the field $K$ which
//is assumed to be Galois.
function RamifiedRationalPrimes(K)
	return PrimeFactors(Discriminant(Integers(AbsoluteField(K))));
end function;

// returns the ramification index of the rational prime $p$ in the GALOIS extension K/QQ
function MyRamificationIndex(K,p)
	return Factorisation(p*RingOfIntegers(AbsoluteField(K)))[1][2];
end function;

///////////////////////////////////////
// Tests for errors catching and general
// admin functions.
///////////////////////////////////////
function TestHilbert(a,b) // 0 for all were 1, else for something was not trivial
	divs:={p : p in PrimeDivisors(a) cat PrimeDivisors(b)}; // dont need to check infinity by reciprocity
	for p in divs do
		if HilbertSymbol(a,b,p) ne 1 then return p; end if;
	end for;
	if a lt 0 and b lt 0 then return -1; end if; //infinite place
	return 0;
end function;

function TestRedeiDefined(a,b,c) // returns 1 if fails Hilbert testing, 2 if fails gcd testing, 0 is pass
	if {TestHilbert(a,b), TestHilbert(b,c), TestHilbert(a,c)} ne {0} then return 1; end if;
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

function RtnValues(triv, additive)
	if triv then
		if additive then return 0;
		else return 1; end if;
	else
		if additive then return 1;
		else return -1; end if;
	end if;
end function;

///////////////////////////////////////
// Redei Symbol specific functions.
///////////////////////////////////////

//Constructs E/K_a/QQ and the beta needed for finding some choice of F
function EFieldBetaConstructor(a,b) 
	if IsSquare(a) or IsSquare(b) then return 9; end if;
	// if IsSquare(a) or IsSquare(b) or IsSquare(a/b) then return 9; end if;
	Ka:=QuadraticField(a);
	isnorm,beta:=NormEquation(Ka, RationalField()!b); //If you do not include the coercion then a=-19,b=5 fails -- magma being magma.
	beta:=beta[1];
	if not isnorm then print "UNDEFINED: a,b does not give a solution to x^2=ay^2+bz^2"; return 9; end if;
	P<x>:=PolynomialRing(RationalField());
	if IsSquare(Ka!b) then
		E:=Ka;
	else
		E:=ext<Ka|x^2-b:Abs:=true>;
	end if;
	return true, E, beta;
end function;

function Is2MinimallyRamified(F,K)
	OK:=RingOfIntegers(K);
	p2:=Factorisation(2*OK)[1][1];
	if not IsSubfield(K,F) then print "ERROR: Catastrophic failure, F was not extension of K."; end if;
	if BaseField(F) ne K then F:=RelativeField(AbsoluteField(K),AbsoluteField(F));end if;//Errors without doing this test first, Magma handles subfields badly!
	if Valuation(Conductor(AbelianExtension(F)), p2) eq 2 then return true; 
	else return false; end if;
end function;

function MinimallyRamifiedFConstructor(a,b,E,beta) // Uses 7.1 of Stevenhagen a lot!
	P<x>:=PolynomialRing(RationalField());
	if IsSquare(a/b) then 
		K:=RationalsAsNumberField();
	else
		K:=NumberField(x^2-a*b);
	end if;
	E:=AbsoluteField(E);
	P<x>:=PolynomialRing(E);
	F:=AbsoluteField(NumberField(x^2-beta));
	Delta_a:=QuadraticDiscriminant(a);
	Delta_b:=QuadraticDiscriminant(b);
	//p-ramification
	ramprimes:=[p: p in RamifiedRationalPrimes(F)| MyRamificationIndex(F,p) ne MyRamificationIndex(K,p)];
	oddprimes:=[p:p in ramprimes| IsOdd(p)];
	for p in oddprimes do
		if Delta_a mod p eq 0 and Delta_b mod p eq 0 then continue; // Forced Ramification
		else beta:=p*beta;end if;
	end for;
	F:=NumberField(x^2-beta);
	ramprimes:=[p: p in RamifiedRationalPrimes(F)| MyRamificationIndex(F,p) ne MyRamificationIndex(K,p)];//Redefining with new F to account for potentially introducing ramification at 2
	avoidablyramifiedat2:= 2 in ramprimes and (IsOdd(Delta_a) or IsOdd(Delta_b));//the second statement is for forced ramification
	//2-minimal ramification
	if avoidablyramifiedat2 then
		if (IsOdd(Delta_a) and IsOdd(Delta_b))//Lemma 7.1 part 1
		or (IsEven(Delta_a) and Delta_b mod 8 eq 1) //Lemma 7.1 part 3
		or (IsEven(Delta_b) and Delta_a mod 8 eq 1) 
		then
		for t in [-1,2,-2] do
			beta_t:=beta*t;
			F_t:=NumberField(x^2-beta_t);
			if not 2 in [p: p in RamifiedRationalPrimes(AbsoluteField(F_t))| MyRamificationIndex(F_t,p) ne MyRamificationIndex(K,p)] then
				beta:=beta_t;
				F:=NumberField(x^2-beta);
			end if;
		end for;
		elif not Is2MinimallyRamified(F,K) then 
			_,sqrta:=IsSquare(E!a);
			tau:=(1+sqrta)^2/2;
			beta:=tau*beta;
			F:=NumberField(x^2-beta);
		end if;
	end if;
	return true, F, beta;
end function;

//Obtains the ideal in K corresponding to the squarefree part of the integer c.
function GetCorrespondingIdeal(c, K)
	OK:=RingOfIntegers(K);
	C:=1*OK;
	for p in PrimeFactors(SquarefreeIntegerPart(c)) do
		C:=C*Factorisation(p*OK)[1][1];
	end for;
	return C;
end function;

//////////////////////////////////////////////
// The RedeiSymbol and minimally Ramified Field functions.
// requires input a,b,c to be integers!
//////////////////////////////////////////////
function RedeiSymbol(a,b,c: Additive:=false)
	T:=TestRedeiDefined(a,b,c);
	if T ne 0 then print ErrorCodesRedei(T); return 2; end if;
	if IsSquare(a) or IsSquare(b) or IsSquare(c) then return RtnValues(true,Additive); end if;
	_, E, beta:=EFieldBetaConstructor(a,b);
	_, F:=MinimallyRamifiedFConstructor(a,b,E,beta);
	P<x>:=PolynomialRing(RationalField());
	if IsSquare(a/b) then 
		K:=RationalsAsNumberField();
	else
		K:=NumberField(x^2-a*b);
	end if;
	if not IsSubfield(K,F) then print "ERROR: Catastrophic failure, F was not extension of K.";end if;
	if IsSquare(a/b) then // abelian extension doesnt like RationalsAsNumberField
		AbF:=AbelianExtension(RelativeField(RationalField(),F));
	else
		AbF:=AbelianExtension(RelativeField(K,F));
	end if;
	Art:=ArtinMap(AbF);
	FF:=NumberField(AbF);
	C:=GetCorrespondingIdeal(c, K);
	if c gt 0 then
		return RtnValues(Art(C)(FF.1) eq FF.1, Additive); // TODO: remove extra return data
	else
		return RtnValues(Art(C)(FF.1) eq ComplexConjugate(FF.1), Additive);
	end if;
end function;


//Implements the above in order to obtain the minimally ramified field $F$ from
//suitable a,b.
function GetMinimallyRamifiedF(a,b)
	if TestHilbert(a,b) ne 0 then
		print "ERROR:  The Hilbert Symbol criterion was not met.";
		return false;
	end if;
	if IsSquare(a) or IsSquare(b) then
		print "ERROR: One of the inputs is square.";
		return false;
	end if;
	_,E,beta:=EFieldBetaConstructor(a,b);
	_,F,beta:=MinimallyRamifiedFConstructor(a,b,E,beta);
	return F,E,beta;
end function;
