load "Redei.mg";

function GenerateTriples(LBound,UBound:WriteOut:=false,Overwrite:=false,IgnoreSquares:=false)
	ValidTriples:=[];
	for X in [LBound..UBound] do
		if X eq 0 or (IgnoreSquares and IsSquare(X)) then continue; end if;
		for Y in [LBound..UBound] do 
			if Y eq 0 or (IgnoreSquares and IsSquare(Y)) then continue; end if;
			for Z in [LBound..UBound] do
		 		if Z eq 0 or (IgnoreSquares and IsSquare(Z)) then continue; end if;
				if TestRedeiDefined(X,Y,Z) eq 0 then 
					Append(~ValidTriples, [X,Y,Z]);
				end if;
			end for;
		end for;
	end for;
	if WriteOut then 
		PrintFile("ValidRedeiTriples.mg", "ValidRedeiTriples:=":Overwrite:=Overwrite);
		PrintFile("ValidRedeiTriples.mg", ValidTriples);
		PrintFile("ValidRedeiTriples.mg", ";");
	end if;
	return ValidTriples;
end function;

function TestRedeiReciprocity(ValidTriples)
	for triple in ValidTriples do
		outcomevalues:={};
		for ordering in Permutations({1,2,3}) do
			outcomevalues:= outcomevalues join {RedeiSymbol(triple[ordering[1]],triple[ordering[2]],triple[ordering[3]])};
		end for;
		if #outcomevalues ne 1 then
			print "FAILURE: Redei symbol at\n", triple, "did not satisfy reciprocity.";
			print "DUMP:";
			print "Outcomevalues: ", outcomevalues;
			print "Triple: ", triple;
			return false;
		end if;
	end for;
	return true;
end function;


function TestSuiteRedei(lbound,ubound)
	ValidTriples:=GenerateTriples(lbound,ubound: WriteOut:=true,Overwrite:=true,IgnoreSquares:=true);
	print "Triples Generated: ", true;
	print "Checking Reciprocity...";
	print "Reciprocity: ", TestRedeiReciprocity(ValidTriples);
	return "Testing Complete.";
end function;
