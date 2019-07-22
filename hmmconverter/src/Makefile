HMMConverter: HMMConverter.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o define.o models.o hmm_state.o hmm.o sequence.o input_from_files.o blastn.o match.o parameters.o evaluation.o info.o scoring_functions.o transitionprob.o emissionprob.o model_parameters.o
	g++ -static -O3 -o HMMConverter HMMConverter.o tinyxml.o tinystr.o tinyxmlerror.o tinyxmlparser.o define.o models.o hmm_state.o hmm.o sequence.o input_from_files.o blastn.o match.o parameters.o evaluation.o info.o scoring_functions.o transitionprob.o emissionprob.o model_parameters.o
HMMConverter.o: HMMConverter.cpp tinyxml.h define.h models.h hmm_state.h multi_array_2.h input_from_files.h blastn.h match.h parameters.h transitionprob.h emissionprob.h model_parameters.h tube.h
	g++ -static -O3 -c -Wno-deprecated HMMConverter.cpp 
tinyxml.o: tinyxml.cpp tinyxml.h
	g++ -static -O3 -c tinyxml.cpp
tinystr.o: tinystr.cpp tinyxml.h
	g++ -static -O3 -c tinystr.cpp
tinyxmlerror.o: tinyxmlerror.cpp tinyxml.h
	g++ -static -O3 -c tinyxmlerror.cpp
tinyxmlparser.o: tinyxmlparser.cpp tinyxml.h
	g++ -static -O3 -c tinyxmlparser.cpp
define.o: define.cpp define.h multi_array_2.h tinyxml.h tube.h parameters.h Stack.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated define.cpp 
transitionprob.o: transitionprob.cpp transitionprob.h define.h tinyxml.h hmm.h model_parameters.h Stack.h
	g++ -static -O3 -c -Wno-deprecated transitionprob.cpp
emissionprob.o: emissionprob.cpp emissionprob.h define.h tinyxml.h hmm.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated emissionprob.cpp
model_parameters.o: model_parameters.cpp model_parameters.h define.h
	g++ -static -O3 -c -Wno-deprecated model_parameters.cpp
models.o: models.cpp models.h define.h multi_array_2.h hmm.h hmm_state.h Stack.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated models.cpp
hmm_state.o: hmm_state.cpp hmm_state.h multi_array_2.h define.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated hmm_state.cpp
hmm.o : hmm.cpp hmm.h hmm_state.h sequence.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated hmm.cpp
sequence.o : sequence.cpp sequence.h define.h match.h evaluation.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated sequence.cpp
input_from_files.o : input_from_files.cpp input_from_files.h define.h multi_array_2.h model_parameters.h sequence.h blastn.h tube.h
	g++ -static -O3 -c -Wno-deprecated input_from_files.cpp
blastn.o : blastn.cpp blastn.h multi_array_2.h define.h sequence.h tube.h match.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated blastn.cpp
match.o : match.cpp match.h
	g++ -static -O3 -c -Wno-deprecated match.cpp
parameters.o : parameters.cpp parameters.h
	g++ -static -O3 -c -Wno-deprecated parameters.cpp
evaluation.o : evaluation.cpp evaluation.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated evaluation.cpp
info.o : info.cpp info.h model_parameters.h
	g++ -static -O3 -c -Wno-deprecated info.cpp
scoring_functions.o : scoring_functions.cpp scoring_functions.h
	g++ -static -O3 -c -Wno-deprecated scoring_functions.cpp