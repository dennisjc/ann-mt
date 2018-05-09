all:
	make evalFitness
	make saveImage
	make showImage
	make testBlobInfluence
	make testRender3d
	make saveModel
	make evalFitness3D
	make showParams3D

evalFitness: evalFitness.cpp functions.cpp functions.h
	g++ -O3 -Wall -o evalFitness evalFitness.cpp functions.cpp `pkg-config --cflags --libs opencv`

saveImage: saveImage.cpp functions.cpp functions.h
	g++ -O3 -Wall -o saveImage saveImage.cpp functions.cpp `pkg-config --cflags --libs opencv`

showImage: showImage.cpp functions.cpp functions.h
	g++ -O3 -Wall -o showImage showImage.cpp functions.cpp `pkg-config --cflags --libs opencv`

testBlobInfluence: testBlobInfluence.cpp functions.cpp functions.h
	g++ -O3 -Wall -o testBlobInfluence testBlobInfluence.cpp functions.cpp `pkg-config --cflags --libs opencv`

testRender3d: testRender3d.cpp functions.cpp functions.h
	g++ -O3 -Wall -o testRender3d testRender3d.cpp functions.cpp  `pkg-config --cflags --libs opencv`

saveModel: saveModel.cpp functions.cpp functions.h
	g++ -O3 -Wall -o saveModel saveModel.cpp functions.cpp  `pkg-config --cflags --libs opencv`

evalFitness3D: evalFitness3D.cpp functions.cpp functions.h
	g++ -O3 -Wall  -o evalFitness3D evalFitness3D.cpp functions.cpp  `pkg-config --cflags --libs opencv`

showParams3D: showParams3D.cpp functions.cpp functions.h
	g++ -O3 -Wall -o  showParams3D showParams3D.cpp functions.cpp `pkg-config --cflags --libs opencv`
